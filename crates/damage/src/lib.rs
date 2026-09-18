#![warn(missing_docs)]
//! Damage and gas-production metrics for irradiated materials: the NRT-dpa
//! and arc-dpa displacement functions, He/H production in appm, and He/dpa
//! ratios — all by spectral folding of a caller-supplied multigroup flux
//! with caller-supplied damage/gas response functions.
//!
//! The crate has three layers, all pure functions over caller inputs:
//!
//! - [`physics`] — the closed-form primary-damage physics: the Lindhard
//!   damage-energy partition (Robinson fit), the NRT displacement function
//!   (Norgett, Robinson & Torrens, Nucl. Eng. Des. 33 (1975) 50–54), and
//!   the arc-dpa efficiency correction (Nordlund et al., Nat. Commun. 9
//!   (2018) 1084, CC BY 4.0). Material constants (`E_d`, `b_arc`, `c_arc`)
//!   and nuclide keys are caller-supplied; the only vendored numbers are
//!   the SPECTER Table VII fallback (see below).
//! - [`fold`] — the multigroup folds `nrt_dpa` / `arc_dpa` / `gas_appm` /
//!   `he_dpa_ratio` over `(flux, response, bounds)` slices with
//!   piecewise-constant-per-group semantics (documented in the module).
//! - [`coil`] — coil fast-fluence and lifetime bookkeeping: the fast-flux
//!   sum above a caller threshold, history accumulation, and the
//!   weakest-link life over caller-supplied limit tables (dpa spectral
//!   weighting stays in [`fold`], consumed as rates here).
//!
//! [`uq`] propagates caller-block uncertainty through the folds with the
//! landed `linalg` MVN machinery (pinned seed, `k`-standard-error gates) —
//! no new sampling code.
//!
//! Provenance stance: every equation is a published fact implemented
//! clean-room; the SPECTER report (Greenwood & Smither, ANL/FPP/TM-197,
//! 1985 — US-government public domain) is the validation oracle
//! (`validation/damage_vs_specter.py` in the repo harness — a handful of
//! transcribed output spots checked at report print precision) *and* the
//! source of one opt-in vendored fallback: [`specter::SpecterTable`]
//! (Table VII spectrum-averaged damage-energy cross sections for 24
//! elements plus the Table II `E_d` column, displacement XS only). The
//! fallback is never consulted implicitly — the [`fold`] kernels only see
//! caller slices. ASTM E693/E521 are paywalled standards referenced by
//! designation string only, never transcribed. PKA-spectra solving (the
//! GPL-3.0 fispact-org PKA evaluator), transport solving, and
//! evaluated-data vendoring (TENDL/JEFF/EAF/IRDFF) are out of scope —
//! callers bring their own response tables.
//!
//! Layering: depends on `nucleide-linalg` and `nucleide-nuclei` only;
//! bindings depend on this crate, never the reverse.

pub mod coil;
pub mod error;
pub mod fold;
pub mod physics;
pub mod specter;
pub mod uq;

pub use coil::{
    accumulate, coil_lifetime, coil_remaining, fast_fluence, fast_flux, CoilLifetime, FPY_SECONDS,
};

pub use error::{Error, Result};
pub use fold::{arc_dpa, gas_appm, he_dpa_ratio, nrt_dpa, BARNS_TO_CM2};
pub use physics::{
    arc_displacements, arc_displacements_for, arc_efficiency, damage_energy, lindhard_partition,
    nrt_displacements, nrt_displacements_for, ArcParams, LINDHARD_E_COEFF, LINDHARD_K_COEFF,
    NRT_EFFICIENCY, ROBINSON_G1, ROBINSON_G16, ROBINSON_G34,
};
pub use specter::{specter_ed_ev, SpecterEntry, SpecterSpectrum, SpecterTable};
pub use uq::{fold_uq, FoldMetric, UqSummary};
