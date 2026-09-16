#![warn(missing_docs)]
//! Tokamak fusion-neutron sources: ring and point geometry (first landing)
//! plus a parametric Miller-geometry plasma (second landing) over the D-D
//! (2.45 MeV) and D-T (14.1 MeV) reactions, sampled to particle vectors and
//! emitted as MCNP `SDEF` / Serpent `src` cards with a drift report.
//!
//! Ring/point sources ([`source`], [`sample`], [`emit_sdef`],
//! [`emit_serpent`]) are the simple full-torus models; the parametric
//! tokamak plasma ([`miller`], [`profile`], [`reactivity`], [`parametric`])
//! takes caller-supplied density/temperature profiles over Miller flux
//! surfaces (Fausser et al., Fus. Eng. Des. **87** (2012) 787) with
//! reactivity-weighted emission (Bosch & Hale, Nucl. Fusion **32** (1992)
//! 611) and Ballabio-broadened spectra. Profiles are caller inputs — nothing
//! computes profiles and no equilibrium is solved. What stays out (loud
//! [`Error::NotYetSupported`], never a guess): reactant distributions beyond
//! the shared-temperature Maxwellian D/T mixture ([`parametric`] — the
//! per-species-temperature Eriksson et al., Comput. Phys. Commun. **199**
//! (2016) 40 generalization), the T-T and D(d,p)T branches, and toroidal
//! sectors.
//!
//! # Physics
//!
//! Thermonuclear neutron spectra are Gaussian with a width proportional to
//! `sqrt(T_i)` (Brysk, Plasma Phys. **15** (1973) 611). [`FusionReaction`]
//! implements the closed-form coefficient fits of Ballabio et al., Nucl.
//! Fusion **38** (1998) 1723, Table III (mean shift and weakly
//! temperature-dependent FWHM), which refine the Brysk scaling; at
//! `T_i = 0` the spectrum is the monoenergetic nominal line. Ring/point
//! spatial moments are closed form (radius/height/azimuth); the parametric
//! source's Miller-map Jacobian and volume weighting carry their own closed
//! forms and finite-difference gates in [`miller`]. The analytic gates in
//! `validation/plasma_source_vs_openmc.py` check sampled particle vectors
//! against both.
//!
//! # Units and dialects
//!
//! Lengths are centimetres, energies MeV, ion temperature keV — the
//! transport-code card convention. MCNP `SDEF` cards render through the
//! typed `nucleide-mcnp-io` reader and round-trip byte-identically (the
//! spectroscopy E9 precedent); Serpent `src` cards follow the public Serpent
//! input-manual spelling and are analytic-by-design in the drift report.
//!
//! # Layering
//!
//! Dependencies are `nucleide-mcnp-io` (the `SDEF` dialect) and
//! `nucleide-nuclei` (`PAR=` designators) only — no `mcpl-io`. MCPL
//! projection stays caller-side, the same rule as `vr-tools` `KdeSampler`:
//! this crate outputs particle vectors and card strings; the caller writes
//! files with `nucleide-mcpl-io` / `nucleide-mcnp-io` when it wants them.
//!
//! # Example
//!
//! ```rust
//! use nucleide_plasma_source::{
//!     emit_sdef, emit_serpent, FusionReaction, SourceSampler, PlasmaSourceConfig,
//! };
//!
//! // D-T tokamak ring, R = 300 cm at the midplane, T_i = 20 keV.
//! let config = PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, 20.0);
//!
//! // Deterministic particle vector (seeded); MCPL stays caller-side.
//! let particles = SourceSampler::new(config, 42).unwrap().sample_n(1000);
//! assert_eq!(particles.len(), 1000);
//!
//! // Cards + drift report.
//! let sdef = emit_sdef(&config, 5, 21).unwrap();
//! sdef.verify_round_trip().unwrap();
//! let serpent = emit_serpent(&config, 21).unwrap();
//! assert!(sdef.text.contains("RAD=D1"));
//! assert!(serpent.text.contains("src 1 rad d1"));
//! ```
//!
//! The parametric plasma builds the same way via
//! [`ParametricPlasmaConfig`] + [`ParametricSampler`] /
//! `emit_sdef_parametric` / `emit_serpent_parametric` (see [`parametric`]).

pub mod emit_sdef;
pub mod emit_serpent;
pub mod error;
pub mod miller;
pub mod parametric;
pub mod profile;
pub mod reaction;
pub mod reactivity;
pub mod report;
pub mod sample;
pub mod source;
pub mod spectrum;

pub use emit_sdef::EmittedCard;
pub use emit_sdef::{emit_sdef, emit_sdef_parametric};
pub use emit_serpent::{emit_serpent, emit_serpent_parametric};
pub use error::{Error, Result};
pub use miller::MillerGeometry;
pub use parametric::{
    BinnedDistribution, EmissionHistograms, FuelMixture, ParametricPlasmaConfig, ParametricSampler,
};
pub use profile::{DensityProfile, ProfileMode, TemperatureProfile};
pub use reaction::FusionReaction;
pub use report::{DriftReport, DriftRow};
pub use sample::{Particle, SourceSampler};
pub use source::{PlasmaSourceConfig, PointSource, RingSource, SourceModel};
pub use spectrum::{SpectrumSpec, SpectrumTable};
