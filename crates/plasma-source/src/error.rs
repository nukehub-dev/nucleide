//! Error type for the plasma-source crate.

use thiserror::Error;

/// Result alias for the `plasma-source` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while validating source configurations, tabulating spectra,
/// or rendering source cards.
///
/// Out-of-scope requests (T-T neutron transport, proton transport, reactant
/// distributions beyond a single deuterium hot-tail fraction on a D/T fuel
/// mixture, per-species ion
/// temperatures without a D/T fuel mixture) are reported through
/// [`Error::NotYetSupported`] — a loud named error, never a
/// panic or a silent fallback. The D(d,p)T proton *rate* is not out of
/// scope: it is accounted alongside the neutron source
/// (`proton_strength_density` / `total_proton_strength`, 50/50 with the D-D
/// neutron branch) while the sampler and the cards stay neutron-only.
/// Invalid fuel-mixture fractions, invalid deuterium-tail parameters,
/// invalid lattice clouds and symmetry folds, and invalid
/// toroidal-sector angles carry their own named errors ([`Error::NonFinite`],
/// [`Error::InvalidFuelMixture`], [`Error::InvalidTail`],
/// [`Error::InvalidLattice`], [`Error::InvalidSymmetry`],
/// [`Error::InvalidSector`]); out-of-range
/// species temperatures reuse [`Error::NonFinite`] and
/// [`Error::NegativeIonTemperature`].

#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    /// A source-model field is `NaN` or infinite: `{0}`.
    #[error("plasma-source: non-finite value in {0}")]
    NonFinite(&'static str),
    /// The ring radius must be finite and strictly positive; got `{0}` cm.
    #[error("plasma-source: ring radius must be > 0 cm, got {0}")]
    NonPositiveRadius(f64),
    /// The ion temperature must be finite and non-negative; got `{0}` keV.
    #[error("plasma-source: ion temperature must be >= 0 keV, got {0}")]
    NegativeIonTemperature(f64),
    /// The reactivity fit leaves its validity domain at this temperature
    /// (the D-D η factor goes non-positive around 965–2720 keV, far above
    /// the published fit range): `{reaction}` at `{ti_kev}` keV.
    #[error("plasma-source: {reaction} reactivity fit out of domain at {ti_kev} keV")]
    FitOutOfDomain {
        /// Reaction whose fit diverged.
        reaction: &'static str,
        /// Ion temperature in keV where the fit broke down.
        ti_kev: f64,
    },
    /// The particle weight must be finite and strictly positive; got `{0}`.
    #[error("plasma-source: particle weight must be > 0, got {0}")]
    NonPositiveWeight(f64),
    /// A spectrum-tabulation request is invalid: `{0}`.
    #[error("plasma-source: invalid spectrum tabulation: {0}")]
    BadTabulation(&'static str),
    /// The MCNP designator version must be 5 or 6; got `{0}`.
    #[error("plasma-source: unsupported MCNP version {0} (supported: 5, 6)")]
    UnsupportedMcnpVersion(u32),
    /// A Miller-geometry field fails its documented range check: `{0}`.
    #[error("plasma-source: invalid Miller geometry: {0}")]
    InvalidGeometry(&'static str),
    /// A profile parameter fails its documented range check: `{0}`.
    #[error("plasma-source: invalid plasma profile: {0}")]
    InvalidProfile(&'static str),
    /// A fuel-mixture fraction is negative, or the pair does not sum to 1
    /// (within the documented 1e-12 tolerance): `{0}`.
    #[error("plasma-source: invalid fuel mixture: {0}")]
    InvalidFuelMixture(&'static str),
    /// A deuterium hot-tail parameter is out of range (fraction outside
    /// `[0, 1]`): `{0}`.
    #[error("plasma-source: invalid deuterium tail: {0}")]
    InvalidTail(&'static str),
    /// A toroidal-sector angle fails its documented range check
    /// (`start_angle` outside `[0, 2π)`, `rotation_angle` outside `(0, 2π]`):
    /// `{0}`.
    #[error("plasma-source: invalid toroidal sector: {0}")]
    InvalidSector(&'static str),
    /// A lattice birth-rate cloud fails its documented range check (empty
    /// cloud, negative point rate, zero total rate): `{0}`.
    #[error("plasma-source: invalid lattice source: {0}")]
    InvalidLattice(&'static str),
    /// A lattice symmetry fold fails its documented precondition (zero field
    /// periods, a point count that is not a multiple of the period count, or
    /// mismatched ion temperatures across a fold orbit): `{0}`.
    #[error("plasma-source: invalid lattice symmetry: {0}")]
    InvalidSymmetry(&'static str),
    /// The emitted card failed to re-parse through the typed `SDEF` reader
    /// (an internal emission invariant; surfaced loudly rather than
    /// delivered unverified).
    #[error("plasma-source: emitted card failed the SDEF reader round trip: {0}")]
    CardRoundTrip(String),
    /// The requested capability is outside the v1 scope: `{0}`.
    #[error("plasma-source: not yet supported: {0}")]
    NotYetSupported(&'static str),
}
