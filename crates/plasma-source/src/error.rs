//! Error type for the plasma-source crate.

use thiserror::Error;

/// Result alias for the `plasma-source` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while validating source configurations, tabulating spectra,
/// or rendering source cards.
///
/// Out-of-scope requests (toroidal sectors, the T-T and D(d,p)T branches,
/// reactant distributions beyond the shared-temperature Maxwellian mixture)
/// are reported through [`Error::NotYetSupported`] — a loud named error,
/// never a panic or a silent fallback. Invalid fuel-mixture fractions carry
/// their own named errors ([`Error::NonFinite`], [`Error::InvalidFuelMixture`]).

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
    /// (the D-D η factor goes non-positive around 300–4700 keV, far above
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
    /// The emitted card failed to re-parse through the typed `SDEF` reader
    /// (an internal emission invariant; surfaced loudly rather than
    /// delivered unverified).
    #[error("plasma-source: emitted card failed the SDEF reader round trip: {0}")]
    CardRoundTrip(String),
    /// The requested capability is outside the v1 scope: `{0}`.
    #[error("plasma-source: not yet supported: {0}")]
    NotYetSupported(&'static str),
}
