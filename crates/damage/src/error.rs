//! Error type for the damage crate.

use thiserror::Error;

/// Result alias for the `damage` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while evaluating damage functions or folding spectra.
///
/// Every rejection names its cause; out-of-scope capabilities (PKA-spectra
/// solving, further displacement tables beyond the vendored SPECTER Table
/// VII fallback, UQ on nonlinear ratios) report through
/// [`Error::NotYetSupported`] — a loud named error, never a panic or
/// a silent fallback.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    /// A required slice is empty.
    #[error("damage: input is empty")]
    Empty,
    /// A length does not match the problem dimension: `{what}` expected
    /// `{expected}`, got `{got}`.
    #[error("damage: {what} length mismatch: expected {expected}, got {got}")]
    DimensionMismatch {
        /// Name of the offending input.
        what: &'static str,
        /// Length required by the leading input.
        expected: usize,
        /// Length actually supplied.
        got: usize,
    },
    /// A non-finite value in the named input (`"flux"`, `"response"`,
    /// `"bounds"`, `"seconds"`, `"pka_energy_ev"`, `"ed_ev"`, `"t_dam_ev"`,
    /// `"b_arc"`, `"c_arc"`, `"k"`, `"mean_delta"`, `"covariance"`).
    #[error("damage: non-finite value in {0}")]
    NonFinite(&'static str),
    /// A physical quantity that must be non-negative holds a negative value.
    #[error("damage: {0} must be non-negative")]
    Negative(&'static str),
    /// A physical quantity that must be strictly positive is not
    /// (`"seconds"`, `"ed_ev"`, `"t_dam_ev"`).
    #[error("damage: {0} must be strictly positive")]
    NonPositive(&'static str),
    /// Group bounds are not strictly increasing: `bounds[{index}]` is not
    /// below `bounds[{index} + 1]`.
    #[error("damage: group bounds must be strictly increasing (violated at index {index})")]
    NonMonotonicBounds {
        /// Index of the lower boundary that fails `bounds[i] < bounds[i+1]`.
        index: usize,
    },
    /// A nucid does not decompose into a valid `(Z, A, state)` triple.
    #[error("damage: invalid nucid {0} (expected a valid Z/A/state packing)")]
    BadNuclide(
        /// Offending raw nucid integer.
        u32,
    ),
    /// The arc-dpa efficiency constants are outside the fitted domain:
    /// `b_arc` must be finite and negative and `c_arc` finite in `(0, 1)`.
    #[error("damage: invalid arc-dpa constants: b_arc={b_arc}, c_arc={c_arc} (need b_arc < 0 and 0 < c_arc < 1)")]
    InvalidArcParams {
        /// Power-law exponent supplied.
        b_arc: f64,
        /// Saturation efficiency supplied.
        c_arc: f64,
    },
    /// The He/dpa ratio was requested where the damage fold is exactly zero
    /// (a zero-flux spectrum, or all-zero damage response). Named error
    /// instead of `f64::INFINITY` — the `ZeroMaxFlux` precedent in
    /// `vr-tools` `magic`.
    #[error("damage: He/dpa ratio at zero dpa is undefined (zero flux or zero damage response)")]
    ZeroDpa,
    /// A line of an embedded data transcription does not parse: TSV line
    /// `line` reads `{msg}`. The committed transcriptions are pinned by
    /// row-count tests, so this fires only on a corrupt edit.
    #[error("damage: TSV line {line}: {msg}")]
    Parse {
        /// 1-based TSV line number of the offending line.
        line: usize,
        /// What was wrong with the line.
        msg: String,
    },
    /// An element symbol outside the vendored SPECTER Table VII set:
    /// `{element}` (the table covers 24 elements, Be through Pb).
    #[error(
        "damage: unknown SPECTER table element `{element}` (vendored Table VII covers 24 elements)"
    )]
    UnknownSpecterElement {
        /// The unrecognized element symbol.
        element: String,
    },
    /// A spectrum name outside the vendored SPECTER Table VII set:
    /// `{spectrum}` (one of `thermal`, `fission`, `14mev`, `hfir`,
    /// `ebr2`, `fftf`, `fusion`).
    #[error("damage: unknown SPECTER table spectrum `{spectrum}` (expected one of thermal, fission, 14mev, hfir, ebr2, fftf, fusion)")]
    UnknownSpecterSpectrum {
        /// The unrecognized spectrum name.
        spectrum: String,
    },
    /// The requested capability is outside the v1 scope: `{0}`.
    #[error("damage: not yet supported: {0}")]
    NotYetSupported(&'static str),
    /// MVN sampling over caller blocks failed (forwarded from `linalg`).
    #[error("damage: {0}")]
    Sampling(
        /// Underlying sampling error.
        #[from]
        nucleide_linalg::SampleError,
    ),
}
