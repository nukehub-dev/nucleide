//! Error type for the spectroscopy crate.

use thiserror::Error;

/// Errors raised while smoothing spectra, counting peaks, calibrating
/// energy/efficiency, parsing `.spe` text, or evaluating X-ray lines.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum Error {
    /// Rectangular smoothing width is below 3 (got `{0}`).
    #[error("spectroscopy: smoothing width {0} is less than 3")]
    SmoothWidthTooSmall(usize),
    /// Rectangular smoothing width is even (got `{0}`); it must be odd.
    #[error("spectroscopy: smoothing width {0} is not odd")]
    SmoothWidthEven(usize),
    /// No counts supplied to a smoothing pass.
    #[error("spectroscopy: no counts to smooth")]
    EmptyCounts,
    /// Five-point smoothing needs at least 4 channels (got `{0}`).
    #[error("spectroscopy: five-point smoothing needs at least 4 channels, got {0}")]
    TooFewChannels(usize),
    /// `c1` must not exceed `c2` (got `c1={c1}`, `c2={c2}`).
    #[error("spectroscopy: c1 ({c1}) must be less than c2 ({c2})")]
    BadChannelRange {
        /// Lower channel bound.
        c1: i64,
        /// Upper channel bound.
        c2: i64,
    },
    /// `c1` must be non-negative (got `{0}`).
    #[error("spectroscopy: c1 ({0}) must be positive number above 0")]
    NegativeChannel(i64),
    /// `c2` exceeds the largest channel label (`{c2}` > `{max}`).
    #[error("spectroscopy: c2 ({c2}) must be less than max number of channels ({max})")]
    ChannelOutOfRange {
        /// Requested upper bound.
        c2: i64,
        /// Largest channel label.
        max: f64,
    },
    /// No channel labels are available for the range check.
    #[error("spectroscopy: no channels available for range check")]
    EmptyChannels,
    /// Background method id is not 1 (got `{0}`); only `m == 1` exists.
    #[error("spectroscopy: m ({0}) is not set to a valid method id (only m == 1 exists)")]
    BadBackgroundMethod(i64),
    /// Detector-efficiency fit selector is not 1 or 2 (got `{0}`).
    #[error("spectroscopy: the selected eff_fit ({0}) is not valid")]
    UnknownEfficiencyFit(i64),
    /// No efficiency coefficients supplied.
    #[error("spectroscopy: no efficiency coefficients supplied")]
    EmptyCoefficients,
    /// Energy calibration needs 3 fit coefficients (got `{0}`).
    #[error("spectroscopy: energy calibration needs 3 fit coefficients, got {0}")]
    MissingCalibration(usize),
    /// `.spe` text is in the other format (magic-line mismatch): `{0}`.
    #[error("spectroscopy: spe file format not supported by this function: {0}")]
    UnsupportedFormat(&'static str),
    /// A required `$TAG:` line is absent from dollar-format text: `{0}`.
    #[error("spectroscopy: missing required tag {0}")]
    MissingTag(&'static str),
    /// A header or data value does not parse: `{0}`.
    #[error("spectroscopy: malformed spe value: {0}")]
    MalformedValue(&'static str),
    /// Dollar-format text ends before all declared data lines are read.
    #[error("spectroscopy: truncated spe data section")]
    TruncatedData,
}
