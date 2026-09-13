#![warn(missing_docs)]
//! Gamma-ray spectroscopy and measurement toolkit.
//!
//! Pure algorithms over caller-supplied spectrum vectors: rectangular and
//! five-point smoothing (E1–E2), background/gross/net counting (E3–E5),
//! quadratic energy calibration (E6), log-polynomial detector-efficiency
//! evaluation (E7) plus the E7-fit log-space coefficient fit, caller-constant
//! X-ray algebra (E8), the two text `.spe` readers (dollar and plain formats),
//! caller-line SDEF decay-source cards (E9), and the runtime decay-lines TSV
//! interchange feeding E9.
//!
//! No tabulated data lives here: atomic constants, calibration fits, decay
//! lines, and efficiency-fit points/weights are all caller inputs (fitted
//! E7 coefficients come out of [`fit_efficiency`]; every other coefficient
//! vector stays an input). No I/O beyond parsing `&str` (callers read files
//! themselves).
//!
//! # Equation labels
//!
//! `E1`–`E9` are pinned in `docs/theory/spectroscopy.mdx` (the
//! coefficient fit is the `E7-fit` extension, not a new equation number);
//! code comments cite those labels.

pub mod calib;
pub mod counts;
pub mod error;
pub mod lines;
pub mod sdef;
pub mod smooth;
pub mod spe;
pub mod spectrum;
pub mod xray;

pub use calib::{detector_efficiency, energy_bins, fit_efficiency};
pub use counts::{calc_bg, gross_count, net_counts};
pub use error::Error;
pub use lines::parse_lines_tsv;
pub use sdef::{normalize_decay_lines, sdef_card, PointSource};
pub use smooth::{five_point_smooth, rect_smooth};
pub use spe::{parse_dollar_spe, parse_plain_spe};
pub use spectrum::{GammaSpectrum, Spectrum};
pub use xray::{xray_lines, AtomicData, XrayLine};
