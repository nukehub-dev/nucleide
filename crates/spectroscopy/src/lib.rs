#![warn(missing_docs)]
//! Gamma-ray spectroscopy and measurement toolkit.
//!
//! Pure algorithms over caller-supplied spectrum vectors: rectangular and
//! five-point smoothing (E1–E2), background/gross/net counting (E3–E5),
//! quadratic energy calibration (E6), log-polynomial detector-efficiency
//! evaluation (E7), caller-constant X-ray algebra (E8), and the two text
//! `.spe` readers (dollar and plain formats).
//!
//! No tabulated data lives here: atomic constants, calibration fits, and
//! efficiency coefficients are all caller inputs. No I/O beyond parsing
//! `&str` (callers read files themselves).
//!
//! # Equation labels
//!
//! `E1`–`E8` are pinned in `docs/theory/spectroscopy.mdx`; code comments
//! cite those labels.

pub mod calib;
pub mod counts;
pub mod error;
pub mod smooth;
pub mod spe;
pub mod spectrum;
pub mod xray;

pub use calib::{detector_efficiency, energy_bins};
pub use counts::{calc_bg, gross_count, net_counts};
pub use error::Error;
pub use smooth::{five_point_smooth, rect_smooth};
pub use spe::{parse_dollar_spe, parse_plain_spe};
pub use spectrum::{GammaSpectrum, Spectrum};
pub use xray::{xray_lines, AtomicData, XrayLine};
