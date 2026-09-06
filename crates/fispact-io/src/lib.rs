//! FISPACT-II output parser producing ALARA-compatible response frames.
//!
//! Reads FISPACT-II `.fis` inventory tables into
//! [`alara_io::output::ResponseFrame`] rows so ALARA and FISPACT-II results
//! share one analysis shape. No activation solving is performed here.

pub mod error;
pub mod inventory;

pub use error::{Error, Result};
pub use inventory::{is_fispact_output, parse_to_frame};
