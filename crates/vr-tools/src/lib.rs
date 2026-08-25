//! Monte Carlo variance-reduction utilities built on [`mcnp_io`] meshtal data.
//!
//! - [`magic`] — MAGIC weight-window generation operating on native
//!   [`mcnp_io::meshtal::MeshTallyData`] instead of MOAB-tagged meshes.
//! - [`sampling`] — Walker/Vose alias-table source sampling plus a
//!   voxel-level `MeshSourceSampler` with ANALOG / UNIFORM / USER bias modes.

pub mod magic;
pub mod sampling;

pub use magic::{magic, magic_with, MagicOutput, MagicParams, MagicSelection};
pub use sampling::{AliasTable, MeshSourceSampler, Mode, SampledVoxel};

/// Errors raised by variance-reduction tools.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Tally carries no volume elements.
    EmptyTally,
    /// A requested array has the wrong length.
    LengthMismatch { expected: usize, got: usize },
    /// Every flux value feeding one energy bin is non-positive, so the
    /// MAGIC normalization `value / (2 * max)` would divide by zero.
    /// (Rather than silently emitting `inf`/`nan`, this is an error.)
    ZeroMaxFlux { energy_group: usize },
    /// PDF input is empty.
    EmptyPdf,
    /// PDF contains a negative entry.
    NegativePdf { index: usize, value: f64 },
    /// PDF contains a non-finite (NaN/infinite) entry.
    NonFinitePdf { index: usize },
    /// PDF sums to zero (or negatively); cannot normalize.
    ZeroSumPdf,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::EmptyTally => write!(f, "tally contains no volume elements"),
            Error::LengthMismatch { expected, got } => {
                write!(f, "length mismatch: expected {expected}, got {got}")
            }
            Error::ZeroMaxFlux { energy_group } => write!(
                f,
                "energy group {energy_group} has no positive flux; \
                 weight-window normalization would divide by zero"
            ),
            Error::EmptyPdf => write!(f, "pdf must contain at least one value"),
            Error::NegativePdf { index, value } => {
                write!(f, "pdf[{index}] = {value} is negative")
            }
            Error::NonFinitePdf { index } => write!(f, "pdf[{index}] is not finite"),
            Error::ZeroSumPdf => write!(f, "pdf sums to zero; cannot normalize"),
        }
    }
}

impl std::error::Error for Error {}
