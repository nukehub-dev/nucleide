#![warn(missing_docs)]
//! Monte Carlo variance-reduction utilities built on [`mcnp_io`] meshtal data.
//!
//! - [`magic`] — MAGIC weight-window generation operating on native
//!   [`nucleide_mcnp_io::meshtal::MeshTallyData`] instead of MOAB-tagged meshes.
//! - [`windows`] — emission of MAGIC weight windows for OpenMC
//!   (`settings.xml` `<mesh>` + `<weight_windows>` text) and Serpent
//!   (`wwin ... wf FILE 2`, the MCNP WWINP spelling Serpent reads).
//! - [`sampling`] — Walker/Vose alias-table source sampling plus a
//!   voxel-level `MeshSourceSampler` with ANALOG / UNIFORM / USER bias modes.
//! - [`kde`] — Gaussian kernel-density source sampling (KDSource-class,
//!   clean-room): fit over caller particle vectors, deterministic resampling.

pub mod kde;
pub mod magic;
pub mod sampling;
pub mod windows;

pub use kde::{Bandwidth, KdeSampler};

pub use magic::{magic, magic_with, MagicOutput, MagicParams, MagicSelection};
pub use sampling::{AliasTable, MeshSourceSampler, Mode, SampledVoxel};
pub use windows::{
    emit_openmc_weight_windows, emit_serpent_wwin, OpenMcOptions, OpenMcWeightWindows, SerpentWwin,
};

/// Result alias for the `vr-tools` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised by variance-reduction tools.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// Tally carries no volume elements.
    EmptyTally,
    /// A requested array has the wrong length.
    LengthMismatch {
        /// Expected element count.
        expected: usize,
        /// Actual element count.
        got: usize,
    },
    /// Every flux value feeding one energy bin is non-positive, so the
    /// MAGIC normalization `value / (2 * max)` would divide by zero.
    /// (Rather than silently emitting `inf`/`nan`, this is an error.)
    ZeroMaxFlux {
        /// Index of the energy bin with no positive flux.
        energy_group: usize,
    },
    /// PDF input is empty.
    EmptyPdf,
    /// PDF contains a negative entry.
    NegativePdf {
        /// Index of the negative entry.
        index: usize,
        /// The offending value.
        value: f64,
    },
    /// PDF contains a non-finite (NaN/infinite) entry.
    NonFinitePdf {
        /// Index of the non-finite entry.
        index: usize,
    },
    /// PDF sums to zero (or negatively); cannot normalize.
    ZeroSumPdf,
    /// Tally or user density contains a negative value.
    NegativeTally {
        /// Index of the negative value.
        index: usize,
        /// The offending value.
        value: f64,
    },
    /// Tally array contains a non-finite (NaN/infinite) value.
    NonFiniteTally {
        /// Name of the tally field holding the value.
        field: &'static str,
        /// Index of the non-finite entry.
        index: usize,
    },
    /// A KDE Silverman dimension has zero variance (a zero bandwidth is a
    /// delta spike, never a density); use an explicit fixed width instead.
    ZeroVarianceDim {
        /// Index of the zero-variance dimension.
        dim: usize,
    },
    /// A KDE adaptive kNN bandwidth has an invalid neighbour rank or
    /// scale: `k` must satisfy `1 <= k <= n - 1` for `n` samples, and
    /// `scale` must be finite and positive.
    BadAdaptiveKnn {
        /// Fitted sample count.
        n: usize,
        /// Caller-supplied 1-based neighbour rank.
        k: usize,
        /// Caller-supplied scale factor.
        scale: f64,
    },
    /// A KDE adaptive kNN sample has a zero `k`-th-neighbour distance
    /// (duplicate points): a zero bandwidth is a delta spike, never a
    /// density. Mirrors [`Error::ZeroVarianceDim`].
    ZeroAdaptiveBandwidth {
        /// Index of the sample with zero `dist_k`.
        sample: usize,
    },
    /// A KDE draw uniform falls outside `[0, 1)`.
    BadDraw {
        /// The offending uniform value.
        value: f64,
    },
    /// A weight-window lower bound is negative and has no spelling in either
    /// target format (non-positive windows are inert; negative never is).
    NegativeWindow {
        /// Flat index into `MagicOutput::lower_bounds_ww`.
        index: usize,
        /// The offending value.
        value: f64,
    },
    /// A weight-window lower bound is NaN or infinite.
    NonFiniteWindow {
        /// Flat index into `MagicOutput::lower_bounds_ww`.
        index: usize,
    },
    /// Energy upper bounds are empty, non-finite, non-positive, or not
    /// strictly increasing, so the window energy grid is ill-defined.
    BadEnergyBounds {
        /// Index of the offending bound.
        index: usize,
    },
    /// Mesh bounds have fewer than two entries, are non-finite, or are not
    /// strictly increasing.
    BadMeshBounds {
        /// Axis index (0 = x, 1 = y, 2 = z).
        axis: usize,
        /// Index of the offending bound.
        index: usize,
    },
    /// An emission tuning parameter or card token falls outside the range the
    /// target code enforces when reading the file, or cannot be interpolated
    /// into the emitted card safely.
    BadEmissionOption {
        /// Parameter name (`upper_bound_ratio`, `survival_ratio`,
        /// `max_split`, `weight_cutoff`, `name`, `file`, `energy_bounds`,
        /// `upper_ww_bounds`).
        option: &'static str,
        /// The offending value.
        value: String,
        /// Why it was rejected.
        detail: &'static str,
    },
    /// The WWINP writer in `nucleide-mcnp-io` refused the assembled file.
    Wwinp(String),
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
            Error::NegativeTally { index, value } => write!(
                f,
                "tally/density value at index {index} = {value} is negative"
            ),
            Error::NonFiniteTally { field, index } => {
                write!(f, "{field}[{index}] is not finite")
            }
            Error::ZeroVarianceDim { dim } => {
                write!(
                    f,
                    "kde dimension {dim} has zero variance; use a fixed bandwidth"
                )
            }
            Error::BadAdaptiveKnn { n, k, scale } => {
                write!(
                    f,
                    "adaptive kNN bandwidth out of range: k = {k}, n = {n}, \
                     scale = {scale}; requires 1 <= k <= n - 1 and a finite \
                     positive scale"
                )
            }
            Error::ZeroAdaptiveBandwidth { sample } => {
                write!(
                    f,
                    "kde sample {sample} has a zero k-th-neighbour distance \
                     (duplicate points); a zero bandwidth is a delta spike"
                )
            }
            Error::BadDraw { value } => {
                write!(f, "kde draw uniform {value} is outside [0, 1)")
            }
            Error::NegativeWindow { index, value } => {
                write!(
                    f,
                    "weight-window lower bound at index {index} = {value} is negative"
                )
            }
            Error::NonFiniteWindow { index } => {
                write!(
                    f,
                    "weight-window lower bound at index {index} is not finite"
                )
            }
            Error::BadEnergyBounds { index } => {
                write!(
                    f,
                    "energy upper bound at index {index} is not positive and increasing"
                )
            }
            Error::BadMeshBounds { axis, index } => {
                write!(
                    f,
                    "mesh bound at axis {axis} index {index} is not finite and increasing"
                )
            }
            Error::BadEmissionOption {
                option,
                value,
                detail,
            } => write!(f, "emission option {option} = {value}: {detail}"),
            Error::Wwinp(m) => write!(f, "wwinp writer error: {m}"),
        }
    }
}

impl std::error::Error for Error {}

#[cfg(test)]
mod tests {
    use super::Error;

    #[test]
    fn error_display_covers_all_variants() {
        let cases: Vec<(Error, &str)> = vec![
            (Error::EmptyTally, "no volume elements"),
            (
                Error::LengthMismatch {
                    expected: 3,
                    got: 5,
                },
                "expected 3, got 5",
            ),
            (Error::ZeroMaxFlux { energy_group: 2 }, "energy group 2"),
            (Error::EmptyPdf, "at least one value"),
            (
                Error::NegativePdf {
                    index: 1,
                    value: -0.5,
                },
                "pdf[1]",
            ),
            (Error::NonFinitePdf { index: 0 }, "pdf[0]"),
            (Error::ZeroSumPdf, "sums to zero"),
            (
                Error::NegativeTally {
                    index: 4,
                    value: -1.0,
                },
                "index 4",
            ),
            (
                Error::NonFiniteTally {
                    field: "flux",
                    index: 7,
                },
                "flux[7]",
            ),
            (Error::ZeroVarianceDim { dim: 1 }, "dimension 1"),
            (
                Error::BadAdaptiveKnn {
                    n: 4,
                    k: 0,
                    scale: 1.0,
                },
                "k = 0",
            ),
            (Error::ZeroAdaptiveBandwidth { sample: 2 }, "sample 2"),
            (Error::BadDraw { value: 1.5 }, "outside [0, 1)"),
            (
                Error::NegativeWindow {
                    index: 2,
                    value: -0.1,
                },
                "index 2",
            ),
            (Error::NonFiniteWindow { index: 3 }, "index 3"),
            (Error::BadEnergyBounds { index: 1 }, "index 1"),
            (Error::BadMeshBounds { axis: 2, index: 4 }, "axis 2"),
            (
                Error::BadEmissionOption {
                    option: "survival_ratio",
                    value: "1".to_string(),
                    detail: "must be greater than 1",
                },
                "survival_ratio",
            ),
            (Error::Wwinp("disk".into()), "wwinp writer error"),
        ];
        for (err, needle) in cases {
            assert!(
                format!("{err}").contains(needle),
                "display of {err:?} should contain {needle:?}"
            );
            // Ensure the std::error::Error impl is linked.
            let _: &dyn std::error::Error = &err;
        }
    }
}
