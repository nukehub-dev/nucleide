//! Error type for the tritium crate.

use thiserror::Error;

/// Result alias for the `tritium` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while validating transport data, trap specs, boundary
/// conditions, grids, or solver options — and when the per-step
/// mobile/trap/face coupling fails to converge.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    /// A transport datum is non-finite or out of range: `{0}`.
    #[error("tritium: invalid transport data: {0}")]
    BadData(&'static str),
    /// A trap spec entry is invalid: `{0}`.
    #[error("tritium: invalid trap spec: {0}")]
    BadTrap(&'static str),
    /// A boundary-condition parameter is invalid: `{0}`.
    #[error("tritium: invalid boundary condition: {0}")]
    BadBoundary(&'static str),
    /// An internal interface condition is not supported (Sieverts and Henry
    /// are): `{0}`.
    #[error("tritium: unsupported internal interface: {0}")]
    UnsupportedInterface(&'static str),
    /// A grid entry is invalid: `{0}`.
    #[error("tritium: invalid grid: {0}")]
    BadGrid(&'static str),
    /// A per-cell query index falls outside the stack's total cell count.
    #[error("tritium: cell index {index} is out of range for a stack of {total} cells")]
    BadCellIndex {
        /// The offending cell index.
        index: usize,
        /// Total cell count across the stack.
        total: usize,
    },
    /// A solver option is invalid: `{0}`.
    #[error("tritium: invalid solver option: {0}")]
    BadOption(&'static str),
    /// An initial-state entry is invalid: `{0}`.
    #[error("tritium: invalid initial state: {0}")]
    BadState(&'static str),
    /// The per-step mobile/trap/face Picard coupling (including the
    /// recombination face Newton) did not converge within its iteration cap.
    #[error("tritium: step coupling did not converge within its iteration cap")]
    NotConverged,
    /// The implicit stepper exceeded `{0}` accepted steps.
    #[error("tritium: step budget exhausted ({0} steps)")]
    StepBudget(usize),
    /// The tridiagonal elimination failed: `{0}`.
    #[error("tritium: tridiagonal solve failed: {0}")]
    Tridiag(String),
}
