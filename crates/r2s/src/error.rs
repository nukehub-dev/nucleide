//! Errors for the R2S shutdown-dose-rate orchestration.

use thiserror::Error;

/// Result alias for the `r2s` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while building or running an R2S workflow.
#[derive(Debug, Error)]
pub enum Error {
    /// Filesystem failure while reading workflow inputs.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// A workflow input failed validation.
    #[error("invalid input: {0}")]
    Invalid(String),
    /// A dangling cross-reference between workflow pieces.
    #[error("dangling cross-reference: {0}")]
    CrossRef(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_carries_context() {
        assert!(Error::Invalid("empty mesh".to_string())
            .to_string()
            .contains("empty mesh"));
        assert!(Error::CrossRef("zone Z pins ghost flux".to_string())
            .to_string()
            .contains("ghost flux"));
    }
}
