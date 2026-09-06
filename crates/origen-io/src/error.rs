//! Errors for the scoped ORIGEN 2.2 TAPE readers.

use thiserror::Error;

/// Result alias for the `origen-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while reading ORIGEN TAPE files.
#[derive(Debug, Error)]
pub enum Error {
    /// Filesystem failure while reading a TAPE file.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// A TAPE line failed to parse; `line` is 1-based.
    #[error("parse error at line {line}: {msg}")]
    Parse {
        /// 1-based line number of the offending line.
        line: usize,
        /// What was expected versus what was found.
        msg: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_carries_line_context() {
        let error = Error::Parse {
            line: 3,
            msg: "expected TAPE card".to_string(),
        };
        assert_eq!(
            error.to_string(),
            "parse error at line 3: expected TAPE card"
        );
    }
}
