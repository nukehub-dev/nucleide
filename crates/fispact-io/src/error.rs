//! Errors for the FISPACT-II output parser.

use thiserror::Error;

/// Result alias for the `fispact-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while parsing FISPACT-II `.fis` output.
#[derive(Debug, Error)]
pub enum Error {
    /// Filesystem failure while reading an output file.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// An output line failed to parse; `line` is 1-based.
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
            line: 9,
            msg: "expected inventory table".to_string(),
        };
        assert_eq!(
            error.to_string(),
            "parse error at line 9: expected inventory table"
        );
    }
}
