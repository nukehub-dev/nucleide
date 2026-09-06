//! Errors for CCCC binary-standard readers and the PARTISN deck writer.

use thiserror::Error;

/// Result alias for the `cccc-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while reading CCCC files or writing PARTISN decks.
#[derive(Debug, Error)]
pub enum Error {
    /// Filesystem failure while reading an input file.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// A record failed to parse; `line` or record index is 1-based where known.
    #[error("parse error at record {line}: {msg}")]
    Parse {
        /// 1-based record or line number.
        line: usize,
        /// What was expected versus what was found.
        msg: String,
    },
    /// A dangling cross-reference (for example a PARTISN zone naming an
    /// unknown ISOTXS nuclide).
    #[error("dangling cross-reference: {0}")]
    CrossRef(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_carries_context() {
        let error = Error::Parse {
            line: 4,
            msg: "short record".to_string(),
        };
        assert_eq!(error.to_string(), "parse error at record 4: short record");
    }
}
