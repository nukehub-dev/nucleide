//! Errors for ALARA activation-code interop.
//!
//! Every variant carries the file/line context available at the failure site
//! so diagnostics can point back at the offending deck or output line.

use thiserror::Error;

/// Result alias for the `alara-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while reading or cross-linking ALARA inputs and outputs.
#[derive(Debug, Error)]
pub enum Error {
    /// Filesystem failure while reading an ALARA input or output file.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// A deck or output line failed to parse; `line` is 1-based.
    #[error("parse error at line {line}: {msg}")]
    Parse {
        /// 1-based line number of the offending line.
        line: usize,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An unrecognized top-level block keyword; `line` is 1-based.
    #[error("unknown block `{block}` at line {line}")]
    UnknownBlock {
        /// The rejected block keyword.
        block: String,
        /// 1-based line number where the keyword appeared.
        line: usize,
    },
    /// A dangling cross-reference between blocks (for example a schedule
    /// item naming an undefined flux, or a mixture naming an undefined
    /// library).
    #[error("dangling cross-reference: {0}")]
    CrossRef(String),
    /// An unrecognized or inconsistent physical unit.
    #[error("bad units: {0}")]
    BadUnits(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_carries_line_context() {
        let error = Error::Parse {
            line: 7,
            msg: "expected float".to_string(),
        };
        assert_eq!(error.to_string(), "parse error at line 7: expected float");

        let error = Error::UnknownBlock {
            block: "bogus".to_string(),
            line: 3,
        };
        assert_eq!(error.to_string(), "unknown block `bogus` at line 3");

        let error = Error::CrossRef("schedule `p1` needs flux `f1`".to_string());
        assert!(error.to_string().contains("f1"));

        let error = Error::BadUnits("fortnight".to_string());
        assert!(error.to_string().contains("fortnight"));
    }

    #[test]
    fn io_converts_from_std() {
        let io = std::io::Error::new(std::io::ErrorKind::NotFound, "gone");
        let error = Error::from(io);
        assert!(matches!(error, Error::Io(_)));
        assert!(error.to_string().contains("gone"));
    }
}
