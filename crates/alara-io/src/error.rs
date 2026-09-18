//! Errors for ALARA activation-code interop.
//!
//! Every variant carries the file/line context available at the failure site
//! so diagnostics can point back at the offending deck or output line.

use thiserror::Error;

/// Result alias for the `alara-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while reading or cross-linking ALARA inputs and outputs.
#[derive(Debug, Error)]
#[non_exhaustive]
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
    /// A clearance computation met an inventory nuclide absent from the
    /// caller-supplied (or vendored) clearance limit table.
    #[error("nuclide {nuclide} has no clearance limit (table holds {table_len} entries)")]
    MissingClearanceLimit {
        /// Name of the unmatched nuclide.
        nuclide: String,
        /// Number of entries in the clearance table consulted.
        table_len: usize,
    },
    /// A clearance input value was rejected: negative or non-finite activity,
    /// a non-positive/non-finite limit (EU table insert, S6 A2, or S7 IAEA
    /// level), or a non-positive/non-finite S7 total mass.
    #[error("bad clearance value for {nuclide}: {msg}")]
    BadClearanceValue {
        /// Name of the offending nuclide.
        nuclide: String,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S1 total-activity input value was rejected: negative or
    /// non-finite activity, or a non-finite overflowed class sum.
    #[error("bad activity value for {nuclide}: {msg}")]
    BadActivityValue {
        /// Name of the offending nuclide (`total` for a sum overflow).
        nuclide: String,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S2 decay-heat input value was rejected: negative or non-finite
    /// activity or decay energy, or a non-finite overflowed class sum.
    #[error("bad heat value for {nuclide}: {msg}")]
    BadHeatValue {
        /// Name of the offending nuclide (`total` for a sum overflow).
        nuclide: String,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S4/S5 hazard input value was rejected: negative or non-finite
    /// activity or dose coefficient, or a non-finite overflowed dose sum.
    #[error("bad hazard value for {nuclide}: {msg}")]
    BadHazardValue {
        /// Name of the offending nuclide (`total` for a sum overflow).
        nuclide: String,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S1 entry carried a Table VI IRT with no pinned α/β/γ mapping
    /// (unused 8, 9; unknown 10; unlisted 0, 5–7, 18, 21–27, 28+).
    #[error("nuclide {nuclide} has unmapped IRT {irt} (no pinned alpha/beta/gamma class)")]
    UnmappedIrt {
        /// Name of the offending nuclide.
        nuclide: String,
        /// The rejected decay-type identifier.
        irt: u8,
    },
    /// An S1 split fraction was missing, superfluous, or out of range:
    /// split IRTs (12, 13, 15) require a finite `alpha_frac` in `[0, 1]`;
    /// wholly-mapped IRTs take none.
    #[error("bad IRT split for {nuclide} (IRT {irt}): {msg}")]
    BadIrtSplit {
        /// Name of the offending nuclide.
        nuclide: String,
        /// The decay-type identifier the fraction was (or was not) for.
        irt: u8,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S3 gamma dose-rate scalar was rejected: negative or non-finite
    /// specific activity or source mass, a non-finite distance, a negative
    /// or non-finite group intensity, or a non-finite overflowed dose sum.
    #[error("bad dose value for {field}: {msg}")]
    BadDoseValue {
        /// Name of the offending field (`total` for a sum overflow).
        field: String,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S3 attenuation input was rejected: negative or non-finite air or
    /// mixture coefficient, a zero mixture coefficient under the slab
    /// kernel (which divides by it), or a malformed
    /// fractions/element-table mixture fold.
    #[error("bad attenuation value at group {group}: {msg}")]
    BadAttenuation {
        /// 0-based gamma group index (`usize::MAX` for whole-input mixture
        /// shape errors with no single group to blame).
        group: usize,
        /// What was expected versus what was found.
        msg: String,
    },
    /// An S6 transport-ratio input value was rejected: negative or
    /// non-finite activity, a non-positive/non-finite A2 limit (zero A2 is
    /// a loud division error, never `inf`), or a non-finite overflowed sum.
    #[error("bad transport value for {nuclide}: {msg}")]
    BadTransportValue {
        /// Name of the offending nuclide (`total` for a sum overflow).
        nuclide: String,
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
    fn clearance_variants_display_key_context() {
        let error = Error::MissingClearanceLimit {
            nuclide: "Mn54".to_string(),
            table_len: 3,
        };
        let text = error.to_string();
        assert!(text.contains("Mn54"), "msg was `{text}`");
        assert!(text.contains("3"), "msg was `{text}`");

        let error = Error::BadClearanceValue {
            nuclide: "Co60".to_string(),
            msg: "activity must be finite and >= 0, got -1".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("Co60"), "msg was `{text}`");
        assert!(text.contains("-1"), "msg was `{text}`");
    }

    #[test]
    fn io_converts_from_std() {
        let io = std::io::Error::new(std::io::ErrorKind::NotFound, "gone");
        let error = Error::from(io);
        assert!(matches!(error, Error::Io(_)));
        assert!(error.to_string().contains("gone"));
    }

    #[test]
    fn sublet_variants_display_key_context() {
        let error = Error::BadActivityValue {
            nuclide: "Co60".to_string(),
            msg: "activity must be finite and >= 0, got -1".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("Co60"), "msg was `{text}`");
        assert!(text.contains("-1"), "msg was `{text}`");

        let error = Error::BadHeatValue {
            nuclide: "H3".to_string(),
            msg: "e_beta_ev must be finite and >= 0, got NaN".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("H3"), "msg was `{text}`");

        let error = Error::UnmappedIrt {
            nuclide: "Co60".to_string(),
            irt: 10,
        };
        let text = error.to_string();
        assert!(text.contains("Co60"), "msg was `{text}`");
        assert!(text.contains("10"), "msg was `{text}`");

        let error = Error::BadIrtSplit {
            nuclide: "U235".to_string(),
            irt: 12,
            msg: "split IRT requires alpha_frac".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("U235"), "msg was `{text}`");
        assert!(text.contains("12"), "msg was `{text}`");

        let error = Error::BadTransportValue {
            nuclide: "Co60".to_string(),
            msg: "A2 limit must be finite and > 0, got 0".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("Co60"), "msg was `{text}`");
        assert!(text.contains("0"), "msg was `{text}`");

        let error = Error::BadHazardValue {
            nuclide: "Co60".to_string(),
            msg: "ingestion: coefficient must be finite and >= 0, got -1".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("Co60"), "msg was `{text}`");
        assert!(text.contains("-1"), "msg was `{text}`");
    }

    #[test]
    fn dose_variants_display_key_context() {
        let error = Error::BadDoseValue {
            field: "activity_bq_per_kg".to_string(),
            msg: "activity_bq_per_kg must be finite and >= 0, got -1".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains("activity_bq_per_kg"), "msg was `{text}`");
        assert!(text.contains("-1"), "msg was `{text}`");

        let error = Error::BadAttenuation {
            group: 2,
            msg: "mu must be finite and > 0 for the slab ratio, got 0".to_string(),
        };
        let text = error.to_string();
        assert!(text.contains('2'), "msg was `{text}`");
        assert!(text.contains("slab ratio"), "msg was `{text}`");
    }
}
