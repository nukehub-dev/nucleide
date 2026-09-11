//! ISOTXS multigroup cross-section library reader.
//!
//! Documented text-analog subset of the CCCC ISOTXS standard:
//!
//! ```text
//! # comment lines (starting with `#`) and blank lines are ignored
//! ISOTXS <ngroups>
//! NUCLIDE <label> <zaid> <ngroups>
//! <total xs group 1> ... <total xs group ngroups>   (may span several lines)
//! ... one NUCLIDE block per nuclide ...
//! ```
//!
//! Validation rules:
//!
//! - The `ISOTXS` header carries the library group count (`> 0`); every
//!   `NUCLIDE` record must repeat exactly that count, otherwise parsing fails
//!   with a [`crate::error::Error::Parse`] error carrying the 1-based physical
//!   line number of the offending record.
//! - Each nuclide contributes exactly `ngroups` total-cross-section values
//!   ([`IsotxsNuclide::total_xs`]), free-format across one or more lines. Short
//!   records (end of file or a new record before enough values), ragged lines
//!   (more values on one line than still needed), and non-numeric tokens all
//!   fail with the 1-based line number where detected.
//! - `zaid` is validated loosely: it must be non-empty.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// One nuclide entry in an ISOTXS library.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct IsotxsNuclide {
    /// Nuclide label as given on the `NUCLIDE` record.
    pub label: String,
    /// Nuclide ZAID (loosely validated: must be non-empty).
    pub zaid: String,
    /// Number of energy groups (always equals the library header count).
    pub groups: usize,
    /// Total cross section per group, in file order.
    pub total_xs: Vec<f64>,
}

/// ISOTXS multigroup cross-section library.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct IsotxsLib {
    /// Nuclide entries in file order.
    pub nuclides: Vec<IsotxsNuclide>,
}

impl IsotxsLib {
    /// Parse ISOTXS text in the documented subset.
    pub fn parse(text: &str) -> Result<Self> {
        let physical: Vec<&str> = text.lines().collect();
        // 1-based physical line count used for short-at-EOF reports; never 0 so
        // the reported line is always valid.
        let total = physical.len().max(1);
        let mut content: Vec<(usize, &str)> = Vec::new();
        for (index, raw) in physical.iter().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            content.push((index + 1, line));
        }
        let Some(&(header_line, header)) = content.first() else {
            return Err(Error::Parse {
                line: 1,
                msg: "empty ISOTXS library".to_string(),
            });
        };
        let head: Vec<&str> = header.split_whitespace().collect();
        if head.len() != 2 || head[0] != "ISOTXS" {
            return Err(Error::Parse {
                line: header_line,
                msg: format!("expected `ISOTXS <ngroups>` header, got `{header}`"),
            });
        }
        let groups: usize = head[1].parse().map_err(|_| Error::Parse {
            line: header_line,
            msg: format!("invalid ISOTXS group count `{}`", head[1]),
        })?;
        if groups == 0 {
            return Err(Error::Parse {
                line: header_line,
                msg: "ISOTXS group count must be positive".to_string(),
            });
        }

        let mut nuclides = Vec::new();
        let mut pos = 1;
        while pos < content.len() {
            let (lineno, line) = content[pos];
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.len() != 4 || tokens[0] != "NUCLIDE" {
                return Err(Error::Parse {
                    line: lineno,
                    msg: format!(
                        "expected `NUCLIDE <label> <zaid> <ngroups>` record, got `{line}`"
                    ),
                });
            }
            let (label, zaid) = (tokens[1], tokens[2]);
            if zaid.is_empty() {
                return Err(Error::Parse {
                    line: lineno,
                    msg: format!("NUCLIDE `{label}` has an empty ZAID"),
                });
            }
            let count: usize = tokens[3].parse().map_err(|_| Error::Parse {
                line: lineno,
                msg: format!("invalid NUCLIDE group count `{}`", tokens[3]),
            })?;
            if count != groups {
                return Err(Error::Parse {
                    line: lineno,
                    msg: format!(
                        "NUCLIDE `{label}` group count {count} does not match ISOTXS header {groups}"
                    ),
                });
            }
            pos += 1;
            let mut xs = Vec::with_capacity(groups);
            while xs.len() < groups {
                let Some(&(xline, xtext)) = content.get(pos) else {
                    return Err(Error::Parse {
                        line: total,
                        msg: format!(
                            "short cross-section record for `{label}`: expected {groups} values, got {}",
                            xs.len()
                        ),
                    });
                };
                let head_token = xtext.split_whitespace().next().unwrap_or_default();
                if head_token == "NUCLIDE" || head_token == "ISOTXS" {
                    return Err(Error::Parse {
                        line: xline,
                        msg: format!(
                            "short cross-section record for `{label}`: expected {groups} values, got {}",
                            xs.len()
                        ),
                    });
                }
                let values: Vec<&str> = xtext.split_whitespace().collect();
                if values.len() > groups - xs.len() {
                    return Err(Error::Parse {
                        line: xline,
                        msg: format!(
                            "ragged cross-section record for `{label}`: expected {} more values, got {} on this line",
                            groups - xs.len(),
                            values.len()
                        ),
                    });
                }
                for token in values {
                    let value: f64 = token.parse().map_err(|_| Error::Parse {
                        line: xline,
                        msg: format!("invalid cross-section value `{token}` for `{label}`"),
                    })?;
                    xs.push(value);
                }
                pos += 1;
            }
            nuclides.push(IsotxsNuclide {
                label: label.to_string(),
                zaid: zaid.to_string(),
                groups,
                total_xs: xs,
            });
        }
        Ok(Self { nuclides })
    }

    /// Read and parse an ISOTXS file from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(&text)
    }

    /// Find a nuclide by label.
    pub fn find(&self, label: &str) -> Option<&IsotxsNuclide> {
        self.nuclides.iter().find(|n| n.label == label)
    }

    /// Number of nuclides in the library.
    pub fn len(&self) -> usize {
        self.nuclides.len()
    }

    /// Whether the library holds no nuclides.
    pub fn is_empty(&self) -> bool {
        self.nuclides.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!("{}/../../fixtures/cccc/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn parses_fixture_two_by_three() {
        let lib = IsotxsLib::from_file(fixture("isotxs_sample")).unwrap();
        assert_eq!(lib.len(), 2);
        assert!(!lib.is_empty());
        let uranium = lib.find("U235").expect("U235 present");
        assert_eq!(uranium.zaid, "92235");
        assert_eq!(uranium.groups, 3);
        assert_eq!(uranium.total_xs, vec![1.1, 2.2, 3.3]);
        // Second nuclide spans multiple XS lines.
        let plutonium = lib.find("PU239").expect("PU239 present");
        assert_eq!(plutonium.zaid, "94239");
        assert_eq!(plutonium.groups, 3);
        assert_eq!(plutonium.total_xs, vec![4.4, 5.5, 6.6]);
        assert!(lib.find("U238").is_none());
        assert!(IsotxsLib::default().is_empty());
    }

    #[test]
    fn rejects_empty() {
        let err = IsotxsLib::parse("  \n").unwrap_err();
        assert!(matches!(err, Error::Parse { line: 1, .. }));
    }

    #[test]
    fn group_count_mismatch_reports_line() {
        let text = "ISOTXS 3\nNUCLIDE U235 92235 2\n1.0 2.0\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 2);
                assert!(msg.contains("2") && msg.contains('3'));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn short_xs_at_eof_reports_line() {
        let text = "ISOTXS 3\nNUCLIDE U235 92235 3\n1.0 2.0\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 3);
                assert!(msg.contains("short"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn short_xs_before_next_nuclide_reports_line() {
        let text = "ISOTXS 3\nNUCLIDE U235 92235 3\n1.0 2.0\nNUCLIDE PU239 94239 3\n4.4 5.5 6.6\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, .. } => assert_eq!(line, 4),
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn ragged_xs_reports_line() {
        let text = "ISOTXS 2\nNUCLIDE U235 92235 2\n1.0 2.0 3.0\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 3);
                assert!(msg.contains("ragged"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn bad_float_reports_line() {
        let text = "ISOTXS 2\nNUCLIDE U235 92235 2\n1.0 banana\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 3);
                assert!(msg.contains("banana"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn header_rejections_report_line() {
        // Wrong token count.
        let err = IsotxsLib::parse("ISOTXS\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("ISOTXS <ngroups>"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Wrong keyword.
        let err = IsotxsLib::parse("BOGUS 2\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("ISOTXS <ngroups>"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Non-numeric group count.
        let err = IsotxsLib::parse("ISOTXS two\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("invalid ISOTXS group count"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Zero group count.
        let err = IsotxsLib::parse("ISOTXS 0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("positive"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn nuclide_group_count_token_rejected() {
        let text = "ISOTXS 2\nNUCLIDE u 92235 x\n1.0 2.0\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 2);
                assert!(msg.contains("invalid NUCLIDE group count"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn isotxs_record_inside_xs_values_is_short() {
        // A second `ISOTXS` header where XS values belong reads as a short record.
        let text = "ISOTXS 2\nNUCLIDE u 92235 2\nISOTXS 2\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 3);
                assert!(msg.contains("short cross-section record"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn from_file_io_error_surfaces() {
        let missing = format!(
            "{}/../../fixtures/cccc/definitely_missing",
            env!("CARGO_MANIFEST_DIR")
        );
        let err = IsotxsLib::from_file(missing).unwrap_err();
        assert!(matches!(err, Error::Io(_)), "got {err:?}");
    }

    #[test]
    fn unexpected_record_reports_line() {
        let text = "ISOTXS 2\nBOGUS U235 92235 2\n";
        let err = IsotxsLib::parse(text).unwrap_err();
        match err {
            Error::Parse { line, .. } => assert_eq!(line, 2),
            err => panic!("expected Parse error, got {err:?}"),
        }
    }
}
