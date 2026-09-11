//! RTFLUX / ATFLUX / RZFLUX flux-file readers.
//!
//! Documented text-analog subset of the CCCC flux-file standards:
//!
//! ```text
//! # comment lines (starting with `#`) and blank lines are ignored
//! RTFLUX|ATFLUX|RZFLUX <npoints> <ngroups>
//! <free-format floats: exactly npoints * ngroups values>
//! ```
//!
//! Validation rules:
//!
//! - The header keyword must match the [`FluxKind`] passed to
//!   [`FluxFile::parse`]; a mismatch fails with a
//!   [`crate::error::Error::Parse`] error at the 1-based header line.
//! - `<npoints>` and `<ngroups>` must both be positive.
//! - The body must hold exactly `npoints * ngroups` floats. A short body fails
//!   at the last physical line; a long body fails at the 1-based line where
//!   the excess value appears. Non-numeric tokens fail at their own line.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Which CCCC flux-file flavor a [`FluxFile`] holds.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum FluxKind {
    /// Regular (rectangular-geometry) flux file.
    #[default]
    Rtflux,
    /// ATFLUX variant.
    Atflux,
    /// RZFLUX (cylindrical-geometry) variant.
    Rzflux,
}

impl FluxKind {
    /// Header keyword for this flavor.
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Rtflux => "RTFLUX",
            Self::Atflux => "ATFLUX",
            Self::Rzflux => "RZFLUX",
        }
    }

    /// Parse a header keyword into a flavor.
    pub fn from_keyword(token: &str) -> Option<Self> {
        match token {
            "RTFLUX" => Some(Self::Rtflux),
            "ATFLUX" => Some(Self::Atflux),
            "RZFLUX" => Some(Self::Rzflux),
            _ => None,
        }
    }
}

/// Regular flux file: `ngroups` flux values per spatial point.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct FluxFile {
    /// Flux-file flavor (must match the header keyword).
    pub kind: FluxKind,
    /// Number of energy groups.
    pub groups: usize,
    /// Flux values in file order (`npoints * ngroups` entries).
    pub values: Vec<f64>,
    /// Values per spatial point (equals `groups` in this subset).
    pub per_point: usize,
}

impl FluxFile {
    /// Parse flux text of the given flavor.
    pub fn parse(kind: FluxKind, text: &str) -> Result<Self> {
        let physical: Vec<&str> = text.lines().collect();
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
                msg: "empty flux file".to_string(),
            });
        };
        let head: Vec<&str> = header.split_whitespace().collect();
        if head.len() != 3 {
            return Err(Error::Parse {
                line: header_line,
                msg: format!(
                    "expected `<RTFLUX|ATFLUX|RZFLUX> <npoints> <ngroups>` header, got `{header}`"
                ),
            });
        }
        let header_kind = FluxKind::from_keyword(head[0]).ok_or_else(|| Error::Parse {
            line: header_line,
            msg: format!(
                "expected `<RTFLUX|ATFLUX|RZFLUX> <npoints> <ngroups>` header, got `{header}`"
            ),
        })?;
        if header_kind != kind {
            return Err(Error::Parse {
                line: header_line,
                msg: format!(
                    "flux kind mismatch: expected {}, found {}",
                    kind.keyword(),
                    head[0]
                ),
            });
        }
        let npoints: usize = head[1].parse().map_err(|_| Error::Parse {
            line: header_line,
            msg: format!("invalid flux point count `{}`", head[1]),
        })?;
        let ngroups: usize = head[2].parse().map_err(|_| Error::Parse {
            line: header_line,
            msg: format!("invalid flux group count `{}`", head[2]),
        })?;
        if npoints == 0 || ngroups == 0 {
            return Err(Error::Parse {
                line: header_line,
                msg: "flux point and group counts must be positive".to_string(),
            });
        }
        let expected = npoints.checked_mul(ngroups).ok_or_else(|| Error::Parse {
            line: header_line,
            msg: "flux point/group count overflows".to_string(),
        })?;

        let mut values = Vec::new();
        for &(lineno, line) in &content[1..] {
            for token in line.split_whitespace() {
                if values.len() >= expected {
                    return Err(Error::Parse {
                        line: lineno,
                        msg: format!("too many flux values: expected {expected}, file holds more"),
                    });
                }
                let value: f64 = token.parse().map_err(|_| Error::Parse {
                    line: lineno,
                    msg: format!("invalid flux value `{token}`"),
                })?;
                values.push(value);
            }
        }
        if values.len() != expected {
            return Err(Error::Parse {
                line: total,
                msg: format!(
                    "short flux record: expected {expected} values, got {}",
                    values.len()
                ),
            });
        }
        Ok(Self {
            kind,
            groups: ngroups,
            values,
            per_point: ngroups,
        })
    }

    /// Read and parse a flux file of the given flavor from disk.
    pub fn from_file(kind: FluxKind, path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(kind, &text)
    }

    /// Number of spatial points in the file.
    pub fn npoints(&self) -> usize {
        self.values.len().checked_div(self.per_point).unwrap_or(0)
    }

    /// Flux vector for point `i`, or `None` when out of range.
    pub fn point(&self, i: usize) -> Option<&[f64]> {
        let start = i.checked_mul(self.per_point)?;
        let end = start.checked_add(self.per_point)?;
        self.values.get(start..end)
    }

    /// Sum of all flux values.
    pub fn total(&self) -> f64 {
        self.values.iter().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!("{}/../../fixtures/cccc/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn parses_fixture_two_points_three_groups() {
        let flux = FluxFile::from_file(FluxKind::Rtflux, fixture("rtflux_sample")).unwrap();
        assert_eq!(flux.kind, FluxKind::Rtflux);
        assert_eq!(flux.groups, 3);
        assert_eq!(flux.per_point, 3);
        assert_eq!(flux.npoints(), 2);
        assert_eq!(flux.values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(flux.point(0), Some([1.0, 2.0, 3.0].as_slice()));
        assert_eq!(flux.point(1), Some([4.0, 5.0, 6.0].as_slice()));
        assert_eq!(flux.point(2), None);
        assert_eq!(flux.total(), 21.0);
    }

    #[test]
    fn all_kinds_parse_their_keyword() {
        for (keyword, kind) in [
            ("RTFLUX", FluxKind::Rtflux),
            ("ATFLUX", FluxKind::Atflux),
            ("RZFLUX", FluxKind::Rzflux),
        ] {
            let text = format!("{keyword} 1 2\n1.0 2.0\n");
            let flux = FluxFile::parse(kind, &text).unwrap();
            assert_eq!(flux.kind, kind);
            assert_eq!(flux.npoints(), 1);
        }
        assert_eq!(FluxKind::Rtflux.keyword(), "RTFLUX");
        assert_eq!(FluxKind::from_keyword("NOPE"), None);
    }

    #[test]
    fn rejects_empty() {
        let err = FluxFile::parse(FluxKind::Rtflux, "  \n").unwrap_err();
        assert!(matches!(err, Error::Parse { line: 1, .. }));
    }

    #[test]
    fn kind_mismatch_reports_header_line() {
        let err = FluxFile::parse(FluxKind::Rtflux, "ATFLUX 1 2\n1.0 2.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("mismatch"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn short_body_reports_line() {
        let err = FluxFile::parse(FluxKind::Rtflux, "RTFLUX 2 3\n1.0 2.0 3.0\n4.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 3);
                assert!(msg.contains("short"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn ragged_body_reports_line() {
        let err = FluxFile::parse(FluxKind::Rtflux, "RTFLUX 1 2\n1.0 2.0 3.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 2);
                assert!(msg.contains("too many"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn keyword_round_trips_for_all_flavors() {
        for kind in [FluxKind::Rtflux, FluxKind::Atflux, FluxKind::Rzflux] {
            assert_eq!(FluxKind::from_keyword(kind.keyword()), Some(kind));
        }
    }

    #[test]
    fn header_shape_errors_report_line() {
        // Too few header tokens.
        let err = FluxFile::parse(FluxKind::Rtflux, "RTFLUX 1\n1.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("RTFLUX|ATFLUX|RZFLUX"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Unknown header keyword.
        let err = FluxFile::parse(FluxKind::Rtflux, "NOPE 1 2\n1.0 2.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("RTFLUX|ATFLUX|RZFLUX"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Non-numeric point count.
        let err = FluxFile::parse(FluxKind::Rtflux, "RTFLUX x 2\n1.0 2.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("point count"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Non-numeric group count.
        let err = FluxFile::parse(FluxKind::Rtflux, "RTFLUX 1 y\n1.0\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("group count"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }

        // Zero counts are rejected.
        for header in ["RTFLUX 0 2", "RTFLUX 1 0"] {
            let err =
                FluxFile::parse(FluxKind::Rtflux, &format!("{header}\n1.0 2.0\n")).unwrap_err();
            match err {
                Error::Parse { line, msg } => {
                    assert_eq!(line, 1);
                    assert!(msg.contains("positive"), "{msg}");
                }
                err => panic!("expected Parse error, got {err:?}"),
            }
        }

        // Point/group product overflow.
        let max = usize::MAX;
        let huge = format!("RTFLUX {max} {max}\n");
        let err = FluxFile::parse(FluxKind::Rtflux, &huge).unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 1);
                assert!(msg.contains("overflow"), "{msg}");
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }

    #[test]
    fn point_math_guards_overflow() {
        // per_point of zero yields no points instead of trapping on division.
        let empty = FluxFile {
            kind: FluxKind::Rtflux,
            groups: 0,
            values: Vec::new(),
            per_point: 0,
        };
        assert_eq!(empty.npoints(), 0);

        // A huge index overflows the offset and reports None.
        let flux = FluxFile {
            kind: FluxKind::Rtflux,
            groups: 2,
            values: vec![1.0, 2.0],
            per_point: 2,
        };
        assert_eq!(flux.point(usize::MAX), None);
    }

    #[test]
    fn from_file_io_error_surfaces() {
        let missing = format!(
            "{}/../../fixtures/cccc/definitely_missing",
            env!("CARGO_MANIFEST_DIR")
        );
        let err = FluxFile::from_file(FluxKind::Rtflux, missing).unwrap_err();
        assert!(matches!(err, Error::Io(_)), "got {err:?}");
    }

    #[test]
    fn bad_float_reports_line() {
        let err = FluxFile::parse(FluxKind::Rzflux, "RZFLUX 1 2\n1.0 banana\n").unwrap_err();
        match err {
            Error::Parse { line, msg } => {
                assert_eq!(line, 2);
                assert!(msg.contains("banana"));
            }
            err => panic!("expected Parse error, got {err:?}"),
        }
    }
}
