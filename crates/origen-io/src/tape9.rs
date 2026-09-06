//! ORIGEN TAPE9 decay-constant reader.
//!
//! Parses the simplified TAPE9 decay table covered by this crate: one entry
//! per line,
//!
//! ```text
//! <nuclide> <lambda-s^-1>
//! ```
//!
//! Blank lines and `#` comments are skipped. Nuclide tokens are kept
//! verbatim (no name validation); decay constants must be finite and
//! non-negative (`0` encodes a stable nuclide).

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// One TAPE9 decay entry.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape9Entry {
    /// Nuclide token as written in the decay line.
    pub nuclide: String,
    /// Decay constant in s^-1 (`>= 0`).
    pub decay_const: f64,
}

impl Tape9Entry {
    /// Parse TAPE9 text into decay entries in file order.
    pub fn parse(text: &str) -> Result<Vec<Self>> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty TAPE9".to_string(),
            });
        }
        let mut entries = Vec::new();
        for (idx, raw) in text.lines().enumerate() {
            let line_no = idx + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 2 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `<nuclide> <lambda>`, got `{line}`"),
                });
            }
            let decay_const = parts[1].parse::<f64>().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("invalid decay constant `{}`", parts[1]),
            })?;
            if !decay_const.is_finite() || decay_const < 0.0 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("invalid decay constant `{}`", parts[1]),
                });
            }
            entries.push(Self {
                nuclide: parts[0].to_string(),
                decay_const,
            });
        }
        Ok(entries)
    }

    /// Read a TAPE9 file into decay entries in file order.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Vec<Self>> {
        Self::parse(&std::fs::read_to_string(path)?)
    }

    /// Find an entry by nuclide name (canonical id match, falling back to
    /// case-insensitive comparison for non-standard tokens).
    pub fn find<'a>(entries: &'a [Self], nuclide: &str) -> Option<&'a Self> {
        entries
            .iter()
            .find(|entry| names_match(&entry.nuclide, nuclide))
    }
}

/// True when two nuclide tokens name the same nuclide.
fn names_match(stored: &str, query: &str) -> bool {
    match (
        nuclei::NuclideId::from_name(stored),
        nuclei::NuclideId::from_name(query),
    ) {
        (Ok(a), Ok(b)) => a == b,
        _ => stored.eq_ignore_ascii_case(query),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        std::fs::read_to_string(format!(
            "{}/../../fixtures/origen/{name}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap()
    }

    fn parse_error(text: &str) -> (usize, String) {
        match Tape9Entry::parse(text).unwrap_err() {
            Error::Parse { line, msg } => (line, msg),
            Error::Io(_) => panic!("expected parse error, got io error"),
        }
    }

    #[test]
    fn rejects_empty() {
        assert!(Tape9Entry::parse("  \n").is_err());
    }

    #[test]
    fn parses_sample_fixture() {
        let entries = Tape9Entry::parse(&fixture("tape9_sample")).unwrap();
        assert_eq!(entries.len(), 5);
        assert_eq!(
            entries[0],
            Tape9Entry {
                nuclide: "U235".to_string(),
                decay_const: 3.1209e-17,
            }
        );
        assert_eq!(entries[4].nuclide, "Co60");
        assert_eq!(entries[4].decay_const, 4.1674e-09);
    }

    #[test]
    fn reads_sample_from_file() {
        let path = format!(
            "{}/../../fixtures/origen/tape9_sample",
            env!("CARGO_MANIFEST_DIR")
        );
        let entries = Tape9Entry::from_file(path).unwrap();
        assert_eq!(entries.len(), 5);
    }

    #[test]
    fn find_matches_by_name_and_misses() {
        let entries = Tape9Entry::parse(&fixture("tape9_sample")).unwrap();
        let hit = Tape9Entry::find(&entries, "Cs137").unwrap();
        assert_eq!(hit.decay_const, 7.3217e-10);
        // Canonical matching tolerates case/dash variants.
        assert!(Tape9Entry::find(&entries, "cs-137").is_some());
        assert!(Tape9Entry::find(&entries, "U233").is_none());
    }

    #[test]
    fn skips_comments_and_blank_lines() {
        let entries =
            Tape9Entry::parse("# header\n\nU235 3.1209e-17\n   # note\n\nU238 4.9161e-18\n")
                .unwrap();
        assert_eq!(entries.len(), 2);
    }

    #[test]
    fn reports_bad_lambda_with_line_number() {
        let (line, msg) = parse_error("U235 3.1209e-17\nU238 fast\n");
        assert_eq!(line, 2);
        assert!(msg.contains("decay constant"), "unexpected message: {msg}");
    }

    #[test]
    fn reports_wrong_field_count_with_line_number() {
        let (line, _) = parse_error("U235 3.1209e-17 extra\n");
        assert_eq!(line, 1);
    }

    #[test]
    fn rejects_negative_lambda_with_line_number() {
        let (line, _) = parse_error("U235 -1.0e-9\n");
        assert_eq!(line, 1);
    }
}
