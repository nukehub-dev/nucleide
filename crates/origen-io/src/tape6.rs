//! ORIGEN TAPE6 output-inventory reader.
//!
//! Parses the simplified TAPE6 inventory covered by this crate: one record
//! per line,
//!
//! ```text
//! <nuclide> <grams> <activity-Bq>
//! ```
//!
//! Blank lines and `#` comments are skipped. Nuclide tokens are validated
//! with [`nuclei::NuclideId`]; unknown names are 1-based parse errors.
//! Nuclide tokens are stored verbatim as written.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// One TAPE6 inventory record.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape6Record {
    /// Nuclide token as written in the inventory line.
    pub nuclide: String,
    /// Inventory mass in grams (`>= 0`).
    pub grams: f64,
    /// Activity in becquerel (`>= 0`).
    pub activity_bq: f64,
}

/// TAPE6 output inventory.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape6 {
    /// Inventory records in file order.
    pub records: Vec<Tape6Record>,
}

impl Tape6 {
    /// Parse TAPE6 text.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty TAPE6".to_string(),
            });
        }
        let mut records = Vec::new();
        for (idx, raw) in text.lines().enumerate() {
            let line_no = idx + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `<nuclide> <grams> <Bq>`, got `{line}`"),
                });
            }
            if nuclei::NuclideId::from_name(parts[0]).is_err() {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("unknown nuclide `{}`", parts[0]),
                });
            }
            let grams = parse_non_negative(parts[1], line_no, "grams")?;
            let activity_bq = parse_non_negative(parts[2], line_no, "activity")?;
            records.push(Tape6Record {
                nuclide: parts[0].to_string(),
                grams,
                activity_bq,
            });
        }
        Ok(Self { records })
    }

    /// Read a TAPE6 file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::parse(&std::fs::read_to_string(path)?)
    }

    /// Find a record by nuclide name (canonical id match, falling back to
    /// case-insensitive comparison for non-standard tokens).
    pub fn find(&self, nuclide: &str) -> Option<&Tape6Record> {
        self.records
            .iter()
            .find(|record| names_match(&record.nuclide, nuclide))
    }

    /// Total inventory activity in becquerel.
    pub fn total_activity(&self) -> f64 {
        self.records.iter().map(|record| record.activity_bq).sum()
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

/// Parse a finite, non-negative inventory value.
fn parse_non_negative(token: &str, line_no: usize, what: &str) -> Result<f64> {
    let value = token.parse::<f64>().map_err(|_| Error::Parse {
        line: line_no,
        msg: format!("invalid {what} value `{token}`"),
    })?;
    if !value.is_finite() || value < 0.0 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("invalid {what} value `{token}`"),
        });
    }
    Ok(value)
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
        match Tape6::parse(text).unwrap_err() {
            Error::Parse { line, msg } => (line, msg),
            Error::Io(_) => panic!("expected parse error, got io error"),
        }
    }

    #[test]
    fn rejects_empty() {
        assert!(Tape6::parse("  \n").is_err());
    }

    #[test]
    fn parses_sample_fixture() {
        let tape = Tape6::parse(&fixture("tape6_sample")).unwrap();
        assert_eq!(tape.records.len(), 4);
        assert_eq!(
            tape.records[0],
            Tape6Record {
                nuclide: "U235".to_string(),
                grams: 10.2,
                activity_bq: 8.16e5,
            }
        );
        assert_eq!(tape.records[3].nuclide, "Cs137");
        assert_eq!(tape.records[3].grams, 0.012);
        assert_eq!(tape.records[3].activity_bq, 3.84e10);
    }

    #[test]
    fn reads_sample_from_file() {
        let path = format!(
            "{}/../../fixtures/origen/tape6_sample",
            env!("CARGO_MANIFEST_DIR")
        );
        let tape = Tape6::from_file(path).unwrap();
        assert_eq!(tape.records.len(), 4);
    }

    #[test]
    fn total_activity_sums_all_records() {
        let tape = Tape6::parse(&fixture("tape6_sample")).unwrap();
        let expected = 8.16e5 + 1.23e4 + 5.63e8 + 3.84e10;
        assert!((tape.total_activity() - expected).abs() < 1.0);
    }

    #[test]
    fn find_matches_by_name_and_misses() {
        let tape = Tape6::parse(&fixture("tape6_sample")).unwrap();
        let hit = tape.find("Pu239").unwrap();
        assert_eq!(hit.grams, 0.245);
        assert_eq!(hit.activity_bq, 5.63e8);
        // Canonical matching tolerates case/dash variants.
        assert!(tape.find("pu-239").is_some());
        assert!(tape.find("U233").is_none());
    }

    #[test]
    fn reports_unknown_nuclide_with_line_number() {
        let (line, msg) = parse_error("U235 10.0 8.16e5\nXx999 1.0 2.0\n");
        assert_eq!(line, 2);
        assert!(msg.contains("Xx999"), "unexpected message: {msg}");
    }

    #[test]
    fn reports_bad_grams_with_line_number() {
        let (line, msg) = parse_error("U235 lots 8.16e5\n");
        assert_eq!(line, 1);
        assert!(msg.contains("grams"), "unexpected message: {msg}");
    }

    #[test]
    fn reports_wrong_field_count_with_line_number() {
        let (line, _) = parse_error("U235 10.0\n");
        assert_eq!(line, 1);
    }
}
