//! ORIGEN TAPE5 input-echo reader.
//!
//! Parses the simplified TAPE5 input echo covered by this crate:
//!
//! ```text
//! <title cards: leading free-text lines>
//! FLUX= <n/cm2/s> DAYS= <days>   (one per irradiation step)
//! MAT <name>                      (starts a material block)
//! <nuclide> <grams>               (material entries)
//! ```
//!
//! Blank lines and `#` comments are skipped. Title cards are the leading
//! free-text lines before the first `FLUX=`/`DAYS=` or `MAT` card. Nuclide
//! tokens in material entries are kept verbatim (no name validation). This
//! is an input summary only; no burnup driving is performed here.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// One TAPE5 irradiation step.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape5Step {
    /// Neutron flux in n/cm²/s (`>= 0`; `0` encodes a decay-only step).
    pub flux: f64,
    /// Step length in days (`> 0`).
    pub days: f64,
}

/// One TAPE5 material block: a name plus `(nuclide, grams)` entries.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape5Material {
    /// Material name from the `MAT` card.
    pub name: String,
    /// `(nuclide token, grams)` pairs in file order.
    pub grams: Vec<(String, f64)>,
}

/// TAPE5 input summary.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tape5 {
    /// Title cards in file order (leading free-text lines).
    pub titles: Vec<String>,
    /// `FLUX=`/`DAYS=` irradiation steps in file order.
    pub irradiation_steps: Vec<Tape5Step>,
    /// `MAT` material blocks in file order.
    pub materials: Vec<Tape5Material>,
}

impl Tape5 {
    /// Parse TAPE5 text.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty TAPE5".to_string(),
            });
        }
        let mut out = Self::default();
        let mut seen_directive = false;
        let mut current: Option<usize> = None;
        for (idx, raw) in text.lines().enumerate() {
            let line_no = idx + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let upper = line.to_ascii_uppercase();
            if upper.contains("FLUX=") || upper.contains("DAYS=") {
                out.irradiation_steps.push(parse_step(line, line_no)?);
                seen_directive = true;
                continue;
            }
            if is_mat_card(line) {
                let name = line[3..].trim().to_string();
                if name.is_empty() {
                    return Err(Error::Parse {
                        line: line_no,
                        msg: "MAT card missing material name".to_string(),
                    });
                }
                out.materials.push(Tape5Material {
                    name,
                    grams: Vec::new(),
                });
                current = Some(out.materials.len() - 1);
                seen_directive = true;
                continue;
            }
            if !seen_directive {
                out.titles.push(line.to_string());
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                let grams = parse_non_negative(parts[1], line_no, "grams")?;
                match current {
                    Some(i) => out.materials[i].grams.push((parts[0].to_string(), grams)),
                    None => {
                        return Err(Error::Parse {
                            line: line_no,
                            msg: "material entry outside MAT block".to_string(),
                        });
                    }
                }
                continue;
            }
            return Err(Error::Parse {
                line: line_no,
                msg: format!("unrecognized TAPE5 card `{line}`"),
            });
        }
        Ok(out)
    }

    /// Read a TAPE5 file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::parse(&std::fs::read_to_string(path)?)
    }
}

/// True when the first whitespace-separated token is `MAT` (any case).
fn is_mat_card(line: &str) -> bool {
    line.split_whitespace()
        .next()
        .is_some_and(|token| token.eq_ignore_ascii_case("MAT"))
}

/// Parse one `FLUX= <v> DAYS= <v>` step line.
fn parse_step(line: &str, line_no: usize) -> Result<Tape5Step> {
    let flux_token = card_value(line, "FLUX=").ok_or_else(|| Error::Parse {
        line: line_no,
        msg: "FLUX step missing FLUX= value".to_string(),
    })?;
    let days_token = card_value(line, "DAYS=").ok_or_else(|| Error::Parse {
        line: line_no,
        msg: "FLUX step missing DAYS= value".to_string(),
    })?;
    let flux = flux_token.parse::<f64>().map_err(|_| Error::Parse {
        line: line_no,
        msg: format!("invalid FLUX value `{flux_token}`"),
    })?;
    if !flux.is_finite() || flux < 0.0 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("invalid FLUX value `{flux_token}`"),
        });
    }
    let days = days_token.parse::<f64>().map_err(|_| Error::Parse {
        line: line_no,
        msg: format!("invalid DAYS value `{days_token}`"),
    })?;
    if !days.is_finite() || days <= 0.0 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("invalid DAYS value `{days_token}`"),
        });
    }
    Ok(Tape5Step { flux, days })
}

/// First whitespace-delimited token after a case-insensitive `KEY=` marker.
fn card_value<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let pos = line.to_ascii_uppercase().find(key)?;
    line[pos + key.len()..].split_whitespace().next()
}

/// Parse a finite, non-negative value such as a material `grams` entry.
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
        match Tape5::parse(text).unwrap_err() {
            Error::Parse { line, msg } => (line, msg),
            Error::Io(_) => panic!("expected parse error, got io error"),
        }
    }

    #[test]
    fn rejects_empty() {
        assert!(Tape5::parse("  \n").is_err());
    }

    #[test]
    fn parses_sample_fixture() {
        let tape = Tape5::parse(&fixture("tape5_sample")).unwrap();
        assert_eq!(
            tape.titles,
            vec![
                "SYNTHETIC PWR PIN - TAPE5 SAMPLE".to_string(),
                "CASE 1 - BASE DEPLETION".to_string(),
            ]
        );
        assert_eq!(
            tape.irradiation_steps,
            vec![
                Tape5Step {
                    flux: 3.0e13,
                    days: 100.0,
                },
                Tape5Step {
                    flux: 0.0,
                    days: 30.0,
                },
            ]
        );
        assert_eq!(tape.materials.len(), 2);
        assert_eq!(tape.materials[0].name, "fuel");
        assert_eq!(
            tape.materials[0].grams,
            vec![
                ("U235".to_string(), 10.5),
                ("U238".to_string(), 1000.0),
                ("Pu239".to_string(), 0.25),
            ]
        );
        assert_eq!(tape.materials[1].name, "clad");
        assert_eq!(tape.materials[1].grams, vec![("Zr90".to_string(), 50.0)]);
    }

    #[test]
    fn reads_sample_from_file() {
        let path = format!(
            "{}/../../fixtures/origen/tape5_sample",
            env!("CARGO_MANIFEST_DIR")
        );
        let tape = Tape5::from_file(path).unwrap();
        assert_eq!(tape.titles.len(), 2);
        assert_eq!(tape.irradiation_steps.len(), 2);
        assert_eq!(tape.materials.len(), 2);
    }

    #[test]
    fn reports_bad_flux_value_with_line_number() {
        let (line, msg) = parse_error("TITLE\nFLUX= hot DAYS= 10.0\n");
        assert_eq!(line, 2);
        assert!(msg.contains("FLUX"), "unexpected message: {msg}");
    }

    #[test]
    fn reports_step_missing_days_with_line_number() {
        let (line, msg) = parse_error("TITLE\nFLUX= 3.0e13\n");
        assert_eq!(line, 2);
        assert!(msg.contains("DAYS"), "unexpected message: {msg}");
    }

    #[test]
    fn reports_entry_outside_mat_block_with_line_number() {
        let (line, msg) = parse_error("TITLE\nFLUX= 3.0e13 DAYS= 10.0\nU235 1.0\n");
        assert_eq!(line, 3);
        assert!(
            msg.contains("outside MAT block"),
            "unexpected message: {msg}"
        );
    }

    #[test]
    fn reports_unknown_card_with_line_number() {
        let (line, msg) = parse_error("TITLE\nMAT fuel\nU235 1.0\nBOGUS CARD HERE\n");
        assert_eq!(line, 4);
        assert!(
            msg.contains("unrecognized TAPE5 card"),
            "unexpected message: {msg}"
        );
    }

    #[test]
    fn reports_nameless_mat_with_line_number() {
        let (line, msg) = parse_error("TITLE\nMAT\n");
        assert_eq!(line, 2);
        assert!(msg.contains("MAT"), "unexpected message: {msg}");
    }

    #[test]
    fn rejects_negative_grams_with_line_number() {
        let (line, _) = parse_error("TITLE\nMAT fuel\nU235 -1.0\n");
        assert_eq!(line, 3);
    }
}
