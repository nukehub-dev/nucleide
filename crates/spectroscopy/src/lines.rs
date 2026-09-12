//! Runtime decay-line tables (SDEF option b): one documented interchange TSV.
//!
//! [`parse_lines_tsv`] reads caller-supplied decay lines from text — one
//! `energy_MeV intensity` pair per line, `#` comments and blank lines
//! skipped — into the `(energy, intensity)` pairs that
//! [`crate::sdef::normalize_decay_lines`] (E9) consumes. No new equation is
//! involved: this module is syntax only (the E9 caller-constants posture,
//! fed from a file instead of a list). Nothing is vendored and no decay
//! data ships with the crate; fixtures stay hand-built synthetic rows.
//!
//! No legacy reader ships alongside: PyNE's gamma/alpha/beta line families
//! query its compiled HDF5 store (no text line-list format exists locally
//! to mirror), a full ENSDF-text grammar is out of proportion, ICRP-107
//! exports are license-blocked under AD-7, and ENDF MF8/MT457 tapes carry
//! mean energies rather than lines. Callers holding such exports convert
//! them to this TSV once, outside the crate.
//!
//! Interchange format (version 1, documented here — the only TSV this
//! crate reads):
//!
//! ```text
//! # energy_MeV intensity   (comment lines start with `#`)
//! 0.662 2.0
//! 1.170 1.0
//! ```
//!
//! - Whitespace-separated, exactly two columns per data line; energies in
//!   MeV, intensities in arbitrary (non-negative) units — E9 normalizes.
//! - Malformed rows are [`Error::MalformedLinesRow`] with the 1-based line
//!   number; text with no data rows is [`Error::EmptyLines`].

use crate::Error;

/// Parse interchange-TSV decay lines into `(energy_MeV, intensity)` pairs.
///
/// Row order is preserved (E9 sorts after merging); use
/// [`crate::sdef::normalize_decay_lines`] for the `(energy, probability)`
/// bins. See the module docs for the format.
pub fn parse_lines_tsv(text: &str) -> Result<Vec<(f64, f64)>, Error> {
    let mut lines = Vec::new();
    for (index, raw) in text.lines().enumerate() {
        let line_no = index + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        let [energy, intensity] = cols.as_slice() else {
            return Err(Error::MalformedLinesRow {
                line: line_no,
                text: raw.to_string(),
            });
        };
        let energy: f64 = energy.parse().map_err(|_| Error::MalformedLinesRow {
            line: line_no,
            text: raw.to_string(),
        })?;
        let intensity: f64 = intensity.parse().map_err(|_| Error::MalformedLinesRow {
            line: line_no,
            text: raw.to_string(),
        })?;
        lines.push((energy, intensity));
    }
    if lines.is_empty() {
        return Err(Error::EmptyLines);
    }
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "# synthetic Cs-137-like pair (hand values, MeV)\n\
                          0.662 2.0\n\
                          \n\
                          1.170 1.0\n";

    #[test]
    fn parses_pairs_skipping_comments_and_blanks() {
        assert_eq!(
            parse_lines_tsv(SAMPLE).unwrap(),
            vec![(0.662, 2.0), (1.17, 1.0)]
        );
    }

    #[test]
    fn empty_and_comment_only_is_empty_lines() {
        assert_eq!(parse_lines_tsv(""), Err(Error::EmptyLines));
        assert_eq!(
            parse_lines_tsv("# nothing here\n\n"),
            Err(Error::EmptyLines)
        );
    }

    #[test]
    fn malformed_rows_carry_line_numbers() {
        for (text, line) in [
            ("0.662\n", 1),
            ("0.662 1.0 3.0\n", 1),
            ("0.662 lots\n", 1),
            ("# ok\n0.662 1.0\nbogus row here\n", 3),
        ] {
            assert_eq!(
                parse_lines_tsv(text),
                Err(Error::MalformedLinesRow {
                    line,
                    text: text.lines().nth(line - 1).unwrap().to_string(),
                }),
                "input {text:?}"
            );
        }
    }

    #[test]
    fn parsed_rows_feed_e9_normalization() {
        let bins = crate::sdef::normalize_decay_lines(&parse_lines_tsv(SAMPLE).unwrap()).unwrap();
        assert_eq!(bins, vec![(0.662, 2.0 / 3.0), (1.17, 1.0 / 3.0)]);
    }
}
