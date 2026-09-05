//! ALARA photon-source (`.photonSrc`) files parsed into [`PhotonSource`].
//!
//! Each row has the form `<nuclide> <time> <strength...>` where `<time>` is
//! `shutdown` (0 s), `pre-irrad` (-1 s), or a `<value> <unit>` pair (`5000 s`,
//! `1 h`, `0.5 y`, ...), followed by one photon strength per energy group in
//! ALARA group order. `TOTAL` aggregate rows are kept verbatim as regular
//! groups; filter by [`PhotonGroup::nuclide`] when per-nuclide data is needed.

use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::error::{Error, Result};
use crate::output::time_unit_to_seconds;

/// Photon strengths of one nuclide at one cooling time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhotonGroup {
    /// Nuclide name in ALARA dialect, kept verbatim (`mn-56`, `TOTAL`).
    pub nuclide: String,
    /// Cooling time in seconds (`-1` = pre-irradiation, `0` = shutdown).
    pub time_s: f64,
    /// Photon strengths in ALARA group order.
    pub strengths: Vec<f64>,
}

impl PhotonGroup {
    /// Sum over groups; `0.0` when empty.
    pub fn total_strength(&self) -> f64 {
        self.strengths.iter().sum()
    }
}

/// Photon-source spectra parsed from a `.photonSrc` file.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PhotonSource {
    /// Groups in file order.
    pub groups: Vec<PhotonGroup>,
}

impl PhotonSource {
    /// Parse `<nuclide> <time> <strength...>` rows from in-memory text.
    ///
    /// Blank lines and `#` comments are skipped. Empty input is an error.
    // Inherent (not `std::str::FromStr`) so the crate-`Result` error type and
    // the `from_file` symmetry survive without trait imports at call sites.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty photon source".to_string(),
            });
        }
        let mut groups = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line_no = index + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            groups.push(parse_group_row(line, line_no)?);
        }
        Ok(Self { groups })
    }

    /// Read and parse a photon source from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::from_str(&text)
    }

    /// Sum over every group and strength (including `TOTAL` rows).
    pub fn total_strength(&self) -> f64 {
        self.groups.iter().map(PhotonGroup::total_strength).sum()
    }

    /// Number of groups.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// True when no groups were parsed.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }
}

/// Parse one `<nuclide> <time> <strength...>` row.
fn parse_group_row(line: &str, line_no: usize) -> Result<PhotonGroup> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let [nuclide, first, ..] = tokens.as_slice() else {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("expected `<nuclide> <time> <strengths...>`, found `{line}`"),
        });
    };
    let (time_s, strengths) = if first.eq_ignore_ascii_case("shutdown") {
        (0.0, &tokens[2..])
    } else if first.eq_ignore_ascii_case("pre-irrad") {
        (-1.0, &tokens[2..])
    } else {
        let value: f64 = first.parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected cooling time, found `{first}`"),
        })?;
        let unit = tokens.get(2).ok_or_else(|| Error::Parse {
            line: line_no,
            msg: format!("expected time unit after `{first}` in `{line}`"),
        })?;
        let factor = time_unit_to_seconds(unit).ok_or_else(|| Error::Parse {
            line: line_no,
            msg: format!("unknown cooling-time unit `{unit}`"),
        })?;
        (value * factor, &tokens[3..])
    };
    if strengths.is_empty() {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("no group strengths in `{line}`"),
        });
    }
    let mut parsed = Vec::with_capacity(strengths.len());
    for token in strengths {
        let strength: f64 = token.parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected group strength, found `{token}`"),
        })?;
        parsed.push(strength);
    }
    Ok(PhotonGroup {
        nuclide: (*nuclide).to_string(),
        time_s,
        strengths: parsed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rows in the style of `sample9.photonSrc` (2-group file).
    const SAMPLE9_STYLE: &str = "mn-56 \tshutdown\t0\t6.35458e+18\t2.83827e+18\n\
                                 mn-56 \t   5000 s   \t0\t4.37493e+18\t1.95406e+18\n\
                                 fe-56 \tshutdown\t0\t0\t0\n\
                                 TOTAL\tshutdown\t0\t6.35458e+18\t2.83827e+18\n";

    #[test]
    fn parse_collects_sample9_style_rows() {
        let source = PhotonSource::from_str(SAMPLE9_STYLE).unwrap();
        assert_eq!(source.len(), 4);
        assert!(!source.is_empty());

        assert_eq!(source.groups[0].nuclide, "mn-56");
        assert_eq!(source.groups[0].time_s, 0.0);
        assert_eq!(source.groups[0].strengths.len(), 3);
        assert_eq!(source.groups[0].strengths[0], 0.0);

        assert_eq!(source.groups[1].time_s, 5_000.0);
        assert_eq!(source.groups[1].strengths[1], 4.37493e18);

        assert_eq!(source.groups[3].nuclide, "TOTAL");

        let expected = 6.35458e18 + 2.83827e18 + 4.37493e18 + 1.95406e18 + 6.35458e18 + 2.83827e18;
        assert!((source.total_strength() - expected).abs() <= 1e4);
    }

    #[test]
    fn parse_handles_named_cooling_times() {
        let source = PhotonSource::from_str(
            "h-1 \t      1 h   \t0\t0\n\
             h-1 \t    0.5 y   \t1.5\n\
             h-1 \tpre-irrad\t2.5\n",
        )
        .unwrap();
        assert_eq!(source.groups[0].time_s, 3_600.0);
        assert_eq!(source.groups[1].time_s, 0.5 * 31_536_000.0);
        assert_eq!(source.groups[2].time_s, -1.0);
    }

    #[test]
    fn parse_rejects_bad_rows_with_line_numbers() {
        assert!(matches!(
            PhotonSource::from_str("mn-56 shutdown 1.0\nnope\n"),
            Err(Error::Parse { line: 2, .. })
        ));
        assert!(matches!(
            PhotonSource::from_str("mn-56 soon s 1.0\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            PhotonSource::from_str("mn-56 1 fortnight 1.0\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            PhotonSource::from_str("mn-56 shutdown\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            PhotonSource::from_str("  \n"),
            Err(Error::Parse { .. })
        ));
    }
}
