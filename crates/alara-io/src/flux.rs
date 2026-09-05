//! ALARA default-format group-flux files.
//!
//! The `flux` deck block points at a flux file plus a scalar normalization
//! and an integer skip value. The file itself holds one free-format list of
//! group fluxes per fine-mesh interval; blank lines may separate intervals
//! but comments are **not** permitted inside flux files. [`FluxSpectra`]
//! parses that layout; the older flat [`FluxSpec`] is kept for compatibility.

use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::error::{Error, Result};

/// One ALARA multigroup flux spectrum.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct FluxSpec {
    /// Spectrum name, as referenced by schedule items.
    pub name: String,
    /// Group fluxes in ALARA group order.
    pub groups: Vec<f64>,
}

impl FluxSpec {
    /// Parse whitespace-separated group fluxes from in-memory text.
    ///
    /// Blank lines and `#` comments are skipped; any other token must parse
    /// as `f64`. Empty input is an error.
    pub fn parse(name: &str, text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: format!("empty flux spectrum `{name}`"),
            });
        }
        let mut groups = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            for token in line.split_whitespace() {
                let value: f64 = token.parse().map_err(|_| Error::Parse {
                    line: index + 1,
                    msg: format!("expected group flux, found `{token}`"),
                })?;
                groups.push(value);
            }
        }
        Ok(Self {
            name: name.to_string(),
            groups,
        })
    }

    /// Read and parse a flux spectrum from disk.
    pub fn from_file(name: &str, path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(name, &text)
    }

    /// Sum over groups; `0.0` when empty.
    pub fn total(&self) -> f64 {
        self.groups.iter().sum()
    }

    /// Number of groups.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// True when no groups are stored.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }
}

/// One ALARA default-format flux file: N-group lists, one per fine-mesh interval.
///
/// Blank lines delimit intervals; every other token must parse as `f64`
/// (comments are not permitted in flux files, so a `#` token is a parse
/// error). All intervals must hold the same number of groups.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct FluxSpectra {
    /// Spectrum name, as referenced by schedule items.
    pub name: String,
    /// Groups per interval (uniform across [`Self::intervals`]).
    pub groups_per_interval: usize,
    /// Group fluxes per fine-mesh interval, in file order.
    pub intervals: Vec<Vec<f64>>,
}

impl FluxSpectra {
    /// Parse default-format flux text: blank-line-separated intervals of
    /// whitespace-separated floats, in a single pass with no regex.
    ///
    /// Empty input and ragged intervals (unequal group counts) are errors;
    /// every error carries a 1-based line number.
    pub fn parse(name: &str, text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: format!("empty flux spectra `{name}`"),
            });
        }
        // Heuristic pre-allocation: formatted values average ~12 B each.
        let mut intervals: Vec<Vec<f64>> = Vec::new();
        let mut end_lines: Vec<usize> = Vec::new();
        let mut current: Vec<f64> = Vec::with_capacity((text.len() / 12).clamp(16, 4096));
        let mut last_data_line: usize = 0;

        for (index, raw) in text.lines().enumerate() {
            let line_no = index + 1;
            let line = raw.trim();
            if line.is_empty() {
                if !current.is_empty() {
                    let cap = current.len();
                    intervals.push(std::mem::replace(&mut current, Vec::with_capacity(cap)));
                    end_lines.push(last_data_line);
                }
                continue;
            }
            last_data_line = line_no;
            for token in line.split_whitespace() {
                let value: f64 = token.parse().map_err(|_| Error::Parse {
                    line: line_no,
                    msg: format!("expected group flux, found `{token}`"),
                })?;
                current.push(value);
            }
        }
        if !current.is_empty() {
            end_lines.push(last_data_line);
            intervals.push(current);
        }
        if intervals.is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: format!("empty flux spectra `{name}`"),
            });
        }
        let groups = intervals[0].len();
        for (position, interval) in intervals.iter().enumerate().skip(1) {
            if interval.len() != groups {
                return Err(Error::Parse {
                    line: end_lines[position],
                    msg: format!(
                        "ragged flux intervals: interval {} has {} groups, expected {groups}",
                        position + 1,
                        interval.len()
                    ),
                });
            }
        }
        Ok(Self {
            name: name.to_string(),
            groups_per_interval: groups,
            intervals,
        })
    }

    /// Read and parse a default-format flux file from disk.
    pub fn from_file(name: &str, path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(name, &text)
    }

    /// Consistency check: uniform group count across intervals.
    ///
    /// Returns the shared count, or a [`Error::Parse`] error naming the
    /// first ragged interval (1-based interval index in the message).
    pub fn groups_per_interval(&self) -> Result<usize> {
        let first = self.intervals.first().map_or(0, Vec::len);
        for (position, interval) in self.intervals.iter().enumerate() {
            if interval.len() != first {
                return Err(Error::Parse {
                    line: 1,
                    msg: format!(
                        "ragged flux intervals: interval {} has {} groups, expected {first}",
                        position + 1,
                        interval.len()
                    ),
                });
            }
        }
        Ok(first)
    }

    /// Number of fine-mesh intervals stored.
    pub fn num_intervals(&self) -> usize {
        self.intervals.len()
    }

    /// True when no intervals are stored.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Interval selected by an ALARA `flux`-block skip value: the
    /// `skip`-th N-group entry in the file. `None` when out of range.
    pub fn interval_after_skip(&self, skip: usize) -> Option<&[f64]> {
        self.intervals.get(skip).map(Vec::as_slice)
    }

    /// Copy with every group flux multiplied by `scale` (the ALARA scalar
    /// normalization); name and group structure are preserved.
    pub fn scaled(&self, scale: f64) -> Self {
        Self {
            name: self.name.clone(),
            groups_per_interval: self.groups_per_interval,
            intervals: self
                .intervals
                .iter()
                .map(|interval| interval.iter().map(|v| v * scale).collect())
                .collect(),
        }
    }

    /// Sum over every group of every interval; `0.0` when empty.
    pub fn total(&self) -> f64 {
        self.intervals.iter().flatten().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        let path = format!(
            "{}/../../fixtures/alara/flux/{name}",
            env!("CARGO_MANIFEST_DIR")
        );
        std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read fixture {path}: {error}"))
    }

    #[test]
    fn parse_collects_group_fluxes() {
        let flux = FluxSpec::parse("core", "# comment\n1.0 2.0\n0.5\n").unwrap();
        assert_eq!(flux.name, "core");
        assert_eq!(flux.groups, [1.0, 2.0, 0.5]);
        assert_eq!(flux.len(), 3);
        assert!(!flux.is_empty());
        assert_eq!(flux.total(), 3.5);
    }

    #[test]
    fn parse_rejects_bad_tokens_with_line_numbers() {
        let error = FluxSpec::parse("core", "1.0\nnope\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }));
    }

    #[test]
    fn parse_rejects_empty_input() {
        assert!(matches!(
            FluxSpec::parse("core", "  \n"),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn spectra_splits_intervals_on_blank_lines() {
        let spectra = FluxSpectra::parse("core", "1.0 2.0\n3.0\n\n4.0 5.0\n6.0\n").unwrap();
        assert_eq!(spectra.name, "core");
        assert_eq!(spectra.groups_per_interval, 3);
        assert_eq!(spectra.num_intervals(), 2);
        assert_eq!(spectra.intervals[0], [1.0, 2.0, 3.0]);
        assert_eq!(spectra.intervals[1], [4.0, 5.0, 6.0]);
        assert_eq!(spectra.groups_per_interval().unwrap(), 3);
        assert_eq!(spectra.total(), 21.0);
    }

    #[test]
    fn spectra_tolerates_leading_and_repeat_blanks() {
        let spectra = FluxSpectra::parse("core", "\n\n1.0\n\n\n2.0\n\n").unwrap();
        assert_eq!(spectra.num_intervals(), 2);
        assert_eq!(spectra.groups_per_interval, 1);
    }

    #[test]
    fn spectra_rejects_comments() {
        // ALARA flux files permit no comments: `#` must fail as a float.
        let error = FluxSpectra::parse("core", "1.0\n# comment\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }));
    }

    #[test]
    fn spectra_rejects_bad_floats_with_line_numbers() {
        let error = FluxSpectra::parse("core", "1.0 2.0\nnope 3.0\n").unwrap_err();
        match error {
            Error::Parse { line, msg } => {
                assert_eq!(line, 2);
                assert!(msg.contains("nope"));
            }
            other => panic!("expected Parse, found {other:?}"),
        }
    }

    #[test]
    fn spectra_rejects_ragged_intervals() {
        let error = FluxSpectra::parse("core", "1.0 2.0\n\n3.0\n").unwrap_err();
        assert!(matches!(error, Error::Parse { .. }));
    }

    #[test]
    fn spectra_rejects_empty_input() {
        assert!(matches!(
            FluxSpectra::parse("core", "  \n \n"),
            Err(Error::Parse { line: 1, .. })
        ));
    }

    #[test]
    fn spectra_skip_selects_interval() {
        let spectra = FluxSpectra::parse("core", "1.0\n\n2.0\n\n3.0\n").unwrap();
        assert_eq!(spectra.interval_after_skip(0), Some([1.0].as_slice()));
        assert_eq!(spectra.interval_after_skip(2), Some([3.0].as_slice()));
        assert_eq!(spectra.interval_after_skip(3), None);
    }

    #[test]
    fn spectra_scaled_multiplies_all_groups() {
        let spectra = FluxSpectra::parse("core", "1.0 2.0\n\n3.0 4.0\n").unwrap();
        let scaled = spectra.scaled(2.0);
        assert_eq!(scaled.name, "core");
        assert_eq!(scaled.groups_per_interval, 2);
        assert_eq!(scaled.intervals[0], [2.0, 4.0]);
        assert_eq!(scaled.total(), 2.0 * spectra.total());
    }

    #[test]
    fn fixture_fluxin2_holds_three_intervals_of_175_groups() {
        let spectra = FluxSpectra::parse("fluxin2", &fixture("fluxin2")).unwrap();
        assert_eq!(spectra.num_intervals(), 3);
        assert_eq!(spectra.groups_per_interval, 175);
        assert_eq!(spectra.groups_per_interval().unwrap(), 175);
        assert_eq!(spectra.intervals.iter().map(Vec::len).sum::<usize>(), 525);
        assert!(spectra.intervals.iter().all(|iv| iv.len() == 175));
        assert!(spectra.total() > 0.0);
        assert_eq!(spectra.interval_after_skip(1).unwrap().len(), 175);
        let twice = spectra.scaled(2.0);
        assert_eq!(twice.total(), 2.0 * spectra.total());
    }

    #[test]
    fn fixture_fluxin_zeros_holds_175_zeros() {
        let spectra = FluxSpectra::parse("zeros", &fixture("fluxin_zeros")).unwrap();
        assert_eq!(spectra.num_intervals(), 1);
        assert_eq!(spectra.groups_per_interval, 175);
        assert!(spectra.intervals[0].iter().all(|v| *v == 0.0));
        assert_eq!(spectra.total(), 0.0);
    }
}
