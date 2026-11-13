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

// Neumaier accumulation retains small contributions beside a large sum.
// Callers check the corrected total for non-finite arithmetic.
#[derive(Clone, Copy, Default)]
struct RegroupSum {
    sum: f64,
    correction: f64,
}

impl RegroupSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    fn total(&self) -> f64 {
        self.sum + self.correction
    }
}

/// Regrouped integrated flux and a report of the assumed spectral model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FluxRegrouping {
    /// Integrated flux in each target interval, in ascending energy order.
    pub flux: Vec<f64>,
    /// Differential-flux assumption; always `"flat_energy"`.
    pub model: String,
    /// Sum of the supplied source group fluxes.
    pub source_total: f64,
    /// Sum of the returned target group fluxes.
    pub target_total: f64,
    /// Zero-based source intervals crossed by a target interior edge.
    /// Includes intervals with zero flux, in source order without duplicates.
    pub split_source_groups: Vec<usize>,
}

/// Regroup integrated flux assuming constant differential flux per eV in
/// each source interval.
///
/// Both grids must contain at least two finite, nonnegative, strictly
/// ascending eV edges, and their endpoints must match exactly. `flux` must
/// contain one finite, nonnegative value per source interval. Each overlap
/// contributes `source_flux * (overlap_width / source_width)`; no clipping,
/// extrapolation, edge merging, or unit conversion is performed. Non-finite
/// widths, totals and output values are errors.
///
/// The overlap sweep uses O(source groups + target groups) time. Identity
/// grids preserve input values exactly. Other grids conserve integrated
/// flux up to floating-point roundoff, using compensated accumulation for
/// totals and contributions to each target group. Contributions smaller
/// than the minimum representable subnormal can round to zero; relative
/// conservation cannot be guaranteed for such underflowing sources.
/// Regrouping does not in general preserve
/// reaction rates or recover structure inside a source group. Photon energy
/// centroids cannot supply the required bin edges.
pub fn regroup_flux_flat_energy(
    source_bounds_ev: &[f64],
    flux: &[f64],
    target_bounds_ev: &[f64],
) -> Result<FluxRegrouping> {
    validate_regroup_grid("source", source_bounds_ev)?;
    validate_regroup_grid("target", target_bounds_ev)?;
    if flux.len() != source_bounds_ev.len() - 1 {
        return Err(Error::BadFluxRegrouping(
            "flux length must equal the number of source intervals".into(),
        ));
    }
    if source_bounds_ev[0] != target_bounds_ev[0]
        || source_bounds_ev.last() != target_bounds_ev.last()
    {
        return Err(Error::BadFluxRegrouping(
            "source and target endpoints must match exactly".into(),
        ));
    }
    let mut source_sum = RegroupSum::default();
    for (i, &value) in flux.iter().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(Error::BadFluxRegrouping(format!(
                "source flux at group {i} must be finite and nonnegative"
            )));
        }
        source_sum.add(value);
        if !source_sum.total().is_finite() {
            return Err(Error::BadFluxRegrouping("source total overflowed".into()));
        }
    }
    let source_total = source_sum.total();
    if source_bounds_ev == target_bounds_ev {
        return Ok(FluxRegrouping {
            flux: flux.to_vec(),
            model: "flat_energy".into(),
            source_total,
            target_total: source_total,
            split_source_groups: Vec::new(),
        });
    }

    let mut target_sums = vec![RegroupSum::default(); target_bounds_ev.len() - 1];
    let mut split_source_groups = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < flux.len() && j < target_sums.len() {
        let source_hi = source_bounds_ev[i + 1];
        let target_hi = target_bounds_ev[j + 1];
        let overlap_lo = source_bounds_ev[i].max(target_bounds_ev[j]);
        let overlap_hi = source_hi.min(target_hi);
        let width = source_hi - source_bounds_ev[i];
        let overlap = overlap_hi - overlap_lo;
        if !overlap.is_finite() || overlap < 0.0 {
            return Err(Error::BadFluxRegrouping("invalid overlap width".into()));
        }
        // Divide widths first: an intermediate flux density can overflow
        // on a narrow interval even when the integrated answer is finite.
        let fraction = overlap / width;
        let contribution = if fraction < f64::MIN_POSITIVE && overlap > 0.0 {
            // Subnormal ratios lose relative precision even before they
            // round to zero. Reorder the arithmetic using a normal density
            // when possible; otherwise the product-first order is safe
            // because overlap / width is smaller than MIN_POSITIVE.
            let density = flux[i] / width;
            if density.is_normal() {
                density * overlap
            } else {
                (flux[i] * overlap) / width
            }
        } else {
            flux[i] * fraction
        };
        target_sums[j].add(contribution);
        if !contribution.is_finite() || !target_sums[j].total().is_finite() {
            return Err(Error::BadFluxRegrouping(format!(
                "target flux at group {j} overflowed"
            )));
        }
        if target_hi < source_hi && split_source_groups.last() != Some(&i) {
            split_source_groups.push(i);
        }
        if source_hi <= target_hi {
            i += 1;
        }
        if target_hi <= source_hi {
            j += 1;
        }
    }
    let target_flux: Vec<f64> = target_sums.iter().map(RegroupSum::total).collect();
    let mut target_sum = RegroupSum::default();
    for &value in &target_flux {
        target_sum.add(value);
    }
    let target_total = target_sum.total();
    if !target_total.is_finite() {
        return Err(Error::BadFluxRegrouping("target total overflowed".into()));
    }
    Ok(FluxRegrouping {
        flux: target_flux,
        model: "flat_energy".into(),
        source_total,
        target_total,
        split_source_groups,
    })
}

fn validate_regroup_grid(name: &str, edges: &[f64]) -> Result<()> {
    if edges.len() < 2 {
        return Err(Error::BadFluxRegrouping(format!(
            "{name} grid must contain at least two edges"
        )));
    }
    for (i, &edge) in edges.iter().enumerate() {
        if !edge.is_finite() || edge < 0.0 {
            return Err(Error::BadFluxRegrouping(format!(
                "{name} edge {i} must be finite and nonnegative"
            )));
        }
    }
    for (i, pair) in edges.windows(2).enumerate() {
        let width = pair[1] - pair[0];
        if !width.is_finite() || width <= 0.0 {
            return Err(Error::BadFluxRegrouping(format!(
                "{name} interval {i} must have a finite positive width"
            )));
        }
    }
    Ok(())
}

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

    fn assert_regroup_close(actual: f64, expected: f64) {
        if expected == 0.0 {
            assert_eq!(actual, 0.0);
        } else {
            assert!((actual - expected).abs() <= 1e-12 * expected.abs());
        }
    }

    #[test]
    fn regroup_exact_overlap_vector_and_report() {
        let result =
            regroup_flux_flat_energy(&[0.0, 2.0, 5.0], &[4.0, 9.0], &[0.0, 1.0, 3.0, 5.0]).unwrap();
        assert_eq!(result.flux, [2.0, 5.0, 6.0]);
        assert_eq!(result.model, "flat_energy");
        assert_eq!(result.source_total, 13.0);
        assert_eq!(result.target_total, 13.0);
        assert_eq!(result.split_source_groups, [0, 1]);
    }

    #[test]
    fn regroup_identity_preserves_bits_and_has_no_splits() {
        let edges = [0.0, 0.25, 1.0, 10.0];
        let values = [-0.0, 1e-200, 1e200];
        let result = regroup_flux_flat_energy(&edges, &values, &edges).unwrap();
        for (actual, expected) in result.flux.iter().zip(values) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert!(result.split_source_groups.is_empty());
        assert_eq!(result.target_total, result.source_total);
    }

    #[test]
    fn regroup_aligned_coarsening_is_group_addition() {
        let values = [2.0, 5.0, 11.0, 23.0, 47.0];
        let result =
            regroup_flux_flat_energy(&[0.0, 1.0, 2.0, 4.0, 8.0, 16.0], &values, &[0.0, 4.0, 16.0])
                .unwrap();
        let expected = [values[..3].iter().sum::<f64>(), values[3..].iter().sum()];
        assert_eq!(result.flux, expected);
        assert_eq!(result.source_total, result.target_total);
        assert!(result.split_source_groups.is_empty());
    }

    #[test]
    fn regroup_compensates_many_small_groups_beside_a_large_group() {
        let small = f64::EPSILON / 4.0;
        let mut values = vec![small; 100_001];
        values[0] = 1.0;
        let edges: Vec<f64> = (0..=100_001).map(f64::from).collect();
        // Each small value is an exact power of two; multiplication by
        // 100000 gives an independent, exactly representable aggregate.
        let tail = 100_000.0 * small;
        let expected = 1.0 + tail;
        let split = regroup_flux_flat_energy(&edges, &values, &[0.0, 1.0, 100_001.0]).unwrap();
        assert_eq!(split.flux, [1.0, tail]);
        assert_eq!(split.source_total, expected);
        assert_eq!(split.target_total, expected);
        let merged = regroup_flux_flat_energy(&edges, &values, &[0.0, 100_001.0]).unwrap();
        assert_eq!(merged.flux, [expected]);
        assert_eq!(merged.source_total, expected);
        assert_eq!(merged.target_total, expected);
        let identity = regroup_flux_flat_energy(&edges, &values, &edges).unwrap();
        assert_eq!(identity.flux, values);
        assert_eq!(identity.source_total, expected);
        assert_eq!(identity.target_total, expected);
    }

    #[test]
    fn regroup_constant_differential_flux_matches_analytic_integrals() {
        let source = [0.0, 0.13, 0.7, 3.1, 19.0];
        let target = [0.0, 0.01, 0.2, 1.0, 7.0, 18.0, 19.0];
        for density in [0.0, 1e-200, 3.5, 1e200] {
            let values: Vec<f64> = source.windows(2).map(|p| density * (p[1] - p[0])).collect();
            let result = regroup_flux_flat_energy(&source, &values, &target).unwrap();
            for (actual, edges) in result.flux.iter().zip(target.windows(2)) {
                assert_regroup_close(*actual, density * (edges[1] - edges[0]));
            }
            assert_regroup_close(result.target_total, result.source_total);
        }
    }

    #[test]
    fn regroup_split_then_recombine_recovers_flat_model() {
        let source = [0.0, 2.0, 5.0];
        let fine = [0.0, 0.5, 1.5, 2.0, 3.0, 4.5, 5.0];
        let values = [7.0, 13.0];
        let split = regroup_flux_flat_energy(&source, &values, &fine).unwrap();
        let restored = regroup_flux_flat_energy(&fine, &split.flux, &source).unwrap();
        for (actual, expected) in restored.flux.iter().zip(values) {
            assert_regroup_close(*actual, expected);
        }
    }

    #[test]
    fn regroup_reports_zero_flux_splits_once_each() {
        let result = regroup_flux_flat_energy(
            &[0.0, 3.0, 4.0],
            &[0.0, 0.0],
            &[0.0, 1.0, 2.0, 3.0, 3.5, 4.0],
        )
        .unwrap();
        assert_eq!(result.flux, [0.0; 5]);
        assert_eq!(result.split_source_groups, [0, 1]);
        assert_eq!(result.source_total, 0.0);
        assert_eq!(result.target_total, 0.0);
    }

    #[test]
    fn regroup_conserves_flux_without_preserving_threshold_response() {
        let result = regroup_flux_flat_energy(&[0.0, 2.0], &[2.0], &[0.0, 1.0, 2.0]).unwrap();
        // A true fine spectrum [0, 2] has threshold response 2 above 1 eV;
        // the flat-energy reconstruction has response 1 despite equal totals.
        assert_eq!(result.target_total, 2.0);
        assert_eq!(result.flux[1], 1.0);
    }

    #[test]
    fn regroup_handles_narrow_intervals_without_flux_density_overflow() {
        let tiny = f64::MIN_POSITIVE;
        let result =
            regroup_flux_flat_energy(&[0.0, tiny], &[1e100], &[0.0, tiny / 2.0, tiny]).unwrap();
        assert_eq!(result.flux, [5e99, 5e99]);
        assert_eq!(result.target_total, result.source_total);
    }

    #[test]
    fn regroup_preserves_representable_flux_when_overlap_fraction_underflows() {
        let tiny = f64::MIN_POSITIVE;
        let result =
            regroup_flux_flat_energy(&[0.0, f64::MAX], &[f64::MAX], &[0.0, tiny, f64::MAX])
                .unwrap();
        assert_eq!(result.flux[0], tiny);
        assert_eq!(result.flux[1], f64::MAX);
        assert_eq!(result.target_total, result.source_total);
    }

    #[test]
    fn regroup_preserves_accuracy_for_nonzero_subnormal_overlap_fraction() {
        let result =
            regroup_flux_flat_energy(&[0.0, f64::MAX], &[f64::MAX], &[0.0, 5e-16, f64::MAX])
                .unwrap();
        // Unit differential flux integrates to the interval width.
        assert_eq!(result.flux[0], 5e-16);
        assert_eq!(result.flux[1], f64::MAX);
    }

    #[test]
    fn regroup_subnormal_fraction_handles_non_normal_density() {
        let tiny = f64::from_bits(1);
        let narrow =
            regroup_flux_flat_energy(&[0.0, 1e-10], &[f64::MAX], &[0.0, tiny, 1e-10]).unwrap();
        // The density overflows, but this integrated overlap is finite.
        assert_regroup_close(narrow.flux[0], (f64::MAX * tiny) / 1e-10);
        let wide =
            regroup_flux_flat_energy(&[0.0, f64::MAX], &[5e-16], &[0.0, 3.0, f64::MAX]).unwrap();
        // A subnormal density must not be rounded before integration.
        assert_eq!(wide.flux[0], f64::from_bits(2));
    }

    #[test]
    fn regroup_rejects_invalid_grids_on_either_side() {
        let invalid: &[&[f64]] = &[
            &[],
            &[0.0],
            &[0.0, 0.0, 2.0],
            &[0.0, 2.0, 1.0],
            &[-1.0, 2.0],
            &[0.0, f64::NAN],
            &[0.0, f64::INFINITY],
            &[f64::NEG_INFINITY, 2.0],
        ];
        for grid in invalid {
            assert!(matches!(
                regroup_flux_flat_energy(grid, &[1.0], &[0.0, 2.0]),
                Err(Error::BadFluxRegrouping(_))
            ));
            assert!(matches!(
                regroup_flux_flat_energy(&[0.0, 2.0], &[1.0], grid),
                Err(Error::BadFluxRegrouping(_))
            ));
        }
    }

    #[test]
    fn regroup_rejects_mismatched_support_without_tolerance_or_clipping() {
        for target in [
            [2.0, 3.0],
            [0.0, 3.0],
            [0.5, 2.0],
            [0.0, 2.0 + f64::EPSILON * 2.0],
        ] {
            assert!(regroup_flux_flat_energy(&[0.0, 2.0], &[1.0], &target).is_err());
        }
    }

    #[test]
    fn regroup_rejects_invalid_flux_and_overflow_even_for_identity() {
        let edges = [0.0, 1.0, 2.0];
        for values in [
            vec![],
            vec![1.0],
            vec![1.0, 2.0, 3.0],
            vec![-1.0, 0.0],
            vec![f64::NAN, 0.0],
            vec![f64::INFINITY, 0.0],
            vec![f64::NEG_INFINITY, 0.0],
            vec![f64::MAX, f64::MAX],
        ] {
            assert!(regroup_flux_flat_energy(&edges, &values, &edges).is_err());
            assert!(regroup_flux_flat_energy(&edges, &values, &[0.0, 2.0]).is_err());
        }
    }

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
