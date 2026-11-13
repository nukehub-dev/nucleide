//! Per-voxel photon-source tags for R2S (no MOAB, no HDF5).
//!
//! [`VoxelTags`] maps zone-indexed photon data onto voxels of the native
//! structured mesh: each voxel names its zone ([`zone_of_voxel`]), carries
//! one total source strength and one decay time, and supports aggregation
//! over zones. [`tag_zone_totals`] copies a zone total onto every
//! voxel of that zone (tag-as-attribute); [`split_zone_totals`]
//! distributes it conservatively instead (the voxel strengths in a zone sum
//! back to the zone total). [`split_zone_totals_by_volume`] instead assumes
//! uniform source density within each zone using explicit voxel volumes.
//! [`photon_groups_at`] selects `.photonSrc`
//! group spectra for caller-named nuclides at one cooling time and
//! [`sum_group_strengths`] adds them element-wise, so callers can feed real
//! group spectra where [`crate::photon::assemble`] leaves its uniform-split
//! placeholder. Sub-voxel discretization stays out of scope by design (it
//! needs DAGMC/MOAB geometry data this crate deliberately does not read).
//!
//! [`zone_of_voxel`]: VoxelTags::zone_of_voxel

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::photon::ZonePhotonSource;

/// Per-voxel photon-source tags over `n_voxels` voxels.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct VoxelTags {
    /// Zone count the indices in [`zone_of_voxel`](Self::zone_of_voxel)
    /// refer to.
    pub n_zones: usize,
    /// Zone index per voxel (length `n_voxels`, each `< n_zones`).
    pub zone_of_voxel: Vec<usize>,
    /// Total source strength per voxel (length `n_voxels`).
    pub source_strength: Vec<f64>,
    /// Decay (cooling) time per voxel in seconds (length `n_voxels`).
    pub decay_time_s: Vec<f64>,
}

impl VoxelTags {
    /// Voxel count.
    pub fn n_voxels(&self) -> usize {
        self.zone_of_voxel.len()
    }

    /// Build tags from parallel arrays, checking lengths, zone-index bounds,
    /// and finiteness (non-finite strengths or times are rejected; sign
    /// conventions stay caller-side).
    pub fn new(
        n_zones: usize,
        zone_of_voxel: Vec<usize>,
        source_strength: Vec<f64>,
        decay_time_s: Vec<f64>,
    ) -> Result<Self> {
        let n = zone_of_voxel.len();
        if source_strength.len() != n || decay_time_s.len() != n {
            return Err(Error::Invalid(format!(
                "voxel tag arrays must share one length, found zone_of_voxel={}, \
                 source_strength={}, decay_time_s={}",
                n,
                source_strength.len(),
                decay_time_s.len()
            )));
        }
        if let Some(bad) = zone_of_voxel.iter().position(|z| *z >= n_zones) {
            return Err(Error::CrossRef(format!(
                "voxel {bad} names zone {} of only {n_zones}",
                zone_of_voxel[bad]
            )));
        }
        if source_strength.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid(
                "voxel source strengths must be finite".to_string(),
            ));
        }
        if decay_time_s.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid(
                "voxel decay times must be finite".to_string(),
            ));
        }
        Ok(VoxelTags {
            n_zones,
            zone_of_voxel,
            source_strength,
            decay_time_s,
        })
    }

    /// Sum of all voxel strengths.
    pub fn total_strength(&self) -> f64 {
        self.source_strength.iter().sum()
    }

    /// Sum voxel strengths with compensated accumulation, rejecting
    /// non-finite inputs or aggregate arithmetic. Signed finite strengths
    /// are supported. This additive checked API leaves the ordinary
    /// summation behavior of [`Self::total_strength`] unchanged.
    ///
    /// Compensation recovers small representable contributions beside
    /// large ones, but cannot recover strengths already rounded to zero
    /// by an earlier split below the representable subnormal range.
    pub fn checked_total_strength(&self) -> Result<f64> {
        let mut sum = CheckedTagSum::default();
        for &value in &self.source_strength {
            if !value.is_finite() {
                return Err(Error::Invalid(
                    "voxel source strengths must be finite".into(),
                ));
            }
            sum.add(value);
            if !sum.total().is_finite() {
                return Err(Error::Invalid(
                    "aggregate voxel source strength must be finite".into(),
                ));
            }
        }
        Ok(sum.total())
    }

    /// Sum of strengths over voxels naming `zone`.
    pub fn zone_total(&self, zone: usize) -> f64 {
        self.zone_of_voxel
            .iter()
            .zip(&self.source_strength)
            .filter(|(z, _)| **z == zone)
            .map(|(_, v)| *v)
            .sum()
    }
}

/// Copy each zone total onto every voxel of that zone (tag-as-attribute:
/// voxel strengths in a zone sum to `count * total`, not `total`).
/// Decay times are `0.0` (shutdown); voxel count is `zone_of_voxel.len()`.
pub fn tag_zone_totals(zones: &[ZonePhotonSource], zone_of_voxel: &[usize]) -> Result<VoxelTags> {
    let totals: Vec<f64> = zones.iter().map(ZonePhotonSource::total).collect();
    let strengths: Vec<f64> = zone_of_voxel
        .iter()
        .map(|z| {
            totals.get(*z).copied().ok_or_else(|| {
                Error::CrossRef(format!(
                    "voxel names zone {} of only {} zones",
                    z,
                    totals.len()
                ))
            })
        })
        .collect::<Result<_>>()?;
    VoxelTags::new(
        zones.len(),
        zone_of_voxel.to_vec(),
        strengths,
        vec![0.0; zone_of_voxel.len()],
    )
}

/// Distribute each zone total conservatively over its voxels (voxel =
/// `total / voxel count in zone`; zones with no voxels contribute nothing).
/// The voxel strengths in a zone sum back to the zone total exactly when
/// the division is exact, and [`VoxelTags::total_strength`] equals the sum
/// of zone totals over zones owning at least one voxel.
pub fn split_zone_totals(zones: &[ZonePhotonSource], zone_of_voxel: &[usize]) -> Result<VoxelTags> {
    let totals: Vec<f64> = zones.iter().map(ZonePhotonSource::total).collect();
    let mut counts = vec![0usize; zones.len()];
    for z in zone_of_voxel {
        counts.get_mut(*z).ok_or_else(|| {
            Error::CrossRef(format!(
                "voxel names zone {z} of only {} zones",
                zones.len()
            ))
        })?;
        counts[*z] += 1;
    }
    let strengths: Vec<f64> = zone_of_voxel
        .iter()
        .map(|z| totals[*z] / counts[*z] as f64)
        .collect();
    VoxelTags::new(
        zones.len(),
        zone_of_voxel.to_vec(),
        strengths,
        vec![0.0; zone_of_voxel.len()],
    )
}

/// Distribute integrated zone strengths in photons/s in proportion to
/// caller-supplied voxel volumes in cm³, assuming uniform source density
/// within each zone. Returned strengths are integrated photons/s, not
/// densities; caller voxel order is preserved and decay times are `0.0`.
///
/// Volumes must match the mapping length and be finite and positive. Zone
/// totals must be finite and nonnegative, and every positive-strength zone
/// must own a voxel. Per-zone volume sums and output strengths must be
/// finite. An empty mesh is allowed when all zone strengths are zero.
/// Conservation holds up to floating-point roundoff. Unlike
/// [`split_zone_totals`], this path rejects unmapped positive sources and
/// negative zone totals; the legacy count-based path retains its semantics.
/// Zone strengths and volume sums use compensated accumulation. Individual
/// strengths below the representable subnormal range can round to zero,
/// so relative conservation cannot be guaranteed for underflowing sources.
pub fn split_zone_totals_by_volume(
    zones: &[ZonePhotonSource],
    zone_of_voxel: &[usize],
    voxel_volumes_cm3: &[f64],
) -> Result<VoxelTags> {
    if zone_of_voxel.len() != voxel_volumes_cm3.len() {
        return Err(Error::Invalid(format!(
            "voxel mapping and volumes must share one length, found {} and {}",
            zone_of_voxel.len(),
            voxel_volumes_cm3.len()
        )));
    }
    let totals: Vec<f64> = zones
        .iter()
        .map(|zone| {
            let mut sum = CheckedTagSum::default();
            for &group in &zone.groups {
                sum.add(group);
            }
            sum.total()
        })
        .collect();
    for (zone, total) in totals.iter().enumerate() {
        if !total.is_finite() || *total < 0.0 {
            return Err(Error::Invalid(format!(
                "zone {zone} source strength must be finite and nonnegative"
            )));
        }
    }
    let mut zone_volume_sums = vec![CheckedTagSum::default(); zones.len()];
    for (voxel, (&zone, &volume)) in zone_of_voxel.iter().zip(voxel_volumes_cm3).enumerate() {
        let sum = zone_volume_sums.get_mut(zone).ok_or_else(|| {
            Error::CrossRef(format!(
                "voxel {voxel} names zone {zone} of only {} zones",
                zones.len()
            ))
        })?;
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::Invalid(format!(
                "voxel {voxel} volume must be finite and positive"
            )));
        }
        sum.add(volume);
        if !sum.total().is_finite() {
            return Err(Error::Invalid(format!(
                "zone {zone} voxel volume sum must be finite"
            )));
        }
    }
    let zone_volumes: Vec<f64> = zone_volume_sums.iter().map(CheckedTagSum::total).collect();
    for (zone, (&total, &volume)) in totals.iter().zip(&zone_volumes).enumerate() {
        if total > 0.0 && volume == 0.0 {
            return Err(Error::Invalid(format!(
                "positive-strength zone {zone} has no voxels"
            )));
        }
    }
    let strengths = zone_of_voxel
        .iter()
        .zip(voxel_volumes_cm3)
        .map(|(&zone, &volume)| {
            // Form the bounded volume fraction first to avoid overflow in
            // total * volume for otherwise representable source strengths.
            let fraction = volume / zone_volumes[zone];
            if fraction < f64::MIN_POSITIVE {
                // Subnormal fractions lose relative precision before
                // reaching zero. Avoid rounding them before multiplication.
                let density = totals[zone] / zone_volumes[zone];
                if density.is_normal() {
                    density * volume
                } else {
                    // Under this small-fraction condition, product-first
                    // arithmetic is safe when density over/underflows.
                    (totals[zone] * volume) / zone_volumes[zone]
                }
            } else {
                totals[zone] * fraction
            }
        })
        .collect();
    VoxelTags::new(
        zones.len(),
        zone_of_voxel.to_vec(),
        strengths,
        vec![0.0; zone_of_voxel.len()],
    )
}

// Used only by the additive checked aggregate and volume path; legacy
// tagging and ordinary summation semantics remain unchanged.
#[derive(Clone, Copy, Default)]
struct CheckedTagSum {
    sum: f64,
    correction: f64,
}

impl CheckedTagSum {
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

/// Select `.photonSrc` group spectra for caller-named `nuclides` at exactly
/// `time_s` seconds (exact match, the same convention
/// [`crate::photon::assemble`] uses for shutdown `0.0`). Unknown nuclides
/// select nothing; `TOTAL` aggregates are ordinary rows — pass `"TOTAL"` to
/// select them. No rescaling is applied: strengths keep the file's own
/// normalization and stay caller-side to interpret.
pub fn photon_groups_at<'a>(
    photon: &'a nucleide_alara_io::photon::PhotonSource,
    nuclides: &[&str],
    time_s: f64,
) -> Vec<&'a nucleide_alara_io::photon::PhotonGroup> {
    photon
        .groups
        .iter()
        .filter(|g| g.time_s == time_s && nuclides.contains(&g.nuclide.as_str()))
        .collect()
}

/// Add selected group spectra element-wise (ALARA group order preserved).
/// Empty selection yields an empty spectrum; ragged group counts are an
/// error. The sums conserve the input total exactly up to float rounding:
/// `sums.iter().sum()` equals the sum of the inputs' totals.
pub fn sum_group_strengths(groups: &[&nucleide_alara_io::photon::PhotonGroup]) -> Result<Vec<f64>> {
    let mut sums: Vec<f64> = Vec::new();
    for g in groups {
        if sums.is_empty() {
            sums = g.strengths.clone();
        } else {
            if sums.len() != g.strengths.len() {
                return Err(Error::Invalid(format!(
                    "photon group spectra have ragged group counts ({} vs {})",
                    sums.len(),
                    g.strengths.len()
                )));
            }
            for (s, v) in sums.iter_mut().zip(&g.strengths) {
                *s += *v;
            }
        }
    }
    Ok(sums)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zones() -> Vec<ZonePhotonSource> {
        vec![
            ZonePhotonSource {
                zone: "a".to_string(),
                groups: vec![100.0, 300.0],
            },
            ZonePhotonSource {
                zone: "b".to_string(),
                groups: vec![50.0],
            },
        ]
    }

    #[test]
    fn tag_copies_zone_totals() {
        let tags = tag_zone_totals(&zones(), &[0, 0, 1]).unwrap();
        assert_eq!(tags.n_voxels(), 3);
        assert_eq!(tags.source_strength, vec![400.0, 400.0, 50.0]);
        assert_eq!(tags.decay_time_s, vec![0.0, 0.0, 0.0]);
        assert_eq!(tags.zone_total(0), 800.0);
        assert_eq!(tags.zone_total(1), 50.0);
        assert_eq!(tags.total_strength(), 850.0);
    }

    #[test]
    fn split_conserves_zone_totals() {
        let tags = split_zone_totals(&zones(), &[0, 0, 1]).unwrap();
        assert_eq!(tags.source_strength, vec![200.0, 200.0, 50.0]);
        assert_eq!(tags.zone_total(0), 400.0);
        assert_eq!(tags.zone_total(1), 50.0);
        assert_eq!(tags.total_strength(), 450.0);
    }

    #[test]
    fn split_zone_without_voxels_contributes_nothing() {
        let tags = split_zone_totals(&zones(), &[1, 1]).unwrap();
        assert_eq!(tags.source_strength, vec![25.0, 25.0]);
        assert_eq!(tags.total_strength(), 50.0);
    }

    fn sources(totals: &[f64]) -> Vec<ZonePhotonSource> {
        totals
            .iter()
            .enumerate()
            .map(|(i, &total)| ZonePhotonSource {
                zone: format!("zone-{i}"),
                groups: vec![total],
            })
            .collect()
    }

    fn assert_close(actual: f64, expected: f64) {
        if expected == 0.0 {
            assert_eq!(actual, 0.0);
        } else {
            assert!(
                (actual - expected).abs() <= 1e-12 * expected.abs(),
                "{actual} != {expected}"
            );
        }
    }

    #[test]
    fn checked_total_compensates_many_small_strengths() {
        let small = f64::EPSILON / 4.0;
        let mut strengths = vec![small; 100_001];
        strengths[0] = 1.0;
        let tags =
            VoxelTags::new(1, vec![0; strengths.len()], strengths, vec![0.0; 100_001]).unwrap();
        let expected = 1.0 + 100_000.0 * small;
        assert_eq!(tags.checked_total_strength().unwrap(), expected);
        // The additive API preserves legacy ordinary summation.
        assert_eq!(tags.total_strength(), 1.0);
        assert_eq!(tags.zone_total(0), 1.0);
    }

    #[test]
    fn checked_total_supports_signed_tags_and_empty_meshes() {
        let tags = VoxelTags::new(1, vec![0; 3], vec![1e16, 1.0, -1e16], vec![0.0; 3]).unwrap();
        assert_eq!(tags.checked_total_strength().unwrap(), 1.0);
        assert_eq!(tags.total_strength(), 0.0);
        assert_eq!(VoxelTags::default().checked_total_strength().unwrap(), 0.0);
    }

    #[test]
    fn checked_total_rejects_overflow_and_publicly_mutated_nonfinite_values() {
        let mut tags = VoxelTags::new(1, vec![0; 2], vec![f64::MAX; 2], vec![0.0; 2]).unwrap();
        assert!(tags.checked_total_strength().is_err());
        assert!(tags.total_strength().is_infinite());
        tags.source_strength = vec![f64::NAN, 1.0];
        assert!(tags.checked_total_strength().is_err());
    }

    #[test]
    fn volume_split_analytic_conservation_and_density() {
        let tags =
            split_zone_totals_by_volume(&sources(&[120.0, 40.0]), &[0, 0, 1], &[1.0, 3.0, 2.0])
                .unwrap();
        assert_eq!(tags.source_strength, vec![30.0, 90.0, 40.0]);
        assert_eq!(tags.zone_of_voxel, vec![0, 0, 1]);
        assert_eq!(tags.n_zones, 2);
        assert_eq!(tags.decay_time_s, vec![0.0; 3]);
        assert_eq!(tags.zone_total(0), 120.0);
        assert_eq!(tags.zone_total(1), 40.0);
        assert_eq!(tags.total_strength(), 160.0);
        assert_eq!(tags.source_strength[0] / 1.0, tags.source_strength[1] / 3.0);
    }

    #[test]
    fn volume_split_preserves_voxel_permutation() {
        let tags =
            split_zone_totals_by_volume(&sources(&[120.0, 40.0]), &[1, 0, 0], &[2.0, 3.0, 1.0])
                .unwrap();
        assert_eq!(tags.source_strength, vec![40.0, 90.0, 30.0]);
        assert_eq!(tags.zone_total(0), 120.0);
        assert_eq!(tags.zone_total(1), 40.0);
    }

    #[test]
    fn volume_split_equal_volumes_recovers_count_split() {
        let mapping = [1, 0, 0, 1, 0];
        for scale in [1e-100, 1.0, 1e100] {
            let zones = sources(&[17.0 * scale, 13.0 * scale]);
            let legacy = split_zone_totals(&zones, &mapping).unwrap();
            let tags = split_zone_totals_by_volume(&zones, &mapping, &[3.0; 5]).unwrap();
            for (&actual, &expected) in tags.source_strength.iter().zip(&legacy.source_strength) {
                assert_close(actual, expected);
            }
            assert_close(tags.zone_total(0), 17.0 * scale);
            assert_close(tags.zone_total(1), 13.0 * scale);
        }
    }

    #[test]
    fn volume_split_compensates_many_small_source_groups() {
        let small = f64::EPSILON / 4.0;
        let mut groups = vec![small; 100_001];
        groups[0] = 1.0;
        let zone = ZonePhotonSource {
            zone: "synthetic".into(),
            groups,
        };
        // Exact power-of-two terms have an independent arithmetic total.
        let expected = 1.0 + 100_000.0 * small;
        let tags = split_zone_totals_by_volume(&[zone], &[0], &[7.0]).unwrap();
        assert_eq!(tags.source_strength, [expected]);
    }

    #[test]
    fn volume_split_compensates_many_small_volumes_and_preserves_order() {
        let small = f64::EPSILON / 4.0;
        let mut volumes = vec![small; 100_001];
        volumes[0] = 1.0;
        let mapping = vec![0; volumes.len()];
        let zones = sources(&[1.0]);
        let tags = split_zone_totals_by_volume(&zones, &mapping, &volumes).unwrap();
        let exact_volume = 1.0 + 100_000.0 * small;
        let expected_large = 1.0 / exact_volume;
        let expected_small = small / exact_volume;
        assert_eq!(tags.source_strength[0], expected_large);
        assert!(tags.source_strength[1..]
            .iter()
            .all(|&v| v == expected_small));
        // Aggregate identical small outputs analytically, independently of
        // the legacy tag-total accessors' ordinary floating-point sum.
        let recovered = tags.source_strength[0] + 100_000.0 * tags.source_strength[1];
        assert!((recovered - 1.0).abs() <= f64::EPSILON);
        assert!((tags.checked_total_strength().unwrap() - 1.0).abs() <= f64::EPSILON);
        volumes.reverse();
        let permuted = split_zone_totals_by_volume(&zones, &mapping, &volumes).unwrap();
        assert!(tags
            .source_strength
            .iter()
            .rev()
            .eq(permuted.source_strength.iter()));
    }

    #[test]
    fn volume_split_one_voxel_per_zone_recovers_totals() {
        let tags =
            split_zone_totals_by_volume(&sources(&[120.0, 40.0]), &[1, 0], &[0.1, 7.0]).unwrap();
        assert_eq!(tags.source_strength, vec![40.0, 120.0]);
    }

    #[test]
    fn volume_split_zero_and_empty_mesh() {
        let zero = sources(&[0.0, 0.0]);
        let tags = split_zone_totals_by_volume(&zero, &[0, 0], &[1.0, 9.0]).unwrap();
        assert_eq!(tags.source_strength, vec![0.0, 0.0]);
        assert_eq!(tags.total_strength(), 0.0);
        let empty = split_zone_totals_by_volume(&zero, &[], &[]).unwrap();
        assert_eq!(empty.n_zones, 2);
        assert_eq!(empty.n_voxels(), 0);
        assert_eq!(empty.total_strength(), 0.0);
        assert!(split_zone_totals_by_volume(&[], &[], &[]).is_ok());
        assert!(split_zone_totals_by_volume(&sources(&[1.0]), &[], &[]).is_err());
        assert!(split_zone_totals_by_volume(&sources(&[1.0, 0.0]), &[1], &[1.0]).is_err());
    }

    #[test]
    fn volume_split_rejects_invalid_inputs() {
        let zones = sources(&[1.0]);
        for volumes in [vec![], vec![1.0, 2.0]] {
            assert!(split_zone_totals_by_volume(&zones, &[0], &volumes).is_err());
        }
        for volume in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(split_zone_totals_by_volume(&zones, &[0], &[volume]).is_err());
        }
        assert!(matches!(
            split_zone_totals_by_volume(&zones, &[1], &[1.0]),
            Err(Error::CrossRef(_))
        ));
        for total in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(split_zone_totals_by_volume(&sources(&[total]), &[0], &[1.0]).is_err());
        }
        // Validate even unreferenced zone strengths.
        assert!(split_zone_totals_by_volume(&sources(&[f64::NAN]), &[], &[]).is_err());
        assert!(split_zone_totals_by_volume(&zones, &[0, 0], &[f64::MAX, f64::MAX]).is_err());
        let overflow_total = vec![ZonePhotonSource {
            zone: "overflow".to_string(),
            groups: vec![f64::MAX, f64::MAX],
        }];
        assert!(split_zone_totals_by_volume(&overflow_total, &[0], &[1.0]).is_err());
    }

    #[test]
    fn volume_split_avoids_spurious_intermediate_overflow_and_underflow() {
        let tags = split_zone_totals_by_volume(&sources(&[f64::MAX]), &[0], &[f64::MAX]).unwrap();
        assert_eq!(tags.source_strength, vec![f64::MAX]);
        let tiny = f64::from_bits(1);
        let tags =
            split_zone_totals_by_volume(&sources(&[f64::MAX]), &[0, 0], &[tiny, f64::MAX]).unwrap();
        assert_eq!(tags.source_strength, vec![tiny, f64::MAX]);
    }

    #[test]
    fn volume_split_unrepresentable_halves_round_to_zero() {
        let minimum_subnormal = f64::from_bits(1);
        let tags =
            split_zone_totals_by_volume(&sources(&[minimum_subnormal]), &[0, 0], &[1.0, 1.0])
                .unwrap();
        // Each exact half lies halfway between zero and the minimum
        // subnormal. Round-to-nearest, ties-to-even selects zero. Relative
        // conservation is unachievable for these independently rounded
        // halves; compensation cannot recreate an unrepresentable value.
        assert_eq!(tags.source_strength, [0.0, 0.0]);
        assert_eq!(tags.checked_total_strength().unwrap(), 0.0);
        assert!(minimum_subnormal > 0.0);
    }

    #[test]
    fn volume_split_preserves_accuracy_for_nonzero_subnormal_fraction() {
        let tags = split_zone_totals_by_volume(&sources(&[f64::MAX]), &[0, 0], &[5e-16, f64::MAX])
            .unwrap();
        // Unit source density gives strengths equal to voxel volumes.
        assert_eq!(tags.source_strength, vec![5e-16, f64::MAX]);
    }

    #[test]
    fn volume_split_subnormal_fraction_handles_non_normal_density() {
        let tiny = f64::from_bits(1);
        let narrow =
            split_zone_totals_by_volume(&sources(&[f64::MAX]), &[0, 0], &[tiny, 1e-10]).unwrap();
        assert_close(narrow.source_strength[0], (f64::MAX * tiny) / 1e-10);
        let wide =
            split_zone_totals_by_volume(&sources(&[5e-16]), &[0, 0], &[3.0, f64::MAX]).unwrap();
        assert_eq!(wide.source_strength[0], f64::from_bits(2));
    }

    #[test]
    fn legacy_split_still_allows_signed_and_unmapped_sources() {
        let tags = split_zone_totals(&sources(&[-12.0, 40.0]), &[0, 0]).unwrap();
        assert_eq!(tags.source_strength, vec![-6.0, -6.0]);
    }

    #[test]
    fn dangling_zone_index_errors() {
        assert!(matches!(
            tag_zone_totals(&zones(), &[0, 7]),
            Err(Error::CrossRef(_))
        ));
        assert!(matches!(
            split_zone_totals(&zones(), &[2]),
            Err(Error::CrossRef(_))
        ));
        assert!(matches!(
            VoxelTags::new(1, vec![0], vec![1.0, 2.0], vec![0.0, 0.0]),
            Err(Error::Invalid(_))
        ));
        assert!(matches!(
            VoxelTags::new(1, vec![0], vec![f64::NAN], vec![0.0]),
            Err(Error::Invalid(_))
        ));
    }

    fn photon_source() -> nucleide_alara_io::photon::PhotonSource {
        nucleide_alara_io::photon::PhotonSource::from_str(
            "mn-56 shutdown 6.0 2.0\n\
             co-60 shutdown 1.0 3.0\n\
             mn-56 1 h 5.0 1.0\n\
             TOTAL shutdown 7.0 5.0\n",
        )
        .unwrap()
    }

    #[test]
    fn groups_select_nuclides_at_time() {
        let photon = photon_source();
        let at = photon_groups_at(&photon, &["mn-56", "co-60"], 0.0);
        assert_eq!(at.len(), 2);
        assert_eq!(at[0].nuclide, "mn-56");
        let sums = sum_group_strengths(&at).unwrap();
        assert_eq!(sums, vec![7.0, 5.0]);
        // Sums conserve the input total.
        let inputs: f64 = at.iter().map(|g| g.total_strength()).sum();
        assert_eq!(sums.iter().sum::<f64>(), inputs);
        // Cooling-time rows and TOTAL aggregates select independently.
        assert_eq!(photon_groups_at(&photon, &["mn-56"], 3600.0).len(), 1);
        assert_eq!(photon_groups_at(&photon, &["TOTAL"], 0.0).len(), 1);
        assert!(photon_groups_at(&photon, &["ghost"], 0.0).is_empty());
        assert!(sum_group_strengths(&[]).unwrap().is_empty());
    }

    #[test]
    fn ragged_groups_error() {
        let photon = nucleide_alara_io::photon::PhotonSource::from_str(
            "mn-56 shutdown 6.0 2.0\nco-60 shutdown 1.0\n",
        )
        .unwrap();
        let at = photon_groups_at(&photon, &["mn-56", "co-60"], 0.0);
        assert!(matches!(sum_group_strengths(&at), Err(Error::Invalid(_))));
    }
}
