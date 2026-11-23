//! Allocate activation-zone photon spectra by zone/voxel intersection volume.

use std::collections::HashSet;

use crate::source::{validate_source_axes, CompensatedSum};
use crate::tags::share_by_volume;
use crate::{Error, PhotonSource, Result};

/// Relative tolerance for comparing declared and mapped volumes of emitting
/// zones, and for comparing supplied full voxel volumes to regular-mesh
/// geometry. This tolerance also bounds any excess in voxel intersection sums.
pub const MIXED_VOXEL_VOLUME_RELATIVE_TOLERANCE: f64 = 1e-12;

/// One activation-zone spectrum on a selected energy grid and cooling time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixedVoxelZoneSpectrum<'a> {
    /// Independently declared total activation-zone volume in cm³.
    pub volume_cm3: f64,
    /// Full ascending energy-bin edges in eV.
    pub edges_ev: &'a [f64],
    /// Cooling offset in seconds after shutdown.
    pub cooling_s: f64,
    /// Integrated photon strength per group in photons/s.
    pub strengths_photons_per_s: &'a [f64],
}

/// Observed overlap coverage and energy-group conservation for an allocation.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedVoxelAllocation {
    /// Relative tolerance applied to volume comparisons.
    pub volume_relative_tolerance: f64,
    /// Caller-declared zone volumes in input zone order.
    pub zone_declared_volumes_cm3: Vec<f64>,
    /// Sum of intersection volumes by zone.
    pub zone_mapped_volumes_cm3: Vec<f64>,
    /// Mapped minus declared volume by zone, in cm³.
    pub zone_volume_deltas_cm3: Vec<f64>,
    /// Full volume minus mapped volume in flat C-order, in cm³.
    /// A negative residual within tolerance is reported as zero.
    pub voxel_unassigned_volumes_cm3: Vec<f64>,
    /// Input group strengths in `[zone][group]` order.
    pub zone_input_strengths_photons_per_s: Vec<Vec<f64>>,
    /// Allocated group strengths in `[zone][group]` order.
    pub zone_allocated_strengths_photons_per_s: Vec<Vec<f64>>,
    /// Input photon strength per group across all zones.
    pub input_group_totals_photons_per_s: Vec<f64>,
    /// Output photon strength per group across all voxels.
    pub allocated_group_totals_photons_per_s: Vec<f64>,
}

/// Validated photon source and its allocation diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedVoxelPhotonSource {
    /// Dense transport source in flat C-order.
    pub source: PhotonSource,
    /// Coverage and conservation of the zone-to-voxel mapping.
    pub allocation: MixedVoxelAllocation,
}

/// Distribute each zone's photon spectrum by caller-supplied intersection
/// volumes and sum the contributions within each voxel.
///
/// `overlaps` contains unique `(voxel_index, zone_index, overlap_cm3)`
/// triples with zero-based indices. Voxels are in flat C-order and zones
/// are in caller order. Callers supply intersections of disjoint zones.
/// The sum of intersections in a voxel cannot exceed its full volume
/// beyond the documented relative tolerance.
/// Emitting zones must be fully mapped within the documented relative
/// tolerance, and every zone's grid and cooling time must match exactly.
/// Zero-strength zones may have incomplete coverage.
///
/// Unassigned volume contributes no photons before homogenization. The
/// output uses uniform position within each whole voxel and uniform
/// density per eV within each energy group. Independently rounded or
/// underflowing contributions may differ from the input photon strength.
#[allow(clippy::too_many_arguments)]
pub fn mixed_voxel_photon_source(
    mesh_shape: [usize; 3],
    lower_left_cm: [f64; 3],
    upper_right_cm: [f64; 3],
    voxel_volumes_cm3: Vec<f64>,
    edges_ev: Vec<f64>,
    cooling_s: f64,
    zones: &[MixedVoxelZoneSpectrum<'_>],
    overlaps: &[(usize, usize, f64)],
) -> Result<MixedVoxelPhotonSource> {
    validate_source_axes(
        mesh_shape,
        lower_left_cm,
        upper_right_cm,
        &voxel_volumes_cm3,
        &edges_ev,
        cooling_s,
    )?;
    validate_voxel_geometry(
        mesh_shape,
        lower_left_cm,
        upper_right_cm,
        &voxel_volumes_cm3,
    )?;
    let groups = edges_ev.len() - 1;
    let voxels = voxel_volumes_cm3.len();
    let mut input_totals = vec![CompensatedSum::default(); groups];
    for (zone_idx, zone) in zones.iter().enumerate() {
        validate_zone(zone_idx, zone, &edges_ev, cooling_s)?;
        for (sum, &strength) in input_totals.iter_mut().zip(zone.strengths_photons_per_s) {
            sum.add(strength);
        }
    }
    let input_group_totals_photons_per_s: Vec<_> =
        input_totals.iter().map(CompensatedSum::total).collect();
    if input_group_totals_photons_per_s
        .iter()
        .any(|sum| !sum.is_finite())
    {
        return Err(Error::Invalid("an input group total is not finite".into()));
    }

    let mut seen_pairs = HashSet::with_capacity(overlaps.len());
    let mut voxel_volume_sums = vec![CompensatedSum::default(); voxels];
    let mut zone_volume_sums = vec![CompensatedSum::default(); zones.len()];
    for (row, &(voxel, zone, volume)) in overlaps.iter().enumerate() {
        if voxel >= voxels || zone >= zones.len() {
            return Err(Error::CrossRef(format!(
                "overlap row {row} references voxel {voxel}, zone {zone}, outside the mesh or zone list"
            )));
        }
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::Invalid(format!(
                "overlap row {row} volume must be finite and positive"
            )));
        }
        if !seen_pairs.insert((voxel, zone)) {
            return Err(Error::Invalid(format!(
                "duplicate overlap for voxel {voxel}, zone {zone}"
            )));
        }
        voxel_volume_sums[voxel].add(volume);
        zone_volume_sums[zone].add(volume);
    }
    let voxel_unassigned_volumes_cm3 = voxel_volume_sums
        .iter()
        .zip(&voxel_volumes_cm3)
        .enumerate()
        .map(|(voxel, (sum, &volume))| {
            let mapped = sum.total();
            if !mapped.is_finite() || (mapped > volume && !volumes_match(mapped, volume)) {
                return Err(Error::Invalid(format!(
                    "voxel {voxel} overlaps fill {mapped} cm3; full volume is {volume} cm3"
                )));
            }
            Ok((volume - mapped).max(0.0))
        })
        .collect::<Result<Vec<_>>>()?;
    let zone_mapped_volumes_cm3 = zone_volume_sums
        .iter()
        .zip(zones)
        .enumerate()
        .map(|(zone_idx, (sum, zone))| {
            let mapped = sum.total();
            let emitting = zone.strengths_photons_per_s.iter().any(|&s| s > 0.0);
            if !mapped.is_finite()
                || ((emitting || mapped > zone.volume_cm3)
                    && !volumes_match(mapped, zone.volume_cm3))
            {
                return Err(Error::Invalid(format!(
                    "zone {zone_idx} coverage is {mapped} cm3 but its declared volume is {} cm3",
                    zone.volume_cm3
                )));
            }
            Ok(mapped)
        })
        .collect::<Result<Vec<_>>>()?;

    let mut voxel_strengths = new_sums(voxels, groups)?;
    let mut zone_strengths = new_sums(zones.len(), groups)?;
    for &(voxel, zone, volume) in overlaps {
        let voxel_groups = &mut voxel_strengths[voxel * groups..][..groups];
        let zone_groups = &mut zone_strengths[zone * groups..][..groups];
        for ((voxel_sum, zone_sum), &strength) in voxel_groups
            .iter_mut()
            .zip(zone_groups)
            .zip(zones[zone].strengths_photons_per_s)
        {
            let contribution = share_by_volume(strength, volume, zone_mapped_volumes_cm3[zone]);
            voxel_sum.add(contribution);
            zone_sum.add(contribution);
        }
    }
    let source = PhotonSource::new(
        mesh_shape,
        lower_left_cm,
        upper_right_cm,
        voxel_volumes_cm3,
        edges_ev,
        cooling_s,
        rows_from_sums(voxel_strengths, groups),
    )?;
    let zone_allocated_strengths_photons_per_s = rows_from_sums(zone_strengths, groups);
    if zone_allocated_strengths_photons_per_s
        .iter()
        .flatten()
        .any(|strength| !strength.is_finite())
    {
        return Err(Error::Invalid(
            "an allocated zone group strength is not finite".into(),
        ));
    }
    let allocation = MixedVoxelAllocation {
        volume_relative_tolerance: MIXED_VOXEL_VOLUME_RELATIVE_TOLERANCE,
        zone_declared_volumes_cm3: zones.iter().map(|zone| zone.volume_cm3).collect(),
        zone_volume_deltas_cm3: zone_mapped_volumes_cm3
            .iter()
            .zip(zones)
            .map(|(&mapped, zone)| mapped - zone.volume_cm3)
            .collect(),
        zone_mapped_volumes_cm3,
        voxel_unassigned_volumes_cm3,
        zone_input_strengths_photons_per_s: zones
            .iter()
            .map(|zone| zone.strengths_photons_per_s.to_vec())
            .collect(),
        zone_allocated_strengths_photons_per_s,
        input_group_totals_photons_per_s,
        allocated_group_totals_photons_per_s: source.group_totals()?,
    };
    Ok(MixedVoxelPhotonSource { source, allocation })
}

fn validate_voxel_geometry(
    mesh_shape: [usize; 3],
    lower_left_cm: [f64; 3],
    upper_right_cm: [f64; 3],
    voxel_volumes_cm3: &[f64],
) -> Result<()> {
    let widths = [0, 1, 2]
        .map(|axis| (upper_right_cm[axis] - lower_left_cm[axis]) / mesh_shape[axis] as f64);
    let geometric_volume = widths[0] * widths[1] * widths[2];
    if !geometric_volume.is_finite() || geometric_volume <= 0.0 {
        return Err(Error::Invalid(
            "regular-mesh cell volume must be finite and positive".into(),
        ));
    }
    for (voxel, &volume) in voxel_volumes_cm3.iter().enumerate() {
        if !volumes_match(volume, geometric_volume) {
            return Err(Error::Invalid(format!(
                "voxel {voxel} full volume {volume} cm3 disagrees with the regular-mesh cell volume {geometric_volume} cm3"
            )));
        }
    }
    Ok(())
}

fn validate_zone(
    index: usize,
    zone: &MixedVoxelZoneSpectrum<'_>,
    edges_ev: &[f64],
    cooling_s: f64,
) -> Result<()> {
    if !zone.volume_cm3.is_finite() || zone.volume_cm3 <= 0.0 {
        return Err(Error::Invalid(format!(
            "zone {index} declared volume must be finite and positive"
        )));
    }
    if zone.edges_ev != edges_ev {
        return Err(Error::Invalid(format!(
            "zone {index} energy grid does not match the requested energy grid"
        )));
    }
    if zone.cooling_s != cooling_s {
        return Err(Error::Invalid(format!(
            "zone {index} cooling time {} s does not match the requested cooling time {cooling_s} s",
            zone.cooling_s
        )));
    }
    let groups = edges_ev.len() - 1;
    if zone.strengths_photons_per_s.len() != groups {
        return Err(Error::Invalid(format!(
            "zone {index} spectrum must contain {groups} group strengths"
        )));
    }
    for (group, &strength) in zone.strengths_photons_per_s.iter().enumerate() {
        if !strength.is_finite() || strength < 0.0 {
            return Err(Error::Invalid(format!(
                "zone {index} group {group} strength must be finite and nonnegative"
            )));
        }
    }
    Ok(())
}

fn volumes_match(lhs: f64, rhs: f64) -> bool {
    (lhs - rhs).abs() / lhs.max(rhs) <= MIXED_VOXEL_VOLUME_RELATIVE_TOLERANCE
}

fn new_sums(rows: usize, groups: usize) -> Result<Vec<CompensatedSum>> {
    let cells = rows
        .checked_mul(groups)
        .ok_or_else(|| Error::Invalid("photon source dimensions overflow".into()))?;
    let mut sums = Vec::new();
    sums.try_reserve_exact(cells)
        .map_err(|_| Error::Invalid("cannot allocate the dense photon group table".into()))?;
    sums.resize(cells, CompensatedSum::default());
    Ok(sums)
}

fn rows_from_sums(sums: Vec<CompensatedSum>, groups: usize) -> Vec<Vec<f64>> {
    sums.chunks_exact(groups)
        .map(|row| row.iter().map(CompensatedSum::total).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{split_zone_totals_by_volume, ZonePhotonSource};

    const EDGES: [f64; 3] = [0.0, 1.0e6, 3.0e6];
    const COOLING_S: f64 = 3600.0;

    fn zone(volume_cm3: f64, strengths_photons_per_s: &[f64]) -> MixedVoxelZoneSpectrum<'_> {
        MixedVoxelZoneSpectrum {
            volume_cm3,
            edges_ev: &EDGES,
            cooling_s: COOLING_S,
            strengths_photons_per_s,
        }
    }

    fn source(
        zones: &[MixedVoxelZoneSpectrum<'_>],
        overlaps: &[(usize, usize, f64)],
    ) -> Result<MixedVoxelPhotonSource> {
        mixed_voxel_photon_source(
            [2, 1, 1],
            [0.0; 3],
            [4.0, 1.0, 1.0],
            vec![2.0; 2],
            EDGES.to_vec(),
            COOLING_S,
            zones,
            overlaps,
        )
    }

    #[test]
    fn mixed_zone_spectra_conserve_groups_and_report_unassigned_volume() {
        let result = source(
            &[zone(2.0, &[40.0, 0.0]), zone(1.0, &[0.0, 60.0])],
            &[(1, 0, 1.5), (0, 1, 1.0), (0, 0, 0.5)],
        )
        .unwrap();
        assert_eq!(result.source.cooling_s, COOLING_S);
        assert_eq!(result.source.edges_ev, EDGES);
        assert_eq!(
            result.source.strengths_photons_per_s,
            [vec![10.0, 60.0], vec![30.0, 0.0]]
        );
        assert_eq!(
            result.source.density_photons_per_s_cm3().unwrap(),
            [vec![5.0, 30.0], vec![15.0, 0.0]]
        );
        assert_eq!(result.allocation.zone_declared_volumes_cm3, [2.0, 1.0]);
        assert_eq!(result.allocation.zone_mapped_volumes_cm3, [2.0, 1.0]);
        assert_eq!(result.allocation.zone_volume_deltas_cm3, [0.0, 0.0]);
        assert_eq!(result.allocation.voxel_unassigned_volumes_cm3, [0.5, 0.5]);
        assert_eq!(
            result.allocation.zone_allocated_strengths_photons_per_s,
            result.allocation.zone_input_strengths_photons_per_s
        );
        assert_eq!(
            result.allocation.input_group_totals_photons_per_s,
            [40.0, 60.0]
        );
        assert_eq!(
            result.allocation.allocated_group_totals_photons_per_s,
            [40.0, 60.0]
        );
    }

    #[test]
    fn single_zone_per_voxel_matches_legacy_volume_split() {
        let groups = [40.0, 60.0];
        let result = source(&[zone(4.0, &groups)], &[(0, 0, 2.0), (1, 0, 2.0)]).unwrap();
        let tags = split_zone_totals_by_volume(
            &[ZonePhotonSource {
                groups: groups.to_vec(),
                ..ZonePhotonSource::default()
            }],
            &[0, 0],
            &[2.0, 2.0],
        )
        .unwrap();
        assert_eq!(result.source.voxel_totals().unwrap(), tags.source_strength);
    }

    #[test]
    fn sparse_overlaps_preserve_asymmetric_c_order() {
        let result = mixed_voxel_photon_source(
            [2, 3, 2],
            [0.0; 3],
            [4.0, 3.0, 2.0],
            vec![2.0; 12],
            EDGES.to_vec(),
            COOLING_S,
            &[zone(1.0, &[0.0, 8.0])],
            &[(9, 0, 1.0)],
        )
        .unwrap();
        let mut expected = vec![vec![0.0; 2]; 12];
        // Voxel (1, 1, 1) on a [2, 3, 2] mesh.
        expected[9][1] = 8.0;
        assert_eq!(result.source.strengths_photons_per_s, expected);
    }

    #[test]
    fn only_emitting_zones_require_full_coverage() {
        let result = source(&[zone(2.0, &[0.0, 0.0])], &[(1, 0, 1.0)]).unwrap();
        assert_eq!(result.source.checked_total_photons_per_s().unwrap(), 0.0);
        assert_eq!(result.allocation.zone_volume_deltas_cm3, [-1.0]);
        assert_eq!(result.allocation.voxel_unassigned_volumes_cm3, [2.0, 1.0]);

        let result = source(&[], &[]).unwrap();
        assert_eq!(result.source.checked_total_photons_per_s().unwrap(), 0.0);
        assert_eq!(result.allocation.voxel_unassigned_volumes_cm3, [2.0, 2.0]);
    }

    #[test]
    fn coverage_tolerance_scales_and_reports_the_volume_difference() {
        // The same relative discrepancy is accepted at very different scales.
        for declared in [1.0e-50, 1.0] {
            let mapped = declared * (1.0 - 5e-13);
            let result = source(&[zone(declared, &[4.0, 0.0])], &[(0, 0, mapped)]).unwrap();
            assert_eq!(
                result.allocation.zone_volume_deltas_cm3,
                [mapped - declared]
            );
            assert_eq!(
                result.allocation.zone_allocated_strengths_photons_per_s[0][0],
                4.0
            );
        }
    }

    #[test]
    fn accepts_roundoff_overfill_and_preserves_positive_residual() {
        let result = mixed_voxel_photon_source(
            [2, 1, 1],
            [0.0; 3],
            [0.6, 1.0, 1.0],
            vec![0.3; 2],
            EDGES.to_vec(),
            COOLING_S,
            &[zone(0.1, &[4.0, 0.0]), zone(0.2, &[0.0, 9.0])],
            &[(0, 0, 0.1), (0, 1, 0.2)],
        )
        .unwrap();
        assert_eq!(result.allocation.voxel_unassigned_volumes_cm3, [0.0, 0.3]);
        assert_eq!(
            result.source.strengths_photons_per_s,
            [vec![4.0, 9.0], vec![0.0, 0.0]]
        );

        let mapped = 2.0 - 1e-13;
        let result = source(&[zone(mapped, &[4.0, 0.0])], &[(0, 0, mapped)]).unwrap();
        assert_eq!(
            result.allocation.voxel_unassigned_volumes_cm3[0],
            2.0 - mapped
        );
    }

    #[test]
    fn rejects_incomplete_coverage_duplicate_pairs_and_overfill() {
        let zones = [zone(2.0, &[40.0, 0.0]), zone(1.0, &[0.0, 60.0])];
        let error = source(&zones, &[(0, 0, 1.0)]).unwrap_err();
        assert!(error.to_string().contains("zone 0 coverage"));
        let error = source(&zones, &[(0, 0, 1.0), (0, 0, 1.0)]).unwrap_err();
        assert!(error.to_string().contains("duplicate overlap"));
        let error = source(&zones, &[(0, 0, 1.5), (0, 1, 1.0)]).unwrap_err();
        assert!(error.to_string().contains("voxel 0 overlaps fill"));
        let error = source(&zones, &[(0, 0, 2.0)]).unwrap_err();
        assert!(error.to_string().contains("zone 1 coverage is 0"));
        let error = source(&zones, &[(2, 0, 1.0)]).unwrap_err();
        assert!(matches!(error, Error::CrossRef(_)));
        let error = source(&zones, &[(0, 0, f64::NAN)]).unwrap_err();
        assert!(error.to_string().contains("overlap row 0 volume"));
    }

    #[test]
    fn rejects_incompatible_zone_grid_time_or_strength() {
        let mut bad = zone(2.0, &[40.0, 0.0]);
        bad.edges_ev = &[0.0, 1.5e6, 3.0e6];
        let error = source(&[bad], &[]).unwrap_err();
        assert!(error.to_string().contains("zone 0 energy grid"));

        bad.edges_ev = &EDGES;
        bad.cooling_s = COOLING_S + 1e-11;
        let error = source(&[bad], &[]).unwrap_err();
        assert!(error.to_string().contains("zone 0 cooling time"));

        bad.cooling_s = COOLING_S;
        bad.strengths_photons_per_s = &[-1.0, 0.0];
        let error = source(&[bad], &[]).unwrap_err();
        assert!(error.to_string().contains("zone 0 group 0 strength"));
    }

    #[test]
    fn rejects_voxel_volumes_incompatible_with_regular_mesh() {
        let error = mixed_voxel_photon_source(
            [2, 1, 1],
            [0.0; 3],
            [4.0, 1.0, 1.0],
            vec![1.0, 3.0],
            EDGES.to_vec(),
            COOLING_S,
            &[],
            &[],
        )
        .unwrap_err();
        assert!(error.to_string().contains("regular-mesh cell volume"));
    }

    #[test]
    fn compensates_small_zone_contributions_in_the_same_group() {
        let small = 2.0_f64.powi(-54);
        let tiny_strength = [small, 0.0];
        let zones = [
            zone(0.5, &[1.0, 0.0]),
            zone(0.5, &tiny_strength),
            zone(0.5, &tiny_strength),
            zone(0.5, &tiny_strength),
        ];
        let result = source(
            &zones,
            &[(0, 0, 0.5), (0, 1, 0.5), (0, 2, 0.5), (0, 3, 0.5)],
        )
        .unwrap();
        assert_eq!(
            result.source.strengths_photons_per_s[0][0],
            1.0 + 3.0 * small
        );
        assert_eq!(
            result.allocation.allocated_group_totals_photons_per_s[0],
            1.0 + 3.0 * small
        );
    }
}
