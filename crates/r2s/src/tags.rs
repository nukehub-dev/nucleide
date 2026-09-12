//! Per-voxel photon-source tags for R2S (no MOAB, no HDF5).
//!
//! [`VoxelTags`] maps zone-indexed photon data onto voxels of the native
//! structured mesh: each voxel names its zone ([`zone_of_voxel`]), carries
//! one total source strength and one decay time, and aggregates back to
//! zone totals exactly. [`tag_zone_totals`] copies a zone total onto every
//! voxel of that zone (tag-as-attribute); [`split_zone_totals`]
//! distributes it conservatively instead (the voxel strengths in a zone sum
//! back to the zone total). [`photon_groups_at`] selects `.photonSrc`
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
