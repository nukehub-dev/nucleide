//! Decay photon-source assembly for R2S.
//!
//! [`assemble`] builds a per-zone photon source from an ALARA activation
//! listing ([`nucleide_alara_io::output::ResponseFrame`]); [`from_photon_file`] loads
//! the group-wise spectra ALARA writes to `.photonSrc` files
//! ([`nucleide_alara_io::photon::PhotonSource`]).

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Photon source per mesh zone.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ZonePhotonSource {
    /// Zone name.
    pub zone: String,
    /// Group strengths.
    pub groups: Vec<f64>,
}

impl ZonePhotonSource {
    /// Total photon strength.
    pub fn total(&self) -> f64 {
        self.groups.iter().sum()
    }
}

/// Assemble the shutdown photon source for `zone` from an activation listing.
///
/// Sums [`nucleide_alara_io::output::ResponseVar::SpecificActivity`] values over every
/// nuclide row whose block name equals `zone` (any block resolution: `Zone`
/// and `Interval` listings both carry the zone name) at shutdown (`time_s ==
/// 0.0`). `total` aggregate rows are skipped so nuclides are counted once.
/// The sum is split uniformly over `groups` energy groups (`groups == 0`
/// yields an empty spectrum).
///
/// # Approximation
///
/// The uniform split preserves only the total shutdown strength: real decay
/// photons follow the nuclide- and energy-dependent emission lines tabulated
/// in the `.photonSrc` spectra (see [`from_photon_file`]). Treat the uniform
/// spectrum as a strength-conserving placeholder until group-wise emission
/// data is wired through per nuclide.
pub fn assemble(
    frame: &nucleide_alara_io::output::ResponseFrame,
    zone: &str,
    groups: usize,
) -> ZonePhotonSource {
    let total: f64 = frame
        .rows
        .iter()
        .filter(|row| row.variable == nucleide_alara_io::output::ResponseVar::SpecificActivity)
        .filter(|row| row.time_s == 0.0)
        .filter(|row| row.block_name == zone)
        .filter(|row| !row.is_total())
        .map(|row| row.value)
        .sum();
    let spectrum = if groups == 0 {
        Vec::new()
    } else {
        vec![total / groups as f64; groups]
    };
    ZonePhotonSource {
        zone: zone.to_string(),
        groups: spectrum,
    }
}

/// Load group-wise photon spectra from an ALARA `.photonSrc` file.
///
/// Delegates to [`nucleide_alara_io::photon::PhotonSource::from_file`]; unlike
/// [`assemble`] this preserves the per-group, per-nuclide, per-cooling-time
/// structure ALARA computed.
pub fn from_photon_file(
    path: impl AsRef<std::path::Path>,
) -> Result<nucleide_alara_io::photon::PhotonSource> {
    nucleide_alara_io::photon::PhotonSource::from_file(path).map_err(|error| match error {
        nucleide_alara_io::Error::Io(inner) => Error::Io(inner),
        nucleide_alara_io::Error::CrossRef(message) => Error::CrossRef(message),
        other => Error::Invalid(other.to_string()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nucleide_alara_io::output::{BlockKind, ResponseRow, ResponseVar};

    fn row(
        nuclide: &str,
        variable: ResponseVar,
        block_name: &str,
        time_s: f64,
        time_label: &str,
        value: f64,
    ) -> ResponseRow {
        ResponseRow {
            time_s,
            time_label: time_label.to_string(),
            nuclide: nuclide.to_string(),
            half_life_s: 1.0,
            run_lbl: "test".to_string(),
            block: BlockKind::Zone,
            block_name: block_name.to_string(),
            block_num: 1,
            variable,
            var_unit: "Bq/cm3".to_string(),
            value,
        }
    }

    fn frame() -> nucleide_alara_io::output::ResponseFrame {
        nucleide_alara_io::output::ResponseFrame {
            rows: vec![
                row(
                    "mn-56",
                    ResponseVar::SpecificActivity,
                    "zone_a",
                    0.0,
                    "shutdown",
                    100.0,
                ),
                row(
                    "co-60",
                    ResponseVar::SpecificActivity,
                    "zone_a",
                    0.0,
                    "shutdown",
                    300.0,
                ),
                // `total` aggregate: skipped (would double-count).
                row(
                    "total",
                    ResponseVar::SpecificActivity,
                    "zone_a",
                    0.0,
                    "shutdown",
                    400.0,
                ),
                // Wrong zone, wrong time, wrong variable: all skipped.
                row(
                    "mn-56",
                    ResponseVar::SpecificActivity,
                    "zone_b",
                    0.0,
                    "shutdown",
                    50.0,
                ),
                row(
                    "mn-56",
                    ResponseVar::SpecificActivity,
                    "zone_a",
                    86_400.0,
                    "1 d",
                    60.0,
                ),
                row(
                    "mn-56",
                    ResponseVar::NumberDensity,
                    "zone_a",
                    0.0,
                    "shutdown",
                    1e20,
                ),
            ],
        }
    }

    #[test]
    fn total_sums_groups() {
        let src = ZonePhotonSource {
            zone: "a".to_string(),
            groups: vec![1.0, 2.0],
        };
        assert_eq!(src.total(), 3.0);
    }

    #[test]
    fn assemble_sums_shutdown_activity_uniformly() {
        let source = assemble(&frame(), "zone_a", 4);
        assert_eq!(source.zone, "zone_a");
        assert_eq!(source.groups.len(), 4);
        // 100 + 300; total/total-label rows excluded.
        assert!(source.groups.iter().all(|value| *value == 100.0));
        assert_eq!(source.total(), 400.0);
    }

    #[test]
    fn assemble_matches_other_zones_exactly() {
        let source = assemble(&frame(), "zone_b", 2);
        assert_eq!(source.groups, [25.0, 25.0]);
    }

    #[test]
    fn assemble_empty_zone_and_zero_groups() {
        let source = assemble(&frame(), "missing", 3);
        assert_eq!(source.groups, [0.0, 0.0, 0.0]);
        assert_eq!(source.total(), 0.0);

        let source = assemble(&frame(), "zone_a", 0);
        assert!(source.groups.is_empty());
    }

    #[test]
    fn from_photon_file_delegates_to_alara_photon() {
        let dir = std::env::temp_dir().join(format!("r2s_photon_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.photonSrc");
        std::fs::write(
            &path,
            "mn-56 \tshutdown\t0\t6.0\t2.0\nTOTAL\tshutdown\t0\t6.0\t2.0\n",
        )
        .unwrap();
        let source = from_photon_file(&path).unwrap();
        assert_eq!(source.len(), 2);
        assert_eq!(source.total_strength(), 16.0);
        assert!(matches!(
            from_photon_file(dir.join("missing.photonSrc")),
            Err(Error::Io(_))
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
