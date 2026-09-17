//! MCNP `SDEF` source-card emission through the typed reader.
//!
//! The card is built as a typed [`SdefProblem`] and rendered with its
//! canonical emitter, so every emitted card round-trips byte-identically
//! through `nucleide-mcnp-io`'s `SDEF` reader — the same verification the
//! spectroscopy decay-source emitter (E9) uses. What the reader cannot
//! express (a radial distribution) was added to its accepted subset as the
//! `AXS`/`RAD`/`EXT` keywords for this crate; anything still outside the
//! subset stays a loud reader error, never a silently misread card.
//!
//! Card shape:
//!
//! - point source: `SDEF POS=x y z`, `ERG=<e>` (monoenergetic) or
//!   `ERG=Dn` with the tabulated spectrum on `SIn L` / `SPn D`;
//! - ring source: `SDEF POS=0 0 z`, `AXS=0 0 1`, `RAD=D1` with the delta
//!   ring on `SI1 L R R` / `SP1 D 0 1` (equal `SI` endpoints pin every
//!   particle to radius `R`), plus the energy field as above.
//!
//! The Gaussian spectrum is tabulated at bin centers over a symmetric
//! `±4 sigma` window (see [`SpectrumSpec::tabulate`]); the truncated tail
//! mass is reported as drift. Sampling semantics of the `SDEF` distributions
//! themselves are MCNP's — the workspace verifies surface syntax only, the
//! same posture as the E9 emitter.

use nucleide_mcnp_io::sdef::{SdefCard, SdefDist, SdefProblem, SdefRef};
use nucleide_nuclei::particles::ParticleId;

use crate::parametric::{emission_histograms, ParametricPlasmaConfig};
use crate::{DriftReport, DriftRow, Error, PlasmaSourceConfig, Result, SourceModel, SpectrumSpec};

/// Half-width of the tabulated spectrum window, in sigma.
const TABLE_WIDTH_SIGMA: f64 = 4.0;

/// One emitted source card plus its drift report.
#[derive(Debug, Clone, PartialEq)]
pub struct EmittedCard {
    /// Canonical card text.
    pub text: String,
    /// Drift report for this emission.
    pub drift: DriftReport,
}

impl EmittedCard {
    /// Verify the card round-trips through the typed SDEF reader
    /// (`parse(text).emit() == text`); a cheap emission-time assertion.
    pub fn verify_round_trip(&self) -> Result<()> {
        let problem = nucleide_mcnp_io::sdef::parse_sdef_text(&self.text)
            .map_err(|e| Error::CardRoundTrip(e.to_string()))?;
        if problem.emit() != self.text {
            return Err(Error::CardRoundTrip(
                "re-emission is not byte-identical to the emitted card".to_string(),
            ));
        }
        Ok(())
    }
}

/// Neutron `PAR=` designator for the requested MCNP version.
fn neutron_designator(version: u32) -> Result<&'static str> {
    let designator = match version {
        5 => ParticleId::Neutron.mcnp(),
        6 => ParticleId::Neutron.mcnp6(),
        _ => return Err(Error::UnsupportedMcnpVersion(version)),
    };
    designator.ok_or(Error::NotYetSupported("MCNP neutron designator"))
}

/// Render a fusion source as an MCNP `SDEF` card plus drift report.
///
/// `version` selects the `PAR=` designator dialect (5 or 6); `n_bins` sets
/// the Gaussian tabulation bin count (ignored for monoenergetic spectra and
/// values `< 2` fall back to the 21-bin default).
pub fn emit_sdef(config: &PlasmaSourceConfig, version: u32, n_bins: usize) -> Result<EmittedCard> {
    config.validate()?;
    let spectrum = config.spectrum()?;
    let table = spectrum.tabulate(if n_bins >= 2 { n_bins } else { 21 }, TABLE_WIDTH_SIGMA)?;
    let par = neutron_designator(version)?;

    // Distribution numbering: the ring radial delta takes D1, so the
    // spectrum table starts at D2 for rings and D1 for points.
    let (spectrum_number, mut dists) = match config.model {
        SourceModel::Ring(r) => {
            let radial = SdefDist {
                number: 1,
                si: vec![r.radius_cm, r.radius_cm],
                sp: Some(vec![0.0, 1.0]),
                sb: None,
                line: 0,
            };
            (2, vec![radial])
        }
        SourceModel::Point(_) => (1, Vec::new()),
    };

    let erg = if spectrum.is_mono() {
        SdefRef::Literal(table.energies_mev[0])
    } else {
        dists.push(SdefDist {
            number: spectrum_number,
            si: table.energies_mev.clone(),
            sp: Some(table.probabilities.clone()),
            sb: None,
            line: 0,
        });
        SdefRef::Dist(spectrum_number)
    };

    let card = match config.model {
        SourceModel::Point(p) => SdefCard {
            pos: Some(SdefRef::Literal([p.x_cm, p.y_cm, p.z_cm])),
            erg: Some(erg),
            wgt: Some(SdefRef::Literal(config.weight)),
            par: Some(SdefRef::Literal(par.to_string())),
            ..SdefCard::default()
        },
        SourceModel::Ring(r) => SdefCard {
            pos: Some(SdefRef::Literal([0.0, 0.0, r.height_cm])),
            axs: Some(SdefRef::Literal([0.0, 0.0, 1.0])),
            rad: Some(SdefRef::Dist(1)),
            erg: Some(erg),
            wgt: Some(SdefRef::Literal(config.weight)),
            par: Some(SdefRef::Literal(par.to_string())),
            ..SdefCard::default()
        },
    };

    let text = SdefProblem { card, dists }.emit();
    let drift = drift_report(&spectrum, &table, true);
    Ok(EmittedCard { text, drift })
}

/// Shared report rows for one emission (SDEF: reparsed; Serpent: analytic).
pub(crate) fn drift_report(
    spectrum: &SpectrumSpec,
    table: &crate::SpectrumTable,
    reparsed: bool,
) -> DriftReport {
    let mut report = DriftReport::new();
    let note = match spectrum {
        SpectrumSpec::Mono { .. } => "monoenergetic line; no probability mass dropped".to_string(),
        SpectrumSpec::Gaussian { .. } => format!(
            "Gaussian tabulated at bin centers over +/-{:.0} sigma ({} lines); \
             tail mass beyond the table is dropped",
            table.width_sigma,
            table.energies_mev.len()
        ),
    };
    report.push(DriftRow::new(
        "emission probability",
        table.coverage,
        reparsed,
        note,
    ));
    report.push(DriftRow::new(
        "spatial distribution",
        1.0,
        reparsed,
        "all particles accounted at the requested position/ring",
    ));
    report
}

/// Render a parametric plasma source as an MCNP `SDEF` card plus drift
/// report.
///
/// The card carries the source's *marginals* as discrete histograms —
/// `RAD=D1` (birth minor-radius profile), `EXT=D2` (birth Z profile), and
/// `ERG=D3` (global marginal energy spectrum) — which is the strongest
/// product-form representation the discrete `SI`/`SP` card subset allows. A
/// partial toroidal sector adds `PHI=D4` (uniform toroidal-angle bins over
/// `[start, start + rotation)` radians); a full rotation (or no sector)
/// keeps the landed three-marginal card bit-for-bit. The drift report
/// quantifies the truncation and the joint-correlation information the card
/// cannot carry (see [`crate::parametric`]); the card still round-trips
/// byte-identically through the typed reader.
pub fn emit_sdef_parametric(
    config: &ParametricPlasmaConfig,
    version: u32,
    n_bins: usize,
) -> Result<EmittedCard> {
    config.validate()?;
    let bins = if n_bins >= 2 { n_bins } else { 21 };
    let hist = emission_histograms(config, bins)?;
    let par = neutron_designator(version)?;

    let dist = |number: u32, dist: &crate::BinnedDistribution| SdefDist {
        number,
        si: dist.centers.clone(),
        sp: Some(dist.masses.clone()),
        sb: None,
        line: 0,
    };
    let mut dists = vec![
        dist(1, &hist.radial),
        dist(2, &hist.vertical),
        dist(3, &hist.energy),
    ];
    let mut card = SdefCard {
        pos: Some(SdefRef::Literal([0.0, 0.0, 0.0])),
        axs: Some(SdefRef::Literal([0.0, 0.0, 1.0])),
        rad: Some(SdefRef::Dist(1)),
        ext: Some(SdefRef::Dist(2)),
        erg: Some(SdefRef::Dist(3)),
        wgt: Some(SdefRef::Literal(config.weight)),
        par: Some(SdefRef::Literal(par.to_string())),
        ..SdefCard::default()
    };
    if let Some(angle) = &hist.toroidal {
        dists.push(dist(4, angle));
        card.phi = Some(SdefRef::Dist(4));
    }
    let text = SdefProblem { card, dists }.emit();

    let mut report = DriftReport::new();
    report.push(DriftRow::new(
        "emission probability",
        hist.energy_coverage,
        true,
        format!(
            "marginal energy spectrum tabulated ({} lines, +/-4 sigma window); \
             tail mass is dropped and the local T_i correlation with birth \
             position is not representable on the card",
            hist.energy.centers.len(),
        ),
    ));
    report.push(DriftRow::new(
        "spatial marginals",
        1.0,
        true,
        "radial and vertical birth-profile marginals preserved as discrete histograms",
    ));
    report.push(DriftRow::new(
        "joint correlation",
        1.0 - hist.joint_correlation,
        true,
        "product-form card: half the L1 distance between the true (r, z) birth \
         joint and the product of its marginals is lost (0 = independent)",
    ));
    if let (Some(sector), Some(angle)) = (&config.sector, &hist.toroidal) {
        report.push(DriftRow::new(
            "toroidal sector",
            1.0,
            true,
            format!(
                "toroidal births uniform over [{:.6}, {:.6}) rad preserved as {} \
                 uniform angle bins (PHI=D4); emission totals scale by \
                 rotation/2π = {:.6}",
                sector.start_angle,
                sector.start_angle + sector.rotation_angle,
                angle.centers.len(),
                sector.fraction(),
            ),
        ));
    }
    Ok(EmittedCard {
        text,
        drift: report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FusionReaction;

    #[test]
    fn point_mono_card_is_the_e9_shape() {
        let config =
            PlasmaSourceConfig::point(1.0, 2.0, 3.0, FusionReaction::Dd, 0.0).with_weight(0.5);
        let card = emit_sdef(&config, 5, 21).unwrap();
        assert_eq!(
            card.text,
            "SDEF POS=1 2 3\
             \n     ERG=2.4495\
             \n     WGT=0.5\
             \n     PAR=n"
        );
        card.verify_round_trip().unwrap();
        assert_eq!(card.drift.worst_rel_drift(), 0.0);
    }

    #[test]
    fn point_gaussian_card_tabulates_spectrum() {
        let config = PlasmaSourceConfig::point(0.0, 0.0, 0.0, FusionReaction::Dt, 20.0);
        let card = emit_sdef(&config, 5, 21).unwrap();
        assert!(card
            .text
            .starts_with("SDEF POS=0 0 0\n     ERG=D1\n     WGT=1\n     PAR=n\nSI1 L "));
        assert!(card.text.contains("\nSP1 D "));
        // The reader sees the same 21-line table we tabulated, whatever the
        // 80-column wrapping did to the card text.
        let table = config.spectrum().unwrap().tabulate(21, 4.0).unwrap();
        let parsed = nucleide_mcnp_io::sdef::parse_sdef_text(&card.text).unwrap();
        assert_eq!(parsed.dists.len(), 1);
        assert_eq!(parsed.dists[0].si.len(), 21);
        // Card text is 6-significant-digit, so compare at card precision.
        for (a, b) in parsed.dists[0].si.iter().zip(table.energies_mev.iter()) {
            assert!(
                (a - b).abs() < 1e-5 * b.abs() + 1e-9,
                "SI energy {a} vs {b}"
            );
        }
        let sp = parsed.dists[0].sp.as_ref().unwrap();
        let sum: f64 = sp.iter().sum();
        assert!((sum - table.coverage).abs() < 1e-6);
        assert!((sum - table.coverage).abs() / table.coverage < 1e-5);
        for line in card.text.lines() {
            assert!(line.len() <= 80, "line over 80 columns: {line:?}");
        }
        card.verify_round_trip().unwrap();
        // Drift is the +/-4 sigma tail mass.
        let tail = 1.0 - table.coverage;
        assert!((card.drift.rows[0].rel_drift - tail).abs() < 1e-12);
        assert!(tail > 0.0 && tail < 1e-4);
    }

    #[test]
    fn ring_card_uses_radial_delta_and_round_trips() {
        let config = PlasmaSourceConfig::ring(300.0, 25.0, FusionReaction::Dt, 0.0);
        let card = emit_sdef(&config, 5, 21).unwrap();
        assert_eq!(
            card.text,
            "SDEF POS=0 0 25\
             \n     AXS=0 0 1\
             \n     RAD=D1\
             \n     ERG=14.021\
             \n     WGT=1\
             \n     PAR=n\
             \nSI1 L 300 300\
             \nSP1 D 0 1"
        );
        card.verify_round_trip().unwrap();
    }

    #[test]
    fn ring_gaussian_uses_distribution_numbers_without_collision() {
        let config = PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, 20.0);
        let card = emit_sdef(&config, 5, 21).unwrap();
        assert!(card.text.contains("RAD=D1\n"));
        assert!(card.text.contains("ERG=D2\n"));
        assert!(card.text.contains("\nSI1 L 300 300\nSP1 D 0 1\n"));
        assert!(card.text.contains("\nSI2 L "));
        assert!(card.text.contains("\nSP2 D "));
        card.verify_round_trip().unwrap();
    }

    #[test]
    fn version_six_designator_and_bad_version() {
        let config = PlasmaSourceConfig::point(0.0, 0.0, 0.0, FusionReaction::Dt, 0.0);
        assert!(emit_sdef(&config, 6, 21)
            .unwrap()
            .text
            .ends_with("\n     PAR=n"));
        assert_eq!(
            emit_sdef(&config, 4, 21),
            Err(Error::UnsupportedMcnpVersion(4))
        );
    }

    #[test]
    fn invalid_configs_are_loud() {
        let bad = PlasmaSourceConfig::ring(0.0, 0.0, FusionReaction::Dt, 10.0);
        assert!(emit_sdef(&bad, 5, 21).is_err());
    }

    #[test]
    fn parametric_card_carries_three_marginals_and_round_trips() {
        let config = crate::parametric::tests::iter_h_mode();
        let card = emit_sdef_parametric(&config, 5, 15).unwrap();
        let head = "SDEF POS=0 0 0\n     AXS=0 0 1\n     RAD=D1\n     EXT=D2\n     ERG=D3\n     WGT=1\n     PAR=n";
        assert!(card.text.starts_with(head), "card head:\n{}", card.text);
        assert!(card.text.contains("\nSI1 L "));
        assert!(card.text.contains("\nSP1 D "));
        assert!(card.text.contains("\nSI2 L "));
        assert!(card.text.contains("\nSP2 D "));
        assert!(card.text.contains("\nSI3 L "));
        assert!(card.text.contains("\nSP3 D "));
        for line in card.text.lines() {
            assert!(line.len() <= 80, "line over 80 columns: {line:?}");
        }
        card.verify_round_trip().unwrap();
        // Three drift rows: truncation, marginals, joint correlation.
        assert_eq!(card.drift.rows.len(), 3);
        assert!(card.drift.rows[2].rel_drift > 0.0);
        assert!(card.drift.rows[2].rel_drift < 1.0);
        // Re-parse and confirm the marginals survived at card precision.
        let parsed = nucleide_mcnp_io::sdef::parse_sdef_text(&card.text).unwrap();
        assert_eq!(parsed.dists.len(), 3);
        let sp1: f64 = parsed.dists[0].sp.as_ref().unwrap().iter().sum();
        assert!((sp1 - 1.0).abs() < 1e-4, "radial masses {sp1}");
    }

    #[test]
    fn parametric_bad_version_and_config_are_loud() {
        let config = crate::parametric::tests::iter_h_mode();
        assert_eq!(
            emit_sdef_parametric(&config, 4, 15),
            Err(Error::UnsupportedMcnpVersion(4))
        );
        let mut bad = config;
        bad.ion_temperature.centre_kev = -1.0;
        assert!(emit_sdef_parametric(&bad, 5, 15).is_err());
    }

    #[test]
    fn parametric_sector_card_carries_phi_marginal_and_round_trips() {
        use crate::parametric::ToroidalSector;
        use std::f64::consts::PI;
        let mut config = crate::parametric::tests::iter_h_mode();
        config.sector = Some(ToroidalSector::new(0.5, PI).unwrap());
        let card = emit_sdef_parametric(&config, 5, 15).unwrap();
        let head = "SDEF POS=0 0 0\n     AXS=0 0 1\n     RAD=D1\n     EXT=D2\n     PHI=D4\n     ERG=D3\n     WGT=1\n     PAR=n";
        assert!(card.text.starts_with(head), "card head:\n{}", card.text);
        assert!(card.text.contains("\nSI4 L "));
        assert!(card.text.contains("\nSP4 D "));
        for line in card.text.lines() {
            assert!(line.len() <= 80, "line over 80 columns: {line:?}");
        }
        card.verify_round_trip().unwrap();
        // Four drift rows: truncation, marginals, joint correlation, sector.
        assert_eq!(card.drift.rows.len(), 4);
        let sector = &card.drift.rows[3];
        assert_eq!(sector.quantity, "toroidal sector");
        assert_eq!(sector.accounted, 1.0);
        assert!(sector.reparsed);
        assert!(sector.note.contains("rotation/2π = 0.5"));
        // Re-parse: the angle marginal survived at card precision (15
        // uniform bins over [0.5, 0.5 + π)).
        let parsed = nucleide_mcnp_io::sdef::parse_sdef_text(&card.text).unwrap();
        assert_eq!(parsed.dists.len(), 4);
        assert_eq!(
            parsed.card.phi,
            Some(nucleide_mcnp_io::sdef::SdefRef::Dist(4))
        );
        let sp4: f64 = parsed.dists[3].sp.as_ref().unwrap().iter().sum();
        assert!((sp4 - 1.0).abs() < 1e-4, "angle masses {sp4}");
        assert_eq!(parsed.dists[3].si.len(), 15);
        assert!((parsed.dists[3].si[0] - (0.5 + PI / 15.0 / 2.0)).abs() < 1e-5);
    }

    #[test]
    fn parametric_full_rotation_card_matches_full_torus_bit_for_bit() {
        use crate::parametric::ToroidalSector;
        use std::f64::consts::PI;
        let full = crate::parametric::tests::iter_h_mode();
        let mut sector = full;
        sector.sector = Some(ToroidalSector::new(0.0, 2.0 * PI).unwrap());
        let a = emit_sdef_parametric(&full, 5, 15).unwrap();
        let b = emit_sdef_parametric(&sector, 5, 15).unwrap();
        assert_eq!(a.text, b.text);
        assert_eq!(a.drift, b.drift);
    }
}
