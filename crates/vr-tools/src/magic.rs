//! MAGIC weight-window generation.
//!
//! Derives MCNP weight-window lower bounds from an MCNP mesh tally:
//!
//! ```text
//! max_val[g]  = max over all volume elements of flux[ve][g]
//! ww[ve][g]   = null_value                       if rel_error[ve][g] > tolerance
//!             = flux[ve][g] / (2 * max_val[g])   otherwise
//! ```
//!
//! Defaults: `tolerance = 0.5`, `null_value = 0.0`. The comparison against
//! the tolerance is strict (`>`), so a cell whose error equals the tolerance
//! keeps its scaled value.
//!
//! Following the legacy mesh-tag convention, the output also carries the
//! single maximum energy bound for total/single-bin tallies (or `e_bounds[1:]`
//! for multi-group tallies) in [`MagicOutput::e_upper_bounds`], plus the tag
//! names the data would be stored under (`ww_n`, `n_e_upper_bounds`, ...) in
//! [`MagicOutput::ww_tag_name`] / [`MagicOutput::e_upper_bounds_tag_name`].
//!
//! If every flux feeding one energy bin is non-positive, the normalization
//! would divide by zero; [`Error::ZeroMaxFlux`] is returned rather than
//! emitting `inf`/`nan`.

use mcnp_io::meshtal::{MeshTallyData, ParticleKind};

use crate::Error;

/// Which tally arrays feed the MAGIC algorithm.
///
/// Legacy tooling selected this implicitly through the mesh tag passed in
/// (`n_total_result` versus `n_result`); here it is explicit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MagicSelection {
    /// Energy-integrated totals: uses [`MeshTallyData::total_result`] /
    /// [`MeshTallyData::total_rel_error`]. One lower bound per volume element.
    Total,
    /// Per-energy-group values: uses [`MeshTallyData::result`] /
    /// [`MeshTallyData::rel_error`]. One lower bound per element per group,
    /// flattened `[ve][group]`.
    PerGroup,
}

/// MAGIC tuning parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagicParams {
    /// Maximum relative error for which a weight-window lower bound is
    /// generated. Elements above it receive `null_value`.
    pub tolerance: f64,
    /// Lower-bound value assigned where relative error exceeds `tolerance`.
    pub null_value: f64,
}

impl Default for MagicParams {
    fn default() -> Self {
        Self {
            tolerance: 0.5,
            null_value: 0.0,
        }
    }
}

/// Output of the MAGIC algorithm — the data written to mesh tags in legacy
/// workflows.
#[derive(Debug, Clone, PartialEq)]
pub struct MagicOutput {
    /// Weight-window lower bounds, row-major over volume elements:
    /// total mode → `lower_bounds_ww[ve]`; per-group mode →
    /// `lower_bounds_ww[ve * groups + g]`.
    pub lower_bounds_ww: Vec<f64>,
    /// Number of bounds per volume element (1 in total mode).
    pub groups_per_ve: usize,
    /// Per-energy-group normalization scale factors `1 / (2 * max_val[g])`
    /// that non-nulled fluxes were multiplied by (length `groups_per_ve`).
    pub scale_factors: Vec<f64>,
    /// Content of the `{particle}_e_upper_bounds` tag: the maximum
    /// energy bound in total mode, `e_bounds[1..]` otherwise.
    pub e_upper_bounds: Vec<f64>,
    /// Mesh tag name for the lower bounds
    /// (e.g. `"ww_n"` for neutron tallies).
    pub ww_tag_name: String,
    /// Mesh tag name for [`MagicOutput::e_upper_bounds`]
    /// (e.g. `"n_e_upper_bounds"`).
    pub e_upper_bounds_tag_name: String,
}

/// Run MAGIC on the energy-integrated totals of a tally with default
/// parameters (`tolerance = 0.5`, `null_value = 0.0`).
pub fn magic(tally: &MeshTallyData) -> Result<MagicOutput, Error> {
    magic_with(tally, MagicSelection::Total, MagicParams::default())
}

/// Run MAGIC on a tally with explicit array selection and parameters.
pub fn magic_with(
    tally: &MeshTallyData,
    selection: MagicSelection,
    params: MagicParams,
) -> Result<MagicOutput, Error> {
    if tally.num_ves() == 0 || tally.e_bounds.len() < 2 {
        return Err(Error::EmptyTally);
    }
    let groups_per_ve = match selection {
        MagicSelection::Total => 1,
        MagicSelection::PerGroup => tally.num_e_groups(),
    };

    // Flatten to ve-major [ve * groups_per_ve + g], matching the legacy
    // `vals[:]` / `errors[:]` tag reads.
    let (vals, errs): (Vec<f64>, Vec<f64>) = match selection {
        MagicSelection::Total => (tally.total_result.clone(), tally.total_rel_error.clone()),
        MagicSelection::PerGroup => {
            let mut v = Vec::with_capacity(tally.result.len() * groups_per_ve);
            let mut e = Vec::with_capacity(tally.rel_error.len() * groups_per_ve);
            for (rv, re) in tally.result.iter().zip(tally.rel_error.iter()) {
                v.extend_from_slice(rv);
                e.extend_from_slice(re);
            }
            (v, e)
        }
    };
    let expected = tally.num_ves() * groups_per_ve;
    if vals.len() != expected || errs.len() != expected {
        return Err(Error::LengthMismatch {
            expected,
            got: vals.len().max(errs.len()),
        });
    }
    for (i, &v) in vals.iter().enumerate() {
        if !v.is_finite() {
            return Err(Error::NonFiniteTally {
                field: "flux",
                index: i,
            });
        }
    }
    for (i, &e) in errs.iter().enumerate() {
        if !e.is_finite() {
            return Err(Error::NonFiniteTally {
                field: "error",
                index: i,
            });
        }
    }

    // max_val[i] = np.max over all ves for each energy bin.
    let mut max_val = vec![f64::NEG_INFINITY; groups_per_ve];
    for (idx, &v) in vals.iter().enumerate() {
        let g = idx % groups_per_ve;
        if v > max_val[g] {
            max_val[g] = v;
        }
    }
    for (g, &m) in max_val.iter().enumerate() {
        if m <= 0.0 {
            return Err(Error::ZeroMaxFlux { energy_group: g });
        }
    }

    // ww[ve][i] = null_value if error > tolerance else value / (2 * max_val[i]).
    let scale_factors: Vec<f64> = max_val.iter().map(|m| 1.0 / (2.0 * m)).collect();
    let mut lower_bounds_ww = Vec::with_capacity(vals.len());
    for (idx, (&v, &e)) in vals.iter().zip(errs.iter()).enumerate() {
        let g = idx % groups_per_ve;
        if e > params.tolerance {
            lower_bounds_ww.push(params.null_value);
        } else {
            lower_bounds_ww.push(v * scale_factors[g]);
        }
    }

    let letter = match tally.particle {
        ParticleKind::Neutron => 'n',
        ParticleKind::Photon => 'p',
    };
    let e_upper_bounds = match selection {
        MagicSelection::Total => vec![tally.e_bounds.last().copied().unwrap_or(0.0)],
        MagicSelection::PerGroup => tally.e_bounds[1..].to_vec(),
    };

    Ok(MagicOutput {
        lower_bounds_ww,
        groups_per_ve,
        scale_factors,
        e_upper_bounds,
        ww_tag_name: format!("ww_{letter}"),
        e_upper_bounds_tag_name: format!("{letter}_e_upper_bounds"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal hand-built neutron tally mirroring reference test meshes:
    /// 4 ves, e_bounds [0, 0.5, 1].
    fn sample_tally(
        result: Vec<Vec<f64>>,
        rel_error: Vec<Vec<f64>>,
        total_result: Vec<f64>,
        total_rel_error: Vec<f64>,
    ) -> MeshTallyData {
        MeshTallyData {
            tally_number: 4,
            particle: ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0, 1.0, 2.0],
            y_bounds: vec![-1.0, 3.0, 4.0],
            z_bounds: vec![10.0, 12.0],
            e_bounds: vec![0.0, 0.5, 1.0],
            column_idx: Default::default(),
            result,
            rel_error,
            total_result,
            total_rel_error,
        }
    }

    fn approx(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-12 * b.abs().max(1e-30),
            "expected {b}, got {a}"
        );
    }

    fn approx_vec(got: &[f64], want: &[f64]) {
        assert_eq!(got.len(), want.len());
        for (a, b) in got.iter().zip(want.iter()) {
            approx(*a, *b);
        }
    }

    #[test]
    fn oracle_total_below_default_tolerance() {
        // Replica of test_magic_below_tolerance: all errors < 0.5, defaults.
        let t = sample_tally(
            vec![vec![0.0]; 4],
            vec![vec![0.0]; 4],
            vec![1.2, 3.3, 1.6, 1.7],
            vec![0.11, 0.013, 0.14, 0.19],
        );
        let out = magic(&t).unwrap();
        assert_eq!(out.groups_per_ve, 1);
        assert_eq!(out.ww_tag_name, "ww_n");
        assert_eq!(out.e_upper_bounds_tag_name, "n_e_upper_bounds");
        assert_eq!(out.scale_factors, vec![1.0 / (2.0 * 3.3)]);
        approx_vec(
            &out.lower_bounds_ww,
            &[1.2 / 6.6, 0.5, 1.6 / 6.6, 1.7 / 6.6],
        );
    }

    #[test]
    fn oracle_total_nulling_custom_params() {
        // Replica of test_magic_e_total: error 0.19 > 0.15 nulls to 0.001.
        let t = sample_tally(
            vec![vec![0.0]; 4],
            vec![vec![0.0]; 4],
            vec![1.2, 3.3, 1.6, 1.7],
            vec![0.11, 0.013, 0.14, 0.19],
        );
        let out = magic_with(
            &t,
            MagicSelection::Total,
            MagicParams {
                tolerance: 0.15,
                null_value: 0.001,
            },
        )
        .unwrap();
        approx_vec(&out.lower_bounds_ww, &[1.2 / 6.6, 0.5, 1.6 / 6.6, 0.001]);
    }

    #[test]
    fn oracle_multi_bin_strict_tolerance_compare() {
        // Replica of test_magic_multi_bins. Group-1 error 0.16 > 0.15 nulls;
        // the boundary case itself is covered by the fixture tests below via
        // strict inequality in the implementation.
        let t = sample_tally(
            vec![
                vec![1.2, 3.3],
                vec![1.6, 1.7],
                vec![1.5, 1.4],
                vec![2.6, 1.0],
            ],
            vec![
                vec![0.11, 0.013],
                vec![0.14, 0.19],
                vec![0.02, 0.16],
                vec![0.04, 0.09],
            ],
            vec![0.0; 4],
            vec![0.0; 4],
        );
        let out = magic_with(
            &t,
            MagicSelection::PerGroup,
            MagicParams {
                tolerance: 0.15,
                null_value: 0.001,
            },
        )
        .unwrap();
        assert_eq!(out.groups_per_ve, 2);
        assert_eq!(out.e_upper_bounds, vec![0.5, 1.0]);
        let want = [
            [1.2 / 5.2, 0.5],
            [1.6 / 5.2, 0.001],
            [1.5 / 5.2, 0.001],
            [0.5, 1.0 / 6.6],
        ];
        for (ve, row) in want.iter().enumerate() {
            approx_vec(&out.lower_bounds_ww[ve * 2..ve * 2 + 2], row);
        }
    }

    #[test]
    fn fixture_total_mode_matches_hand_computed_magic() {
        let m = mcnp_io::meshtal::Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt"))
            .unwrap();
        let t = &m.tallies[&4];
        let out = magic(t).unwrap();
        assert_eq!(out.groups_per_ve, 1);
        assert_eq!(out.ww_tag_name, "ww_n");
        // max(total_result) = 1.31488E-06, no rel error exceeds 0.5.
        assert_eq!(out.e_upper_bounds, vec![1.0]);
        assert_eq!(out.scale_factors, vec![1.0 / (2.0 * 1.31488e-06)]);
        let want_first = [
            0.07277089924555853,
            0.09303472560233633,
            0.059732447067413,
            0.1160508943781942,
            0.1535896811876369,
        ];
        approx_vec(&out.lower_bounds_ww[..5], &want_first);
        let n = t.num_ves();
        let want_last = [
            0.07414554940374787,
            0.09211715137503042,
            0.06157482051594062,
        ];
        approx_vec(&out.lower_bounds_ww[n - 3..], &want_last);
    }

    #[test]
    fn fixture_group_mode_matches_hand_computed_magic() {
        let m = mcnp_io::meshtal::Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt"))
            .unwrap();
        let t = &m.tallies[&4];
        let out = magic_with(t, MagicSelection::PerGroup, MagicParams::default()).unwrap();
        assert_eq!(out.groups_per_ve, 3);
        assert_eq!(out.e_upper_bounds, vec![0.1, 0.2, 1.0]);
        // Per-group flux maxima across all 45 ves (hand-computed from the
        // fixture): 4.45445E-08, 1.04704E-07, 1.16563E-06.
        assert_eq!(out.scale_factors.len(), 3);
        approx(out.scale_factors[0], 1.0 / (2.0 * 4.45445e-08));
        approx(out.scale_factors[1], 1.0 / (2.0 * 1.04704e-07));
        approx(out.scale_factors[2], 1.0 / (2.0 * 1.16563e-06));
        // Cell 22 holds the maximum of every group simultaneously.
        approx_vec(&out.lower_bounds_ww[22 * 3..22 * 3 + 3], &[0.5, 0.5, 0.5]);
        approx_vec(
            &out.lower_bounds_ww[..3],
            &[
                0.05572753089607023,
                0.08699954156479218,
                0.07214424817480675,
            ],
        );
        let last = t.num_ves() - 1;
        approx_vec(
            &out.lower_bounds_ww[last * 3..],
            &[
                0.04513037524273479,
                0.06796731738997555,
                0.06162933349347562,
            ],
        );
    }

    #[test]
    fn fixture_tight_tolerance_produces_null_values() {
        let m = mcnp_io::meshtal::Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt"))
            .unwrap();
        let t = &m.tallies[&4];
        let out = magic_with(
            t,
            MagicSelection::PerGroup,
            MagicParams {
                tolerance: 0.06,
                null_value: 1e-3,
            },
        )
        .unwrap();
        // Hand-computed: 87 of the 135 (ve, group) pairs exceed tol = 0.06.
        let nulled = out.lower_bounds_ww.iter().filter(|&&w| w == 1e-3).count();
        assert_eq!(nulled, 87);
        approx_vec(
            &out.lower_bounds_ww[..3],
            &[1e-3, 1e-3, 0.07214424817480675],
        );
    }

    #[test]
    fn zero_max_flux_group_is_rejected() {
        let t = sample_tally(
            vec![vec![0.0, 1.0]; 4],
            vec![vec![0.0, 0.0]; 4],
            vec![0.0; 4],
            vec![0.0; 4],
        );
        assert_eq!(
            magic_with(&t, MagicSelection::PerGroup, MagicParams::default()),
            Err(Error::ZeroMaxFlux { energy_group: 0 })
        );
        // Totals all zero too.
        assert_eq!(magic(&t), Err(Error::ZeroMaxFlux { energy_group: 0 }));
    }

    #[test]
    fn empty_tally_rejected() {
        let t = MeshTallyData {
            tally_number: 1,
            particle: ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0],
            y_bounds: vec![0.0],
            z_bounds: vec![0.0],
            e_bounds: vec![],
            column_idx: Default::default(),
            result: Vec::new(),
            rel_error: Vec::new(),
            total_result: Vec::new(),
            total_rel_error: Vec::new(),
        };
        assert_eq!(magic(&t), Err(Error::EmptyTally));
    }

    #[test]
    fn length_mismatch_detected() {
        let mut t = sample_tally(
            vec![vec![1.0], vec![1.0]],
            vec![vec![0.0], vec![0.0]],
            vec![1.0],
            vec![0.0],
        );
        t.total_rel_error = vec![0.0, 0.0];
        assert!(matches!(magic(&t), Err(Error::LengthMismatch { .. })));
    }

    #[test]
    fn empty_e_bounds_with_volume_elements_does_not_panic() {
        // num_ves() > 0 but e_bounds is empty used to underflow in
        // num_e_groups(). The guard must check e_bounds.len() < 2 directly.
        let t = MeshTallyData {
            tally_number: 1,
            particle: ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0, 1.0],
            y_bounds: vec![0.0, 1.0],
            z_bounds: vec![0.0, 1.0],
            e_bounds: vec![],
            column_idx: Default::default(),
            result: vec![vec![1.0]],
            rel_error: vec![vec![0.0]],
            total_result: vec![1.0],
            total_rel_error: vec![0.0],
        };
        assert_eq!(magic(&t), Err(Error::EmptyTally));
    }

    #[test]
    fn non_finite_flux_rejected() {
        let mut t = sample_tally(
            vec![vec![0.0]; 4],
            vec![vec![0.0]; 4],
            vec![1.0, f64::NAN, 1.0, 1.0],
            vec![0.0; 4],
        );
        assert!(matches!(
            magic(&t),
            Err(Error::NonFiniteTally {
                field: "flux",
                index: 1
            })
        ));

        t.total_result[1] = f64::INFINITY;
        assert!(matches!(
            magic(&t),
            Err(Error::NonFiniteTally {
                field: "flux",
                index: 1
            })
        ));
    }

    #[test]
    fn non_finite_error_rejected() {
        let t = sample_tally(
            vec![vec![0.0]; 4],
            vec![vec![0.0]; 4],
            vec![1.0; 4],
            vec![0.0, f64::NAN, 0.0, 0.0],
        );
        assert!(matches!(
            magic(&t),
            Err(Error::NonFiniteTally {
                field: "error",
                index: 1
            })
        ));
    }

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/meshtal/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }
}
