//! Sublet S3 gamma dose-rate kernel: slab and point dose over caller data.
//!
//! Pure arithmetic over caller-supplied gamma groups, pinned to two open
//! sources (never paywalled journal pages):
//!
//! - the FISPACT-II `DOSE` keyword documentation (the open
//!   `fispact.github.io` user documentation): the slab dose
//!   `D = C·B/2·Σᵢ μa(Eᵢ)/μm(Eᵢ)·Sγ(Eᵢ)` with buildup factor `B = 2` and the
//!   point dose `D = C·Σᵢ μa(Eᵢ)/(4πr²)·e^(−μ(Eᵢ)r)·mₛ·Sγ(Eᵢ)` with the
//!   gamma source rate `Sγ(Eᵢ) = Iᵢ·A(t)` summed over gamma energy groups
//!   (no nuclide split key), the mixture rule `μm = Σⱼ fⱼ·μmⱼ`, and the
//!   conversion constant `C = 3.6e9·|e|` (MeV/kg/s → Sv/h) with the exact SI
//!   elementary charge `|e| = 1.602176634e-19`;
//! - the CCFE-PR(16)53 preprint (Sublet, Eastwood, Morgan, Gilbert,
//!   Fleming, Arter — the author preprint of the FISPACT-II system paper):
//!   the dose-response context behind the same slab/point kernel.
//!
//! Every input is caller-supplied: the specific activity `A(t)` in Bq/kg,
//! the per-group intensities `Iᵢ`, the air and mixture attenuation
//! coefficients (directly per group, or folded from elemental values with
//! [`mixture_mu`]), the source mass `mₛ` in kg (point), and the distance `r`
//! in m (point). No attenuation table is vendored.
//!
//! Dose rates are Sv/h. Distances below [`MIN_DISTANCE_M`] (0.3 m) are
//! clamped up to it and reported LOUD via [`PointDose::clamped`] (never
//! silent). Explicitly out of scope: the >20% no-spectral-data warning and
//! the bremsstrahlung correction (both stay unimplemented by design).

use std::f64::consts::PI;

use crate::error::{Error, Result};

/// Exact SI elementary charge `|e|` in C (J/eV), pinning `C` below.
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// MeV/kg/s → Sv/h conversion constant `C = 3.6e9·|e|`.
pub const DOSE_CONVERSION_C: f64 = 3.6e9 * ELEMENTARY_CHARGE;

/// Slab buildup factor `B = 2` (pinned by the `DOSE` keyword documentation).
pub const SLAB_BUILDUP_B: f64 = 2.0;

/// Minimum point-source distance in m: smaller `r` clamps here, LOUD.
pub const MIN_DISTANCE_M: f64 = 0.3;

/// One gamma energy group: caller-supplied intensity and attenuation.
///
/// `intensity` is the group yield `Iᵢ` (photons per decay, finite, `>= 0`);
/// `mu_air` is `μa(Eᵢ)` and `mu` is the mixture `μm(Eᵢ)` (both finite,
/// `>= 0`, sharing one unit basis with the caller's distance convention:
/// `mu` in m⁻¹ when `r` is in m, so the `e^(−μr)` exponent is
/// dimensionless). The slab kernel additionally requires `mu > 0` (it divides
/// by the mixture coefficient); the point kernel allows `mu == 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoseGroup {
    /// Group yield `Iᵢ` in photons per decay (finite, `>= 0`).
    pub intensity: f64,
    /// Air attenuation coefficient `μa(Eᵢ)` (finite, `>= 0`).
    pub mu_air: f64,
    /// Mixture attenuation coefficient `μm(Eᵢ)` (finite, `>= 0`; `> 0` for
    /// the slab kernel).
    pub mu: f64,
}

/// S3 slab dose rate (Sv/h).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlabDose {
    /// `D = C·B/2·Σᵢ μa(Eᵢ)/μm(Eᵢ)·Sγ(Eᵢ)`, Sv/h.
    pub dose_sv_per_h: f64,
}

/// S3 point dose rate (Sv/h) with the distance-clamp record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointDose {
    /// `D = C·Σᵢ μa(Eᵢ)/(4πr²)·e^(−μ(Eᵢ)r)·mₛ·Sγ(Eᵢ)`, Sv/h.
    pub dose_sv_per_h: f64,
    /// Distance actually used in m (`max(r, 0.3)`).
    pub distance_used_m: f64,
    /// `true` when the caller `r` was below 0.3 m and got clamped (LOUD,
    /// never silent).
    pub clamped: bool,
}

/// Fold elemental mass attenuation coefficients into mixture values.
///
/// `μm(Eᵢ) = Σⱼ fⱼ·μmⱼ(Eᵢ)`: `fractions[j]` weights every group of
/// `element_mus[j]`. Fractions must be finite, `>= 0`, and sum to 1 within
/// `1e-9` (absolute); every coefficient must be finite and `>= 0`; every
/// element vector must hold the same group count. Violations are loud
/// [`Error`]s. Returns one mixture coefficient per group.
pub fn mixture_mu(fractions: &[f64], element_mus: &[Vec<f64>]) -> Result<Vec<f64>> {
    if fractions.is_empty() || element_mus.is_empty() {
        return Err(Error::BadAttenuation {
            group: usize::MAX,
            msg: "mixture needs at least one element fraction and one coefficient row".to_string(),
        });
    }
    if fractions.len() != element_mus.len() {
        return Err(Error::BadAttenuation {
            group: usize::MAX,
            msg: format!(
                "fraction count {} != element row count {}",
                fractions.len(),
                element_mus.len()
            ),
        });
    }
    let groups = element_mus[0].len();
    if groups == 0 {
        return Err(Error::BadAttenuation {
            group: usize::MAX,
            msg: "element coefficient rows must hold at least one group".to_string(),
        });
    }
    let mut frac_sum = 0.0f64;
    for (j, frac) in fractions.iter().enumerate() {
        if !frac.is_finite() || *frac < 0.0 {
            return Err(Error::BadAttenuation {
                group: usize::MAX,
                msg: format!("fraction[{j}] must be finite and >= 0, got {frac}"),
            });
        }
        frac_sum += frac;
    }
    if (frac_sum - 1.0).abs() > 1e-9 {
        return Err(Error::BadAttenuation {
            group: usize::MAX,
            msg: format!("fractions must sum to 1 within 1e-9, got {frac_sum}"),
        });
    }
    let mut out = vec![0.0f64; groups];
    for (j, row) in element_mus.iter().enumerate() {
        if row.len() != groups {
            return Err(Error::BadAttenuation {
                group: usize::MAX,
                msg: format!(
                    "element row {j} holds {} groups, expected {groups}",
                    row.len()
                ),
            });
        }
        for (i, mu) in row.iter().enumerate() {
            if !mu.is_finite() || *mu < 0.0 {
                return Err(Error::BadAttenuation {
                    group: i,
                    msg: format!("element {j} mu must be finite and >= 0, got {mu}"),
                });
            }
            out[i] += fractions[j] * mu;
        }
    }
    if !out.iter().all(|v| v.is_finite()) {
        return Err(Error::BadAttenuation {
            group: usize::MAX,
            msg: "mixture sum overflowed to non-finite".to_string(),
        });
    }
    Ok(out)
}

/// Check one shared group field pair, returning the gamma source rate.
fn check_group(
    group: &DoseGroup,
    index: usize,
    activity: f64,
    need_positive_mu: bool,
) -> Result<f64> {
    if !group.intensity.is_finite() || group.intensity < 0.0 {
        return Err(Error::BadDoseValue {
            field: format!("groups[{index}].intensity"),
            msg: format!("intensity must be finite and >= 0, got {}", group.intensity),
        });
    }
    if !group.mu_air.is_finite() || group.mu_air < 0.0 {
        return Err(Error::BadAttenuation {
            group: index,
            msg: format!("mu_air must be finite and >= 0, got {}", group.mu_air),
        });
    }
    if !group.mu.is_finite() || group.mu < 0.0 || (need_positive_mu && group.mu <= 0.0) {
        return Err(Error::BadAttenuation {
            group: index,
            msg: if need_positive_mu {
                format!(
                    "mu must be finite and > 0 for the slab ratio, got {}",
                    group.mu
                )
            } else {
                format!("mu must be finite and >= 0, got {}", group.mu)
            },
        });
    }
    Ok(group.intensity * activity)
}

/// Reject a bad specific activity or source mass.
fn check_scalar(value: f64, field: &str) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(Error::BadDoseValue {
            field: field.to_string(),
            msg: format!("{field} must be finite and >= 0, got {value}"),
        });
    }
    Ok(())
}

/// S3 slab dose `D = C·B/2·Σᵢ μa(Eᵢ)/μm(Eᵢ)·Sγ(Eᵢ)` (Sv/h).
///
/// `activity_bq_per_kg` is the specific activity `A(t)`; `Sγ(Eᵢ) =
/// Iᵢ·A(t)` per group. An empty group list doses `0.0`. Bad values and a
/// non-finite sum are loud [`Error`]s.
pub fn dose_slab(activity_bq_per_kg: f64, groups: &[DoseGroup]) -> Result<SlabDose> {
    check_scalar(activity_bq_per_kg, "activity_bq_per_kg")?;
    let mut sum = 0.0f64;
    for (index, group) in groups.iter().enumerate() {
        let source = check_group(group, index, activity_bq_per_kg, true)?;
        sum += (group.mu_air / group.mu) * source;
    }
    let dose = DOSE_CONVERSION_C * (SLAB_BUILDUP_B / 2.0) * sum;
    if !dose.is_finite() {
        return Err(Error::BadDoseValue {
            field: "total".to_string(),
            msg: format!(
                "slab dose sum overflowed to {dose} over {} groups",
                groups.len()
            ),
        });
    }
    Ok(SlabDose {
        dose_sv_per_h: dose,
    })
}

/// S3 point dose `D = C·Σᵢ μa(Eᵢ)/(4πr²)·e^(−μ(Eᵢ)r)·mₛ·Sγ(Eᵢ)` (Sv/h).
///
/// `source_mass_kg` is `mₛ`; `distance_m` is `r` and must be finite. A finite
/// `r` below [`MIN_DISTANCE_M`] clamps to it with [`PointDose::clamped`]
/// set (LOUD, never silent). An empty group list doses `0.0`. Bad values and
/// a non-finite sum are loud [`Error`]s.
pub fn dose_point(
    activity_bq_per_kg: f64,
    source_mass_kg: f64,
    distance_m: f64,
    groups: &[DoseGroup],
) -> Result<PointDose> {
    check_scalar(activity_bq_per_kg, "activity_bq_per_kg")?;
    check_scalar(source_mass_kg, "source_mass_kg")?;
    if !distance_m.is_finite() {
        return Err(Error::BadDoseValue {
            field: "distance_m".to_string(),
            msg: format!("distance_m must be finite, got {distance_m}"),
        });
    }
    let clamped = distance_m < MIN_DISTANCE_M;
    let r = if clamped { MIN_DISTANCE_M } else { distance_m };
    let denom = 4.0 * PI * r * r;
    let mut sum = 0.0f64;
    for (index, group) in groups.iter().enumerate() {
        let source = check_group(group, index, activity_bq_per_kg, false)?;
        sum += (group.mu_air / denom) * (-group.mu * r).exp() * source_mass_kg * source;
    }
    let dose = DOSE_CONVERSION_C * sum;
    if !dose.is_finite() {
        return Err(Error::BadDoseValue {
            field: "total".to_string(),
            msg: format!(
                "point dose sum overflowed to {dose} over {} groups",
                groups.len()
            ),
        });
    }
    Ok(PointDose {
        dose_sv_per_h: dose,
        distance_used_m: r,
        clamped,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(intensity: f64, mu_air: f64, mu: f64) -> DoseGroup {
        DoseGroup {
            intensity,
            mu_air,
            mu,
        }
    }

    #[test]
    fn constants_pin_the_pinned_definition() {
        assert_eq!(ELEMENTARY_CHARGE, 1.602_176_634e-19);
        assert_eq!(DOSE_CONVERSION_C, 3.6e9 * 1.602_176_634e-19);
        assert_eq!(SLAB_BUILDUP_B, 2.0);
        assert_eq!(MIN_DISTANCE_M, 0.3);
    }

    #[test]
    fn hand_computed_slab_vectors_at_exact_equality() {
        // A=2, I=0.5 -> S=1; ratio 1/2=0.5; B/2=1; D=C*0.5 exactly.
        let out = dose_slab(2.0, &[group(0.5, 1.0, 2.0)]).unwrap();
        assert_eq!(
            out.dose_sv_per_h,
            DOSE_CONVERSION_C * (SLAB_BUILDUP_B / 2.0) * ((1.0 / 2.0) * (0.5 * 2.0))
        );
        assert_eq!(out.dose_sv_per_h, DOSE_CONVERSION_C * 0.5);
        // Two groups: (A=4, I=0.25 -> S=1, ratio 2) + (A=4, I=0.5 -> S=2, ratio 0.5).
        let out = dose_slab(4.0, &[group(0.25, 2.0, 1.0), group(0.5, 1.0, 2.0)]).unwrap();
        let first = (2.0 / 1.0) * (0.25 * 4.0);
        let second = (1.0 / 2.0) * (0.5 * 4.0);
        assert_eq!(
            out.dose_sv_per_h,
            DOSE_CONVERSION_C * (SLAB_BUILDUP_B / 2.0) * (first + second)
        );
        assert_eq!(out.dose_sv_per_h, DOSE_CONVERSION_C * 3.0);
        // Empty group list doses zero.
        let out = dose_slab(4.0, &[]).unwrap();
        assert_eq!(out.dose_sv_per_h, 0.0);
        // Zero activity doses zero.
        let out = dose_slab(0.0, &[group(0.5, 1.0, 2.0)]).unwrap();
        assert_eq!(out.dose_sv_per_h, 0.0);
    }

    #[test]
    fn hand_computed_point_vectors_at_exact_equality() {
        // mu=0 kills the exponential (e^0=1): mirror the implementation order.
        let out = dose_point(4.0, 3.0, 1.0, &[group(0.5, 2.0, 0.0)]).unwrap();
        let denom = 4.0 * PI * 1.0 * 1.0;
        assert_eq!(
            out.dose_sv_per_h,
            DOSE_CONVERSION_C * ((2.0 / denom) * (-0.0f64 * 1.0).exp() * 3.0 * (0.5 * 4.0))
        );
        assert!(!out.clamped);
        assert_eq!(out.distance_used_m, 1.0);
        // Attenuated hand vector: mirror the implementation order exactly.
        let out = dose_point(2.0, 1.0, 2.0, &[group(1.0, 1.0, 0.5)]).unwrap();
        let denom = 4.0 * PI * 2.0 * 2.0;
        assert_eq!(
            out.dose_sv_per_h,
            DOSE_CONVERSION_C * ((1.0 / denom) * (-0.5f64 * 2.0).exp() * 1.0 * (1.0 * 2.0))
        );
        // Empty group list doses zero without clamping at r >= 0.3.
        let out = dose_point(2.0, 1.0, 1.0, &[]).unwrap();
        assert_eq!(out.dose_sv_per_h, 0.0);
        assert!(!out.clamped);
    }

    #[test]
    fn short_distance_clamps_loud_never_silent() {
        // r=0.1 clamps to 0.3 and says so; the dose equals the r=0.3 dose.
        let near = dose_point(2.0, 1.0, 0.1, &[group(1.0, 1.0, 0.5)]).unwrap();
        let at_floor = dose_point(2.0, 1.0, 0.3, &[group(1.0, 1.0, 0.5)]).unwrap();
        assert!(near.clamped);
        assert_eq!(near.distance_used_m, 0.3);
        assert!(!at_floor.clamped);
        assert_eq!(near.dose_sv_per_h, at_floor.dose_sv_per_h);
        // Zero and negative distances clamp too (still finite, still LOUD).
        for bad_r in [0.0, -1.0] {
            let out = dose_point(2.0, 1.0, bad_r, &[group(1.0, 1.0, 0.5)]).unwrap();
            assert!(out.clamped, "r={bad_r} should clamp");
            assert_eq!(out.distance_used_m, 0.3);
            assert_eq!(out.dose_sv_per_h, at_floor.dose_sv_per_h);
        }
    }

    #[test]
    fn mixture_folds_fractions_per_group() {
        // 0.25*Fe + 0.75*Al per group: [0.25*4+0.75*8, 0.25*2+0.75*6] = [7, 5].
        let out = mixture_mu(&[0.25, 0.75], &[vec![4.0, 2.0], vec![8.0, 6.0]]).unwrap();
        assert_eq!(out, vec![7.0, 5.0]);
        // Single element at fraction 1 reproduces its row.
        let out = mixture_mu(&[1.0], &[vec![3.0, 1.5]]).unwrap();
        assert_eq!(out, vec![3.0, 1.5]);
    }

    #[test]
    fn malformed_dose_inputs_are_loud() {
        // Bad specific activity / source mass.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    dose_slab(bad, &[group(0.5, 1.0, 2.0)]),
                    Err(Error::BadDoseValue { .. })
                ),
                "activity {bad} should fail"
            );
            assert!(
                matches!(
                    dose_point(bad, 1.0, 1.0, &[group(0.5, 1.0, 0.0)]),
                    Err(Error::BadDoseValue { .. })
                ),
                "activity {bad} should fail"
            );
            assert!(
                matches!(
                    dose_point(1.0, bad, 1.0, &[group(0.5, 1.0, 0.0)]),
                    Err(Error::BadDoseValue { .. })
                ),
                "mass {bad} should fail"
            );
        }
        // Non-finite distance (NaN/inf never clamps silently).
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                matches!(
                    dose_point(1.0, 1.0, bad, &[group(0.5, 1.0, 0.0)]),
                    Err(Error::BadDoseValue { .. })
                ),
                "distance {bad} should fail"
            );
        }
        // Bad group fields: intensity -> BadDoseValue; mus -> BadAttenuation.
        assert!(matches!(
            dose_slab(1.0, &[group(-0.5, 1.0, 2.0)]),
            Err(Error::BadDoseValue { .. })
        ));
        assert!(matches!(
            dose_slab(1.0, &[group(f64::NAN, 1.0, 2.0)]),
            Err(Error::BadDoseValue { .. })
        ));
        assert!(matches!(
            dose_slab(1.0, &[group(0.5, -1.0, 2.0)]),
            Err(Error::BadAttenuation { group: 0, .. })
        ));
        assert!(matches!(
            dose_slab(1.0, &[group(0.5, 1.0, f64::INFINITY)]),
            Err(Error::BadAttenuation { group: 0, .. })
        ));
        // Slab divides by mu: zero mu is loud (point allows it).
        assert!(matches!(
            dose_slab(1.0, &[group(0.5, 1.0, 0.0)]),
            Err(Error::BadAttenuation { group: 0, .. })
        ));
        dose_point(1.0, 1.0, 1.0, &[group(0.5, 1.0, 0.0)]).unwrap();
        // The bad group index is reported.
        assert!(matches!(
            dose_slab(1.0, &[group(0.5, 1.0, 2.0), group(0.5, 1.0, -3.0)]),
            Err(Error::BadAttenuation { group: 1, .. })
        ));
    }

    #[test]
    fn malformed_mixture_inputs_are_loud() {
        // Empty, ragged, or count-mismatched inputs.
        assert!(matches!(
            mixture_mu(&[], &[]),
            Err(Error::BadAttenuation { .. })
        ));
        assert!(matches!(
            mixture_mu(&[1.0], &[vec![1.0], vec![2.0]]),
            Err(Error::BadAttenuation { .. })
        ));
        assert!(matches!(
            mixture_mu(&[1.0], &[vec![]]),
            Err(Error::BadAttenuation { .. })
        ));
        assert!(matches!(
            mixture_mu(&[0.5, 0.5], &[vec![1.0, 2.0], vec![3.0]]),
            Err(Error::BadAttenuation { .. })
        ));
        // Fractions must be finite, >= 0, and sum to 1.
        assert!(matches!(
            mixture_mu(&[0.5, 0.25], &[vec![1.0], vec![2.0]]),
            Err(Error::BadAttenuation { .. })
        ));
        assert!(matches!(
            mixture_mu(&[-0.5, 1.5], &[vec![1.0], vec![2.0]]),
            Err(Error::BadAttenuation { .. })
        ));
        assert!(matches!(
            mixture_mu(&[f64::NAN], &[vec![1.0]]),
            Err(Error::BadAttenuation { .. })
        ));
        // Negative or non-finite coefficients.
        assert!(matches!(
            mixture_mu(&[1.0], &[vec![-1.0]]),
            Err(Error::BadAttenuation { group: 0, .. })
        ));
        assert!(matches!(
            mixture_mu(&[1.0], &[vec![f64::INFINITY]]),
            Err(Error::BadAttenuation { group: 0, .. })
        ));
    }

    #[test]
    fn overflowing_dose_sum_is_a_loud_error() {
        let out = dose_slab(f64::MAX, &[group(f64::MAX, f64::MAX, 1.0)]);
        match out {
            Err(Error::BadDoseValue { field, msg }) => {
                assert_eq!(field, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
        let out = dose_point(f64::MAX, f64::MAX, 1.0, &[group(1.0, f64::MAX, 0.0)]);
        assert!(
            matches!(out, Err(Error::BadDoseValue { .. })),
            "point overflow should fail, got {out:?}"
        );
    }
}
