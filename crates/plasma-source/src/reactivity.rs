//! Thermonuclear reactivity fits ⟨σv⟩(Tᵢ) — Bosch & Hale, Nucl. Fusion
//! **32** (1992) 611, in the compact parametrization reproduced by
//! Atzeni & Meyer-ter-Vehn, *The Physics of Inertial Fusion* (and the form
//! the MIT `openmc-plasma-source`/`NeSST` oracles implement):
//!
//! ```text
//! ⟨σv⟩ = C1 · η^(−5/6) · ξ² · exp(−3·η^(1/3)·ξ)        [m³/s], T in keV
//! ξ    = A · T^(−1/3)
//! η    = 1 − num(T)/den(T)        (rational correction, per-reaction coeffs)
//! ```
//!
//! Coefficients are the published fits (plain published facts). Reactivity
//! vanishes as `T → 0` (the exponential beats the power-law growth), which
//! the sampler relies on at cold separatrix temperatures.
//!
//! Gates pinned: hand recomputation at 1/10/20 keV against an independent
//! float64 evaluation of the published formula, the `T → 0` limit, and a
//! container-oracle leg against NeSST's `reac_DT`/`reac_DD` (agreement to
//! ~5e-9 — the published coefficients' own precision).

use crate::{Error, FusionReaction, Result};

/// Published Bosch–Hale/Atzeni coefficients for one reaction.
#[derive(Debug, Clone, Copy)]
struct ReactivityFit {
    /// Overall scale `C1` \[m³/s\].
    c1: f64,
    /// `ξ = xi_coeff · T^(−1/3)` coefficient `A`.
    xi_coeff: f64,
    /// `η` numerator polynomial `n3·T³ + n2·T² + n1·T + n0`.
    num: [f64; 4],
    /// `η` denominator polynomial `d3·T³ + d2·T² + d1·T + d0`.
    den: [f64; 4],
}

/// T(d,n)⁴He fit (Atzeni & Meyer-ter-Vehn Table; = Bosch–Hale DT).
const DT_REACTIVITY: ReactivityFit = ReactivityFit {
    c1: 643.41e-22,
    xi_coeff: 6.6610,
    num: [-1.0675e-4, 4.6064e-3, 1.5136e-2, 0.0],
    den: [1.366e-5, 1.35e-2, 7.5189e-2, 1.0],
};

/// D-D fit (the neutron-relevant D(d,n)³He branch shares the published D-D
/// fit to the precision of the coefficients; see the module docs).
const DD_REACTIVITY: ReactivityFit = ReactivityFit {
    c1: 3.5741e-22,
    xi_coeff: 6.2696,
    num: [0.0, 0.0, 5.8577e-3, 0.0],
    den: [0.0, -2.964e-6, 7.6822e-3, 1.0],
};

fn eval_poly(c: &[f64; 4], t: f64) -> f64 {
    c[3] + t * (c[2] + t * (c[1] + t * c[0]))
}

impl FusionReaction {
    /// The published reactivity fit for this reaction.
    fn reactivity_fit(self) -> &'static ReactivityFit {
        match self {
            FusionReaction::Dt => &DT_REACTIVITY,
            FusionReaction::Dd => &DD_REACTIVITY,
        }
    }

    /// Thermonuclear reactivity ⟨σv⟩ \[m³/s\] at ion temperature `ti_kev`
    /// \[keV\] (Bosch & Hale 1992; Atzeni–Meyer-ter-Vehn parametrization).
    ///
    /// Returns 0 at `T = 0` and errors on negative/non-finite temperature.
    /// Outside the fit's validity domain the result is a loud
    /// [`Error::FitOutOfDomain`], never `NaN` (the D-D η factor goes
    /// non-positive on roughly 965–2720 keV, far above the published range).
    pub fn reactivity_m3_per_s(self, ti_kev: f64) -> Result<f64> {
        if !ti_kev.is_finite() {
            return Err(Error::NonFinite("ion temperature"));
        }
        if ti_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(ti_kev));
        }
        if ti_kev == 0.0 {
            return Ok(0.0);
        }
        let fit = self.reactivity_fit();
        let xi = fit.xi_coeff * ti_kev.powf(-1.0 / 3.0);
        let eta = 1.0 - eval_poly(&fit.num, ti_kev) / eval_poly(&fit.den, ti_kev);
        if eta <= 0.0 {
            return Err(Error::FitOutOfDomain {
                reaction: self.label(),
                ti_kev,
            });
        }
        let value =
            fit.c1 * eta.powf(-5.0 / 6.0) * xi * xi * (-3.0 * eta.powf(1.0 / 3.0) * xi).exp();
        if !value.is_finite() {
            return Err(Error::FitOutOfDomain {
                reaction: self.label(),
                ti_kev,
            });
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden values recomputed independently from the published
    /// coefficients (float64; documented in
    /// `validation/plasma_source_vs_openmc.py` gate P6). They pin the
    /// coefficient transcription and the unit convention (keV in, m³/s out).
    #[test]
    fn reactivity_matches_independent_recomputation() {
        let gold: &[(FusionReaction, f64, f64)] = &[
            (FusionReaction::Dt, 1.0, 6.857_069_187_960_237e-27),
            (FusionReaction::Dt, 10.0, 1.136_173_488_699_605_8e-22),
            (FusionReaction::Dt, 20.0, 4.330_220_397_689_175e-22),
            (FusionReaction::Dd, 1.0, 9.932_467_364_810_979e-29),
            (FusionReaction::Dd, 10.0, 6.022_478_475_561_915e-25),
            (FusionReaction::Dd, 20.0, 2.602_582_958_721_524e-24),
        ];
        for &(reaction, ti, want) in gold {
            let got = reaction.reactivity_m3_per_s(ti).unwrap();
            assert!(
                (got - want).abs() < 1e-15 * want,
                "{reaction:?} T={ti}: {got} vs {want}"
            );
        }
    }

    #[test]
    fn reactivity_vanishes_at_zero_and_rises_through_the_operating_range() {
        for reaction in [FusionReaction::Dt, FusionReaction::Dd] {
            assert_eq!(reaction.reactivity_m3_per_s(0.0).unwrap(), 0.0);
            let mut prev = 0.0;
            for ti in [1.0, 5.0, 10.0, 20.0, 40.0] {
                let v = reaction.reactivity_m3_per_s(ti).unwrap();
                assert!(v > prev, "{reaction:?} rising to 40 keV");
                prev = v;
            }
        }
        // DT reactivity peaks near 64 keV and turns over by 80 keV (the
        // Bosch–Hale fit shape), another independent transcription pin.
        let at64 = FusionReaction::Dt.reactivity_m3_per_s(64.0).unwrap();
        let at80 = FusionReaction::Dt.reactivity_m3_per_s(80.0).unwrap();
        assert!((at64 - 8.940_872_172_398_506e-22).abs() < 1e-15 * at64);
        assert!(at80 < at64);
        assert_eq!(
            FusionReaction::Dt.reactivity_m3_per_s(-1.0),
            Err(Error::NegativeIonTemperature(-1.0))
        );
    }

    #[test]
    fn dd_reactivity_is_loud_outside_the_fit_domain() {
        // The D-D η factor goes non-positive on roughly 965–2720 keV (far above
        // the published fit range); the old code returned Ok(NaN) there.
        match FusionReaction::Dd.reactivity_m3_per_s(1000.0) {
            Err(Error::FitOutOfDomain { reaction, ti_kev }) => {
                assert_eq!(reaction, "D-D");
                assert_eq!(ti_kev, 1000.0);
            }
            other => panic!("expected FitOutOfDomain, got {other:?}"),
        }
        // In-domain values stay finite and positive.
        assert!(FusionReaction::Dd.reactivity_m3_per_s(100.0).unwrap() > 0.0);
    }
}
