//! Calibration: quadratic energy bins (E6), efficiency evaluation (E7), and
//! the E7-fit log-space efficiency-coefficient fit.
//!
//! FWHM coefficients are parsed by the `.spe` readers but never evaluated
//! here — the upstream module carries them without any evaluation routine.

use crate::Error;

/// Energy per channel (E6): `ebin[ch] = a0 + a1*ch + a2*ch^2`.
pub fn energy_bins(channels: &[f64], fit: &[f64]) -> Result<Vec<f64>, Error> {
    if fit.len() < 3 {
        return Err(Error::MissingCalibration(fit.len()));
    }
    Ok(channels
        .iter()
        .map(|ch| fit[0] + fit[1] * ch + fit[2] * ch * ch)
        .collect())
}

/// Detector efficiency at `energy_mev` (E7, energy in MeV).
///
/// `fit == 1`: `eff = exp(sum a_i * (ln E)^i)`;
/// `fit == 2`: `eff = exp(sum a_i * (1/E)^i)`;
/// anything else is an error. Coefficients are caller-supplied or come out
/// of [`fit_efficiency`]: the upstream module has no coefficient-fitting
/// routine (so the fit has no upstream parity), but this crate fits the E7
/// log-polynomial form to caller points through the workspace least-squares
/// kernel.
pub fn detector_efficiency(energy_mev: f64, coeff: &[f64], fit: i64) -> Result<f64, Error> {
    if coeff.is_empty() {
        return Err(Error::EmptyCoefficients);
    }
    let log_eff = match fit {
        1 => {
            let l = energy_mev.ln();
            let mut acc = coeff[0];
            let mut power = 1.0;
            for &a in &coeff[1..] {
                power *= l;
                acc += a * power;
            }
            acc
        }
        2 => {
            let inv = 1.0 / energy_mev;
            let mut acc = coeff[0];
            let mut power = 1.0;
            for &a in &coeff[1..] {
                power *= inv;
                acc += a * power;
            }
            acc
        }
        other => return Err(Error::UnknownEfficiencyFit(other)),
    };
    Ok(log_eff.exp())
}

/// Efficiency-coefficient fit (E7-fit, energies in MeV).
///
/// Log-space linear least squares over caller points: targets `y_i = ln eff_i`
/// with design rows `X[i][j] = (ln E_i)^j` (`fit == 1`) or `(1/E_i)^j`
/// (`fit == 2`), minimizing `sum w_i (X_i β − y_i)^2` through the workspace
/// least-squares kernel (`nucleide-linalg` `lstsq`, faer column-pivoted QR).
/// Weights are caller-supplied (a zero weight skips its row; all-zero is an
/// error). Returns `order + 1` coefficients feeding [`detector_efficiency`].
/// This is the one coefficient-fitting path in the crate — every other
/// coefficient vector stays a caller input.
pub fn fit_efficiency(
    energies_mev: &[f64],
    effs: &[f64],
    weights: &[f64],
    order: usize,
    fit: i64,
) -> Result<Vec<f64>, Error> {
    if energies_mev.is_empty() {
        return Err(Error::EmptyFitInput);
    }
    if energies_mev.len() != effs.len() || energies_mev.len() != weights.len() {
        return Err(Error::MismatchedFitLengths {
            energies: energies_mev.len(),
            effs: effs.len(),
            weights: weights.len(),
        });
    }
    if fit != 1 && fit != 2 {
        return Err(Error::UnknownEfficiencyFit(fit));
    }
    let ncoeff = order + 1;
    if energies_mev.len() < ncoeff {
        return Err(Error::UnderdeterminedFit {
            points: energies_mev.len(),
            coeffs: ncoeff,
        });
    }
    for &e in energies_mev {
        if !e.is_finite() {
            return Err(Error::NonFiniteFitValue("energies"));
        }
        if e <= 0.0 {
            return Err(Error::NonPositiveEnergy(e));
        }
    }
    for &v in effs {
        if !v.is_finite() {
            return Err(Error::NonFiniteFitValue("effs"));
        }
        if v <= 0.0 {
            return Err(Error::NonPositiveEfficiency(v));
        }
    }
    for &w in weights {
        if !w.is_finite() {
            return Err(Error::NonFiniteFitValue("weights"));
        }
        if w < 0.0 {
            return Err(Error::NegativeWeight(w));
        }
    }
    let mut design = Vec::with_capacity(energies_mev.len());
    for &e in energies_mev {
        let basis = if fit == 1 { e.ln() } else { 1.0 / e };
        let mut row = Vec::with_capacity(ncoeff);
        let mut power = 1.0;
        for _ in 0..ncoeff {
            row.push(power);
            power *= basis;
        }
        design.push(row);
    }
    let target: Vec<f64> = effs.iter().map(|v| v.ln()).collect();
    nucleide_linalg::lstsq::weighted_lstsq(&design, &target, weights)
        .map_err(|e| Error::LstsqFit(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-computed: fit [1.5, 2.0, 0.5] on channels [0,1,2].
    #[test]
    fn energy_bins_quadratic() {
        let out = energy_bins(&[0.0, 1.0, 2.0], &[1.5, 2.0, 0.5]).unwrap();
        assert_eq!(out, vec![1.5, 4.0, 7.5]);
    }

    // Upstream golden vector (test oracle): at E = 1 MeV, ln E = 0, so
    // eff = exp(a0) = exp(-2.81861504261204) = 0.059688551591347033.
    // (Oracle literals trimmed of trailing zeros; identical f64 values.)
    #[test]
    fn efficiency_golden_fit1_at_1mev() {
        let coeff = [
            -2.81861504261204,
            -0.727352820018942,
            -0.0395798886481904,
            -0.0592305254664096,
            0.023772637347443,
            0.0325306475072671,
        ];
        let eff = detector_efficiency(1.0, &coeff, 1).unwrap();
        let rel = (eff - 0.059_688_551_591_347_03).abs() / 0.059_688_551_591_347_03;
        assert!(rel < 1e-12, "rel err {rel:.3e}");
    }

    // Hand-computed fit2: coeff [0.5, -1.0], E = 2 MeV:
    // log = 0.5 - 1/2 = 0, eff = 1.
    #[test]
    fn efficiency_fit2_hand_value() {
        let eff = detector_efficiency(2.0, &[0.5, -1.0], 2).unwrap();
        assert!((eff - 1.0).abs() < 1e-15);
    }

    #[test]
    fn efficiency_rejects() {
        let coeff = [0.5];
        assert_eq!(
            detector_efficiency(1.0, &coeff, 10),
            Err(Error::UnknownEfficiencyFit(10))
        );
        assert_eq!(
            detector_efficiency(1.0, &[], 1),
            Err(Error::EmptyCoefficients)
        );
        assert_eq!(
            energy_bins(&[0.0], &[1.0]),
            Err(Error::MissingCalibration(1))
        );
    }

    // E7-fit: points exactly on a degree-2 ln-basis curve recover the
    // coefficients (closed-form synthetic gate, tolerance 1e-9).
    #[test]
    fn fit_efficiency_fit1_recovers_closed_form() {
        let coeff = [-2.0, -0.8, -0.05];
        let energies = [0.3, 0.5, 0.662, 1.0, 1.33, 2.0];
        let effs: Vec<f64> = energies
            .iter()
            .map(|e| detector_efficiency(*e, &coeff, 1).unwrap())
            .collect();
        let weights = vec![1.0, 1.0, 2.0, 1.0, 1.0, 1.0];
        let got = fit_efficiency(&energies, &effs, &weights, 2, 1).unwrap();
        for (g, w) in got.iter().zip(coeff.iter()) {
            assert!((g - w).abs() < 1e-9, "got {got:?} want {coeff:?}");
        }
        // Round-trip: the fitted coefficients re-evaluate to the inputs.
        for (e, v) in energies.iter().zip(effs.iter()) {
            let back = detector_efficiency(*e, &got, 1).unwrap();
            let rel = (back - v).abs() / v;
            assert!(rel < 1e-9, "E={e} rel err {rel:.3e}");
        }
    }

    // E7-fit: points exactly on a degree-1 inverse-E curve recover the
    // coefficients (closed-form synthetic gate, tolerance 1e-9).
    #[test]
    fn fit_efficiency_fit2_recovers_closed_form() {
        let coeff = [0.5, -1.0];
        let energies = [0.5, 1.0, 2.0, 3.0];
        let effs: Vec<f64> = energies
            .iter()
            .map(|e| detector_efficiency(*e, &coeff, 2).unwrap())
            .collect();
        let weights = vec![1.0; 4];
        let got = fit_efficiency(&energies, &effs, &weights, 1, 2).unwrap();
        for (g, w) in got.iter().zip(coeff.iter()) {
            assert!((g - w).abs() < 1e-9, "got {got:?} want {coeff:?}");
        }
    }

    // E7-fit: a zero-weight outlier is skipped (weights are caller-supplied).
    #[test]
    fn fit_efficiency_zero_weight_skips_row() {
        // y = ln eff on the fit-1 line a0 + a1 ln E with a = [0.0, 1.0]:
        // eff = E, so (1,1) and (2,2) pin the line; (3,300) is skipped.
        let got =
            fit_efficiency(&[1.0, 2.0, 3.0], &[1.0, 2.0, 300.0], &[1.0, 1.0, 0.0], 1, 1).unwrap();
        assert!(got[0].abs() < 1e-9, "a0 = {}", got[0]);
        assert!((got[1] - 1.0).abs() < 1e-9, "a1 = {}", got[1]);
    }

    #[test]
    fn fit_efficiency_rejects() {
        let energies = [0.5, 1.0, 2.0];
        let effs = [0.1, 0.2, 0.3];
        let weights = [1.0, 1.0, 1.0];
        assert_eq!(
            fit_efficiency(&[], &[], &[], 1, 1),
            Err(Error::EmptyFitInput)
        );
        assert_eq!(
            fit_efficiency(&energies, &effs[..2], &weights, 1, 1),
            Err(Error::MismatchedFitLengths {
                energies: 3,
                effs: 2,
                weights: 3
            })
        );
        assert_eq!(
            fit_efficiency(&energies, &effs, &weights, 1, 7),
            Err(Error::UnknownEfficiencyFit(7))
        );
        assert_eq!(
            fit_efficiency(&energies, &effs, &weights, 5, 1),
            Err(Error::UnderdeterminedFit {
                points: 3,
                coeffs: 6
            })
        );
        assert_eq!(
            fit_efficiency(&[0.0, 1.0], &[0.1, 0.2], &[1.0, 1.0], 1, 1),
            Err(Error::NonPositiveEnergy(0.0))
        );
        assert_eq!(
            fit_efficiency(&[0.5, 1.0], &[0.0, 0.2], &[1.0, 1.0], 1, 1),
            Err(Error::NonPositiveEfficiency(0.0))
        );
        assert_eq!(
            fit_efficiency(&[f64::NAN, 1.0], &[0.1, 0.2], &[1.0, 1.0], 1, 1),
            Err(Error::NonFiniteFitValue("energies"))
        );
        assert_eq!(
            fit_efficiency(&[0.5, 1.0], &[0.1, 0.2], &[1.0, -1.0], 1, 1),
            Err(Error::NegativeWeight(-1.0))
        );
        // All-zero weights reach the shared kernel, which names the cause.
        assert!(matches!(
            fit_efficiency(&[0.5, 1.0], &[0.1, 0.2], &[0.0, 0.0], 0, 1),
            Err(Error::LstsqFit(_))
        ));
    }
}
