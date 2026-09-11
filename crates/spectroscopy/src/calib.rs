//! Calibration: quadratic energy bins (E6) and efficiency evaluation (E7).
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
/// anything else is an error. Coefficients are caller-supplied: the
/// upstream module has no coefficient-fitting routine, so neither does
/// this crate.
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
}
