//! Kinetic data (`betas`, `lambdas`, `Lambda`) and equilibrium initials (E2).
//!
//! All data are caller-supplied in SI-adjacent units: `betas` dimensionless
//! delayed-neutron fractions, `lambdas` precursor decay constants \[1/s\],
//! `lambda_gen` (`Λ`) the prompt-neutron generation time \[s\], reactivity in
//! Δk. No tabulated nuclear data ship with this crate.

use crate::error::Error;

/// Validated point-kinetics data for `G >= 1` precursor groups.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticParams {
    /// Per-group delayed-neutron fractions (dimensionless, `>= 0`).
    betas: Vec<f64>,
    /// Per-group precursor decay constants \[1/s\] (`> 0`).
    lambdas: Vec<f64>,
    /// Prompt-neutron generation time Λ \[s\] (`> 0`).
    lambda_gen: f64,
    /// Cached `sum(betas)`.
    beta_total: f64,
}

impl KineticParams {
    /// Validate raw kinetic data.
    ///
    /// Rejects: empty groups, length mismatch, non-finite entries,
    /// `lambda_i <= 0`, `beta_i < 0`, `sum(beta) <= 0` or `>= 1` (a total
    /// delayed fraction at or above unity is not a delayed-supercritical
    /// kinetics problem), and non-positive/non-finite `Λ`.
    pub fn new(betas: Vec<f64>, lambdas: Vec<f64>, lambda_gen: f64) -> Result<Self, Error> {
        if betas.is_empty() {
            return Err(Error::EmptyGroups);
        }
        if betas.len() != lambdas.len() {
            return Err(Error::LengthMismatch {
                betas: betas.len(),
                lambdas: lambdas.len(),
            });
        }
        if !lambda_gen.is_finite() || lambda_gen <= 0.0 {
            return Err(Error::BadData("Lambda must be finite and > 0"));
        }
        for (b, l) in betas.iter().zip(&lambdas) {
            if !b.is_finite() || !l.is_finite() {
                return Err(Error::BadData("betas and lambdas must be finite"));
            }
            if *b < 0.0 {
                return Err(Error::BadData("betas must be >= 0"));
            }
            if *l <= 0.0 {
                return Err(Error::BadData("lambdas must be > 0"));
            }
        }
        let beta_total: f64 = betas.iter().sum();
        if beta_total <= 0.0 || beta_total >= 1.0 {
            return Err(Error::BadData("sum(betas) must satisfy 0 < sum < 1"));
        }
        Ok(Self {
            betas,
            lambdas,
            lambda_gen,
            beta_total,
        })
    }

    /// Build from OpenMC iterated-fission-probability (IFP) kinetics data.
    ///
    /// OpenMC's IFP estimator reports effective delayed fractions (`betas`)
    /// and the generation time (`lambda_gen`) but **no** precursor decay
    /// constants: the caller supplies `lambdas` from the same delayed-neutron
    /// data library the IFP run used. Validation is identical to [`Self::new`].
    pub fn from_ifp(betas: Vec<f64>, lambda_gen: f64, lambdas: Vec<f64>) -> Result<Self, Error> {
        Self::new(betas, lambdas, lambda_gen)
    }

    /// Number of precursor groups `G`.
    pub fn groups(&self) -> usize {
        self.betas.len()
    }

    /// Per-group delayed-neutron fractions.
    pub fn betas(&self) -> &[f64] {
        &self.betas
    }

    /// Per-group precursor decay constants \[1/s\].
    pub fn lambdas(&self) -> &[f64] {
        &self.lambdas
    }

    /// Prompt-neutron generation time Λ \[s\].
    pub fn lambda_gen(&self) -> f64 {
        self.lambda_gen
    }

    /// Total delayed-neutron fraction `β = sum(beta_i)`.
    pub fn beta_total(&self) -> f64 {
        self.beta_total
    }

    /// Equilibrium precursor populations for a steady state at `n0` (E2):
    ///
    /// ```text
    /// C_i(0) = beta_i / (lambda_i * Lambda) * n0.
    /// ```
    ///
    /// Setting `dC_i/dt = 0` in (E1b) at constant reactivity gives
    /// `C_i = beta_i n / (lambda_i Lambda)`; this is the discrete form of
    /// the upstream driver initial conditions (power normalized, precursors
    /// at the steady-state ratio). Rejects non-finite or negative `n0`.
    pub fn equilibrium_precursors(&self, n0: f64) -> Result<Vec<f64>, Error> {
        if !n0.is_finite() || n0 < 0.0 {
            return Err(Error::BadState("n0 must be finite and >= 0"));
        }
        Ok(self
            .betas
            .iter()
            .zip(&self.lambdas)
            .map(|(b, l)| b / (l * self.lambda_gen) * n0)
            .collect())
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn one_group() -> KineticParams {
        KineticParams::new(vec![0.0065], vec![0.08], 1e-4).unwrap()
    }

    pub(crate) fn six_group() -> KineticParams {
        // Synthetic 6-group set (hand-picked round numbers, NOT a data
        // library evaluation): betas sum to 0.0065.
        KineticParams::new(
            vec![0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032],
            vec![0.01, 0.03, 0.1, 0.3, 1.0, 3.0],
            1e-5,
        )
        .unwrap()
    }

    #[test]
    fn accepts_one_six_eight_groups() {
        assert_eq!(one_group().groups(), 1);
        assert_eq!(six_group().groups(), 6);
        let eight = KineticParams::new(vec![0.0008; 8], vec![0.1; 8], 1e-5).unwrap();
        assert_eq!(eight.groups(), 8);
        assert!((eight.beta_total() - 0.0064).abs() < 1e-15);
    }

    #[test]
    fn equilibrium_matches_e2_ratio() {
        let p = one_group();
        let c = p.equilibrium_precursors(1.0).unwrap();
        assert_eq!(c.len(), 1);
        let want = 0.0065 / (0.08 * 1e-4);
        assert!((c[0] - want).abs() / want < 1e-12);
    }

    #[test]
    fn rejects_bad_data() {
        assert_eq!(
            KineticParams::new(vec![], vec![], 1e-5),
            Err(Error::EmptyGroups)
        );
        assert_eq!(
            KineticParams::new(vec![0.001], vec![0.1, 0.2], 1e-5),
            Err(Error::LengthMismatch {
                betas: 1,
                lambdas: 2
            })
        );
        assert!(KineticParams::new(vec![0.001], vec![0.1], 0.0).is_err());
        assert!(KineticParams::new(vec![0.001], vec![0.0], 1e-5).is_err());
        assert!(KineticParams::new(vec![-0.001], vec![0.1], 1e-5).is_err());
        assert!(KineticParams::new(vec![0.6, 0.6], vec![0.1, 0.2], 1e-5).is_err());
        assert!(KineticParams::new(vec![0.0], vec![0.1], 1e-5).is_err());
        assert!(KineticParams::new(vec![f64::NAN], vec![0.1], 1e-5).is_err());
        assert!(six_group().equilibrium_precursors(-1.0).is_err());
    }

    #[test]
    fn from_ifp_matches_new() {
        let a = KineticParams::new(vec![0.0065], vec![0.08], 1e-4).unwrap();
        let b = KineticParams::from_ifp(vec![0.0065], 1e-4, vec![0.08]).unwrap();
        assert_eq!(a, b);
    }
}
