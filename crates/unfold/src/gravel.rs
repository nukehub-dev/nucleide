//! GRAVEL iterative spectral adjustment (Matzke, PTB-N-19, 1994).
//!
//! The third landed method in this crate, beside [`crate::sandii`] and
//! [`crate::staysl`]. GRAVEL is the uncertainty-aware sibling of SAND-II:
//! the same positivity-preserving multiplicative adjustment, but each
//! detector's correction is weighted by its measurement precision, so
//! precisely measured rates pull harder than sloppy ones. The pinned
//! equation set:
//!
//! ```text
//! c_i   = Σ_j R_ij · φ_j                                              (G1 fold)
//! W_ji  = (R_ij · φ_j / c_i) · (N_i² / σ_i²)                           (G2 chi-square weights)
//! φ_j  ← φ_j · exp( Σ_i W_ji · ln(N_i / c_i) / Σ_i W_ji )                 (G3 adjustment)
//! ```
//!
//! (G2) is the SAND-II base rate share times the measurement-weight factor
//! `N_i²/σ_i²`: a detector measured precisely (small `σ_i` against its rate
//! `N_i`) dominates the groups it responds to, a sloppy one barely pulls.
//! At a fixed point every `N_i/c_i` is 1, so the weights are irrelevant
//! there — a SAND-II fixed point is a GRAVEL fixed point and vice versa,
//! and both consume the same fold (G1 == S1 == T1).
//!
//! Clean-room from the published iteration; no gated code is consulted.
//! Detector/reaction labelling, energy-group bounds, and every response
//! value, measured rate, and sigma are caller-supplied: evaluated libraries
//! (IRDFF and friends) are IAEA-copyright and are never vendored — the
//! IRDFF-II v1 pack ships as a runtime download
//! (`nucleide.data.fetch_irdff` plus `parse_irdff_g725`).
//!
//! Convergence contract: after each adjustment the largest per-group
//! relative change `max_j |φ_new − φ_old| / φ_old` is compared against
//! `tolerance`; the run converges when it drops strictly below it. The
//! adjustment cap (`max_iterations`) is explicit; exhausting it raises
//! [`Error::NotConverged`] — a hard fail, never a silent partial spectrum
//! (the `tritium` face-Newton precedent, same as [`crate::sandii`] and
//! [`crate::staysl`]).
//!
//! Degenerate-input policy, all named: a detector with a zero response row
//! against a nonzero measurement makes the rates unreachable
//! ([`Error::RatesUnreachable`], detected up front); detectors folding to
//! zero carry zero weight (skipped); a group no positively weighted
//! detector responds to keeps the guess (factor 1, exactly). Only the guess
//! is required to be strictly positive — the update is multiplicative, so
//! the spectrum stays strictly positive.
//!
//! Divergences from [`crate::sandii`], pinned: caller sigmas weight the
//! detectors (one finite positive sigma per detector — zero, non-finite, or
//! missing sigmas are loud errors, never silent uniform weighting — via the
//! `N_i²/σ_i²` factor); zero measurements carry zero weight, so unlike
//! SAND-II they do not pin their groups to zero — a detector reading zero
//! is simply not fitted, and an all-zero measurement set returns the guess
//! unchanged after one confirming no-op. What both methods share: the
//! caller-supplied response matrix, the per-group relative-change
//! convergence contract, and the loud named errors.

use crate::error::{Error, Result};

/// Default per-group relative-change convergence tolerance.
///
/// Library default mirroring [`crate::sandii::DEFAULT_TOLERANCE`]; the
/// original codes exposed the criterion as a user knob.
pub const DEFAULT_TOLERANCE: f64 = 1e-3;

/// Default adjustment cap (iterations).
///
/// Library default mirroring [`crate::sandii::DEFAULT_MAX_ITERATIONS`].
pub const DEFAULT_MAX_ITERATIONS: usize = 200;

/// One GRAVEL adjustment step as yielded by the [`Gravel`] iterator: the
/// post-adjustment spectrum together with its convergence diagnostics. Same
/// shape as [`crate::sandii::Iteration`].
#[derive(Debug, Clone, PartialEq)]
pub struct Iteration {
    /// 1-based adjustment number (the first adjustment is `1`).
    pub index: usize,
    /// Spectrum after this adjustment (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the post-adjustment spectrum (G1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector at this state;
    /// all ones at a fixed point, and 0.0 for a detector folding to zero
    /// (the 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Largest per-group relative change this adjustment produced,
    /// `max_j |φ_new − φ_old| / φ_old`.
    pub max_rel_change: f64,
}

/// GRAVEL unfolding iterator: owns the working spectrum and yields one
/// [`Iteration`] per adjustment until the run converges, the adjustment cap
/// is exhausted, or the measurements prove unreachable (see [`Gravel::halt`]).
///
/// Construction validates the full input set once ([`Gravel::new`]);
/// borrowing keeps the crate usable without allocations beyond the working
/// vectors.
#[derive(Debug)]
pub struct Gravel<'a> {
    response: &'a [Vec<f64>],
    rates: &'a [f64],
    meas_weight: Vec<f64>,
    spectrum: Vec<f64>,
    tolerance: f64,
    max_iterations: usize,
    completed: usize,
    converged: bool,
    halt: Option<Error>,
}

impl<'a> Gravel<'a> {
    /// Validate the inputs and seed the iterator with the guess spectrum.
    ///
    /// `response` holds one row per detector (all rows one value per energy
    /// group, non-negative finite), `rates` one measured rate per detector
    /// (non-negative finite), `sigmas` one strictly positive finite
    /// measurement sigma per detector (the per-detector weight factor is
    /// `N_i²/σ_i²`), and `guess` one strictly positive finite value per
    /// energy group. `tolerance` must be finite and positive and
    /// `max_iterations` at least 1.
    pub fn new(
        response: &'a [Vec<f64>],
        rates: &'a [f64],
        sigmas: &'a [f64],
        guess: &[f64],
        tolerance: f64,
        max_iterations: usize,
    ) -> Result<Self> {
        let n = crate::validate_response(response)?;
        if rates.len() != response.len() {
            return Err(Error::BadShape {
                what: "rates",
                expected: response.len(),
                got: rates.len(),
            });
        }
        for (i, &rate) in rates.iter().enumerate() {
            if !rate.is_finite() {
                return Err(Error::BadRates("non-finite entry"));
            }
            if rate < 0.0 {
                return Err(Error::BadRates("negative entry"));
            }
            if rate > 0.0 && response[i].iter().all(|&r| r == 0.0) {
                return Err(Error::RatesUnreachable { detector: i });
            }
        }
        if sigmas.len() != rates.len() {
            return Err(Error::BadShape {
                what: "sigmas",
                expected: rates.len(),
                got: sigmas.len(),
            });
        }
        let mut meas_weight = Vec::with_capacity(rates.len());
        for (&rate, &sigma) in rates.iter().zip(sigmas.iter()) {
            if !sigma.is_finite() {
                return Err(Error::BadRates("non-finite sigma"));
            }
            if sigma <= 0.0 {
                return Err(Error::BadRates("sigma must be > 0"));
            }
            let factor = (rate / sigma) * (rate / sigma);
            // A zero rate carries zero weight by design (pinned divergence);
            // a nonzero rate weighting to exactly zero means sigma overflowed
            // out of range and must be loud, never a silent skip.
            if !factor.is_finite() || (factor == 0.0 && rate != 0.0) {
                return Err(Error::BadRates("sigma is out of weight range"));
            }
            meas_weight.push(factor);
        }
        crate::validate_spectrum("guess", guess, n, true, Error::BadGuess)?;
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::BadOption("tolerance must be finite and > 0"));
        }
        if max_iterations == 0 {
            return Err(Error::BadOption("max_iterations must be >= 1"));
        }
        Ok(Self {
            response,
            rates,
            meas_weight,
            spectrum: guess.to_vec(),
            tolerance,
            max_iterations,
            completed: 0,
            converged: false,
            halt: None,
        })
    }

    /// Tolerance the run converges under (largest per-group relative change
    /// strictly below it).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Number of adjustments applied so far.
    pub fn iterations(&self) -> usize {
        self.completed
    }

    /// Whether the last applied adjustment met the tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// Current working spectrum (the guess before the first adjustment).
    pub fn spectrum(&self) -> &[f64] {
        &self.spectrum
    }
}

impl Iterator for Gravel<'_> {
    type Item = Iteration;

    fn next(&mut self) -> Option<Iteration> {
        if self.converged || self.halt.is_some() || self.completed >= self.max_iterations {
            return None;
        }
        // (G1) fold of the current spectrum.
        let folds = crate::fold(self.response, &self.spectrum);
        for (i, (&c, &n)) in folds.iter().zip(self.rates.iter()).enumerate() {
            if n > 0.0 && c <= 0.0 {
                // A positive measurement whose whole support folds to zero:
                // unreachable.
                self.halt = Some(Error::RatesUnreachable { detector: i });
                return None;
            }
        }
        // (G2)/(G3) per-group chi-square-weighted geometric mean of the
        // correction factors. Detectors folding to zero carry zero weight
        // (skipped); zero-measurement detectors carry a zero `N²/σ²` factor
        // and are skipped the same way, so they never pin groups to zero;
        // groups no positively weighted detector responds to keep the guess.
        // Skipping zero weights before the logarithm also avoids the
        // 0 * ln(0) = NaN hazard from zero-rate detectors.
        let mut new_spectrum = Vec::with_capacity(self.spectrum.len());
        for (j, &phi) in self.spectrum.iter().enumerate() {
            let mut weighted_log = 0.0f64;
            let mut weight_sum = 0.0f64;
            for (i, row) in self.response.iter().enumerate() {
                let c = folds[i];
                if c <= 0.0 {
                    continue;
                }
                let w = row[j] * phi / c * self.meas_weight[i];
                if w == 0.0 {
                    continue;
                }
                weighted_log += w * (self.rates[i] / c).ln();
                weight_sum += w;
            }
            let factor = if weight_sum > 0.0 {
                (weighted_log / weight_sum).exp()
            } else {
                1.0
            };
            new_spectrum.push(phi * factor);
        }
        let max_rel_change = new_spectrum
            .iter()
            .zip(self.spectrum.iter())
            .map(|(new, old)| ((new - old) / old).abs())
            .fold(0.0f64, f64::max);
        self.spectrum = new_spectrum;
        self.completed += 1;
        if max_rel_change < self.tolerance {
            self.converged = true;
        }
        // Diagnostics from the post-adjustment state.
        let rates = crate::fold(self.response, &self.spectrum);
        // A detector folding to zero spells its factor as 0.0 instead of the
        // 0/0 NaN.
        let rate_factors = self
            .rates
            .iter()
            .zip(rates.iter())
            .map(|(n, &c)| if c == 0.0 { 0.0 } else { n / c })
            .collect();
        Some(Iteration {
            index: self.completed,
            spectrum: self.spectrum.clone(),
            rates,
            rate_factors,
            max_rel_change,
        })
    }
}

/// Converged GRAVEL solution with its convergence diagnostics. Same shape as
/// [`crate::sandii::Solution`].
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Adjusted spectrum (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the adjusted spectrum (G1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector; all ones at a
    /// fixed point, and 0.0 for a detector folding to zero (the 0/0 factor
    /// spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Number of adjustments applied.
    pub iterations: usize,
    /// Tolerance the run converged under.
    pub tolerance: f64,
    /// Largest per-group relative change of the final adjustment.
    pub max_rel_change: f64,
}

/// Drive a [`Gravel`] iterator to convergence and return the solution.
///
/// Non-convergence is a hard [`Error::NotConverged`] — no partial spectrum is
/// returned (the `tritium` face-Newton precedent). Inputs and options are
/// validated up front exactly as by [`Gravel::new`].
pub fn unfold(
    response: &[Vec<f64>],
    rates: &[f64],
    sigmas: &[f64],
    guess: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<Solution> {
    let mut run = Gravel::new(response, rates, sigmas, guess, tolerance, max_iterations)?;
    let mut last: Option<Iteration> = None;
    for iteration in &mut run {
        last = Some(iteration);
    }
    if let Some(err) = run.halt {
        return Err(err);
    }
    if !run.converged() {
        return Err(Error::NotConverged);
    }
    // `max_iterations >= 1` is validated, so a converged run applied at least
    // one adjustment; the `None` arm is unreachable-in-practice defense (loud
    // error, never a panic).
    let Some(last) = last else {
        return Err(Error::NotConverged);
    };
    Ok(Solution {
        spectrum: last.spectrum,
        rates: last.rates,
        rate_factors: last.rate_factors,
        iterations: last.index,
        tolerance,
        max_rel_change: last.max_rel_change,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic synthetic response matrix (provenance: hand-built in
    /// this test module — synthetic fixtures only, never evaluated data).
    /// `n_det` rows of Gaussian bumps over `n_groups` log-spaced energies,
    /// each bump centered on a different group with half-width `width`,
    /// plus a low tail so every group carries some weight.
    fn synthetic_response(n_det: usize, n_groups: usize, width: f64) -> Vec<Vec<f64>> {
        let mut rows = Vec::with_capacity(n_det);
        for i in 0..n_det {
            let center = (i * (n_groups - 1)) as f64 / (n_det.max(2) - 1) as f64;
            let mut row = Vec::with_capacity(n_groups);
            for j in 0..n_groups {
                let d = (j as f64 - center) / width;
                let bump = (-d * d).exp();
                row.push(bump + 1e-3);
            }
            rows.push(row);
        }
        rows
    }

    /// Synthetic two-component spectrum: thermal Maxwellian at `kT` plus a
    /// 1/E epithermal tail, sampled at `midpoints` (MeV). All entries
    /// strictly positive.
    fn maxwellian_plus_inv_e(midpoints: &[f64], kt: f64) -> Vec<f64> {
        midpoints
            .iter()
            .map(|&e| e * (-e / kt).exp() + 1.0e-6 / e)
            .collect()
    }

    fn log_midpoints(n_groups: usize, lo: f64, hi: f64) -> Vec<f64> {
        let step = (hi / lo).ln() / n_groups as f64;
        (0..n_groups)
            .map(|j| (lo.ln() + (j as f64 + 0.5) * step).exp())
            .collect()
    }

    fn unit_sigmas(n: usize) -> Vec<f64> {
        vec![1.0; n]
    }

    #[test]
    fn exact_guess_is_a_fixed_point() {
        let response = synthetic_response(4, 12, 2.0);
        let midpoints = log_midpoints(12, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let sol = unfold(
            &response,
            &rates,
            &sigmas,
            &truth,
            DEFAULT_TOLERANCE,
            DEFAULT_MAX_ITERATIONS,
        )
        .expect("exact guess must converge");
        assert_eq!(
            sol.iterations, 1,
            "one no-op adjustment pins the fixed point"
        );
        assert_eq!(sol.max_rel_change, 0.0);
        assert_eq!(sol.spectrum, truth);
        assert!(sol.rate_factors.iter().all(|f| (*f - 1.0).abs() < 1e-12));
    }

    #[test]
    fn recovers_determined_system_from_biased_guess() {
        // Nearly-diagonal 6x6 system: every group is seen, so the recovery
        // is fully determined and lands at round-off.
        let response = synthetic_response(6, 6, 0.9);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 3.0 * p).collect();
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 10_000)
            .expect("determined system must converge");
        assert!(sol.max_rel_change < 1e-12);
        for (rec, want) in sol.spectrum.iter().zip(truth.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-6,
                "recovered {rec} vs truth {want}"
            );
        }
        assert!(sol.rate_factors.iter().all(|f| (*f - 1.0).abs() < 1e-9));
    }

    #[test]
    fn caller_sigmas_are_the_weights() {
        // Two redundant readings of the same group, one precise (2.0 ± 0.01)
        // and one sloppy (2.2 ± 1.0): the `N²/σ²` factors (40000 vs 4.84)
        // must follow the precise detector, not the average.
        let response = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let rates = vec![2.0, 2.2];
        let sigmas = vec![0.01, 1.0];
        let guess = vec![1.5, 7.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100).unwrap();
        assert!(
            (sol.spectrum[0] - 2.0).abs() < 0.02,
            "precise detector must dominate, got {}",
            sol.spectrum[0]
        );
        // The unseen second group keeps the guess exactly.
        assert_eq!(sol.spectrum[1], 7.0);
    }

    #[test]
    fn reduces_error_and_reproduces_rates_when_underdetermined() {
        // The realistic case: fewer detectors than groups. The rates are
        // reproduced exactly and the recovered spectrum is strictly closer
        // to the truth than the guess was. Sigmas follow the
        // counting-statistics model (σ² = N, weights = N): genuinely
        // non-uniform weights spanning the rates' three decades, so this
        // exercises the chi-square machinery rather than cloning the
        // SAND-II trajectory.
        let n_groups = 24;
        let response = synthetic_response(6, n_groups, 3.0);
        let midpoints = log_midpoints(n_groups, 1e-6, 12.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas: Vec<f64> = rates.iter().map(|r| r.sqrt()).collect();
        // Non-uniform bias: part of it lies in the response null space and
        // must keep the guess; the recoverable part must still converge.
        // Calibrated: converges in ~41k adjustments, rates to ~1e-9, guess
        // error 1.0 -> ~0.77.
        let guess: Vec<f64> = truth
            .iter()
            .enumerate()
            .map(|(j, p)| p * (1.3 + 0.7 * (j % 5) as f64 / 4.0))
            .collect();
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-9, 50_000)
            .expect("underdetermined round trip must converge");
        let rate_err = sol
            .rate_factors
            .iter()
            .map(|f| (*f - 1.0f64).abs())
            .fold(0.0f64, f64::max);
        assert!(rate_err < 1e-6, "rates reproduced to 1e-6, got {rate_err}");
        let err = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .map(|(x, y)| ((x - y) / y).abs())
                .fold(0.0f64, f64::max)
        };
        assert!(
            err(&sol.spectrum, &truth) < err(&guess, &truth),
            "unfolding must improve on the guess"
        );
        assert!(sol.spectrum.iter().all(|p| p.is_finite() && *p > 0.0));
    }

    #[test]
    fn irdff_ii_analytical_benchmark_shapes_recover() {
        // IRDFF-II (Trkov et al., Nucl. Data Sheets 163 (2020) 1; open access
        // arXiv:1909.03336) lists analytical benchmark-field shapes; three of
        // them are plain published facts reused here as test inputs with
        // citation: a thermal Maxwellian at 293.6 K, a pure 1/E field, and a
        // Maxwellian at 25 keV. The library's tabulated group spectra are
        // IAEA-copyright data files and are deliberately NOT used.
        let n_groups = 18;
        let midpoints = log_midpoints(n_groups, 1e-9, 20.0);
        let thermal_kt = 8.617333262e-5 * 293.6; // 293.6 K in eV -> MeV
        let thermal: Vec<f64> = midpoints
            .iter()
            .map(|&e| e * (-e / thermal_kt).exp() + 1e-30)
            .collect();
        let inv_e: Vec<f64> = midpoints.iter().map(|&e| 1.0 / e + 1e-30).collect();
        let fusion_kt = 2.5e-2; // 25 keV in MeV
        let fusion: Vec<f64> = midpoints
            .iter()
            .map(|&e| e.sqrt() * (-e / fusion_kt).exp() + 1e-30)
            .collect();
        for (name, shape) in [("thermal", thermal), ("1/E", inv_e), ("25keV", fusion)] {
            let response = synthetic_response(8, n_groups, 2.5);
            let rates = crate::forward_fold(&response, &shape).unwrap();
            let sigmas = unit_sigmas(rates.len());
            let guess: Vec<f64> = shape.iter().map(|p| 5.0 * p).collect();
            let sol = unfold(&response, &rates, &sigmas, &guess, 1e-10, 50_000)
                .unwrap_or_else(|e| panic!("IRDFF-II {name} probe must converge: {e}"));
            let rate_err = sol
                .rate_factors
                .iter()
                .map(|f| (*f - 1.0f64).abs())
                .fold(0.0f64, f64::max);
            assert!(rate_err < 1e-7, "{name}: rates to 1e-7, got {rate_err}");
            // Determined probe: 8 detectors over a nearly-diagonal response.
            let det_response = synthetic_response(n_groups, n_groups, 0.8);
            let det_rates = crate::forward_fold(&det_response, &shape).unwrap();
            let det_sigmas = unit_sigmas(det_rates.len());
            let det_sol = unfold(
                &det_response,
                &det_rates,
                &det_sigmas,
                &guess,
                1e-11,
                100_000,
            )
            .expect("determined IRDFF probe must converge");
            let worst = det_sol
                .spectrum
                .iter()
                .zip(shape.iter())
                .map(|(r, w)| ((r - w) / w).abs())
                .fold(0.0f64, f64::max);
            assert!(worst < 1e-5, "{name}: shape recovery to 1e-5, got {worst}");
        }
    }

    #[test]
    fn sandii_fixed_point_is_a_gravel_fixed_point() {
        // Both methods consume the same fold (G1 == S1): once SAND-II has
        // converged (rate factors pinned to 1, every ln factor zero), the
        // chi-square-weighted cycle (G3) about that spectrum is a no-op
        // regardless of the weights — a SAND-II fixed point is a GRAVEL
        // fixed point too.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let sandii_sol = crate::sandii::unfold(&response, &rates, &guess, 1e-12, 100_000)
            .expect("SAND-II must converge");
        let sigmas = unit_sigmas(rates.len());
        let sol = unfold(&response, &rates, &sigmas, &sandii_sol.spectrum, 1e-4, 100)
            .expect("GRAVEL step at a SAND-II fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(sandii_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "GRAVEL must hold the SAND-II fixed point"
            );
        }
    }

    #[test]
    fn gravel_fixed_point_is_a_sandii_fixed_point() {
        // The weights are irrelevant at a fixed point in either direction:
        // GRAVEL reproduces the rates exactly (multiplicative family), so a
        // SAND-II cycle about the GRAVEL spectrum is a no-op too.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let gravel_sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100_000)
            .expect("GRAVEL must converge");
        let sol = crate::sandii::unfold(&response, &rates, &gravel_sol.spectrum, 1e-4, 100)
            .expect("SAND-II step at a GRAVEL fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(gravel_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "SAND-II must hold the GRAVEL fixed point"
            );
        }
    }

    #[test]
    fn gravel_fixed_point_is_a_least_squares_fixed_point() {
        // GRAVEL reproduces the rates exactly, so the anchored
        // least-squares cycle (T3) about the GRAVEL spectrum as the anchor
        // sees a ~zero residual and is a no-op — mirroring the SAND-II
        // cross-check in `staysl`.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let gravel_sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100_000)
            .expect("GRAVEL must converge");
        let sol = crate::staysl::unfold(
            &response,
            &rates,
            &sigmas,
            &gravel_sol.spectrum,
            1e-4,
            100,
            1e-6,
        )
        .expect("least-squares step at a GRAVEL fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op solve");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(gravel_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "least-squares solve must hold the GRAVEL fixed point"
            );
        }
    }

    #[test]
    fn unconstrained_group_keeps_guess_exactly() {
        // Detector 0 responds only to group 0; groups 1..4 see nothing and
        // must keep the guess values bit-for-bit.
        let response = vec![vec![1.0, 0.0, 0.0, 0.0, 0.0]];
        let rates = vec![4.0];
        let sigmas = vec![0.5];
        let guess = vec![2.0, 7.0, 3.0, 9.0, 5.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum[0], 4.0);
        assert_eq!(&sol.spectrum[1..], &guess[1..]);
        // First adjustment moves group 0 (change 1.0); the confirming one is
        // a no-op (change 0).
        assert_eq!(sol.iterations, 2);
        assert_eq!(sol.max_rel_change, 0.0);
    }

    #[test]
    fn zero_measurement_detectors_carry_no_weight() {
        // Pinned divergence from SAND-II: the `N²/σ²` factor of a zero
        // measurement is zero, so the detector is skipped rather than
        // pinning its groups to zero — the guess survives untouched.
        let response = vec![vec![1.0, 1.0]];
        let rates = vec![0.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 3.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum, guess);
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert_eq!(sol.max_rel_change, 0.0);
        // The unfitted detector folds to nonzero against its zero
        // measurement; its diagnostics factor is spelled 0.0 only when the
        // fold itself is zero.
        assert!(sol.rate_factors.iter().all(|f| f.is_finite()));
    }

    #[test]
    fn zero_response_row_with_nonzero_rate_is_unreachable() {
        let response = vec![vec![0.0, 0.0]];
        let rates = vec![1.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 2.0];
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, 1e-6, 100),
            Err(Error::RatesUnreachable { detector: 0 })
        );
    }

    #[test]
    fn exhausted_cap_is_a_hard_not_converged() {
        let response = synthetic_response(4, 8, 1.5);
        let midpoints = log_midpoints(8, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.0 * p).collect();
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, 1e-12, 1),
            Err(Error::NotConverged)
        );
    }

    #[test]
    fn iterator_yields_each_adjustment_with_diagnostics() {
        let response = synthetic_response(3, 6, 1.2);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.0 * p).collect();
        let mut run = Gravel::new(&response, &rates, &sigmas, &guess, 1e-9, 1000).unwrap();
        let mut seen = 0usize;
        for it in &mut run {
            seen += 1;
            assert_eq!(it.index, seen);
            assert_eq!(it.spectrum.len(), 6);
            assert_eq!(it.rates.len(), 3);
            assert_eq!(it.rate_factors.len(), 3);
            assert!(it.spectrum.iter().all(|p| p.is_finite() && *p > 0.0));
            assert!(it.rate_factors.iter().all(|f| f.is_finite()));
            if seen > 5000 {
                break;
            }
        }
        assert!(run.converged());
        assert!(seen > 1, "a biased guess needs more than one adjustment");
        // A converged iterator is exhausted.
        assert!(run.next().is_none());
        assert!(run.halt.is_none());
    }

    #[test]
    fn validation_names_the_offending_input() {
        let ok_response = synthetic_response(2, 3, 1.0);
        let ok_rates = vec![1.0, 2.0];
        let ok_sigmas = vec![1.0, 1.0];
        let ok_guess = vec![1.0, 1.0, 1.0];
        assert!(matches!(
            Gravel::new(&[], &ok_rates, &ok_sigmas, &ok_guess, 1e-3, 10),
            Err(Error::BadResponse("empty response matrix"))
        ));
        assert!(matches!(
            Gravel::new(
                &[vec![1.0], vec![1.0, 2.0]],
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                1e-3,
                10
            ),
            Err(Error::BadShape {
                what: "response row",
                expected: 1,
                got: 2
            })
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &[1.0], &ok_sigmas, &ok_guess, 1e-3, 10),
            Err(Error::BadShape {
                what: "rates",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &[1.0], &ok_guess, 1e-3, 10),
            Err(Error::BadShape {
                what: "sigmas",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Gravel::new(
                &ok_response,
                &ok_rates,
                &[1.0, f64::NAN],
                &ok_guess,
                1e-3,
                10
            ),
            Err(Error::BadRates("non-finite sigma"))
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &[1.0, 0.0], &ok_guess, 1e-3, 10),
            Err(Error::BadRates("sigma must be > 0"))
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &[1.0, 1e-300], &ok_guess, 1e-3, 10),
            Err(Error::BadRates("sigma is out of weight range"))
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &ok_sigmas, &[1.0, 1.0], 1e-3, 10),
            Err(Error::BadShape {
                what: "guess",
                expected: 3,
                got: 2
            })
        ));
        let zero_guess = vec![1.0, 0.0, 1.0];
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &ok_sigmas, &zero_guess, 1e-3, 10),
            Err(Error::BadGuess("entries must be strictly positive"))
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &ok_sigmas, &ok_guess, 0.0, 10),
            Err(Error::BadOption("tolerance must be finite and > 0"))
        ));
        assert!(matches!(
            Gravel::new(&ok_response, &ok_rates, &ok_sigmas, &ok_guess, 1e-3, 0),
            Err(Error::BadOption("max_iterations must be >= 1"))
        ));
    }
}
