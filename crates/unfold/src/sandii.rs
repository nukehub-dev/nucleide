//! SAND-II iterative spectral adjustment (McElroy et al., AFWL-TR-67-41, 1967).
//!
//! The forward problem folds a spectrum `phi` (one value per energy group)
//! through a response matrix `R` (one row per detector/reaction) into
//! calculated rates `c` (S1); unfolding adjusts a caller-supplied guess
//! spectrum against the measured rates `N` until the fold reproduces them.
//! SAND-II applies one multiplicative adjustment per iteration (S2–S3), each
//! a positivity-preserving weighted geometric mean of per-detector
//! correction factors, so the spectrum stays positive and a guess that
//! already fits is a fixed point:
//!
//! ```text
//! c_i   = Σ_j R_ij · φ_j                                (S1 fold)
//! W_ji  = R_ij · φ_j / c_i                              (S2 detector-i rate share from group j)
//! φ_j  ← φ_j · exp( Σ_i W_ji · ln(N_i / c_i) / Σ_i W_ji )   (S3 adjustment)
//! ```
//!
//! `W_ij` is the canonical base SAND-II weight; the original
//! implementation's optional per-detector statistics factor (`N_i²/σ_i`)
//! is not applied, so adjustments are not Poisson-weighted by detector
//! uncertainty.
//!
//! Clean-room from the public-domain report (US government work); the
//! adjustment form is the one reproduced across the open unfolding
//! literature. Detector/detector-reaction labelling, energy-group bounds,
//! and every response value are caller-supplied: evaluated libraries
//! (IRDFF and friends) are IAEA-copyright and are never vendored — a
//! runtime-download pack stays a later decision under the FGR-15 precedent.
//!
//! Convergence contract: after each adjustment the largest per-group
//! relative change `max_j |φ_new − φ_old| / φ_old` is compared against
//! `tolerance`; the run converges when it drops strictly below it. The
//! adjustment cap (`max_iterations`) is explicit; exhausting it raises
//! [`Error::NotConverged`] — a hard fail, never a silent partial spectrum
//! (the `tritium` face-Newton precedent). Degenerate-input policy, all
//! named: a detector with a zero fold against a nonzero measurement makes
//! the rates unreachable ([`Error::RatesUnreachable`], detected up front
//! for zero response rows and mid-iteration if zero rates pin away a
//! detector's whole support); detectors folding to zero with zero
//! measurement contribute no weight; a group no detector responds to keeps
//! the guess (factor 1, exactly); a zero measurement pins its groups to
//! zero. Only the guess is required to be strictly positive — the update is
//! multiplicative, so the spectrum stays non-negative (strictly positive
//! when all measured rates are positive).

use crate::error::{Error, Result};

/// Default per-group relative-change convergence tolerance.
///
/// This is a library default, not a value from the report (the original
/// code exposed the criterion as a user knob).
pub const DEFAULT_TOLERANCE: f64 = 1e-3;

/// Default adjustment cap (iterations).
///
/// Library default; the report's cap was likewise a user knob.
pub const DEFAULT_MAX_ITERATIONS: usize = 200;

/// One SAND-II adjustment step as yielded by the [`SandII`] iterator: the
/// post-adjustment spectrum together with its convergence diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct Iteration {
    /// 1-based adjustment number (the first adjustment is `1`).
    pub index: usize,
    /// Spectrum after this adjustment (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the post-adjustment spectrum (S1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector at this state;
    /// all ones at a fixed point, and 0.0 for a pinned detector (zero
    /// measurement folding to zero — the 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Largest per-group relative change this adjustment produced,
    /// `max_j |φ_new − φ_old| / φ_old`.
    pub max_rel_change: f64,
}

/// SAND-II unfolding iterator: owns the working spectrum and yields one
/// [`Iteration`] per adjustment until the run converges, the adjustment cap
/// is exhausted, or the measurements prove unreachable (see [`SandII::halt`]).
///
/// Construction validates the full input set once
/// ([`SandII::new`]); borrowing keeps the crate usable without allocations
/// beyond the working vectors.
#[derive(Debug)]
pub struct SandII<'a> {
    response: &'a [Vec<f64>],
    rates: &'a [f64],
    spectrum: Vec<f64>,
    tolerance: f64,
    max_iterations: usize,
    completed: usize,
    converged: bool,
    halt: Option<Error>,
}

impl<'a> SandII<'a> {
    /// Validate the inputs and seed the iterator with the guess spectrum.
    ///
    /// `response` holds one row per detector (all rows one value per energy
    /// group, non-negative finite), `rates` one measured rate per detector
    /// (non-negative finite), and `guess` one strictly positive finite value
    /// per energy group. `tolerance` must be finite and positive and
    /// `max_iterations` at least 1.
    pub fn new(
        response: &'a [Vec<f64>],
        rates: &'a [f64],
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

impl Iterator for SandII<'_> {
    type Item = Iteration;

    fn next(&mut self) -> Option<Iteration> {
        if self.converged || self.halt.is_some() || self.completed >= self.max_iterations {
            return None;
        }
        // (S1) fold of the current spectrum.
        let folds = crate::fold(self.response, &self.spectrum);
        for (i, (&c, &n)) in folds.iter().zip(self.rates.iter()).enumerate() {
            if n > 0.0 && c <= 0.0 {
                // A positive measurement whose whole support was pinned to
                // zero by other (zero) measurements: unreachable.
                self.halt = Some(Error::RatesUnreachable { detector: i });
                return None;
            }
        }
        // (S2)/(S3) per-group weighted geometric mean of the correction
        // factors. Detectors folding to zero carry zero weight (skipped);
        // groups no active detector responds to keep the guess.
        let mut new_spectrum = Vec::with_capacity(self.spectrum.len());
        for (j, &phi) in self.spectrum.iter().enumerate() {
            let mut weighted_log = 0.0f64;
            let mut weight_sum = 0.0f64;
            for (i, row) in self.response.iter().enumerate() {
                let c = folds[i];
                if c <= 0.0 {
                    continue;
                }
                let w = row[j] * phi / c;
                if w == 0.0 {
                    // Zero response in this group: skipping also avoids the
                    // 0 * ln(0) = NaN hazard from zero-rate detectors.
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
        // A pinned detector (zero measurement, support pinned to zero) folds
        // to zero; spell its factor as 0.0 instead of the 0/0 NaN.
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

/// Converged SAND-II solution with its convergence diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Adjusted spectrum (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the adjusted spectrum (S1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector; all ones at a
    /// fixed point, and 0.0 for a pinned detector (zero measurement folding
    /// to zero — the 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Number of adjustments applied.
    pub iterations: usize,
    /// Tolerance the run converged under.
    pub tolerance: f64,
    /// Largest per-group relative change of the final adjustment.
    pub max_rel_change: f64,
}

/// Drive a [`SandII`] iterator to convergence and return the solution.
///
/// Non-convergence is a hard [`Error::NotConverged`] — no partial spectrum is
/// returned (the `tritium` face-Newton precedent). Inputs and options are
/// validated up front exactly as by [`SandII::new`].
pub fn unfold(
    response: &[Vec<f64>],
    rates: &[f64],
    guess: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<Solution> {
    let mut run = SandII::new(response, rates, guess, tolerance, max_iterations)?;
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

    #[test]
    fn forward_fold_matches_hand_matvec() {
        let response = synthetic_response(3, 4, 1.5);
        let spectrum = vec![1.0, 2.0, 3.0, 4.0];
        let rates = crate::forward_fold(&response, &spectrum).unwrap();
        for (i, row) in response.iter().enumerate() {
            let hand: f64 = row.iter().zip(&spectrum).map(|(r, p)| r * p).sum();
            assert_eq!(rates[i], hand);
        }
    }

    #[test]
    fn exact_guess_is_a_fixed_point() {
        let response = synthetic_response(4, 12, 2.0);
        let midpoints = log_midpoints(12, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sol = unfold(
            &response,
            &rates,
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
        let guess: Vec<f64> = truth.iter().map(|p| 3.0 * p).collect();
        let sol = unfold(&response, &rates, &guess, 1e-12, 10_000)
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
    fn reduces_error_and_reproduces_rates_when_underdetermined() {
        // The realistic case: fewer detectors than groups. The rates are
        // reproduced exactly and the recovered spectrum is strictly closer
        // to the truth than the guess was.
        let n_groups = 24;
        let response = synthetic_response(6, n_groups, 3.0);
        let midpoints = log_midpoints(n_groups, 1e-6, 12.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        // Non-uniform bias: part of it lies in the response null space and
        // must keep the guess; the recoverable part must still converge
        // (deterministic residual: guess error 1.0 -> ~0.40, a 2.5x cut).
        let guess: Vec<f64> = truth
            .iter()
            .enumerate()
            .map(|(j, p)| p * (1.3 + 0.7 * (j % 5) as f64 / 4.0))
            .collect();
        let sol = unfold(&response, &rates, &guess, 1e-9, 50_000)
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
            err(&sol.spectrum, &truth) < err(&guess, &truth) / 2.0,
            "unfolding must cut the guess error by 2x"
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
            let guess: Vec<f64> = shape.iter().map(|p| 5.0 * p).collect();
            let sol = unfold(&response, &rates, &guess, 1e-10, 50_000)
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
            let det_sol = unfold(&det_response, &det_rates, &guess, 1e-11, 100_000)
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
    fn unconstrained_group_keeps_guess_exactly() {
        // Detector 0 responds only to group 0; groups 1..4 see nothing and
        // must keep the guess values bit-for-bit.
        let response = vec![vec![1.0, 0.0, 0.0, 0.0, 0.0]];
        let rates = vec![4.0];
        let guess = vec![2.0, 7.0, 3.0, 9.0, 5.0];
        let sol = unfold(&response, &rates, &guess, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum[0], 4.0);
        assert_eq!(&sol.spectrum[1..], &guess[1..]);
        // First adjustment moves group 0 (change 1.0); the confirming one is
        // a no-op (change 0).
        assert_eq!(sol.iterations, 2);
        assert_eq!(sol.max_rel_change, 0.0);
    }

    #[test]
    fn zero_measurement_pins_its_groups_to_zero() {
        let response = vec![vec![1.0, 1.0]];
        let rates = vec![0.0];
        let guess = vec![2.0, 3.0];
        let sol = unfold(&response, &rates, &guess, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum, vec![0.0, 0.0]);
        // The pinned detector folds to zero against its zero measurement;
        // its diagnostics factor is spelled 0.0, never the 0/0 NaN.
        assert_eq!(sol.rate_factors, vec![0.0]);
        assert!(sol.rate_factors.iter().all(|f| f.is_finite()));
    }

    #[test]
    fn zero_response_row_with_zero_rate_is_skipped() {
        let response = vec![vec![1.0, 0.5], vec![0.0, 0.0]];
        let rates = vec![3.0, 0.0];
        let guess = vec![2.0, 2.0];
        let sol = unfold(&response, &rates, &guess, 1e-12, 100).unwrap();
        assert!((sol.rate_factors[0] - 1.0).abs() < 1e-12);
        assert!(sol.spectrum.iter().all(|p| p.is_finite() && *p > 0.0));
    }

    #[test]
    fn zero_response_row_with_nonzero_rate_is_unreachable() {
        let response = vec![vec![0.0, 0.0]];
        let rates = vec![1.0];
        let guess = vec![2.0, 2.0];
        assert_eq!(
            unfold(&response, &rates, &guess, 1e-6, 100),
            Err(Error::RatesUnreachable { detector: 0 })
        );
    }

    #[test]
    fn inconsistent_rates_halt_mid_iteration() {
        // Detector 0 (rate 0) pins group 0 to zero; detector 1 then has a
        // zero fold against a nonzero rate -> unreachable, named, hard fail.
        let response = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let rates = vec![0.0, 1.0];
        let guess = vec![2.0, 2.0];
        assert_eq!(
            unfold(&response, &rates, &guess, 1e-6, 100),
            Err(Error::RatesUnreachable { detector: 1 })
        );
    }

    #[test]
    fn exhausted_cap_is_a_hard_not_converged() {
        let response = synthetic_response(4, 8, 1.5);
        let midpoints = log_midpoints(8, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let guess: Vec<f64> = truth.iter().map(|p| 2.0 * p).collect();
        assert_eq!(
            unfold(&response, &rates, &guess, 1e-12, 1),
            Err(Error::NotConverged)
        );
    }

    #[test]
    fn iterator_yields_each_adjustment_with_diagnostics() {
        let response = synthetic_response(3, 6, 1.2);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let guess: Vec<f64> = truth.iter().map(|p| 2.0 * p).collect();
        let mut run = SandII::new(&response, &rates, &guess, 1e-9, 1000).unwrap();
        let mut seen = 0usize;
        let mut last_change = f64::INFINITY;
        for it in &mut run {
            seen += 1;
            assert_eq!(it.index, seen);
            assert_eq!(it.spectrum.len(), 6);
            assert_eq!(it.rates.len(), 3);
            assert_eq!(it.rate_factors.len(), 3);
            assert!(
                it.max_rel_change <= last_change,
                "changes must not increase"
            );
            last_change = it.max_rel_change;
            if seen > 200 {
                break;
            }
        }
        assert!(run.converged());
        assert!(seen > 1, "a biased guess needs more than one adjustment");
        assert!(last_change < 1e-9);
        // A converged iterator is exhausted.
        assert!(run.next().is_none());
        assert!(run.halt.is_none());
    }

    #[test]
    fn validation_names_the_offending_input() {
        let ok_response = synthetic_response(2, 3, 1.0);
        let ok_rates = vec![1.0, 2.0];
        let ok_guess = vec![1.0, 1.0, 1.0];
        assert!(matches!(
            SandII::new(&[], &ok_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadResponse("empty response matrix"))
        ));
        assert!(matches!(
            SandII::new(&[vec![1.0], vec![1.0, 2.0]], &ok_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadShape {
                what: "response row",
                expected: 1,
                got: 2
            })
        ));
        let bad = vec![vec![1.0, f64::NAN, 1.0], vec![1.0; 3]];
        assert!(matches!(
            SandII::new(&bad, &ok_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadResponse("non-finite entry"))
        ));
        let neg = vec![vec![1.0, -1.0, 1.0], vec![1.0; 3]];
        assert!(matches!(
            SandII::new(&neg, &ok_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadResponse("negative entry"))
        ));
        assert!(matches!(
            SandII::new(&ok_response, &[1.0], &ok_guess, 1e-3, 10),
            Err(Error::BadShape {
                what: "rates",
                expected: 2,
                got: 1
            })
        ));
        let nan_rates = vec![1.0, f64::INFINITY];
        assert!(matches!(
            SandII::new(&ok_response, &nan_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadRates("non-finite entry"))
        ));
        let neg_rates = vec![1.0, -1.0];
        assert!(matches!(
            SandII::new(&ok_response, &neg_rates, &ok_guess, 1e-3, 10),
            Err(Error::BadRates("negative entry"))
        ));
        assert!(matches!(
            SandII::new(&ok_response, &ok_rates, &[1.0, 1.0], 1e-3, 10),
            Err(Error::BadShape {
                what: "guess",
                expected: 3,
                got: 2
            })
        ));
        let zero_guess = vec![1.0, 0.0, 1.0];
        assert!(matches!(
            SandII::new(&ok_response, &ok_rates, &zero_guess, 1e-3, 10),
            Err(Error::BadGuess("entries must be strictly positive"))
        ));
        assert!(matches!(
            SandII::new(&ok_response, &ok_rates, &ok_guess, 0.0, 10),
            Err(Error::BadOption("tolerance must be finite and > 0"))
        ));
        assert!(matches!(
            SandII::new(&ok_response, &ok_rates, &ok_guess, 1e-3, 0),
            Err(Error::BadOption("max_iterations must be >= 1"))
        ));
    }

    #[test]
    fn error_display_strings() {
        assert_eq!(
            Error::BadResponse("negative entry").to_string(),
            "unfold: invalid response matrix: negative entry"
        );
        assert_eq!(
            Error::BadShape {
                what: "rates",
                expected: 2,
                got: 1
            }
            .to_string(),
            "unfold: shape mismatch: rates holds 1 entries, expected 2"
        );
        assert_eq!(
            Error::RatesUnreachable { detector: 3 }.to_string(),
            "unfold: measured rates are unreachable: detector 3 cannot be satisfied by any positive spectrum"
        );
        assert_eq!(
            Error::NotConverged.to_string(),
            "unfold: SAND-II adjustment did not converge within its iteration cap"
        );
    }

    #[test]
    fn error_implements_std_error_trait() {
        let err: &dyn std::error::Error = &Error::NotConverged;
        assert!(err.source().is_none());
    }
}
