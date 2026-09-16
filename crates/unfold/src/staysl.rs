//! STAYSL-class damped least-squares spectral adjustment (Perey, ORNL/TM-6062,
//! 1977).
//!
//! The second landed method in this crate, beside [`crate::sandii`]. Where
//! SAND-II multiplies the guess by per-detector correction factors,
//! STAYSL-class methods adjust it toward the least-squares fit of the
//! measured rates (Perey, ORNL/TM-6062, 1977 — US government work, public
//! domain; the convention is reproduced across the open unfolding
//! literature, e.g. Greenwood, EPJ Web Conf. 106 (2016) 07001, CC-BY, and
//! the PNNL-22253 user guide, which describe the modern code as a one-step
//! least-squares adjustment of a starting spectrum). Clean-room from those
//! published equations; the UMG/STAYSL program package itself (NEA/RSICC
//! gated) is never consulted.
//!
//! The forward fold is linear (T1), so the "linearization" of the
//! least-squares formulation is exact and the objective being minimized is
//! the quadratic anchored-damped form (T2): the adjusted spectrum stays
//! near the caller's starting spectrum (the prior) with a pull strength set
//! by the damping `λ`, while fitting the rates weighted by the inverse
//! squared measurement sigmas. The pinned equation set:
//!
//! ```text
//! c_i  = Σ_j R_ij · φ_j                                              (T1 fold)
//! φ†   = argmin_φ Σ_i w_i (c_i(φ) − N_i)² + λ Σ_j (φ_j − φ_guess,j)²  (T2 anchored damped objective)
//! w_i  = 1 / σ_i²
//! ```
//!
//! Each iterator cycle re-solves (T2) on the current active set as one
//! weighted least-squares system through the shared
//! [`nucleide_linalg::lstsq`] kernel (per the crate rule; there is no
//! second least-squares implementation anywhere). The cycle solves for the
//! increment `Δ = φ† − φ_k` about the current iterate, which keeps the
//! kernel call uniform across cycles:
//!
//! ```text
//! Δ_A  = argmin_Δ Σ_i w_i (Σ_{j∈A} R_ij Δ_j − r_i)²
//!                  + λ Σ_{j∈A} (Δ_j − (φ_guess,j − φ_k,j))²           (T3 increment form)
//! r_i  = N_i − c_i(φ_k)
//! φ_j ← φ_j + Δ_j for j ∈ A, clipped to zero and dropped from A
//!        when the update lands non-positive (T4 clip + active-set pin)
//! ```
//!
//! (T3) is built as the augmented design rows `[R_·A; I]` with targets
//! `[r; φ_guess,A − φ_k,A]` and weights `[w; λ]` — the identity rows are
//! the damping, the weight column carries `λ`. Because the anchor in (T2)
//! is fixed at the guess for the whole run, every cycle solves the same
//! anchored problem: once a cycle produces no clip and no movement above
//! the tolerance, the next cycle is a bit-for-bit no-op and the run is
//! converged. The iteration that remains — and the reason this is an
//! iterator rather than a one-shot solve — is the published active-set
//! convention: groups the solve drives non-positive are pinned to zero and
//! the system is re-solved on the reduced set (Greenwood 2016: the
//! least-squares process is a one-step solve, repeated while constraints
//! change), so a run lasts as many cycles as the clipping needs plus one
//! confirming no-op.
//!
//! Weighting convention (v1, pinned): the caller supplies one strictly
//! positive sigma per detector and the weight is `w_i = 1/σ_i²`, the
//! classical measurement-uncertainty weight. These are the same sigmas the
//! landed UQ-lite conventions attach to the rates; this method consumes
//! them deterministically as weights only — no covariance propagation.
//!
//! Damping stance (v1, pinned): `λ` is absolute Tikhonov damping on the
//! distance to the guess — the fixed point of the unconstrained problem is
//! the ridge-regression solution `φ† = (Rᵀ W R + λ I)⁻¹ (Rᵀ W N + λ
//! φ_guess)`, so energy-group directions the data constrains much stronger
//! than `λ` are fit to the rates, directions constrained much weaker than
//! `λ` keep the guess, and the crossover is gradual. Smaller `λ` fits the
//! rates harder on weakly constrained (possibly ill-conditioned)
//! directions at the cost of a wilder solve; larger `λ` holds the spectrum
//! closer to the prior. Tune `λ` against the diagonal scale of `Rᵀ W R`.
//!
//! Convergence contract: after each cycle the largest per-group relative
//! change `max_j |φ_new − φ_old| / |φ_old|` (the 0/0 case of a pinned
//! group spelling 0.0, a group moving away from an exact zero spelling
//! infinity) is compared against `tolerance`; the run converges when it
//! drops strictly below it. The cycle cap (`max_iterations`) is explicit;
//! exhausting it raises [`Error::NotConverged`] — a hard fail, never a
//! silent partial spectrum (the `tritium` face-Newton precedent, same as
//! [`crate::sandii`]).
//!
//! Degenerate-input policy, all named: a detector with a zero response row
//! against a nonzero measurement is [`Error::RatesUnreachable`], detected
//! up front; an active set emptied by clipping against nonzero residual
//! rates is likewise [`Error::RatesUnreachable`]; a kernel failure surfaces
//! as [`Error::Lstsq`]; a group no active detector responds to keeps the
//! guess exactly (every solve returns its increment to zero — the damping
//! ties it to the prior and the data never pulls); zero measurements pull
//! their groups toward zero, and a group is pinned to zero exactly when
//! that pull drives it non-positive (T4) — the least-squares analogue of
//! SAND-II's zero-measurement pinning. Only the guess is required to be
//! strictly positive; additive updates may legitimately cross zero
//! mid-cycle, which is what the active set exists to catch. Pinned groups
//! stay pinned for the run (NNLS-style gradient re-activation stays a
//! later decision), and [`crate::sandii`] remains the strictly
//! positivity-preserving method.
//!
//! Divergences from [`crate::sandii`], pinned: the update is additive
//! (never multiplicative); the fixed point reproduces the rates up to the
//! damping slack (`λ`-scaled) rather than exactly; caller sigmas weight
//! the detectors (base SAND-II carries no per-detector statistics factor);
//! clipping, not multiplication, keeps groups non-negative. What both
//! methods share: the caller-supplied response matrix (evaluated libraries
//! such as IRDFF are IAEA-copyright and never vendored — the IRDFF-II v1
//! pack ships as a runtime download), the per-group
//! relative-change convergence contract, and the loud named errors.

use crate::error::{Error, Result};
use nucleide_linalg::lstsq;

/// Default per-group relative-change convergence tolerance.
///
/// Library default mirroring [`crate::sandii::DEFAULT_TOLERANCE`]; the
/// original codes exposed the criterion as a user knob.
pub const DEFAULT_TOLERANCE: f64 = 1e-3;

/// Default cycle cap (iterations).
///
/// Library default mirroring [`crate::sandii::DEFAULT_MAX_ITERATIONS`].
pub const DEFAULT_MAX_ITERATIONS: usize = 200;

/// Default damping `λ` for the anchored objective (T2).
///
/// Library default, not a value from the reports: it damps lightly against
/// response matrices whose sigma-weighted normal diagonals are `O(1)`. Tune
/// it against the diagonal scale of `Rᵀ W R` — smaller fits the rates
/// harder on weakly constrained directions, larger holds the spectrum
/// closer to the prior (and the run closer to a no-op).
pub const DEFAULT_DAMPING: f64 = 1e-3;

/// One STAYSL-class cycle as yielded by the [`Staysl`] iterator: the
/// post-solve spectrum together with its convergence diagnostics. Same
/// shape as [`crate::sandii::Iteration`].
#[derive(Debug, Clone, PartialEq)]
pub struct Iteration {
    /// 1-based cycle number (the first solve is `1`).
    pub index: usize,
    /// Spectrum after this cycle (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the post-cycle spectrum (T1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector at this state;
    /// all ones at an undamped fixed point, and 0.0 for a detector folding
    /// to zero (the 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Largest per-group relative change this cycle produced,
    /// `max_j |φ_new − φ_old| / |φ_old|`.
    pub max_rel_change: f64,
}

/// STAYSL-class least-squares unfolding iterator: owns the working spectrum,
/// the guess anchor, and the active set, and yields one [`Iteration`] per
/// damped solve until the run converges, the cycle cap is exhausted, or the
/// measurements prove unreachable (see [`Staysl::halt`]).
///
/// Construction validates the full input set once ([`Staysl::new`]);
/// borrowing keeps the crate usable without allocations beyond the working
/// vectors and the per-cycle augmented design matrix.
#[derive(Debug)]
pub struct Staysl<'a> {
    response: &'a [Vec<f64>],
    rates: &'a [f64],
    sigmas: &'a [f64],
    guess: Vec<f64>,
    spectrum: Vec<f64>,
    active: Vec<bool>,
    tolerance: f64,
    max_iterations: usize,
    damping: f64,
    completed: usize,
    converged: bool,
    halt: Option<Error>,
}

impl<'a> Staysl<'a> {
    /// Validate the inputs and seed the iterator with the guess spectrum.
    ///
    /// `response` holds one row per detector (all rows one value per energy
    /// group, non-negative finite), `rates` one measured rate per detector
    /// (non-negative finite), `sigmas` one strictly positive finite
    /// measurement sigma per detector (the weight is `1/σ²`), and `guess`
    /// one strictly positive finite value per energy group. `tolerance`
    /// must be finite and positive, `max_iterations` at least 1, and
    /// `damping` finite and positive.
    pub fn new(
        response: &'a [Vec<f64>],
        rates: &'a [f64],
        sigmas: &'a [f64],
        guess: &[f64],
        tolerance: f64,
        max_iterations: usize,
        damping: f64,
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
        for &sigma in sigmas {
            if !sigma.is_finite() {
                return Err(Error::BadRates("non-finite sigma"));
            }
            if sigma <= 0.0 {
                return Err(Error::BadRates("sigma must be > 0"));
            }
            let weight = 1.0 / (sigma * sigma);
            if !weight.is_finite() || weight <= 0.0 {
                return Err(Error::BadRates("sigma is out of weight range"));
            }
        }
        crate::validate_spectrum("guess", guess, n, true, Error::BadGuess)?;
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::BadOption("tolerance must be finite and > 0"));
        }
        if max_iterations == 0 {
            return Err(Error::BadOption("max_iterations must be >= 1"));
        }
        if !damping.is_finite() || damping <= 0.0 {
            return Err(Error::BadOption("damping must be finite and > 0"));
        }
        Ok(Self {
            response,
            rates,
            sigmas,
            guess: guess.to_vec(),
            spectrum: guess.to_vec(),
            active: vec![true; n],
            tolerance,
            max_iterations,
            damping,
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

    /// Damping `λ` applied to the anchored objective (T2).
    pub fn damping(&self) -> f64 {
        self.damping
    }

    /// Number of cycles applied so far.
    pub fn iterations(&self) -> usize {
        self.completed
    }

    /// Whether the last applied cycle met the tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// Current working spectrum (the guess before the first cycle).
    pub fn spectrum(&self) -> &[f64] {
        &self.spectrum
    }
}

impl Iterator for Staysl<'_> {
    type Item = Iteration;

    fn next(&mut self) -> Option<Iteration> {
        if self.converged || self.halt.is_some() || self.completed >= self.max_iterations {
            return None;
        }
        let n_active = self.active.iter().filter(|a| **a).count();
        if n_active == 0 {
            // Every group is pinned to zero; the zero spectrum is the only
            // spectrum left, and it satisfies the measurements only when all
            // rates vanish.
            if self.rates.iter().all(|&n| n == 0.0) {
                self.completed += 1;
                self.converged = true;
                let rates = crate::fold(self.response, &self.spectrum);
                return Some(Iteration {
                    index: self.completed,
                    spectrum: self.spectrum.clone(),
                    rates,
                    rate_factors: vec![0.0; self.rates.len()],
                    max_rel_change: 0.0,
                });
            }
            let detector = self.rates.iter().position(|&n| n > 0.0).unwrap_or(0);
            self.halt = Some(Error::RatesUnreachable { detector });
            return None;
        }
        // (T1)/(T2) residual of the current iterate against the rates.
        let folds = crate::fold(self.response, &self.spectrum);
        let mut residual: Vec<f64> = self.rates.iter().zip(&folds).map(|(n, c)| n - c).collect();
        // (T3) damped least-squares step on the active set, anchored at the
        // guess, as the augmented system [R_·A; I] Δ = [r; φ_guess − φ_k]
        // with weights [w; λ] through the shared kernel.
        let m = self.response.len();
        let mut design: Vec<Vec<f64>> = Vec::with_capacity(m + n_active);
        for row in self.response {
            design.push(
                row.iter()
                    .zip(&self.active)
                    .filter_map(|(r, &a)| a.then_some(*r))
                    .collect(),
            );
        }
        let mut k = 0usize;
        for j in 0..self.spectrum.len() {
            if !self.active[j] {
                continue;
            }
            let mut damp_row = vec![0.0; n_active];
            damp_row[k] = 1.0;
            design.push(damp_row);
            // Damping targets ride along in the target vector below.
            residual.push(self.guess[j] - self.spectrum[j]);
            k += 1;
        }
        let mut weights: Vec<f64> = self.sigmas.iter().map(|s| 1.0 / (s * s)).collect();
        weights.resize(m + n_active, self.damping);
        let step = match lstsq::weighted_lstsq(&design, &residual, &weights) {
            Ok(step) => step,
            Err(err) => {
                self.halt = Some(Error::Lstsq(err.to_string()));
                return None;
            }
        };
        // (T4) additive update; groups landing non-positive are pinned to
        // zero and leave the active set (the published clip-and-re-solve
        // convention — the re-solve on the reduced set is the next cycle).
        let previous = self.spectrum.clone();
        let mut k = 0usize;
        for j in 0..self.spectrum.len() {
            if !self.active[j] {
                continue;
            }
            let updated = self.spectrum[j] + step[k];
            k += 1;
            if updated <= 0.0 {
                self.spectrum[j] = 0.0;
                self.active[j] = false;
            } else {
                self.spectrum[j] = updated;
            }
        }
        let max_rel_change = previous
            .iter()
            .zip(&self.spectrum)
            .map(|(old, new)| {
                if *old != 0.0 {
                    ((new - old) / old).abs()
                } else if *new == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            })
            .fold(0.0f64, f64::max);
        self.completed += 1;
        if max_rel_change < self.tolerance {
            self.converged = true;
        }
        // Diagnostics from the post-cycle state.
        let rates = crate::fold(self.response, &self.spectrum);
        // A detector folding to zero against a zero measurement spells its
        // factor as 0.0 instead of the 0/0 NaN.
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

/// Converged STAYSL-class solution with its convergence diagnostics. Same
/// shape as [`crate::sandii::Solution`].
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Adjusted spectrum (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the adjusted spectrum (T1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector; all ones at an
    /// undamped fixed point, and 0.0 for a detector folding to zero (the
    /// 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Number of cycles applied.
    pub iterations: usize,
    /// Tolerance the run converged under.
    pub tolerance: f64,
    /// Largest per-group relative change of the final cycle.
    pub max_rel_change: f64,
}

/// Drive a [`Staysl`] iterator to convergence and return the solution.
///
/// Non-convergence is a hard [`Error::NotConverged`] — no partial spectrum
/// is returned (the `tritium` face-Newton precedent). Inputs and options are
/// validated up front exactly as by [`Staysl::new`].
pub fn unfold(
    response: &[Vec<f64>],
    rates: &[f64],
    sigmas: &[f64],
    guess: &[f64],
    tolerance: f64,
    max_iterations: usize,
    damping: f64,
) -> Result<Solution> {
    let mut run = Staysl::new(
        response,
        rates,
        sigmas,
        guess,
        tolerance,
        max_iterations,
        damping,
    )?;
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
    // `max_iterations >= 1` is validated, so a converged run applied at
    // least one cycle; the `None` arm is unreachable-in-practice defense
    // (loud error, never a panic).
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
            DEFAULT_DAMPING,
        )
        .expect("exact guess must converge");
        assert_eq!(sol.iterations, 1, "one no-op solve pins the fixed point");
        assert_eq!(sol.max_rel_change, 0.0);
        assert_eq!(sol.spectrum, truth);
        assert!(sol.rate_factors.iter().all(|f| (*f - 1.0).abs() < 1e-12));
    }

    #[test]
    fn recovers_determined_system_from_biased_guess() {
        // Nearly-diagonal 6x6 system: every group is seen, and with light
        // damping the single solve lands at the ridge solution, round-off
        // away from the truth.
        let response = synthetic_response(6, 6, 0.9);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 3.0 * p).collect();
        // Calibrated: at damping 1e-12 the ridge fixed point lands 2.1e-8 off
        // the truth (rate factors to 1.5e-11), two cycles.
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-9, 10_000, 1e-12)
            .expect("determined system must converge");
        for (rec, want) in sol.spectrum.iter().zip(truth.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-6,
                "recovered {rec} vs truth {want}"
            );
        }
        assert!(sol.rate_factors.iter().all(|f| (*f - 1.0).abs() < 1e-6));
    }

    #[test]
    fn caller_sigmas_are_the_weights() {
        // Two redundant readings of the same group, one precise (2.0 ± 0.01)
        // and one sloppy (2.2 ± 1.0): the weighted fit must follow the
        // precise detector, not the average.
        let response = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let rates = vec![2.0, 2.2];
        let sigmas = vec![0.01, 1.0];
        let guess = vec![1.5, 7.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-12, 100, 1e-8).unwrap();
        assert!(
            (sol.spectrum[0] - 2.0).abs() < 0.02,
            "precise detector must dominate, got {}",
            sol.spectrum[0]
        );
        // The unseen second group keeps the guess exactly.
        assert_eq!(sol.spectrum[1], 7.0);
    }

    #[test]
    fn underdetermined_run_reproduces_the_rates() {
        // The realistic case: fewer detectors than groups. The least-squares
        // fixed point interpolates the measured rates (the physically
        // measured quantities) up to the damping slack. Unlike SAND-II's
        // multiplicative spreading — whose weighted geometric mean happened
        // to cut this particular bias by 2.5x — the arithmetic fixed point
        // does NOT promise a spectrum closer to the truth: group directions
        // the data constrains only weakly take their values from the fit
        // (soft-direction min-norm content), not from the prior, so no
        // error-cut assertion is pinned here. That divergence from SAND-II
        // is deliberate and documented in the module rustdoc.
        let n_groups = 24;
        let response = synthetic_response(6, n_groups, 3.0);
        let midpoints = log_midpoints(n_groups, 1e-6, 12.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth
            .iter()
            .enumerate()
            .map(|(j, p)| p * (1.3 + 0.7 * (j % 5) as f64 / 4.0))
            .collect();
        // Calibrated: rates reproduced to 3.9e-6 in six cycles (one clipping
        // wave plus the confirming no-op).
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-9, 100, 1e-6)
            .expect("underdetermined round trip must converge");
        let rate_err = sol
            .rate_factors
            .iter()
            .map(|f| (*f - 1.0f64).abs())
            .fold(0.0f64, f64::max);
        assert!(rate_err < 1e-4, "rates reproduced to 1e-4, got {rate_err}");
        assert!(sol.spectrum.iter().all(|p| p.is_finite()));
    }

    #[test]
    fn irdff_ii_analytical_benchmark_shapes_recover() {
        // IRDFF-II (Trkov et al., Nucl. Data Sheets 163 (2020) 1; open access
        // arXiv:1909.03336) lists analytical benchmark-field shapes; three of
        // them are plain published facts reused here as test inputs with
        // citation: a thermal Maxwellian at 293.6 K, a pure 1/E field, and a
        // Maxwellian at 25 keV. The library's tabulated group spectra are
        // IAEA-copyright data files and are deliberately NOT used.
        // Group bounds bracket each shape's support so every group sees a
        // significant flux: with bounds extending decades past the shape's
        // decay the trailing groups sit at the 1e-30 positivity floor, and a
        // pinned-to-zero recovered group would read as a spurious relative
        // error of exactly 1.0 against that floor.
        let n_groups = 18;
        let thermal_kt = 8.617333262e-5 * 293.6; // 293.6 K in eV -> MeV
        let thermal: Vec<f64> = log_midpoints(n_groups, 1e-9, 1e-3)
            .iter()
            .map(|&e| e * (-e / thermal_kt).exp() + 1e-30)
            .collect();
        let inv_e: Vec<f64> = log_midpoints(n_groups, 1e-6, 10.0)
            .iter()
            .map(|&e| 1.0 / e + 1e-30)
            .collect();
        let fusion_kt = 2.5e-2; // 25 keV in MeV
        let fusion: Vec<f64> = log_midpoints(n_groups, 1e-3, 0.3)
            .iter()
            .map(|&e| e.sqrt() * (-e / fusion_kt).exp() + 1e-30)
            .collect();
        for (name, shape) in [("thermal", thermal), ("1/E", inv_e), ("25keV", fusion)] {
            let det_response = synthetic_response(n_groups, n_groups, 0.8);
            let det_rates = crate::forward_fold(&det_response, &shape).unwrap();
            let sigmas = unit_sigmas(det_rates.len());
            let guess: Vec<f64> = shape.iter().map(|p| 5.0 * p).collect();
            // Calibrated: at damping 1e-10 the determined ridge recovery is
            // 3.8e-7 (thermal), 3.1e-6 (1/E), 1.5e-8 (25 keV) — recovery
            // error scales with the damping, the pinned tolerance sits a
            // factor ~3-30 below the measured values.
            let det_sol = unfold(&det_response, &det_rates, &sigmas, &guess, 1e-9, 100, 1e-10)
                .unwrap_or_else(|e| panic!("determined IRDFF probe {name} must converge: {e}"));
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
    fn sandii_fixed_point_is_a_least_squares_fixed_point() {
        // Both methods consume the same fold (T1 == S1): once SAND-II has
        // converged (rate factors pinned to 1, residual ~ 0), the anchored
        // least-squares cycle (T3) against that residual is a no-op about
        // the SAND-II spectrum as the anchor, so the SAND-II fixed point is
        // a STAYSL-class fixed point too.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let sandii_sol = crate::sandii::unfold(&response, &rates, &guess, 1e-12, 100_000)
            .expect("SAND-II must converge");
        let sigmas = unit_sigmas(rates.len());
        let sol = unfold(
            &response,
            &rates,
            &sigmas,
            &sandii_sol.spectrum,
            1e-4,
            100,
            1e-6,
        )
        .expect("least-squares step at a SAND-II fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op solve");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(sandii_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "least-squares solve must hold the SAND-II fixed point"
            );
        }
    }

    #[test]
    fn unconstrained_group_keeps_guess_exactly() {
        // Detector 0 responds only to group 0; groups 1..4 see nothing and
        // must keep the guess values bit-for-bit. With damping the fixed
        // point is the damped optimum (w·N + λ·guess)/(w + λ), not N.
        let response = vec![vec![1.0, 0.0, 0.0, 0.0, 0.0]];
        let rates = vec![4.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 7.0, 3.0, 9.0, 5.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-3, 100, 1e-3).unwrap();
        let want0 = (1.0 * 4.0 + 1e-3 * 2.0) / (1.0 + 1e-3);
        assert!(
            (sol.spectrum[0] - want0).abs() < 1e-6,
            "damped optimum {want0}, got {}",
            sol.spectrum[0]
        );
        assert_eq!(&sol.spectrum[1..], &guess[1..]);
        // First cycle moves group 0 (change ~1.0); the confirming one is a
        // no-op (change below tolerance).
        assert_eq!(sol.iterations, 2);
    }

    #[test]
    fn zero_measurement_pulls_its_groups_toward_zero_and_clips() {
        // The zero rate pulls both groups toward zero; the pull overshoots
        // group 0 into the active set (T4), and the re-solve on the reduced
        // set settles at the clipped damped optimum — the least-squares
        // analogue of SAND-II's zero-measurement pinning.
        let response = vec![vec![1.0, 1.0]];
        let rates = vec![0.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 3.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-3, 100, 1e-3).unwrap();
        assert_eq!(sol.spectrum[0], 0.0, "overshot group pinned to zero");
        assert!(
            sol.spectrum[1] > 0.0 && sol.spectrum[1] < 1e-2,
            "remaining group settles near zero, got {}",
            sol.spectrum[1]
        );
        assert!(
            sol.rates[0].abs() < 1e-2,
            "fold decays toward zero, got {}",
            sol.rates[0]
        );
        assert_eq!(sol.rate_factors, vec![0.0]);
        assert!(sol.spectrum.iter().all(|p| p.is_finite()));
    }

    #[test]
    fn zero_response_row_with_zero_rate_is_skipped() {
        let response = vec![vec![1.0, 0.5], vec![0.0, 0.0]];
        let rates = vec![3.0, 0.0];
        let sigmas = vec![1.0, 1.0];
        let guess = vec![2.0, 2.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, 1e-9, 100, 1e-6).unwrap();
        assert!((sol.rate_factors[0] - 1.0).abs() < 1e-3);
        assert!(sol.spectrum.iter().all(|p| p.is_finite()));
    }

    #[test]
    fn zero_response_row_with_nonzero_rate_is_unreachable() {
        let response = vec![vec![0.0, 0.0]];
        let rates = vec![1.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 2.0];
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, 1e-6, 100, 1e-3),
            Err(Error::RatesUnreachable { detector: 0 })
        );
    }

    #[test]
    fn exhausted_cap_is_a_hard_not_converged() {
        // The unconstrained-group case needs two cycles at tolerance 1e-3;
        // a cap of one exhausts mid-run — a hard fail, never a partial.
        let response = vec![vec![1.0, 0.0, 0.0]];
        let rates = vec![4.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 7.0, 3.0];
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, 1e-3, 1, 1e-3),
            Err(Error::NotConverged)
        );
    }

    #[test]
    fn iterator_yields_each_cycle_with_diagnostics() {
        let response = synthetic_response(3, 6, 1.2);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.0 * p).collect();
        let mut run = Staysl::new(&response, &rates, &sigmas, &guess, 1e-9, 1000, 1e-6).unwrap();
        let mut seen = 0usize;
        for it in &mut run {
            seen += 1;
            assert_eq!(it.index, seen);
            assert_eq!(it.spectrum.len(), 6);
            assert_eq!(it.rates.len(), 3);
            assert_eq!(it.rate_factors.len(), 3);
            if seen > 200 {
                break;
            }
        }
        assert!(run.converged());
        assert!(seen > 1, "a biased guess needs more than one cycle");
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
            Staysl::new(&[], &ok_rates, &ok_sigmas, &ok_guess, 1e-3, 10, 1e-3),
            Err(Error::BadResponse("empty response matrix"))
        ));
        assert!(matches!(
            Staysl::new(
                &[vec![1.0], vec![1.0, 2.0]],
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadShape {
                what: "response row",
                expected: 1,
                got: 2
            })
        ));
        let neg = vec![vec![1.0, -1.0, 1.0], vec![1.0; 3]];
        assert!(matches!(
            Staysl::new(&neg, &ok_rates, &ok_sigmas, &ok_guess, 1e-3, 10, 1e-3),
            Err(Error::BadResponse("negative entry"))
        ));
        assert!(matches!(
            Staysl::new(&ok_response, &[1.0], &ok_sigmas, &ok_guess, 1e-3, 10, 1e-3),
            Err(Error::BadShape {
                what: "rates",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Staysl::new(&ok_response, &ok_rates, &[1.0], &ok_guess, 1e-3, 10, 1e-3),
            Err(Error::BadShape {
                what: "sigmas",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &[1.0, f64::NAN],
                &ok_guess,
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadRates("non-finite sigma"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &[1.0, 0.0],
                &ok_guess,
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadRates("sigma must be > 0"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &[1.0, 1e-300],
                &ok_guess,
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadRates("sigma is out of weight range"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &[1.0, 1.0],
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadShape {
                what: "guess",
                expected: 3,
                got: 2
            })
        ));
        let zero_guess = vec![1.0, 0.0, 1.0];
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &zero_guess,
                1e-3,
                10,
                1e-3
            ),
            Err(Error::BadGuess("entries must be strictly positive"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                0.0,
                10,
                1e-3
            ),
            Err(Error::BadOption("tolerance must be finite and > 0"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                1e-3,
                0,
                1e-3
            ),
            Err(Error::BadOption("max_iterations must be >= 1"))
        ));
        assert!(matches!(
            Staysl::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                1e-3,
                10,
                0.0
            ),
            Err(Error::BadOption("damping must be finite and > 0"))
        ));
    }

    #[test]
    fn error_display_strings() {
        assert_eq!(
            Error::Lstsq("boom".to_string()).to_string(),
            "unfold: least-squares solve failed: boom"
        );
    }
}
