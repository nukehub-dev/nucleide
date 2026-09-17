//! MAXED maximum-entropy spectral adjustment (Reginatto & Goldhagen,
//! Health Phys. 77 (1999) 579; NIM A 476 (2002) 242).
//!
//! The fourth and last landed method in this crate, beside
//! [`crate::sandii`], [`crate::staysl`], and [`crate::gravel`]: one method
//! per adjustment family — geometric-mean (SAND-II), least-squares
//! (STAYSL-class), chi-square-weighted (GRAVEL), maximum entropy (MAXED).
//! Where SAND-II/GRAVEL spread each detector's correction geometrically and
//! STAYSL-class solves a damped least-squares system in the spectrum, MAXED
//! maximizes the Shannon relative entropy against the guess-as-prior
//! subject to the chi-square fit, so among the spectra that reproduce the
//! rates it keeps the one closest to the guess in the information sense.
//! The pinned equation set:
//!
//! ```text
//! c_i   = Σ_j R_ij · φ_j                                              (M1 fold)
//! S(φ)  = Σ_j (φ_j − φ⁰_j − φ_j · ln(φ_j / φ⁰_j))                        (M2 relative entropy)
//! χ²(φ) = Σ_i ((c_i − N_i) / σ_i)²                                     (M3 chi-square)
//! φ_j(λ) = φ⁰_j · exp(−Σ_i R_ij · λ_i)                                   (M4 exponential family)
//! J_ik  = Σ_j R_ij · φ_j · R_kj;                                         (M5 dual curvature)
//!         min_Δ Σ_{i ∈ fitted} (Σ_k J_ik · Δ_k − (c_i − N_i))² / σ_i²
//!         min_δ Σ_{i ∈ fitted} (Σ_k H_ik · δ_k − v_i)² / σ'ᵢ² + μ·‖Δ‖²        (M6 damped step)
//!             H_ik = J_ik / (p_i·q_k), v_i = (c_i − N_i) / p_i,
//!             1/σ'ᵢ² = p_i² / σ_i², Δ_k = δ_k / q_k,
//!             accept iff χ² drops and max_j|ln(φ_new_j/φ_old_j)| ≤ 1      (M6 cap)
//! ```
//!
//! (M2) is the relative entropy with the guess spectrum as the default
//! model `φ⁰`: `S ≤ 0` with equality only at `φ = φ⁰`. Its Lagrange
//! stationary form under rate constraints is the exponential family (M4) in
//! the per-detector multipliers `λ`; each cycle linearizes the fold about
//! the current multipliers (`∂c_i/∂λ_k = −J_ik`) and takes the dual step
//! (M5/M6) — one weighted least-squares solve through the shared
//! [`nucleide_linalg::lstsq`] kernel (per the crate rule), over the fitted
//! detectors only, with the caller `1/σ_i²` weights carried through
//! untouched. (M6) equilibrates that solve: max-abs row scales `p` and
//! column scales `q` bring every entry to order unity before a
//! Levenberg-Marquardt damping term `μ` (adaptive: a step is accepted only
//! when it strictly lowers the chi-square (M3), shrinking `μ` toward
//! Newton; rejection grows `μ` toward scaled steepest descent). The row
//! scales are compensated into the weights (`p_i²/σ_i²`), so equilibration
//! changes neither the weights nor the fixed points (a zero step still
//! means the dual gradient `JᵀWr` vanishes bit-for-bit); it is needed
//! because the raw curvature (M5) inherits the spectrum's dynamic range —
//! a thermal guess spans thirty orders across groups — and redundant
//! detectors make it exactly singular, so the raw Gauss-Newton step is
//! both inaccurate and, at rank deficiency, stalled. For consistent rates
//! the multipliers converge to the maximum-entropy exact-fit spectrum; for
//! inconsistent (e.g. redundant disagreeing) readings they converge to the
//! minimum-chi-square member of the family.
//!
//! Clean-room from the published iteration; no gated code is consulted (the
//! UMG package itself is closed and is never touched — journal equations
//! only). Detector/reaction labelling, energy-group bounds, and every
//! response value, measured rate, and sigma are caller-supplied: evaluated
//! libraries (IRDFF and friends) are IAEA-copyright and are never
//! vendored — the IRDFF-II v1 pack ships as a runtime download
//! (`nucleide.data.fetch_irdff` plus `parse_irdff_g725`).
//!
//! Convergence contract: after each adjustment the largest per-group
//! relative change `max_j |φ_new − φ_old| / φ_old` is compared against
//! `tolerance`, and the chi-square (M3) against the caller `target_chi2`
//! (default: the detector count, the chi-square expectation). The run
//! converges when the change drops strictly below the tolerance while the
//! chi-square sits at or below its target. The adjustment cap
//! (`max_iterations`) is explicit; exhausting it — or converging in
//! relative change to a fit whose chi-square still exceeds the target —
//! raises [`Error::NotConverged`]: a hard fail, never a silent partial
//! spectrum (the `tritium` face-Newton precedent, same as [`crate::sandii`],
//! [`crate::staysl`], and [`crate::gravel`]).
//!
//! Degenerate-input policy, all named: a detector with a zero response row
//! against a nonzero measurement makes the rates unreachable
//! ([`Error::RatesUnreachable`], detected up front); detectors folding to
//! zero carry no step (they drop out of the dual solve); a group no fitted
//! detector responds to keeps the guess bit-for-bit (its exponent is
//! exactly zero). Only the guess is required to be strictly positive — the
//! update is exponential, so the spectrum stays strictly positive.
//!
//! Divergences from [`crate::sandii`], pinned: caller sigmas weight the
//! detectors through the chi-square (one finite positive sigma per
//! detector — zero, non-finite, or out-of-weight-range sigmas are loud
//! errors, never silent uniform weighting — via the `1/σ_i²` factor);
//! zero measurements carry zero weight, so like GRAVEL (and unlike
//! SAND-II) they are simply not fitted — an all-zero measurement set
//! returns the guess unchanged after one confirming no-op. What all four
//! methods share: the caller-supplied response matrix, the per-group
//! relative-change convergence contract, and the loud named errors.

use crate::error::{Error, Result};
use nucleide_linalg::lstsq;

/// Default per-group relative-change convergence tolerance.
///
/// Library default mirroring [`crate::sandii::DEFAULT_TOLERANCE`]; the
/// original codes exposed the criterion as a user knob.
pub const DEFAULT_TOLERANCE: f64 = 1e-3;

/// Default adjustment cap (iterations).
///
/// Library default mirroring [`crate::sandii::DEFAULT_MAX_ITERATIONS`].
pub const DEFAULT_MAX_ITERATIONS: usize = 200;

/// Initial Levenberg-Marquardt damping (M6) for the equilibrated dual step.
///
/// Starts mild (near-Newton on order-unity equilibrated curvature, heavy
/// only against wild dynamic ranges): the first attempt of every
/// adjustment is close to the raw equilibrated Newton step, so
/// well-conditioned systems polish from full steps; damping only grows
/// through rejection, and only to tame overshoot. Fixed points are
/// independent of it.
const INITIAL_DAMPING: f64 = 1e-3;

/// Floor for the adaptive damping: stays strictly positive so the augmented
/// rows always restore full column rank (a zero damping would hand the
/// shared kernel the raw, possibly singular, curvature). At 1e-20 the
/// floor penalty is far below any data weight, so converged directions
/// finish unfreezed.
const MIN_DAMPING: f64 = 1e-20;

/// Ceiling for the adaptive damping; its square root still fits comfortably
/// in finite floating point as an augmented-row entry.
const MAX_DAMPING: f64 = 1e15;

/// Trust-region cap (M6): the largest accepted per-group log movement
/// `max_j |ln(φ_new_j / φ_old_j)|` per adjustment (a factor of `e`).
///
/// A full dual-Newton step from a biased guess can overshoot along
/// weakly-constrained multiplier directions — collapsing a low-signal
/// group by dozens of orders in one accepted step (the chi-square only
/// sees the well-fitted groups improve) — and the exponential family
/// cannot resurrect it afterwards. Capping each adjustment to a factor of
/// `e` keeps the path synchronous across groups (the GRAVEL precedent:
/// small multiplicative moves every cycle); the cap goes quiet near the
/// fit, where Newton steps are tiny, so the quadratic finish is
/// unaffected. Overshoot is absorbed by halving the step length, never by
/// inflating the damping: only a direction that cannot improve the fit at
/// any length grows `μ`.
///
/// Floor for the trust-region step-length halving: below it the trial is
/// the current point to working precision, so the direction is exhausted
/// and the damping grows instead.
const MAX_LOG_STEP: f64 = 1.0;
const ALPHA_MIN: f64 = 1e-15;

/// Dual-step retries per adjustment before the cycle yields a no-op (whose
/// zero relative change then trips the tolerance leg when the fit is
/// acceptable, or idles toward the cap — a hard [`Error::NotConverged`] —
/// when the chi-square target stays out of reach).
const DAMPING_RETRIES: usize = 25;

/// Default chi-square target: the detector count.
///
/// The chi-square expectation for `m` fitted detectors is `m`, so an
/// unspecified target accepts any fit at or below it. Pass an explicit
/// smaller target (down to 0.0) to demand a tighter fit; exact-fit
/// round-trips drive the chi-square to roundoff, far below any positive
/// target.
pub fn default_target_chi2(n_detectors: usize) -> f64 {
    n_detectors as f64
}

/// Relative entropy (M2) of `spectrum` against the `prior` default model.
///
/// `S = Σ_j (φ_j − φ⁰_j − φ_j · ln(φ_j / φ⁰_j)) ≤ 0`, with equality only at
/// `φ = φ⁰`. Both inputs must be strictly positive and finite with equal
/// length; violations are [`Error::BadResponse`]. Exposed so callers can
/// score unfold solutions in the information sense (the quantity MAXED
/// maximizes subject to the chi-square target).
pub fn relative_entropy(spectrum: &[f64], prior: &[f64]) -> Result<f64> {
    if spectrum.len() != prior.len() {
        return Err(Error::BadShape {
            what: "spectrum",
            expected: prior.len(),
            got: spectrum.len(),
        });
    }
    let mut entropy = 0.0f64;
    for (&phi, &phi0) in spectrum.iter().zip(prior.iter()) {
        if !phi.is_finite() || !phi0.is_finite() || phi <= 0.0 || phi0 <= 0.0 {
            return Err(Error::BadResponse(
                "entropy needs strictly positive spectra",
            ));
        }
        entropy += phi - phi0 - phi * (phi / phi0).ln();
    }
    Ok(entropy)
}

/// Chi-square (M3) of a folded rate vector against the measurements.
///
/// `χ² = Σ_{i ∈ fitted} ((c_i − N_i) / σ_i)²`; zero measurements carry
/// zero weight (skipped, the GRAVEL house rule). `inv_var` holds the
/// `1/σ_i²` weights.
pub(crate) fn chi_square(folds: &[f64], rates: &[f64], inv_var: &[f64]) -> f64 {
    folds
        .iter()
        .zip(rates.iter())
        .zip(inv_var.iter())
        .map(|((&c, &n), &w)| {
            if n == 0.0 {
                0.0
            } else {
                let r = c - n;
                r * r * w
            }
        })
        .sum()
}

/// One MAXED adjustment step as yielded by the [`Maxed`] iterator: the
/// post-adjustment spectrum together with its convergence diagnostics. Same
/// shape as [`crate::sandii::Iteration`].
#[derive(Debug, Clone, PartialEq)]
pub struct Iteration {
    /// 1-based adjustment number (the first adjustment is `1`).
    pub index: usize,
    /// Spectrum after this adjustment (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the post-adjustment spectrum (M1).
    pub rates: Vec<f64>,
    /// Measured-rate over folded-rate ratios per detector at this state;
    /// all ones at a fixed point, and 0.0 for a detector folding to zero
    /// (the 0/0 factor spelled as 0.0).
    pub rate_factors: Vec<f64>,
    /// Largest per-group relative change this adjustment produced,
    /// `max_j |φ_new − φ_old| / φ_old`.
    pub max_rel_change: f64,
}

/// MAXED maximum-entropy unfolding iterator: owns the dual multipliers and
/// the working spectrum, and yields one [`Iteration`] per adjustment until
/// the run converges, the adjustment cap is exhausted, or the measurements
/// prove unreachable (see the module docs for the acceptance leg).
///
/// Construction validates the full input set once ([`Maxed::new`]);
/// borrowing keeps the crate usable without allocations beyond the working
/// vectors.
#[derive(Debug)]
pub struct Maxed<'a> {
    response: &'a [Vec<f64>],
    rates: &'a [f64],
    inv_var: Vec<f64>,
    prior: Vec<f64>,
    lambda: Vec<f64>,
    spectrum: Vec<f64>,
    target_chi2: f64,
    tolerance: f64,
    max_iterations: usize,
    damping: f64,
    completed: usize,
    converged: bool,
    halt: Option<Error>,
}

impl<'a> Maxed<'a> {
    /// Validate the inputs and seed the iterator with the guess spectrum.
    ///
    /// `response` holds one row per detector (all rows one value per energy
    /// group, non-negative finite), `rates` one measured rate per detector
    /// (non-negative finite), `sigmas` one strictly positive finite
    /// measurement sigma per detector (the chi-square weight is `1/σ_i²`),
    /// and `guess` one strictly positive finite value per energy group —
    /// the maximum-entropy default model. `target_chi2` must be finite and
    /// non-negative (`None` selects [`default_target_chi2`]); `tolerance`
    /// must be finite and positive and `max_iterations` at least 1.
    pub fn new(
        response: &'a [Vec<f64>],
        rates: &'a [f64],
        sigmas: &'a [f64],
        guess: &[f64],
        target_chi2: Option<f64>,
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
        let mut inv_var = Vec::with_capacity(rates.len());
        for &sigma in sigmas.iter() {
            if !sigma.is_finite() {
                return Err(Error::BadRates("non-finite sigma"));
            }
            if sigma <= 0.0 {
                return Err(Error::BadRates("sigma must be > 0"));
            }
            let weight = 1.0 / (sigma * sigma);
            if !weight.is_finite() || weight == 0.0 {
                return Err(Error::BadRates("sigma is out of weight range"));
            }
            inv_var.push(weight);
        }
        crate::validate_spectrum("guess", guess, n, true, Error::BadGuess)?;
        let target_chi2 = target_chi2.unwrap_or_else(|| default_target_chi2(response.len()));
        if !target_chi2.is_finite() || target_chi2 < 0.0 {
            return Err(Error::BadOption("target_chi2 must be finite and >= 0"));
        }
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::BadOption("tolerance must be finite and > 0"));
        }
        if max_iterations == 0 {
            return Err(Error::BadOption("max_iterations must be >= 1"));
        }
        Ok(Self {
            response,
            rates,
            inv_var,
            prior: guess.to_vec(),
            lambda: vec![0.0; response.len()],
            spectrum: guess.to_vec(),
            target_chi2,
            tolerance,
            max_iterations,
            damping: INITIAL_DAMPING,
            completed: 0,
            converged: false,
            halt: None,
        })
    }

    /// Tolerance the run converges under (largest per-group relative change
    /// strictly below it, with the chi-square at or below its target).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Chi-square target the run accepts at or below.
    pub fn target_chi2(&self) -> f64 {
        self.target_chi2
    }

    /// Current Levenberg-Marquardt damping (M6): adapts per adjustment
    /// (shrinks on acceptance, grows on rejection).
    pub fn damping(&self) -> f64 {
        self.damping
    }

    /// Number of adjustments applied so far.
    pub fn iterations(&self) -> usize {
        self.completed
    }

    /// Whether the last applied adjustment met the tolerance and the target.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// Current working spectrum (the guess before the first adjustment).
    pub fn spectrum(&self) -> &[f64] {
        &self.spectrum
    }

    /// Spectrum of the exponential family (M4) at the given multipliers.
    ///
    /// The exponent is clamped to `±700` as overflow defense (the
    /// trust-region cap keeps the multipliers bounded in practice); an
    /// all-zero exponent row reproduces the prior bit-for-bit.
    fn family_spectrum(&self, lambda: &[f64]) -> Vec<f64> {
        self.prior
            .iter()
            .enumerate()
            .map(|(j, &phi0)| {
                let mut exponent = 0.0f64;
                for (i, row) in self.response.iter().enumerate() {
                    exponent -= row[j] * lambda[i];
                }
                phi0 * exponent.clamp(-700.0, 700.0).exp()
            })
            .collect()
    }

    /// One equilibrated damped dual step (M6): returns the trial multipliers,
    /// or `None` for a no-op (zero residual, or no damping level improves
    /// the chi-square). Records a loud halt and returns `None` when the
    /// curvature itself is degenerate or the shared kernel refuses the
    /// solve — the caller then yields nothing further.
    fn dual_step(
        &mut self,
        fitted: &[usize],
        curvature: &[Vec<f64>],
        residual: &[f64],
        chi2: f64,
    ) -> Option<Vec<f64>> {
        let dim = fitted.len();
        // Equilibration scales (M6): max-abs row scales over the raw
        // curvature, then max-abs column scales over the row-scaled rows,
        // so the QR route inside the shared kernel sees order-unity entries
        // no matter how many orders the spectrum spans. The row scales are
        // compensated into the least-squares weights (`p_i²/σ_i²`), so the
        // equilibration changes neither the caller weights nor the fixed
        // points (a zero scaled step still means the dual gradient `JᵀWr`
        // vanishes bit-for-bit). A zero scale is a vacuous dual direction —
        // unreachable under the strictly positive update, so loud, never
        // silent.
        let mut row_scale = vec![0.0f64; dim];
        for (a, row) in curvature.iter().enumerate() {
            row_scale[a] = row.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            if row_scale[a] == 0.0 {
                self.halt = Some(Error::Lstsq("dual curvature row is vacuous".to_string()));
                return None;
            }
        }
        let mut scaled: Vec<Vec<f64>> = Vec::with_capacity(dim);
        let mut target: Vec<f64> = Vec::with_capacity(dim);
        for (a, row) in curvature.iter().enumerate() {
            scaled.push(row.iter().map(|&v| v / row_scale[a]).collect());
            target.push(residual[a] / row_scale[a]);
        }
        let mut col_scale = vec![0.0f64; dim];
        for k in 0..dim {
            col_scale[k] = scaled.iter().map(|row| row[k].abs()).fold(0.0f64, f64::max);
            if col_scale[k] == 0.0 {
                self.halt = Some(Error::Lstsq("dual curvature column is vacuous".to_string()));
                return None;
            }
            for row in scaled.iter_mut() {
                row[k] /= col_scale[k];
            }
        }
        if residual.iter().all(|&v| v == 0.0) {
            // Already at the family optimum (an exact-fit fixed point, or
            // the minimum-chi-square member): any step is a no-op.
            return None;
        }
        // Adaptive damping loop: accept the first trial that strictly
        // lowers the chi-square (M3), with the caller weights on the data
        // rows and unit weight on the damping rows. Data weights are
        // max-normalized so the dimensionless damping starts mild
        // (near-Newton first step) on every problem scale; the
        // normalization rescales the objective only, never its minimizer.
        // Both ratios stay in (0, 1] — no overflow for extreme sigmas, no
        // division by zero (sigmas are finite positive, scales guarded
        // positive above).
        let inv_peak = fitted
            .iter()
            .map(|&i| self.inv_var[i])
            .fold(0.0f64, f64::max);
        let row_peak = row_scale.iter().cloned().fold(0.0f64, f64::max);
        let peak = fitted
            .iter()
            .enumerate()
            .map(|(a, &i)| (self.inv_var[i] / inv_peak) * (row_scale[a] / row_peak).powi(2))
            .fold(0.0f64, f64::max);
        let mut design = Vec::with_capacity(2 * dim);
        let mut augmented = Vec::with_capacity(2 * dim);
        let mut weights = Vec::with_capacity(2 * dim);
        for _ in 0..DAMPING_RETRIES {
            design.clear();
            augmented.clear();
            weights.clear();
            for (a, row) in scaled.iter().enumerate() {
                design.push(row.clone());
                augmented.push(target[a]);
                // (M6) weight compensation: dividing the equation by `p_i`
                // divides its squared residual by `p_i²`, so multiplying the
                // caller weight by `p_i²` keeps the objective — and the
                // minimum-chi-square member — exact.
                weights.push(
                    (self.inv_var[fitted[a]] / inv_peak) * (row_scale[a] / row_peak).powi(2) / peak,
                );
            }
            let root = self.damping.sqrt();
            for k in 0..dim {
                // Damping penalizes the physical multiplier step
                // `Δ_k = δ_k / q_k` (MINPACK-style D-scaling with the
                // equilibration scales): a weakly-constrained direction
                // (tiny `q_k`) may not buy a visible fit improvement with a
                // huge multiplier — huge multipliers annihilate their own
                // groups through the diagonal response and the exponential
                // family cannot resurrect them. Strong directions solve
                // essentially undamped.
                let mut row = vec![0.0; dim];
                row[k] = root / col_scale[k];
                design.push(row);
                augmented.push(0.0);
                weights.push(1.0);
            }
            let delta = match lstsq::weighted_lstsq(&design, &augmented, &weights) {
                Ok(delta) => delta,
                Err(err) => {
                    self.halt = Some(Error::Lstsq(err.to_string()));
                    return None;
                }
            };
            // Trust region (M6): halve the step toward the current
            // multipliers until the trial both strictly improves the fit
            // and moves no group by more than a factor of `e`. Shrinking
            // toward the current point always enters the region, so only a
            // trial that cannot improve at any length grows the damping;
            // overshoot alone never does.
            let mut scale = 1.0f64;
            let accepted = loop {
                let mut trial_lambda = self.lambda.clone();
                for (pos, &i) in fitted.iter().enumerate() {
                    trial_lambda[i] += scale * delta[pos] / col_scale[pos];
                }
                let trial_spectrum = self.family_spectrum(&trial_lambda);
                let trial_chi2 = chi_square(
                    &crate::fold(self.response, &trial_spectrum),
                    self.rates,
                    &self.inv_var,
                );
                let log_step = trial_spectrum
                    .iter()
                    .zip(self.spectrum.iter())
                    .map(|(new, old)| (new / old).ln().abs())
                    .fold(0.0f64, f64::max);
                if trial_chi2 < chi2 && log_step <= MAX_LOG_STEP {
                    break Some(trial_lambda);
                }
                scale *= 0.5;
                if scale < ALPHA_MIN {
                    break None;
                }
            };
            if let Some(trial_lambda) = accepted {
                self.damping = (self.damping / 5.0).max(MIN_DAMPING);
                return Some(trial_lambda);
            }
            self.damping = (self.damping * 10.0).min(MAX_DAMPING);
            if self.damping >= MAX_DAMPING {
                break;
            }
        }
        // No damping level improves the fit (already at the family optimum
        // to working precision): a no-op, whose zero relative change then
        // trips the tolerance leg.
        None
    }
}

impl Iterator for Maxed<'_> {
    type Item = Iteration;

    fn next(&mut self) -> Option<Iteration> {
        if self.converged || self.halt.is_some() || self.completed >= self.max_iterations {
            return None;
        }
        // (M1) fold of the current spectrum.
        let folds = crate::fold(self.response, &self.spectrum);
        for (i, (&c, &n)) in folds.iter().zip(self.rates.iter()).enumerate() {
            if n > 0.0 && c <= 0.0 {
                // A positive measurement whose whole support folds to zero:
                // unreachable (only a zero response row folds to zero under
                // the exponential update, already rejected up front).
                self.halt = Some(Error::RatesUnreachable { detector: i });
                return None;
            }
        }
        let chi2 = chi_square(&folds, self.rates, &self.inv_var);
        // Fitted detectors only: zero measurements carry zero weight (the
        // GRAVEL house rule — skipped, never pinned to zero). An all-zero
        // measurement set is one confirming no-op returning the guess.
        let fitted: Vec<usize> = self
            .rates
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n > 0.0).then_some(i))
            .collect();
        let new_spectrum = if fitted.is_empty() {
            self.spectrum.clone()
        } else {
            // (M5) dual curvature over the fitted detectors plus the fold
            // residual it must drive to zero.
            let dim = fitted.len();
            let mut curvature = Vec::with_capacity(dim);
            let mut residual = Vec::with_capacity(dim);
            for &i in &fitted {
                let mut row = Vec::with_capacity(dim);
                for &k in &fitted {
                    let jacobian: f64 = self.response[i]
                        .iter()
                        .zip(self.spectrum.iter())
                        .zip(self.response[k].iter())
                        .map(|((&r_ij, &phi), &r_kj)| r_ij * phi * r_kj)
                        .sum();
                    row.push(jacobian);
                }
                curvature.push(row);
                residual.push(folds[i] - self.rates[i]);
            }
            match self.dual_step(&fitted, &curvature, &residual, chi2) {
                None => {
                    if self.halt.is_some() {
                        return None;
                    }
                    self.spectrum.clone()
                }
                Some(trial_lambda) => {
                    let trial_spectrum = self.family_spectrum(&trial_lambda);
                    self.lambda = trial_lambda;
                    trial_spectrum
                }
            }
        };
        let max_rel_change = new_spectrum
            .iter()
            .zip(self.spectrum.iter())
            .map(|(new, old)| ((new - old) / old).abs())
            .fold(0.0f64, f64::max);
        self.spectrum = new_spectrum;
        self.completed += 1;
        let new_chi2 = chi_square(
            &crate::fold(self.response, &self.spectrum),
            self.rates,
            &self.inv_var,
        );
        if max_rel_change < self.tolerance && new_chi2 <= self.target_chi2 {
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

/// Converged MAXED solution with its convergence diagnostics. Same shape as
/// [`crate::sandii::Solution`].
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Adjusted spectrum (one value per energy group).
    pub spectrum: Vec<f64>,
    /// Rates folded from the adjusted spectrum (M1).
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

/// Drive a [`Maxed`] iterator to convergence and return the solution.
///
/// Non-convergence is a hard [`Error::NotConverged`] — no partial spectrum is
/// returned (the `tritium` face-Newton precedent). That covers both the
/// exhausted cap and a relatively converged spectrum whose chi-square still
/// exceeds `target_chi2`. Inputs and options are validated up front exactly
/// as by [`Maxed::new`].
pub fn unfold(
    response: &[Vec<f64>],
    rates: &[f64],
    sigmas: &[f64],
    guess: &[f64],
    target_chi2: Option<f64>,
    tolerance: f64,
    max_iterations: usize,
) -> Result<Solution> {
    let mut run = Maxed::new(
        response,
        rates,
        sigmas,
        guess,
        target_chi2,
        tolerance,
        max_iterations,
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
    fn entropy_is_maximal_at_the_prior() {
        let prior = vec![1.0, 2.0, 3.0];
        let at_prior = relative_entropy(&prior, &prior).unwrap();
        assert_eq!(at_prior, 0.0);
        let moved = vec![1.5, 1.0, 4.0];
        assert!(relative_entropy(&moved, &prior).unwrap() < 0.0);
        assert!(relative_entropy(&prior, &moved).unwrap() < 0.0);
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
            None,
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
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 10_000)
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
        // and one sloppy (2.2 ± 1.0): the `1/σ²` chi-square weights (10000
        // vs 1) must follow the precise detector, not the average.
        let response = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let rates = vec![2.0, 2.2];
        let sigmas = vec![0.01, 1.0];
        let guess = vec![1.5, 7.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 1000).unwrap();
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
        // to the truth than the guess was — the maximum-entropy payoff:
        // among exact-fit spectra the family member closest to the prior
        // wins. Sigmas follow the counting-statistics model (σ² = N).
        let n_groups = 24;
        let response = synthetic_response(6, n_groups, 3.0);
        let midpoints = log_midpoints(n_groups, 1e-6, 12.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas: Vec<f64> = rates.iter().map(|r| r.sqrt()).collect();
        let guess: Vec<f64> = truth
            .iter()
            .enumerate()
            .map(|(j, p)| p * (1.3 + 0.7 * (j % 5) as f64 / 4.0))
            .collect();
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-9, 50_000)
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
        // Group bounds bracket each shape's support so every group sees a
        // significant flux (the STAYSL-class house rule for solve-based
        // methods): with bounds extending decades past the shape's decay
        // the trailing groups sit at the 1e-30 positivity floor, and a
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
            let response = synthetic_response(n_groups, n_groups, 0.8);
            let rates = crate::forward_fold(&response, &shape).unwrap();
            let sigmas = unit_sigmas(rates.len());
            let guess: Vec<f64> = shape.iter().map(|p| 5.0 * p).collect();
            let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-11, 100_000)
                .unwrap_or_else(|e| panic!("IRDFF-II {name} probe must converge: {e}"));
            let worst = sol
                .spectrum
                .iter()
                .zip(shape.iter())
                .map(|(r, w)| ((r - w) / w).abs())
                .fold(0.0f64, f64::max);
            assert!(worst < 1e-5, "{name}: shape recovery to 1e-5, got {worst}");
        }
    }

    #[test]
    fn sandii_fixed_point_is_a_maxed_fixed_point() {
        // Both methods consume the same fold (M1 == S1): once SAND-II has
        // converged (rate factors pinned to 1), the dual residual is zero,
        // so the Gauss-Newton cycle about that spectrum as its own prior is
        // a no-op — a SAND-II fixed point is a MAXED fixed point too.
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
            None,
            1e-4,
            100,
        )
        .expect("MAXED step at a SAND-II fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(sandii_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "MAXED must hold the SAND-II fixed point"
            );
        }
    }

    #[test]
    fn gravel_fixed_point_is_a_maxed_fixed_point() {
        // GRAVEL reproduces the rates exactly (multiplicative family), so a
        // MAXED cycle about the GRAVEL spectrum as its own prior sees a zero
        // dual residual and is a no-op too.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let gravel_sol = crate::gravel::unfold(&response, &rates, &sigmas, &guess, 1e-12, 100_000)
            .expect("GRAVEL must converge");
        let sol = unfold(
            &response,
            &rates,
            &sigmas,
            &gravel_sol.spectrum,
            None,
            1e-4,
            100,
        )
        .expect("MAXED step at a GRAVEL fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(gravel_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "MAXED must hold the GRAVEL fixed point"
            );
        }
    }

    #[test]
    fn staysl_fixed_point_is_held_by_maxed() {
        // The anchored least-squares fixed point reproduces the rates up to
        // the damping slack, so the MAXED dual residual about it is small:
        // the run must settle back onto it within the shared tolerance.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let staysl_sol =
            crate::staysl::unfold(&response, &rates, &sigmas, &guess, 1e-12, 100_000, 1e-6)
                .expect("STAYSL-class must converge");
        let sol = unfold(
            &response,
            &rates,
            &sigmas,
            &staysl_sol.spectrum,
            None,
            1e-4,
            100,
        )
        .expect("MAXED must settle at a least-squares fixed point");
        let worst = sol
            .spectrum
            .iter()
            .zip(staysl_sol.spectrum.iter())
            .map(|(r, w)| ((r - w) / w).abs())
            .fold(0.0f64, f64::max);
        assert!(
            worst < 1e-3,
            "MAXED must hold the least-squares fixed point, got {worst}"
        );
    }

    #[test]
    fn maxed_fixed_point_is_a_sandii_fixed_point() {
        // MAXED reproduces the rates exactly, so a SAND-II cycle about the
        // MAXED spectrum is a no-op too.
        let response = synthetic_response(6, 6, 1.0);
        let midpoints = log_midpoints(6, 1e-6, 10.0);
        let truth = maxwellian_plus_inv_e(&midpoints, 2.53e-5);
        let rates = crate::forward_fold(&response, &truth).unwrap();
        let sigmas = unit_sigmas(rates.len());
        let guess: Vec<f64> = truth.iter().map(|p| 2.5 * p).collect();
        let maxed_sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 100_000)
            .expect("MAXED must converge");
        let sol = crate::sandii::unfold(&response, &rates, &maxed_sol.spectrum, 1e-4, 100)
            .expect("SAND-II step at a MAXED fixed point must be a no-op");
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert!(
            sol.max_rel_change < 1e-4,
            "no-op step, got change {}",
            sol.max_rel_change
        );
        for (rec, want) in sol.spectrum.iter().zip(maxed_sol.spectrum.iter()) {
            assert!(
                ((rec - want) / want).abs() < 1e-5,
                "SAND-II must hold the MAXED fixed point"
            );
        }
    }

    #[test]
    fn unconstrained_group_keeps_guess_exactly() {
        // Detector 0 responds only to group 0; groups 1..4 see nothing and
        // must keep the guess values bit-for-bit (their exponents are
        // exactly zero).
        let response = vec![vec![1.0, 0.0, 0.0, 0.0, 0.0]];
        let rates = vec![4.0];
        let sigmas = vec![0.5];
        let guess = vec![2.0, 7.0, 3.0, 9.0, 5.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 100).unwrap();
        assert!(
            (sol.spectrum[0] - 4.0).abs() < 1e-9,
            "group 0 must converge to 4.0, got {}",
            sol.spectrum[0]
        );
        assert_eq!(&sol.spectrum[1..], &guess[1..]);
        assert!(sol.iterations >= 2, "needs more than one dual step");
        assert!(sol.max_rel_change < 1e-12);
    }

    #[test]
    fn zero_measurement_detectors_carry_no_weight() {
        // Pinned house rule shared with GRAVEL: a zero measurement is
        // skipped rather than fitted, so the guess survives untouched.
        let response = vec![vec![1.0, 1.0]];
        let rates = vec![0.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 3.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum, guess);
        assert_eq!(sol.iterations, 1, "one confirming no-op adjustment");
        assert_eq!(sol.max_rel_change, 0.0);
        assert!(sol.rate_factors.iter().all(|f| f.is_finite()));
    }

    #[test]
    fn mixed_zero_measurements_fit_the_positive_subset() {
        // Detector 0 reads zero (skipped); detector 1 pins group 1 to 4.0.
        // Group 0 sees only the skipped detector and keeps the guess.
        let response = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let rates = vec![0.0, 4.0];
        let sigmas = vec![1.0, 1.0];
        let guess = vec![2.0, 7.0];
        let sol = unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 100).unwrap();
        assert_eq!(sol.spectrum[0], 2.0);
        assert!(
            (sol.spectrum[1] - 4.0).abs() < 1e-9,
            "group 1 must converge to 4.0, got {}",
            sol.spectrum[1]
        );
    }

    #[test]
    fn chi2_target_violation_is_a_hard_not_converged() {
        // Two redundant readings of one group that disagree well beyond
        // their sigmas: no member of the exponential family fits both, so
        // the chi-square bottoms out at 50 while the target demands 1. The
        // run stalls in relative change at that minimum yet is rejected
        // loudly, never returned as a partial spectrum.
        let response = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let rates = vec![2.0, 2.1];
        let sigmas = vec![0.01, 0.01];
        let guess = vec![1.5, 7.0];
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, Some(1.0), 1e-3, 200),
            Err(Error::NotConverged)
        );
    }

    #[test]
    fn huge_sigma_is_out_of_weight_range() {
        // σ² overflowing to inf would zero-weight the detector silently;
        // out-of-weight-range sigmas are loud errors instead.
        let response = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let rates = vec![2.0, 4.0];
        let sigmas = vec![1.0, 1e200];
        let guess = vec![2.0, 7.0];
        assert!(matches!(
            unfold(&response, &rates, &sigmas, &guess, None, 1e-6, 100),
            Err(Error::BadRates(_))
        ));
    }

    #[test]
    fn zero_response_row_with_nonzero_rate_is_unreachable() {
        let response = vec![vec![0.0, 0.0]];
        let rates = vec![1.0];
        let sigmas = vec![1.0];
        let guess = vec![2.0, 2.0];
        assert_eq!(
            unfold(&response, &rates, &sigmas, &guess, None, 1e-6, 100),
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
            unfold(&response, &rates, &sigmas, &guess, None, 1e-12, 1),
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
        let mut run = Maxed::new(&response, &rates, &sigmas, &guess, None, 1e-9, 1000).unwrap();
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
            Maxed::new(&[], &ok_rates, &ok_sigmas, &ok_guess, None, 1e-3, 10),
            Err(Error::BadResponse("empty response matrix"))
        ));
        assert!(matches!(
            Maxed::new(
                &[vec![1.0], vec![1.0, 2.0]],
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                None,
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
            Maxed::new(&ok_response, &[1.0], &ok_sigmas, &ok_guess, None, 1e-3, 10),
            Err(Error::BadShape {
                what: "rates",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Maxed::new(&ok_response, &ok_rates, &[1.0], &ok_guess, None, 1e-3, 10),
            Err(Error::BadShape {
                what: "sigmas",
                expected: 2,
                got: 1
            })
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &[1.0, f64::NAN],
                &ok_guess,
                None,
                1e-3,
                10
            ),
            Err(Error::BadRates("non-finite sigma"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &[1.0, 0.0],
                &ok_guess,
                None,
                1e-3,
                10
            ),
            Err(Error::BadRates("sigma must be > 0"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &[1.0, 1e-300],
                &ok_guess,
                None,
                1e-3,
                10
            ),
            Err(Error::BadRates("sigma is out of weight range"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &[1.0, 1.0],
                None,
                1e-3,
                10
            ),
            Err(Error::BadShape {
                what: "guess",
                expected: 3,
                got: 2
            })
        ));
        let zero_guess = vec![1.0, 0.0, 1.0];
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &zero_guess,
                None,
                1e-3,
                10
            ),
            Err(Error::BadGuess("entries must be strictly positive"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                Some(f64::NAN),
                1e-3,
                10
            ),
            Err(Error::BadOption("target_chi2 must be finite and >= 0"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                Some(-1.0),
                1e-3,
                10
            ),
            Err(Error::BadOption("target_chi2 must be finite and >= 0"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                None,
                0.0,
                10
            ),
            Err(Error::BadOption("tolerance must be finite and > 0"))
        ));
        assert!(matches!(
            Maxed::new(
                &ok_response,
                &ok_rates,
                &ok_sigmas,
                &ok_guess,
                None,
                1e-3,
                0
            ),
            Err(Error::BadOption("max_iterations must be >= 1"))
        ));
    }
}
