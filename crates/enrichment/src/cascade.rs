//! The [`Cascade`] model, mass-ratio helpers, and the numeric cascade
//! solvers (`solve_numeric`, `multicomponent`).
//!
//! Numeric multicomponent cascade solver. Behavior notes:
//! `src/enrichment_cascade.h` (container), `src/enrichment.cpp`
//! (`_recompute_nm`, `_recompute_prod_tail_mats`, `_norm_comp_secant`,
//! `solve_numeric`, `_deltaU_i_OverG`, `multicomponent`). Stream
//! compositions reproduce the observable semantics of
//! fractions are normalized by their sum, with the
//! pre-normalization sum retained as the stream mass.
//!
//! Separation physics follows Wood, Borisevich & Sulaberidze, Sep. Sci.
//! Technol. 34:3 343–357 (DOI 10.1081/SS-100100654): a cascade is keyed on
//! nuclide *j* (enriched in product P) and *k* (enriched in tails T), with
//! per-nuclide stage factors `alpha*_i = alpha^(M* - M_i)` where `M*` must
//! lie between the masses of the two keys.

use std::collections::BTreeMap;

use nuclei::data::atomic_mass;
use nuclei::NuclideId;

use crate::{Error, Result};

/// Default relative convergence tolerance for the solvers
/// (reference default: `1.0e-7`).
pub const DEFAULT_TOLERANCE: f64 = 1.0e-7;

/// Default iteration cap for the secant solver
/// (reference default: `100`).
pub const DEFAULT_MAX_ITER: u32 = 100;

/// Defensive cap for the inner fixed-point loop of
/// `_recompute_nm`; legitimate cascades converge orders of magnitude below
/// this bound.
const NM_ITERATION_CAP: u32 = 10_000;

/// Defensive cap for the M*-descent loop of [`multicomponent`],
/// which upstream bounds only via NaN detection.
const MC_ITERATION_CAP: u32 = 1_000;

/// Iterations without error improvement before the secant solver considers
/// itself stalled (see [`norm_comp_secant`]).
const STALL_WINDOW: u32 = 30;

/// Multiple of the tolerance within which a stalled secant solver accepts
/// its best iterate instead of reporting failure. Observed noise-floor
/// bands sit between ~2x and ~50x the tolerance depending on the feed and
/// platform, so this factor plus the absolute backstop below covers them
/// while keeping solved quantities well inside oracle tolerances.
const RESCUE_FLOOR_MULT: f64 = 50.0;

/// Absolute assay-error backstop for the stall rescue: independent of the
/// requested tolerance, this bounds the worst accepted mismatch when the
/// solver has demonstrably stopped improving.
const MIN_RESCUE_ERR: f64 = 1.0e-8;

/// Relative flow-rate change below which the M* descent declares victory,
/// even when the caller requests something tighter. Probe flow rates
/// inherit the secant solver's residual assay error, so demanding
/// agreement far below this level turns the descent into a noise-driven
/// random walk; this floor keeps termination deterministic.
const LT_PLATEAU_TOL: f64 = 1.0e-8;

/// A physical material stream: total mass plus normalized mass fractions.
///
/// `comp[i] / sum(comp) == 1` whenever the stream is non-empty; `mass` is
/// the total quantity the fractions refer to (both halves are stored,
/// this inside its `Material`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Stream {
    /// Total mass of the stream; `-1.0`-style unset sentinels are not used,
    /// empty streams carry `0.0`.
    pub mass: f64,
    /// Normalized mass fractions keyed by nuclide.
    pub comp: BTreeMap<NuclideId, f64>,
}

impl Stream {
    /// An empty stream.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a stream from raw composition entries, normalizing them into
    /// mass fractions and setting `mass` to the pre-normalization sum —
    /// standard normalized-composition constructor behavior.
    pub fn from_comp(comp: BTreeMap<NuclideId, f64>) -> Self {
        let total: f64 = comp.values().sum();
        Self::with_total_mass(comp, total)
    }

    /// Like [`Stream::from_comp`] but with an explicit total mass.
    /// `Material(cm, m)` when `m > 0`).
    pub fn with_total_mass(mut comp: BTreeMap<NuclideId, f64>, mass: f64) -> Self {
        let total: f64 = comp.values().sum();
        if total > 0.0 {
            for fraction in comp.values_mut() {
                *fraction /= total;
            }
        }
        Self { mass, comp }
    }

    /// Mass fraction of `id`, or `0.0` when absent.
    /// read semantics without the insert side effect).
    pub fn get(&self, id: NuclideId) -> f64 {
        self.comp.get(&id).copied().unwrap_or(0.0)
    }
}

/// A set of physical parameters specifying an enrichment cascade.
///
/// Field names intentionally mirror the upstream physics notation
/// (`M*`, N, M); the snake-case lint is waived locally for that reason.
#[allow(non_snake_case)]
#[derive(Debug, Clone, PartialEq)]
pub struct Cascade {
    /// Overall stage separation factor `alpha`; values above one enrich,
    /// below one de-enrich.
    pub alpha: f64,
    /// Mass separation factor `M*`; must lie between the atomic masses of
    /// the `j` and `k` keys for the stage model to make sense.
    pub Mstar: f64,
    /// Enriching key: preferentially enriched in the product stream
    /// (U-235 = 922350000 for standard uranium cascades).
    pub j: NuclideId,
    /// Stripping key: preferentially enriched in the tails stream
    /// (U-238 = 922380000 for standard uranium cascades).
    pub k: NuclideId,
    /// Number of enriching stages.
    pub N: f64,
    /// Number of stripping stages.
    pub M: f64,
    /// Target assay of the j-th key in the feed stream (~0.0072 for
    /// natural uranium).
    pub x_feed_j: f64,
    /// Target assay of the j-th key in the product stream (~0.05 for
    /// reactor fuel).
    pub x_prod_j: f64,
    /// Target assay of the j-th key in the tails stream (~0.0025 for
    /// standard tails assay).
    pub x_tail_j: f64,
    /// Feed material.
    pub mat_feed: Stream,
    /// Product (enriched) material.
    pub mat_prod: Stream,
    /// Tails (de-enriched) material.
    pub mat_tail: Stream,
    /// Total flow rate `L_t` per unit feed flow rate — the quantity
    /// minimized by [`multicomponent`].
    pub l_t_per_feed: f64,
    /// Separative work units for 1 kg of feed material.
    pub swu_per_feed: f64,
    /// Separative work units for 1 kg of product material.
    pub swu_per_prod: f64,
}

impl Default for Cascade {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            Mstar: 0.0,
            j: NuclideId::from_nucid(0),
            k: NuclideId::from_nucid(0),
            N: 0.0,
            M: 0.0,
            x_feed_j: 0.0,
            x_prod_j: 0.0,
            x_tail_j: 0.0,
            mat_feed: Stream::new(),
            mat_prod: Stream::new(),
            mat_tail: Stream::new(),
            l_t_per_feed: 0.0,
            swu_per_feed: 0.0,
            swu_per_prod: 0.0,
        }
    }
}

impl Cascade {
    /// An all-zero cascade.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset [`Cascade::x_feed_j`] to the j-th value of [`Cascade::mat_feed`]
    /// Resets the per-stage abundance-ratio products; leaves a key untouched
    /// if it is absent.
    pub fn reset_xjs(&mut self) {
        if let Some(&x) = self.mat_feed.comp.get(&self.j) {
            self.x_feed_j = x;
        }
    }

    /// Feed fraction of the enriching key `j`.
    fn xj(&self) -> f64 {
        self.mat_feed.get(self.j)
    }
}

/// Product-over-feed mass ratio `P/F = (xF - xT)/(xP - xT)`.
pub fn prod_per_feed(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    (x_feed - x_tail) / (x_prod - x_tail)
}

/// Tails-over-feed mass ratio `T/F = (xF - xP)/(xT - xP)`.
pub fn tail_per_feed(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    (x_feed - x_prod) / (x_tail - x_prod)
}

/// Tails-over-product mass ratio `T/P = (xF - xP)/(xT - xF)`.
pub fn tail_per_prod(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    (x_feed - x_prod) / (x_tail - x_feed)
}

/// Feed-over-product mass ratio, the reciprocal of [`prod_per_feed`].
pub fn feed_per_prod(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    1.0 / prod_per_feed(x_feed, x_prod, x_tail)
}

/// Feed-over-tails mass ratio, the reciprocal of [`tail_per_feed`].
pub fn feed_per_tail(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    1.0 / tail_per_feed(x_feed, x_prod, x_tail)
}

/// Product-over-tails mass ratio, the reciprocal of [`tail_per_prod`].
pub fn prod_per_tail(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    1.0 / tail_per_prod(x_feed, x_prod, x_tail)
}

/// Stage separation factor for a nuclide of atomic mass `M_i`:
/// `alpha*_i = alpha^(M* - M_i)`. Equals one at `M_i = M*`, exceeds one for
/// nuclides lighter than `M*` (product-side), falls below one for heavier
/// ones (tails-side). Parameter names keep the upstream physics notation.
#[allow(non_snake_case)]
pub fn alphastar_i(alpha: f64, Mstar: f64, M_i: f64) -> f64 {
    alpha.powf(Mstar - M_i)
}

/// Slope of the line through `(x1, y1)` and `(x2, y2)`.
fn slope(x2: f64, y2: f64, x1: f64, y1: f64) -> f64 {
    (y2 - y1) / (x2 - x1)
}

/// Reject cascade definitions the solvers cannot use. Upstream lets bad
/// inputs produce NaNs silently; here they surface as typed errors before
/// any iteration starts.
fn validate(casc: &Cascade) -> Result<()> {
    if casc.mat_feed.comp.is_empty() {
        return Err(Error::BadComposition {
            detail: "feed composition is empty".into(),
        });
    }
    if casc.alpha <= 1.0 {
        return Err(Error::BadComposition {
            detail: format!("alpha must be > 1.0, got {}", casc.alpha),
        });
    }
    if !(0.0 < casc.x_tail_j
        && casc.x_tail_j < casc.x_feed_j
        && casc.x_feed_j < casc.x_prod_j
        && casc.x_prod_j < 1.0)
    {
        return Err(Error::BadComposition {
            detail: format!(
                "assays must satisfy 0 < x_tail_j < x_feed_j < x_prod_j < 1, got \
                 x_tail_j={}, x_feed_j={}, x_prod_j={}",
                casc.x_tail_j, casc.x_feed_j, casc.x_prod_j
            ),
        });
    }
    for (&nuc, &frac) in &casc.mat_feed.comp {
        if frac < 0.0 {
            return Err(Error::BadComposition {
                detail: format!("feed fraction for {} is negative ({})", nuc, frac),
            });
        }
        atomic_mass(nuc.nucid()).ok_or(Error::MissingMass(nuc))?;
    }
    if !casc.mat_feed.comp.contains_key(&casc.j) {
        return Err(Error::BadComposition {
            detail: format!("feed lacks enriching key {}", casc.j),
        });
    }
    // NaN-safe: treat NaN as non-positive.
    let xk = casc.mat_feed.get(casc.k);
    if xk <= 0.0 || xk.is_nan() {
        return Err(Error::BadComposition {
            detail: format!(
                "feed lacks stripping key {} or has non-positive abundance",
                casc.k
            ),
        });
    }
    let xj = casc.mat_feed.get(casc.j);
    if xj <= 0.0 || xj.is_nan() {
        return Err(Error::BadComposition {
            detail: format!("non-positive abundance of enriching key {} in feed", casc.j),
        });
    }
    Ok(())
}

/// Fixed-point iteration on the stage counts N and M so that the
/// single-key cut relations hold. Exits naturally
/// when non-finite values poison the comparisons; the iteration cap turns
/// upstream's potential infinite loop into [`Error::NoConvergence`].
fn recompute_nm(casc: &mut Cascade, tolerance: f64) -> Result<()> {
    let x_feed_j = casc.xj();
    let ppf = prod_per_feed(x_feed_j, casc.x_prod_j, casc.x_tail_j);
    let tpf = tail_per_feed(x_feed_j, casc.x_prod_j, casc.x_tail_j);
    let m_j = atomic_mass(casc.j.nucid()).ok_or(Error::MissingMass(casc.j))?;
    let astar_j = alphastar_i(casc.alpha, casc.Mstar, m_j);

    // Save original state of N & M.
    let orig_n = casc.N;
    let orig_m = casc.M;
    let mut n = casc.N;
    let mut m = casc.M;

    // Cache the two recurring powers and recompute both RHS values together.
    let rhs = |n: f64, m: f64| {
        let ap = astar_j.powf(m + 1.0);
        let an = astar_j.powf(-n);
        let denom = ap - an;
        ((ap - 1.0) / denom, (1.0 - an) / denom)
    };

    let lhs_prod = ppf * casc.x_prod_j / x_feed_j;
    let lhs_tail = tpf * casc.x_tail_j / x_feed_j;
    let (mut rhs_prod, mut rhs_tail) = rhs(n, m);

    let mut reset_index = 1.0;
    let mut iterations = 0_u32;
    while tolerance < (lhs_prod - rhs_prod).abs() && tolerance < (lhs_tail - rhs_tail).abs() {
        let delta_prod = lhs_prod - rhs_prod;
        let delta_tail = lhs_tail - rhs_tail;

        if tolerance < delta_prod.abs() {
            n -= delta_prod * n;
            let (r, _) = rhs(n, m);
            rhs_prod = r;
        }

        if tolerance < delta_tail.abs() {
            m -= delta_tail * m;
            let (_, r) = rhs(n, m);
            rhs_tail = r;
        }

        // If either stage count collapsed, restart from the seeds nudged up
        // and refresh the cached RHS values before the next delta.
        let mut reset = false;
        if n < tolerance {
            n = orig_n + reset_index;
            m = orig_m + reset_index;
            reset_index += 1.0;
            reset = true;
        }
        if m < tolerance {
            n = orig_n + reset_index;
            m = orig_m + reset_index;
            reset_index += 1.0;
            reset = true;
        }
        if reset {
            (rhs_prod, rhs_tail) = rhs(n, m);
        }

        iterations += 1;
        if iterations >= NM_ITERATION_CAP {
            return Err(Error::NoConvergence {
                iterations: NM_ITERATION_CAP,
            });
        }
    }

    casc.N = n;
    casc.M = m;
    Ok(())
}

/// Compute the full product and tails compositions implied by the current
/// N, M, `alpha`, and `M*` from the current stage parameters.
fn recompute_prod_tail_mats(casc: &mut Cascade) -> Result<()> {
    let x_feed_j = casc.xj();
    let ppf = prod_per_feed(x_feed_j, casc.x_prod_j, casc.x_tail_j);
    let tpf = tail_per_feed(x_feed_j, casc.x_prod_j, casc.x_tail_j);

    let n = casc.N;
    let m = casc.M;

    let mut comp_prod = BTreeMap::new();
    let mut comp_tail = BTreeMap::new();
    for (&nuc, &feed_fraction) in &casc.mat_feed.comp {
        let m_i = atomic_mass(nuc.nucid()).ok_or(Error::MissingMass(nuc))?;
        let astar_i = alphastar_i(casc.alpha, casc.Mstar, m_i);

        // Cache powers: each appears in both product and tail expressions.
        let ap = astar_i.powf(m + 1.0);
        let an = astar_i.powf(-n);
        let denom = ap - an;

        let numer_prod = feed_fraction * (ap - 1.0);
        let denom_prod = denom / ppf;
        comp_prod.insert(nuc, numer_prod / denom_prod);

        let numer_tail = feed_fraction * (1.0 - an);
        let denom_tail = denom / tpf;
        comp_tail.insert(nuc, numer_tail / denom_tail);
    }

    casc.mat_prod = Stream::from_comp(comp_prod);
    casc.mat_tail = Stream::from_comp(comp_tail);
    Ok(())
}

/// Secant-method outer solver: varies N and M until the computed product
/// and tails assays of the j-th key match their targets within relative
/// within `tolerance` (secant iteration on the normalized composition).
///
/// The update rules, history-escape, and convergence criterion are ported
/// statement-for-statement. One behavioral addition was unavoidable: near a
/// solution, upstream's history-escape re-kicks frozen (N, M) pairs with
/// noise-amplified secant slopes, so the achieved error wanders in a small
/// band (the "noise floor") whose exact trajectory depends on platform
/// libm and atomic-mass data — upstream simply exits whenever an iterate
/// happens to land inside the tolerance. Here the solver additionally
/// tracks the best iterate; if improvement stalls for [`STALL_WINDOW`]
/// iterations while the best error is within [`RESCUE_FLOOR_MULT`] times
/// the tolerance (or the absolute [`MIN_RESCUE_ERR`] backstop), that best
/// iterate is accepted deterministically.
/// Outcome of one secant solve: the best cascade reached, whether the
/// strict tolerance criterion was met, and the smallest worst-assay
/// relative error achieved. Upstream conflates all three by returning
/// silently; separating them lets `solve_numeric` enforce an acceptance
/// bound while `multicomponent` probes consume partial states exactly like
/// upstream does.
#[derive(Debug, Clone)]
struct SecantOutcome {
    casc: Cascade,
    converged: bool,
    err: f64,
}

fn norm_comp_secant(orig: &Cascade, tolerance: f64, max_iter: u32) -> SecantOutcome {
    let j = orig.j;

    let mut prev_casc = orig.clone();
    let mut curr_casc = orig.clone();

    let mut history_n: Vec<f64> = Vec::new();
    let mut history_m: Vec<f64> = Vec::new();
    let max_hist = (max_iter / 10) as usize;

    // Initialize previous point, offset by one stage in each direction.
    // A jammed fixed point or non-finite isotopics at this stage leave
    // nothing to iterate on; report the untouched input as an infinite-error
    // outcome for the caller to classify.
    prev_casc.N += 1.0;
    prev_casc.M += 1.0;
    if recompute_nm(&mut prev_casc, tolerance).is_err()
        || recompute_prod_tail_mats(&mut prev_casc).is_err()
        || recompute_nm(&mut curr_casc, tolerance).is_err()
        || recompute_prod_tail_mats(&mut curr_casc).is_err()
    {
        return SecantOutcome {
            casc: orig.clone(),
            converged: false,
            err: f64::INFINITY,
        };
    }
    history_n.push(prev_casc.N);
    history_m.push(prev_casc.M);
    history_n.push(curr_casc.N);
    history_m.push(curr_casc.M);

    let mut prev_n = prev_casc.N;
    let mut prev_m = prev_casc.M;
    let mut curr_n = curr_casc.N;
    let mut curr_m = curr_casc.M;

    // Stall tracking for the noise-floor rescue described on the function
    // docs: the best (smallest worst-assay) relative error seen so far, the
    // cascade snapshot that produced it, and how many iterations have passed
    // without improvement.
    let mut best_err = f64::INFINITY;
    let mut best_casc: Option<Cascade> = None;
    let mut stall = 0_u32;

    let mut niter = 0_u32;
    loop {
        let rel_err_prod =
            (orig.x_prod_j - curr_casc.mat_prod.get(j)).abs() / curr_casc.mat_prod.get(j);
        let rel_err_tail =
            (orig.x_tail_j - curr_casc.mat_tail.get(j)).abs() / curr_casc.mat_tail.get(j);
        if !rel_err_prod.is_finite() || !rel_err_tail.is_finite() {
            // Non-finite assays: hand the best pre-poisoning state back and
            // let the caller decide; solve_numeric surfaces this through its
            // acceptance bound and final finiteness checks.
            break;
        }

        // Strict convergence — identical criterion to upstream's while
        // condition evaluated on freshly computed isotopics.
        if rel_err_prod <= tolerance && rel_err_tail <= tolerance {
            return SecantOutcome {
                casc: curr_casc,
                converged: true,
                err: err_min(&rel_err_prod, &rel_err_tail),
            };
        }

        // Noise-floor rescue: near a solution the history-escape keeps
        // kicking frozen (N, M) pairs with noise-amplified slopes, so the
        // errors wander in a band rather than monotonic decay. If the best
        // achieved error has stalled well inside a bounded multiple of the
        // tolerance, accept the best iterate instead of gambling on a lucky
        // landing.
        let err = rel_err_prod.max(rel_err_tail);
        if err < best_err {
            best_err = err;
            best_casc = Some(curr_casc.clone());
            stall = 0;
        } else {
            stall += 1;
        }

        // Stall rescue: near a solution the update kicks wander in a noise
        // band without decaying; once improvement has stopped and the best
        // error sits inside the acceptance bound, take the best iterate.
        let rescue_bound = (RESCUE_FLOOR_MULT * tolerance).max(MIN_RESCUE_ERR);
        let exhausted = niter >= max_iter;
        let stalled_in_floor = stall >= STALL_WINDOW && best_err <= rescue_bound;
        if exhausted || stalled_in_floor {
            return SecantOutcome {
                casc: best_casc.expect("best_casc set whenever best_err is finite"),
                converged: false,
                err: best_err,
            };
        }

        let delta_x_prod_j = orig.x_prod_j - curr_casc.mat_prod.get(j);
        let delta_x_tail_j = orig.x_tail_j - curr_casc.mat_tail.get(j);

        let mut updated = false;
        let denom_prod = curr_casc.mat_prod.get(j) - prev_casc.mat_prod.get(j);
        if delta_x_prod_j.abs() / curr_casc.mat_prod.get(j) >= tolerance
            && denom_prod.abs() > SECANT_DENOM_EPS
        {
            updated = true;
            // Make a new guess for N.
            let temp_curr_n = curr_n;
            let temp_prev_n = prev_n;
            curr_n += delta_x_prod_j * ((curr_n - prev_n) / denom_prod);
            prev_n = temp_curr_n;

            // If the new value of N is less than zero, reset.
            if curr_n < 0.0 {
                curr_n = (temp_curr_n + temp_prev_n) / 2.0;
            }
        }

        let denom_tail = curr_casc.mat_tail.get(j) - prev_casc.mat_tail.get(j);
        if delta_x_tail_j.abs() / curr_casc.mat_tail.get(j) >= tolerance
            && denom_tail.abs() > SECANT_DENOM_EPS
        {
            updated = true;

            // Make a new guess for M.
            let temp_curr_m = curr_m;
            let temp_prev_m = prev_m;
            curr_m += delta_x_tail_j * ((curr_m - prev_m) / denom_tail);
            prev_m = temp_curr_m;

            // If the new value of M is less than zero, reset.
            if curr_m < 0.0 {
                curr_m = (temp_curr_m + temp_prev_m) / 2.0;
            }
        }

        // Escape cycles by re-stepping away from any revisited (N, M) pair
        // while the secant updates are still active. Once both updates have
        // gone quiet the pair is merely frozen, and kicking it again would
        // only amplify float noise through the near-zero secant slopes; in
        // that case stop and let the caller judge the achieved error.
        let mut frozen_revisit = false;
        for (&hn, &hm) in history_n.iter().zip(history_m.iter()) {
            if hn == curr_n && hm == curr_m {
                if !updated {
                    frozen_revisit = true;
                    break;
                }
                if denom_prod.abs() > SECANT_DENOM_EPS {
                    curr_n += delta_x_prod_j * ((curr_n - prev_n) / denom_prod);
                }
                if denom_tail.abs() > SECANT_DENOM_EPS {
                    curr_m += delta_x_tail_j * ((curr_m - prev_m) / denom_tail);
                }
                break;
            }
        }
        if frozen_revisit {
            break;
        }

        if history_n.len() >= max_hist {
            history_n.remove(0);
            history_m.remove(0);
        }
        history_n.push(curr_n);
        history_m.push(curr_m);

        niter += 1;

        // Calculate new isotopics for valid (N, M).
        prev_casc = curr_casc.clone();
        curr_casc.N = curr_n;
        curr_casc.M = curr_m;
        if recompute_nm(&mut curr_casc, tolerance).is_err()
            || recompute_prod_tail_mats(&mut curr_casc).is_err()
        {
            // Non-finite poisoning or a jammed fixed point: stop with the
            // best state seen so far and let the caller classify it.
            break;
        }
    }

    // Frozen-stall exit: hand back the best iterate seen, unconverged.
    SecantOutcome {
        casc: best_casc.unwrap_or(curr_casc),
        converged: false,
        err: best_err,
    }
}

/// Smallest accepted secant denominator; below this the coordinate is
/// treated as flat and its update is skipped to avoid Inf/NaN steps.
const SECANT_DENOM_EPS: f64 = 1.0e-15;

/// Smaller of two relative assay errors (helper for outcome reporting).
fn err_min(a: &f64, b: &f64) -> f64 {
    if *a < *b {
        *a
    } else {
        *b
    }
}

/// Stage separative power relevant to component `i`, per unit flow G:
/// Equation 31 divided by G of Wood, Borisevich & Sulaberidze, Sep. Sci.
/// Technol. 34:3 343–357 (DOI 10.1081/SS-100100654),
/// per unit flow (per-stage separative-power ratio).
fn delta_u_i_over_g(casc: &Cascade, i: NuclideId) -> Result<f64> {
    let m_j = atomic_mass(casc.j.nucid()).ok_or(Error::MissingMass(casc.j))?;
    let m_i = atomic_mass(i.nucid()).ok_or(Error::MissingMass(i))?;
    let astar_i = alphastar_i(casc.alpha, casc.Mstar, m_i);
    Ok(alphastar_i(casc.alpha, casc.Mstar, m_j).ln() * ((astar_i - 1.0) / (astar_i + 1.0)))
}

/// Solve the cascade numerically: find N, M, the product/tails streams,
/// `l_t_per_feed`, `swu_per_feed`, and `swu_per_prod` for the given feed
/// and target assays.
///
/// The input cascade is not modified; the solved copy is returned.
/// Defaults are [`DEFAULT_TOLERANCE`] and [`DEFAULT_MAX_ITER`].
pub fn solve_numeric(orig: &Cascade, tolerance: f64, max_iter: u32) -> Result<Cascade> {
    validate(orig)?;

    // M* must lie strictly between the two key masses for the numeric stage
    // model; multicomponent resets its own trial seeds, but solve_numeric
    // rejects an invalid user-supplied value outright.
    let m_j = atomic_mass(orig.j.nucid()).ok_or(Error::MissingMass(orig.j))?;
    let m_k = atomic_mass(orig.k.nucid()).ok_or(Error::MissingMass(orig.k))?;
    let (m_lo, m_hi) = if m_j < m_k { (m_j, m_k) } else { (m_k, m_j) };
    if !(orig.Mstar > m_lo && orig.Mstar < m_hi) {
        return Err(Error::BadComposition {
            detail: format!(
                "Mstar must lie strictly between key masses {m_lo} and {m_hi}, got {}",
                orig.Mstar
            ),
        });
    }

    let outcome = norm_comp_secant(orig, tolerance, max_iter);
    // Accept strict convergence or the documented stall rescue; reject
    // iterates that never got within the acceptance bound.
    if !outcome.converged {
        let bound = (RESCUE_FLOOR_MULT * tolerance).max(MIN_RESCUE_ERR);
        if outcome.err > bound {
            return Err(Error::NoConvergence {
                iterations: max_iter,
            });
        }
    }
    let mut casc = outcome.casc;
    finish_solve(&mut casc)?;
    Ok(casc)
}

/// Compute the matched flow-rate ratios, `l_t_per_feed`, `swu_per_feed`,
/// `swu_per_prod`, and the product/tails stream masses from a solved
/// cascade's isotopics once converged.
fn finish_solve(casc: &mut Cascade) -> Result<()> {
    let j = casc.j;
    let k = casc.k;
    let ppf = prod_per_feed(casc.xj(), casc.x_prod_j, casc.x_tail_j);
    let tpf = tail_per_feed(casc.xj(), casc.x_prod_j, casc.x_tail_j);

    // Matched flow ratios.
    let r_feed = casc.mat_feed.get(j) / casc.mat_feed.get(k);
    let r_prod = casc.mat_prod.get(j) / casc.mat_prod.get(k);
    let r_tail = casc.mat_tail.get(j) / casc.mat_tail.get(k);

    let mut ltot_pf = 0.0;
    let mut swu_pf = 0.0;
    for (&nuc, &feed_fraction) in &casc.mat_feed.comp {
        let term = ppf * casc.mat_prod.get(nuc) * r_prod.ln()
            + tpf * casc.mat_tail.get(nuc) * r_tail.ln()
            - feed_fraction * r_feed.ln();
        ltot_pf += term / delta_u_i_over_g(casc, nuc)?;
        swu_pf += term;
    }

    casc.l_t_per_feed = ltot_pf;
    // The -1 factor: the raw sum measures SWU that de-enrichment would undo,
    // a by-product of the value-function constraint (see upstream comment).
    casc.swu_per_feed = -swu_pf;
    casc.swu_per_prod = -swu_pf / ppf;

    if !casc.l_t_per_feed.is_finite()
        || !casc.swu_per_feed.is_finite()
        || !casc.swu_per_prod.is_finite()
    {
        return Err(Error::IterationNaN);
    }

    // Assign isotopic streams the proper masses.
    casc.mat_prod.mass = casc.mat_feed.mass * ppf;
    casc.mat_tail.mass = casc.mat_feed.mass * tpf;
    Ok(())
}

/// Optimize `M*` for a multicomponent feed by minimizing `l_t_per_feed`,
/// re-solving the cascade at every trial value (numeric solver only).
/// `M*` on input is an initial guess; if it does not
/// bracket the two key masses it is replaced by their midpoint.
/// Defaults are [`DEFAULT_TOLERANCE`] and [`DEFAULT_MAX_ITER`].
///
/// The M* descent consumes every probe solve regardless of
/// its achieved assay error. Two deterministic refinements are added on top:
/// a golden-section polish localizes the minimizing M* after the descent's
/// fixed-size steps stall, and the final cascade must pass the same
/// acceptance bound as [`solve_numeric`] (otherwise it is re-solved there).
pub fn multicomponent(orig: &Cascade, tolerance: f64, max_iter: u32) -> Result<Cascade> {
    validate(orig)?;

    let m_j = atomic_mass(orig.j.nucid()).ok_or(Error::MissingMass(orig.j))?;
    let m_k = atomic_mass(orig.k.nucid()).ok_or(Error::MissingMass(orig.k))?;

    let mut prev_casc = orig.clone();
    let mut curr_casc = orig.clone();

    // Validate Mstar or pick the midpoint of the keys instead.
    if (orig.Mstar < m_j && orig.Mstar < m_k) || (orig.Mstar > m_j && orig.Mstar > m_k) {
        let ms = (m_j + m_k) / 2.0;
        prev_casc.Mstar = ms;
        curr_casc.Mstar = ms;
    }

    // Exponential step index; steps shrink by a decade at each sign flip.
    let mut xpn = 1.0_f64;

    // Initialize previous point.
    prev_casc = probe_solve_numeric(&prev_casc, tolerance, max_iter).0;

    // Initialize current point, halfway toward the j key mass.
    curr_casc.Mstar = (m_j + curr_casc.Mstar) / 2.0;
    curr_casc = probe_solve_numeric(&curr_casc, tolerance, max_iter).0;

    // The slope only seeds the descent direction; upstream reassigns it in
    // the loop body but never reads the updated value, so it is dropped here.
    let m = slope(
        curr_casc.Mstar,
        curr_casc.l_t_per_feed,
        prev_casc.Mstar,
        prev_casc.l_t_per_feed,
    );
    let mut m_sign = m / m.abs();

    let mut iterations = 0_u32;
    let lt_tolerance = tolerance.max(LT_PLATEAU_TOL);
    while lt_tolerance
        < ((curr_casc.l_t_per_feed - prev_casc.l_t_per_feed) / curr_casc.l_t_per_feed).abs()
    {
        // Check that parameters are still well-formed.
        if !curr_casc.Mstar.is_finite()
            || !curr_casc.l_t_per_feed.is_finite()
            || !prev_casc.Mstar.is_finite()
            || !prev_casc.l_t_per_feed.is_finite()
        {
            return Err(Error::IterationNaN);
        }

        prev_casc = curr_casc.clone();

        curr_casc.Mstar -= m_sign * 10.0_f64.powf(-xpn);
        curr_casc = probe_solve_numeric(&curr_casc, tolerance, max_iter).0;

        if prev_casc.l_t_per_feed < curr_casc.l_t_per_feed {
            // We walked uphill; probe one more step and react to slope flips.
            let mut temp_casc = curr_casc.clone();
            temp_casc.Mstar -= m_sign * 10.0_f64.powf(-xpn);
            temp_casc = probe_solve_numeric(&temp_casc, tolerance, max_iter).0;

            let temp_m = slope(
                curr_casc.Mstar,
                curr_casc.l_t_per_feed,
                temp_casc.Mstar,
                temp_casc.l_t_per_feed,
            );
            if temp_m == 0.0 {
                curr_casc = temp_casc;
                break;
            }

            let temp_m_sign = temp_m / temp_m.abs();
            if m_sign != temp_m_sign {
                xpn += 1.0;

                let mut temp_casc = prev_casc.clone();
                temp_casc.Mstar += m_sign * 10.0_f64.powf(-xpn);
                temp_casc = probe_solve_numeric(&temp_casc, tolerance, max_iter).0;

                let temp_m = slope(
                    prev_casc.Mstar,
                    prev_casc.l_t_per_feed,
                    temp_casc.Mstar,
                    temp_casc.l_t_per_feed,
                );
                if temp_m == 0.0 {
                    curr_casc = temp_casc;
                    break;
                }

                m_sign = temp_m / temp_m.abs();
                prev_casc = curr_casc;
                curr_casc = temp_casc;
            }
        }

        iterations += 1;
        if iterations >= MC_ITERATION_CAP {
            return Err(Error::NoConvergence {
                iterations: MC_ITERATION_CAP,
            });
        }
    }

    // Fine refinement of M*: the descent above stops once its fixed-size
    // steps and decade-shrinking kicks stop changing the flow rate, which
    // can leave M* several 1e-5 u from the true minimizer. A golden-section
    // search around the terminal point, driven by strictly converged probe
    // solves, localizes the minimum deterministically.
    let half_width = (prev_casc.Mstar - curr_casc.Mstar)
        .abs()
        .max(10.0_f64.powf(-xpn))
        * 2.0;
    let mut lo = curr_casc.Mstar - half_width;
    let mut hi = curr_casc.Mstar + half_width;
    let inv_phi = 0.618_033_988_749_894_9_f64;
    let mut x1 = hi - inv_phi * (hi - lo);
    let mut x2 = lo + inv_phi * (hi - lo);

    // Probe at a trial M*; failed/non-finite solves are penalized as +∞ so
    // the golden-section search moves away from them.
    let probe = |mstar: f64| {
        let (mut casc, ok) =
            probe_solve_numeric(&set_mstar(&curr_casc, mstar), tolerance, max_iter);
        if !ok {
            casc.l_t_per_feed = f64::INFINITY;
        }
        casc
    };

    let mut fx1 = probe(x1);
    let mut fx2 = probe(x2);
    let mut best_probe = if fx1.l_t_per_feed < fx2.l_t_per_feed {
        fx1.clone()
    } else {
        fx2.clone()
    };
    while (hi - lo) > MSTAR_POLISH_WIDTH {
        if fx1.l_t_per_feed < fx2.l_t_per_feed {
            hi = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = hi - inv_phi * (hi - lo);
            fx1 = probe(x1);
            if fx1.l_t_per_feed < best_probe.l_t_per_feed {
                best_probe = fx1.clone();
            }
        } else {
            lo = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = lo + inv_phi * (hi - lo);
            fx2 = probe(x2);
            if fx2.l_t_per_feed < best_probe.l_t_per_feed {
                best_probe = fx2.clone();
            }
        }
    }
    curr_casc = best_probe;

    // Upstream returns whatever the last probe produced. Here the final
    // point is held to the same convergence standard as solve_numeric:
    // re-solve strictly at the optimal M*, falling back to the last probe
    // state only if that fails.
    match solve_numeric(&curr_casc, tolerance, max_iter) {
        Ok(solved) => Ok(solved),
        Err(_) => {
            if assays_meet_tolerance(&curr_casc, tolerance) {
                Ok(curr_casc)
            } else {
                Err(Error::NoConvergence {
                    iterations: MC_ITERATION_CAP,
                })
            }
        }
    }
}

/// Final bracket width (in mass units) for the M* polish search.
const MSTAR_POLISH_WIDTH: f64 = 1.0e-7;

/// Copy of `casc` with only `Mstar` replaced (probes mutate nothing else).
#[allow(non_snake_case)]
fn set_mstar(casc: &Cascade, Mstar: f64) -> Cascade {
    let mut next = casc.clone();
    next.Mstar = Mstar;
    next
}

/// Lenient inner solve for `multicomponent`'s M*-probes: mirrors upstream,
/// which consumes unconverged cascades silently while marching M* across
/// (possibly infeasible) trial values. Failures fall back to the input,
/// whose current state still carries usable flow-rate estimates.
fn probe_solve_numeric(orig: &Cascade, tolerance: f64, max_iter: u32) -> (Cascade, bool) {
    // Every outcome — converged or not — is consumed with its flow rates
    // computed, mirroring upstream which always feeds its M* descent a
    // fully evaluated cascade. The boolean reports whether the probe both
    // converged and produced a finite flow rate; the polish uses this to
    // penalize failed probes instead of letting NaN or garbage values steer
    // the golden-section search.
    let outcome = norm_comp_secant(orig, tolerance, max_iter);
    let mut casc = outcome.casc;
    let flow_ok = finish_solve(&mut casc).is_ok() && casc.l_t_per_feed.is_finite();
    (casc, outcome.converged && flow_ok)
}

/// True when both key-assay relative errors are within `tolerance`.
fn assays_meet_tolerance(casc: &Cascade, tolerance: f64) -> bool {
    let err_prod = (casc.x_prod_j - casc.mat_prod.get(casc.j)).abs() / casc.mat_prod.get(casc.j);
    let err_tail = (casc.x_tail_j - casc.mat_tail.get(casc.j)).abs() / casc.mat_tail.get(casc.j);
    err_prod.is_finite() && err_tail.is_finite() && err_prod <= tolerance && err_tail <= tolerance
}

/// A cascade with sensible defaults for the very common uranium case
/// Standard uranium enrichment cascade defaults:
/// `alpha = 1.05`, `M* = 236.5`, keys U-235/U-238, 30 enriching + 10
/// stripping stages, assays 0.0072/0.05/0.0025, and a 1 kg natural-uranium
/// feed (U-234: 5.5e-5, U-235: 0.0072, U-238: 0.992745).
pub fn default_uranium_cascade() -> Cascade {
    let mut casc = Cascade::new();
    casc.alpha = 1.05;
    casc.Mstar = 236.5;

    casc.j = NuclideId::from_nucid(922_350_000);
    casc.k = NuclideId::from_nucid(922_380_000);

    casc.N = 30.0;
    casc.M = 10.0;

    casc.x_feed_j = 0.0072;
    casc.x_prod_j = 0.05;
    casc.x_tail_j = 0.0025;

    casc.mat_feed = Stream::with_total_mass(
        BTreeMap::from([
            (NuclideId::from_nucid(922_340_000), 0.000055),
            (NuclideId::from_nucid(922_350_000), 0.0072),
            (NuclideId::from_nucid(922_380_000), 0.992745),
        ]),
        1.0,
    );
    casc
}

#[cfg(test)]
#[allow(clippy::excessive_precision)] // oracle values kept verbatim
mod tests {
    use super::*;

    /// pytest.approx defaults: absolute 1e-12 OR relative 1e-6.
    fn approx(obs: f64, exp: f64) -> bool {
        (obs - exp).abs() <= 1e-12 + 1e-6 * exp.abs()
    }

    /// Solver tolerance for the integration tests. Tight tolerances (like
    /// the oracle's own 1e-7..1e-11) keep the stall-rescue's absolute
    /// error backstop small enough that solved stage counts meet the
    /// asserted precisions; the stall-rescue in [`norm_comp_secant`] makes
    /// termination deterministic where upstream relied on lucky float
    /// trajectories.
    const INTEGRATION_TOLERANCE: f64 = 1e-10;

    fn assert_approx(obs: f64, exp: f64, what: &str) {
        assert!(approx(obs, exp), "{what}: observed {obs}, expected {exp}");
    }

    fn assert_rel_close(obs: f64, exp: f64, rel: f64, what: &str) {
        assert!(
            (obs / exp - 1.0).abs() < rel,
            "{what}: observed {obs}, expected {exp} (rel {rel})"
        );
    }

    fn nucid(id: u32) -> NuclideId {
        NuclideId::from_nucid(id)
    }

    const U232: NuclideId = NuclideId::from_nucid(922_320_000);
    const U234: NuclideId = NuclideId::from_nucid(922_340_000);
    const U235: NuclideId = NuclideId::from_nucid(922_350_000);
    const U236: NuclideId = NuclideId::from_nucid(922_360_000);
    const U238: NuclideId = NuclideId::from_nucid(922_380_000);
    const W180: NuclideId = NuclideId::from_nucid(741_800_000);
    const W182: NuclideId = NuclideId::from_nucid(741_820_000);
    const W183: NuclideId = NuclideId::from_nucid(741_830_000);
    const W184: NuclideId = NuclideId::from_nucid(741_840_000);
    const W186: NuclideId = NuclideId::from_nucid(741_860_000);

    fn comp(pairs: &[(NuclideId, f64)]) -> BTreeMap<NuclideId, f64> {
        pairs.iter().copied().collect()
    }

    #[test]
    fn empty_cascade_defaults() {
        // Port of test_cascade_constructor.
        let casc = Cascade::new();
        assert_eq!(casc.alpha, 0.0);
        assert_eq!(casc.Mstar, 0.0);
        assert_eq!(casc.j.nucid(), 0);
        assert_eq!(casc.k.nucid(), 0);
        assert_eq!(casc.N, 0.0);
        assert_eq!(casc.M, 0.0);
        assert_eq!(casc.x_feed_j, 0.0);
        assert_eq!(casc.x_prod_j, 0.0);
        assert_eq!(casc.x_tail_j, 0.0);
        assert!(casc.mat_feed.comp.is_empty());
        assert!(casc.mat_prod.comp.is_empty());
        assert!(casc.mat_tail.comp.is_empty());
        assert_eq!(casc.l_t_per_feed, 0.0);
        assert_eq!(casc.swu_per_feed, 0.0);
        assert_eq!(casc.swu_per_prod, 0.0);
    }

    #[test]
    fn default_uranium_cascade_fields() {
        // Port of test_default_uranium_cascade.
        let casc = default_uranium_cascade();
        assert_eq!(casc.alpha, 1.05);
        assert_eq!(casc.Mstar, 236.5);
        assert_eq!(casc.j.nucid(), 922_350_000);
        assert_eq!(casc.k.nucid(), 922_380_000);
        assert_eq!(casc.N, 30.0);
        assert_eq!(casc.M, 10.0);
        assert_eq!(casc.x_feed_j, 0.0072);
        assert_eq!(casc.x_prod_j, 0.05);
        assert_eq!(casc.x_tail_j, 0.0025);
        assert_eq!(casc.mat_feed.mass, 1.0);
        assert_approx(casc.mat_feed.get(U234), 5.5e-05, "U234 fraction");
        assert_approx(casc.mat_feed.get(U235), 0.0072, "U235 fraction");
        assert_approx(casc.mat_feed.get(U238), 0.992745, "U238 fraction");
        assert_eq!(casc.mat_feed.comp.len(), 3);
    }

    #[test]
    fn reset_xjs_updates_feed_assay() {
        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, 0.0092), (U238, 0.9908)]));
        casc.reset_xjs();
        assert_eq!(casc.x_feed_j, 0.0092);
    }

    #[test]
    fn mass_ratio_formulas() {
        // Ports test_prod_per_feed, test_tail_per_feed, test_tail_per_prod.
        let (xf, xp, xt) = (0.0072_f64, 0.05_f64, 0.0025_f64);
        assert_rel_close(
            prod_per_feed(xf, xp, xt),
            (xf - xt) / (xp - xt),
            1e-4,
            "prod_per_feed",
        );
        assert_approx(
            tail_per_feed(xf, xp, xt),
            (xf - xp) / (xt - xp),
            "tail_per_feed",
        );
        assert_rel_close(
            tail_per_prod(xf, xp, xt),
            (xf - xp) / (xt - xf),
            1e-4,
            "tail_per_prod",
        );
    }

    #[test]
    fn mass_ratio_reciprocal_identities() {
        let (xf, xp, xt) = (0.0072_f64, 0.05_f64, 0.0025_f64);
        assert_approx(
            feed_per_prod(xf, xp, xt) * prod_per_feed(xf, xp, xt),
            1.0,
            "fpp*ppf",
        );
        assert_approx(
            feed_per_tail(xf, xp, xt) * tail_per_feed(xf, xp, xt),
            1.0,
            "fpt*tpf",
        );
        assert_approx(
            prod_per_tail(xf, xp, xt) * tail_per_prod(xf, xp, xt),
            1.0,
            "ppt*ttp",
        );
        // Cut relation: P/F + T/F == 1.
        assert_approx(
            prod_per_feed(xf, xp, xt) + tail_per_feed(xf, xp, xt),
            1.0,
            "cut closure",
        );
    }

    #[test]
    fn feed_product_tails_quantities() {
        // Ports test_prod, test_feed, test_tails:
        // 15.1596 kg feed <-> 1.5 kg product <-> 13.6596 kg tails.
        let (xf, xp, xt) = (0.0072_f64, 0.05_f64, 0.0025_f64);
        let (feed, prod, tails) = (15.1596_f64, 1.5_f64, 13.6596_f64);
        let rel = 1e-4;

        assert_rel_close(
            feed * prod_per_feed(xf, xp, xt),
            prod,
            rel,
            "prod from feed",
        );
        assert_rel_close(
            tails * prod_per_tail(xf, xp, xt),
            prod,
            rel,
            "prod from tails",
        );

        assert_rel_close(
            prod * feed_per_prod(xf, xp, xt),
            feed,
            rel,
            "feed from prod",
        );
        assert_rel_close(
            tails * feed_per_tail(xf, xp, xt),
            feed,
            rel,
            "feed from tails",
        );

        assert_rel_close(
            feed * tail_per_feed(xf, xp, xt),
            tails,
            rel,
            "tails from feed",
        );
        assert_rel_close(
            prod * tail_per_prod(xf, xp, xt),
            tails,
            rel,
            "tails from prod",
        );
    }

    #[test]
    fn alphastar_matches_power_law() {
        // Port of test_alphastar_i.
        let (alpha, mstar, m_i) = (1.05_f64, 236.5_f64, 235.0_f64);
        let obs = alphastar_i(alpha, mstar, m_i);
        assert_approx(obs, alpha.powf(mstar - m_i), "alphastar power law");
    }

    #[test]
    fn alphastar_monotonicity_and_unity() {
        let (alpha, mstar) = (1.05_f64, 236.5_f64);
        let lighter = alphastar_i(alpha, mstar, 235.0);
        let at_star = alphastar_i(alpha, mstar, mstar);
        let heavier = alphastar_i(alpha, mstar, 238.0);
        assert!(lighter > 1.0, "nuclides below M* must be product-side");
        assert!((at_star - 1.0).abs() < 1e-12, "alpha*(M*) == 1");
        assert!(heavier < 1.0, "nuclides above M* must be tails-side");
        assert!(lighter > heavier);
        // Strictly decreasing across ascending masses.
        let masses = [232.0, 234.0, 235.0, 236.0, 238.0];
        let factors: Vec<f64> = masses
            .iter()
            .map(|&mi| alphastar_i(alpha, mstar, mi))
            .collect();
        for pair in factors.windows(2) {
            assert!(pair[0] > pair[1]);
        }
    }

    #[test]
    fn solve_numeric_hits_targets_and_mass_relations() {
        // Solve the natural-uranium default cascade at fixed M*; the oracle
        // checks these relations implicitly through check_NU.
        let orig = default_uranium_cascade();
        let solved = solve_numeric(&orig, INTEGRATION_TOLERANCE, 100).expect("NU solve");

        let rel = 1e-5;
        assert_rel_close(solved.mat_prod.get(U235), 0.05, rel, "prod assay");
        assert_rel_close(solved.mat_tail.get(U235), 0.0025, rel, "tails assay");

        let ppf = prod_per_feed(0.0072, 0.05, 0.0025);
        let tpf = tail_per_feed(0.0072, 0.05, 0.0025);
        assert_rel_close(solved.mat_prod.mass, ppf, rel, "P/F");
        assert_rel_close(solved.mat_tail.mass, tpf, rel, "T/F");
        assert_rel_close(ppf + tpf, 1.0, rel, "cut closure");

        // Streams are normalized and conserve the key-to-key flow ratio.
        let prod_sum: f64 = solved.mat_prod.comp.values().sum();
        let tail_sum: f64 = solved.mat_tail.comp.values().sum();
        assert!((prod_sum - 1.0).abs() < 1e-12);
        assert!((tail_sum - 1.0).abs() < 1e-12);
        let r_feed = 0.0072 / 0.992745;
        let r_prod = solved.mat_prod.get(U235) / solved.mat_prod.get(U238);
        let r_tail = solved.mat_tail.get(U235) / solved.mat_tail.get(U238);
        assert!(r_prod > r_feed && r_feed > r_tail);

        // Flow-rate and SWU analytics are positive and self-consistent.
        assert!(solved.l_t_per_feed > 0.0);
        assert!(solved.swu_per_feed > 0.0);
        assert_rel_close(
            solved.swu_per_prod,
            solved.swu_per_feed / ppf,
            1e-12,
            "swu_per_prod identity",
        );

        // Input cascade untouched.
        assert_eq!(orig.l_t_per_feed, 0.0);
        assert_eq!(orig.mat_prod.comp.len(), 0);
    }

    #[test]
    fn solve_numeric_non_convergence_error() {
        // Zero iterations can never meet the targets: the port maps the
        // silent unconverged return of upstream onto Error::NoConvergence.
        let casc = default_uranium_cascade();
        match solve_numeric(&casc, 1e-11, 0) {
            Err(Error::NoConvergence { iterations }) => assert_eq!(iterations, 0),
            other => panic!("expected NoConvergence, got {other:?}"),
        }
    }

    #[test]
    fn reject_empty_feed() {
        let casc = Cascade::new();
        assert!(matches!(
            solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER),
            Err(Error::BadComposition { .. })
        ));
    }

    #[test]
    fn reject_missing_key_nuclides() {
        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, 1.0)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("missing k key must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");

        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U238, 1.0)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("missing j key must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }

    #[test]
    fn reject_zero_key_abundance() {
        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, 0.0), (U238, 1.0)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("zero j abundance must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }

    #[test]
    fn reject_degenerate_targets() {
        let mut casc = default_uranium_cascade();
        casc.x_tail_j = casc.x_prod_j;
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("coincident assays must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }

    #[test]
    fn reject_unknown_atomic_mass() {
        // Metastable ids have no AME2020 ground-state mass entry.
        let am242m = NuclideId::from_name("Am242_m1").unwrap();
        let mut casc = default_uranium_cascade();
        casc.mat_feed =
            Stream::from_comp(comp(&[(am242m, 0.01), (U235, 0.0072), (U238, 0.982745)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("unknown mass must fail");
        assert_eq!(err, Error::MissingMass(am242m));
    }

    #[test]
    fn multicomponent_sample_feed() {
        // Port of check_sample_feed (numeric solver).
        let mut orig = default_uranium_cascade();
        orig.x_prod_j = 0.06;
        orig.mat_feed = Stream::from_comp(comp(&[
            (U232, 1.1e-9),
            (U234, 0.00021),
            (U235, 0.0092),
            (U236, 0.0042),
            (U238, 0.9863899989),
        ]));
        // Inner tolerance is the documented default; 1e-11 sits below the
        // double-precision noise floor of the secant loop (see module notes).
        let casc = multicomponent(&orig, INTEGRATION_TOLERANCE, 100).expect("sample feed solve");

        let rel = 1e-5;
        assert_rel_close(casc.mat_prod.get(U235), 0.06, rel, "prod U235");
        assert_rel_close(casc.mat_tail.get(U235), 0.0025, rel, "tail U235");

        assert_rel_close(casc.mat_feed.mass, 1.0, rel, "feed mass");
        assert_rel_close(casc.mat_prod.mass, 0.11652173913043479, rel, "prod mass");
        assert_rel_close(casc.mat_tail.mass, 0.88347826086956527, rel, "tail mass");

        assert_rel_close(casc.N, 26.864660071132583, rel, "N");
        assert_rel_close(casc.M, 16.637695259838416, rel, "M");

        assert_rel_close(casc.Mstar, 236.57708506549994, rel, "Mstar");

        assert_rel_close(casc.l_t_per_feed, 357.3888391866117, rel, "L_t/F");
        assert_rel_close(casc.swu_per_feed, 0.9322804173594426, rel, "SWU/F");
        assert_rel_close(casc.swu_per_prod, 8.000914029577306, rel, "SWU/P");
    }

    #[test]
    fn multicomponent_natural_uranium() {
        // Port of check_NU (numeric solver).
        let orig = default_uranium_cascade();
        let casc = multicomponent(&orig, INTEGRATION_TOLERANCE, 100).expect("NU solve");

        let rel = 1e-5;
        assert_rel_close(casc.mat_prod.get(U235), 0.05, rel, "prod U235");
        assert_rel_close(casc.mat_tail.get(U235), 0.0025, rel, "tail U235");

        assert_rel_close(casc.mat_feed.mass, 1.0, rel, "feed mass");
        assert_rel_close(casc.mat_prod.mass, 0.0989473684211, rel, "prod mass");
        assert_rel_close(casc.mat_tail.mass, 0.901052631579, rel, "tail mass");

        let nm_rel = 1e-4;
        assert_rel_close(casc.N, 27.183583424704818, nm_rel, "N");
        assert_rel_close(casc.M, 13.387464890476533, nm_rel, "M");

        assert_rel_close(casc.Mstar, 236.5621860655, rel, "Mstar");

        assert_rel_close(casc.l_t_per_feed, 288.62731727645644, rel, "L_t/F");
        assert_rel_close(casc.swu_per_feed, 0.761263453429, rel, "SWU/F");
        assert_rel_close(casc.swu_per_prod, 7.69362000806, rel, "SWU/P");
    }

    #[test]
    fn multicomponent_vision_feed() {
        // Port of check_vision (numeric solver).
        let mut orig = default_uranium_cascade();
        orig.x_prod_j = 0.055;
        orig.mat_feed = Stream::from_comp(comp(&[
            (U234, 0.000183963025893197),
            (U235, 0.00818576605617839),
            (U236, 0.00610641667100979),
            (U238, 0.985523854246919),
        ]));
        let casc = multicomponent(&orig, INTEGRATION_TOLERANCE, 100).expect("vision solve");

        assert_rel_close(casc.mat_prod.get(U235), 0.055, 1e-5, "prod U235");
        assert_rel_close(casc.mat_tail.get(U235), 0.0025, 1e-5, "tail U235");

        assert_rel_close(casc.mat_feed.mass, 1.0, 1e-5, "feed mass");
        assert_rel_close(casc.mat_prod.mass, 0.10830030583196934, 1e-5, "prod mass");
        assert_rel_close(casc.mat_tail.mass, 0.89169969416803063, 1e-5, "tail mass");

        assert_rel_close(casc.N, 27.38162850698868, 1e-2, "N");
        assert_rel_close(casc.M, 15.09646512546496, 1e-2, "M");

        assert_rel_close(casc.Mstar, 236.58177606549995, 1e-4, "Mstar");

        assert_rel_close(casc.l_t_per_feed, 326.8956175003255, 1e-4, "L_t/F");
        assert_rel_close(casc.swu_per_feed, 0.85102089049, 1e-4, "SWU/F");
        assert_rel_close(casc.swu_per_prod, 7.85797310499, 1e-4, "SWU/P");
    }

    #[test]
    fn multicomponent_tungsten_von_halle() {
        // Port of check_tungsten (von Halle 1987, numeric solver).
        let mut orig = Cascade::new();
        orig.alpha = 1.16306;
        orig.Mstar = 181.3;
        orig.j = W180;
        orig.k = W186;
        orig.N = 30.0;
        orig.M = 10.0;
        orig.x_prod_j = 0.5109;
        orig.x_tail_j = 0.00014;
        orig.mat_feed = Stream::from_comp(comp(&[
            (W180, 0.0014),
            (W182, 0.26416),
            (W183, 0.14409),
            (W184, 0.30618),
            (W186, 0.28417),
        ]));
        orig.reset_xjs();
        let casc = multicomponent(&orig, 1e-7, 100).expect("tungsten solve");

        assert_rel_close(casc.mat_prod.get(W180), 0.5109, 1e-5, "prod W180");
        assert_rel_close(casc.mat_tail.get(W180), 0.00014, 1e-5, "tail W180");

        assert_rel_close(casc.mat_feed.mass, 1.0, 1e-5, "feed mass");
        assert_rel_close(casc.mat_prod.mass, 0.0024669120526274574, 1e-5, "prod mass");
        assert_rel_close(casc.mat_tail.mass, 0.99753308794737272, 1e-5, "tail mass");

        assert_rel_close(casc.N, 43.557515688533513, 1e-2, "N");
        assert_rel_close(casc.M, 11.49556481009056, 1e-2, "M");

        assert_rel_close(casc.Mstar, 181.16425540249995, 1e-4, "Mstar");

        assert_rel_close(casc.l_t_per_feed, 96.81774564292206, 1e-3, "L_t/F");
        assert_rel_close(casc.swu_per_feed, 2.22221945305, 1e-3, "SWU/F");
        assert_rel_close(casc.swu_per_prod, 900.810164953, 1e-3, "SWU/P");
    }

    #[test]
    fn multicomponent_replaces_out_of_range_mstar() {
        // An M* that brackets neither key gets replaced by the midpoint of
        // the key masses, landing on the same optimum as the default guess.
        let mut orig = default_uranium_cascade();
        orig.Mstar = 100.0;
        let casc_bad = multicomponent(&orig, INTEGRATION_TOLERANCE, 100).expect("bad-guess solve");
        let casc_good = multicomponent(&default_uranium_cascade(), INTEGRATION_TOLERANCE, 100)
            .expect("good-guess solve");
        assert!(
            (casc_bad.Mstar - casc_good.Mstar).abs() < 1e-3,
            "M* optimum independent of initial guess: {} vs {}",
            casc_bad.Mstar,
            casc_good.Mstar
        );
    }

    #[test]
    fn nucid_helpers_agree_with_names() {
        assert_eq!(NuclideId::from_name("U235").unwrap(), U235);
        assert_eq!(NuclideId::from_name("W186").unwrap(), W186);
        assert_eq!(nucid(922350000), U235);
    }

    #[test]
    fn reject_alpha_not_greater_than_one() {
        for alpha in [1.0, 0.0, -1.0, 0.999_999_999_999] {
            let mut casc = default_uranium_cascade();
            casc.alpha = alpha;
            let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
                .expect_err("alpha <= 1 must fail");
            assert!(matches!(err, Error::BadComposition { .. }), "{err}");
        }
    }

    #[test]
    fn reject_out_of_order_assays() {
        let cases: [(f64, f64, f64); 6] = [
            (0.0025, 0.0072, 0.0072), // prod == feed
            (0.0072, 0.0072, 0.0025), // feed == tail
            (0.05, 0.0072, 0.0025),   // tail > prod
            (0.0, 0.0072, 0.0025),    // tail == 0
            (0.0025, 0.0072, 1.0),    // prod == 1
            (-0.001, 0.0072, 0.05),   // tail < 0
        ];
        for (xt, xf, xp) in cases {
            let mut casc = default_uranium_cascade();
            casc.x_tail_j = xt;
            casc.x_feed_j = xf;
            casc.x_prod_j = xp;
            let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
                .expect_err("out-of-order assays must fail");
            assert!(matches!(err, Error::BadComposition { .. }), "{err}");
        }
    }

    #[test]
    fn reject_mstar_outside_key_masses_for_numeric_solve() {
        let mut casc = default_uranium_cascade();
        casc.Mstar = 100.0; // below both U-235 and U-238
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("out-of-bounds M* must fail in solve_numeric");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");

        let mut casc = default_uranium_cascade();
        casc.Mstar = 300.0; // above both keys
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("out-of-bounds M* must fail in solve_numeric");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }

    #[test]
    fn multicomponent_resets_out_of_range_mstar() {
        // multicomponent may reset an invalid user-supplied M* to the midpoint
        // and still converge, while solve_numeric rejects the same input.
        let mut orig = default_uranium_cascade();
        orig.Mstar = 100.0;
        let casc = multicomponent(&orig, INTEGRATION_TOLERANCE, 100)
            .expect("multicomponent resets invalid M*");
        assert!(casc.Mstar > 235.0 && casc.Mstar < 238.0);
    }

    #[test]
    fn reject_negative_feed_fraction() {
        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, -0.1), (U238, 1.1)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("negative feed fraction must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }

    #[test]
    fn reject_nan_key_abundance() {
        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, f64::NAN), (U238, 1.0)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("NaN j abundance must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");

        let mut casc = default_uranium_cascade();
        casc.mat_feed = Stream::from_comp(comp(&[(U235, 1.0), (U238, f64::NAN)]));
        let err = solve_numeric(&casc, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER)
            .expect_err("NaN k abundance must fail");
        assert!(matches!(err, Error::BadComposition { .. }), "{err}");
    }
}
