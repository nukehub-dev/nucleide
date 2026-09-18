//! Uncertainty propagation through the spectral folds over caller-supplied
//! MVN blocks — the UQ-lite hook, gated like the landed U1–U4/U7 pattern.
//!
//! [`fold_uq`] perturbs the stacked `[flux, response]` vector with seeded
//! MVN draws from `linalg::sample::sample_mvn` (Cholesky primary,
//! eigen-clip fallback, pinned `seed` → bit-identical streams) and refolds
//! the metric per draw. Deltas are **relative-unit** perturbations applied
//! as `nominal·(1 + δ)` (the SANDY convention); a draw leaving the physical
//! domain (negative flux or response) surfaces as the fold's named error,
//! never a silent clip — callers keep covariances in the small-perturbation
//! regime where this cannot trigger.
//!
//! The linear folds (dpa, appm) admit moment analytics the gate compares
//! at `k` standard errors:
//!
//! - `E[m] = m₀ + scale·(fᵀμr + rᵀμf + μfᵀμr + Σ_g cov[g][G+g])` — exact
//!   for the bilinear fold of a Gaussian (`E[δf_g·δr_g] = μf_g·μr_g +
//!   cov[g][G+g]`).
//! - `std[m] = scale·sqrt(J·C·Jᵀ)` with `J = [(f⊙r)⊙(1+μr) | (f⊙r)⊙(1+μf)]`
//!   — the first-order propagation about the block mean (for relative
//!   perturbations the per-group sensitivity is the nominal product
//!   `f_g·r_g` times the other block's mean shift); for small relative
//!   blocks the neglected quartic term sits orders of magnitude inside the
//!   gate width (the in-crate tests size it at ~3e-4 of the linear term).
//!
//! No new sampling machinery lives here: `sample_mvn` and the perturbation
//! conventions are `linalg`'s.

use nucleide_linalg::sample::{apply_perturbation, sample_mvn, PerturbConvention};

use crate::error::{Error, Result};
use crate::fold::{fold_scaled, BARNS_TO_CM2};

/// Which spectral fold the UQ sweep refolds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FoldMetric {
    /// [`crate::fold::nrt_dpa`].
    NrtDpa,
    /// [`crate::fold::arc_dpa`].
    ArcDpa,
    /// [`crate::fold::gas_appm`].
    GasAppm,
    /// [`crate::fold::he_dpa_ratio`] — not a linear fold, so [`fold_uq`]
    /// rejects it loudly; ratio UQ lives in [`he_dpa_ratio_uq`], which
    /// propagates the same seeded joint-MVN block through both folds and
    /// forms the ratio per draw.
    HeDpaRatio,
}

impl FoldMetric {
    /// Stable short name for reports and the Python facade.
    pub fn name(&self) -> &'static str {
        match self {
            FoldMetric::NrtDpa => "nrt_dpa",
            FoldMetric::ArcDpa => "arc_dpa",
            FoldMetric::GasAppm => "gas_appm",
            FoldMetric::HeDpaRatio => "he_dpa_ratio",
        }
    }
}

/// Convergence report from [`fold_uq`].
#[derive(Debug, Clone, PartialEq)]
pub struct UqSummary {
    /// Metric being propagated.
    pub metric: FoldMetric,
    /// Nominal (unperturbed) metric value.
    pub nominal: f64,
    /// Sample mean over the `n` perturbed folds.
    pub mean: f64,
    /// Unbiased sample standard deviation over the draws.
    pub std: f64,
    /// Exact expectation of the bilinear fold under the Gaussian block.
    pub expected: f64,
    /// First-order propagated standard deviation.
    pub analytic_std: f64,
    /// Gate width multiplier supplied by the caller.
    pub k: f64,
    /// Draw count.
    pub n: usize,
    /// Pinned seed that produced the draws.
    pub seed: u64,
    /// Both moment gates passed at `k` standard errors.
    pub passed: bool,
}

/// Scale factor mapping the raw fold sum to the metric's units.
fn metric_scale(metric: FoldMetric, seconds: f64) -> f64 {
    match metric {
        FoldMetric::NrtDpa | FoldMetric::ArcDpa => BARNS_TO_CM2 * seconds,
        FoldMetric::GasAppm => 1.0e-18 * seconds,
        FoldMetric::HeDpaRatio => 1.0,
    }
}

/// Propagate caller-block uncertainty through one spectral fold.
///
/// `mean_delta`/`cov` describe relative perturbations of the stacked
/// `[flux, response]` vector (dimension `2G`). Draws are seeded and
/// reproducible; the summary carries both the sample moments and the
/// analytic propagation, with the `k`-standard-error gates evaluated and
/// echoed back. `FoldMetric::HeDpaRatio` is a loud
/// [`Error::NotYetSupported`] naming [`he_dpa_ratio_uq`]: the ratio needs
/// both the He and the damage response, which the single-response `fold_uq`
/// shape cannot carry.
#[allow(clippy::too_many_arguments)] // the fold tuple plus the UQ block is the natural call shape
pub fn fold_uq(
    metric: FoldMetric,
    flux: &[f64],
    response: &[f64],
    bounds: &[f64],
    seconds: f64,
    mean_delta: &[f64],
    cov: &[Vec<f64>],
    n: usize,
    seed: u64,
    k: f64,
) -> Result<UqSummary> {
    if matches!(metric, FoldMetric::HeDpaRatio) {
        return Err(Error::NotYetSupported(
            "UQ on the He/dpa ratio needs both responses; use he_dpa_ratio_uq",
        ));
    }
    if !k.is_finite() || k <= 0.0 {
        return Err(Error::NonPositive("k"));
    }
    let g = flux.len();
    if g == 0 {
        return Err(Error::Empty);
    }
    if mean_delta.len() != 2 * g {
        return Err(Error::DimensionMismatch {
            what: "mean_delta",
            expected: 2 * g,
            got: mean_delta.len(),
        });
    }
    // Validate the nominal fold up front (bounds, signs, seconds).
    let is_appm = matches!(metric, FoldMetric::GasAppm);
    let nominal = fold_scaled(flux, response, bounds, seconds, is_appm)?;
    // Draw once up front: linalg validates the block (dimension, symmetry,
    // finiteness, PSD-adjacency) and pins the factorisation path before
    // any per-draw work happens.
    let set = sample_mvn(mean_delta, cov, n, seed)?;
    let scale = metric_scale(metric, seconds);

    // Exact expectation of the bilinear fold under δ ~ N(μ, C) with
    // relative perturbations: E[(f(1+δf))(r(1+δr))] = f·r·(1 + μf + μr +
    // μf·μr + cov_fr), since E[δf·δr] = μf·μr + cov_fr. The block shape is
    // guaranteed (2G)×(2G) by sample_mvn's validation above.
    let (mu_f, mu_r) = mean_delta.split_at(g);
    let mut expected = 0.0;
    for i in 0..g {
        expected +=
            flux[i] * response[i] * (1.0 + mu_f[i] + mu_r[i] + mu_f[i] * mu_r[i] + cov[i][g + i]);
    }
    expected *= scale;

    // First-order variance about the block mean. For relative perturbations
    // the sensitivity of the group product f(1+δf)·r(1+δr) to δf is
    // f·r·(1+μr) (and symmetrically for δr), so
    // J = scale·[(f⊙r)⊙(1+μr) | (f⊙r)⊙(1+μf)] and Var ≈ scale²·J·C·Jᵀ.
    let j_at = |i: usize| -> f64 {
        let g_i = i % g;
        let other_mean = if i < g { mu_r[g_i] } else { mu_f[g_i] };
        scale * flux[g_i] * response[g_i] * (1.0 + other_mean)
    };
    let mut jcj = 0.0;
    for (i, row) in cov.iter().enumerate() {
        let j_i = j_at(i);
        for (j, c_ij) in row.iter().enumerate() {
            jcj += j_i * c_ij * j_at(j);
        }
    }
    let analytic_std = jcj.max(0.0).sqrt();

    let mut draws: Vec<f64> = Vec::with_capacity(n);
    for delta in &set.samples {
        let (d_f, d_r) = delta.split_at(g);
        let pert_f = apply_perturbation(flux, d_f, PerturbConvention::Relative)?;
        let pert_r = apply_perturbation(response, d_r, PerturbConvention::Relative)?;
        draws.push(fold_scaled(&pert_f, &pert_r, bounds, seconds, is_appm)?);
    }

    // Unbiased sample moments (1/(n-1) variance), like linalg's estimators.
    let n_f = n as f64;
    let mean = draws.iter().sum::<f64>() / n_f;
    let std = if n < 2 {
        0.0
    } else {
        let var = draws.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n_f - 1.0);
        var.max(0.0).sqrt()
    };

    // Honest k-SE gates: the sample mean against the exact expectation, and
    // the sample standard deviation against the first-order propagation
    // (whose own normal-theory standard error is σ/√(2(n−1))).
    let mean_se = analytic_std / n_f.sqrt();
    let std_se = if n >= 2 {
        analytic_std / (2.0 * (n as f64 - 1.0)).sqrt()
    } else {
        0.0
    };
    let passed = (mean - expected).abs() <= k * mean_se && (std - analytic_std).abs() <= k * std_se;

    Ok(UqSummary {
        metric,
        nominal,
        mean,
        std,
        expected,
        analytic_std,
        k,
        n,
        seed,
        passed,
    })
}

/// Standard normal CDF via Abramowitz & Stegun 7.1.26.
///
/// `Φ(x) = 1 − φ(x)·(b₁·t + … + b₅·t⁵)` with `t = 1/(1 + 0.2316419·|x|)`,
/// `|ε(x)| ≤ 7.5e−8`, and `Φ(−x) = 1 − Φ(x)` by construction. Universal
/// mathematics, pinned by reference values in the test module — plenty
/// accurate for quantile gates whose widths sit at 1e−3 relative.
fn normal_cdf(x: f64) -> f64 {
    const A: f64 = 0.2316419;
    const B: [f64; 5] = [
        0.319381530,
        -0.356563782,
        1.781477937,
        -1.821255978,
        1.330274429,
    ];
    const INV_SQRT_2PI: f64 = 0.3989422804014327; // 1/√(2π)
    let ax = x.abs();
    let t = 1.0 / (1.0 + A * ax);
    let poly = ((((B[4] * t + B[3]) * t + B[2]) * t + B[1]) * t + B[0]) * t;
    let tail = INV_SQRT_2PI * (-0.5 * ax * ax).exp() * poly;
    if x >= 0.0 {
        1.0 - tail
    } else {
        tail
    }
}

/// Implied joint-normal moments of `(He, dpa)` for the ratio interval:
///
/// exact bilinear expectations with the first-order (co)variances — the
/// same ingredients as the delta-method summary above.
struct RatioMoments {
    mu_h: f64,
    mu_d: f64,
    var_h: f64,
    var_d: f64,
    cov: f64,
}

/// CDF of the Fieller ratio construction at `t`.
///
/// Under the joint-normal approximation, `H − t·D` is normal with mean
/// `μH − t·μD` and variance `V(t) = VH − 2·t·C + t²·VD`, so
/// `F(t) = Φ((t·μD − μH)/√V(t))`. Where `V(t) ≤ 0` (a degenerate,
/// perfectly-correlated direction — measure-zero in t) the CDF jumps and
/// the function reports the jump side so bisection still brackets. This
/// construction assumes negligible non-positive-dpa mass — the regime the
/// loud per-draw [`Error::ZeroDpa`] check enforces (any offending draw
/// fails the whole summary instead of entering it).
fn fieller_cdf(t: f64, m: &RatioMoments) -> f64 {
    let num = t * m.mu_d - m.mu_h;
    let v = m.var_h - 2.0 * t * m.cov + t * t * m.var_d;
    if v <= 0.0 {
        if num > 0.0 {
            return 1.0;
        }
        if num < 0.0 {
            return 0.0;
        }
        return 0.5;
    }
    normal_cdf(num / v.sqrt())
}

/// Density of the Fieller construction at `t` (quotient rule on the CDF).
///
/// `f(t) = φ(z)·[μD·V(t) − num·(t·VD − C)] / V(t)^{3/2}` with
/// `z = num/√V(t)`. Drives the quantile standard error
/// `√(p·(1−p)/(n·f²))`; zero at degenerate directions (the gate then goes
/// infinitely wide there — honest about the mass point). Like the CDF, this
/// is the construction's formal density under the negligible-non-positive-dpa-mass
/// assumption — bisection only ever evaluates it inside the bracketing
/// regime, where the construction is monotone.
fn fieller_pdf(t: f64, m: &RatioMoments) -> f64 {
    const INV_SQRT_2PI: f64 = 0.3989422804014327; // 1/√(2π)
    let num = t * m.mu_d - m.mu_h;
    let v = m.var_h - 2.0 * t * m.cov + t * t * m.var_d;
    if v <= 0.0 {
        return 0.0;
    }
    let z = num / v.sqrt();
    let dz = (m.mu_d * v - num * (t * m.var_d - m.cov)) / (v * v.sqrt());
    INV_SQRT_2PI * (-0.5 * z * z).exp() * dz
}

/// Invert the Fieller CDF by bisection around `center`.
///
/// Starts at `±width` and doubles geometrically until `p` brackets (up to
/// 1000 doublings), then bisects 200 times. A fully degenerate block
/// (zero spread) returns `center` directly; an unbracketable `p` — only
/// reachable with overflowed moments or genuine non-positive-dpa mass —
/// fails loudly with [`Error::NonFinite`].
fn fieller_quantile(p: f64, m: &RatioMoments, center: f64, width: f64) -> Result<f64> {
    if m.var_h == 0.0 && m.var_d == 0.0 {
        return Ok(center);
    }
    let mut w = if width > 0.0 {
        width
    } else {
        center.abs().max(1.0) * 1e-6
    };
    let (mut lo, mut hi) = (center - w, center + w);
    let (mut f_lo, mut f_hi) = (fieller_cdf(lo, m), fieller_cdf(hi, m));
    for _ in 0..1000 {
        if f_lo <= p && p <= f_hi {
            break;
        }
        w *= 2.0;
        if !w.is_finite() {
            return Err(Error::NonFinite("propagated moments"));
        }
        lo = center - w;
        hi = center + w;
        f_lo = fieller_cdf(lo, m);
        f_hi = fieller_cdf(hi, m);
    }
    if !(f_lo <= p && p <= f_hi) {
        return Err(Error::NonFinite("propagated moments"));
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if fieller_cdf(mid, m) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Linear-interpolating order-statistic quantile of sorted draws.
///
/// Position `p·(n−1)` lerped between adjacent order statistics
/// (deterministic, documented method — not a statistical claim).
fn draw_quantile(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    debug_assert!(n >= 1);
    let pos = p * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
    }
}

/// Interval-augmented convergence report from [`he_dpa_ratio_uq`].
///
/// Everything [`UqSummary`] carries for the ratio (nominal, draw moments,
/// bias-corrected expectation, first-order spread, the small-perturbation
/// `passed` verdict) plus the distribution-free 68% interval: the draw
/// 16th/50th/84th percentiles against the Fieller-construction quantiles
/// with their own `k`-SE gates (`quantiles_passed`). The interval verdict
/// holds at honestly-sized blocks where the symmetric-spread gate cannot —
/// a skewed ratio needs asymmetric margins, not a wider symmetric one.
#[derive(Debug, Clone, PartialEq)]
pub struct RatioUqSummary {
    /// Always [`FoldMetric::HeDpaRatio`].
    pub metric: FoldMetric,
    /// Nominal (unperturbed) He/dpa ratio.
    pub nominal: f64,
    /// Sample mean over the `n` per-draw ratios.
    pub mean: f64,
    /// Unbiased sample standard deviation over the draws.
    pub std: f64,
    /// Second-order bias-corrected expectation of the mean.
    pub expected: f64,
    /// First-order delta-propagated standard deviation.
    pub analytic_std: f64,
    /// Draw 16th percentile (linear-interpolating order statistic).
    pub q16: f64,
    /// Draw median.
    pub q50: f64,
    /// Draw 84th percentile.
    pub q84: f64,
    /// Fieller-construction 16th percentile.
    pub expected_q16: f64,
    /// Fieller-construction median.
    pub expected_q50: f64,
    /// Fieller-construction 84th percentile.
    pub expected_q84: f64,
    /// All three quantile gates passed at `k` standard errors.
    pub quantiles_passed: bool,
    /// Gate width multiplier supplied by the caller.
    pub k: f64,
    /// Draw count.
    pub n: usize,
    /// Pinned seed that produced the draws.
    pub seed: u64,
    /// Both moment gates passed at `k` standard errors (small-perturbation
    /// regime verdict; see `quantiles_passed` for the all-regime verdict).
    pub passed: bool,
}

/// Propagate caller-block uncertainty through the He/dpa ratio.
///
/// `mean_delta`/`cov` describe relative perturbations of the stacked
/// `[flux, he_response, damage_response]` vector (dimension `3G`). Draws are
/// seeded and reproducible via `linalg::sample::sample_mvn` (the same
/// machinery [`fold_uq`] uses — no new sampling code); each draw refolds He
/// (`1e-18` scale) and dpa (`1e-24` scale) and forms the ratio per draw, and
/// the summary reports the draw mean ± draw standard deviation with the
/// landed `k`-standard-error gate style.
///
/// Analytic cross-check (delta method through the two bilinear folds, with
/// `E[H]` / `E[D]` the exact bilinear expectations and `Var(H)` / `Var(D)` /
/// `Cov(H, D)` the first-order Jacobian propagations — the [`fold_uq`]
/// constructions applied per fold, the cross term as `J_H·C·J_Dᵀ`):
///
/// ```text
/// r0           = E[H] / E[D]
/// expected     = r0 · (1 + Var(D)/E[D]² − Cov(H,D)/(E[H]·E[D]))
/// analytic_std = r0 · √(Var(H)/E[H]² + Var(D)/E[D]² − 2·Cov(H,D)/(E[H]·E[D]))
/// ```
///
/// The mean of ratios carries the second-order `Var(D)` bias, so the mean
/// gate compares against the bias-corrected `expected`, not the bare ratio
/// of means — gating against `r0` would fail honestly-sized blocks.
///
/// The nominal ratio comes from [`crate::fold::he_dpa_ratio`] (which
/// validates the fold shape and fails loudly with [`Error::ZeroDpa`] at zero
/// nominal dpa). A draw landing at `dpa <= 0` fails loudly with
/// [`Error::ZeroDpa`] too — never `inf`/`NaN` with a spread, never silent
/// dropping of draws.
///
/// Besides the moment summary, every call reports the distribution-free 68%
/// interval (draw 16th/50th/84th percentiles) gated against the
/// Fieller-construction quantiles — the all-regime verdict
/// (`quantiles_passed`) that holds at honestly-sized blocks where the
/// symmetric-spread gate cannot.
#[allow(clippy::too_many_arguments)] // the two-response fold tuple plus the UQ block
pub fn he_dpa_ratio_uq(
    flux: &[f64],
    he_response: &[f64],
    damage_response: &[f64],
    bounds: &[f64],
    seconds: f64,
    mean_delta: &[f64],
    cov: &[Vec<f64>],
    n: usize,
    seed: u64,
    k: f64,
) -> Result<RatioUqSummary> {
    use crate::fold::he_dpa_ratio;

    if !k.is_finite() || k <= 0.0 {
        return Err(Error::NonPositive("k"));
    }
    let g = flux.len();
    if g == 0 {
        return Err(Error::Empty);
    }
    if mean_delta.len() != 3 * g {
        return Err(Error::DimensionMismatch {
            what: "mean_delta",
            expected: 3 * g,
            got: mean_delta.len(),
        });
    }
    // Validate the nominal ratio up front (response shapes, bounds, signs,
    // seconds, and the zero-dpa stance).
    let nominal = he_dpa_ratio(flux, he_response, damage_response, bounds, seconds)?;
    // Draw once up front: linalg validates the block before any per-draw
    // work happens.
    let set = sample_mvn(mean_delta, cov, n, seed)?;
    let scale_he = 1.0e-18 * seconds;
    let scale_dpa = BARNS_TO_CM2 * seconds;

    // Block thirds: [flux | he_response | damage_response].
    let (mu_f, rest) = mean_delta.split_at(g);
    let (mu_h, mu_d) = rest.split_at(g);
    // Exact bilinear expectation per fold: E = scale·Σ f·r·(1 + μf + μr +
    // μf·μr + cov_fr) (the `fold_uq` formula).
    let mut exp_he = 0.0;
    let mut exp_dpa = 0.0;
    for i in 0..g {
        exp_he += flux[i]
            * he_response[i]
            * (1.0 + mu_f[i] + mu_h[i] + mu_f[i] * mu_h[i] + cov[i][g + i]);
        exp_dpa += flux[i]
            * damage_response[i]
            * (1.0 + mu_f[i] + mu_d[i] + mu_f[i] * mu_d[i] + cov[i][2 * g + i]);
    }
    exp_he *= scale_he;
    exp_dpa *= scale_dpa;
    if exp_dpa == 0.0 {
        return Err(Error::ZeroDpa);
    }

    // First-order Jacobians per fold (the `fold_uq` construction): the
    // sensitivity of `scale·f(1+δf)·r(1+δr)` to δf is
    // `scale·f·r·(1+μr)`, symmetrically for δr.
    let j_he = |i: usize| -> f64 {
        let group = i % g;
        if i < g {
            scale_he * flux[group] * he_response[group] * (1.0 + mu_h[group])
        } else if i < 2 * g {
            scale_he * flux[group] * he_response[group] * (1.0 + mu_f[group])
        } else {
            0.0
        }
    };
    let j_dpa = |i: usize| -> f64 {
        let group = i % g;
        if i < g {
            scale_dpa * flux[group] * damage_response[group] * (1.0 + mu_d[group])
        } else if i < 2 * g {
            0.0
        } else {
            scale_dpa * flux[group] * damage_response[group] * (1.0 + mu_f[group])
        }
    };
    let mut var_he = 0.0;
    let mut var_dpa = 0.0;
    let mut cov_hd = 0.0;
    for (i, row) in cov.iter().enumerate() {
        for (j, c_ij) in row.iter().enumerate() {
            var_he += j_he(i) * c_ij * j_he(j);
            var_dpa += j_dpa(i) * c_ij * j_dpa(j);
            cov_hd += j_he(i) * c_ij * j_dpa(j);
        }
    }
    let r0 = exp_he / exp_dpa;
    let expected = r0 * (1.0 + var_dpa / (exp_dpa * exp_dpa) - cov_hd / (exp_he * exp_dpa));
    let analytic_var = r0
        * r0
        * (var_he / (exp_he * exp_he) + var_dpa / (exp_dpa * exp_dpa)
            - 2.0 * cov_hd / (exp_he * exp_dpa));
    let analytic_std = analytic_var.max(0.0).sqrt();
    // Absurd-magnitude blocks can overflow the analytic moments themselves;
    // a non-finite gate width would pass any draw (or none) meaninglessly,
    // so it fails loudly instead of reporting `passed` on garbage.
    if !expected.is_finite() || !analytic_std.is_finite() {
        return Err(Error::NonFinite("propagated moments"));
    }

    let mut draws: Vec<f64> = Vec::with_capacity(n);
    for delta in &set.samples {
        let (d_f, rest) = delta.split_at(g);
        let (d_h, d_d) = rest.split_at(g);
        let pert_f = apply_perturbation(flux, d_f, PerturbConvention::Relative)?;
        let pert_h = apply_perturbation(he_response, d_h, PerturbConvention::Relative)?;
        let pert_d = apply_perturbation(damage_response, d_d, PerturbConvention::Relative)?;
        let he = fold_scaled(&pert_f, &pert_h, bounds, seconds, true)?;
        let dpa = fold_scaled(&pert_f, &pert_d, bounds, seconds, false)?;
        if dpa <= 0.0 {
            return Err(Error::ZeroDpa);
        }
        draws.push(he / dpa);
    }

    // Unbiased sample moments (1/(n-1) variance), like linalg's estimators.
    let n_f = n as f64;
    let mean = draws.iter().sum::<f64>() / n_f;
    let std = if n < 2 {
        0.0
    } else {
        let var = draws.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n_f - 1.0);
        var.max(0.0).sqrt()
    };

    // Honest k-SE gates: the sample mean against the bias-corrected
    // expectation, and the sample standard deviation against the first-order
    // propagation (whose own normal-theory standard error is σ/√(2(n−1))).
    let mean_se = analytic_std / n_f.sqrt();
    let std_se = if n >= 2 {
        analytic_std / (2.0 * (n as f64 - 1.0)).sqrt()
    } else {
        0.0
    };
    let passed = (mean - expected).abs() <= k * mean_se && (std - analytic_std).abs() <= k * std_se;

    // Distribution-free 68% interval: linear-interpolating order statistics
    // of the draws, gated against the Fieller-construction quantiles with
    // the quantile standard error √(p(1−p)/(n·f²)) at the analytic density.
    let mut sorted = draws.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let (q16, q50, q84) = (
        draw_quantile(&sorted, 0.16),
        draw_quantile(&sorted, 0.5),
        draw_quantile(&sorted, 0.84),
    );
    let moments = RatioMoments {
        mu_h: exp_he,
        mu_d: exp_dpa,
        var_h: var_he,
        var_d: var_dpa,
        cov: cov_hd,
    };
    let center = r0;
    let width = 10.0 * analytic_std;
    let expected_q16 = fieller_quantile(0.16, &moments, center, width)?;
    let expected_q50 = fieller_quantile(0.5, &moments, center, width)?;
    let expected_q84 = fieller_quantile(0.84, &moments, center, width)?;
    if !expected_q16.is_finite() || !expected_q50.is_finite() || !expected_q84.is_finite() {
        return Err(Error::NonFinite("propagated moments"));
    }
    let q_se = |p: f64, q: f64| {
        let f = fieller_pdf(q, &moments);
        (p * (1.0 - p) / (n_f * f * f)).sqrt()
    };
    let quantiles_passed = [
        (0.16, q16, expected_q16),
        (0.5, q50, expected_q50),
        (0.84, q84, expected_q84),
    ]
    .into_iter()
    .all(|(p, got, want)| (got - want).abs() <= k * q_se(p, want));

    Ok(RatioUqSummary {
        metric: FoldMetric::HeDpaRatio,
        nominal,
        mean,
        std,
        expected,
        analytic_std,
        q16,
        q50,
        q84,
        expected_q16,
        expected_q50,
        expected_q84,
        quantiles_passed,
        k,
        n,
        seed,
        passed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEED: u64 = 20260915;

    /// Synthetic two-group problem tuple: (flux, response, bounds, mean, cov).
    type Problem = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<Vec<f64>>);

    /// Synthetic two-group fold and a small relative block (1–5%
    /// standard deviations — the small-perturbation regime the gate
    /// documents). No evaluated data anywhere.
    fn problem() -> Problem {
        let flux = vec![1.0e13, 3.0e12];
        let resp = vec![80.0, 160.0];
        let bounds = vec![0.0, 0.5, 20.0];
        let mean = vec![0.0, 0.0, 0.0, 0.0];
        let cov = vec![
            vec![0.0004, 0.0, 0.0, 0.0],
            vec![0.0, 0.0009, 0.0, 0.0],
            vec![0.0, 0.0, 0.0025, 0.0],
            vec![0.0, 0.0, 0.0, 0.0001],
        ];
        (flux, resp, bounds, mean, cov)
    }

    #[test]
    fn mvn_uq_recovers_propagated_moments_within_k_se() {
        let (flux, resp, bounds, mean, cov) = problem();
        let n = 20_000usize;
        let k = 5.0;
        for metric in [FoldMetric::NrtDpa, FoldMetric::ArcDpa, FoldMetric::GasAppm] {
            let s = fold_uq(metric, &flux, &resp, &bounds, 2.0, &mean, &cov, n, SEED, k).unwrap();
            assert!(s.passed, "{metric:?}: {s:?}");
            assert_eq!(s.n, n);
            assert_eq!(s.seed, SEED);
            // With zero-mean block-diagonal perturbations the exact
            // expectation is the nominal fold itself.
            let nominal = fold_scaled(
                &flux,
                &resp,
                &bounds,
                2.0,
                matches!(metric, FoldMetric::GasAppm),
            )
            .unwrap();
            assert_eq!(s.nominal, nominal);
            assert_eq!(s.expected, nominal);
            // The reported moments sit inside the gate widths the summary
            // itself echoes (mean_se = analytic_std/sqrt(n), etc.).
            let mean_se = s.analytic_std / (n as f64).sqrt();
            assert!((s.mean - s.expected).abs() <= k * mean_se);
            let std_se = s.analytic_std / (2.0 * (n as f64 - 1.0)).sqrt();
            assert!((s.std - s.analytic_std).abs() <= k * std_se);
        }
        // Determinism under the pinned seed.
        let a = fold_uq(
            FoldMetric::NrtDpa,
            &flux,
            &resp,
            &bounds,
            2.0,
            &mean,
            &cov,
            64,
            SEED,
            5.0,
        )
        .unwrap();
        let b = fold_uq(
            FoldMetric::NrtDpa,
            &flux,
            &resp,
            &bounds,
            2.0,
            &mean,
            &cov,
            64,
            SEED,
            5.0,
        )
        .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mean_bias_is_caught_exactly() {
        // A relative mean bias on the response block shifts E[m] by
        // scale·Σ f·μr — the analytic expectation must track it. The bias
        // is sized well above the f64 summation noise of the big terms.
        let (flux, resp, bounds, mut mean, cov) = problem();
        mean[2] = 0.5; // +50% response bias, group 1
        let s = fold_uq(
            FoldMetric::NrtDpa,
            &flux,
            &resp,
            &bounds,
            1.0,
            &mean,
            &cov,
            4_000,
            SEED,
            5.0,
        )
        .unwrap();
        let bias = 1.0e-24 * flux[0] * resp[0] * 0.5;
        assert!((s.expected - (s.nominal + bias)).abs() <= 1e-12 * s.nominal.abs());
        assert!(s.passed);
    }

    #[test]
    fn ratio_uq_names_its_two_response_entry_point() {
        let (flux, resp, bounds, mean, cov) = problem();
        assert!(matches!(
            fold_uq(
                FoldMetric::HeDpaRatio,
                &flux,
                &resp,
                &bounds,
                1.0,
                &mean,
                &cov,
                16,
                SEED,
                5.0
            ),
            Err(Error::NotYetSupported(_))
        ));
        assert_eq!(FoldMetric::HeDpaRatio.name(), "he_dpa_ratio");
    }

    /// Single-group ratio problem tuple: (flux, he, dpa, bounds, mean, cov).
    type RatioProblem = (
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<Vec<f64>>,
    );

    /// Single-group ratio problem: flux 1e12, He response 5e6 → He = 5.0,
    /// damage response 2e12 → dpa = 2.0, ratio of means 2.5.
    fn ratio_problem() -> RatioProblem {
        let flux = vec![1.0e12];
        let he = vec![5.0e6];
        let dpa = vec![2.0e12];
        let bounds = vec![0.0, 20.0];
        // 3G block [flux | he | dpa]: flux exact, ±10% uncorrelated on both
        // responses — independent He = 5.0 ± 0.5, dpa = 2.0 ± 0.2 (1σ).
        let mean = vec![0.0; 3];
        let cov = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.01, 0.0],
            vec![0.0, 0.0, 0.01],
        ];
        (flux, he, dpa, bounds, mean, cov)
    }

    #[test]
    fn normal_cdf_pins_reference_values() {
        // Universal mathematics (Abramowitz & Stegun Table 7.1), not outputs
        // of this implementation: Φ(0) = 1/2, Φ(1), Φ(2), antisymmetry.
        assert!((normal_cdf(0.0) - 0.5).abs() <= 1e-9);
        assert!((normal_cdf(1.0) - 0.8413447460685429).abs() <= 1e-7);
        assert!((normal_cdf(2.0) - 0.9772498680518208).abs() <= 1e-7);
        assert!((normal_cdf(-1.0) - 0.1586552539314571).abs() <= 1e-7);
        assert!((normal_cdf(5.0) - 0.9999997133484281).abs() <= 1e-7);
        assert!((normal_cdf(-5.0) - 0.0000002866515719).abs() <= 1e-7);
        // Far tails keep the absolute (not relative) error bound: the true
        // Φ(±10) = 1 − 7.61985302416047e−24 is far below f64 resolution at
        // 1.0, and the approximant stays within 7.5e−8 of it.
        assert!((normal_cdf(10.0) - (1.0 - 7.61985302416047e-24)).abs() <= 1e-7);
        assert!((normal_cdf(-10.0) - 7.61985302416047e-24).abs() <= 1e-7);
    }

    #[test]
    fn ratio_interval_holds_where_the_spread_gate_cannot() {
        // The ±10% hand-vector block: the symmetric-spread gate honestly
        // fails here (heavy ratio tails), while the Fieller interval gates
        // pass — a skewed ratio needs asymmetric margins.
        let (flux, he, dpa_resp, bounds, mean, cov) = ratio_problem();
        let s = he_dpa_ratio_uq(
            &flux, &he, &dpa_resp, &bounds, 1.0, &mean, &cov, 20_000, SEED, 5.0,
        )
        .unwrap();
        assert!(!s.passed, "spread gate should honestly fail at ±10%");
        assert!(s.quantiles_passed, "interval gates should pass: {s:?}");
        // Ordered interval around the median, asymmetric upward (positive skew).
        assert!(s.q16 < s.q50 && s.q50 < s.q84);
        assert!(s.expected_q16 < s.expected_q50 && s.expected_q50 < s.expected_q84);
        assert!(
            s.q84 - s.q50 > s.q50 - s.q16,
            "draw interval should skew upward: {:?}",
            (s.q16, s.q50, s.q84)
        );
        assert!(
            s.expected_q84 - s.expected_q50 > s.expected_q50 - s.expected_q16,
            "analytic interval should skew upward"
        );
        // The median is nearly unbiased: within a percent of the nominal 2.5.
        assert!((s.q50 - 2.5).abs() <= 0.025);
        assert!((s.expected_q50 - 2.5).abs() <= 0.025);
        // The 68% interval brackets the nominal on both sides.
        assert!(s.q16 < 2.5 && 2.5 < s.q84);
    }

    #[test]
    fn degenerate_block_collapses_the_interval_exactly() {
        // An all-zero block cannot pass `sample_mvn` (no positive
        // eigenvalue), so the degenerate quantile path is pinned directly:
        // zero spread returns the center for every probability.
        let moments = RatioMoments {
            mu_h: 5.0,
            mu_d: 2.0,
            var_h: 0.0,
            var_d: 0.0,
            cov: 0.0,
        };
        for p in [0.16, 0.5, 0.84] {
            assert_eq!(fieller_quantile(p, &moments, 2.5, 0.0).unwrap(), 2.5);
        }
        // Near-degenerate end to end (1e-12 block): the interval collapses
        // onto the nominal ratio and both verdicts pass.
        let (flux, he, dpa_resp, bounds, _, _) = ratio_problem();
        let mean = vec![0.0; 3];
        let tiny = 1e-12;
        let cov = vec![
            vec![tiny, 0.0, 0.0],
            vec![0.0, tiny, 0.0],
            vec![0.0, 0.0, tiny],
        ];
        let s = he_dpa_ratio_uq(
            &flux, &he, &dpa_resp, &bounds, 1.0, &mean, &cov, 64, SEED, 5.0,
        )
        .unwrap();
        assert!((s.nominal - 2.5).abs() <= 1e-12);
        // The 1e-6 relative block spreads draws ~3.5e-6 absolute; the
        // interval still collapses onto the nominal far inside any honest
        // block width.
        for q in [
            s.q16,
            s.q50,
            s.q84,
            s.expected_q16,
            s.expected_q50,
            s.expected_q84,
        ] {
            assert!((q - s.nominal).abs() <= 1e-4, "quantile {q} vs nominal");
        }
        assert!(s.quantiles_passed);
        assert!(s.passed);
    }

    #[test]
    fn ratio_uq_pins_the_hand_vector() {
        let (flux, he, dpa_resp, bounds, mean, cov) = ratio_problem();
        let s = he_dpa_ratio_uq(
            &flux, &he, &dpa_resp, &bounds, 1.0, &mean, &cov, 64, SEED, 5.0,
        )
        .unwrap();
        assert_eq!(s.metric, FoldMetric::HeDpaRatio);
        // Nominal is the point ratio 5/2 up to fold fp rounding.
        assert!(
            (s.nominal - 2.5).abs() <= 1e-12,
            "nominal {} vs 2.5",
            s.nominal
        );
        // First-order propagated sd: 2.5·√(0.01 + 0.01) = √0.125.
        let pinned = 0.3535533905932738f64;
        assert!(
            (s.analytic_std - pinned).abs() <= 1e-12,
            "analytic_std {} vs pinned {pinned}",
            s.analytic_std
        );
        // Mean of ratios carries the second-order Var(D) bias: the
        // expectation corrects the bare ratio 2.5·(1 + 0.01) = 2.525.
        assert!(
            (s.expected - 2.525).abs() <= 1e-12,
            "expected {} vs 2.525",
            s.expected
        );
        assert_eq!((s.n, s.seed), (64, SEED));
        // Determinism under the pinned seed.
        let twice = he_dpa_ratio_uq(
            &flux, &he, &dpa_resp, &bounds, 1.0, &mean, &cov, 64, SEED, 5.0,
        )
        .unwrap();
        assert_eq!(s, twice);
    }

    #[test]
    fn ratio_uq_converges_within_k_se_in_the_small_perturbation_regime() {
        // ±2% uncorrelated on both responses (variances 4e-4): the
        // higher-order ratio moments the first-order gate neglects sit
        // orders of magnitude inside the gate width (the documented
        // small-perturbation regime). At ±10% the heavy 1/(1+δd) tails
        // inflate the sample sd past any honest first-order k-SE width —
        // the pin test above covers that block analytically instead.
        let (flux, he, dpa_resp, bounds, _, _) = ratio_problem();
        let mean = vec![0.0; 3];
        let cov = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0004, 0.0],
            vec![0.0, 0.0, 0.0004],
        ];
        let n = 20_000usize;
        let k = 5.0;
        let s =
            he_dpa_ratio_uq(&flux, &he, &dpa_resp, &bounds, 1.0, &mean, &cov, n, SEED, k).unwrap();
        assert!(s.passed, "ratio UQ did not converge: {s:?}");
        // Analytic spot at ±2%: 2.5·√(0.0008) with the 4e-4 bias correction.
        let analytic = 2.5 * 0.0008f64.sqrt();
        assert!((s.analytic_std - analytic).abs() <= 1e-12);
        assert!((s.expected - 2.5 * 1.0004).abs() <= 1e-12);
    }

    #[test]
    fn ratio_uq_zero_dpa_is_loud() {
        let (flux, he, dpa_resp, bounds, mean, cov) = ratio_problem();
        // Zero nominal dpa: loud before any draw.
        let zero_dpa = vec![0.0];
        assert!(matches!(
            he_dpa_ratio_uq(&flux, &he, &zero_dpa, &bounds, 1.0, &mean, &cov, 8, SEED, 5.0),
            Err(Error::ZeroDpa)
        ));
        // A block that drives a draw to non-positive dpa surfaces loudly
        // (ZeroDpa at exactly zero; the fold's loud Negative below it).
        let mean_neg = vec![0.0, 0.0, -2.0];
        let out = he_dpa_ratio_uq(
            &flux, &he, &dpa_resp, &bounds, 1.0, &mean_neg, &cov, 8, SEED, 5.0,
        );
        assert!(
            matches!(out, Err(Error::ZeroDpa) | Err(Error::Negative(_))),
            "negative-dpa block should fail loud, got {out:?}"
        );
    }

    #[test]
    fn ratio_uq_overflowed_moments_are_loud() {
        // Absurd-magnitude inputs overflow the analytic moments themselves;
        // a non-finite gate width must fail loudly, never report `passed`
        // on garbage.
        let flux = vec![f64::MAX];
        let resp = vec![f64::MAX];
        let bounds = vec![0.0, 20.0];
        let mean = vec![0.0; 3];
        let cov = vec![
            vec![0.01, 0.0, 0.0],
            vec![0.0, 0.01, 0.0],
            vec![0.0, 0.0, 0.01],
        ];
        assert!(matches!(
            he_dpa_ratio_uq(&flux, &resp, &resp, &bounds, 1.0, &mean, &cov, 8, SEED, 5.0),
            Err(Error::NonFinite("propagated moments"))
        ));
    }

    #[test]
    fn bad_blocks_name_their_cause() {
        let (flux, resp, bounds, mean, cov) = problem();
        assert!(matches!(
            fold_uq(
                FoldMetric::NrtDpa,
                &flux,
                &resp,
                &bounds,
                1.0,
                &[0.0; 3],
                &cov,
                8,
                SEED,
                5.0
            ),
            Err(Error::DimensionMismatch {
                what: "mean_delta",
                ..
            })
        ));
        assert!(matches!(
            fold_uq(
                FoldMetric::NrtDpa,
                &flux,
                &resp,
                &bounds,
                1.0,
                &mean,
                &[vec![1.0]],
                8,
                SEED,
                5.0
            ),
            Err(Error::Sampling(
                nucleide_linalg::SampleError::DimensionMismatch { .. },
            ))
        ));
        assert!(matches!(
            fold_uq(
                FoldMetric::NrtDpa,
                &flux,
                &resp,
                &bounds,
                1.0,
                &mean,
                &cov,
                8,
                SEED,
                0.0
            ),
            Err(Error::NonPositive("k"))
        ));
    }
}
