//! Seeded multivariate-normal sampling over caller-supplied covariance blocks.
//!
//! UQ-lite kernel (decay-only sub-scope): draws `n`
//! reproducible samples `x ~ N(mean, cov)` for uncertainty propagation in
//! the SANDY role (seeded MVN draws over caller-supplied blocks), without
//! vendored covariance stores, transport coupling, or any ERRORR/NJOY
//! machinery. All inputs are caller-supplied; this module never reads
//! evaluated data.
//!
//! Engine difference (verified against the SANDY 1.1.0 wheel): SANDY's
//! `CategoryCov.draw_sample` factorises with SVD ("better than QR or
//! cholesky") behind `numpy.random.default_rng`; this module factorises
//! with Cholesky plus an eigen-clipping fallback behind a ChaCha8
//! `StdRng`. Same target distribution, different streams — draws are not
//! interchangeable across implementations. What agrees exactly (Tier-2
//! oracle gate at 1e-9) are the moment estimators: [`sample_mean`](crate::sample::sample_mean) is the
//! row mean like `Samples.get_mean`, and [`sample_cov`](crate::sample::sample_cov) is the unbiased
//! `1/(n-1)` covariance like `Samples.get_cov` (`pandas.DataFrame.cov`
//! default `ddof=1`).
//!
//! Method: Cholesky factorisation of `cov` (via faer 0.20, the workspace
//! backend) is the primary path. When `cov` is positive-semidefinite but
//! numerically singular (Cholesky fails), the fallback is an eigen-clipping
//! reconstruction `B = U sqrt(clip(S))` from the real-symmetric
//! eigendecomposition, with eigenvalues floored at
//! [`EIGEN_FLOOR_REL`](crate::sample::EIGEN_FLOOR_REL) times the largest eigenvalue. Which path was taken is
//! reported in [`SampleSet::method`](crate::sample::SampleSet) — never silent.
//!
//! Perturbation conventions ([`PerturbConvention`](crate::sample::PerturbConvention)):
//!
//! - `Relative`: `perturbed = nominal * (1 + delta)` — deltas are fractional
//!   deviations (relative-unit perturbations as in SANDY's perturbation
//!   model, whose `truncate_normal` assumes samples centered at 1).
//! - `Absolute`: `perturbed = nominal + delta` — deltas carry the
//!   nominal's units.
//! - `LogNormal`: `perturbed = nominal * exp(delta)` — log-space deltas
//!   (multiplicative log-normal perturbations, always non-negative for
//!   non-negative nominals).
//!
//! Log-normal sampling ([`sample_lognormal`](crate::sample::sample_lognormal)):
//! draws `x ~ N(mean_log, cov_log)` with the shared [`sample_mvn`] factor
//! path and RNG, then maps `y = exp(x)` elementwise. `mean_log`/`cov_log`
//! are log-space MVN parameters (never the moments of `y` itself); the
//! closed-form moments are `E[y_i] = exp(mu_i + C_ii/2)` and
//! `Cov(y_i, y_j) = exp(mu_i + mu_j + (C_ii + C_jj)/2) (exp(C_ij) - 1)`.
//!
//! Convergence diagnostics ([`check_convergence`](crate::sample::check_convergence)) compare the sample mean
//! and unbiased sample covariance against the inputs that generated them —
//! the same moment estimators SANDY's `Samples.get_mean` / `get_cov`
//! compute (Tier-2 gate), but with our own gate formulation: caller-supplied
//! tolerances, with the validation oracle (`validation/uq_lite_vs_sandy.py`)
//! pinning statistical ones derived from the MVN sampling variances.
//!
//! Explicitly OUT: transport-coupled UQ, ERRORR/NJOY
//! reimplementation, vendored covariance stores, MF32-resonance machinery,
//! MF40, fission-yield perturbation (named-open in [`crate::decay`]), and
//! Latin-hypercube sampling (needs a gate redesign, stays out).

use faer::{Mat, Side};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

/// Relative floor for eigen-clipping: clipped eigenvalues are at least this
/// fraction of the largest eigenvalue.
pub const EIGEN_FLOOR_REL: f64 = 1e-12;

/// Relative tolerance for covariance symmetry: inputs more asymmetric than
/// this (relative to the largest absolute entry) are rejected instead of
/// being silently symmetrised.
pub const SYMMETRY_TOL_REL: f64 = 1e-8;

/// Errors surfaced by the sampling kernel. Every rejection names its cause;
/// the factorisation path itself is reported, not errored (see
/// [`FactorMethod`]).
#[derive(Debug, Clone, PartialEq)]
pub enum SampleError {
    /// `mean` (or the sample set) is empty.
    Empty,
    /// A length does not match the problem dimension.
    DimensionMismatch {
        /// Dimension required by the leading input.
        expected: usize,
        /// Dimension actually supplied.
        got: usize,
    },
    /// A covariance row is not as long as the problem dimension.
    NonSquare {
        /// Problem dimension (row count).
        dim: usize,
        /// Offending row index.
        row: usize,
        /// Offending row length.
        got: usize,
    },
    /// A non-finite value in the named input (`"mean"` or `"covariance"`).
    NonFinite(&'static str),
    /// `cov` is asymmetric beyond [`SYMMETRY_TOL_REL`]; the value is the
    /// largest `|C[i][j] - C[j][i]|`.
    NonSymmetric {
        /// Largest absolute asymmetry found.
        deviation: f64,
    },
    /// `n == 0`: no samples requested.
    NoSamples,
    /// The covariance has no positive eigenvalue, so even the eigen-clipping
    /// fallback cannot build a factor.
    NoPositiveEigenvalue,
}

impl std::fmt::Display for SampleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SampleError::Empty => write!(f, "sampling input is empty"),
            SampleError::DimensionMismatch { expected, got } => write!(
                f,
                "sampling dimension mismatch: expected {expected}, got {got}"
            ),
            SampleError::NonSquare { dim, row, got } => {
                write!(f, "covariance row {row} has length {got}, expected {dim}")
            }
            SampleError::NonFinite(what) => {
                write!(f, "sampling input `{what}` holds a non-finite value")
            }
            SampleError::NonSymmetric { deviation } => write!(
                f,
                "covariance is asymmetric (max deviation {deviation:.3e}); \
                 symmetrise caller-side or stay within {SYMMETRY_TOL_REL:.0e} relative"
            ),
            SampleError::NoSamples => write!(f, "sampling requested zero draws"),
            SampleError::NoPositiveEigenvalue => write!(
                f,
                "covariance has no positive eigenvalue; no sampling factor exists"
            ),
        }
    }
}

impl std::error::Error for SampleError {}

/// Which factorisation produced a [`SampleSet`].
#[derive(Debug, Clone, PartialEq)]
pub enum FactorMethod {
    /// Cholesky `cov = L Lᵀ` succeeded (positive-definite input).
    Cholesky,
    /// Cholesky failed; samples use the eigen-clipping reconstruction
    /// `cov ≈ B Bᵀ` with `B = U sqrt(clip(S))`. `min_eigen`/`max_eigen`
    /// are the unclipped extremes.
    EigenClip {
        /// Smallest eigenvalue before clipping.
        min_eigen: f64,
        /// Largest eigenvalue before clipping.
        max_eigen: f64,
    },
}

impl FactorMethod {
    /// Stable short name for reports and the Python facade (`"cholesky"` /
    /// `"eigen_clip"`).
    pub fn name(&self) -> &'static str {
        match self {
            FactorMethod::Cholesky => "cholesky",
            FactorMethod::EigenClip { .. } => "eigen_clip",
        }
    }
}

/// How a sampled delta maps onto a nominal value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerturbConvention {
    /// `perturbed = nominal * (1 + delta)`; deltas are fractional.
    Relative,
    /// `perturbed = nominal + delta`; deltas carry nominal units.
    Absolute,
    /// `perturbed = nominal * exp(delta)`; deltas live in log space, so
    /// non-negative nominals stay non-negative without clamping.
    LogNormal,
}

impl PerturbConvention {
    /// Parse `"relative"` / `"absolute"` / `"lognormal"` (case-insensitive,
    /// `-`/`_`/nothing interchangeable, so `"log-normal"`, `"log_normal"`,
    /// and `"lognormal"` all match); anything else names the accepted
    /// spellings.
    pub fn parse(name: &str) -> Result<Self, String> {
        let norm: String = name
            .chars()
            .map(|c| {
                if c == '-' || c == '_' || c == ' ' {
                    '\0'
                } else {
                    c.to_ascii_lowercase()
                }
            })
            .filter(|c| *c != '\0')
            .collect();
        match norm.as_str() {
            "relative" | "rel" => Ok(PerturbConvention::Relative),
            "absolute" | "abs" => Ok(PerturbConvention::Absolute),
            "lognormal" | "lognorm" | "ln" => Ok(PerturbConvention::LogNormal),
            other => Err(format!(
                "unknown perturbation convention `{other}` \
                 (supported: \"relative\", \"absolute\", \"lognormal\")"
            )),
        }
    }

    /// Apply one delta to one nominal value under this convention.
    pub fn apply(self, nominal: f64, delta: f64) -> f64 {
        match self {
            PerturbConvention::Relative => nominal * (1.0 + delta),
            PerturbConvention::Absolute => nominal + delta,
            PerturbConvention::LogNormal => nominal * delta.exp(),
        }
    }
}

/// One seeded draw batch plus the factorisation path that produced it.
#[derive(Debug, Clone)]
pub struct SampleSet {
    /// `n` samples, each of length `dim`, in draw order.
    pub samples: Vec<Vec<f64>>,
    /// Primary (Cholesky) or fallback (eigen-clip) factorisation.
    pub method: FactorMethod,
}

/// Draw `n` samples `x ~ N(mean, cov)` reproducibly from `seed`.
///
/// Implements theory (U1): `x = mean + B z` with `B` the Cholesky factor
/// (primary) or eigen-clip reconstruction (fallback, reported in
/// [`SampleSet::method`]).
///
/// The same `(mean, cov, n, seed)` inputs always yield bit-identical
/// `samples`: `StdRng::seed_from_u64` (ChaCha, no OS entropy) drives a
/// Box–Muller normal generator, and the faer factorisation is deterministic.
/// `cov` is symmetrised by averaging with its transpose when already
/// symmetric within [`SYMMETRY_TOL_REL`]; larger asymmetries are
/// [`SampleError::NonSymmetric`] errors, never silent.
pub fn sample_mvn(
    mean: &[f64],
    cov: &[Vec<f64>],
    n: usize,
    seed: u64,
) -> Result<SampleSet, SampleError> {
    let dim = mean.len();
    if dim == 0 {
        return Err(SampleError::Empty);
    }
    if n == 0 {
        return Err(SampleError::NoSamples);
    }
    if cov.len() != dim {
        return Err(SampleError::DimensionMismatch {
            expected: dim,
            got: cov.len(),
        });
    }
    for (i, row) in cov.iter().enumerate() {
        if row.len() != dim {
            return Err(SampleError::NonSquare {
                dim,
                row: i,
                got: row.len(),
            });
        }
    }
    if mean.iter().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("mean"));
    }
    if cov.iter().flatten().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("covariance"));
    }

    // Symmetrise small asymmetries (SANDY instead requires exact symmetry
    // and raises `TypeError` otherwise); reject gross asymmetry.
    let scale = cov
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(0.0f64, f64::max)
        .max(1.0);
    let mut dev = 0.0f64;
    for (i, row) in cov.iter().enumerate() {
        for (j, v) in row.iter().enumerate().skip(i + 1) {
            dev = dev.max((v - cov[j][i]).abs());
        }
    }
    if dev > SYMMETRY_TOL_REL * scale {
        return Err(SampleError::NonSymmetric { deviation: dev });
    }

    let mat = Mat::<f64>::from_fn(dim, dim, |i, j| 0.5 * (cov[i][j] + cov[j][i]));

    // Primary path: Cholesky. Fallback: eigen-clipping.
    enum Factor {
        LowerTriangular(Vec<Vec<f64>>),
        Dense(Vec<Vec<f64>>),
    }
    let (factor, method) = match mat.as_ref().cholesky(Side::Lower) {
        Ok(chol) => {
            let l = chol.compute_l();
            let mut rows = vec![vec![0.0; dim]; dim];
            for i in 0..dim {
                for j in 0..=i {
                    rows[i][j] = l[(i, j)];
                }
            }
            (Factor::LowerTriangular(rows), FactorMethod::Cholesky)
        }
        Err(_) => {
            let eig = mat.as_ref().selfadjoint_eigendecomposition(Side::Lower);
            let evals: Vec<f64> = (0..dim).map(|i| eig.s().column_vector()[i]).collect();
            let max_eig = evals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let min_eigen = evals.iter().cloned().fold(f64::INFINITY, f64::min);
            if max_eig <= 0.0 {
                return Err(SampleError::NoPositiveEigenvalue);
            }
            let floor = EIGEN_FLOOR_REL * max_eig;
            let u = eig.u();
            let mut rows = vec![vec![0.0; dim]; dim];
            for i in 0..dim {
                for j in 0..dim {
                    rows[i][j] = u[(i, j)] * evals[j].max(floor).sqrt();
                }
            }
            (
                Factor::Dense(rows),
                FactorMethod::EigenClip {
                    min_eigen,
                    max_eigen: max_eig,
                },
            )
        }
    };

    let mut rng = StdRng::seed_from_u64(seed);
    let mut spare: Option<f64> = None;
    let mut normal = move || -> f64 {
        if let Some(z) = spare.take() {
            return z;
        }
        let u1 = rng.random::<f64>().max(f64::MIN_POSITIVE);
        let u2 = rng.random::<f64>();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        spare = Some(r * theta.sin());
        r * theta.cos()
    };

    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        let z: Vec<f64> = (0..dim).map(|_| normal()).collect();
        let mut x = vec![0.0; dim];
        match &factor {
            Factor::LowerTriangular(l) => {
                for i in 0..dim {
                    let mut acc = mean[i];
                    for j in 0..=i {
                        acc += l[i][j] * z[j];
                    }
                    x[i] = acc;
                }
            }
            Factor::Dense(b) => {
                for i in 0..dim {
                    let mut acc = mean[i];
                    for j in 0..dim {
                        acc += b[i][j] * z[j];
                    }
                    x[i] = acc;
                }
            }
        }
        samples.push(x);
    }
    Ok(SampleSet { samples, method })
}

/// Draw `n` log-normal samples reproducibly from `seed`.
///
/// Log-space semantics: draws `x ~ N(mean_log, cov_log)` with the shared
/// [`sample_mvn`] factor path (Cholesky primary, eigen-clip fallback) and
/// RNG stream, then maps `y = exp(x)` elementwise (`f64::exp` is the only
/// addition). `mean_log`/`cov_log` are the MVN parameters in log space —
/// never the mean/covariance of `y` itself. The closed-form moments of `y`
/// are `E[y_i] = exp(mu_i + C_ii/2)` and
/// `Cov(y_i, y_j) = exp(mu_i + mu_j + (C_ii + C_jj)/2) (exp(C_ij) - 1)`.
///
/// The same `(mean_log, cov_log, n, seed)` inputs always yield bit-identical
/// `samples`, and `y[k][i] == exp(x[k][i])` elementwise against
/// [`sample_mvn`] at the same inputs. All draws are strictly positive
/// (up to floating-point range); validation, symmetry handling, and the
/// reported [`FactorMethod`] are exactly the MVN ones.
pub fn sample_lognormal(
    mean_log: &[f64],
    cov_log: &[Vec<f64>],
    n: usize,
    seed: u64,
) -> Result<SampleSet, SampleError> {
    let mut set = sample_mvn(mean_log, cov_log, n, seed)?;
    for row in &mut set.samples {
        for v in row {
            *v = v.exp();
        }
    }
    Ok(set)
}

/// Closed-form mean of the log-normal draw `y = exp(x)`,
/// `x ~ N(mean_log, cov_log)`: `E[y_i] = exp(mu_i + C_ii/2)`.
///
/// Thin helper so callers and tests share one spelling of the theory (U6);
/// inputs are the same log-space parameters [`sample_lognormal`] takes.
pub fn lognormal_mean(mean_log: &[f64], cov_log: &[Vec<f64>]) -> Result<Vec<f64>, SampleError> {
    let dim = mean_log.len();
    if dim == 0 {
        return Err(SampleError::Empty);
    }
    if cov_log.len() != dim {
        return Err(SampleError::DimensionMismatch {
            expected: dim,
            got: cov_log.len(),
        });
    }
    for (i, row) in cov_log.iter().enumerate() {
        if row.len() != dim {
            return Err(SampleError::NonSquare {
                dim,
                row: i,
                got: row.len(),
            });
        }
    }
    if mean_log.iter().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("mean_log"));
    }
    if cov_log.iter().flatten().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("covariance_log"));
    }
    Ok(mean_log
        .iter()
        .enumerate()
        .map(|(i, mu)| (mu + 0.5 * cov_log[i][i]).exp())
        .collect())
}

/// Closed-form covariance of the log-normal draw `y = exp(x)`,
/// `x ~ N(mean_log, cov_log)`:
/// `Cov(y_i, y_j) = exp(mu_i + mu_j + (C_ii + C_jj)/2) (exp(C_ij) - 1)`.
///
/// Same validation as [`lognormal_mean`]; shares one spelling of theory (U6).
pub fn lognormal_cov(mean_log: &[f64], cov_log: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, SampleError> {
    let dim = mean_log.len();
    if dim == 0 {
        return Err(SampleError::Empty);
    }
    if cov_log.len() != dim {
        return Err(SampleError::DimensionMismatch {
            expected: dim,
            got: cov_log.len(),
        });
    }
    for (i, row) in cov_log.iter().enumerate() {
        if row.len() != dim {
            return Err(SampleError::NonSquare {
                dim,
                row: i,
                got: row.len(),
            });
        }
    }
    if mean_log.iter().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("mean_log"));
    }
    if cov_log.iter().flatten().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("covariance_log"));
    }
    let mut out = vec![vec![0.0; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            out[i][j] = (mean_log[i] + mean_log[j] + 0.5 * (cov_log[i][i] + cov_log[j][j])).exp()
                * (cov_log[i][j].exp() - 1.0);
        }
    }
    Ok(out)
}

/// Sample mean over draws (one entry per dimension).
pub fn sample_mean(samples: &[Vec<f64>]) -> Result<Vec<f64>, SampleError> {
    if samples.is_empty() {
        return Err(SampleError::Empty);
    }
    let dim = samples[0].len();
    if dim == 0 {
        return Err(SampleError::Empty);
    }
    let mut mean = vec![0.0; dim];
    for s in samples.iter() {
        if s.len() != dim {
            return Err(SampleError::DimensionMismatch {
                expected: dim,
                got: s.len(),
            });
        }
        if s.iter().any(|v| !v.is_finite()) {
            return Err(SampleError::NonFinite("samples"));
        }
        for (i, v) in s.iter().enumerate() {
            mean[i] += *v;
        }
    }
    let n = samples.len() as f64;
    for m in &mut mean {
        *m /= n;
    }
    Ok(mean)
}

/// Unbiased sample covariance (`1/(n-1)` denominator, matching SANDY
/// `Samples.get_cov` / pandas `DataFrame.cov`), as a `dim × dim` matrix.
pub fn sample_cov(samples: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, SampleError> {
    let mean = sample_mean(samples)?;
    let dim = mean.len();
    let n = samples.len();
    if n < 2 {
        return Err(SampleError::DimensionMismatch {
            expected: 2,
            got: n,
        });
    }
    let mut cov = vec![vec![0.0; dim]; dim];
    for s in samples {
        for i in 0..dim {
            for j in 0..dim {
                cov[i][j] += (s[i] - mean[i]) * (s[j] - mean[j]);
            }
        }
    }
    let denom = (n - 1) as f64;
    for row in &mut cov {
        for v in row {
            *v /= denom;
        }
    }
    Ok(cov)
}

/// Convergence report from [`check_convergence`].
#[derive(Debug, Clone, PartialEq)]
pub struct Convergence {
    /// Largest absolute deviation of the sample mean from `mean`.
    pub mean_err_max: f64,
    /// Frobenius norm of (`sample_cov` − `cov`).
    pub cov_err_fro: f64,
    /// Caller-supplied mean tolerance the report was checked against.
    pub mean_tol: f64,
    /// Caller-supplied covariance tolerance the report was checked against.
    pub cov_tol: f64,
    /// `mean_err_max <= mean_tol && cov_err_fro <= cov_tol`.
    pub passed: bool,
}

/// Sample mean/covariance convergence diagnostics: recompute the
/// moments of `samples` and compare against the `mean`/`cov` inputs.
///
/// Implements theory (U3) with caller-supplied tolerances, echoed back in
/// the report. The validation oracle pins statistical ones: for `n` MVN draws,
/// `std(mean_i) = sqrt(C[i][i]/n)` and
/// `var(S[i][j]) = (C[i][i]*C[j][j] + C[i][j]²)/(n-1)`, so honest gates are
/// `k` standard errors with a stated `k` — never hand-waved round numbers.
pub fn check_convergence(
    mean: &[f64],
    cov: &[Vec<f64>],
    samples: &[Vec<f64>],
    mean_tol: f64,
    cov_tol: f64,
) -> Result<Convergence, SampleError> {
    if !mean_tol.is_finite() || mean_tol < 0.0 {
        return Err(SampleError::NonFinite("mean_tol"));
    }
    if !cov_tol.is_finite() || cov_tol < 0.0 {
        return Err(SampleError::NonFinite("cov_tol"));
    }
    let dim = mean.len();
    if cov.len() != dim {
        return Err(SampleError::DimensionMismatch {
            expected: dim,
            got: cov.len(),
        });
    }
    let sm = sample_mean(samples)?;
    if sm.len() != dim {
        return Err(SampleError::DimensionMismatch {
            expected: dim,
            got: sm.len(),
        });
    }
    let sc = sample_cov(samples)?;
    let mean_err_max = sm
        .iter()
        .zip(mean.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    let mut acc = 0.0;
    for (srow, crow) in sc.iter().zip(cov.iter()) {
        if srow.len() != dim || crow.len() != dim {
            return Err(SampleError::DimensionMismatch {
                expected: dim,
                got: srow.len().max(crow.len()),
            });
        }
        for (a, b) in srow.iter().zip(crow.iter()) {
            acc += (a - b) * (a - b);
        }
    }
    let cov_err_fro = acc.sqrt();
    Ok(Convergence {
        mean_err_max,
        cov_err_fro,
        mean_tol,
        cov_tol,
        passed: mean_err_max <= mean_tol && cov_err_fro <= cov_tol,
    })
}

/// Apply one perturbation vector to one nominal vector under `convention`
/// (lengths must match; all values finite).
pub fn apply_perturbation(
    nominal: &[f64],
    delta: &[f64],
    convention: PerturbConvention,
) -> Result<Vec<f64>, SampleError> {
    if nominal.len() != delta.len() {
        return Err(SampleError::DimensionMismatch {
            expected: nominal.len(),
            got: delta.len(),
        });
    }
    if nominal.iter().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("nominal"));
    }
    if delta.iter().any(|v| !v.is_finite()) {
        return Err(SampleError::NonFinite("delta"));
    }
    Ok(nominal
        .iter()
        .zip(delta.iter())
        .map(|(x, d)| convention.apply(*x, *d))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEED: u64 = 20260913;

    fn cov_2x2() -> (Vec<f64>, Vec<Vec<f64>>) {
        // Synthetic block (fixtures/uq/cov_2x2.json): variances 0.25/0.16,
        // covariance 0.10 (correlation 0.5). No evaluated data anywhere.
        (vec![1.0, 2.0], vec![vec![0.25, 0.10], vec![0.10, 0.16]])
    }

    #[test]
    fn seeded_draws_are_reproducible() {
        let (mean, cov) = cov_2x2();
        let a = sample_mvn(&mean, &cov, 64, SEED).unwrap();
        let b = sample_mvn(&mean, &cov, 64, SEED).unwrap();
        assert_eq!(a.method, FactorMethod::Cholesky);
        assert_eq!(a.samples, b.samples);
        let c = sample_mvn(&mean, &cov, 64, SEED + 1).unwrap();
        assert_ne!(a.samples, c.samples);
    }

    #[test]
    fn moments_recover_closed_form_covariance() {
        // Honest statistical gate: for n MVN draws, std(mean_i) =
        // sqrt(C[i][i]/n) and var(S[i][j]) =
        // (C[i][i]*C[j][j] + C[i][j]^2)/(n-1); assert within k = 5
        // standard errors at the pinned seed (not a round-number guess).
        let (mean, cov) = cov_2x2();
        let n = 20_000usize;
        let set = sample_mvn(&mean, &cov, n, SEED).unwrap();
        let sm = sample_mean(&set.samples).unwrap();
        let sc = sample_cov(&set.samples).unwrap();
        let k = 5.0;
        for i in 0..2 {
            let se = (cov[i][i] / n as f64).sqrt();
            assert!(
                (sm[i] - mean[i]).abs() <= k * se,
                "mean[{i}] off by {:.3e} (> 5 SE = {:.3e})",
                (sm[i] - mean[i]).abs(),
                k * se
            );
            for j in 0..2 {
                let se_cov =
                    ((cov[i][i] * cov[j][j] + cov[i][j] * cov[i][j]) / (n - 1) as f64).sqrt();
                assert!(
                    (sc[i][j] - cov[i][j]).abs() <= k * se_cov,
                    "cov[{i}][{j}] off by {:.3e} (> 5 SE = {:.3e})",
                    (sc[i][j] - cov[i][j]).abs(),
                    k * se_cov
                );
            }
        }
    }

    #[test]
    fn singular_covariance_uses_eigen_clip_fallback() {
        // Rank-1 block: Cholesky must fail, eigen-clip takes over, and the
        // recovered covariance still matches within the same 5-SE gate.
        let mean = vec![0.0, 0.0];
        let cov = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let set = sample_mvn(&mean, &cov, 20_000, SEED).unwrap();
        match &set.method {
            FactorMethod::EigenClip {
                min_eigen,
                max_eigen,
            } => {
                assert!(min_eigen.abs() < 1e-9, "min_eigen = {min_eigen}");
                assert!((max_eigen - 2.0).abs() < 1e-9, "max_eigen = {max_eigen}");
            }
            FactorMethod::Cholesky => panic!("rank-1 block must not take the Cholesky path"),
        }
        let sc = sample_cov(&set.samples).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let se_cov = ((cov[i][i] * cov[j][j] + cov[i][j] * cov[i][j]) / 19_999.0).sqrt();
                assert!((sc[i][j] - cov[i][j]).abs() <= 5.0 * se_cov);
            }
        }
    }

    #[test]
    fn convergence_report_echoes_tolerances() {
        let (mean, cov) = cov_2x2();
        let set = sample_mvn(&mean, &cov, 512, SEED).unwrap();
        let rep = check_convergence(&mean, &cov, &set.samples, 1.0, 10.0).unwrap();
        assert!(rep.passed);
        assert_eq!((rep.mean_tol, rep.cov_tol), (1.0, 10.0));
        let strict = check_convergence(&mean, &cov, &set.samples, 1e-12, 1e-12).unwrap();
        assert!(!strict.passed);
    }

    #[test]
    fn conventions_apply_documented_identities() {
        let nominal = vec![2.0, 4.0];
        let delta = vec![0.5, -0.25];
        assert_eq!(
            apply_perturbation(&nominal, &delta, PerturbConvention::Relative).unwrap(),
            vec![3.0, 3.0]
        );
        assert_eq!(
            apply_perturbation(&nominal, &delta, PerturbConvention::Absolute).unwrap(),
            vec![2.5, 3.75]
        );
        // Log-normal: nominal * exp(delta); exp(0) is the identity.
        assert_eq!(
            apply_perturbation(&nominal, &[0.0, 0.0], PerturbConvention::LogNormal).unwrap(),
            nominal
        );
        let out = apply_perturbation(
            &[2.0],
            &[std::f64::consts::LN_2],
            PerturbConvention::LogNormal,
        )
        .unwrap();
        assert!((out[0] - 4.0).abs() < 1e-12, "lognormal apply = {out:?}");
        assert_eq!(
            PerturbConvention::parse("relative").unwrap(),
            PerturbConvention::Relative
        );
        assert_eq!(
            PerturbConvention::parse("ABS").unwrap(),
            PerturbConvention::Absolute
        );
        for spelling in [
            "lognormal",
            "log_normal",
            "log-normal",
            "LOGNORMAL",
            "lognorm",
        ] {
            assert_eq!(
                PerturbConvention::parse(spelling).unwrap(),
                PerturbConvention::LogNormal,
                "spelling {spelling}"
            );
        }
        assert!(PerturbConvention::parse("lhs").is_err());
    }

    #[test]
    fn lognormal_reuses_mvn_stream_then_exps() {
        // Same factor path + RNG as sample_mvn at the same seed; f64::exp is
        // the only addition. Synthetic log-space block (no evaluated data).
        let mean_log = vec![0.0, 0.0];
        let cov_log = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
        let mvn = sample_mvn(&mean_log, &cov_log, 64, SEED).unwrap();
        let logn = sample_lognormal(&mean_log, &cov_log, 64, SEED).unwrap();
        assert_eq!(mvn.method, logn.method);
        assert_eq!(mvn.method, FactorMethod::Cholesky);
        for (x, y) in mvn.samples.iter().zip(logn.samples.iter()) {
            for (a, b) in x.iter().zip(y.iter()) {
                assert_eq!(*b, a.exp());
                assert!(b.is_finite() && *b > 0.0);
            }
        }
        // Seeded reproducibility.
        let again = sample_lognormal(&mean_log, &cov_log, 64, SEED).unwrap();
        assert_eq!(logn.samples, again.samples);
        let other = sample_lognormal(&mean_log, &cov_log, 64, SEED + 1).unwrap();
        assert_ne!(logn.samples, other.samples);
    }

    #[test]
    fn lognormal_moments_match_closed_form() {
        // Honest gates on the committed synthetic fixture
        // (fixtures/uq/lognormal_2x2.json): the log draws ln(y) recover the
        // log-space inputs within the MVN k-SE gate (exact, since y = exp(x)
        // is a bijection), and the sample mean of y recovers the closed-form
        // E[y_i] = exp(mu_i + C_ii/2) within k standard errors with
        // Var(y_i) = exp(2 mu_i + C_ii)(exp(C_ii) - 1) from theory (U6).
        let text = include_str!("../../../fixtures/uq/lognormal_2x2.json");
        let v: serde_json::Value = serde_json::from_str(text).unwrap();
        let mean_log: Vec<f64> = serde_json::from_value(v["mean_log"].clone()).unwrap();
        let cov_log: Vec<Vec<f64>> = serde_json::from_value(v["cov"].clone()).unwrap();
        let seed = v["seed"].as_u64().unwrap();
        let n = v["n"].as_u64().unwrap() as usize;
        let k = v["k"].as_f64().unwrap();
        let dim = mean_log.len();
        let set = sample_lognormal(&mean_log, &cov_log, n, seed).unwrap();
        assert_eq!(set.method, FactorMethod::Cholesky);
        assert!(set.samples.iter().all(|s| s.iter().all(|x| *x > 0.0)));
        // L1: ln(y) == MVN draws recover the log-space block.
        let ln: Vec<Vec<f64>> = set
            .samples
            .iter()
            .map(|s| s.iter().map(|x| x.ln()).collect())
            .collect();
        let sm = sample_mean(&ln).unwrap();
        let sc = sample_cov(&ln).unwrap();
        for i in 0..dim {
            let se = (cov_log[i][i] / n as f64).sqrt();
            assert!(
                (sm[i] - mean_log[i]).abs() <= k * se,
                "ln mean[{i}] off by {:.3e} (> {k} SE)",
                (sm[i] - mean_log[i]).abs()
            );
            for j in 0..dim {
                let se_cov = ((cov_log[i][i] * cov_log[j][j] + cov_log[i][j] * cov_log[i][j])
                    / (n - 1) as f64)
                    .sqrt();
                assert!(
                    (sc[i][j] - cov_log[i][j]).abs() <= k * se_cov,
                    "ln cov[{i}][{j}] gate failed"
                );
            }
        }
        // L2: sample mean of y against the closed-form log-normal mean.
        let expected = lognormal_mean(&mean_log, &cov_log).unwrap();
        let cov_y = lognormal_cov(&mean_log, &cov_log).unwrap();
        let sm_y = sample_mean(&set.samples).unwrap();
        for i in 0..dim {
            let se = (cov_y[i][i] / n as f64).sqrt();
            assert!(
                (sm_y[i] - expected[i]).abs() <= k * se,
                "lognormal mean[{i}] off by {:.3e} (> {k} SE = {:.3e})",
                (sm_y[i] - expected[i]).abs(),
                k * se
            );
        }
    }

    #[test]
    fn lognormal_closed_form_helpers_match_theory() {
        // Hand-computed spot: mu = [0, 0], C = [[0.04, 0.01], [0.01, 0.09]].
        let mu = vec![0.0, 0.0];
        let cov = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
        let m = lognormal_mean(&mu, &cov).unwrap();
        assert!((m[0] - 0.02f64.exp()).abs() < 1e-15);
        assert!((m[1] - 0.045f64.exp()).abs() < 1e-15);
        let c = lognormal_cov(&mu, &cov).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let want =
                    (mu[i] + mu[j] + 0.5 * (cov[i][i] + cov[j][j])).exp() * (cov[i][j].exp() - 1.0);
                assert!((c[i][j] - want).abs() < 1e-15);
            }
        }
        assert!(lognormal_mean(&[], &[]).is_err());
        assert!(lognormal_cov(&mu, &[vec![1.0]]).is_err());
    }

    #[test]
    fn malformed_inputs_name_their_cause() {
        let (mean, cov) = cov_2x2();
        assert_eq!(
            sample_mvn(&[], &[], 4, SEED).unwrap_err(),
            SampleError::Empty
        );
        assert_eq!(
            sample_mvn(&mean, &cov, 0, SEED).unwrap_err(),
            SampleError::NoSamples
        );
        assert!(matches!(
            sample_mvn(&mean, &[vec![1.0]], 4, SEED).unwrap_err(),
            SampleError::DimensionMismatch { .. }
        ));
        assert!(matches!(
            sample_mvn(&mean, &[vec![1.0], vec![1.0, 2.0]], 4, SEED).unwrap_err(),
            SampleError::NonSquare { .. }
        ));
        assert_eq!(
            sample_mvn(&[f64::NAN, 0.0], &cov, 4, SEED).unwrap_err(),
            SampleError::NonFinite("mean")
        );
        let asym = vec![vec![1.0, 0.0], vec![0.5, 1.0]];
        assert!(matches!(
            sample_mvn(&mean, &asym, 4, SEED).unwrap_err(),
            SampleError::NonSymmetric { .. }
        ));
        let neg = vec![vec![-1.0, 0.0], vec![0.0, -2.0]];
        assert_eq!(
            sample_mvn(&mean, &neg, 4, SEED).unwrap_err(),
            SampleError::NoPositiveEigenvalue
        );
        assert!(sample_cov(&sample_mvn(&mean, &cov, 1, SEED).unwrap().samples).is_err());
    }

    #[test]
    fn fixture_blocks_recover_within_statistical_gates() {
        // Committed synthetic fixtures (no evaluated data): parse the exact
        // files the Python tests and validation oracle read, and apply the
        // same k-SE gates. Fixture schema: {mean, cov, seed, n, k}.
        for text in [
            include_str!("../../../fixtures/uq/cov_2x2.json"),
            include_str!("../../../fixtures/uq/cov_3x3.json"),
        ] {
            let v: serde_json::Value = serde_json::from_str(text).unwrap();
            let mean: Vec<f64> = serde_json::from_value(v["mean"].clone()).unwrap();
            let cov: Vec<Vec<f64>> = serde_json::from_value(v["cov"].clone()).unwrap();
            let seed = v["seed"].as_u64().unwrap();
            let n = v["n"].as_u64().unwrap() as usize;
            let k = v["k"].as_f64().unwrap();
            let dim = mean.len();
            let set = sample_mvn(&mean, &cov, n, seed).unwrap();
            assert_eq!(set.method, FactorMethod::Cholesky);
            let sm = sample_mean(&set.samples).unwrap();
            let sc = sample_cov(&set.samples).unwrap();
            for i in 0..dim {
                let se = (cov[i][i] / n as f64).sqrt();
                assert!((sm[i] - mean[i]).abs() <= k * se, "mean[{i}] gate failed");
                for j in 0..dim {
                    let se_cov =
                        ((cov[i][i] * cov[j][j] + cov[i][j] * cov[i][j]) / (n - 1) as f64).sqrt();
                    assert!(
                        (sc[i][j] - cov[i][j]).abs() <= k * se_cov,
                        "cov[{i}][{j}] gate failed"
                    );
                }
            }
        }
    }
}
