//! Analytic Bateman decay fast path with an f64-careful (HP) variant.
//!
//! Textbook clean-room implementation (Bateman 1910; Amaku–Pascholati–Vanin,
//! Comput. Phys. Commun. 181 (2010) 21–23, eigendecomposition form) over the
//! crate's own chain data. No third-party data or code is involved.
//!
//! # Equations
//!
//! For decay-only systems the burnup matrix is lower-triangular `Λ` with
//!
//! ```text
//! (E1)  dN/dt = Λ N,
//!       Λ_jj = -λ_j,  Λ_ij = b_j→i · λ_j  (i > j),  λ = ln 2 / T½ (stable → 0).
//! ```
//!
//! `Λ` is assembled exactly like [`crate::matrix::DepletionSystem`] decay
//! entries (spontaneous-fission daughters skipped, `He4`/`H1` light-particle
//! secondaries from alpha/proton decay included, zero branches skipped).
//!
//! For distinct diagonal entries on an acyclic (lower-triangular) chain,
//! `Λ = C·D·C⁻¹` with `D = diag(-λ)` and unit-diagonal `C`:
//!
//! ```text
//! (E2)  C_ii = 1,  C_ij = 0 (i < j),
//! (E3)  C_ij = Σ_{k=j}^{i-1} Λ_ik · C_kj / (Λ_jj − Λ_ii)   (i > j),
//! (E4)  C⁻¹_ii = 1,  C⁻¹_ij = 0 (i < j),
//!       C⁻¹_ij = −Σ_{k=j}^{i-1} C_ik · C⁻¹_kj               (i > j),
//! (E5)  N(t) = C · diag(e^{−λt}) · C⁻¹ · N₀.
//! ```
//!
//! `C`/`C⁻¹` are time-independent, so they are built once per chain
//! ([`BatemanCache`]) and reused across timesteps. The solve itself factors
//! as two triangular matvecs (`y = C⁻¹·N₀`, `z = e∘y`, `N = C·z`).
//!
//! Observables reuse the crate-wide channels:
//!
//! ```text
//! (E8)  A_i = λ_i · N_i,
//!       ∫₀ᵗ e^{−λs} ds = (1 − e^{−λt}) / λ,  λ = 0 limit → t.
//! ```
//!
//! Activity is [`crate::integrate::activity_vec`]; the diagonal cumulative
//! integral is [`crate::inventory::cumulative_decays`]
//! (`exp_m1` form, stable nuclides report `0.0`).
//!
//! # Precision variants
//!
//! [`Method::Bateman`] evaluates (E5) with plain index-order accumulation.
//! [`Method::BatemanHp`] evaluates the same closed form with f64-careful
//! summation only (no `f128`, no new dependencies): per-row terms sorted by
//! magnitude and summed with Neumaier compensation, plus `exp_m1` forms for
//! small `λt` (same precedent as `inventory.rs` cumulative decays).
//!
//! # Fallback rules (D1–D4)
//!
//! The fast path is attempted only for decay-only, lower-triangular systems
//! with usable eigenvalue gaps. Anything else falls back to CRAM-48 inside
//! [`solve_with_method`] (never an error — the caller cannot tell other
//! than by timing):
//!
//! - (D1) near-equal decay constants with a live coupling path
//!   (`|λⱼ − λᵢ| ≤ tol · max(λᵢ, λⱼ)`, `tol = 1e-12`, numerator `≠ 0`)
//!   → CRAM-48. Pairs whose numerator vanishes (e.g. two stable nuclides
//!   with no decay path between them) take the limit form `C_ij = 0`
//!   instead of falling back, so ordinary multi-stable chains stay on the
//!   fast path.
//! - (D2) cyclic or out-of-order decay topology (any decay gain with
//!   `row < col`, i.e. not lower-triangular) → CRAM-48.
//! - (D3) stable nuclides (`λ = 0`) → limit forms handled inline
//!   (`e⁰ = 1`, `C_ij = 0` for uncoupled pairs); no fallback.
//! - (D4) reactions on (assembled matrix differs from the pure-decay
//!   reconstruction: fission yields, capture/transmutation gains or losses)
//!   → CRAM-48, since the system is no longer triangular decay-only.
//!
//! A zero timestep (`dt == 0.0`) returns the input vector exactly.

use std::collections::BTreeMap;
use std::str::FromStr;

use crate::cram::{Error as CramError, Order};
use crate::matrix::DepletionSystem;

/// Solver selector for single-step and series depletion.
///
/// `Order` itself is intentionally untouched (it is matched as `16 | 48` in
/// the Python/WASM binding layers); the CRAM order rides inside
/// [`Method::Cram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// IPF CRAM at the given order (production default is `Cram(Order48)`).
    Cram(Order),
    /// Analytic Bateman closed form (E5), plain f64 accumulation.
    Bateman,
    /// Same closed form with f64-careful summation (sorted terms, Neumaier
    /// compensation, `exp_m1` for small `λt`).
    BatemanHp,
}

impl Method {
    /// Production default: CRAM-48.
    pub fn default_cram() -> Self {
        Method::Cram(Order::Order48)
    }

    /// Whether this method uses the Bateman fast path (standard or HP).
    pub fn is_bateman(self) -> bool {
        matches!(self, Method::Bateman | Method::BatemanHp)
    }
}

impl FromStr for Method {
    type Err = String;

    /// Parse method spellings (case-insensitive, `-`/`_` interchangeable):
    /// `cram`, `cram16`, `cram48`, `bateman`, `bateman_hp`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let norm: String = s.trim().to_ascii_lowercase().replace(['-', ' '], "_");
        match norm.as_str() {
            "cram" | "cram48" | "cram_48" => Ok(Method::Cram(Order::Order48)),
            "cram16" | "cram_16" => Ok(Method::Cram(Order::Order16)),
            "bateman" | "bateman_standard" => Ok(Method::Bateman),
            "bateman_hp" | "batemanhp" | "bateman_high_precision" | "hp" => Ok(Method::BatemanHp),
            other => Err(format!(
                "unsupported method `{s}` (supported: cram16, cram48, bateman, bateman_hp); got `{other}`"
            )),
        }
    }
}

impl std::fmt::Display for Method {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Method::Cram(Order::Order16) => write!(f, "cram16"),
            Method::Cram(Order::Order48) => write!(f, "cram48"),
            Method::Bateman => write!(f, "bateman"),
            Method::BatemanHp => write!(f, "bateman_hp"),
        }
    }
}

/// Why a system cannot use the Bateman fast path (routes to CRAM-48).
#[derive(Debug, Clone, PartialEq)]
pub enum FallbackReason {
    /// Near-equal decay constants with a live coupling path (D1).
    NearDegenerate {
        /// Row index (daughter side).
        i: usize,
        /// Column index (ancestor side).
        j: usize,
        /// `|λⱼ − λᵢ|`.
        gap: f64,
    },
    /// Decay gain above the diagonal: cyclic or out-of-order topology (D2).
    NonTriangular {
        /// Offending gain row.
        row: usize,
        /// Offending gain column.
        col: usize,
    },
    /// Assembled matrix differs from the pure-decay reconstruction:
    /// reactions/fission active (D4).
    ReactionsPresent,
}

impl std::fmt::Display for FallbackReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FallbackReason::NearDegenerate { i, j, gap } => write!(
                f,
                "near-equal decay constants between chain indices {j} and {i} (gap {gap:e}); CRAM fallback"
            ),
            FallbackReason::NonTriangular { row, col } => write!(
                f,
                "non-triangular decay edge ({col} -> {row}); CRAM fallback"
            ),
            FallbackReason::ReactionsPresent => {
                write!(f, "reaction/fission contributions present; CRAM fallback")
            }
        }
    }
}

/// Relative gap tolerance for (D1) near-degenerate detection.
pub const DEGENERATE_TOL: f64 = 1e-12;

/// Time-independent Bateman eigendecomposition for one chain.
///
/// Built once via [`BatemanCache::build`] (which enforces D1/D2/D4), then
/// reused across timesteps via [`BatemanCache::solve`]. The cache holds the
/// same `C`/`C⁻¹` for both precision variants; `hp` only changes the
/// summation in the triangular matvecs.
#[derive(Debug, Clone)]
pub struct BatemanCache {
    n: usize,
    lambdas: Vec<f64>,
    c: Vec<Vec<f64>>,
    cinv: Vec<Vec<f64>>,
}

impl BatemanCache {
    /// Build `C`/`C⁻¹` for a decay-only lower-triangular system.
    ///
    /// Returns [`FallbackReason`] when the fast path does not apply (the
    /// caller falls back to CRAM-48). Half-life defects cannot occur here:
    /// [`DepletionSystem::build`] already rejects them, so a stable `0.0`
    /// is used if one ever slipped through.
    ///
    /// Index loops below are intentional: every recurrence addresses two or
    /// more dense triangles by `(row, col)`, which iterators obscure.
    #[allow(clippy::needless_range_loop)]
    pub fn build(sys: &DepletionSystem) -> Result<Self, FallbackReason> {
        let n = sys.chain.len();
        let mut lambdas = Vec::with_capacity(n);
        for nuc in &sys.chain.nuclides {
            lambdas.push(nuc.decay_constant().unwrap_or(0.0));
        }

        // Pure-decay reconstruction (mirrors matrix.rs decay section).
        let decay = decay_entries(&sys.chain);
        // D4: the assembled matrix must equal the decay-only reconstruction.
        if !matrix_matches_decay(sys, &decay) {
            return Err(FallbackReason::ReactionsPresent);
        }
        // D2: every decay gain must sit below the diagonal.
        for &(r, c) in decay.keys() {
            if r < c {
                return Err(FallbackReason::NonTriangular { row: r, col: c });
            }
        }

        // Dense lower-triangular Λ in chain order.
        let mut lam = vec![vec![0.0; n]; n];
        for (&(r, c), &v) in &decay {
            lam[r][c] += v;
        }

        // (E2–E3): unit-diagonal C, column by column, rows ascending so every
        // C_kj (k < i) needed by the numerator is already computed.
        let mut c = vec![vec![0.0; n]; n];
        for i in 0..n {
            c[i][i] = 1.0;
        }
        for j in 0..n {
            for i in (j + 1)..n {
                let denom = lam[j][j] - lam[i][i]; // Λ_jj − Λ_ii = λ_i − λ_j
                let mut num = 0.0;
                for k in j..i {
                    num += lam[i][k] * c[k][j];
                }
                if num == 0.0 {
                    // Uncoupled pair: limit form C_ij = 0 (covers D3
                    // stable–stable exact degeneracy with no path).
                    c[i][j] = 0.0;
                    continue;
                }
                let scale = DEGENERATE_TOL * lambdas[i].max(lambdas[j]);
                if denom.abs() <= scale {
                    return Err(FallbackReason::NearDegenerate {
                        i,
                        j,
                        gap: denom.abs(),
                    });
                }
                c[i][j] = num / denom;
            }
        }

        // (E4): C⁻¹ needs no division (unit diagonal).
        let mut cinv = vec![vec![0.0; n]; n];
        for i in 0..n {
            cinv[i][i] = 1.0;
        }
        for j in 0..n {
            for i in (j + 1)..n {
                let mut acc = 0.0;
                for k in j..i {
                    acc += c[i][k] * cinv[k][j];
                }
                cinv[i][j] = -acc;
            }
        }

        Ok(Self {
            n,
            lambdas,
            c,
            cinv,
        })
    }

    /// Number of nuclides.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Whether the cache holds no nuclides.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Decay constants in chain order [1/s] (stable → 0.0).
    pub fn lambdas(&self) -> &[f64] {
        &self.lambdas
    }

    /// Solve `N(dt) = C·diag(e^{−λdt})·C⁻¹·N₀` (E5).
    ///
    /// `dt == 0.0` returns the input exactly. Non-positive (other than
    /// exactly zero), non-finite `dt`, length mismatches, and non-finite or
    /// negative `n0` entries are errors like the CRAM entry points.
    ///
    /// Index loops below are intentional (paired dense-triangle addressing).
    #[allow(clippy::needless_range_loop)]
    pub fn solve(&self, n0: &[f64], dt: f64, hp: bool) -> Result<Vec<f64>, CramError> {
        if n0.len() != self.n {
            return Err(CramError::Linalg(format!(
                "n0 length {} != chain size {}",
                n0.len(),
                self.n
            )));
        }
        if dt == 0.0 {
            return Ok(n0.to_vec());
        }
        if dt <= 0.0 || !dt.is_finite() {
            return Err(CramError::Linalg(format!("invalid timestep dt: {dt}")));
        }
        if n0.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(CramError::Linalg(
                "n0 must hold finite atom counts >= 0".to_string(),
            ));
        }
        // y = C⁻¹·N₀ (forward substitution; C⁻¹_ii = 1).
        let mut y = vec![0.0; self.n];
        for i in 0..self.n {
            y[i] = if hp {
                neumaier_row(&self.cinv[i][..=i], n0, i)
            } else {
                let mut acc = n0[i];
                for k in 0..i {
                    acc += self.cinv[i][k] * n0[k];
                }
                acc
            };
        }
        // z = e^{−λdt} ∘ y (HP uses exp_m1 for small arguments).
        for i in 0..self.n {
            let x = -self.lambdas[i] * dt;
            let e = if hp && x.abs() < 0.5 {
                // 1 + expm1(x): keeps the O(x) part for tiny x.
                1.0 + x.exp_m1()
            } else {
                x.exp()
            };
            y[i] *= e;
        }
        // N = C·z (C_ii = 1).
        let mut out = vec![0.0; self.n];
        for i in 0..self.n {
            out[i] = if hp {
                neumaier_row(&self.c[i][..=i], &y, i)
            } else {
                let mut acc = y[i];
                for k in 0..i {
                    acc += self.c[i][k] * y[k];
                }
                acc
            };
        }
        Ok(out)
    }
}

/// Neumaier compensated sum of `row[k] * x[k]` over `k = 0..=i`, with terms
/// accumulated largest-magnitude-first (the HP summation kernel).
fn neumaier_row(row: &[f64], x: &[f64], i: usize) -> f64 {
    let mut terms: Vec<f64> = (0..=i).map(|k| row[k] * x[k]).collect();
    terms.sort_by(|a, b| b.abs().total_cmp(&a.abs()));
    let mut sum = 0.0;
    let mut comp = 0.0;
    for t in terms {
        let u = sum + t;
        comp += if sum.abs() >= t.abs() {
            (sum - u) + t
        } else {
            (t - u) + sum
        };
        sum = u;
    }
    sum + comp
}

/// Pure-decay `(row, col) → value` entries, replicating the decay section of
/// `DepletionSystem::build` (matrix.rs): diagonal loss, non-sf daughter
/// gains, `He4`/`H1` secondaries from alpha/proton decay, zero branches
/// skipped, every diagonal present.
fn decay_entries(chain: &crate::chain::Chain) -> BTreeMap<(usize, usize), f64> {
    let mut acc: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut add = |r: usize, c: usize, v: f64| {
        *acc.entry((r, c)).or_insert(0.0) += v;
    };
    for (i, nuc) in chain.nuclides.iter().enumerate() {
        let lambda = nuc.decay_constant().unwrap_or(0.0);
        if lambda > 0.0 {
            add(i, i, -lambda);
            for mode in &nuc.decay_modes {
                let branch_val = lambda * mode.branching_ratio;
                if branch_val == 0.0 {
                    continue;
                }
                if !mode.kind.contains("sf") {
                    if let Some(j) = chain.index_of(&mode.target) {
                        add(j, i, branch_val);
                    }
                }
                if mode.kind.contains("alpha") {
                    if let Some(j) = chain.index_of("He4") {
                        let count = mode.kind.matches("alpha").count();
                        add(j, i, count as f64 * branch_val);
                    }
                } else if mode.kind.contains('p') {
                    if let Some(j) = chain.index_of("H1") {
                        let count = mode.kind.matches('p').count();
                        add(j, i, count as f64 * branch_val);
                    }
                }
            }
        }
    }
    for i in 0..chain.len() {
        acc.entry((i, i)).or_insert(0.0);
    }
    acc
}

/// Whether the assembled system matrix equals the pure-decay reconstruction
/// (D4 check): identical sparsity pattern plus matching values.
///
/// Values use a tight relative tolerance (`1e-12`, absolute floor `1e-30`)
/// rather than bit equality so the check is robust to accumulation order;
/// any genuine reaction/fission contribution dwarfs this band.
fn matrix_matches_decay(sys: &DepletionSystem, decay: &BTreeMap<(usize, usize), f64>) -> bool {
    if sys.entries.len() != decay.len() {
        return false;
    }
    for (entry, base) in sys.entries.iter().zip(sys.base_values.iter()) {
        let key = (entry.row, entry.col);
        let Some(&want) = decay.get(&key) else {
            return false;
        };
        let got = base.re;
        let scale = want.abs().max(got.abs()).max(1e-30);
        if (got - want).abs() > 1e-12 * scale {
            return false;
        }
    }
    true
}

/// Solve one depletion step with an explicit [`Method`].
///
/// `Cram(order)` runs the IPF CRAM path; `Bateman`/`BatemanHp` build the
/// cached closed form and solve (E5), falling back to CRAM-48 on D1/D2/D4
/// systems (including any system whose matrix carries reaction/fission
/// contributions).
pub fn solve_with_method(
    sys: &DepletionSystem,
    method: Method,
    n0: &[f64],
    dt: f64,
) -> Result<Vec<f64>, CramError> {
    match method {
        Method::Cram(order) => crate::cram(sys, order, n0, dt),
        Method::Bateman => match BatemanCache::build(sys) {
            Ok(cache) => cache.solve(n0, dt, false),
            Err(_) => crate::cram(sys, Order::Order48, n0, dt),
        },
        Method::BatemanHp => match BatemanCache::build(sys) {
            Ok(cache) => cache.solve(n0, dt, true),
            Err(_) => crate::cram(sys, Order::Order48, n0, dt),
        },
    }
}

/// [`solve_with_method`] with a caller-managed symbolic factorization (used
/// by the CRAM arms; the Bateman arms bypass LU entirely).
pub fn solve_with_method_symbolic(
    sys: &DepletionSystem,
    sym: &nucleide_linalg::SymbolicLu,
    method: Method,
    n0: &[f64],
    dt: f64,
) -> Result<Vec<f64>, CramError> {
    match method {
        Method::Cram(order) => crate::cram_with_symbolic(sys, sym, order, n0, dt),
        Method::Bateman => match BatemanCache::build(sys) {
            Ok(cache) => cache.solve(n0, dt, false),
            Err(_) => crate::cram_with_symbolic(sys, sym, Order::Order48, n0, dt),
        },
        Method::BatemanHp => match BatemanCache::build(sys) {
            Ok(cache) => cache.solve(n0, dt, true),
            Err(_) => crate::cram_with_symbolic(sys, sym, Order::Order48, n0, dt),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::{Chain, ChainNuclide, DecayMode};
    use crate::matrix::ReactionRates;

    const L1: f64 = 1.0e-6;
    const L2: f64 = 1.0e-5;

    fn hl(lam: f64) -> f64 {
        std::f64::consts::LN_2 / lam
    }

    fn abc_chain() -> Chain {
        Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(L1)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                half_life: Some(hl(L2)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "C".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "C".into(),
                ..Default::default()
            },
        ])
        .unwrap()
    }

    fn bateman_abc(n0: f64, t: f64) -> [f64; 3] {
        let na = n0 * (-L1 * t).exp();
        let nb = n0 * L1 / (L2 - L1) * ((-L1 * t).exp() - (-L2 * t).exp());
        [na, nb, n0 - na - nb]
    }

    #[test]
    fn method_spellings() {
        assert_eq!("cram48".parse(), Ok(Method::Cram(Order::Order48)));
        assert_eq!("cram16".parse(), Ok(Method::Cram(Order::Order16)));
        assert_eq!("cram".parse(), Ok(Method::Cram(Order::Order48)));
        assert_eq!("bateman".parse(), Ok(Method::Bateman));
        assert_eq!("bateman_hp".parse(), Ok(Method::BatemanHp));
        assert_eq!("Bateman-HP".parse(), Ok(Method::BatemanHp));
        assert!("cram32".parse::<Method>().is_err());
        assert_eq!(Method::Bateman.to_string(), "bateman");
        assert!(Method::BatemanHp.is_bateman());
        assert!(!Method::Cram(Order::Order48).is_bateman());
    }

    #[test]
    fn standard_matches_analytic() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let n0 = vec![1.0e15, 0.0, 0.0];
        let dt = 1.0e5;
        for hp in [false, true] {
            let got = BatemanCache::build(&sys)
                .unwrap()
                .solve(&n0, dt, hp)
                .unwrap();
            let want = bateman_abc(n0[0], dt);
            assert!(
                (got[0] - want[0]).abs() / want[0] < 1e-8,
                "{got:?} {want:?}"
            );
            assert!(
                (got[1] - want[1]).abs() / want[1] < 1e-7,
                "{got:?} {want:?}"
            );
            assert!(
                (got[2] - want[2]).abs() / want[2] < 1e-9,
                "{got:?} {want:?}"
            );
            let total: f64 = got.iter().sum();
            assert!((total - n0[0]).abs() / n0[0] < 1e-8);
        }
    }

    #[test]
    fn cache_reuse_matches_fresh() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        for dt in [1.0, 2.0e4, 1.0e5, 5.0e5] {
            let n0 = vec![1.0e14, 5e13, 1e10];
            let a = cache.solve(&n0, dt, false).unwrap();
            let b = BatemanCache::build(&sys)
                .unwrap()
                .solve(&n0, dt, false)
                .unwrap();
            assert_eq!(a, b);
            let h = cache.solve(&n0, dt, true).unwrap();
            for (x, y) in a.iter().zip(&h) {
                assert!((x - y).abs() / x.max(1e-30) < 1e-12, "{x} vs {y}");
            }
        }
    }

    #[test]
    fn zero_dt_returns_input_exactly() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        let n0 = vec![1.0e14, 5e13, 1e10];
        assert_eq!(cache.solve(&n0, 0.0, false).unwrap(), n0);
        assert_eq!(cache.solve(&n0, 0.0, true).unwrap(), n0);
    }

    #[test]
    fn solve_rejects_bad_dt_and_bad_n0() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        let good = vec![1.0e14, 5e13, 1e10];
        // dt=0 echoes; every other non-positive or non-finite dt errors
        // (inf/NaN already errored before this change).
        assert_eq!(cache.solve(&good, 0.0, false).unwrap(), good);
        for dt in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(cache.solve(&good, dt, false).is_err(), "dt={dt}");
            assert!(cache.solve(&good, dt, true).is_err(), "dt={dt}");
        }
        // NaN/infinite/negative n0 entries error with the Linalg style.
        for bad in [
            vec![f64::NAN, 0.0, 0.0],
            vec![1.0e14, f64::INFINITY, 0.0],
            vec![1.0e14, 5e13, -1.0],
        ] {
            let err = cache.solve(&bad, 1.0e5, false).unwrap_err();
            assert!(err.to_string().contains("n0"), "{err}");
            assert!(cache.solve(&bad, 1.0e5, true).is_err());
        }
        // Length mismatch still errors; valid input still solves.
        assert!(cache.solve(&[1.0], 1.0e5, false).is_err());
        assert!(cache.solve(&good, 1.0e5, false).is_ok());
    }

    #[test]
    fn wide_spread_chain_hp_matches_and_conserves() {
        // Es-254-style spread: By-scale parent, µs-scale progeny (10 orders).
        let lam_p = std::f64::consts::LN_2 / 3.15576e7; // ~1 y parent
        let lam_d = std::f64::consts::LN_2 / 1.0e-6; // 1 µs daughter
        let chain = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "P".into(),
                half_life: Some(std::f64::consts::LN_2 / lam_p),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "D".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "D".into(),
                half_life: Some(std::f64::consts::LN_2 / lam_d),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "S".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "S".into(),
                ..Default::default()
            },
        ])
        .unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        let n0 = vec![1.0e15, 0.0, 0.0];
        let dt = 1.0e6;
        let std = cache.solve(&n0, dt, false).unwrap();
        let hp = cache.solve(&n0, dt, true).unwrap();
        // Analytic reference (daughter in secular equilibrium, then decayed).
        let exp_p = (-lam_p * dt).exp();
        let exp_d = (-lam_d * dt).exp();
        let want_p = n0[0] * exp_p;
        let want_d = n0[0] * lam_p / (lam_d - lam_p) * (exp_p - exp_d);
        let want_s = n0[0] - want_p - want_d;
        for (got, want) in hp.iter().zip([want_p, want_d, want_s]) {
            let scale = want.abs().max(1e-30);
            assert!((got - want).abs() / scale < 1e-8, "{got} vs {want}");
        }
        // HP beats-or-matches standard against the analytic reference.
        let err = |v: &[f64]| {
            v.iter()
                .zip([want_p, want_d, want_s])
                .map(|(g, w)| (g - w).abs() / w.abs().max(1e-30))
                .fold(0.0_f64, f64::max)
        };
        assert!(err(&hp) <= err(&std), "{} vs {}", err(&hp), err(&std));
        for v in [&std, &hp] {
            let total: f64 = v.iter().sum();
            assert!((total - n0[0]).abs() / n0[0] < 1e-8, "{v:?}");
        }
    }

    #[test]
    fn bateman_vs_cram48_within_series_band() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let n0 = vec![1.0e15, 0.0, 0.0];
        for dt in [1.0e4, 1.0e5, 5.0e5] {
            let cram = crate::cram(&sys, Order::Order48, &n0, dt).unwrap();
            for method in [Method::Bateman, Method::BatemanHp] {
                let got = solve_with_method(&sys, method, &n0, dt).unwrap();
                for (g, c) in got.iter().zip(&cram) {
                    assert!((g - c).abs() / c.max(1e-30) < 1e-6, "{g} vs {c}");
                }
            }
        }
    }

    #[test]
    fn d1_near_equal_falls_back_to_cram() {
        // Relative λ gap 1e-13 < tol 1e-12 with a live A→B edge.
        let lam = 1.0e-6;
        let chain = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(lam)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                half_life: Some(hl(lam * (1.0 + 1e-13))),
                ..Default::default()
            },
        ])
        .unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        assert!(matches!(
            BatemanCache::build(&sys),
            Err(FallbackReason::NearDegenerate { .. })
        ));
        // Dispatch still succeeds via CRAM fallback.
        let n0 = vec![1.0e12, 0.0];
        let got = solve_with_method(&sys, Method::Bateman, &n0, 1.0e5).unwrap();
        let want = crate::cram(&sys, Order::Order48, &n0, 1.0e5).unwrap();
        assert_eq!(got, want);
    }

    #[test]
    fn d2_cyclic_falls_back() {
        let chain = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(L1)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                half_life: Some(hl(L2)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "A".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
        ])
        .unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        assert!(matches!(
            BatemanCache::build(&sys),
            Err(FallbackReason::NonTriangular { .. })
        ));
    }

    #[test]
    fn d3_stable_limit_forms_stay_fast() {
        // Two stables + one decaying parent feeding one of them: no D1 trip.
        let chain = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(L1)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                ..Default::default()
            },
            ChainNuclide {
                name: "S".into(),
                ..Default::default()
            },
        ])
        .unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        let n0 = vec![1.0e12, 0.0, 3.0e11];
        let dt = 1.0e5;
        let got = cache.solve(&n0, dt, false).unwrap();
        assert!((got[0] - n0[0] * (-L1 * dt).exp()).abs() / got[0] < 1e-12);
        assert!((got[1] - (n0[0] - got[0])).abs() / got[1].max(1e-30) < 1e-9);
        assert_eq!(got[2], n0[2]);
    }

    #[test]
    fn d4_reactions_fall_back() {
        let sys_decay = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        assert!(BatemanCache::build(&sys_decay).is_ok());
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-7);
        // Give A a capture channel so the rate is live.
        let chain = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(L1)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                reactions: vec![crate::chain::Reaction {
                    kind: "(n,gamma)".into(),
                    target: Some("B".into()),
                    q: 0.0,
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                half_life: Some(hl(L2)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "C".into(),
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "C".into(),
                ..Default::default()
            },
        ])
        .unwrap();
        let sys_rx = DepletionSystem::build(chain, &rates).unwrap();
        assert!(matches!(
            BatemanCache::build(&sys_rx),
            Err(FallbackReason::ReactionsPresent)
        ));
        let n0 = vec![1.0e12, 0.0, 0.0];
        let got = solve_with_method(&sys_rx, Method::BatemanHp, &n0, 1.0e5).unwrap();
        let want = crate::cram(&sys_rx, Order::Order48, &n0, 1.0e5).unwrap();
        assert_eq!(got, want);
    }

    #[test]
    fn invalid_inputs_error() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let cache = BatemanCache::build(&sys).unwrap();
        assert!(cache.solve(&[1.0], 1.0, false).is_err());
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(cache.solve(&[1.0, 0.0, 0.0], bad, false).is_err());
        }
    }

    #[test]
    fn symbolic_dispatch_matches_plain() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let sym = nucleide_linalg::SymbolicLu::try_new(&sys.pattern).unwrap();
        let n0 = vec![1.0e15, 0.0, 0.0];
        let dt = 1.0e5;
        for method in [Method::Bateman, Method::BatemanHp] {
            let a = solve_with_method(&sys, method, &n0, dt).unwrap();
            let b = solve_with_method_symbolic(&sys, &sym, method, &n0, dt).unwrap();
            assert_eq!(a, b);
        }
        // Fallback arm through the symbolic entry point is bit-identical
        // (live capture rate on A makes the matrix non-decay).
        let chain_rx = Chain::from_nuclides(vec![
            ChainNuclide {
                name: "A".into(),
                half_life: Some(hl(L1)),
                decay_modes: vec![DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 1.0,
                }],
                reactions: vec![crate::chain::Reaction {
                    kind: "(n,gamma)".into(),
                    target: Some("B".into()),
                    q: 0.0,
                    branching_ratio: 1.0,
                }],
                ..Default::default()
            },
            ChainNuclide {
                name: "B".into(),
                half_life: Some(hl(L2)),
                ..Default::default()
            },
        ])
        .unwrap();
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-7);
        let sys_rx = DepletionSystem::build(chain_rx, &rates).unwrap();
        assert!(BatemanCache::build(&sys_rx).is_err());
        let sym_rx = nucleide_linalg::SymbolicLu::try_new(&sys_rx.pattern).unwrap();
        let n0_rx = vec![1.0e12, 0.0];
        let a = solve_with_method(&sys_rx, Method::Bateman, &n0_rx, dt).unwrap();
        let b = solve_with_method_symbolic(&sys_rx, &sym_rx, Method::Bateman, &n0_rx, dt).unwrap();
        assert_eq!(a, b);
        // Fallback reasons display without panicking.
        for r in [
            FallbackReason::NearDegenerate {
                i: 1,
                j: 0,
                gap: 1e-19,
            },
            FallbackReason::NonTriangular { row: 0, col: 1 },
            FallbackReason::ReactionsPresent,
        ] {
            assert!(format!("{r}").contains("CRAM fallback"));
        }
        assert!(!BatemanCache::build(&sys).unwrap().is_empty());
        assert_eq!(BatemanCache::build(&sys).unwrap().len(), 3);
        assert_eq!(BatemanCache::build(&sys).unwrap().lambdas().len(), 3);
    }
}
