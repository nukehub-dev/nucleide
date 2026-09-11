//! Stiff-aware PKE integrator, prompt-jump approximation (E4), initial rate.
//!
//! The system `dy/dt = A(t) y` with `y = [n, C_1..C_G]` is stiff: the
//! prompt-neutron timescale `Λ/(β−ρ)` is microseconds against precursor
//! times of seconds. [`solve`] therefore uses an **adaptive implicit**
//! scheme — trapezoidal (second order, A-stable; default) or backward Euler
//! (first order, L-stable) — with step-doubling error control, a hard step
//! budget, and exact stepping onto reactivity knots and output times. No
//! explicit `dopri5`-style parity is attempted (see the crate docs).
//!
//! The prompt-jump approximation (E4) follows from (E1) in the `Λ → 0`
//! limit with precursors frozen across a step `ρ⁻ → ρ⁺`: prompt equilibrium
//! `0 = (ρ−β)n + ΛS` before and after (same source `S = ΣλC`) gives
//!
//! ```text
//! n⁺ = n⁻ (β − ρ⁻) / (β − ρ⁺),
//! ```
//!
//! valid for `ρ⁺ < β`. It is an approximation (documented ~percent-level
//! gate), not an algebraic identity like (E2)/(E3).

use crate::error::Error;
use crate::params::KineticParams;
use crate::reactivity::Reactivity;

/// Initial state: neutron level plus precursor populations.
///
/// `None` precursors default to the equilibrium (E2) populations for `n0`.
#[derive(Debug, Clone, PartialEq)]
pub struct State {
    /// Neutron population (arbitrary units, `>= 0`).
    pub n0: f64,
    /// Precursor populations (length `G`, `>= 0`, finite).
    pub c0: Vec<f64>,
}

impl State {
    /// Build a state, defaulting precursors to equilibrium (E2) when `c0`
    /// is `None`. Rejects non-finite/negative `n0` and mismatched or
    /// non-finite/negative `c0`.
    pub fn new(params: &KineticParams, n0: f64, c0: Option<Vec<f64>>) -> Result<Self, Error> {
        if !n0.is_finite() || n0 < 0.0 {
            return Err(Error::BadState("n0 must be finite and >= 0"));
        }
        let c0 = match c0 {
            Some(c) => {
                if c.len() != params.groups() {
                    return Err(Error::BadState("C0 length must match group count"));
                }
                if c.iter().any(|x| !x.is_finite() || *x < 0.0) {
                    return Err(Error::BadState("C0 must be finite and >= 0"));
                }
                c
            }
            None => params.equilibrium_precursors(n0)?,
        };
        Ok(Self { n0, c0 })
    }
}

/// Output time grid: strictly increasing times `> 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeGrid {
    /// Output times \[s\].
    pub times: Vec<f64>,
}

impl TimeGrid {
    /// Validate: at least one time, all finite and `> 0`, strictly
    /// increasing.
    pub fn new(times: Vec<f64>) -> Result<Self, Error> {
        if times.is_empty() {
            return Err(Error::BadGrid("time grid must not be empty"));
        }
        let mut prev = 0.0_f64;
        for t in &times {
            if !t.is_finite() || *t <= prev {
                return Err(Error::BadGrid(
                    "times must be finite, > 0, strictly increasing",
                ));
            }
            prev = *t;
        }
        Ok(Self { times })
    }
}

/// Implicit integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    /// Second-order A-stable trapezoidal rule (default).
    #[default]
    Trapezoidal,
    /// First-order L-stable backward Euler (extra damping on the prompt
    /// transient; larger phase error on the tail).
    BackwardEuler,
}

/// Solver tolerances and budgets.
#[derive(Debug, Clone, PartialEq)]
pub struct SolverOptions {
    /// Integration method.
    pub method: Method,
    /// Relative tolerance for step-doubling control.
    pub rtol: f64,
    /// Absolute tolerance (in `n`/`C` units) for step-doubling control.
    pub atol: f64,
    /// Smallest allowed internal step \[s\].
    pub dt_min: f64,
    /// Largest allowed internal step \[s\].
    pub dt_max: f64,
    /// Maximum accepted internal steps over the whole grid.
    pub max_steps: usize,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            method: Method::Trapezoidal,
            rtol: 1e-9,
            atol: 1e-12,
            dt_min: 1e-14,
            dt_max: f64::INFINITY,
            max_steps: 1_000_000,
        }
    }
}

impl SolverOptions {
    /// Validate tolerances and budgets (finite, positive, ordered).
    pub fn validate(&self) -> Result<(), Error> {
        if !self.rtol.is_finite() || self.rtol <= 0.0 {
            return Err(Error::BadOption("rtol must be finite and > 0"));
        }
        if !self.atol.is_finite() || self.atol <= 0.0 {
            return Err(Error::BadOption("atol must be finite and > 0"));
        }
        if !self.dt_min.is_finite() || self.dt_min <= 0.0 {
            return Err(Error::BadOption("dt_min must be finite and > 0"));
        }
        if self.dt_max <= 0.0 || self.dt_max < self.dt_min {
            return Err(Error::BadOption("need 0 < dt_min <= dt_max"));
        }
        if self.max_steps == 0 {
            return Err(Error::BadOption("max_steps must be > 0"));
        }
        Ok(())
    }
}

/// Full transient: one row per output time (the `t = 0` state is the
/// caller-supplied [`State`], echoed here as `initial`).
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Output times \[s\] (echo of the grid).
    pub times: Vec<f64>,
    /// Neutron population at each output time.
    pub n: Vec<f64>,
    /// Precursor populations at each output time (`[time][group]`).
    pub c: Vec<Vec<f64>>,
    /// The `t = 0` state.
    pub initial: State,
}

/// Right-hand side of (E1) at reactivity `rho`: returns `(dn/dt, dC/dt)`.
///
/// This is the discrete form of the upstream `dpdt`/`dzetadt` pair; the O1
/// gate asserts `dn/dt = n0·ρ(0)/Λ` at equilibrium initials to 1e-12.
pub fn rhs(params: &KineticParams, rho: f64, n: f64, c: &[f64]) -> (f64, Vec<f64>) {
    let beta = params.beta_total();
    let gen = params.lambda_gen();
    let mut precursors = 0.0;
    for (l, ci) in params.lambdas().iter().zip(c) {
        precursors += l * ci;
    }
    let dn = n * (rho - beta) / gen + precursors;
    let dc = params
        .betas()
        .iter()
        .zip(params.lambdas())
        .zip(c)
        .map(|((b, l), ci)| b * n / gen - l * ci)
        .collect();
    (dn, dc)
}

/// Initial rate `dn/dt|₀` at equilibrium initials.
///
/// Algebraically `n0·ρ(0)/Λ` (the precursor source cancels the
/// `(ρ−β)` term exactly); computed here through [`rhs`] so the O1 gate
/// pins the implementation, not the formula.
pub fn initial_rate(params: &KineticParams, rho: &Reactivity, state: &State) -> f64 {
    rhs(params, rho.eval(0.0), state.n0, &state.c0).0
}

/// Prompt-jump estimate (E4) across a step `rho_before → rho_after`:
///
/// ```text
/// n⁺ = n⁻ (β − ρ⁻) / (β − ρ⁺).
/// ```
///
/// Errors when `ρ⁺ >= β` (no finite prompt equilibrium past prompt
/// critical). Approximate: exact only in the `Λ → 0` limit with frozen
/// precursors.
pub fn prompt_jump(
    n_before: f64,
    rho_before: f64,
    rho_after: f64,
    beta_total: f64,
) -> Result<f64, Error> {
    if !n_before.is_finite() || n_before < 0.0 {
        return Err(Error::BadState("n_before must be finite and >= 0"));
    }
    if !rho_before.is_finite() || !rho_after.is_finite() || !beta_total.is_finite() {
        return Err(Error::BadData("reactivities and beta must be finite"));
    }
    if rho_after >= beta_total {
        return Err(Error::PromptSupercritical {
            rho_after,
            beta_total,
        });
    }
    Ok(n_before * (beta_total - rho_before) / (beta_total - rho_after))
}

/// Solve the PKE system over `grid` under `rho` from `state`.
///
/// Knots of `rho` inside `(0, t_end]` and every output time are hit
/// exactly; between consecutive event times the adaptive implicit stepper
/// advances with step-doubling control under `opts`. Returns [`Solution`]
/// with one row per grid time.
pub fn solve(
    params: &KineticParams,
    rho: &Reactivity,
    grid: &TimeGrid,
    state: &State,
    opts: &SolverOptions,
) -> Result<Solution, Error> {
    rho.validate()?;
    opts.validate()?;
    if state.c0.len() != params.groups() {
        return Err(Error::BadState("C0 length must match group count"));
    }
    let g = params.groups();
    let dim = g + 1;
    let t_end = grid.times[grid.times.len() - 1];

    // Event times: knots inside (0, t_end] plus output times; dedupe exact
    // duplicates (a step edge exactly on an output time is one event).
    let mut events: Vec<f64> = Vec::new();
    for k in rho.knots() {
        if k > 0.0 && k <= t_end {
            events.push(k);
        }
    }
    events.extend(grid.times.iter().cloned());
    events.sort_by(|a, b| a.total_cmp(b));
    events.dedup();

    let mut y = vec![0.0; dim];
    y[0] = state.n0;
    y[1..].copy_from_slice(&state.c0);

    let mut n_out = Vec::with_capacity(grid.times.len());
    let mut c_out = Vec::with_capacity(grid.times.len());
    let mut gi = 0;
    let mut t = 0.0_f64;
    let mut h = (events[0]).min(opts.dt_max).max(opts.dt_min);
    let mut steps = 0_usize;

    for &t_ev in &events {
        while t < t_ev {
            let h_try = (h.min(opts.dt_max).max(opts.dt_min)).min(t_ev - t);
            if h_try < opts.dt_min && t_ev - t > 0.0 {
                // Final sliver smaller than dt_min: take it exactly rather
                // than stalling (error already controlled on prior steps).
                y = implicit_step(params, rho, t, t_ev - t, &y, opts.method);
                t = t_ev;
                steps += 1;
                h = h_try;
                break;
            }
            let (y_new, err) = doubled_step(params, rho, t, h_try, &y, opts);
            if err <= 1.0 {
                y = y_new;
                t += h_try;
                steps += 1;
                if steps > opts.max_steps {
                    return Err(Error::StepBudget(opts.max_steps));
                }
            }
            let order = match opts.method {
                Method::Trapezoidal => 2.0,
                Method::BackwardEuler => 1.0,
            };
            let factor = 0.9 * err.max(1e-16).powf(-1.0 / (order + 1.0));
            h = h_try * factor.clamp(0.2, 5.0);
            if err > 1.0 {
                continue;
            }
            if t >= t_ev {
                break;
            }
        }
        // Record output rows whose grid time is (numerically) this event.
        while gi < grid.times.len() && (grid.times[gi] - t_ev).abs() <= 1e-12 * t_ev.max(1.0) {
            n_out.push(y[0]);
            c_out.push(y[1..].to_vec());
            gi += 1;
        }
    }

    debug_assert_eq!(gi, grid.times.len());
    Ok(Solution {
        times: grid.times.clone(),
        n: n_out,
        c: c_out,
        initial: state.clone(),
    })
}

/// One implicit step of size `h` from `(t, y)` with `rho` sampled at the
/// interval ends. Solves `(I − h·θ·A₁) y₁ = y₀ + h·(1−θ)·A₀y₀ + h·θ·... `
/// in the standard θ-form (`θ = 1/2` trapezoidal, `θ = 1` backward Euler)
/// via small-dense partial-pivot elimination (`dim = G+1 ≤ 9`, no
/// factorization reuse needed).
fn implicit_step(
    params: &KineticParams,
    rho: &Reactivity,
    t: f64,
    h: f64,
    y: &[f64],
    method: Method,
) -> Vec<f64> {
    let theta = match method {
        Method::Trapezoidal => 0.5,
        Method::BackwardEuler => 1.0,
    };
    let dim = y.len();
    let a0 = system_matrix(params, rho.eval(t));
    let a1 = system_matrix(params, rho.eval(t + h));
    // M = I − h·θ·A₁; b = (I + h·(1−θ)·A₀) y.
    let mut m = vec![vec![0.0; dim]; dim];
    let mut b = vec![0.0; dim];
    for i in 0..dim {
        for j in 0..dim {
            m[i][j] = -h * theta * a1[i][j];
        }
        m[i][i] += 1.0;
        let mut bi = y[i];
        if theta < 1.0 {
            let mut ay = 0.0;
            for j in 0..dim {
                ay += a0[i][j] * y[j];
            }
            bi += h * (1.0 - theta) * ay;
        }
        b[i] = bi;
    }
    solve_dense(&mut m, &mut b)
}

/// Step-doubling pair: one full step vs two half steps; returns the
/// higher-accuracy (half-step) state plus the scaled max-norm error.
fn doubled_step(
    params: &KineticParams,
    rho: &Reactivity,
    t: f64,
    h: f64,
    y: &[f64],
    opts: &SolverOptions,
) -> (Vec<f64>, f64) {
    let full = implicit_step(params, rho, t, h, y, opts.method);
    let half = implicit_step(params, rho, t, h * 0.5, y, opts.method);
    let two = implicit_step(params, rho, t + h * 0.5, h * 0.5, &half, opts.method);
    let mut err = 0.0_f64;
    for i in 0..y.len() {
        let scale = opts.atol + opts.rtol * full[i].abs().max(two[i].abs());
        err = err.max((two[i] - full[i]).abs() / scale);
    }
    (two, err)
}

/// PKE system matrix `A(ρ)` in `[n, C]` ordering.
fn system_matrix(params: &KineticParams, rho: f64) -> Vec<Vec<f64>> {
    let g = params.groups();
    let dim = g + 1;
    let mut a = vec![vec![0.0; dim]; dim];
    let beta = params.beta_total();
    let gen = params.lambda_gen();
    a[0][0] = (rho - beta) / gen;
    for (i, (b, l)) in params.betas().iter().zip(params.lambdas()).enumerate() {
        a[0][1 + i] = *l;
        a[1 + i][0] = b / gen;
        a[1 + i][1 + i] = -*l;
    }
    a
}

/// Dense solve `M x = b` by Gaussian elimination with partial pivoting.
/// `M` is `(G+1)² ≤ 81` entries; singularity cannot occur for the
/// well-posed implicit step (`I − hθA` with `h > 0`), so near-zero pivots
/// fall back to a tiny regularization instead of an error.
fn solve_dense(m: &mut [Vec<f64>], b: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let mut piv = k;
        for i in (k + 1)..n {
            if m[i][k].abs() > m[piv][k].abs() {
                piv = i;
            }
        }
        m.swap(k, piv);
        b.swap(k, piv);
        if m[k][k].abs() < 1e-300 {
            m[k][k] = 1e-300;
        }
        for i in (k + 1)..n {
            let f = m[i][k] / m[k][k];
            m[i][k] = 0.0;
            let pivot_row = m[k].clone();
            for (j, mij) in m[i].iter_mut().enumerate().skip(k + 1) {
                *mij -= f * pivot_row[j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::tests::six_group;

    fn grid(times: Vec<f64>) -> TimeGrid {
        TimeGrid::new(times).unwrap()
    }

    fn opts() -> SolverOptions {
        SolverOptions {
            rtol: 1e-10,
            atol: 1e-14,
            ..Default::default()
        }
    }

    #[test]
    fn o1_initial_rate_matches_n0_rho_over_lambda() {
        // O1 algebraic gate @1e-12: at equilibrium initials the precursor
        // source cancels the (ρ−β) term, leaving n0·ρ/Λ.
        let p = six_group();
        let s = State::new(&p, 1.0, None).unwrap();
        // ρ(0) = ρ_final via a constant schedule (a step at t > 0 reads
        // rho_init at t = 0; covered by the stationary test below).
        let rc = Reactivity::Constant { rho: 0.002 };
        let got = initial_rate(&p, &rc, &s);
        let want = 1.0 * 0.002 / 1e-5;
        assert!((got - want).abs() / want < 1e-12, "got {got} want {want}");
    }

    #[test]
    fn o1_zero_reactivity_is_stationary() {
        // Invariant: ρ = 0 at equilibrium initials stays put.
        let p = six_group();
        let rho = Reactivity::zero();
        let s = State::new(&p, 2.5, None).unwrap();
        let sol = solve(&p, &rho, &grid(vec![0.1, 1.0, 10.0, 100.0]), &s, &opts()).unwrap();
        for n in &sol.n {
            assert!((n - 2.5).abs() / 2.5 < 1e-9, "n={n}");
        }
    }

    /// Closed-form 1-group step response (independent path from the
    /// solver): n'' + [λ + (β−ρ)/Λ] n' − λρ/Λ n = 0 with n(0) = n0,
    /// n'(0) = n0 ρ/Λ.
    fn analytic_one_group(beta: f64, lam: f64, gen: f64, rho: f64, n0: f64, t: f64) -> f64 {
        let a = lam + (beta - rho) / gen;
        let disc = (a * a + 4.0 * lam * rho / gen).sqrt();
        let s1 = 0.5 * (-a + disc);
        let s2 = 0.5 * (-a - disc);
        let np = n0 * rho / gen;
        let big_a = (np - s2 * n0) / (s1 - s2);
        let big_b = n0 - big_a;
        big_a * (s1 * t).exp() + big_b * (s2 * t).exp()
    }

    #[test]
    fn o3_one_group_step_matches_analytic_1e_6() {
        let p = KineticParams::new(vec![0.0065], vec![0.08], 1e-4).unwrap();
        let rho = Reactivity::Step {
            t_step: 1.0,
            rho_init: 0.0,
            rho_final: 0.002,
        };
        let s = State::new(&p, 1.0, None).unwrap();
        let times: Vec<f64> = vec![0.5, 1.0, 1.0001, 1.01, 1.1, 2.0, 5.0, 20.0, 100.0];
        let sol = solve(&p, &rho, &grid(times.clone()), &s, &opts()).unwrap();
        for (t, n) in times.iter().zip(&sol.n) {
            let r = if *t < 1.0 { 0.0 } else { 0.002 };
            // Before the step the equilibrium state is stationary.
            let want = if *t < 1.0 {
                1.0
            } else {
                analytic_one_group(0.0065, 0.08, 1e-4, r, 1.0, t - 1.0)
            };
            let err = (n - want).abs() / want.max(1e-300);
            assert!(err < 1e-6, "t={t} got={n} want={want} err={err}");
        }
    }

    #[test]
    fn o2_prompt_jump_within_two_percent() {
        // Step 0 → 0.003 (β = 0.0065): prompt plateau ≈ β/(β−ρ) = 1.857.
        let p = six_group();
        let rho = Reactivity::Step {
            t_step: 1.0,
            rho_init: 0.0,
            rho_final: 0.003,
        };
        let s = State::new(&p, 1.0, None).unwrap();
        let pj = prompt_jump(1.0, 0.0, 0.003, p.beta_total()).unwrap();
        assert!((pj - 0.0065 / (0.0065 - 0.003)).abs() < 1e-12);
        // Sample just after the prompt transient (Λ/(β−ρ) ≈ 2.9 ms) but
        // long before precursors move: solver must sit on the plateau.
        let sol = solve(&p, &rho, &grid(vec![1.05]), &s, &opts()).unwrap();
        let err = (sol.n[0] - pj).abs() / pj;
        assert!(err < 0.02, "plateau={} jump={pj} err={err}", sol.n[0]);
    }

    #[test]
    fn o4_six_group_tail_matches_stable_period() {
        // Late-time log slope over [50, 100] s must equal the dominant
        // inhour root within 2%.
        let p = six_group();
        let rho = Reactivity::Step {
            t_step: 0.0,
            rho_init: 0.002,
            rho_final: 0.002,
        };
        let s = State::new(&p, 1.0, None).unwrap();
        let sol = solve(&p, &rho, &grid(vec![50.0, 60.0, 80.0, 100.0]), &s, &opts()).unwrap();
        let slope = (sol.n[3] / sol.n[0]).ln() / 50.0;
        let t_period = crate::stable_period(&p, 0.002).unwrap();
        let err = (slope - 1.0 / t_period).abs() / (1.0 / t_period);
        assert!(
            err < 0.02,
            "slope={slope} omega={} err={err}",
            1.0 / t_period
        );
    }

    #[test]
    fn positivity_and_monotone_step_up() {
        let p = six_group();
        let rho = Reactivity::Step {
            t_step: 1.0,
            rho_init: 0.0,
            rho_final: 0.002,
        };
        let s = State::new(&p, 1.0, None).unwrap();
        let times: Vec<f64> = (0..50).map(|k| 0.2 * (k + 1) as f64).collect();
        let sol = solve(&p, &rho, &grid(times), &s, &opts()).unwrap();
        for (n, c) in sol.n.iter().zip(&sol.c) {
            assert!(*n > 0.0);
            assert!(c.iter().all(|x| *x > 0.0));
        }
        for w in sol.n.windows(2) {
            // Tolerance-aware: the pre-step stationary segment holds only
            // to stepper tolerance.
            assert!(w[1] >= w[0] - 1e-9, "positive step must not decrease");
        }
    }

    #[test]
    fn backward_euler_agrees_on_tail() {
        let p = six_group();
        let rho = Reactivity::Constant { rho: 0.002 };
        let s = State::new(&p, 1.0, None).unwrap();
        let g = grid(vec![10.0, 50.0, 100.0]);
        let a = solve(&p, &rho, &g, &s, &opts()).unwrap();
        let b = solve(
            &p,
            &rho,
            &g,
            &s,
            &SolverOptions {
                method: Method::BackwardEuler,
                ..opts()
            },
        )
        .unwrap();
        for (x, y) in a.n.iter().zip(&b.n) {
            assert!((x - y).abs() / x < 1e-3, "{x} vs {y}");
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        let p = six_group();
        let s = State::new(&p, 1.0, None).unwrap();
        assert!(TimeGrid::new(vec![]).is_err());
        assert!(TimeGrid::new(vec![1.0, 1.0]).is_err());
        assert!(TimeGrid::new(vec![-1.0]).is_err());
        assert!(State::new(&p, -1.0, None).is_err());
        assert!(State::new(&p, 1.0, Some(vec![1.0])).is_err());
        assert!(prompt_jump(1.0, 0.0, 0.0065, 0.0065).is_err());
        assert!(prompt_jump(1.0, 0.0, 0.01, 0.0065).is_err());
        let bad_opts = SolverOptions {
            rtol: -1.0,
            ..opts()
        };
        assert!(solve(&p, &Reactivity::zero(), &grid(vec![1.0]), &s, &bad_opts).is_err());
        let tight = SolverOptions {
            dt_min: 1.0,
            dt_max: 0.5,
            ..opts()
        };
        assert!(solve(&p, &Reactivity::zero(), &grid(vec![1.0]), &s, &tight).is_err());
    }

    #[test]
    fn rejects_bad_state_options_and_jumps() {
        let p = six_group();
        // Precursor vectors must match the group count and stay non-negative.
        let explicit = State::new(&p, 1.0, Some(vec![1.0; 6])).unwrap();
        assert_eq!(explicit.c0, vec![1.0; 6]);
        assert!(State::new(
            &p,
            1.0,
            Some(vec![1.0; 6].into_iter().map(|_| f64::NAN).collect())
        )
        .is_err());
        assert!(State::new(&p, 1.0, Some(vec![-1.0; 6])).is_err());
        assert!(State::new(&p, f64::INFINITY, None).is_err());
        // Every solver budget must be finite and positive.
        for bad in [
            SolverOptions {
                atol: -1.0,
                ..opts()
            },
            SolverOptions {
                atol: f64::NAN,
                ..opts()
            },
            SolverOptions {
                dt_min: 0.0,
                ..opts()
            },
            SolverOptions {
                dt_min: f64::NAN,
                ..opts()
            },
            SolverOptions {
                max_steps: 0,
                ..opts()
            },
        ] {
            assert!(bad.validate().is_err());
        }
        // Prompt-jump inputs must be finite with a non-negative population.
        assert!(prompt_jump(-1.0, 0.0, 0.001, 0.0065).is_err());
        assert!(prompt_jump(f64::NAN, 0.0, 0.001, 0.0065).is_err());
        assert!(prompt_jump(1.0, f64::NAN, 0.001, 0.0065).is_err());
        assert!(prompt_jump(1.0, 0.0, 0.001, f64::INFINITY).is_err());
        // A hand-built state whose precursors miss the group count is
        // rejected at solve time even though it type-checks.
        let s = State::new(&p, 1.0, None).unwrap();
        let mismatched = State {
            n0: s.n0,
            c0: vec![1.0],
        };
        assert!(solve(
            &p,
            &Reactivity::zero(),
            &grid(vec![1.0]),
            &mismatched,
            &opts()
        )
        .is_err());
    }

    #[test]
    fn solver_takes_final_sliver_and_enforces_budget() {
        let p = six_group();
        let s = State::new(&p, 1.0, None).unwrap();
        // A dt_min larger than the output span takes the sliver path.
        let sliver = SolverOptions {
            dt_min: 1.0,
            ..opts()
        };
        let sol = solve(&p, &Reactivity::zero(), &grid(vec![0.5]), &s, &sliver).unwrap();
        assert!((sol.n[0] - 1.0).abs() < 1e-6);
        // One accepted step per output time: a budget of one fails on two.
        let budgeted = SolverOptions {
            max_steps: 1,
            ..opts()
        };
        assert!(solve(
            &p,
            &Reactivity::zero(),
            &grid(vec![1.0, 2.0]),
            &s,
            &budgeted
        )
        .is_err());
    }
}
