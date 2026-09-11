//! Error type for the kinetics crate.

use thiserror::Error;

/// Errors raised while validating kinetic data, reactivity schedules, time
/// grids, or solver options — and when a requested analysis has no solution
/// (prompt-supercritical prompt jump, non-delayed-supercritical period).
#[derive(Debug, Clone, PartialEq, Error)]
pub enum Error {
    /// No precursor groups supplied (`betas`/`lambdas` empty).
    #[error("kinetics: at least one precursor group is required")]
    EmptyGroups,
    /// `betas` and `lambdas` lengths differ (`{betas}` vs `{lambdas}`).
    #[error("kinetics: betas ({betas}) and lambdas ({lambdas}) lengths differ")]
    LengthMismatch {
        /// Length of `betas`.
        betas: usize,
        /// Length of `lambdas`.
        lambdas: usize,
    },
    /// A kinetic datum is non-finite or out of range: `{0}`.
    #[error("kinetics: invalid kinetic data: {0}")]
    BadData(&'static str),
    /// A reactivity schedule parameter is invalid: `{0}`.
    #[error("kinetics: invalid reactivity schedule: {0}")]
    BadReactivity(&'static str),
    /// A time-grid entry is invalid: `{0}`.
    #[error("kinetics: invalid time grid: {0}")]
    BadGrid(&'static str),
    /// A solver option is invalid: `{0}`.
    #[error("kinetics: invalid solver option: {0}")]
    BadOption(&'static str),
    /// An initial state entry is invalid: `{0}`.
    #[error("kinetics: invalid initial state: {0}")]
    BadState(&'static str),
    /// The prompt-jump formula needs `rho_after < beta_total`; got
    /// `rho_after={rho_after}` with `beta_total={beta_total}`.
    #[error("kinetics: prompt jump needs rho_after ({rho_after}) < beta_total ({beta_total})")]
    PromptSupercritical {
        /// Reactivity after the step \[Δk\].
        rho_after: f64,
        /// Total delayed-neutron fraction.
        beta_total: f64,
    },
    /// The stable-period solve needs `0 < rho < beta_total`; got
    /// `rho={rho}` with `beta_total={beta_total}`.
    #[error("kinetics: stable period needs 0 < rho ({rho}) < beta_total ({beta_total})")]
    NoStablePeriod {
        /// Requested reactivity \[Δk\].
        rho: f64,
        /// Total delayed-neutron fraction.
        beta_total: f64,
    },
    /// The inhour evaluation point is at or beyond a pole
    /// (`omega={omega}` with slowest decay `lambda_min={lambda_min}`).
    #[error("kinetics: inhour pole crossed at omega={omega} (lambda_min={lambda_min})")]
    InhourPole {
        /// Angular root candidate [1/s].
        omega: f64,
        /// Slowest precursor decay constant [1/s].
        lambda_min: f64,
    },
    /// The adaptive integrator exceeded `{0}` accepted steps.
    #[error("kinetics: step budget exhausted ({0} steps)")]
    StepBudget(usize),
}
