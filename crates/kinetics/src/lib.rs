#![warn(missing_docs)]
//! Point-kinetics transients with prescribed reactivity.
//!
//! Solves the point-kinetics equations (E1) for the neutron population `n`
//! and `G` delayed-neutron precursor groups `C_i` under a caller-supplied
//! reactivity schedule `rho(t)` in units of Δk:
//!
//! ```text
//! dn/dt   = (rho - beta) / Lambda * n + sum_i lambda_i * C_i      (E1a)
//! dC_i/dt = beta_i / Lambda * n - lambda_i * C_i                  (E1b)
//! ```
//!
//! Equation set (E1) reproduces the PyRK neutronics block
//! (`Neutronics::dpdt` / `Neutronics::dzetadt` in the upstream BSD-3
//! `pyrk/pyrk` package): same state vector, same right-hand sides, same Δk
//! reactivity units. The equilibrium initial conditions (E2),
//!
//! ```text
//! C_i(0) = beta_i / (lambda_i * Lambda) * n0,
//! ```
//!
//! reproduce the PyRK driver initial conditions (`y0` in the upstream
//! `driver` module: power normalized, precursors at the steady-state ratio).
//!
//! What is intentionally different from PyRK:
//!
//! - The time integrator is an adaptive implicit trapezoidal scheme
//!   (A-stable, exact stepping to reactivity knots) instead of the explicit
//!   `dopri5` driver loop. The PKE system is stiff (prompt-neutron lifetime
//!   `Lambda ~ 1e-5 s` against precursor decay times `1/lambda_i ~ 0.1-100
//!   s`); an explicit scheme needs unphysically small steps on the prompt
//!   timescale, so no `dopri5` parity is attempted.
//! - Thermal-hydraulic feedback, the decay-heat block, tabulated precursor
//!   data, and input-file plumbing are out of scope. All kinetic data
//!   (`betas`, `lambdas`, `Lambda`) are caller-supplied; see
//!   [`KineticParams::from_ifp`] for the OpenMC-IFP provenance note.
//!
//! Modules: [`params`] (data + E2), [`reactivity`] (insertion taxonomy),
//! [`mod@solve`] (stiff-aware integrator + E4 prompt jump), [`inhour`] (E3),
//! [`error`] (error type).

pub mod error;
pub mod inhour;
pub mod params;
pub mod reactivity;
pub mod solve;

pub use error::Error;
pub use inhour::{residual as inhour_residual, rho_of_omega, stable_period};
pub use params::KineticParams;
pub use reactivity::Reactivity;
pub use solve::{prompt_jump, solve, Method, Solution, SolverOptions, State, TimeGrid};

#[cfg(test)]
mod fixture_tests;
