//! Inhour equation (E3) and the delayed-supercritical stable period.
//!
//! Assuming exponential solutions `n, C_i ~ e^{ωt}` at constant reactivity
//! in (E1) gives `C_i = β_i n / (Λ(ω + λ_i))`; substituting into (E1a) and
//! dividing by `n` yields the inhour relation
//!
//! ```text
//! rho(ω) = Λ ω + Σ_i β_i ω / (ω + λ_i).                        (E3)
//! ```
//!
//! This is a textbook algebraic consequence of (E1) — no new physics is
//! introduced, and [`rho_of_omega`] is verified by direct substitution (the
//! residual gate below). For `0 < ρ < β` the right-hand side is strictly
//! increasing on `ω > 0` (derivative `Λ + Σ β_i λ_i/(ω+λ_i)² > 0`), so there
//! is exactly one positive root: the stable (asymptotic) period
//! `T = 1/ω`.

use crate::error::Error;
use crate::params::KineticParams;

/// Evaluate the inhour right-hand side `ρ(ω)` (E3) \[Δk\].
///
/// Errors on non-finite `ω` or `ω` at/below the slowest-decay pole
/// (`ω <= -min λ_i`), where a denominator vanishes.
pub fn rho_of_omega(params: &KineticParams, omega: f64) -> Result<f64, Error> {
    if !omega.is_finite() {
        return Err(Error::BadData("omega must be finite"));
    }
    let lambda_min = params
        .lambdas()
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    if omega <= -lambda_min {
        return Err(Error::InhourPole { omega, lambda_min });
    }
    let mut rho = params.lambda_gen() * omega;
    for (b, l) in params.betas().iter().zip(params.lambdas()) {
        rho += b * omega / (omega + l);
    }
    Ok(rho)
}

/// Inhour residual `ρ − ρ(ω)` \[Δk\].
///
/// The oracle gate asserts `|residual| < 1e-10` at a known root: for the
/// [`stable_period`] root this is a self-consistency check of the
/// bisection, and for hand-built `(ρ, ω)` fixture pairs it pins (E3)
/// against values computed outside this crate.
pub fn residual(params: &KineticParams, rho: f64, omega: f64) -> Result<f64, Error> {
    Ok(rho - rho_of_omega(params, omega)?)
}

/// Asymptotic (stable) period `T = 1/ω` \[s\] for `0 < ρ < β`.
///
/// Finds the unique positive inhour root by bisection: brackets `[tiny, hi]`
/// with `hi` doubled from `min λ_i` until `ρ(hi) > ρ` (monotonicity on
/// `ω > 0` guarantees a sign change), then bisects to a relative width of
/// `1e-13`. Errors when `ρ` is outside `(0, β)`.
pub fn stable_period(params: &KineticParams, rho: f64) -> Result<f64, Error> {
    let beta = params.beta_total();
    if !rho.is_finite() || rho <= 0.0 || rho >= beta {
        return Err(Error::NoStablePeriod {
            rho,
            beta_total: beta,
        });
    }
    let lambda_min = params
        .lambdas()
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let f = |w: f64| rho_of_omega(params, w).map(|r| r - rho).unwrap_or(f64::NAN);
    let mut lo = 0.0_f64;
    let mut flo = -rho;
    let mut hi = lambda_min;
    let mut fhi = f(hi);
    for _ in 0..1024 {
        if fhi.is_nan() {
            return Err(Error::InhourPole {
                omega: hi,
                lambda_min,
            });
        }
        if fhi > 0.0 {
            break;
        }
        hi *= 2.0;
        if !hi.is_finite() {
            return Err(Error::NoStablePeriod {
                rho,
                beta_total: beta,
            });
        }
        fhi = f(hi);
    }
    if fhi.is_nan() || fhi <= 0.0 {
        return Err(Error::NoStablePeriod {
            rho,
            beta_total: beta,
        });
    }
    for _ in 0..1000 {
        let mid = 0.5 * (lo + hi);
        if hi - lo <= 1e-13 * hi.max(1e-300) {
            break;
        }
        let fm = f(mid);
        if fm > 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    let _ = flo;
    let omega = 0.5 * (lo + hi);
    Ok(1.0 / omega)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::KineticParams;

    fn one_group() -> KineticParams {
        KineticParams::new(vec![0.0065], vec![0.08], 1e-4).unwrap()
    }

    #[test]
    fn rho_of_zero_is_zero() {
        assert_eq!(rho_of_omega(&one_group(), 0.0).unwrap(), 0.0);
    }

    #[test]
    fn residual_at_bisected_root_below_1e_10() {
        // Self-consistency gate: the bisection root must satisfy (E3).
        let p = crate::params::tests::six_group();
        let rho = 0.002;
        let t = stable_period(&p, rho).unwrap();
        let res = residual(&p, rho, 1.0 / t).unwrap().abs();
        assert!(res < 1e-10, "residual {res}");
    }

    #[test]
    fn one_group_period_matches_quadratic() {
        // Independent check: 1-group inhour  ρ = Λω + βω/(ω+λ) rearranges to
        // Λω² + (Λλ + β − ρ)ω − ρλ = 0; take the positive root directly.
        let p = one_group();
        let (beta, lam, gen) = (0.0065_f64, 0.08_f64, 1e-4_f64);
        for rho in [0.0005, 0.002, 0.005] {
            let t = stable_period(&p, rho).unwrap();
            let (a, b, c) = (gen, gen * lam + beta - rho, -rho * lam);
            let w = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
            assert!((t - 1.0 / w).abs() / t < 1e-9, "rho={rho} t={t}");
        }
    }

    #[test]
    fn rejects_out_of_range() {
        let p = one_group();
        assert!(stable_period(&p, 0.0).is_err());
        assert!(stable_period(&p, -0.001).is_err());
        assert!(stable_period(&p, 0.0065).is_err());
        assert!(stable_period(&p, 0.01).is_err());
        assert!(rho_of_omega(&p, f64::NAN).is_err());
        assert!(rho_of_omega(&p, -0.08).is_err());
    }
}
