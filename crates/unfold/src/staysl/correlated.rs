//! Correlated weighting for the anchored STAYSL-class adjustment.

use crate::error::{Error, Result};
use nucleide_linalg::{covariance, lstsq};

/// Default relative first-order optimality tolerance for a correlated fit.
pub const DEFAULT_CORRELATED_TOLERANCE: f64 = 1e-10;

/// Spectrum and diagnostics for a nonnegative correlated adjustment.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelatedSolution {
    /// Adjusted spectrum in original group order.
    pub spectrum: Vec<f64>,
    /// Rates folded through the original response matrix, in detector order.
    pub rates: Vec<f64>,
    /// Squared norm of the whitened measurement residual.
    pub measurement_chi2: f64,
    /// Damping coefficient times squared distance from the guess.
    pub damping_penalty: f64,
    /// Number of active-set iterations performed.
    pub iterations: usize,
    /// Requested relative first-order optimality tolerance.
    pub tolerance: f64,
    /// Largest scaled violation of the bound-constrained gradient conditions.
    pub optimality: f64,
}

/// Adjust a spectrum with caller-supplied measurement covariance.
///
/// For rates `y`, response `R`, guess `x₀`, and covariance `C`, minimize
/// `(Rx-y)ᵀ C⁻¹ (Rx-y) + damping * ||x-x₀||²` subject to `x >= 0`.
/// The covariance rows and columns match the detector order and their entries
/// carry products of rate units. The rates and responses must be nonnegative
/// and finite, the guess strictly positive and finite, and damping positive
/// and finite. Zero measured rates and zero response rows participate in the
/// residual, and groups with no response keep the guess exactly.
///
/// The solver whitens the response and measurements, then uses a feasible
/// active-set solve on the damped design. Bound variables can re-enter the
/// free set. Each unconstrained subproblem uses the shared QR least-squares
/// kernel. A successful fit satisfies the relative first-order gradient
/// check for the free and bound groups. `tolerance` must be in `(0, 1)`.
/// Exhausting the positive iteration cap returns [`Error::NotConverged`].
pub fn unfold_correlated(
    response: &[Vec<f64>],
    rates: &[f64],
    covariance: &[Vec<f64>],
    guess: &[f64],
    tolerance: f64,
    max_iterations: usize,
    damping: f64,
) -> Result<CorrelatedSolution> {
    let groups = crate::validate_response(response)?;
    crate::validate_spectrum("guess", guess, groups, true, Error::BadGuess)?;
    if rates.len() != response.len() {
        return Err(Error::BadShape {
            what: "rates",
            expected: response.len(),
            got: rates.len(),
        });
    }
    if rates.iter().any(|rate| !rate.is_finite() || *rate < 0.0) {
        return Err(Error::BadRates("rates must be finite and nonnegative"));
    }
    if !(0.0 < tolerance && tolerance < 1.0) {
        return Err(Error::BadOption("tolerance must be in (0, 1)"));
    }
    if max_iterations == 0 {
        return Err(Error::BadOption("max_iterations must be at least 1"));
    }
    if !damping.is_finite() || damping <= 0.0 {
        return Err(Error::BadOption("damping must be finite and positive"));
    }

    let observations: Vec<Vec<f64>> = response
        .iter()
        .zip(rates)
        .map(|(row, &rate)| {
            let mut augmented = row.clone();
            augmented.push(rate);
            augmented
        })
        .collect();
    let mut design = covariance::whiten(covariance, &observations)?;
    let mut target = Vec::with_capacity(response.len() + groups);
    for row in &mut design {
        target.push(row[groups]);
        row.truncate(groups);
    }
    let sqrt_damping = damping.sqrt();
    for (j, prior) in guess.iter().enumerate() {
        let mut row = vec![0.0; groups];
        row[j] = sqrt_damping;
        design.push(row);
        target.push(sqrt_damping * prior);
    }
    let scales: Vec<f64> = (0..groups)
        .map(|j| design.iter().fold(0.0f64, |norm, row| norm.hypot(row[j])))
        .collect();
    let initial: Vec<f64> = guess
        .iter()
        .zip(&scales)
        .map(|(x, scale)| x * scale)
        .collect();
    finite(scales.iter().chain(&target).chain(&initial).copied())?;
    for row in &mut design {
        for (value, scale) in row.iter_mut().zip(&scales) {
            *value /= scale;
        }
    }

    let (mut scaled_spectrum, iterations) =
        nonnegative_fit(&design, &target, initial, tolerance, max_iterations)?;
    let mut spectrum = Vec::with_capacity(groups);
    for (j, &scale) in scales.iter().enumerate() {
        let flux = if response.iter().all(|row| row[j] == 0.0) {
            guess[j]
        } else {
            scaled_spectrum[j] / scale
        };
        if scaled_spectrum[j] > 0.0 && flux == 0.0 {
            return Err(Error::CorrelatedNumerical("positive spectrum underflowed"));
        }
        spectrum.push(flux);
        scaled_spectrum[j] = flux * scale;
    }
    finite(spectrum.iter().chain(&scaled_spectrum).copied())?;
    let folded_rates = crate::fold(response, &spectrum);
    let fitted = crate::fold(&design, &scaled_spectrum);
    let (optimality, _) = check_gradient(&design, &target, &fitted, &scaled_spectrum)?;
    if optimality > tolerance {
        return Err(Error::CorrelatedNumerical(
            "gradient check failed after scaling",
        ));
    }
    let (measurement, prior) = fitted.split_at(response.len());
    let (observed_rates, observed_prior) = target.split_at(response.len());
    let measurement_chi2 = squared_error(measurement, observed_rates);
    let damping_penalty = squared_error(prior, observed_prior);
    finite(
        folded_rates
            .iter()
            .copied()
            .chain([measurement_chi2, damping_penalty]),
    )?;
    Ok(CorrelatedSolution {
        spectrum,
        rates: folded_rates,
        measurement_chi2,
        damping_penalty,
        iterations,
        tolerance,
        optimality,
    })
}

fn finite(values: impl IntoIterator<Item = f64>) -> Result<()> {
    if values.into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(Error::CorrelatedNumerical("non-finite fit arithmetic"))
    }
}

fn squared_error(fitted: &[f64], target: &[f64]) -> f64 {
    fitted
        .iter()
        .zip(target)
        .map(|(&a, &b)| (a - b).powi(2))
        .sum()
}

/// Solve the current free subproblem, then move along a feasible segment
/// to its optimum or the first blocking bound. At a feasible subproblem
/// optimum, release the bound with the largest negative scaled gradient.
fn nonnegative_fit(
    design: &[Vec<f64>],
    target: &[f64],
    mut spectrum: Vec<f64>,
    tolerance: f64,
    max_iterations: usize,
) -> Result<(Vec<f64>, usize)> {
    let mut is_free = vec![true; spectrum.len()];
    let weights = vec![1.0; target.len()];
    for iteration in 1..=max_iterations {
        let free: Vec<usize> = (0..spectrum.len()).filter(|&j| is_free[j]).collect();
        if !free.is_empty() {
            let reduced: Vec<Vec<f64>> = design
                .iter()
                .map(|row| free.iter().map(|&j| row[j]).collect())
                .collect();
            let candidate = lstsq::weighted_lstsq(&reduced, target, &weights)
                .map_err(|err| Error::Lstsq(err.to_string()))?;
            if move_to_candidate(&free, &candidate, &mut spectrum, &mut is_free) {
                continue;
            }
        }
        let fitted = crate::fold(design, &spectrum);
        let (optimality, release) = check_gradient(design, target, &fitted, &spectrum)?;
        if optimality <= tolerance {
            return Ok((spectrum, iteration));
        }
        let Some(index) = release else {
            return Err(Error::CorrelatedNumerical(
                "free-set gradient did not meet tolerance",
            ));
        };
        is_free[index] = true;
    }
    Err(Error::NotConverged)
}

fn move_to_candidate(
    free: &[usize],
    candidate: &[f64],
    spectrum: &mut [f64],
    is_free: &mut [bool],
) -> bool {
    let boundary = |current: f64, trial: f64| {
        if current == 0.0 {
            return 0.0;
        }
        let scale = current.max(-trial);
        (current / scale) / (current / scale - trial / scale)
    };
    let step = free
        .iter()
        .zip(candidate)
        .filter(|(_, trial)| **trial <= 0.0)
        .map(|(&j, &trial)| boundary(spectrum[j], trial))
        .reduce(f64::min);
    let Some(step) = step else {
        for (&j, &trial) in free.iter().zip(candidate) {
            spectrum[j] = trial;
        }
        return false;
    };
    for (&j, &trial) in free.iter().zip(candidate) {
        let current = spectrum[j];
        let value = (1.0 - step) * current + step * trial;
        if value <= 0.0 || (trial <= 0.0 && boundary(current, trial) == step) {
            spectrum[j] = 0.0;
            is_free[j] = false;
        } else {
            spectrum[j] = value;
        }
    }
    true
}

// Compute the half-objective gradient. Positive groups require zero
// gradient and zero groups require nonnegative gradient. Scale each
// component by the sum of magnitudes of its terms before cancellation.
fn check_gradient(
    design: &[Vec<f64>],
    target: &[f64],
    fitted: &[f64],
    spectrum: &[f64],
) -> Result<(f64, Option<usize>)> {
    finite(fitted.iter().copied())?;
    let mut worst = 0.0f64;
    let mut worst_bound_violation = 0.0;
    let mut release = None;
    for (j, &value) in spectrum.iter().enumerate() {
        let (gradient, scale) = design.iter().zip(target).zip(fitted).fold(
            (0.0, 0.0),
            |(gradient, scale), ((row, &b), &prediction)| {
                let a = row[j];
                (
                    gradient + a * (prediction - b),
                    scale + a.abs() * prediction.abs() + a.abs() * b.abs(),
                )
            },
        );
        finite([gradient, scale])?;
        let relative = if scale == 0.0 { 0.0 } else { gradient / scale };
        let violation = if value > 0.0 {
            relative.abs()
        } else {
            (-relative).max(0.0)
        };
        worst = worst.max(violation);
        if value == 0.0 && violation > worst_bound_violation {
            worst_bound_violation = violation;
            release = Some(j);
        }
    }
    Ok((worst, release))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = DEFAULT_CORRELATED_TOLERANCE;

    #[test]
    fn correlated_two_detector_fit_matches_closed_form() {
        let fit = unfold_correlated(
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
            &[3.0, 2.0],
            &[vec![1.0, 0.5], vec![0.5, 1.0]],
            &[1.0, 1.0],
            TOL,
            1,
            1.0,
        )
        .unwrap();
        assert!((fit.spectrum[0] - 29.0 / 15.0).abs() < 1e-12);
        assert!((fit.spectrum[1] - 19.0 / 15.0).abs() < 1e-12);
        assert_eq!(fit.rates, fit.spectrum);
        assert!(fit.optimality <= fit.tolerance);
    }

    #[test]
    fn correlated_bound_solution_respects_iteration_cap() {
        let response = [
            vec![0.0, 3.0, 3.0],
            vec![0.0, 1.0, 0.0],
            vec![3.0, 2.0, 0.0],
        ];
        let covariance = [
            vec![1.0, 0.75, 0.0],
            vec![0.75, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let solve = |cap| {
            unfold_correlated(
                &response,
                &[6.0, 4.0, 1.0],
                &covariance,
                &[1.0; 3],
                TOL,
                cap,
                1.0,
            )
        };
        assert_eq!(solve(1), Err(Error::NotConverged));
        let fit = solve(20).unwrap();
        assert_eq!(fit.spectrum[0], 0.0);
        assert!((fit.spectrum[1] - 151.0 / 141.0).abs() < 1e-12);
        assert!((fit.spectrum[2] - 11.0 / 47.0).abs() < 1e-12);
        assert!(fit.optimality <= TOL);
    }

    #[test]
    fn a_group_can_return_after_hitting_zero() {
        // Group 3 hits zero on the first step and reenters after group 1
        // hits zero. Freezing either group permanently breaks the fit.
        let response = [
            vec![4.0, 4.0, 2.0, 2.0],
            vec![2.0, 0.0, 1.0, 2.0],
            vec![3.0, 1.0, 3.0, 3.0],
        ];
        let fit = unfold_correlated(
            &response,
            &[2.0, 1.0, 3.0],
            &[
                vec![1.0, 0.75, 0.0],
                vec![0.75, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            &[3.0, 2.0, 2.0, 1.0],
            TOL,
            20,
            1.0,
        )
        .unwrap();
        assert_eq!(fit.iterations, 4);
        for (actual, expected) in
            fit.spectrum
                .iter()
                .zip([235.0 / 764.0, 0.0, 519.0 / 764.0, 17.0 / 2292.0])
        {
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn zero_response_and_zero_rates_are_valid_in_the_damped_objective() {
        let fit = unfold_correlated(
            &[vec![1.0, 0.0], vec![0.0, 0.0]],
            &[0.0, 2.0],
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
            &[1.0, 7.2],
            TOL,
            20,
            1.0,
        )
        .unwrap();
        assert!((fit.spectrum[0] - 0.5).abs() < 1e-12);
        assert_eq!(fit.spectrum[1], 7.2);
        assert_eq!(fit.rates[1], 0.0);
        assert!((fit.measurement_chi2 - 4.25).abs() < 1e-12);
        assert!((fit.damping_penalty - 0.25).abs() < 1e-12);
    }
}
