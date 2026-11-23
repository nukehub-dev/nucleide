//! Covariance whitening by diagonal equilibration and a Cholesky solve.

use faer::{Mat, Side};

/// Lower bound on the accepted ratio of minimum to maximum correlation eigenvalue.
pub const MIN_RELATIVE_EIGENVALUE: f64 = 1e-12;

const SYMMETRY_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// A covariance or observation block that cannot be whitened reliably.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CovarianceError {
    /// The covariance is empty, nonsquare, or does not match the observation rows.
    BadShape,
    /// Observation rows must have a common nonzero width.
    BadObservations,
    /// An entry is not finite.
    NonFinite,
    /// The covariance has a nonpositive variance.
    NonPositiveVariance,
    /// The matrix is asymmetric beyond the tolerance on its correlation scale.
    Asymmetric,
    /// The covariance is not positive definite.
    NotPositiveDefinite,
    /// The correlation matrix is too ill-conditioned for whitening.
    IllConditioned,
    /// A scaled input or whitened output is not finite.
    Numerical,
}

impl std::fmt::Display for CovarianceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let detail = match self {
            Self::BadShape => "expected a square matrix matching the observation count",
            Self::BadObservations => "expected equal-length nonempty observation rows",
            Self::NonFinite => "non-finite covariance or observation entry",
            Self::NonPositiveVariance => "diagonal variances must be positive",
            Self::Asymmetric => "matrix must be symmetric on the correlation scale",
            Self::NotPositiveDefinite => "matrix must be positive definite",
            Self::IllConditioned => "correlation matrix is too ill-conditioned to whiten",
            Self::Numerical => "whitening arithmetic produced a non-finite value",
        };
        write!(f, "covariance: {detail}")
    }
}

impl std::error::Error for CovarianceError {}

/// Result alias for covariance whitening.
pub type Result<T> = std::result::Result<T, CovarianceError>;

/// Whiten columns of an observation matrix with a shared covariance factor.
///
/// For covariance `C = D K D`, uses standard deviations on the diagonal of
/// `D`, and factors the correlation matrix as `K = L Lᵀ`. Returns
/// `L⁻¹ D⁻¹ observations` using triangular solves. Correlation pairs may
/// differ by at most `64 * f64::EPSILON` and are averaged within that
/// tolerance. Both positive definiteness and the conditioning of `K` are
/// checked, so changing observation units does not change the condition test.
pub fn whiten(covariance: &[Vec<f64>], observations: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let rows = observations.len();
    if rows == 0 || covariance.len() != rows || covariance.iter().any(|row| row.len() != rows) {
        return Err(CovarianceError::BadShape);
    }
    let columns = observations[0].len();
    if columns == 0 || observations.iter().any(|row| row.len() != columns) {
        return Err(CovarianceError::BadObservations);
    }
    if covariance
        .iter()
        .chain(observations)
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(CovarianceError::NonFinite);
    }
    let mut sigmas = Vec::with_capacity(rows);
    for (i, row) in covariance.iter().enumerate() {
        if row[i] <= 0.0 {
            return Err(CovarianceError::NonPositiveVariance);
        }
        sigmas.push(row[i].sqrt());
    }
    let mut correlation = Mat::<f64>::identity(rows, rows);
    for i in 0..rows {
        for j in 0..i {
            // Divide by the smaller sigma first to retain correlations even
            // when the observation units differ by many orders of magnitude.
            let small = sigmas[i].min(sigmas[j]);
            let large = sigmas[i].max(sigmas[j]);
            let lower = covariance[i][j] / small / large;
            let upper = covariance[j][i] / small / large;
            if lower.abs() > 1.0 || upper.abs() > 1.0 {
                return Err(CovarianceError::NotPositiveDefinite);
            }
            if (lower - upper).abs() > SYMMETRY_TOLERANCE {
                return Err(CovarianceError::Asymmetric);
            }
            correlation[(i, j)] = lower + (upper - lower) * 0.5;
        }
    }
    let factor = correlation
        .as_ref()
        .cholesky(Side::Lower)
        .map_err(|_| CovarianceError::NotPositiveDefinite)?;
    let eigenvalues = correlation.as_ref().selfadjoint_eigenvalues(Side::Lower);
    if eigenvalues.iter().any(|v| !v.is_finite()) {
        return Err(CovarianceError::Numerical);
    }
    let minimum = eigenvalues.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum = eigenvalues.iter().copied().fold(0.0f64, f64::max);
    if minimum <= MIN_RELATIVE_EIGENVALUE * maximum {
        return Err(CovarianceError::IllConditioned);
    }

    let mut whitened = Mat::from_fn(rows, columns, |i, j| observations[i][j] / sigmas[i]);
    if (0..rows).any(|i| (0..columns).any(|j| !whitened[(i, j)].is_finite())) {
        return Err(CovarianceError::Numerical);
    }
    factor
        .compute_l()
        .as_ref()
        .solve_lower_triangular_in_place(whitened.as_mut());
    (0..rows)
        .map(|i| {
            (0..columns)
                .map(|j| {
                    let value = whitened[(i, j)];
                    value
                        .is_finite()
                        .then_some(value)
                        .ok_or(CovarianceError::Numerical)
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn whitening_preserves_geometry_under_unequal_observation_units() {
        let base = whiten(&[vec![4.0, 3.0], vec![3.0, 9.0]], &[vec![2.0], vec![6.0]]).unwrap();
        let scaled = whiten(
            &[vec![4e-200, 3.0], vec![3.0, 9e200]],
            &[vec![2e-100], vec![6e100]],
        )
        .unwrap();
        assert_eq!(base[0][0], 1.0);
        assert!((base[1][0] - 3.0f64.sqrt()).abs() < 1e-14);
        for (a, b) in base.iter().zip(scaled) {
            assert!((a[0] - b[0]).abs() < 1e-14);
        }
    }

    #[test]
    fn conditioning_rejects_a_nearly_singular_positive_definite_covariance() {
        let rho = 1.0 - 1e-13;
        assert_eq!(
            whiten(&[vec![1.0, rho], vec![rho, 1.0]], &[vec![1.0], vec![1.0]]),
            Err(CovarianceError::IllConditioned)
        );
    }
}
