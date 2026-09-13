//! Dense real weighted least-squares over caller-supplied matrices.
//!
//! Shared kernel: solves `min_β Σ w_i (X_i β − y_i)²` for a
//! caller design matrix `X` (row-major `Vec<Vec<f64>>`), targets `y`, and
//! non-negative weights `w`. Rows scale by `sqrt(w)` and the scaled system
//! solves through the workspace backend (faer 0.20 column-pivoted QR
//! [`solve_lstsq`](faer::linalg::solvers::SpSolverLstsq::solve_lstsq)); this
//! crate stays the only workspace member that names `faer`, so no dependency
//! change was needed to land this module.
//!
//! Consumers: the UQ-lite sampler in [`crate::sample`] keeps its own
//! Cholesky/eigen factor path for covariance blocks, while the
//! spectroscopy E7 efficiency-coefficient fit builds its log-space design
//! matrix and calls [`weighted_lstsq`] — one dense-real module, two
//! consumers, no second least-squares route.
//!
//! No evaluated data lives here: every entry of `X`/`y`/`w` is a caller
//! input, validated for shape and finiteness before the backend runs.

use faer::linalg::solvers::{ColPivQr, SpSolverLstsq};
use faer::Mat;

/// Errors surfaced by [`weighted_lstsq`]. Every rejection names its cause.
#[derive(Debug, Clone, PartialEq)]
pub enum LstsqError {
    /// No rows supplied, or the leading row is empty.
    Empty,
    /// `y` (or `w`) length does not match the `X` row count.
    DimensionMismatch {
        /// Row count of `X`.
        expected: usize,
        /// Length actually supplied.
        got: usize,
    },
    /// Row `row` length differs from the leading row width.
    RowLength {
        /// Offending row index.
        row: usize,
        /// Leading row width.
        expected: usize,
        /// Offending row length.
        got: usize,
    },
    /// A non-finite value in the named input (`"x"`, `"y"`, or `"weights"`).
    NonFinite(&'static str),
    /// A negative weight at the given row index.
    NegativeWeight {
        /// Offending row index.
        index: usize,
    },
    /// Every weight is zero: no observation constrains the fit.
    NoPositiveWeight,
    /// Fewer rows than columns (`rows` observations, `cols` unknowns).
    Underdetermined {
        /// Observation count.
        rows: usize,
        /// Unknown count.
        cols: usize,
    },
    /// The backend failed or returned a non-finite solution (for example a
    /// rank-deficient design).
    Backend(String),
}

impl std::fmt::Display for LstsqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LstsqError::Empty => write!(f, "least-squares input is empty"),
            LstsqError::DimensionMismatch { expected, got } => write!(
                f,
                "least-squares dimension mismatch: expected {expected}, got {got}"
            ),
            LstsqError::RowLength { row, expected, got } => write!(
                f,
                "least-squares row {row} has length {got}, expected {expected}"
            ),
            LstsqError::NonFinite(what) => {
                write!(f, "least-squares input `{what}` holds a non-finite value")
            }
            LstsqError::NegativeWeight { index } => {
                write!(f, "least-squares weight at row {index} is negative")
            }
            LstsqError::NoPositiveWeight => write!(
                f,
                "least-squares weights are all zero; no observation constrains the fit"
            ),
            LstsqError::Underdetermined { rows, cols } => write!(
                f,
                "least-squares underdetermined: {rows} rows, {cols} columns"
            ),
            LstsqError::Backend(m) => write!(f, "linalg backend error: {m}"),
        }
    }
}

impl std::error::Error for LstsqError {}

/// Weighted least-squares solve: `min_β Σ w_i (X_i β − y_i)²`.
///
/// `x` holds one row per observation (all rows the same width), `y` and `w`
/// hold one entry per row. Weights must be finite and non-negative with at
/// least one positive entry; a zero weight skips its row. Needs at least as
/// many rows as columns. Returns one coefficient per column.
pub fn weighted_lstsq(x: &[Vec<f64>], y: &[f64], w: &[f64]) -> Result<Vec<f64>, LstsqError> {
    let m = x.len();
    if m == 0 {
        return Err(LstsqError::Empty);
    }
    let n = x[0].len();
    if n == 0 {
        return Err(LstsqError::Empty);
    }
    for (i, row) in x.iter().enumerate() {
        if row.len() != n {
            return Err(LstsqError::RowLength {
                row: i,
                expected: n,
                got: row.len(),
            });
        }
    }
    if y.len() != m {
        return Err(LstsqError::DimensionMismatch {
            expected: m,
            got: y.len(),
        });
    }
    if w.len() != m {
        return Err(LstsqError::DimensionMismatch {
            expected: m,
            got: w.len(),
        });
    }
    if x.iter().flatten().any(|v| !v.is_finite()) {
        return Err(LstsqError::NonFinite("x"));
    }
    if y.iter().any(|v| !v.is_finite()) {
        return Err(LstsqError::NonFinite("y"));
    }
    if w.iter().any(|v| !v.is_finite()) {
        return Err(LstsqError::NonFinite("weights"));
    }
    for (i, wi) in w.iter().enumerate() {
        if *wi < 0.0 {
            return Err(LstsqError::NegativeWeight { index: i });
        }
    }
    if !w.iter().any(|wi| *wi > 0.0) {
        return Err(LstsqError::NoPositiveWeight);
    }
    if m < n {
        return Err(LstsqError::Underdetermined { rows: m, cols: n });
    }

    // Scale rows by sqrt(w), then solve the scaled least-squares system with
    // the column-pivoted QR route (already enabled via the `linalg` feature).
    let a = Mat::<f64>::from_fn(m, n, |i, j| x[i][j] * w[i].sqrt());
    let b = Mat::<f64>::from_fn(m, 1, |i, _| y[i] * w[i].sqrt());
    let qr = ColPivQr::<f64>::new(a.as_ref());
    let sol = qr.solve_lstsq(b.as_ref());
    let mut beta = Vec::with_capacity(n);
    for j in 0..n {
        beta.push(sol[(j, 0)]);
    }
    if beta.iter().any(|v| !v.is_finite()) {
        return Err(LstsqError::Backend(
            "least-squares solution is non-finite (rank-deficient design?)".to_string(),
        ));
    }
    Ok(beta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_line_recovers_closed_form() {
        // y = 1 + 2t on t = 0, 1, 2 (exact, full rank).
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]];
        let y = vec![1.0, 3.0, 5.0];
        let w = vec![1.0, 1.0, 1.0];
        let beta = weighted_lstsq(&x, &y, &w).unwrap();
        assert!((beta[0] - 1.0).abs() < 1e-12, "b0 = {}", beta[0]);
        assert!((beta[1] - 2.0).abs() < 1e-12, "b1 = {}", beta[1]);
    }

    #[test]
    fn nonuniform_weights_keep_exact_line() {
        // Same collinear points: any positive weights recover the line.
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]];
        let y = vec![1.0, 3.0, 5.0];
        let w = vec![1.0, 4.0, 9.0];
        let beta = weighted_lstsq(&x, &y, &w).unwrap();
        assert!((beta[0] - 1.0).abs() < 1e-12, "b0 = {}", beta[0]);
        assert!((beta[1] - 2.0).abs() < 1e-12, "b1 = {}", beta[1]);
    }

    #[test]
    fn zero_weight_row_is_ignored() {
        // (0,0) and (1,1) pin y = x; the (2,100) outlier carries zero weight.
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]];
        let y = vec![0.0, 1.0, 100.0];
        let w = vec![1.0, 1.0, 0.0];
        let beta = weighted_lstsq(&x, &y, &w).unwrap();
        assert!(beta[0].abs() < 1e-12, "b0 = {}", beta[0]);
        assert!((beta[1] - 1.0).abs() < 1e-12, "b1 = {}", beta[1]);
    }

    #[test]
    fn quadratic_recovers_closed_form() {
        // y = 0.5 - t + 0.25 t^2 on t = 0, 1, 2, 3 (exact, overdetermined).
        let ts = [0.0, 1.0, 2.0, 3.0];
        let x: Vec<Vec<f64>> = ts.iter().map(|t| vec![1.0, *t, t * t]).collect();
        let y: Vec<f64> = ts.iter().map(|t| 0.5 - t + 0.25 * t * t).collect();
        let w = vec![1.0; 4];
        let beta = weighted_lstsq(&x, &y, &w).unwrap();
        assert!((beta[0] - 0.5).abs() < 1e-12, "b0 = {}", beta[0]);
        assert!((beta[1] + 1.0).abs() < 1e-12, "b1 = {}", beta[1]);
        assert!((beta[2] - 0.25).abs() < 1e-12, "b2 = {}", beta[2]);
    }

    #[test]
    fn malformed_inputs_name_their_cause() {
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0]];
        let y = vec![1.0, 2.0];
        let w = vec![1.0, 1.0];
        assert_eq!(weighted_lstsq(&[], &y, &w).unwrap_err(), LstsqError::Empty);
        assert_eq!(
            weighted_lstsq(&[vec![]], &[1.0], &[1.0]).unwrap_err(),
            LstsqError::Empty
        );
        assert_eq!(
            weighted_lstsq(&x, &[1.0], &w).unwrap_err(),
            LstsqError::DimensionMismatch {
                expected: 2,
                got: 1
            }
        );
        assert_eq!(
            weighted_lstsq(&x, &y, &[1.0]).unwrap_err(),
            LstsqError::DimensionMismatch {
                expected: 2,
                got: 1
            }
        );
        assert_eq!(
            weighted_lstsq(&[vec![1.0], vec![1.0, 2.0]], &y, &w).unwrap_err(),
            LstsqError::RowLength {
                row: 1,
                expected: 1,
                got: 2
            }
        );
        assert_eq!(
            weighted_lstsq(&x, &[f64::NAN, 1.0], &w).unwrap_err(),
            LstsqError::NonFinite("y")
        );
        assert_eq!(
            weighted_lstsq(&[vec![f64::INFINITY, 0.0], vec![1.0, 1.0]], &y, &w).unwrap_err(),
            LstsqError::NonFinite("x")
        );
        assert_eq!(
            weighted_lstsq(&x, &y, &[1.0, f64::NAN]).unwrap_err(),
            LstsqError::NonFinite("weights")
        );
        assert_eq!(
            weighted_lstsq(&x, &y, &[1.0, -1.0]).unwrap_err(),
            LstsqError::NegativeWeight { index: 1 }
        );
        assert_eq!(
            weighted_lstsq(&x, &y, &[0.0, 0.0]).unwrap_err(),
            LstsqError::NoPositiveWeight
        );
        assert_eq!(
            weighted_lstsq(&[vec![1.0, 2.0, 3.0]], &[1.0], &[1.0]).unwrap_err(),
            LstsqError::Underdetermined { rows: 1, cols: 3 }
        );
    }

    #[test]
    fn error_display_strings() {
        assert_eq!(
            LstsqError::Empty.to_string(),
            "least-squares input is empty"
        );
        assert_eq!(
            LstsqError::DimensionMismatch {
                expected: 2,
                got: 1
            }
            .to_string(),
            "least-squares dimension mismatch: expected 2, got 1"
        );
        assert_eq!(
            LstsqError::RowLength {
                row: 1,
                expected: 1,
                got: 2
            }
            .to_string(),
            "least-squares row 1 has length 2, expected 1"
        );
        assert_eq!(
            LstsqError::NonFinite("y").to_string(),
            "least-squares input `y` holds a non-finite value"
        );
        assert_eq!(
            LstsqError::NegativeWeight { index: 1 }.to_string(),
            "least-squares weight at row 1 is negative"
        );
        assert_eq!(
            LstsqError::NoPositiveWeight.to_string(),
            "least-squares weights are all zero; no observation constrains the fit"
        );
        assert_eq!(
            LstsqError::Underdetermined { rows: 1, cols: 3 }.to_string(),
            "least-squares underdetermined: 1 rows, 3 columns"
        );
        assert_eq!(
            LstsqError::Backend("boom".to_string()).to_string(),
            "linalg backend error: boom"
        );
    }

    #[test]
    fn error_implements_std_error_trait() {
        let err: &dyn std::error::Error = &LstsqError::Empty;
        assert!(err.source().is_none());
    }
}
