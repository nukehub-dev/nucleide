//! Numeric backend isolation layer.
//!
//! The workspace depends on this crate — never on `faer` (or any other LA
//! backend) directly. This keeps a single swap point if the backend API
//! churns pre-1.0.
//!
//! Sparse complex LU with symbolic/numeric split: analyze the sparsity
//! pattern once per chain topology, then reuse it across CRAM theta values
//! and depletion time steps.
//!
//! Current backend: faer 0.20 (pure Rust, native complex support, rayon).

use std::sync::Arc;

/// Backend selected for the implementation (informational).
pub const BACKEND: &str = "faer-0.20";

/// Complex scalar used throughout the facade.
pub use faer::complex_native::c64 as C64;

/// Zero constant for the scalar type.
pub const C64_ZERO: C64 = C64 { re: 0.0, im: 0.0 };

/// Errors surfaced by the numeric backend.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Matrix construction or factorization failed inside the backend.
    Backend(String),
    /// Value count does not match the pattern's nonzeros.
    Shape { expected: usize, got: usize },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Backend(m) => write!(f, "linalg backend error: {m}"),
            Error::Shape { expected, got } => {
                write!(f, "shape mismatch: expected {expected} values, got {got}")
            }
        }
    }
}
impl std::error::Error for Error {}

/// Shared sparsity pattern of a square CSC matrix.
///
/// Build once per chain topology; clones are cheap (Arc).
#[derive(Clone)]
pub struct Pattern {
    symbolic: Arc<faer::sparse::SymbolicSparseColMat<usize>>,
    order: Arc<faer::sparse::ValuesOrder<usize>>,
    pub n: usize,
    pub nnz: usize,
}

impl Pattern {
    /// Build a pattern from `(row, col)` pairs. Duplicate coordinates are
    /// collapsed by position (values must be supplied accordingly — the
    /// depletion matrix builder pre-sums duplicates before calling this).
    pub fn from_entries(n: usize, entries: &[(usize, usize)]) -> Result<Self, Error> {
        let (symbolic, order) =
            faer::sparse::SymbolicSparseColMat::<usize>::try_new_from_indices(n, n, entries)
                .map_err(|e| Error::Backend(e.to_string()))?;
        let nnz = symbolic.compute_nnz();
        Ok(Self {
            symbolic: Arc::new(symbolic),
            order: Arc::new(order),
            n,
            nnz,
        })
    }

    pub fn nrows(&self) -> usize {
        self.n
    }
}

/// A complex CSC matrix bound to a shared [`Pattern`].
#[derive(Clone)]
pub struct ComplexCsc {
    inner: Arc<faer::sparse::SparseColMat<usize, C64>>,
    n: usize,
}

impl ComplexCsc {
    /// Assemble a matrix on `pattern`. `values_in_entry_order[k]` belongs to
    /// entry `k` as given to [`Pattern::from_entries`]; permutation to CSC
    /// order happens here.
    pub fn from_entries(pattern: &Pattern, values_in_entry_order: &[C64]) -> Result<Self, Error> {
        if values_in_entry_order.len() != pattern.nnz {
            return Err(Error::Shape {
                expected: pattern.nnz,
                got: values_in_entry_order.len(),
            });
        }
        let inner = faer::sparse::SparseColMat::new_from_order_and_values(
            (*pattern.symbolic).clone(),
            &pattern.order,
            values_in_entry_order,
        )
        .map_err(|e| Error::Backend(e.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
            n: pattern.n,
        })
    }

    /// Build from full triplets (duplicates allowed upstream but pre-summed
    /// here defensively).
    pub fn from_triplets(n: usize, triplets: &[(usize, usize, C64)]) -> Result<Self, Error> {
        faer::sparse::SparseColMat::<usize, C64>::try_new_from_triplets(n, n, triplets)
            .map(|m| Self {
                inner: Arc::new(m),
                n,
            })
            .map_err(|e| Error::Backend(e.to_string()))
    }

    pub fn nrows(&self) -> usize {
        self.n
    }

    fn as_faer_ref(&self) -> faer::sparse::SparseColMatRef<'_, usize, C64> {
        (*self.inner).as_ref()
    }

    /// Dense reconstruction for tests and small systems.
    pub fn to_dense(&self) -> Vec<Vec<C64>> {
        let n = self.n;
        let mut dense = vec![vec![C64_ZERO; n]; n];
        let sym = self.inner.symbolic();
        let ptr = sym.col_ptrs();
        let row_ind = sym.row_indices();
        let values = self.inner.values();
        for col in 0..n {
            for k in ptr[col]..ptr[col + 1] {
                dense[row_ind[k]][col] = values[k];
            }
        }
        dense
    }
}

/// Symbolic LU factorization of a [`Pattern`] (once per topology).
#[derive(Clone)]
pub struct SymbolicLu {
    inner: faer::sparse::linalg::solvers::SymbolicLu<usize>,
}

impl SymbolicLu {
    pub fn try_new(pattern: &Pattern) -> Result<Self, Error> {
        Ok(Self {
            inner: faer::sparse::linalg::solvers::SymbolicLu::try_new((*pattern.symbolic).as_ref())
                .map_err(|e| Error::Backend(e.to_string()))?,
        })
    }
}

/// Numeric LU factors bound to a [`SymbolicLu`].
pub struct ComplexLu {
    inner: faer::sparse::linalg::solvers::Lu<usize, C64>,
}

impl ComplexLu {
    /// Numeric factorization reusing the symbolic analysis of the pattern.
    pub fn try_new_with_symbolic(sym: &SymbolicLu, mat: &ComplexCsc) -> Result<Self, Error> {
        Ok(Self {
            inner: faer::sparse::linalg::solvers::Lu::<usize, C64>::try_new_with_symbolic(
                sym.inner.clone(),
                mat.as_faer_ref(),
            )
            .map_err(|e| Error::Backend(e.to_string()))?,
        })
    }

    /// Solve `A x = rhs`.
    pub fn solve(&self, rhs: &[C64]) -> Vec<C64> {
        use faer::linalg::solvers::SpSolver as _;
        let mut rhs_mat = faer::mat::Mat::<C64>::from_fn(rhs.len(), 1, |i, _| rhs[i]);
        self.inner.solve_in_place(&mut rhs_mat);
        (0..rhs.len()).map(|i| rhs_mat[(i, 0)]).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_small_system() {
        // [[2, 1], [0, 3]] x = [4, 6]  =>  x = [1, 2]
        let entries = [(0usize, 0usize), (0, 1), (1, 1)];
        let vals = [
            C64 { re: 2.0, im: 0.0 },
            C64 { re: 1.0, im: 0.0 },
            C64 { re: 3.0, im: 0.0 },
        ];
        let pattern = Pattern::from_entries(2, &entries).unwrap();
        let a = ComplexCsc::from_entries(&pattern, &vals).unwrap();
        let sym = SymbolicLu::try_new(&pattern).unwrap();
        let lu = ComplexLu::try_new_with_symbolic(&sym, &a).unwrap();

        let x = lu.solve(&[C64 { re: 4.0, im: 0.0 }, C64 { re: 6.0, im: 0.0 }]);
        assert!((x[0].re - 1.0).abs() < 1e-12);
        assert!((x[1].re - 2.0).abs() < 1e-12);
    }

    #[test]
    fn complex_shifted_solve() {
        // (A - theta I) x = b with A = diag(1,2), theta = 3+4j
        let theta = C64 { re: 3.0, im: 4.0 };
        let entries = [(0usize, 0usize), (1usize, 1usize)];
        let vals = [
            C64 { re: 1.0, im: 0.0 } - theta,
            C64 { re: 2.0, im: 0.0 } - theta,
        ];
        let pattern = Pattern::from_entries(2, &entries).unwrap();
        let a = ComplexCsc::from_entries(&pattern, &vals).unwrap();
        let sym = SymbolicLu::try_new(&pattern).unwrap();
        let lu = ComplexLu::try_new_with_symbolic(&sym, &a).unwrap();

        let x = lu.solve(&[C64 { re: 1.0, im: 0.0 }, C64 { re: 1.0, im: 0.0 }]);
        // Verify by residual: A' x == b
        let ap = [
            (C64 { re: 1.0, im: 0.0 } - theta),
            (C64 { re: 2.0, im: 0.0 } - theta),
        ];
        let r0 = ap[0] * x[0];
        let r1 = ap[1] * x[1];
        assert!((r0.re - 1.0).abs() < 1e-10 && r0.im.abs() < 1e-10);
        assert!((r1.re - 1.0).abs() < 1e-10 && r1.im.abs() < 1e-10);
    }

    #[test]
    fn dense_reconstruction_matches_triplets() {
        let triplets = vec![
            (0usize, 1usize, C64 { re: 5.0, im: 1.0 }),
            (1usize, 0usize, C64 { re: -2.0, im: 0.0 }),
            (1usize, 1usize, C64 { re: 4.0, im: 0.0 }),
        ];
        let m = ComplexCsc::from_triplets(2, &triplets).unwrap();
        let d = m.to_dense();
        assert_eq!(d[0][1], C64 { re: 5.0, im: 1.0 });
        assert_eq!(d[1][0], C64 { re: -2.0, im: 0.0 });
        assert_eq!(d[0][0], C64_ZERO);
    }
}
