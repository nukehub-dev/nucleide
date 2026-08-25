//! Multicomponent enrichment cascade solving and separative-work analytics.
//!
//! The numeric solver is provided; generated symbolic solvers are
//! deliberately out of scope.
//!
//! # Model
//!
//! A [`Cascade`] describes a matched-abundance-ratio cascade (MARC): a train
//! of ideal stages whose overall stage separation factor `alpha` and mass
//! separation factor `M*` split an `n`-component mixture into a product
//! stream (P), enriched in the *j*-th key nuclide, and a tails stream (T),
//! enriched in the *k*-th key nuclide, from a common feed (F). Per-component
//! stage factors follow
//!
//! ```text
//! alpha*_i = alpha^(M* - M_i)
//! ```
//!
//! with `M*` chosen strictly between the atomic masses of the `j` and `k`
//! keys ([`multicomponent`] optimizes it to minimize the total flow rate
//! `L_t / F`). Stage-number relations and the flow-rate/separative-power
//! sums follow Wood, Borisevich & Sulaberidze, "On a Criterion Efficiency
//! for Multi-Isotope Mixtures Separation", Sep. Sci. Technol. 34:3
//! 343–357 (DOI 10.1081/SS-100100654).
//!
//! # Conventions
//!
//! - Feed/product/tails compositions are stored as normalized mass
//!   fractions keyed by [`NuclideId`], with the stream total mass carried
//!   alongside.
//! - This crate is materials-free on purpose: no dependency on the
//!   workspace `material` crate, so enrichment math stays usable standalone.
//! - Atomic masses come from the AME2020 tables in `nuclei::data`.
//!
//! # Error handling deviations from legacy implementations
//!
//! - Solvers report failures as [`Result`]s ([`Error::NoConvergence`],
//!   [`Error::IterationNaN`]) instead of returning silently-unconverged
//!   cascades or throwing `EnrichmentIterationLimit` /
//!   `EnrichmentIterationNaN`; internal fixed-point loops carry generous
//!   iteration caps so pathological inputs cannot hang the solver.
//! - `multicomponent` supports the numeric solver only (the symbolic one
//!   is upstream-generated code and is out of scope).
//! - The thin `feed()`/`product()`/`tails()`/`swu()` quantity wrappers are
//!   replaced by the underlying `*_per_*` mass-ratio primitives, which the
//!   caller multiplies by the known stream quantity.

mod cascade;
mod swu;

pub use cascade::{
    alphastar_i, default_uranium_cascade, feed_per_prod, feed_per_tail, multicomponent,
    prod_per_feed, prod_per_tail, solve_numeric, tail_per_feed, tail_per_prod, Cascade, Stream,
    DEFAULT_MAX_ITER, DEFAULT_TOLERANCE,
};
pub use swu::{swu_per_feed, swu_per_prod, swu_per_tail, value_func};

use nuclei::NuclideId;

/// Errors raised by the cascade solvers.
///
/// Reports solver failures: iteration limits, NaN states, and
/// input-validation failure modes.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// A solver failed to reach its convergence tolerance within the
    /// allotted iterations.
    NoConvergence {
        /// Number of iterations performed when the solver gave up.
        iterations: u32,
    },
    /// A solver iteration produced non-finite intermediate values.
    IterationNaN,
    /// A cascade definition is unusable (empty or degenerate composition,
    /// coincident assay targets, ...).
    BadComposition { detail: String },
    /// A nuclide in the cascade has no entry in the atomic-mass tables.
    MissingMass(NuclideId),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::NoConvergence { iterations } => {
                write!(
                    f,
                    "enrichment solver did not converge in {iterations} iterations"
                )
            }
            Error::IterationNaN => {
                write!(f, "enrichment solver iteration produced non-finite values")
            }
            Error::BadComposition { detail } => write!(f, "bad cascade composition: {detail}"),
            Error::MissingMass(id) => {
                write!(f, "no atomic mass available for nuclide {id}")
            }
        }
    }
}

impl std::error::Error for Error {}

/// Crate-local result alias.
pub type Result<T> = std::result::Result<T, Error>;
