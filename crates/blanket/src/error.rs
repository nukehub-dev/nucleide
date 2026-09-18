//! Error type for the blanket crate.

use thiserror::Error;

/// Result alias for the `blanket` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors raised while validating bookkeeping inputs.
///
/// Every rejection names its cause; capabilities outside the arithmetic
/// scope (transport solving, geometry optimization, evaluated-data
/// consumption) have no spelling here by construction.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    /// A named input holds a non-finite value (`"tritons_bred"`,
    /// `"source_neutrons"`, `"raw_tbr"`, `"port_fractions"`,
    /// `"effective_tbr"`, `"required_tbr"`, `"blanket_power_mw"`,
    /// `"fusion_power_mw"`, `"multiplication"`).
    #[error("blanket: non-finite value in {0}")]
    NonFinite(&'static str),
    /// A quantity that must be non-negative holds a negative value.
    #[error("blanket: {0} must be non-negative")]
    Negative(&'static str),
    /// A quantity that must be strictly positive is not
    /// (`"source_neutrons"`, `"fusion_power_mw"`, `"required_tbr"`).
    #[error("blanket: {0} must be strictly positive")]
    NonPositive(&'static str),
    /// A port coverage fraction sits outside `[0, 1)`: entry `index`
    /// holds `value`.
    #[error("blanket: port fraction at index {index} is {value} (need 0 <= f < 1)")]
    InvalidPortFraction {
        /// Position of the offending fraction in the caller slice.
        index: usize,
        /// Offending fraction value.
        value: f64,
    },
}
