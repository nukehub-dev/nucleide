//! Neutron spectrum unfolding from activation-type measurements.
//!
//! Foil-activation spectrometry measures a handful of detector/reaction
//! rates; unfolding adjusts a caller-supplied guess spectrum until folding
//! it through the caller-supplied response matrix reproduces those rates
//! (JET/ITER diagnostics workhorse). Inverse-problem iteration, not a
//! transport solve (the kinetics θ-method precedent).
//!
//! One method lands per cycle, each in its own module behind the same
//! iteration shape. Shipped so far: SAND-II ([`sandii`], McElroy et al.,
//! AFWL-TR-67-41, 1967 — US government work, public domain, clean-room
//! from the report), STAYSL-class damped least-squares ([`staysl`],
//! Perey, ORNL/TM-6062, 1977 — US government work, clean-room) with
//! caller-supplied per-detector sigmas as weights, every solve routed
//! through the shared [`nucleide_linalg::lstsq`] kernel, and GRAVEL
//! ([`gravel`], Matzke, PTB-N-19, 1994 — the SAND-II update with
//! measurement-error weights), and MAXED ([`maxed`], Reginatto &
//! Goldhagen, Health Phys. 77 (1999) 579 — maximum entropy over the
//! exponential Lagrange family with a chi-square target).
//!
//! Data provenance is pinned: every response value, measured rate, guess
//! entry, and energy-group bound is caller-supplied. Evaluated libraries
//! (IRDFF and other IAEA-copyright data) are never vendored; the IRDFF-II
//! v1 pack ships as a runtime download (`nucleide.data.fetch_irdff` plus
//! `parse_irdff_g725`, parsed into caller-ready response rows). Test gates are synthetic forward-fold-then-recover
//! round-trips plus the IRDFF-II analytical benchmark-field *shapes*
//! (Trkov et al., Nucl. Data Sheets 163 (2020) 1), which are plain
//! published facts used with citation.
//!
//! Convergence contract (all methods): a per-group relative-change
//! tolerance with an explicit iteration cap; exhausting the cap is a hard
//! [`Error::NotConverged`], never a silent partial spectrum (the `tritium`
//! face-Newton precedent).
//!
//! Modules: [`error`] (named error set), [`sandii`] (SAND-II adjustment
//! iterator + driver), [`staysl`] (STAYSL-class damped least-squares
//! iterator + driver), [`gravel`] (GRAVEL chi-square-weighted adjustment
//! iterator + driver), [`maxed`] (MAXED maximum-entropy dual adjustment
//! iterator + driver). [`forward_fold`] is the shared forward operator.
//! [`staysl`], [`gravel`], and [`maxed`] keep their `Iteration`/`Solution`
//! types module-scoped: the names intentionally mirror [`sandii`]'s, which
//! stay re-exported at the crate root for compatibility.

#![warn(missing_docs)]

pub mod error;
pub mod gravel;
pub mod maxed;
pub mod sandii;
pub mod staysl;

pub use error::{Error, Result};
pub use sandii::{Iteration, SandII, Solution, DEFAULT_MAX_ITERATIONS, DEFAULT_TOLERANCE};

use error::Error as UnfoldError;

/// Validate a response matrix: non-empty, ragged-free, finite, non-negative.
/// Returns the group count `n` (the common row width).
pub(crate) fn validate_response(response: &[Vec<f64>]) -> Result<usize> {
    if response.is_empty() {
        return Err(UnfoldError::BadResponse("empty response matrix"));
    }
    let n = response[0].len();
    if n == 0 {
        return Err(UnfoldError::BadResponse("empty response rows"));
    }
    for row in response {
        if row.len() != n {
            return Err(UnfoldError::BadShape {
                what: "response row",
                expected: n,
                got: row.len(),
            });
        }
    }
    if response.iter().flatten().any(|v| !v.is_finite()) {
        return Err(UnfoldError::BadResponse("non-finite entry"));
    }
    if response.iter().flatten().any(|v| *v < 0.0) {
        return Err(UnfoldError::BadResponse("negative entry"));
    }
    Ok(n)
}

/// Validate a per-group spectrum vector against the expected group count.
/// With `strict`, entries must be strictly positive (guess spectra); value
/// problems are wrapped by the caller's `bad_value` constructor so the error
/// names the input being validated.
pub(crate) fn validate_spectrum(
    what: &'static str,
    spectrum: &[f64],
    expected: usize,
    strict: bool,
    bad_value: fn(&'static str) -> UnfoldError,
) -> Result<()> {
    if spectrum.len() != expected {
        return Err(UnfoldError::BadShape {
            what,
            expected,
            got: spectrum.len(),
        });
    }
    if spectrum.iter().any(|v| !v.is_finite()) {
        return Err(bad_value("non-finite entry"));
    }
    if spectrum.iter().any(|v| *v <= 0.0) {
        return Err(bad_value(if strict {
            "entries must be strictly positive"
        } else {
            "entries must be non-negative"
        }));
    }
    Ok(())
}

/// Forward operator: fold a spectrum through the response matrix (S1).
///
/// `response` holds one row per detector (all rows one value per energy
/// group); returns the calculated rate per detector. This is the map the
/// unfolding adjusts against — also the natural way to synthesize
/// round-trip test rates from a known spectrum.
pub fn forward_fold(response: &[Vec<f64>], spectrum: &[f64]) -> Result<Vec<f64>> {
    let n = validate_response(response)?;
    validate_spectrum("spectrum", spectrum, n, false, UnfoldError::BadResponse)?;
    Ok(fold(response, spectrum))
}

/// Unchecked dense fold over validated inputs.
pub(crate) fn fold(response: &[Vec<f64>], spectrum: &[f64]) -> Vec<f64> {
    response
        .iter()
        .map(|row| row.iter().zip(spectrum).map(|(r, p)| r * p).sum())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_fold_validates_its_inputs() {
        let response = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        assert!(matches!(
            forward_fold(&response, &[1.0]),
            Err(UnfoldError::BadShape {
                what: "spectrum",
                expected: 2,
                got: 1
            })
        ));
        let ragged = vec![vec![1.0], vec![1.0, 2.0]];
        assert!(matches!(
            forward_fold(&ragged, &[1.0, 2.0]),
            Err(UnfoldError::BadShape {
                what: "response row",
                expected: 1,
                got: 2
            })
        ));
        let neg = vec![vec![-1.0, 2.0]];
        assert!(matches!(
            forward_fold(&neg, &[1.0, 1.0]),
            Err(UnfoldError::BadResponse("negative entry"))
        ));
    }
}
