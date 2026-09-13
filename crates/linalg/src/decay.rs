//! Decay-only perturbation consumer for the UQ-lite sampling kernel.
//!
//! Applies caller-supplied perturbation vectors to caller-supplied nominal
//! decay data. This module owns no nuclear data: branch fractions and decay
//! energies arrive as plain slices (the evaluated stores live in
//! `nucleide-nuclei` behind `LazyLock` + `include_str!`, and this crate must
//! never depend on them — the layering runs the other way).
//!
//! Branch fractions ([`perturb_branches`](crate::decay::perturb_branches)): the evaluated store drops
//! spontaneous-fission/fission branches at generation time, so kept branches
//! sum to `1 − BR(SF)` (see `nucleide-nuclei` `data.rs`: the `decay_branches`
//! docs and the `decay_branch_table` contract). Perturbed branches are
//! renormalised to preserve exactly that incoming deficit — the perturber
//! never invents or removes the dropped-SF mass, it only redistributes the
//! kept share. Negative raw branches clamp to zero before renormalisation
//! (the negative-clamp half of SANDY `Samples.truncate_normal(mode=1)` —
//! that routine additionally caps values `> 2` at `2` for relative-unit
//! samples centered at 1, which has no analogue here since branches
//! renormalise to the kept share); an all-clamped draw is a [`DecayError`](crate::decay::DecayError),
//! never a silent zero vector.
//!
//! Decay energies ([`perturb_energies`](crate::decay::perturb_energies)) are independent per nuclide (no sum
//! constraint) and follow the [`PerturbConvention`](crate::sample::PerturbConvention)
//! elementwise, with negative results clamped to zero as above.
//!
//! Fission-yield perturbation is explicitly OUT: no ENDF fission-yield
//! tapes are vendored or read anywhere in this workspace, so
//! [`perturb_fission_yields`](crate::decay::perturb_fission_yields) is a named-open hook that always returns
//! [`DecayError::FissionYieldsOpen`](crate::decay::DecayError). Next step is
//! landing a tape reader, then wiring real FY blocks through this hook.

use crate::sample::PerturbConvention;

/// Errors surfaced by the decay perturbation consumer.
#[derive(Debug, Clone, PartialEq)]
pub enum DecayError {
    /// An input slice is empty.
    Empty,
    /// A perturbation vector does not match the nominal length.
    LengthMismatch {
        /// Nominal (base) length.
        expected: usize,
        /// Perturbation length actually supplied.
        got: usize,
    },
    /// A non-finite value in the named input.
    NonFinite(&'static str),
    /// A nominal branch fraction or energy is negative at the given index.
    NegativeBase {
        /// Offending index.
        index: usize,
    },
    /// Renormalisation is impossible: the kept-branch sum is zero, or every
    /// perturbed branch clamped to zero.
    Degenerate,
    /// Fission-yield perturbation is named-open (waits on ENDF
    /// fission-yield tapes); no FY perturbation is performed.
    FissionYieldsOpen,
}

impl std::fmt::Display for DecayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecayError::Empty => write!(f, "decay perturbation input is empty"),
            DecayError::LengthMismatch { expected, got } => write!(
                f,
                "decay perturbation length mismatch: expected {expected}, got {got}"
            ),
            DecayError::NonFinite(what) => {
                write!(
                    f,
                    "decay perturbation input `{what}` holds a non-finite value"
                )
            }
            DecayError::NegativeBase { index } => {
                write!(f, "decay nominal value at index {index} is negative")
            }
            DecayError::Degenerate => write!(
                f,
                "decay branch renormalisation is degenerate (zero kept sum)"
            ),
            DecayError::FissionYieldsOpen => write!(
                f,
                "fission-yield perturbation is named-open (waits on ENDF fission-yield tapes)"
            ),
        }
    }
}

impl std::error::Error for DecayError {}

fn check_pair(base: &[f64], delta: &[f64], what: &'static str) -> Result<(), DecayError> {
    if base.is_empty() {
        return Err(DecayError::Empty);
    }
    if base.len() != delta.len() {
        return Err(DecayError::LengthMismatch {
            expected: base.len(),
            got: delta.len(),
        });
    }
    if base.iter().any(|v| !v.is_finite()) {
        return Err(DecayError::NonFinite("base"));
    }
    if delta.iter().any(|v| !v.is_finite()) {
        return Err(DecayError::NonFinite(what));
    }
    if let Some(index) = base.iter().position(|v| *v < 0.0) {
        return Err(DecayError::NegativeBase { index });
    }
    Ok(())
}

/// Perturb one parent's kept branch fractions with relative deltas.
///
/// Implements theory (U5): `raw[i] = base[i] * (1 + rel[i])`, negatives clamped to zero, then
/// rescaled so `sum(out) == sum(base)`: the incoming `1 − BR(SF)` deficit is
/// preserved exactly. Returns [`DecayError::Degenerate`] when the kept sum
/// is zero or every branch clamps to zero.
pub fn perturb_branches(base: &[f64], rel: &[f64]) -> Result<Vec<f64>, DecayError> {
    check_pair(base, rel, "rel")?;
    let kept: f64 = base.iter().sum();
    if kept <= 0.0 {
        return Err(DecayError::Degenerate);
    }
    let raw: Vec<f64> = base
        .iter()
        .zip(rel.iter())
        .map(|(b, r)| (b * (1.0 + r)).max(0.0))
        .collect();
    let total: f64 = raw.iter().sum();
    if total <= 0.0 {
        return Err(DecayError::Degenerate);
    }
    Ok(raw.iter().map(|v| v * kept / total).collect())
}

/// Perturb decay energies (mean prompt recoverable MeV/decay or any
/// caller-supplied per-nuclide energies) under `convention`.
///
/// Energies carry no sum constraint; each entry is perturbed independently
/// and negative results clamp to zero (the negative-clamp half of SANDY
/// `Samples.truncate_normal(mode=1)`; that routine additionally caps values
/// `> 2` at `2` for relative-unit samples centered at 1, which has no
/// analogue here — energies are absolute and carry no upper bound).
pub fn perturb_energies(
    base: &[f64],
    delta: &[f64],
    convention: PerturbConvention,
) -> Result<Vec<f64>, DecayError> {
    check_pair(base, delta, "delta")?;
    Ok(base
        .iter()
        .zip(delta.iter())
        .map(|(b, d)| convention.apply(*b, *d).max(0.0))
        .collect())
}

/// Passthrough for perturbation vectors that need no decay physics (e.g.
/// sampled group-wise deltas forwarded to a downstream consumer): validates
/// finiteness and returns a copy.
pub fn passthrough(delta: &[f64]) -> Result<Vec<f64>, DecayError> {
    if delta.is_empty() {
        return Err(DecayError::Empty);
    }
    if delta.iter().any(|v| !v.is_finite()) {
        return Err(DecayError::NonFinite("delta"));
    }
    Ok(delta.to_vec())
}

/// Fission-yield perturbation — named-open hook, always errors.
///
/// No ENDF fission-yield tapes exist on disk, so there is no evaluated
/// nominal to perturb and no FY covariance contract to sample from.
/// Next step: land a tape reader, then replace this stub with a real FY
/// perturber taking caller-supplied FY blocks (same deficit discipline as
/// [`perturb_branches`] where the evaluation defines one).
pub fn perturb_fission_yields(_base: &[f64], _rel: &[f64]) -> Result<Vec<f64>, DecayError> {
    Err(DecayError::FissionYieldsOpen)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn branches_preserve_the_sf_deficit() {
        // Synthetic kept branches summing to 0.90 (deficit 0.10 standing in
        // for a dropped SF branch — fixtures/uq/decay_perturb.json; the real
        // convention is documented on nucleide-nuclei `decay_branch_table`).
        let base = vec![0.50, 0.30, 0.10];
        let rel = vec![0.10, -0.20, 0.0];
        let out = perturb_branches(&base, &rel).unwrap();
        let kept: f64 = out.iter().sum();
        assert!((kept - 0.90).abs() < 1e-15, "deficit moved: sum = {kept}");
        // raw = [0.55, 0.24, 0.10], total 0.89, rescaled to 0.90.
        assert!((out[0] - 0.55 * 0.90 / 0.89).abs() < 1e-12);
        assert!((out[1] - 0.24 * 0.90 / 0.89).abs() < 1e-12);
        assert!((out[2] - 0.10 * 0.90 / 0.89).abs() < 1e-12);
    }

    #[test]
    fn branches_clamp_negatives_then_renormalise() {
        let base = vec![0.6, 0.3];
        let out = perturb_branches(&base, &[0.0, -1.5]).unwrap();
        assert!((out[0] - 0.9).abs() < 1e-12 && out[1] == 0.0);
        assert!(perturb_branches(&base, &[-2.0, -3.0]).is_err());
    }

    #[test]
    fn energies_follow_conventions_and_clamp() {
        let base = vec![0.5, 1.0];
        assert_eq!(
            perturb_energies(&base, &[0.2, -0.5], PerturbConvention::Relative).unwrap(),
            vec![0.6, 0.5]
        );
        assert_eq!(
            perturb_energies(&base, &[0.1, -2.0], PerturbConvention::Absolute).unwrap(),
            vec![0.6, 0.0]
        );
    }

    #[test]
    fn passthrough_copies_valid_vectors() {
        let d = vec![0.1, -0.2, 0.0];
        assert_eq!(passthrough(&d).unwrap(), d);
        assert_eq!(passthrough(&[]).unwrap_err(), DecayError::Empty);
        assert_eq!(
            passthrough(&[f64::INFINITY]).unwrap_err(),
            DecayError::NonFinite("delta")
        );
    }

    #[test]
    fn fission_yields_stay_named_open() {
        assert_eq!(
            perturb_fission_yields(&[0.5], &[0.1]).unwrap_err(),
            DecayError::FissionYieldsOpen
        );
    }

    #[test]
    fn malformed_inputs_name_their_cause() {
        assert_eq!(perturb_branches(&[], &[]).unwrap_err(), DecayError::Empty);
        assert_eq!(
            perturb_branches(&[0.5], &[0.1, 0.2]).unwrap_err(),
            DecayError::LengthMismatch {
                expected: 1,
                got: 2
            }
        );
        assert_eq!(
            perturb_branches(&[0.5, -0.1], &[0.0, 0.0]).unwrap_err(),
            DecayError::NegativeBase { index: 1 }
        );
        assert_eq!(
            perturb_branches(&[0.0, 0.0], &[0.1, 0.2]).unwrap_err(),
            DecayError::Degenerate
        );
        assert_eq!(
            perturb_energies(&[f64::NAN], &[0.0], PerturbConvention::Absolute).unwrap_err(),
            DecayError::NonFinite("base")
        );
    }

    #[test]
    fn fixture_decay_inputs_behave() {
        // Same file the Python tests read (fixtures/uq/decay_perturb.json).
        let text = include_str!("../../../fixtures/uq/decay_perturb.json");
        let v: serde_json::Value = serde_json::from_str(text).unwrap();
        let base: Vec<f64> = serde_json::from_value(v["base_branches"].clone()).unwrap();
        let rel: Vec<f64> = serde_json::from_value(v["rel"].clone()).unwrap();
        let out = perturb_branches(&base, &rel).unwrap();
        let (kept, total): (f64, f64) = (base.iter().sum(), out.iter().sum());
        assert!((kept - total).abs() < 1e-15);
        let energies: Vec<f64> = serde_json::from_value(v["base_energies"].clone()).unwrap();
        let delta: Vec<f64> = serde_json::from_value(v["delta"].clone()).unwrap();
        let conv = PerturbConvention::parse(v["convention"].as_str().unwrap()).unwrap();
        let e = perturb_energies(&energies, &delta, conv).unwrap();
        assert!(e.iter().all(|x| *x >= 0.0));
    }
}
