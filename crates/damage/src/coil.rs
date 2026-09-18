//! Coil fast-fluence and lifetime bookkeeping over caller spectra.
//!
//! Magnet lifetime is arithmetic once spectra exist: accumulate fluence
//! (`flux × time`, spectrally weighted through the landed [`crate::fold`]
//! kernels for dpa or through the fast-flux sum here), compare each channel
//! against its caller-supplied limit, and report the weakest link. No new
//! physics, no transport solving, no magnetics/quench/structural analysis.
//!
//! Pinned accumulation rules (one per metric):
//!
//! - **Fast fluence** ([`fast_flux`]/[`fast_fluence`]): the piecewise-constant
//!   group sum `Σ flux[g]` over groups with `bounds[g + 1] > threshold_mev`,
//!   times `seconds` for fluence. The threshold is caller-supplied (a common
//!   fast-neutron choice is 0.1 MeV); when it cuts through a group the whole
//!   group counts — callers align the grid so the threshold sits on a group
//!   boundary when that conservatism matters. Units are caller-consistent:
//!   whatever area normalization `flux` uses (m⁻²s⁻¹ or cm⁻²s⁻¹), the
//!   fluence and the limit table share it.
//! - **dpa**: the landed [`crate::fold`] kernels (`nrt_dpa`/`arc_dpa` at
//!   `seconds = 1` is the dpa rate); this module consumes the resulting
//!   rates without re-weighting. See the `coil_lifetime` doctest.
//! - **History** ([`accumulate`]): piecewise-constant irradiation history,
//!   `Σ rates[i]·durations[i]` — the caller loops its per-interval spectra
//!   through the metric of choice and accumulates here.
//!
//! Pinned life rule ([`coil_lifetime`]/[`coil_remaining`]): weakest-link,
//! `min` over channels of `limit / rate` (or `(limit − accumulated) / rate`
//! for remaining life, clamped at zero once a limit is reached). A zero rate
//! never fails, so it is skipped — unless its limit is already reached, in
//! which case it reports zero. All rates zero (and nothing exceeded) is
//! infinite life (`seconds = INFINITY`, `limiting = None`). Limit tables are
//! always caller-supplied; published design numbers are validation gates,
//! never defaults (see the `coil_lifetime` doctest).
//!
//! [`CoilLifetime`] carries the life in seconds plus the limiting channel
//! index (`None` only for the all-zero-rate infinite case).

use crate::error::{Error, Result};

/// Seconds in one full-power year (365.25 d), the life unit of the coil gates.
pub const FPY_SECONDS: f64 = 365.25 * 24.0 * 3600.0;

/// Lifetime of a coil set: `seconds` to the first limit breach and the
/// `limiting` channel index (`None` only when every rate is zero, i.e.
/// infinite life).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoilLifetime {
    /// Seconds to the first breach (`f64::INFINITY` when no channel ages).
    pub seconds: f64,
    /// Index into the caller's limit/rate tables of the weakest link.
    pub limiting: Option<usize>,
}

/// Validate the shared `(flux, bounds)` shape for the fast-flux sum.
///
/// Returns the group count `G`.
fn validate_fast(flux: &[f64], bounds: &[f64], threshold_mev: f64) -> Result<usize> {
    if flux.is_empty() {
        return Err(Error::Empty);
    }
    if bounds.len() != flux.len() + 1 {
        return Err(Error::DimensionMismatch {
            what: "bounds",
            expected: flux.len() + 1,
            got: bounds.len(),
        });
    }
    if !threshold_mev.is_finite() {
        return Err(Error::NonFinite("threshold_mev"));
    }
    if threshold_mev < 0.0 {
        return Err(Error::Negative("threshold_mev"));
    }
    if !flux.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("flux"));
    }
    if flux.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("flux"));
    }
    if !bounds.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("bounds"));
    }
    for (i, w) in bounds.windows(2).enumerate() {
        if w[0] >= w[1] {
            return Err(Error::NonMonotonicBounds { index: i });
        }
    }
    Ok(flux.len())
}

/// Fast flux above `threshold_mev`: `Σ flux[g]` over groups with
/// `bounds[g + 1] > threshold_mev` (caller area units, e.g. m⁻²s⁻¹).
///
/// A threshold at or above the top boundary sums to exactly `0.0`; a
/// threshold below the bottom boundary sums the whole spectrum. A threshold
/// cutting through a group includes the whole group (documented
/// conservatism — align the grid when it matters).
///
/// ```rust
/// use nucleide_damage::coil::fast_flux;
/// let bounds = vec![0.0, 0.1, 1.0, 20.0];
/// let flux = vec![1.0e12, 2.0e12, 4.0e12];
/// assert_eq!(fast_flux(&flux, &bounds, 0.1).unwrap(), 6.0e12);
/// assert_eq!(fast_flux(&flux, &bounds, 1.0).unwrap(), 4.0e12);
/// assert_eq!(fast_flux(&flux, &bounds, 0.0).unwrap(), 7.0e12);
/// assert_eq!(fast_flux(&flux, &bounds, 20.0).unwrap(), 0.0);
/// ```
pub fn fast_flux(flux: &[f64], bounds: &[f64], threshold_mev: f64) -> Result<f64> {
    validate_fast(flux, bounds, threshold_mev)?;
    Ok(flux
        .iter()
        .zip(bounds.windows(2))
        .filter(|(_, w)| w[1] > threshold_mev)
        .map(|(f, _)| *f)
        .sum())
}

/// Fast fluence above `threshold_mev`: [`fast_flux`] held for `seconds`
/// (`fluence = flux × time`, caller area units).
///
/// ```rust
/// use nucleide_damage::coil::{fast_fluence, fast_flux};
/// let bounds = vec![0.0, 0.1, 1.0, 20.0];
/// let flux = vec![1.0e12, 2.0e12, 4.0e12];
/// let rate = fast_flux(&flux, &bounds, 0.1).unwrap();
/// assert_eq!(fast_fluence(&flux, &bounds, 0.1, 2.0).unwrap(), rate * 2.0);
/// ```
pub fn fast_fluence(flux: &[f64], bounds: &[f64], threshold_mev: f64, seconds: f64) -> Result<f64> {
    if !seconds.is_finite() || seconds <= 0.0 {
        return Err(Error::NonPositive("seconds"));
    }
    Ok(fast_flux(flux, bounds, threshold_mev)? * seconds)
}

/// Accumulate a piecewise-constant irradiation history: `Σ rates[i]·durations[i]`.
///
/// The caller folds each interval's spectrum through its metric of choice
/// ([`fast_flux`] or a [`crate::fold`] dpa rate) and accumulates the rate
/// history here. Zero rates and durations are fine; negatives are loud.
///
/// ```rust
/// use nucleide_damage::coil::accumulate;
/// assert_eq!(accumulate(&[1.0e13, 2.0e13], &[10.0, 5.0]).unwrap(), 2.0e14);
/// ```
pub fn accumulate(rates: &[f64], durations: &[f64]) -> Result<f64> {
    if rates.is_empty() {
        return Err(Error::Empty);
    }
    if durations.len() != rates.len() {
        return Err(Error::DimensionMismatch {
            what: "durations",
            expected: rates.len(),
            got: durations.len(),
        });
    }
    if !rates.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("rates"));
    }
    if !durations.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("durations"));
    }
    if rates.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("rates"));
    }
    if durations.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("durations"));
    }
    Ok(rates
        .iter()
        .zip(durations.iter())
        .map(|(r, d)| *r * *d)
        .sum())
}

/// Validate the shared limit/rate tables; returns the channel count.
fn validate_life(limits: &[f64], rates: &[f64]) -> Result<usize> {
    if limits.is_empty() {
        return Err(Error::Empty);
    }
    if rates.len() != limits.len() {
        return Err(Error::DimensionMismatch {
            what: "rates",
            expected: limits.len(),
            got: rates.len(),
        });
    }
    if !limits.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("limits"));
    }
    if !rates.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("rates"));
    }
    if limits.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("limits"));
    }
    if rates.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("rates"));
    }
    Ok(limits.len())
}

/// Weakest-link coil life: `min` over channels of `limits[i] / rates[i]`.
///
/// Limits are caller-supplied fluences (or dpa — any accumulated metric, as
/// long as each limit shares its rate's units); rates are per-second. A zero
/// rate is skipped (it never breaches); all zero is infinite life
/// (`seconds = INFINITY`, `limiting = None`). A zero limit reports zero life
/// at that channel.
///
/// Design-anchor gate (caller numbers, not defaults): a 99th-percentile coil
/// fast flux of 9.5e13 m⁻²s⁻¹ against fast-fluence limits of 3e22 m⁻² and
/// 1.5e23 m⁻² gives `3e22 / 9.5e13 ≈ 3.16e8 s ≈ 10.0 FPY`, limited by the
/// first channel:
///
/// ```rust
/// use nucleide_damage::coil::{coil_lifetime, FPY_SECONDS};
/// let life = coil_lifetime(&[3.0e22, 1.5e23], &[9.5e13, 9.5e13]).unwrap();
/// assert_eq!(life.limiting, Some(0));
/// assert!((life.seconds / FPY_SECONDS - 10.0).abs() < 0.1);
/// ```
pub fn coil_lifetime(limits: &[f64], rates: &[f64]) -> Result<CoilLifetime> {
    let n = validate_life(limits, rates)?;
    let mut best: Option<(f64, usize)> = None;
    for i in 0..n {
        if rates[i] == 0.0 {
            continue;
        }
        let life = limits[i] / rates[i];
        if best.is_none_or(|(b, _)| life < b) {
            best = Some((life, i));
        }
    }
    Ok(match best {
        Some((seconds, i)) => CoilLifetime {
            seconds,
            limiting: Some(i),
        },
        None => CoilLifetime {
            seconds: f64::INFINITY,
            limiting: None,
        },
    })
}

/// Weakest-link remaining life from accumulated damage:
/// `min` over channels of `(limits[i] − accumulated[i]) / rates[i]`,
/// clamped at zero once a limit is reached.
///
/// `accumulated` shares each limit's units (caller loops its history through
/// [`fast_fluence`]/[`accumulate`] or a [`crate::fold`] dpa fold). Zero-rate
/// channels are skipped unless already breached (then they report zero);
/// nothing breached and no live channel is infinite remaining life.
///
/// ```rust
/// use nucleide_damage::coil::coil_remaining;
/// let life = coil_remaining(&[10.0, 10.0], &[4.0, 9.0], &[2.0, 1.0]).unwrap();
/// assert_eq!((life.seconds, life.limiting), (1.0, Some(1)));
/// let spent = coil_remaining(&[10.0, 10.0], &[11.0, 5.0], &[1.0, 1.0]).unwrap();
/// assert_eq!((spent.seconds, spent.limiting), (0.0, Some(0)));
/// ```
pub fn coil_remaining(limits: &[f64], accumulated: &[f64], rates: &[f64]) -> Result<CoilLifetime> {
    if limits.is_empty() {
        return Err(Error::Empty);
    }
    if accumulated.len() != limits.len() {
        return Err(Error::DimensionMismatch {
            what: "accumulated",
            expected: limits.len(),
            got: accumulated.len(),
        });
    }
    if rates.len() != limits.len() {
        return Err(Error::DimensionMismatch {
            what: "rates",
            expected: limits.len(),
            got: rates.len(),
        });
    }
    if !limits.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("limits"));
    }
    if !accumulated.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("accumulated"));
    }
    if !rates.iter().all(|v| v.is_finite()) {
        return Err(Error::NonFinite("rates"));
    }
    if limits.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("limits"));
    }
    if accumulated.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("accumulated"));
    }
    if rates.iter().any(|v| *v < 0.0) {
        return Err(Error::Negative("rates"));
    }
    let mut best: Option<(f64, usize)> = None;
    for i in 0..limits.len() {
        let headroom = limits[i] - accumulated[i];
        if headroom <= 0.0 {
            if best.is_none_or(|(b, _)| 0.0 < b) {
                best = Some((0.0, i));
            }
            continue;
        }
        if rates[i] == 0.0 {
            continue;
        }
        let life = headroom / rates[i];
        if best.is_none_or(|(b, _)| life < b) {
            best = Some((life, i));
        }
    }
    Ok(match best {
        Some((seconds, i)) => CoilLifetime {
            seconds,
            limiting: Some(i),
        },
        None => CoilLifetime {
            seconds: f64::INFINITY,
            limiting: None,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::nrt_dpa;

    /// Three-group synthetic grid (MeV) — the fold tests' round bounds.
    fn bounds3() -> Vec<f64> {
        vec![0.0, 0.1, 1.0, 20.0]
    }

    #[test]
    fn fast_flux_threshold_pins_groups() {
        let bounds = bounds3();
        let flux = vec![1.0e12, 2.0e12, 4.0e12];
        // Groups count iff their upper edge clears the threshold.
        assert_eq!(fast_flux(&flux, &bounds, 0.0).unwrap(), 7.0e12);
        assert_eq!(fast_flux(&flux, &bounds, 0.1).unwrap(), 6.0e12);
        assert_eq!(fast_flux(&flux, &bounds, 1.0).unwrap(), 4.0e12);
        assert_eq!(fast_flux(&flux, &bounds, 20.0).unwrap(), 0.0);
        assert_eq!(fast_flux(&flux, &bounds, 25.0).unwrap(), 0.0);
        // A threshold cutting through a group includes the whole group
        // (documented conservatism): 0.5 sits inside [0.1, 1.0).
        assert_eq!(fast_flux(&flux, &bounds, 0.5).unwrap(), 6.0e12);
        // Fluence is the rate held for `seconds`, exactly.
        let rate = fast_flux(&flux, &bounds, 0.1).unwrap();
        assert_eq!(fast_fluence(&flux, &bounds, 0.1, 2.0).unwrap(), rate * 2.0);
    }

    #[test]
    fn stellaris_anchor_is_ten_fpy_class() {
        // Caller-supplied design numbers (gates, not defaults): 99th
        // percentile coil fast flux 9.5e13 m⁻²s⁻¹ against fast-fluence
        // limits 3e22 / 1.5e23 m⁻². Hand life: 3e22/9.5e13 s.
        let life = coil_lifetime(&[3.0e22, 1.5e23], &[9.5e13, 9.5e13]).unwrap();
        assert_eq!(life.limiting, Some(0));
        let hand = 3.0e22 / 9.5e13;
        assert!(
            (life.seconds - hand).abs() / hand < 1e-9,
            "life {} vs hand {hand}",
            life.seconds
        );
        let fpy = life.seconds / FPY_SECONDS;
        assert!((fpy - 10.0).abs() < 0.1, "life {fpy} FPY");
        // Ten full-power years at that flux lands on the first limit.
        let fluence = fast_fluence(&[9.5e13], &[0.0, 20.0], 0.1, 10.0 * FPY_SECONDS).unwrap();
        assert!(
            (fluence - 3.0e22).abs() / 3.0e22 < 1e-3,
            "10-FPY fluence {fluence}"
        );
        // The copper-class channel alone is ~50 FPY.
        let cu = coil_lifetime(&[1.5e23], &[9.5e13]).unwrap();
        assert!((cu.seconds / FPY_SECONDS - 50.0).abs() < 0.5);
    }

    #[test]
    fn dpa_rate_feeds_the_life_kernel() {
        // dpa spectral weighting stays in the landed fold: the rate is the
        // NRT fold at one second, the life kernel consumes it unchanged.
        let bounds = bounds3();
        let flux = vec![1.0e12, 2.0e12, 4.0e12];
        let resp = vec![100.0, 200.0, 50.0];
        let rate = nrt_dpa(&flux, &resp, &bounds, 1.0).unwrap();
        assert_eq!(rate, 1.0e-24 * 7.0e14);
        let life = coil_lifetime(&[7.0e-10], &[rate]).unwrap();
        assert_eq!(life.limiting, Some(0));
        assert!((life.seconds - 1.0).abs() < 1e-9, "life {}", life.seconds);
    }

    #[test]
    fn weakest_link_skips_zero_rates() {
        let life = coil_lifetime(&[10.0, 10.0], &[1.0, 2.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (5.0, Some(1)));
        // A zero rate never breaches, so the live channel limits.
        let life = coil_lifetime(&[10.0, 10.0], &[0.0, 2.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (5.0, Some(1)));
        // Nothing aging at all is infinite life with no limiting channel.
        let life = coil_lifetime(&[10.0, 10.0], &[0.0, 0.0]).unwrap();
        assert_eq!(life.seconds, f64::INFINITY);
        assert_eq!(life.limiting, None);
        // A zero limit is spent at once.
        let life = coil_lifetime(&[0.0, 10.0], &[1.0, 1.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (0.0, Some(0)));
    }

    #[test]
    fn remaining_clamps_at_spent_limits() {
        let life = coil_remaining(&[10.0, 10.0], &[4.0, 9.0], &[2.0, 1.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (1.0, Some(1)));
        // Breached channels report zero, first one wins ties.
        let life = coil_remaining(&[10.0, 10.0], &[11.0, 5.0], &[1.0, 1.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (0.0, Some(0)));
        // A breached zero-rate channel still reports zero …
        let life = coil_remaining(&[10.0, 10.0], &[11.0, 5.0], &[0.0, 1.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (0.0, Some(0)));
        // … while an unbreached zero-rate channel is skipped.
        let life = coil_remaining(&[10.0, 10.0], &[5.0, 5.0], &[0.0, 1.0]).unwrap();
        assert_eq!((life.seconds, life.limiting), (5.0, Some(1)));
        // Nothing live and nothing spent is infinite.
        let life = coil_remaining(&[10.0], &[5.0], &[0.0]).unwrap();
        assert_eq!(life.seconds, f64::INFINITY);
        assert_eq!(life.limiting, None);
    }

    #[test]
    fn history_accumulates_rate_slices() {
        assert_eq!(accumulate(&[1.0e13, 2.0e13], &[10.0, 5.0]).unwrap(), 2.0e14);
        assert_eq!(accumulate(&[0.0, 2.0e13], &[10.0, 0.0]).unwrap(), 0.0);
    }

    #[test]
    fn malformed_coil_inputs_name_their_cause() {
        let bounds = bounds3();
        let flux = vec![1.0e12, 2.0e12, 4.0e12];
        assert!(matches!(fast_flux(&[], &[], 0.1), Err(Error::Empty)));
        assert!(matches!(
            fast_flux(&flux, &[0.0, 1.0], 0.1),
            Err(Error::DimensionMismatch { what: "bounds", .. })
        ));
        assert!(matches!(
            fast_flux(&flux, &bounds, f64::NAN),
            Err(Error::NonFinite("threshold_mev"))
        ));
        assert!(matches!(
            fast_flux(&flux, &bounds, -0.1),
            Err(Error::Negative("threshold_mev"))
        ));
        assert!(matches!(
            fast_flux(&[-1.0, 1.0, 1.0], &bounds, 0.1),
            Err(Error::Negative("flux"))
        ));
        assert!(matches!(
            fast_fluence(&flux, &bounds, 0.1, 0.0),
            Err(Error::NonPositive("seconds"))
        ));
        assert!(matches!(
            fast_flux(&flux, &[0.0, 0.5, 0.5, 1.0], 0.1),
            Err(Error::NonMonotonicBounds { index: 1 })
        ));
        assert!(matches!(accumulate(&[], &[]), Err(Error::Empty)));
        assert!(matches!(
            accumulate(&[1.0, 2.0], &[1.0]),
            Err(Error::DimensionMismatch {
                what: "durations",
                ..
            })
        ));
        assert!(matches!(
            accumulate(&[1.0, -1.0], &[1.0, 1.0]),
            Err(Error::Negative("rates"))
        ));
        assert!(matches!(
            accumulate(&[1.0, 1.0], &[1.0, f64::INFINITY]),
            Err(Error::NonFinite("durations"))
        ));
        assert!(matches!(coil_lifetime(&[], &[]), Err(Error::Empty)));
        assert!(matches!(
            coil_lifetime(&[1.0, 2.0], &[1.0]),
            Err(Error::DimensionMismatch { what: "rates", .. })
        ));
        assert!(matches!(
            coil_lifetime(&[1.0], &[-1.0]),
            Err(Error::Negative("rates"))
        ));
        assert!(matches!(
            coil_lifetime(&[f64::NAN], &[1.0]),
            Err(Error::NonFinite("limits"))
        ));
        assert!(matches!(
            coil_remaining(&[1.0], &[2.0, 3.0], &[1.0]),
            Err(Error::DimensionMismatch {
                what: "accumulated",
                ..
            })
        ));
        assert!(matches!(
            coil_remaining(&[1.0], &[0.5], &[1.0, 2.0]),
            Err(Error::DimensionMismatch { what: "rates", .. })
        ));
        assert!(matches!(
            coil_remaining(&[1.0], &[-0.5], &[1.0]),
            Err(Error::Negative("accumulated"))
        ));
    }
}
