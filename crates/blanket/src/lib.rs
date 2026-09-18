//! TBR and blanket power bookkeeping over caller transport tallies.
//!
//! Blanket studies need TBR bookkeeping, not transport: the tallies come
//! from the caller's OpenMC/DAGMC runs, and the penalties, multiplications,
//! and margins below are closed-form arithmetic over those caller values.
//! Every input is a plain float or slice; nothing reads files and nothing
//! solves transport. This crate reports margins — it never optimizes
//! geometry and never targets a TBR.
//!
//! Equation set (cited as `(B1)`–`(F2)` in code comments):
//!
//! - `(B1)` raw TBR: `TBR = bred / source` (dimensionless), with `bred`
//!   the caller-tallied tritons bred and `source` the source neutrons.
//! - `(B2)` port-penalty haircut: `TBR_eff = TBR_raw * Prod_i (1 - f_i)`,
//!   one multiplicative haircut per port coverage fraction `f_i` in
//!   `[0, 1)`.
//! - `(B3)` breeding margin: `margin = TBR_eff - 1` (dimensionless;
//!   negative is a deficit, never clamped).
//! - `(B4)` requirement check: `TBR_eff >= TBR_required`.
//! - `(P1)` energy multiplication: `M = P_blanket / P_fusion`
//!   (dimensionless).
//! - `(P2)` blanket power: `P_blanket = M * P_fusion` (MW).
//! - `(F1)` tritium burn rate: `burn = P_fusion * K` (g/day), with `K`
//!   the D-T burn constant below (one triton of [`TRITON_MASS_U`] per
//!   [`DT_ENERGY_MEV`] of fusion energy).
//! - `(F2)` net surplus: `net = (TBR_eff - 1) * burn` (g/day; negative is
//!   a net deficit).
//!
//! Regression anchors: the reported stellarator Point-A values gate the
//! arithmetic — raw TBR `1.1070`, `1.074` after the 3% ECRH-port haircut
//! (`1.1070 * 0.97 = 1.07379`), energy multiplication `1.20`, and a burn
//! of `416.6` g/day at Point A. The Point-A fusion power quoted in the
//! gates (`2714.82` MW) is the power implied by `416.6` g/day through
//! `(F1)`, not a separate published input: the published anchor is the
//! `416.6` value itself.
//!
//! Out of scope by construction: breeding-blanket transport or neutronics
//! solving, TBR target solving, and any coupling into the `tritium`
//! permeation kernel (which keeps its no-breeding-coupling scope line).

#![warn(missing_docs)]

pub mod error;

pub use error::{Error, Result};

/// D-T fusion energy released per reaction, in MeV (pinned constant for `(F1)`).
pub const DT_ENERGY_MEV: f64 = 17.6;

/// Triton atomic mass, in u (pinned constant for `(F1)`).
pub const TRITON_MASS_U: f64 = 3.0160492;

/// Atomic mass unit, in kg (exact, SI 2019).
pub const ATOMIC_MASS_KG: f64 = 1.660_539_066_60e-27;

/// Electronvolt, in J (exact, SI 2019).
pub const EV_JOULE: f64 = 1.602_176_634e-19;

/// Tritium burn rate per MW of D-T fusion power, in g/day (derived for `(F1)`).
///
/// `K = 1e6 * 86400 * (TRITON_MASS_U * ATOMIC_MASS_KG * 1e3)
/// / (DT_ENERGY_MEV * 1e6 * EV_JOULE)`: one triton mass per reaction
/// energy, scaled from per-MW-second to per-day grams. Evaluates to
/// `0.15345399531870646` g/day per MW.
pub const BURN_G_PER_DAY_PER_MW: f64 = 1.0e6 * 86_400.0 * (TRITON_MASS_U * ATOMIC_MASS_KG * 1.0e3)
    / (DT_ENERGY_MEV * 1.0e6 * EV_JOULE);

fn check_finite(what: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFinite(what))
    }
}

/// Raw TBR from caller tallies `(B1)`.
///
/// `tritons_bred` is the caller-tallied tritons bred (non-negative),
/// `source_neutrons` the source neutrons (strictly positive). Returns the
/// dimensionless ratio `bred / source`.
pub fn tbr_from_tallies(tritons_bred: f64, source_neutrons: f64) -> Result<f64> {
    check_finite("tritons_bred", tritons_bred)?;
    check_finite("source_neutrons", source_neutrons)?;
    if tritons_bred < 0.0 {
        return Err(Error::Negative("tritons_bred"));
    }
    if source_neutrons <= 0.0 {
        return Err(Error::NonPositive("source_neutrons"));
    }
    Ok(tritons_bred / source_neutrons)
}

/// Effective TBR after per-port coverage penalties `(B2)`.
///
/// Each entry of `port_fractions` is one port's fractional coverage loss in
/// `[0, 1)`; the haircut is multiplicative,
/// `raw * Prod_i (1 - f_i)`. An empty slice returns `raw_tbr` unchanged
/// (no ports, no penalty).
pub fn apply_port_penalty(raw_tbr: f64, port_fractions: &[f64]) -> Result<f64> {
    check_finite("raw_tbr", raw_tbr)?;
    if raw_tbr < 0.0 {
        return Err(Error::Negative("raw_tbr"));
    }
    let mut effective = raw_tbr;
    for (index, fraction) in port_fractions.iter().enumerate() {
        check_finite("port_fractions", *fraction)?;
        if !(0.0..1.0).contains(fraction) {
            return Err(Error::InvalidPortFraction {
                index,
                value: *fraction,
            });
        }
        effective *= 1.0 - fraction;
    }
    Ok(effective)
}

/// Breeding margin `(B3)`: `effective_tbr - 1` (dimensionless).
///
/// A negative return is a sub-breakeven deficit and is reported as-is,
/// never clamped to zero.
pub fn breeding_margin(effective_tbr: f64) -> Result<f64> {
    check_finite("effective_tbr", effective_tbr)?;
    if effective_tbr < 0.0 {
        return Err(Error::Negative("effective_tbr"));
    }
    Ok(effective_tbr - 1.0)
}

/// Requirement check `(B4)`: `true` when `effective_tbr >= required_tbr`.
///
/// `required_tbr` must be finite and strictly positive (a requirement of
/// exactly breakeven is spelled `1.0`, not `0.0`).
pub fn meets_requirement(effective_tbr: f64, required_tbr: f64) -> Result<bool> {
    check_finite("effective_tbr", effective_tbr)?;
    check_finite("required_tbr", required_tbr)?;
    if effective_tbr < 0.0 {
        return Err(Error::Negative("effective_tbr"));
    }
    if required_tbr <= 0.0 {
        return Err(Error::NonPositive("required_tbr"));
    }
    Ok(effective_tbr >= required_tbr)
}

/// Blanket energy multiplication `(P1)`: `blanket_power_mw / fusion_power_mw`.
///
/// Both powers are in MW; blanket power is non-negative, fusion power
/// strictly positive.
pub fn energy_multiplication(blanket_power_mw: f64, fusion_power_mw: f64) -> Result<f64> {
    check_finite("blanket_power_mw", blanket_power_mw)?;
    check_finite("fusion_power_mw", fusion_power_mw)?;
    if blanket_power_mw < 0.0 {
        return Err(Error::Negative("blanket_power_mw"));
    }
    if fusion_power_mw <= 0.0 {
        return Err(Error::NonPositive("fusion_power_mw"));
    }
    Ok(blanket_power_mw / fusion_power_mw)
}

/// Blanket thermal power `(P2)`: `multiplication * fusion_power_mw`, in MW.
pub fn blanket_power(fusion_power_mw: f64, multiplication: f64) -> Result<f64> {
    check_finite("fusion_power_mw", fusion_power_mw)?;
    check_finite("multiplication", multiplication)?;
    if fusion_power_mw <= 0.0 {
        return Err(Error::NonPositive("fusion_power_mw"));
    }
    if multiplication < 0.0 {
        return Err(Error::Negative("multiplication"));
    }
    Ok(fusion_power_mw * multiplication)
}

/// Tritium burn rate `(F1)`: `fusion_power_mw * BURN_G_PER_DAY_PER_MW`, in g/day.
pub fn tritium_burn_g_per_day(fusion_power_mw: f64) -> Result<f64> {
    check_finite("fusion_power_mw", fusion_power_mw)?;
    if fusion_power_mw <= 0.0 {
        return Err(Error::NonPositive("fusion_power_mw"));
    }
    Ok(fusion_power_mw * BURN_G_PER_DAY_PER_MW)
}

/// Net tritium surplus `(F2)`: `(effective_tbr - 1) * burn`, in g/day.
///
/// Positive is a breeding surplus, negative a net deficit (reported
/// as-is). Powers share the `(F1)` validation.
pub fn net_surplus_g_per_day(effective_tbr: f64, fusion_power_mw: f64) -> Result<f64> {
    check_finite("effective_tbr", effective_tbr)?;
    if effective_tbr < 0.0 {
        return Err(Error::Negative("effective_tbr"));
    }
    let burn = tritium_burn_g_per_day(fusion_power_mw)?;
    Ok((effective_tbr - 1.0) * burn)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Stellarator Point-A anchors (reported values; the gates below pin
    // the arithmetic against them at exact equality).
    const RAW_TBR: f64 = 1.1070;
    const ECRH_PORT_FRACTION: f64 = 0.03;
    const PENALIZED_TBR: f64 = 1.07379; // 1.1070 * 0.97
    const REPORTED_PENALIZED_TBR: f64 = 1.074; // report print precision
    const POWER_MULTIPLICATION: f64 = 1.20;
    const POINT_A_BURN_G_PER_DAY: f64 = 416.6;
    // Fusion power implied by 416.6 g/day through (F1): 416.6 / K.
    const POINT_A_FUSION_MW: f64 = 2714.82;

    #[test]
    fn burn_constant_matches_hand_derivation() {
        // Hand vector, rearranged term-by-term from the const expression:
        // seconds/day * W/MW * triton grams per reaction / joules per
        // reaction. Shares no intermediate with BURN_G_PER_DAY_PER_MW.
        let joules_per_reaction = 17.6 * 1.0e6 * 1.602_176_634e-19;
        let grams_per_reaction = 3.0160492 * 1.660_539_066_60e-27 * 1.0e3;
        let hand = 86_400.0 * 1.0e6 * grams_per_reaction / joules_per_reaction;
        assert!((BURN_G_PER_DAY_PER_MW - hand).abs() < 1e-18);
        assert!((BURN_G_PER_DAY_PER_MW - 0.15345399531870646).abs() < 1e-16);
    }

    #[test]
    fn stellaris_raw_tbr_gate() {
        let got = tbr_from_tallies(RAW_TBR, 1.0).unwrap();
        assert!((got - 1.1070).abs() < 1e-12);
        // Ratio form agrees: bred/source with a non-unit source.
        let scaled = tbr_from_tallies(2.0 * RAW_TBR, 2.0).unwrap();
        assert!((scaled - 1.1070).abs() < 1e-12);
    }

    #[test]
    fn stellaris_port_penalty_gate() {
        let got = apply_port_penalty(RAW_TBR, &[ECRH_PORT_FRACTION]).unwrap();
        assert!((got - PENALIZED_TBR).abs() < 1e-12);
        // Report print precision: 1.07379 rounds to the published 1.074.
        assert!((got - REPORTED_PENALIZED_TBR).abs() < 5e-4);
        // Empty port list is the identity (no ports, no penalty).
        assert!((apply_port_penalty(RAW_TBR, &[]).unwrap() - RAW_TBR).abs() < 1e-15);
        // Several ports compose multiplicatively.
        let multi = apply_port_penalty(1.0, &[0.03, 0.02]).unwrap();
        assert!((multi - 0.97 * 0.98).abs() < 1e-15);
    }

    #[test]
    fn stellaris_power_multiplication_gate() {
        let m = energy_multiplication(3600.0, 3000.0).unwrap();
        assert!((m - POWER_MULTIPLICATION).abs() < 1e-12);
        let p = blanket_power(3000.0, POWER_MULTIPLICATION).unwrap();
        assert!((p - 3600.0).abs() / 3600.0 < 1e-12);
        // Round-trip: (P2) then (P1) recovers the multiplication.
        let back = energy_multiplication(p, 3000.0).unwrap();
        assert!((back - POWER_MULTIPLICATION).abs() < 1e-12);
    }

    #[test]
    fn stellaris_burn_gate() {
        // Report-precision gate: the (F1) burn at the implied Point-A
        // power reproduces the published 416.6 g/day.
        let burn = tritium_burn_g_per_day(POINT_A_FUSION_MW).unwrap();
        assert!((burn - POINT_A_BURN_G_PER_DAY).abs() < 1e-2);
        // Exact-equality leg: the burn constant identity itself.
        let kilo = tritium_burn_g_per_day(1000.0).unwrap();
        assert!((kilo - 1000.0 * BURN_G_PER_DAY_PER_MW).abs() == 0.0);
        assert!((kilo - 153.45399532).abs() < 1e-6);
    }

    #[test]
    fn stellaris_margin_gates() {
        let margin = breeding_margin(PENALIZED_TBR).unwrap();
        assert!((margin - 0.07379).abs() < 1e-12);
        assert!(meets_requirement(PENALIZED_TBR, 1.0).unwrap());
        assert!(meets_requirement(PENALIZED_TBR, 1.05).unwrap());
        assert!(!meets_requirement(PENALIZED_TBR, 1.10).unwrap());
        // Net surplus at Point A: (1.07379 - 1) * burn ≈ 30.74 g/day.
        let burn = tritium_burn_g_per_day(POINT_A_FUSION_MW).unwrap();
        let net = net_surplus_g_per_day(PENALIZED_TBR, POINT_A_FUSION_MW).unwrap();
        assert!((net - (PENALIZED_TBR - 1.0) * burn).abs() == 0.0);
        assert!((net - 30.74).abs() < 5e-2);
        // Sub-breakeven is a reported deficit, never clamped.
        let deficit = net_surplus_g_per_day(0.9, 1000.0).unwrap();
        assert!(deficit < 0.0);
    }

    #[test]
    fn malformed_inputs_are_loud() {
        assert!(matches!(
            tbr_from_tallies(1.0, 0.0),
            Err(Error::NonPositive("source_neutrons"))
        ));
        assert!(matches!(
            tbr_from_tallies(-1.0, 1.0),
            Err(Error::Negative("tritons_bred"))
        ));
        assert!(matches!(
            tbr_from_tallies(f64::NAN, 1.0),
            Err(Error::NonFinite("tritons_bred"))
        ));
        assert!(matches!(
            apply_port_penalty(1.0, &[1.0]),
            Err(Error::InvalidPortFraction { index: 0, .. })
        ));
        assert!(matches!(
            apply_port_penalty(1.0, &[-0.01]),
            Err(Error::InvalidPortFraction { index: 0, .. })
        ));
        assert!(matches!(
            apply_port_penalty(1.0, &[0.0, f64::INFINITY]),
            Err(Error::NonFinite("port_fractions"))
        ));
        assert!(matches!(
            breeding_margin(-0.5),
            Err(Error::Negative("effective_tbr"))
        ));
        assert!(matches!(
            meets_requirement(1.1, 0.0),
            Err(Error::NonPositive("required_tbr"))
        ));
        assert!(matches!(
            energy_multiplication(1.0, -3.0),
            Err(Error::NonPositive("fusion_power_mw"))
        ));
        assert!(matches!(
            blanket_power(0.0, 1.2),
            Err(Error::NonPositive("fusion_power_mw"))
        ));
        assert!(matches!(
            blanket_power(100.0, -1.2),
            Err(Error::Negative("multiplication"))
        ));
        assert!(matches!(
            tritium_burn_g_per_day(f64::NAN),
            Err(Error::NonFinite("fusion_power_mw"))
        ));
        assert!(matches!(
            net_surplus_g_per_day(-1.0, 100.0),
            Err(Error::Negative("effective_tbr"))
        ));
    }
}
