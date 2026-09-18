//! Closed-form ECRH accessibility kernel: cold electron-cyclotron resonance,
//! relativistic (Maxwell–Jüttner thermal) resonance shift, and the O1/X1
//! cut-off densities, located along a caller-supplied beamline.
//!
//! The kernel answers heating-access questions without any ray or beam
//! tracing: given a wave frequency plus the caller's magnetic-field and
//! density profiles tabulated along a ray path, it reports where the wave
//! meets the cyclotron resonance and where it meets a cut-off. The
//! per-port absorbed-power penalty itself lives downstream (blanket
//! bookkeeping consumes these positions as plain floats); this module
//! reports accessibility, never designs launchers.
//!
//! # Equation set (published facts, cited — never derived here)
//!
//! - **ECRH-1 — cold gyrofrequency.** `f_ce = e·B / (2π·m_e)` (Hutchinson,
//!   *Principles of Plasma Diagnostics*, 2nd ed. (2002), §7.2; the
//!   electron-cyclotron frequency; ≈ 27.992 GHz/T).
//! - **ECRH-2 — cold resonance.** Fundamental (harmonic `n`) absorption sits
//!   where `f = n·f_ce`, i.e. `B_res = 2π·m_e·f / (n·e)` (Bornatici et al.,
//!   Nucl. Fusion **23** (1983) 1153 — the ECRH review pinning the cold
//!   resonance condition alongside the cut-offs below).
//! - **ECRH-3 — relativistic shift.** At finite electron temperature the
//!   resonance condition carries the Lorentz factor, `f = n·f_ce/γ`, so the
//!   resonant field moves out to `B_res(T_e) = γ·B_res(cold)` with the
//!   thermal `γ = 1 + T_e/m_e·c²` (`T_e` in keV, `m_e·c² = 510.999 keV;
//!   Bornatici et al. 1983 — the weakly-relativistic Maxwell–Jüttner
//!   downshift of the absorption layer).
//! - **ECRH-4 — O-mode cut-off.** `P = 0` of the cold Appleton–Hartree
//!   dispersion gives `f = f_pe`, i.e. `n_c,O1 = ε₀·m_e·(2π·f)²/e²`
//!   (Stix, *Waves in Plasmas* (1992), Ch. 1).
//! - **ECRH-5 — X-mode (right-hand) cut-off.** `f = f_R = f_ce/2 +
//!   √((f_ce/2)² + f_pe²)`, inverted to
//!   `n_c,X1 = n_c,O1·(f − f_ce)/f`, which exists only for `f > f_ce`
//!   (Stix 1992, Ch. 1 — the low-field-side launch cut-off). Where
//!   `f ≤ f_ce` the X1 wave is evanescent and [`x1_cutoff_density_m3`]
//!   returns [`None`].
//!
//! # Inputs, outputs, units
//!
//! Scalar kernels take the wave frequency in GHz, fields in tesla, densities
//! in m⁻³, electron temperatures in keV. Beamline locators take caller
//! profiles tabulated at strictly increasing positions `s_m` in metres and
//! return crossing positions in metres as plain `f64` vectors (linear
//! interpolation between bracketing nodes; exact nodal touches are
//! reported). No files, no HDF5, no transport.
//!
//! # Layering
//!
//! This module lives in `nucleide-plasma-source` because ECRH access is a
//! property of the wave plus the caller's magnetic geometry — the same
//! caller-profile-in / plain-floats-out shape as the parametric source's
//! caller-supplied profiles — and it needs no dependency beyond the crate's
//! existing ones (in fact it needs none at all). It performs no sampling,
//! emits no cards, and knows nothing about fuel mixtures or toroidal
//! sectors.

use thiserror::Error;

/// SI elementary charge in coulombs (exact by definition).
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
/// Electron mass in kilograms.
pub const ELECTRON_MASS_KG: f64 = 9.109_383_701_5e-31;
/// Vacuum permittivity in farads per metre.
pub const VACUUM_PERMITTIVITY_F_PER_M: f64 = 8.854_187_812_8e-12;
/// Electron rest energy in keV.
pub const ELECTRON_REST_ENERGY_KEV: f64 = 510.998_95;

/// Cold electron gyrofrequency per tesla in GHz/T (ECRH-1 at `B = 1 T`).
///
/// Pinned by gate: ≈ 27.99249 GHz/T from the constants above.
pub const GYROFREQUENCY_GHZ_PER_T: f64 =
    ELEMENTARY_CHARGE_C / (2.0 * core::f64::consts::PI * ELECTRON_MASS_KG) / 1.0e9;

/// O1 cut-off coefficient: `n_c,O1 = O1_COEFFICIENT_M3_PER_GHZ2 · f_GHz²`
/// (ECRH-4 with `f` in GHz).
pub const O1_COEFFICIENT_M3_PER_GHZ2: f64 = VACUUM_PERMITTIVITY_F_PER_M
    * ELECTRON_MASS_KG
    * (2.0 * core::f64::consts::PI * 1.0e9)
    * (2.0 * core::f64::consts::PI * 1.0e9)
    / (ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C);

/// Errors raised by the ECRH accessibility kernel.
///
/// Every malformed input is loud; out-of-range physics (X1 evanescence) is
/// [`None`], never a silent substitution.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum EcrhError {
    /// A scalar input is `NaN` or infinite: `{0}`.
    #[error("ecrh: non-finite value in {0}")]
    NonFinite(&'static str),
    /// The wave frequency must be finite and strictly positive; got `{0}` GHz.
    #[error("ecrh: frequency must be > 0 GHz, got {0}")]
    InvalidFrequency(f64),
    /// The cyclotron harmonic must be ≥ 1; got `{0}`.
    #[error("ecrh: harmonic must be >= 1, got {0}")]
    InvalidHarmonic(u32),
    /// The electron temperature must be finite and non-negative; got `{0}` keV.
    #[error("ecrh: electron temperature must be >= 0 keV, got {0}")]
    NegativeTemperature(f64),
    /// A beamline profile is empty.
    #[error("ecrh: beamline profile is empty")]
    EmptyBeamline,
    /// Beamline slices disagree in length: `{0}`.
    #[error("ecrh: beamline length mismatch: {0}")]
    LengthMismatch(&'static str),
    /// Beamline positions must be strictly increasing in metres.
    #[error("ecrh: beamline positions must be strictly increasing")]
    NonIncreasingPositions,
    /// A beamline entry is `NaN` or infinite: `{0}`.
    #[error("ecrh: non-finite beamline entry in {0}")]
    NonFiniteBeamline(&'static str),
    /// An electron density is negative (`{0}` m⁻³).
    #[error("ecrh: electron density must be >= 0 m^-3, got {0}")]
    NegativeDensity(f64),
}

fn check_frequency(frequency_ghz: f64) -> Result<f64, EcrhError> {
    if !frequency_ghz.is_finite() {
        return Err(EcrhError::NonFinite("frequency_ghz"));
    }
    if frequency_ghz <= 0.0 {
        return Err(EcrhError::InvalidFrequency(frequency_ghz));
    }
    Ok(frequency_ghz)
}

fn check_harmonic(harmonic: u32) -> Result<f64, EcrhError> {
    if harmonic < 1 {
        return Err(EcrhError::InvalidHarmonic(harmonic));
    }
    Ok(f64::from(harmonic))
}

fn check_field(b_t: f64, name: &'static str) -> Result<f64, EcrhError> {
    if !b_t.is_finite() {
        return Err(EcrhError::NonFinite(name));
    }
    Ok(b_t)
}

/// Cold electron gyrofrequency in GHz at field `b_t` in tesla (ECRH-1).
///
/// # Errors
///
/// Returns [`EcrhError::NonFinite`] for a non-finite field.
pub fn gyrofrequency_ghz(b_t: f64) -> Result<f64, EcrhError> {
    Ok(check_field(b_t, "b_t")? * GYROFREQUENCY_GHZ_PER_T)
}

/// Cold resonant field in tesla for wave frequency `frequency_ghz` at
/// harmonic `harmonic` (ECRH-2).
///
/// # Errors
///
/// Returns [`EcrhError::InvalidFrequency`] for a non-positive frequency,
/// [`EcrhError::InvalidHarmonic`] for harmonic 0, [`EcrhError::NonFinite`]
/// for non-finite inputs.
pub fn cold_resonant_field_t(frequency_ghz: f64, harmonic: u32) -> Result<f64, EcrhError> {
    let f = check_frequency(frequency_ghz)?;
    let n = check_harmonic(harmonic)?;
    Ok(f / (n * GYROFREQUENCY_GHZ_PER_T))
}

/// Thermal Lorentz factor `γ = 1 + T_e/m_e·c²` for electron temperature
/// `electron_temperature_kev` in keV (ECRH-3).
///
/// # Errors
///
/// Returns [`EcrhError::NegativeTemperature`] for a negative temperature,
/// [`EcrhError::NonFinite`] for a non-finite one.
pub fn relativistic_gamma(electron_temperature_kev: f64) -> Result<f64, EcrhError> {
    if !electron_temperature_kev.is_finite() {
        return Err(EcrhError::NonFinite("electron_temperature_kev"));
    }
    if electron_temperature_kev < 0.0 {
        return Err(EcrhError::NegativeTemperature(electron_temperature_kev));
    }
    Ok(1.0 + electron_temperature_kev / ELECTRON_REST_ENERGY_KEV)
}

/// Relativistically shifted resonant field in tesla (ECRH-3):
/// `B_res(T_e) = γ·B_res(cold)`.
///
/// # Errors
///
/// Same as [`cold_resonant_field_t`] plus [`EcrhError::NegativeTemperature`]
/// for a negative electron temperature.
pub fn relativistic_resonant_field_t(
    frequency_ghz: f64,
    harmonic: u32,
    electron_temperature_kev: f64,
) -> Result<f64, EcrhError> {
    Ok(cold_resonant_field_t(frequency_ghz, harmonic)?
        * relativistic_gamma(electron_temperature_kev)?)
}

/// O-mode (ordinary) cut-off density in m⁻³ at `frequency_ghz` in GHz
/// (ECRH-4).
///
/// # Errors
///
/// Returns [`EcrhError::InvalidFrequency`] for a non-positive frequency,
/// [`EcrhError::NonFinite`] for a non-finite one.
pub fn o1_cutoff_density_m3(frequency_ghz: f64) -> Result<f64, EcrhError> {
    let f = check_frequency(frequency_ghz)?;
    Ok(O1_COEFFICIENT_M3_PER_GHZ2 * f * f)
}

/// X-mode (right-hand) cut-off density in m⁻³ at `frequency_ghz` in GHz and
/// local field `b_t` in tesla (ECRH-5).
///
/// Returns [`None`] where `f ≤ f_ce(B)` — the X1 wave is evanescent there
/// and no cut-off exists.
///
/// # Errors
///
/// Returns [`EcrhError::InvalidFrequency`] for a non-positive frequency,
/// [`EcrhError::NonFinite`] for non-finite inputs.
pub fn x1_cutoff_density_m3(frequency_ghz: f64, b_t: f64) -> Result<Option<f64>, EcrhError> {
    let f = check_frequency(frequency_ghz)?;
    let b = check_field(b_t, "b_t")?;
    let f_ce = b * GYROFREQUENCY_GHZ_PER_T;
    if f <= f_ce {
        return Ok(None);
    }
    Ok(Some(O1_COEFFICIENT_M3_PER_GHZ2 * f * (f - f_ce)))
}

fn check_beamline(s_m: &[f64], name_s: &'static str) -> Result<(), EcrhError> {
    if s_m.is_empty() {
        return Err(EcrhError::EmptyBeamline);
    }
    if !s_m.iter().all(|v| v.is_finite()) {
        return Err(EcrhError::NonFiniteBeamline(name_s));
    }
    if !s_m.windows(2).all(|w| w[1] > w[0]) {
        return Err(EcrhError::NonIncreasingPositions);
    }
    Ok(())
}

fn check_profile(values: &[f64], expected: usize, name: &'static str) -> Result<(), EcrhError> {
    if values.len() != expected {
        return Err(EcrhError::LengthMismatch(name));
    }
    if !values.iter().all(|v| v.is_finite()) {
        return Err(EcrhError::NonFiniteBeamline(name));
    }
    Ok(())
}

/// Linear-interpolation zero crossings of `a − b` along `s_m`.
///
/// Exact nodal touches are reported once (at the node); strict sign changes
/// interpolate. `s_m` must be strictly increasing (checked by callers).
fn crossings(s_m: &[f64], a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for w in 0..s_m.len().saturating_sub(1) {
        let d0 = a[w] - b[w];
        let d1 = a[w + 1] - b[w + 1];
        if d0 == 0.0 {
            out.push(s_m[w]);
        } else if d1 != 0.0 && d0.signum() != d1.signum() {
            let t = d0 / (d0 - d1);
            out.push(s_m[w] + t * (s_m[w + 1] - s_m[w]));
        }
    }
    let (&s_last, &a_last, &b_last) = (
        s_m.last().expect("non-empty"),
        a.last().expect("non-empty"),
        b.last().expect("non-empty"),
    );
    if a_last == b_last {
        out.push(s_last);
    }
    out
}

/// Cold/relativistic resonance positions in metres along the beamline: where
/// the caller field `b_t` crosses `resonant_field_t` (ECRH-2/ECRH-3 located
/// on caller geometry).
///
/// # Errors
///
/// Loud on an empty or ragged beamline, non-increasing positions,
/// non-finite entries, or a non-finite resonant field.
pub fn resonance_positions_m(
    s_m: &[f64],
    b_t: &[f64],
    resonant_field_t: f64,
) -> Result<Vec<f64>, EcrhError> {
    check_beamline(s_m, "s_m")?;
    check_profile(b_t, s_m.len(), "b_t vs s_m")?;
    let target = check_field(resonant_field_t, "resonant_field_t")?;
    Ok(crossings(s_m, b_t, &vec![target; s_m.len()]))
}

/// O1 cut-off positions in metres along the beamline: where the caller
/// density `ne_m3` crosses `cutoff_density_m3` (ECRH-4 located on caller
/// geometry).
///
/// # Errors
///
/// Loud on an empty or ragged beamline, non-increasing positions,
/// non-finite entries, negative densities, or a non-finite cut-off density.
pub fn o1_cutoff_positions_m(
    s_m: &[f64],
    ne_m3: &[f64],
    cutoff_density_m3: f64,
) -> Result<Vec<f64>, EcrhError> {
    check_beamline(s_m, "s_m")?;
    check_profile(ne_m3, s_m.len(), "ne_m3 vs s_m")?;
    if !cutoff_density_m3.is_finite() {
        return Err(EcrhError::NonFinite("cutoff_density_m3"));
    }
    if let Some(&bad) = ne_m3.iter().find(|v| **v < 0.0) {
        return Err(EcrhError::NegativeDensity(bad));
    }
    Ok(crossings(s_m, ne_m3, &vec![cutoff_density_m3; s_m.len()]))
}

/// X1 cut-off positions in metres along the beamline: where the caller
/// density `ne_m3` crosses the local right-hand cut-off density
/// `n_c,X1(B(s))` (ECRH-5 located on caller geometry).
///
/// Segments with an evanescent endpoint (`f ≤ f_ce(B)`, where no cut-off
/// exists) report no crossing — the wave never meets a cut-off there.
///
/// # Errors
///
/// Loud on an empty or ragged beamline, non-increasing positions,
/// non-finite entries, negative densities, or a non-positive frequency.
pub fn x1_cutoff_positions_m(
    s_m: &[f64],
    b_t: &[f64],
    ne_m3: &[f64],
    frequency_ghz: f64,
) -> Result<Vec<f64>, EcrhError> {
    check_beamline(s_m, "s_m")?;
    check_profile(b_t, s_m.len(), "b_t vs s_m")?;
    check_profile(ne_m3, s_m.len(), "ne_m3 vs s_m")?;
    let f = check_frequency(frequency_ghz)?;
    if let Some(&bad) = ne_m3.iter().find(|v| **v < 0.0) {
        return Err(EcrhError::NegativeDensity(bad));
    }
    let local_cutoff: Vec<Option<f64>> = b_t
        .iter()
        .map(|b| {
            let f_ce = b * GYROFREQUENCY_GHZ_PER_T;
            if f <= f_ce {
                None
            } else {
                Some(O1_COEFFICIENT_M3_PER_GHZ2 * f * (f - f_ce))
            }
        })
        .collect();
    let mut out = Vec::new();
    for w in 0..s_m.len().saturating_sub(1) {
        let (Some(c0), Some(c1)) = (local_cutoff[w], local_cutoff[w + 1]) else {
            continue;
        };
        let d0 = ne_m3[w] - c0;
        let d1 = ne_m3[w + 1] - c1;
        if d0 == 0.0 {
            out.push(s_m[w]);
        } else if d1 != 0.0 && d0.signum() != d1.signum() {
            let t = d0 / (d0 - d1);
            out.push(s_m[w] + t * (s_m[w + 1] - s_m[w]));
        }
    }
    if let (Some(&s_last), Some(&ne_last), Some(&Some(c_last))) =
        (s_m.last(), ne_m3.last(), local_cutoff.last())
    {
        if ne_last == c_last {
            out.push(s_last);
        }
    }
    Ok(out)
}

/// Whole-beamline ECRH accessibility report: scalar resonance/cut-off
/// quantities plus their located positions in metres.
///
/// This is the per-port penalty input shape the blanket bookkeeping
/// consumes: plain floats, no files.
#[derive(Debug, Clone, PartialEq)]
pub struct EcrhAccessibility {
    /// Cold resonant field in tesla (ECRH-2).
    pub resonant_field_t: f64,
    /// Relativistically shifted resonant field in tesla (ECRH-3), or `None`
    /// when no electron temperature was supplied.
    pub relativistic_field_t: Option<f64>,
    /// O1 cut-off density in m⁻³ (ECRH-4).
    pub o1_cutoff_density_m3: f64,
    /// Resonance positions in metres (field crossings of the located
    /// resonant field: relativistic when `T_e` is given, cold otherwise).
    pub resonance_m: Vec<f64>,
    /// O1 cut-off positions in metres.
    pub o1_cutoff_m: Vec<f64>,
    /// X1 cut-off positions in metres (empty where the wave is evanescent).
    pub x1_cutoff_m: Vec<f64>,
}

/// Locate cold resonance (or the relativistic shift when
/// `electron_temperature_kev` is `Some`) plus O1/X1 cut-offs along the
/// caller beamline `(s_m, b_t, ne_m3)` at `frequency_ghz`/`harmonic`.
///
/// # Errors
///
/// Any malformed scalar or beamline input is a loud [`EcrhError`]; see the
/// individual kernels for the exact conditions.
pub fn accessibility(
    s_m: &[f64],
    b_t: &[f64],
    ne_m3: &[f64],
    frequency_ghz: f64,
    harmonic: u32,
    electron_temperature_kev: Option<f64>,
) -> Result<EcrhAccessibility, EcrhError> {
    let cold = cold_resonant_field_t(frequency_ghz, harmonic)?;
    let rel = electron_temperature_kev
        .map(relativistic_gamma)
        .transpose()?
        .map(|g| cold * g);
    let located = rel.unwrap_or(cold);
    let o1 = o1_cutoff_density_m3(frequency_ghz)?;
    Ok(EcrhAccessibility {
        resonant_field_t: cold,
        relativistic_field_t: rel,
        o1_cutoff_density_m3: o1,
        resonance_m: resonance_positions_m(s_m, b_t, located)?,
        o1_cutoff_m: o1_cutoff_positions_m(s_m, ne_m3, o1)?,
        x1_cutoff_m: x1_cutoff_positions_m(s_m, b_t, ne_m3, frequency_ghz)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const REL_TOL: f64 = 1e-9;

    fn approx(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= REL_TOL * expected.abs().max(1e-300),
            "actual {actual} != expected {expected}"
        );
    }

    #[test]
    fn constants_pin_ghz_per_t() {
        approx(GYROFREQUENCY_GHZ_PER_T, 27.992_489_872_333_04);
        approx(O1_COEFFICIENT_M3_PER_GHZ2, 1.240_442_606_115_044_2e16);
    }

    #[test]
    fn gyrofrequency_hand_vector() {
        approx(gyrofrequency_ghz(5.0).unwrap(), 139.962_449_361_665_2);
        approx(gyrofrequency_ghz(3.0).unwrap(), 83.977_469_616_999_12);
    }

    #[test]
    fn cold_resonance_hand_vector() {
        // 170 GHz fundamental: ITER-class cold resonance.
        approx(
            cold_resonant_field_t(170.0, 1).unwrap(),
            6.073_057_479_892_957,
        );
        approx(
            cold_resonant_field_t(170.0, 2).unwrap(),
            6.073_057_479_892_957 / 2.0,
        );
    }

    #[test]
    fn relativistic_shift_hand_vector() {
        approx(relativistic_gamma(0.0).unwrap(), 1.0);
        approx(relativistic_gamma(20.0).unwrap(), 1.039_139_023_671_183_7);
        approx(
            relativistic_resonant_field_t(170.0, 1, 20.0).unwrap(),
            6.310_751_020_354_947,
        );
    }

    #[test]
    fn cutoff_hand_vectors() {
        approx(
            o1_cutoff_density_m3(170.0).unwrap(),
            3.584_879_131_672_478e20,
        );
        approx(
            o1_cutoff_density_m3(140.0).unwrap(),
            2.431_267_507_985_486_6e20,
        );
        approx(
            x1_cutoff_density_m3(170.0, 5.0).unwrap().unwrap(),
            6.334_175_791_171_868e19,
        );
        // 100 GHz at 5 T: f < f_ce = 139.96 GHz, evanescent, no cut-off.
        assert_eq!(x1_cutoff_density_m3(100.0, 5.0).unwrap(), None);
    }

    #[test]
    fn resonance_positions_hand_vector() {
        // B falls 7 -> 5 T over s = 0..2 m; 170 GHz cold resonance at
        // 6.073 T sits 0.9269 m along.
        let s = [0.0, 1.0, 2.0];
        let b = [7.0, 6.0, 5.0];
        let pos = resonance_positions_m(&s, &b, 6.073_057_479_892_957).unwrap();
        assert_eq!(pos.len(), 1);
        approx(pos[0], 0.926_942_520_107_043_1);
        // No crossing when the profile never reaches the target.
        assert!(resonance_positions_m(&s, &b, 8.0).unwrap().is_empty());
        // Exact nodal touch is reported.
        let touch = resonance_positions_m(&s, &b, 6.0).unwrap();
        assert_eq!(touch, vec![1.0]);
    }

    #[test]
    fn o1_cutoff_positions_hand_vector() {
        // ne rises 1e20 -> 3e20 -> 5e20; O1 cut-off at 140 GHz (2.431e20)
        // sits at s = 0.7156 m in the first segment.
        let s = [0.0, 1.0, 2.0];
        let ne = [1.0e20, 3.0e20, 5.0e20];
        let pos = o1_cutoff_positions_m(&s, &ne, 2.431_267_507_985_486_6e20).unwrap();
        assert_eq!(pos.len(), 1);
        approx(pos[0], 0.715_633_753_992_743_3);
    }

    #[test]
    fn x1_cutoff_positions_hand_vector() {
        // Uniform 5 T, ne 1e19 -> 1e20 against n_c,X1(170 GHz, 5 T).
        let s = [0.0, 1.0];
        let b = [5.0, 5.0];
        let ne = [1.0e19, 1.0e20];
        let pos = x1_cutoff_positions_m(&s, &b, &ne, 170.0).unwrap();
        assert_eq!(pos.len(), 1);
        approx(pos[0], 0.592_686_199_019_096_5);
        // Evanescent everywhere: 100 GHz at 5 T reports no cut-off.
        assert!(x1_cutoff_positions_m(&s, &b, &ne, 100.0)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn accessibility_report_is_consistent() {
        let s = [0.0, 1.0, 2.0];
        let b = [7.0, 6.0, 5.0];
        let ne = [1.0e20, 3.0e20, 5.0e20];
        let rep = accessibility(&s, &b, &ne, 170.0, 1, Some(20.0)).unwrap();
        approx(rep.resonant_field_t, 6.073_057_479_892_957);
        approx(rep.relativistic_field_t.unwrap(), 6.310_751_020_354_947);
        approx(rep.o1_cutoff_density_m3, 3.584_879_131_672_478e20);
        // Relativistic resonance at 6.311 T sits 0.689 m along.
        assert_eq!(rep.resonance_m.len(), 1);
        approx(rep.resonance_m[0], 0.689_248_979_645_053);
        assert!(rep.relativistic_field_t.is_some());
        let cold_rep = accessibility(&s, &b, &ne, 170.0, 1, None).unwrap();
        assert_eq!(cold_rep.relativistic_field_t, None);
        approx(cold_rep.resonance_m[0], 0.926_942_520_107_043_1);
    }

    #[test]
    fn malformed_inputs_are_loud() {
        assert!(matches!(
            cold_resonant_field_t(0.0, 1),
            Err(EcrhError::InvalidFrequency(_))
        ));
        assert!(matches!(
            cold_resonant_field_t(f64::NAN, 1),
            Err(EcrhError::NonFinite(_))
        ));
        assert!(matches!(
            cold_resonant_field_t(170.0, 0),
            Err(EcrhError::InvalidHarmonic(0))
        ));
        assert!(matches!(
            relativistic_gamma(-1.0),
            Err(EcrhError::NegativeTemperature(_))
        ));
        assert!(matches!(
            gyrofrequency_ghz(f64::INFINITY),
            Err(EcrhError::NonFinite(_))
        ));
        let s = [0.0, 1.0];
        let b = [5.0, 5.0];
        let ne = [1.0e19, 1.0e20];
        assert!(matches!(
            resonance_positions_m(&[], &[], 5.0),
            Err(EcrhError::EmptyBeamline)
        ));
        assert!(matches!(
            resonance_positions_m(&s, &[5.0], 5.0),
            Err(EcrhError::LengthMismatch(_))
        ));
        assert!(matches!(
            resonance_positions_m(&[1.0, 0.0], &b, 5.0),
            Err(EcrhError::NonIncreasingPositions)
        ));
        assert!(matches!(
            resonance_positions_m(&[0.0, f64::NAN], &b, 5.0),
            Err(EcrhError::NonFiniteBeamline(_))
        ));
        assert!(matches!(
            o1_cutoff_positions_m(&s, &[1.0e19, -1.0], 1.0e20),
            Err(EcrhError::NegativeDensity(_))
        ));
        assert!(matches!(
            x1_cutoff_positions_m(&s, &b, &ne, -170.0),
            Err(EcrhError::InvalidFrequency(_))
        ));
        assert!(matches!(
            accessibility(&s, &b, &ne, 170.0, 1, Some(-5.0)),
            Err(EcrhError::NegativeTemperature(_))
        ));
    }
}
