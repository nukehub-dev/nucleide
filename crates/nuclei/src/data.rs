//! Static nuclear reference data: atomic masses (AME2020), natural
//! abundances, and radioactive half-lives.
//!
//! # Provenance
//!
//! - Atomic masses: the [AME2020 atomic mass evaluation][ame] (Huang et al.,
//!   *Chinese Physics C* **45**, 030002/030003, 2021; data courtesy the
//!   IAEA-supported AMDC), condensed into the compact [`crate::data`]
//!   tables below.
//! - Natural abundances: standard isotopic compositions from the
//!   ENDF/B-VIII.0 evaluation, expressed as fractions in 0..1.
//! - Half-lives: ENDF/B-VIII.0 decay evaluations (distributed via the IAEA
//!   and BNL/NNDC), expressed in seconds.
//!
//! All tables are vendored as tab-separated text under `src/data/` and
//! embedded with `include_str!`; no runtime dependencies beyond `std`.
//! Parsing happens lazily on first lookup into fixed static maps.
//!
//! Table contents:
//!
//! - `data/ame2020.tsv`: one row per ground-state nuclide (`nucid`,
//!   `mass_u`, `uncertainty_u`), covering all 3 557 nuclides with `Z >= 1`.
//!   The free-neutron row of the source file is dropped, and excited levels
//!   are absent because `mass.mas20` lists one entry per (Z, A).
//! - `data/natural_abundance.tsv`: `GNDS name` → `fraction`, including the
//!   lone naturally occurring isomer Ta180_m1.
//! - `data/half_life.tsv`: `GNDS name` → `half_life_seconds`, for every
//!   radionuclide in the evaluation (stable nuclides are simply absent).
//!
//! Masses are keyed by the canonical [`NuclideId`] nucid layout; metastable
//! states have no mass entries (ground state and isomers share an atomic
//! mass — query the state-0 nucid).
//!
//! [ame]: https://doi.org/10.1088/1674-1137/abddaf
//!
//! Ownership note: this module + its data files only — do NOT edit `lib.rs`.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::NuclideId;

const AME2020_TSV: &str = include_str!("data/ame2020.tsv");
const NATURAL_ABUNDANCE_TSV: &str = include_str!("data/natural_abundance.tsv");
const HALF_LIFE_TSV: &str = include_str!("data/half_life.tsv");

static MASSES: OnceLock<BTreeMap<u32, f64>> = OnceLock::new();
static ABUNDANCES: OnceLock<BTreeMap<u32, f64>> = OnceLock::new();
static HALF_LIVES: OnceLock<BTreeMap<u32, f64>> = OnceLock::new();

/// MeV per unified atomic mass unit `c²` (2022 CODATA consistent with
/// AME2020 usage).
pub const MEV_PER_U: f64 = 931.494_102_42;

/// Free-neutron atomic mass in u (AME2020: 1.00866491595 u).
pub const NEUTRON_MASS_U: f64 = 1.008_664_915_95;

/// Helium-4 atomic mass in u (AME2020), for alpha-decay Q-values.
pub const HELIUM4_MASS_U: f64 = 4.002_603_254_13;

/// The free neutron's mass, in u.
///
/// See [`NEUTRON_MASS_U`].
pub const fn neutron_mass_u() -> f64 {
    NEUTRON_MASS_U
}

/// Parse `nucid \t mass_u [\t uncertainty_u]` rows, skipping comment lines.
fn parse_masses(tsv: &str) -> BTreeMap<u32, f64> {
    tsv.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .filter_map(|line| {
            let mut cols = line.split('\t');
            let nucid = cols.next()?.parse().ok()?;
            let mass = cols.next()?.parse().ok()?;
            Some((nucid, mass))
        })
        .collect()
}

/// Parse `GNDS name \t fraction` rows into nucid-keyed fractions.
fn parse_abundances(tsv: &str) -> BTreeMap<u32, f64> {
    tsv.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let mut cols = line.split('\t');
            let name = cols.next().expect("non-empty line");
            let fraction: f64 = cols
                .next()
                .expect("fraction column")
                .parse()
                .expect("valid fraction");
            let nucid = NuclideId::from_name(name).unwrap_or_else(|e| {
                panic!("invalid nuclide name `{name}` in abundance table: {e}")
            });
            (nucid.nucid(), fraction)
        })
        .collect()
}

/// Parse `GNDS name \t half_life_seconds` rows into nucid-keyed seconds.
fn parse_half_lives(tsv: &str) -> BTreeMap<u32, f64> {
    tsv.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let mut cols = line.split('\t');
            let name = cols.next().expect("non-empty line");
            let seconds: f64 = cols
                .next()
                .expect("half-life column")
                .parse()
                .expect("valid half-life");
            let nucid = NuclideId::from_name(name).unwrap_or_else(|e| {
                panic!("invalid nuclide name `{name}` in half-life table: {e}")
            });
            (nucid.nucid(), seconds)
        })
        .collect()
}

fn masses() -> &'static BTreeMap<u32, f64> {
    MASSES.get_or_init(|| parse_masses(AME2020_TSV))
}

fn abundances() -> &'static BTreeMap<u32, f64> {
    ABUNDANCES.get_or_init(|| parse_abundances(NATURAL_ABUNDANCE_TSV))
}

fn half_lives() -> &'static BTreeMap<u32, f64> {
    HALF_LIVES.get_or_init(|| parse_half_lives(HALF_LIFE_TSV))
}

/// The full AME2020 ground-state mass table, keyed by nucid (in u).
pub fn mass_table() -> &'static BTreeMap<u32, f64> {
    masses()
}

/// The full natural-abundance table, keyed by nucid (fraction in 0..1).
pub fn abundance_table() -> &'static BTreeMap<u32, f64> {
    abundances()
}

/// The full half-life table in seconds, keyed by nucid.
pub fn half_life_table() -> &'static BTreeMap<u32, f64> {
    half_lives()
}

/// Atomic mass of the nuclide with the given ground-state `nucid`, in u.
///
/// Returns `None` for unknown nuclides, non-nuclide ids (`Z = 0`), and
/// metastable state ids (query the ground state instead).
pub fn atomic_mass(nucid: u32) -> Option<f64> {
    masses().get(&nucid).copied()
}

/// Atomic mass of a named nuclide (see [`NuclideId::from_name`]), in u.
pub fn atomic_mass_by_name(name: &str) -> Option<f64> {
    atomic_mass(NuclideId::from_name(name).ok()?.nucid())
}

/// Natural abundance of the nuclide with the given `nucid`, as a fraction
/// in 0..1.
///
/// Only naturally occurring nuclides have entries, including the isomer
/// Ta180_m1; everything else returns `None`.
pub fn natural_abundance(nucid: u32) -> Option<f64> {
    abundances().get(&nucid).copied()
}

/// Natural abundance of a named nuclide, as a fraction in 0..1.
pub fn natural_abundance_by_name(name: &str) -> Option<f64> {
    natural_abundance(NuclideId::from_name(name).ok()?.nucid())
}

/// Radioactive half-life of the nuclide with the given `nucid`, in seconds.
///
/// Values come from the ENDF/B-VIII.0-derived evaluation; stable
/// nuclides and unknown ids have no entry and return `None`.
pub fn half_life(nucid: u32) -> Option<f64> {
    half_lives().get(&nucid).copied()
}

/// Half-life in seconds of a named nuclide (see [`NuclideId::from_name`]).
pub fn half_life_by_name(name: &str) -> Option<f64> {
    half_life(NuclideId::from_name(name).ok()?.nucid())
}

/// Decay constant λ = ln(2) / t½ of the given `nucid`, in inverse seconds.
///
/// Computed from [`half_life`]; `None` wherever the half-life is unknown.
pub fn decay_constant(nucid: u32) -> Option<f64> {
    half_life(nucid).map(|t_half| std::f64::consts::LN_2 / t_half)
}

/// Decay constant λ = ln(2) / t½ (inverse seconds) of a named nuclide.
pub fn decay_constant_by_name(name: &str) -> Option<f64> {
    decay_constant(NuclideId::from_name(name).ok()?.nucid())
}

/// Q-value of neutron radiative capture X(n,γ)X', in MeV.
///
/// From atomic masses (AME2020, in u):
///
/// ```text
/// Q = [m(X + A) + m_n − m(X' = Z, A+1)] · c²   (MeV/u × 931.49410242)
/// ```
///
/// using the AME2020 free-neutron mass ([`NEUTRON_MASS_U`]). Electron
/// binding differences are neglected (standard for capture Q-values).
/// Returns `None` when either ground-state mass is missing, or for
/// metastable-state or non-nuclide inputs.
pub fn q_value_neutron_capture(nucid: u32) -> Option<f64> {
    let id = NuclideId::from_nucid(nucid);
    if id.state() != 0 || id.z() == 0 {
        return None;
    }
    let product = (id.z() * 1000 + id.a() + 1) * 10_000;
    let q_u = atomic_mass(nucid)? + NEUTRON_MASS_U - atomic_mass(product)?;
    Some(q_u * MEV_PER_U)
}

/// Q-value of neutron capture on a named nuclide, in MeV.
pub fn q_value_neutron_capture_by_name(name: &str) -> Option<f64> {
    q_value_neutron_capture(NuclideId::from_name(name).ok()?.nucid())
}

/// Q-value of alpha decay X → Y + He4, in MeV.
///
/// From atomic masses (AME2020, in u):
///
/// ```text
/// Q = [m(Z,A) − m(Z−2,A−4) − m(He4)] · c²   (MeV/u × 931.49410242)
/// ```
///
/// with [`HELIUM4_MASS_U`] = 4.00260325413 u. Using neutral atomic masses,
/// the two released electrons cancel exactly to first order. Returns `None`
/// when the daughter mass is absent from the table, or for metastable or
/// non-nuclide inputs.
pub fn q_value_alpha(nucid: u32) -> Option<f64> {
    let id = NuclideId::from_nucid(nucid);
    if id.state() != 0 || id.z() <= 2 || id.a() <= 4 {
        return None;
    }
    let daughter = ((id.z() - 2) * 1000 + (id.a() - 4)) * 10_000;
    let q_u = atomic_mass(nucid)? - atomic_mass(daughter)? - HELIUM4_MASS_U;
    Some(q_u * MEV_PER_U)
}

/// Q-value of alpha decay on a named nuclide, in MeV.
pub fn q_value_alpha_by_name(name: &str) -> Option<f64> {
    q_value_alpha(NuclideId::from_name(name).ok()?.nucid())
}

/// Zero-sized façade over this module's lookups.
///
/// Standalone today; intended to satisfy the material crate's future
/// `MassProvider` trait once that integration lands (kept free of cross-crate
/// coupling for now).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AmeMasses;

impl AmeMasses {
    /// See [`atomic_mass`].
    pub fn atomic_mass(&self, nucid: u32) -> Option<f64> {
        atomic_mass(nucid)
    }

    /// See [`atomic_mass_by_name`].
    pub fn atomic_mass_by_name(&self, name: &str) -> Option<f64> {
        atomic_mass_by_name(name)
    }

    /// See [`natural_abundance`].
    pub fn natural_abundance(&self, nucid: u32) -> Option<f64> {
        natural_abundance(nucid)
    }

    /// See [`natural_abundance_by_name`].
    pub fn natural_abundance_by_name(&self, name: &str) -> Option<f64> {
        natural_abundance_by_name(name)
    }
}

/// Zero-sized façade over the decay-data lookups ([`half_life`],
/// [`decay_constant`]).
///
/// Mirrors [`AmeMasses`]; standalone today, kept free of cross-crate
/// coupling for future provider-trait integration.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DecayData;

impl DecayData {
    /// See [`half_life`].
    pub fn half_life(&self, nucid: u32) -> Option<f64> {
        half_life(nucid)
    }

    /// See [`half_life_by_name`].
    pub fn half_life_by_name(&self, name: &str) -> Option<f64> {
        half_life_by_name(name)
    }

    /// See [`decay_constant`].
    pub fn decay_constant(&self, nucid: u32) -> Option<f64> {
        decay_constant(nucid)
    }

    /// See [`decay_constant_by_name`].
    pub fn decay_constant_by_name(&self, name: &str) -> Option<f64> {
        decay_constant_by_name(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H1: u32 = 10_010_000;
    const O16: u32 = 80_160_000;
    const FE56: u32 = 260_560_000;
    const U235: u32 = 922_350_000;

    #[test]
    fn h1_exact_ame2020_value() {
        assert_eq!(atomic_mass(H1), Some(1.007825031_898));
        assert_eq!(atomic_mass_by_name("H1"), Some(1.007825031_898));
    }

    #[test]
    fn heavy_nuclide_spot_values() {
        assert_eq!(atomic_mass(U235), Some(235.043_928_117));
        assert_eq!(atomic_mass(FE56), Some(55.934_935_537));
        assert_eq!(atomic_mass(O16), Some(15.994_914_619_26));
    }

    #[test]
    fn non_nuclides_and_unknown_ids_return_none() {
        // Free neutron row (Z = 0) is excluded from the vendored table.
        assert_eq!(atomic_mass(10_000), None);
        assert_eq!(atomic_mass(999_999_999), None);
        // Ground states only: metastable ids have no mass entries.
        assert_eq!(atomic_mass(922_350_001), None);
        assert_eq!(atomic_mass_by_name("U235_m1"), None);
    }

    #[test]
    fn by_name_agrees_with_nucid_lookup() {
        for name in ["H1", "O16", "Fe56", "U235", "Og294"] {
            let nucid = NuclideId::from_name(name).unwrap().nucid();
            assert_eq!(atomic_mass_by_name(name), atomic_mass(nucid), "{name}");
        }
    }

    #[test]
    fn natural_abundance_spot_values() {
        assert_eq!(natural_abundance(U235), Some(0.007_204));
        assert_eq!(natural_abundance_by_name("U235"), Some(0.007_204));
        assert_eq!(natural_abundance(O16), Some(0.997_620_6));
        assert_eq!(natural_abundance_by_name("H1"), Some(0.999_844_26));
        // The lone naturally occurring isomer.
        assert_eq!(natural_abundance_by_name("Ta180_m1"), Some(0.000_120_1));
    }

    #[test]
    fn natural_abundance_unknown_returns_none() {
        assert_eq!(natural_abundance(999_999_999), None);
        assert_eq!(natural_abundance_by_name("C14"), None);
        assert_eq!(natural_abundance_by_name("Xx999"), None);
    }

    #[test]
    fn abundances_sum_to_one_per_element() {
        let mut totals = [0.0_f64; 119];
        for (nucid, frac) in abundance_table() {
            totals[NuclideId::from_nucid(*nucid).z() as usize] += frac;
        }
        for (z, total) in totals.iter().enumerate() {
            if *total > 0.0 {
                assert!(
                    (total - 1.0).abs() < 1e-6,
                    "Z={z} abundances sum to {total}"
                );
            }
        }
    }

    #[test]
    fn mass_sanity_sweep() {
        let table = mass_table();
        for (nucid, mass) in table {
            let id = NuclideId::from_nucid(*nucid);
            assert_eq!(id.state(), 0, "only ground states expected");
            let (lo, hi) = (0.9 * f64::from(id.a()), 1.2 * f64::from(id.a()));
            assert!(*mass > lo && *mass < hi, "{} mass {mass}", id.to_name());
            assert!(*mass > 0.0);
        }
    }

    #[test]
    fn vendored_row_counts_match_tables() {
        let mass_rows = AME2020_TSV
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .count();
        let abundance_rows = NATURAL_ABUNDANCE_TSV
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .count();
        assert_eq!(mass_table().len(), mass_rows);
        assert_eq!(mass_rows, 3557);
        assert_eq!(abundance_table().len(), abundance_rows);
        assert_eq!(abundance_rows, 289);
    }

    #[test]
    fn ame_masses_facade_delegates() {
        let provider = AmeMasses;
        assert_eq!(provider.atomic_mass(U235), atomic_mass(U235));
        assert_eq!(provider.atomic_mass_by_name("Fe56"), Some(55.934_935_537));
        assert_eq!(provider.natural_abundance(O16), Some(0.997_620_6));
        assert_eq!(provider.natural_abundance_by_name("Nope1"), None);
    }

    const U238: u32 = 922_380_000;
    const I135: u32 = 531_350_000;
    const CS137: u32 = 551_370_000;

    #[test]
    fn half_life_spot_values() {
        // U-238: 4.468e9 yr × 3.1557e7 s/yr ≈ 1.4100e17 s.
        let t_u238 = half_life(U238).unwrap();
        assert!((t_u238 - 1.409_99e17).abs() / t_u238 < 1e-9);
        // I-135: chain fixture value 2.3652e4 s.
        assert_eq!(half_life(I135), Some(23_652.0));
        // Cs-135: 7.25825e13 s (≈ 2.3 Myr).
        let t_cs135 = half_life_by_name("Cs135").unwrap();
        assert!((t_cs135 - 7.258_25e13).abs() / t_cs135 < 1e-12);
        // Cs-137: 30.08 yr ≈ 9.49e8 s.
        let t_cs137 = half_life(CS137).unwrap();
        assert!((t_cs137 / 3.155_695_2e7 - 30.08).abs() < 0.01, "{t_cs137}");
        // Isomers carry their own rows.
        assert_eq!(half_life_by_name("Am242_m1"), Some(4_449_622_000.0));
    }

    #[test]
    fn stable_and_unknown_nuclides_have_no_half_life() {
        // Stable nuclides are absent from the vendored table.
        assert_eq!(half_life(O16), None);
        assert_eq!(half_life_by_name("Fe56"), None);
        assert_eq!(decay_constant(H1), None);
        assert_eq!(half_life(999_999_999), None);
        assert_eq!(half_life_by_name("Xx999"), None);
    }

    #[test]
    fn decay_constant_is_ln2_over_half_life() {
        for nucid in [
            U238,
            I135,
            CS137,
            NuclideId::from_name("Te132").unwrap().nucid(),
        ] {
            let t_half = half_life(nucid).unwrap();
            let lambda = decay_constant(nucid).unwrap();
            let rel = (lambda * t_half - std::f64::consts::LN_2).abs() / std::f64::consts::LN_2;
            assert!(rel < 1e-12, "nucid {nucid}: rel err {rel}");
        }
        let lam = decay_constant_by_name("I135").unwrap();
        assert!((lam - std::f64::consts::LN_2 / 23_652.0).abs() < 1e-18);
    }

    #[test]
    fn shorter_half_life_gives_larger_decay_constant() {
        // Te132 (3.2 d) vs I135 (6.57 h) vs Xe135 (9.14 h) chain ordering.
        let te = NuclideId::from_name("Te132").unwrap().nucid();
        let xe = NuclideId::from_name("Xe135").unwrap().nucid();
        let (t_te, t_i, t_xe) = (
            half_life(te).unwrap(),
            half_life(I135).unwrap(),
            half_life(xe).unwrap(),
        );
        assert!(t_te > t_xe && t_i < t_xe);
        assert!(decay_constant(te).unwrap() < decay_constant(xe).unwrap());
        assert!(decay_constant(xe).unwrap() < decay_constant(I135).unwrap());
    }

    #[test]
    fn half_lives_are_positive_and_finite() {
        for (nucid, t_half) in half_life_table() {
            assert!(*t_half > 0.0 && t_half.is_finite(), "{nucid}: {t_half}");
            let id = NuclideId::from_nucid(*nucid);
            assert!(id.a() >= id.z(), "{}", id.to_name());
        }
    }

    #[test]
    fn vendored_half_life_row_count_matches_table() {
        let rows = HALF_LIFE_TSV
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .count();
        assert_eq!(half_life_table().len(), rows);
        assert_eq!(rows, 3561);
    }

    #[test]
    fn neutron_capture_q_value_anchors() {
        // H1(n,γ): textbook 2.224566 MeV.
        let q_h = q_value_neutron_capture(H1).unwrap();
        assert!((q_h - 2.224_566).abs() < 1e-3, "{q_h}");
        // U238(n,γ): ≈ 4.8 MeV.
        let q_u238 = q_value_neutron_capture(U238).unwrap();
        assert!((q_u238 - 4.806_382).abs() < 1e-3, "{q_u238}");
        assert!(q_u238 > 4.79 && q_u238 < 4.81);
        // O16(n,γ): 4.143 MeV.
        let q_o16 = q_value_neutron_capture(O16).unwrap();
        assert!((q_o16 - 4.143_080).abs() < 1e-3, "{q_o16}");
    }

    #[test]
    fn capture_q_value_matches_manual_formula() {
        let expected = (atomic_mass(U238).unwrap() + NEUTRON_MASS_U
            - atomic_mass(922_390_000).unwrap())
            * MEV_PER_U;
        let q = q_value_neutron_capture(U238).unwrap();
        assert!((q - expected).abs() < 1e-9);
        assert_eq!(q_value_neutron_capture_by_name("U238"), Some(q));
        assert_eq!(neutron_mass_u(), NEUTRON_MASS_U);
        assert_eq!(neutron_mass_u(), 1.008_664_915_95);
    }

    #[test]
    fn capture_q_value_missing_or_bad_targets_return_none() {
        // Product nuclide beyond the chart of nuclides (no A+1 mass row).
        let he10 = NuclideId::from_name("He10").unwrap().nucid();
        assert_eq!(atomic_mass(he10 + 10_000), None);
        assert_eq!(q_value_neutron_capture(he10), None);
        // Metastable states and non-nuclides are rejected outright.
        assert_eq!(q_value_neutron_capture(922_350_001), None);
        assert_eq!(q_value_neutron_capture(10_000), None);
        assert_eq!(q_value_neutron_capture_by_name("U235_m1"), None);
        assert_eq!(q_value_neutron_capture_by_name("Nope1"), None);
    }

    #[test]
    fn alpha_q_value_anchors() {
        // Literature: U238 α → 4.2698 MeV, Po210 → 5.4075 MeV,
        // Ra226 → 4.8706 MeV; atomic-mass formula reproduces all three.
        for (name, lit) in [
            ("U238", 4.269_858),
            ("Po210", 5.407_530),
            ("Ra226", 4.870_703),
        ] {
            let q = q_value_alpha_by_name(name).unwrap();
            assert!((q - lit).abs() < 1e-3, "{name}: {q} vs {lit}");
        }
        assert_eq!(q_value_alpha(U238), q_value_alpha_by_name("U238"));
    }

    #[test]
    fn alpha_q_value_rejects_light_and_metastable() {
        // Too light to alpha-decay within the table's Z/A domain.
        assert_eq!(q_value_alpha(H1), None);
        assert_eq!(q_value_alpha_by_name("He4"), None);
        assert_eq!(q_value_alpha(922_350_001), None);
        assert_eq!(q_value_alpha_by_name("Am242_m1"), None);
        // Endothermic "decay" still yields the (negative) Q-value:
        // O16 → C12 + He4 costs ≈ 7.162 MeV.
        let q_o16 = q_value_alpha_by_name("O16").unwrap();
        assert!((q_o16 + 7.162).abs() < 1e-3, "{q_o16}");
    }

    #[test]
    fn decay_data_facade_delegates() {
        let provider = DecayData;
        assert_eq!(provider.half_life(I135), Some(23_652.0));
        assert_eq!(provider.half_life_by_name("I135"), Some(23_652.0));
        assert_eq!(provider.decay_constant(U238), decay_constant(U238));
        assert_eq!(provider.decay_constant_by_name("Fe56"), None);
    }
}
