//! Unit-aware decay-inventory facade plus cumulative/progeny observables.
//!
//! [`DecayInventory`] is a thin Tier-1 UX layer over atom counts: it converts
//! between activity / mass / mole / atom-count units, decays an inventory with
//! a single CRAM solve (decay-only, empty reaction rates), and formats
//! human-readable half-lives. [`cumulative_decays`] reports time-integrated
//! decays per nuclide, and the `progeny` / `branching_fraction` /
//! `decay_mode` / `chain_edges` helpers expose chain-daughter topology.
//!
//! Key decisions:
//! - Decay constants resolve **chain half-life first, then the nuclei table**
//!   (mirroring [`crate::integrate::decay_constants`] for chain members, with
//!   a `nucleide_nuclei::data::half_life_by_name` fallback so chain-stable
//!   entries with tabulated half-lives still convert). Truly unknown nuclides
//!   are [`Error::UnknownNuclide`], never silent zero; a stable nuclide with
//!   zero activity is legitimately zero.
//! - Masses come from `nucleide_nuclei::data::atomic_mass_by_name` (u ==
//!   g/mol numerically). This crate must not depend on `material`.
//! - [`cumulative_decays`] is **diagonal**: nuclide `i` reports
//!   `N0[i] * (1 - exp(-λ_i t))`, ignoring in-growth from parents during the
//!   interval. Stable nuclides report `0.0`.
//! - A year is exactly 365.25 days (31557600 s).

use std::collections::BTreeMap;
use std::f64::consts::LN_2;
use std::str::FromStr;

use crate::chain::{Chain, Error};
use crate::matrix::{DepletionSystem, ReactionRates};

/// Avogadro constant [1/mol] (exact SI definition).
pub const AVOGADRO: f64 = 6.02214076e23;

/// Curie in becquerel (exact definition).
pub const CURIE_TO_BQ: f64 = 3.7e10;

/// Seconds per Julian year (365.25 d).
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

/// Seconds per day.
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Time units accepted by [`DecayInventory::decay`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeUnit {
    /// Seconds.
    Seconds,
    /// Minutes (60 s).
    Minutes,
    /// Hours (3600 s).
    Hours,
    /// Days (86400 s).
    Days,
    /// Julian years (365.25 d).
    Years,
}

impl TimeUnit {
    /// Length of one unit in seconds.
    pub fn as_seconds(&self) -> f64 {
        match self {
            TimeUnit::Seconds => 1.0,
            TimeUnit::Minutes => 60.0,
            TimeUnit::Hours => 3_600.0,
            TimeUnit::Days => SECONDS_PER_DAY,
            TimeUnit::Years => SECONDS_PER_YEAR,
        }
    }
}

/// Parse a time-unit string (case-insensitive, surrounding whitespace
/// ignored). Accepts `s`/`sec`/`second(s)`, `m`/`min`/`minute(s)`,
/// `h`/`hr`/`hour(s)`, `d`/`day(s)`, `y`/`yr`/`year(s)` and plurals.
pub fn time_unit_from_str(s: &str) -> Result<TimeUnit, Error> {
    match s.trim().to_ascii_lowercase().as_str() {
        "s" | "sec" | "secs" | "second" | "seconds" => Ok(TimeUnit::Seconds),
        "m" | "min" | "mins" | "minute" | "minutes" => Ok(TimeUnit::Minutes),
        "h" | "hr" | "hrs" | "hour" | "hours" => Ok(TimeUnit::Hours),
        "d" | "day" | "days" => Ok(TimeUnit::Days),
        "y" | "yr" | "yrs" | "year" | "years" => Ok(TimeUnit::Years),
        other => Err(Error::BadStructure(format!("unknown time unit `{other}`"))),
    }
}

/// Quantity kinds behind [`QuantityUnit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QuantityKind {
    Activity,
    Mass,
    Moles,
    Atoms,
}

/// Units accepted for inventory quantities.
///
/// Activity: `Bq`, `kBq`–`PBq`, `Ci`, `mCi`, `uCi`, `dpm`.
/// Mass: `g`, `kg`, `mg`, `ug`, `t`.
/// Amount: `mol`, `mmol`.
/// Counts: `atoms`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantityUnit {
    /// Becquerel.
    Bq,
    /// Kilobecquerel.
    KBq,
    /// Megabecquerel.
    MBq,
    /// Gigabecquerel.
    GBq,
    /// Terabecquerel.
    TBq,
    /// Petabecquerel.
    PBq,
    /// Curie (3.7e10 Bq).
    Ci,
    /// Millicurie.
    MCi,
    /// Microcurie.
    UCi,
    /// Decays per minute (1/60 Bq).
    Dpm,
    /// Gram.
    Gram,
    /// Kilogram.
    Kilogram,
    /// Milligram.
    Milligram,
    /// Microgram.
    Microgram,
    /// Metric tonne (1e6 g).
    Tonne,
    /// Mole.
    Mol,
    /// Millimole.
    Mmol,
    /// Plain atom counts.
    Atoms,
}

impl std::str::FromStr for QuantityUnit {
    type Err = Error;

    /// Parse via [`QuantityUnit::parse`]; also reachable as
    /// `QuantityUnit::from_str` with `std::str::FromStr` in scope.
    fn from_str(s: &str) -> Result<Self, Error> {
        Self::parse(s)
    }
}

impl QuantityUnit {
    /// Parse a quantity-unit string (case-insensitive, whitespace ignored).
    ///
    /// Accepts SI prefixes on becquerel (`kbq`…`pbq`, plus `becquerel(s)`
    /// and `kbecquerel(s)`-style long forms), `ci`/`curie(s)`,
    /// `mci`/`millicurie(s)`, `uci`/`µci`/`microcurie(s)`, `dpm`, mass
    /// forms (`g`/`gram(s)`, `kg`…, `t`/`tonne(s)`/`ton(s)`), amount forms
    /// (`mol`/`mole(s)`, `mmol`/`millimole(s)`), and count forms
    /// (`atoms`/`atom`/`counts`/`count`/`number(s)`).
    fn parse(s: &str) -> Result<Self, Error> {
        match s.trim().to_ascii_lowercase().as_str() {
            "bq" | "becquerel" | "becquerels" => Ok(QuantityUnit::Bq),
            "kbq" | "kbecquerel" | "kbecquerels" => Ok(QuantityUnit::KBq),
            "mbq" | "mbecquerel" | "mbecquerels" => Ok(QuantityUnit::MBq),
            "gbq" | "gbecquerel" | "gbecquerels" => Ok(QuantityUnit::GBq),
            "tbq" | "tbecquerel" | "tbecquerels" => Ok(QuantityUnit::TBq),
            "pbq" | "pbecquerel" | "pbecquerels" => Ok(QuantityUnit::PBq),
            "ci" | "curie" | "curies" => Ok(QuantityUnit::Ci),
            "mci" | "millicurie" | "millicuries" => Ok(QuantityUnit::MCi),
            "uci" | "µci" | "μci" | "microcurie" | "microcuries" => Ok(QuantityUnit::UCi),
            "dpm" => Ok(QuantityUnit::Dpm),
            "g" | "gram" | "grams" => Ok(QuantityUnit::Gram),
            "kg" | "kilogram" | "kilograms" => Ok(QuantityUnit::Kilogram),
            "mg" | "milligram" | "milligrams" => Ok(QuantityUnit::Milligram),
            "ug" | "µg" | "μg" | "microgram" | "micrograms" => Ok(QuantityUnit::Microgram),
            "t" | "tonne" | "tonnes" | "ton" | "tons" => Ok(QuantityUnit::Tonne),
            "mol" | "mole" | "moles" => Ok(QuantityUnit::Mol),
            "mmol" | "millimole" | "millimoles" => Ok(QuantityUnit::Mmol),
            "atoms" | "atom" | "counts" | "count" | "number" | "numbers" => Ok(QuantityUnit::Atoms),
            other => Err(Error::BadStructure(format!(
                "unknown quantity unit `{other}`"
            ))),
        }
    }

    fn kind(&self) -> QuantityKind {
        match self {
            QuantityUnit::Bq
            | QuantityUnit::KBq
            | QuantityUnit::MBq
            | QuantityUnit::GBq
            | QuantityUnit::TBq
            | QuantityUnit::PBq
            | QuantityUnit::Ci
            | QuantityUnit::MCi
            | QuantityUnit::UCi
            | QuantityUnit::Dpm => QuantityKind::Activity,
            QuantityUnit::Gram
            | QuantityUnit::Kilogram
            | QuantityUnit::Milligram
            | QuantityUnit::Microgram
            | QuantityUnit::Tonne => QuantityKind::Mass,
            QuantityUnit::Mol | QuantityUnit::Mmol => QuantityKind::Moles,
            QuantityUnit::Atoms => QuantityKind::Atoms,
        }
    }

    /// Multiplicative factor to the kind's base unit
    /// (Bq, g, mol, atoms).
    fn base_factor(self) -> f64 {
        match self {
            QuantityUnit::Bq => 1.0,
            QuantityUnit::KBq => 1.0e3,
            QuantityUnit::MBq => 1.0e6,
            QuantityUnit::GBq => 1.0e9,
            QuantityUnit::TBq => 1.0e12,
            QuantityUnit::PBq => 1.0e15,
            QuantityUnit::Ci => CURIE_TO_BQ,
            QuantityUnit::MCi => CURIE_TO_BQ * 1.0e-3,
            QuantityUnit::UCi => CURIE_TO_BQ * 1.0e-6,
            QuantityUnit::Dpm => 1.0 / 60.0,
            QuantityUnit::Gram => 1.0,
            QuantityUnit::Kilogram => 1.0e3,
            QuantityUnit::Milligram => 1.0e-3,
            QuantityUnit::Microgram => 1.0e-6,
            QuantityUnit::Tonne => 1.0e6,
            QuantityUnit::Mol => 1.0,
            QuantityUnit::Mmol => 1.0e-3,
            QuantityUnit::Atoms => 1.0,
        }
    }

    /// Canonical short label used by [`DecayInventory::to_csv`].
    fn label(&self) -> &'static str {
        match self {
            QuantityUnit::Bq => "Bq",
            QuantityUnit::KBq => "kBq",
            QuantityUnit::MBq => "MBq",
            QuantityUnit::GBq => "GBq",
            QuantityUnit::TBq => "TBq",
            QuantityUnit::PBq => "PBq",
            QuantityUnit::Ci => "Ci",
            QuantityUnit::MCi => "mCi",
            QuantityUnit::UCi => "uCi",
            QuantityUnit::Dpm => "dpm",
            QuantityUnit::Gram => "g",
            QuantityUnit::Kilogram => "kg",
            QuantityUnit::Milligram => "mg",
            QuantityUnit::Microgram => "ug",
            QuantityUnit::Tonne => "t",
            QuantityUnit::Mol => "mol",
            QuantityUnit::Mmol => "mmol",
            QuantityUnit::Atoms => "atoms",
        }
    }
}

/// Resolve the decay constant for `name` [1/s]: chain half-life first, then
/// the nuclei half-life table.
///
/// Chain members with a positive chain half-life use it directly. Chain
/// members without one (chain-stable) fall back to the nuclei table so
/// tabulated half-lives still convert; chain-stable names absent from the
/// table are legitimately `0.0`. Names outside the chain resolve through the
/// nuclei table, and names unknown to both are [`Error::UnknownNuclide`].
fn decay_constant_for(sys: &DepletionSystem, name: &str) -> Result<f64, Error> {
    if let Some(idx) = sys.chain.index_of(name) {
        let lam = sys.chain.nuclides[idx].decay_constant()?;
        if lam > 0.0 {
            return Ok(lam);
        }
        if let Some(t) = nucleide_nuclei::data::half_life_by_name(name) {
            if t > 0.0 && t.is_finite() {
                return Ok(LN_2 / t);
            }
        }
        return Ok(0.0);
    }
    match nucleide_nuclei::data::half_life_by_name(name) {
        Some(t) if t > 0.0 && t.is_finite() => Ok(LN_2 / t),
        // A parseable name with no tabulated half-life is treated as stable.
        Some(_) => Ok(0.0),
        None => Err(Error::UnknownNuclide {
            name: name.to_string(),
            context: "inventory",
        }),
    }
}

/// Resolve the molar mass for `name` [g/mol, numerically equal to u].
fn molar_mass_for(name: &str) -> Result<f64, Error> {
    nucleide_nuclei::data::atomic_mass_by_name(name).ok_or_else(|| Error::UnknownNuclide {
        name: name.to_string(),
        context: "mass table",
    })
}

fn require_finite_nonnegative(name: &str, value: f64) -> Result<(), Error> {
    if !value.is_finite() || value < 0.0 {
        return Err(Error::BadStructure(format!(
            "invalid quantity for `{name}`: {value}"
        )));
    }
    Ok(())
}

/// Unit-aware decay inventory over atom counts.
///
/// The map holds atom counts keyed by nuclide name; all other units are
/// converted on entry ([`from_units`](Self::from_units)) and on exit
/// ([`activities`](Self::activities), [`masses`](Self::masses),
/// [`moles`](Self::moles)).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DecayInventory {
    /// Atom counts by nuclide name.
    pub atoms: BTreeMap<String, f64>,
}

impl DecayInventory {
    /// Empty inventory.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build directly from atom counts (must be finite and non-negative).
    pub fn from_atoms(map: &BTreeMap<String, f64>) -> Result<Self, Error> {
        for (name, value) in map {
            require_finite_nonnegative(name, *value)?;
        }
        Ok(Self { atoms: map.clone() })
    }

    /// Build from quantities in `unit`, converting to atom counts.
    ///
    /// `sys` supplies chain membership and chain-first decay constants (see
    /// [`decay_constant_for`]); masses and fallback half-lives come from the
    /// nuclei tables. Unknown nuclides and missing masses are errors, never
    /// silent zero — except a stable nuclide at zero activity, which is
    /// legitimately zero atoms. A non-zero activity for a stable nuclide is
    /// an error.
    pub fn from_units(
        map: &BTreeMap<String, f64>,
        unit: QuantityUnit,
        sys: &DepletionSystem,
    ) -> Result<Self, Error> {
        let mut atoms = BTreeMap::new();
        for (name, value) in map {
            require_finite_nonnegative(name, *value)?;
            let n = match unit.kind() {
                QuantityKind::Atoms => *value,
                QuantityKind::Activity => {
                    let a_bq = value * unit.base_factor();
                    let lam = decay_constant_for(sys, name)?;
                    if lam == 0.0 {
                        if a_bq != 0.0 {
                            return Err(Error::BadStructure(format!(
                                "non-zero activity for stable nuclide `{name}`"
                            )));
                        }
                        0.0
                    } else {
                        a_bq / lam
                    }
                }
                QuantityKind::Mass => {
                    let grams = value * unit.base_factor();
                    let molar = molar_mass_for(name)?;
                    grams / molar * AVOGADRO
                }
                QuantityKind::Moles => value * unit.base_factor() * AVOGADRO,
            };
            atoms.insert(name.clone(), n);
        }
        Ok(Self { atoms })
    }

    /// Decay the inventory by `dt` in `unit` with an explicit [`Method`].
    ///
    /// The solve is decay-only: the system is rebuilt from `sys.chain` with
    /// empty reaction rates, so [`Method::Bateman`]/[`Method::BatemanHp`]
    /// stay on the fast path here. Result covers all chain nuclides (missing
    /// inputs start at zero); inventory names absent from the chain are an
    /// error.
    pub fn decay_with_method(
        &self,
        sys: &DepletionSystem,
        dt: f64,
        unit: TimeUnit,
        method: crate::bateman::Method,
    ) -> Result<Self, Error> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(Error::BadStructure(format!("invalid timestep dt: {dt}")));
        }
        let dt_s = dt * unit.as_seconds();
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(Error::BadStructure(format!("invalid timestep dt: {dt}")));
        }
        let decay_sys = DepletionSystem::build(sys.chain.clone(), &ReactionRates::new())?;
        let n = decay_sys.chain.len();
        let mut n0 = vec![0.0; n];
        for (name, value) in &self.atoms {
            let idx = decay_sys
                .chain
                .index_of(name)
                .ok_or_else(|| Error::UnknownNuclide {
                    name: name.clone(),
                    context: "inventory",
                })?;
            require_finite_nonnegative(name, *value)?;
            n0[idx] = *value;
        }
        let solved = crate::solve_with_method(&decay_sys, method, &n0, dt_s)
            .map_err(|e| Error::BadStructure(e.to_string()))?;
        Ok(Self {
            atoms: decay_sys
                .chain
                .nuclides
                .iter()
                .zip(solved)
                .map(|(nuc, v)| (nuc.name.clone(), v))
                .collect(),
        })
    }

    /// Decay the inventory by `dt` in `unit` with a single CRAM-48 solve.
    ///
    /// Thin shim over [`DecayInventory::decay_with_method`] keeping the
    /// pre-`Method` call shape; new code should pass an explicit [`Method`].
    ///
    /// The solve is decay-only: the system is rebuilt from `sys.chain` with
    /// empty reaction rates. Result covers all chain nuclides (missing
    /// inputs start at zero); inventory names absent from the chain are an
    /// error.
    pub fn decay(&self, sys: &DepletionSystem, dt: f64, unit: TimeUnit) -> Result<Self, Error> {
        self.decay_with_method(
            sys,
            dt,
            unit,
            crate::bateman::Method::Cram(crate::Order::Order48),
        )
    }

    /// Activity per nuclide in `unit` (must be an activity unit).
    pub fn activities(
        &self,
        sys: &DepletionSystem,
        unit: QuantityUnit,
    ) -> Result<BTreeMap<String, f64>, Error> {
        if unit.kind() != QuantityKind::Activity {
            return Err(Error::BadStructure(format!(
                "unit `{}` is not an activity unit",
                unit.label()
            )));
        }
        let mut out = BTreeMap::new();
        for (name, n) in &self.atoms {
            let lam = decay_constant_for(sys, name)?;
            out.insert(name.clone(), (lam * n).max(0.0) / unit.base_factor());
        }
        Ok(out)
    }

    /// Mass per nuclide in `unit` (must be a mass unit).
    pub fn masses(&self, unit: QuantityUnit) -> Result<BTreeMap<String, f64>, Error> {
        if unit.kind() != QuantityKind::Mass {
            return Err(Error::BadStructure(format!(
                "unit `{}` is not a mass unit",
                unit.label()
            )));
        }
        let mut out = BTreeMap::new();
        for (name, n) in &self.atoms {
            let molar = molar_mass_for(name)?;
            out.insert(name.clone(), n / AVOGADRO * molar / unit.base_factor());
        }
        Ok(out)
    }

    /// Amount per nuclide in `unit` (must be a mole unit).
    pub fn moles(&self, unit: QuantityUnit) -> Result<BTreeMap<String, f64>, Error> {
        if unit.kind() != QuantityKind::Moles {
            return Err(Error::BadStructure(format!(
                "unit `{}` is not a mole unit",
                unit.label()
            )));
        }
        let mut out = BTreeMap::new();
        for (name, n) in &self.atoms {
            out.insert(name.clone(), n / AVOGADRO / unit.base_factor());
        }
        Ok(out)
    }

    /// Atom counts (copy of the underlying map).
    pub fn numbers(&self) -> BTreeMap<String, f64> {
        self.atoms.clone()
    }

    /// Fraction of total activity per nuclide (all-stable inventories give
    /// all-zero fractions rather than NaN).
    pub fn activity_fractions(
        &self,
        sys: &DepletionSystem,
    ) -> Result<BTreeMap<String, f64>, Error> {
        let mut acts = BTreeMap::new();
        for (name, n) in &self.atoms {
            let lam = decay_constant_for(sys, name)?;
            acts.insert(name.clone(), (lam * n).max(0.0));
        }
        Ok(normalize(acts))
    }

    /// Fraction of total mass per nuclide (empty/zero-mass inventories give
    /// all-zero fractions rather than NaN).
    pub fn mass_fractions(&self) -> Result<BTreeMap<String, f64>, Error> {
        let mut grams = BTreeMap::new();
        for (name, n) in &self.atoms {
            let molar = molar_mass_for(name)?;
            grams.insert(name.clone(), (n / AVOGADRO * molar).max(0.0));
        }
        Ok(normalize(grams))
    }

    /// Fraction of total atoms per nuclide (mole fractions; empty
    /// inventories give all-zero fractions rather than NaN).
    pub fn mole_fractions(&self) -> BTreeMap<String, f64> {
        normalize(self.atoms.clone())
    }

    /// Human-readable half-lives (`"5.27 y"`, `"3.2 d"`, `"stable"` for
    /// tabulated stable nuclides, `"unknown"` for names outside the nuclei
    /// table).
    pub fn half_lives_readable(&self) -> BTreeMap<String, String> {
        self.atoms
            .keys()
            .map(|name| {
                let text = match nucleide_nuclei::data::half_life_by_name(name) {
                    Some(t) if t > 0.0 && t.is_finite() => format_duration(t),
                    Some(_) => "stable".to_string(),
                    None => {
                        if nucleide_nuclei::NuclideId::from_name(name).is_ok() {
                            "stable".to_string()
                        } else {
                            "unknown".to_string()
                        }
                    }
                };
                (name.clone(), text)
            })
            .collect()
    }

    /// Element-wise sum (union of keys, missing entries count as zero).
    pub fn add(&self, other: &Self) -> Self {
        let mut atoms = self.atoms.clone();
        for (name, value) in &other.atoms {
            *atoms.entry(name.clone()).or_insert(0.0) += value;
        }
        Self { atoms }
    }

    /// Element-wise difference, clamped at zero (inventories stay physical).
    pub fn sub(&self, other: &Self) -> Self {
        let mut atoms = self.atoms.clone();
        for (name, value) in &other.atoms {
            let entry = atoms.entry(name.clone()).or_insert(0.0);
            *entry = (*entry - value).max(0.0);
        }
        Self { atoms }
    }

    /// Scale all counts by `s`.
    pub fn mul(&self, s: f64) -> Self {
        Self {
            atoms: self
                .atoms
                .iter()
                .map(|(name, value)| (name.clone(), value * s))
                .collect(),
        }
    }

    /// Divide all counts by `s` (IEEE semantics: division by zero yields
    /// infinities, matching plain `f64` arithmetic).
    pub fn div(&self, s: f64) -> Self {
        Self {
            atoms: self
                .atoms
                .iter()
                .map(|(name, value)| (name.clone(), value / s))
                .collect(),
        }
    }

    /// Serialize as CSV with a `nuclide,quantity,unit` header; quantities
    /// are atom counts (`atoms` unit).
    pub fn to_csv(&self) -> String {
        let mut out = String::from("nuclide,quantity,unit\n");
        for (name, value) in &self.atoms {
            out.push_str(&format!("{name},{value},atoms\n"));
        }
        out
    }

    /// Parse [`to_csv`](Self::to_csv) output (or hand-written
    /// `nuclide,quantity[,unit]` rows with an optional header line).
    ///
    /// Only atom-count rows round-trip without provider context: a missing
    /// unit means `atoms`, and any other unit is an error (converting
    /// activity/mass/mole rows needs decay constants and masses — use
    /// [`from_units`](Self::from_units) instead). Quantities must be finite
    /// and non-negative; duplicate nuclides are an error.
    pub fn from_csv(s: &str) -> Result<Self, Error> {
        let mut atoms = BTreeMap::new();
        for (lineno, raw) in s.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split(',').map(str::trim).collect();
            if lineno == 0
                && fields
                    .first()
                    .is_some_and(|h| h.eq_ignore_ascii_case("nuclide"))
            {
                continue;
            }
            if fields.len() < 2 || fields.len() > 3 {
                return Err(Error::BadStructure(format!(
                    "bad CSV row {lineno}: `{raw}`"
                )));
            }
            let name = fields[0].to_string();
            if name.is_empty() {
                return Err(Error::BadStructure(format!(
                    "bad CSV row {lineno}: empty nuclide name"
                )));
            }
            if atoms.contains_key(&name) {
                return Err(Error::BadStructure(format!(
                    "duplicate nuclide `{name}` in CSV"
                )));
            }
            let value: f64 = fields[1].parse().map_err(|_| {
                Error::BadStructure(format!("bad CSV quantity in row {lineno}: `{raw}`"))
            })?;
            require_finite_nonnegative(&name, value)?;
            if fields.len() == 3 {
                let unit = QuantityUnit::from_str(fields[2])?;
                if unit.kind() != QuantityKind::Atoms {
                    return Err(Error::BadStructure(format!(
                        "CSV unit `{}` needs provider context; use from_units",
                        fields[2].trim()
                    )));
                }
            }
            atoms.insert(name, value);
        }
        Ok(Self { atoms })
    }
}

fn normalize(map: BTreeMap<String, f64>) -> BTreeMap<String, f64> {
    let total: f64 = map.values().sum();
    if total > 0.0 && total.is_finite() {
        map.into_iter()
            .map(|(name, value)| (name, value / total))
            .collect()
    } else {
        map.into_keys().map(|name| (name, 0.0)).collect()
    }
}

fn format_duration(seconds: f64) -> String {
    let (value, unit) = if seconds < 60.0 {
        (seconds, "s")
    } else if seconds < 3_600.0 {
        (seconds / 60.0, "m")
    } else if seconds < SECONDS_PER_DAY {
        (seconds / 3_600.0, "h")
    } else if seconds < SECONDS_PER_YEAR {
        (seconds / SECONDS_PER_DAY, "d")
    } else {
        (seconds / SECONDS_PER_YEAR, "y")
    };
    format!("{} {unit}", format_sig(value, 3))
}

/// Format with `sig` significant figures (a minimal `{:.3g}` analogue:
/// fixed notation in a sane exponent range, scientific outside it, no
/// trailing zeros).
fn format_sig(value: f64, sig: u32) -> String {
    if value == 0.0 || !value.is_finite() {
        return format!("{value}");
    }
    let exp = value.abs().log10().floor() as i32;
    if !(-4..21).contains(&exp) {
        return format!("{:.prec$e}", value, prec = sig as usize - 1);
    }
    let decimals = (sig as i32 - 1 - exp).max(0) as usize;
    let raw = format!("{value:.decimals$}");
    if raw.contains('.') {
        raw.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        raw
    }
}

/// Time-integrated decays per nuclide over `[0, dt]` in chain order.
///
/// Diagonal Bateman integral: nuclide `i` reports
/// `n0[i] * (1 - exp(-λ_i dt))` with chain decay constants
/// ([`crate::integrate::decay_constants`]); in-growth from parents during
/// the interval is ignored. Stable nuclides report `0.0`.
pub fn cumulative_decays(sys: &DepletionSystem, n0: &[f64], dt: f64) -> Result<Vec<f64>, Error> {
    if n0.len() != sys.chain.len() {
        return Err(Error::BadStructure(format!(
            "n0 length {} != chain size {}",
            n0.len(),
            sys.chain.len()
        )));
    }
    if dt <= 0.0 || !dt.is_finite() {
        return Err(Error::BadStructure(format!("invalid timestep dt: {dt}")));
    }
    for (i, v) in n0.iter().enumerate() {
        if !v.is_finite() {
            return Err(Error::BadStructure(format!(
                "non-finite n0 at index {i}: {v}"
            )));
        }
    }
    Ok(crate::integrate::decay_constants(sys)
        .iter()
        .zip(n0)
        .map(|(lam, n)| {
            if *lam <= 0.0 {
                0.0
            } else {
                n * (-(-lam * dt).exp_m1())
            }
        })
        .collect())
}

/// Direct daughters of `name`: `(child, branching_ratio, mode kind)` in
/// chain order. Unknown names give an empty vector.
pub fn progeny(chain: &Chain, name: &str) -> Vec<(String, f64, String)> {
    let Some(idx) = chain.index_of(name) else {
        return Vec::new();
    };
    chain.nuclides[idx]
        .decay_modes
        .iter()
        .map(|m| (m.target.clone(), m.branching_ratio, m.kind.clone()))
        .collect()
}

/// Summed branching fraction from `parent` to `child` (`None` when the
/// parent is unknown or has no such decay edge).
pub fn branching_fraction(chain: &Chain, parent: &str, child: &str) -> Option<f64> {
    let idx = chain.index_of(parent)?;
    let sum: f64 = chain.nuclides[idx]
        .decay_modes
        .iter()
        .filter(|m| m.target == child)
        .map(|m| m.branching_ratio)
        .sum();
    // Distinguish "no edge" (sum of nothing) from a genuine zero branch.
    if chain.nuclides[idx]
        .decay_modes
        .iter()
        .any(|m| m.target == child)
    {
        Some(sum)
    } else {
        None
    }
}

/// Mode kind of the first `parent` → `child` decay edge (`None` when
/// absent).
pub fn decay_mode(chain: &Chain, parent: &str, child: &str) -> Option<String> {
    let idx = chain.index_of(parent)?;
    chain.nuclides[idx]
        .decay_modes
        .iter()
        .find(|m| m.target == child)
        .map(|m| m.kind.clone())
}

/// All decay edges as `(parent, child, branching_ratio, mode kind)` in
/// chain order.
pub fn chain_edges(chain: &Chain) -> Vec<(String, String, f64, String)> {
    let mut edges = Vec::new();
    for nuc in &chain.nuclides {
        for mode in &nuc.decay_modes {
            edges.push((
                nuc.name.clone(),
                mode.target.clone(),
                mode.branching_ratio,
                mode.kind.clone(),
            ));
        }
    }
    edges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::{ChainNuclide, DecayMode};

    const L1: f64 = 1.0e-6;
    const L2: f64 = 1.0e-5;

    fn hl(lam: f64) -> f64 {
        LN_2 / lam
    }

    fn abc_chain() -> Chain {
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(hl(L1)),
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "B".into(),
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            half_life: Some(hl(L2)),
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "C".into(),
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let c = ChainNuclide {
            name: "C".into(),
            ..Default::default()
        };
        Chain::from_nuclides(vec![a, b, c]).unwrap()
    }

    fn co60_sys() -> DepletionSystem {
        let t_half = nucleide_nuclei::data::half_life_by_name("Co60").unwrap();
        let co60 = ChainNuclide {
            name: "Co60".into(),
            half_life: Some(t_half),
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "Ni60".into(),
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let ni60 = ChainNuclide {
            name: "Ni60".into(),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![co60, ni60]).unwrap();
        DepletionSystem::build(chain, &ReactionRates::new()).unwrap()
    }

    fn map(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
        pairs.iter().map(|(k, v)| ((*k).to_string(), *v)).collect()
    }

    #[test]
    fn time_unit_strings() {
        assert_eq!(time_unit_from_str("s").unwrap(), TimeUnit::Seconds);
        assert_eq!(time_unit_from_str("SEC").unwrap(), TimeUnit::Seconds);
        assert_eq!(time_unit_from_str("seconds").unwrap(), TimeUnit::Seconds);
        assert_eq!(time_unit_from_str("m").unwrap(), TimeUnit::Minutes);
        assert_eq!(time_unit_from_str("Min").unwrap(), TimeUnit::Minutes);
        assert_eq!(time_unit_from_str("H").unwrap(), TimeUnit::Hours);
        assert_eq!(time_unit_from_str("hr").unwrap(), TimeUnit::Hours);
        assert_eq!(time_unit_from_str("d").unwrap(), TimeUnit::Days);
        assert_eq!(time_unit_from_str("Day").unwrap(), TimeUnit::Days);
        assert_eq!(time_unit_from_str("y").unwrap(), TimeUnit::Years);
        assert_eq!(time_unit_from_str("yr").unwrap(), TimeUnit::Years);
        assert_eq!(time_unit_from_str("YEARS").unwrap(), TimeUnit::Years);
        assert_eq!(TimeUnit::Years.as_seconds(), SECONDS_PER_YEAR);
        assert_eq!(TimeUnit::Days.as_seconds(), SECONDS_PER_DAY);
        assert!(time_unit_from_str("fortnight").is_err());
        assert!(time_unit_from_str("").is_err());
    }

    #[test]
    fn quantity_unit_strings() {
        assert_eq!(QuantityUnit::from_str("Bq").unwrap(), QuantityUnit::Bq);
        assert_eq!(QuantityUnit::from_str("kbq").unwrap(), QuantityUnit::KBq);
        assert_eq!(QuantityUnit::from_str("MBq").unwrap(), QuantityUnit::MBq);
        assert_eq!(QuantityUnit::from_str("gbq").unwrap(), QuantityUnit::GBq);
        assert_eq!(QuantityUnit::from_str("tbq").unwrap(), QuantityUnit::TBq);
        assert_eq!(QuantityUnit::from_str("PBq").unwrap(), QuantityUnit::PBq);
        assert_eq!(QuantityUnit::from_str("Ci").unwrap(), QuantityUnit::Ci);
        assert_eq!(QuantityUnit::from_str("curies").unwrap(), QuantityUnit::Ci);
        assert_eq!(QuantityUnit::from_str("mCi").unwrap(), QuantityUnit::MCi);
        assert_eq!(QuantityUnit::from_str("uCi").unwrap(), QuantityUnit::UCi);
        assert_eq!(QuantityUnit::from_str("µCi").unwrap(), QuantityUnit::UCi);
        assert_eq!(QuantityUnit::from_str("dpm").unwrap(), QuantityUnit::Dpm);
        assert_eq!(QuantityUnit::from_str("g").unwrap(), QuantityUnit::Gram);
        assert_eq!(
            QuantityUnit::from_str("KG").unwrap(),
            QuantityUnit::Kilogram
        );
        assert_eq!(
            QuantityUnit::from_str("mg").unwrap(),
            QuantityUnit::Milligram
        );
        assert_eq!(
            QuantityUnit::from_str("ug").unwrap(),
            QuantityUnit::Microgram
        );
        assert_eq!(QuantityUnit::from_str("t").unwrap(), QuantityUnit::Tonne);
        assert_eq!(QuantityUnit::from_str("mol").unwrap(), QuantityUnit::Mol);
        assert_eq!(QuantityUnit::from_str("mmol").unwrap(), QuantityUnit::Mmol);
        assert_eq!(
            QuantityUnit::from_str("atoms").unwrap(),
            QuantityUnit::Atoms
        );
        assert!(QuantityUnit::from_str("poodle").is_err());
        assert!(QuantityUnit::from_str("").is_err());
    }

    #[test]
    fn co60_unit_round_trip() {
        let sys = co60_sys();
        let inv =
            DecayInventory::from_units(&map(&[("Co60", 1.0)]), QuantityUnit::Ci, &sys).unwrap();
        // 1 Ci == 3.7e10 Bq.
        let bq = inv.activities(&sys, QuantityUnit::Bq).unwrap();
        let rel = (bq["Co60"] - 3.7e10).abs() / 3.7e10;
        assert!(rel < 1e-9, "{}", bq["Co60"]);

        // Bq -> Ci round trip.
        let back =
            DecayInventory::from_units(&map(&[("Co60", bq["Co60"])]), QuantityUnit::Bq, &sys)
                .unwrap();
        let ci = back.activities(&sys, QuantityUnit::Ci).unwrap();
        assert!((ci["Co60"] - 1.0).abs() < 1e-9);

        // Mass and moles agree with N = A / λ.
        let lam = LN_2 / nucleide_nuclei::data::half_life_by_name("Co60").unwrap();
        let expect_n = 3.7e10 / lam;
        let molar = nucleide_nuclei::data::atomic_mass_by_name("Co60").unwrap();
        let grams = inv.masses(QuantityUnit::Gram).unwrap();
        assert!((grams["Co60"] - expect_n / AVOGADRO * molar).abs() / grams["Co60"] < 1e-9);
        let moles = inv.moles(QuantityUnit::Mol).unwrap();
        assert!((moles["Co60"] - expect_n / AVOGADRO).abs() / moles["Co60"] < 1e-9);

        // g -> atoms -> mol closes the loop.
        let from_g =
            DecayInventory::from_units(&map(&[("Co60", grams["Co60"])]), QuantityUnit::Gram, &sys)
                .unwrap();
        assert!((from_g.atoms["Co60"] - expect_n).abs() / expect_n < 1e-9);
        let mmol = inv.moles(QuantityUnit::Mmol).unwrap();
        assert!((mmol["Co60"] - moles["Co60"] * 1.0e3).abs() / mmol["Co60"] < 1e-12);
        assert_eq!(inv.numbers()["Co60"], inv.atoms["Co60"]);
    }

    #[test]
    fn decay_matches_bateman() {
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let inv = DecayInventory::from_atoms(&map(&[("A", 1.0e15)])).unwrap();
        let dt = 1.0e5;
        let out = inv.decay(&sys, dt, TimeUnit::Seconds).unwrap();
        let na = 1.0e15 * (-L1 * dt).exp();
        let nb = 1.0e15 * L1 / (L2 - L1) * ((-L1 * dt).exp() - (-L2 * dt).exp());
        let nc = 1.0e15 - na - nb;
        assert!((out.atoms["A"] - na).abs() / na < 1e-8);
        assert!((out.atoms["B"] - nb).abs() / nb < 1e-7);
        assert!((out.atoms["C"] - nc).abs() / nc < 1e-8);
        // Time-unit plumbing: 1 day == 86400 s.
        let out_d = inv.decay(&sys, dt / 86_400.0, TimeUnit::Days).unwrap();
        assert!((out_d.atoms["A"] - na).abs() / na < 1e-8);
    }

    #[test]
    fn decay_with_method_matches_cram() {
        use crate::bateman::Method;
        let sys = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let inv = DecayInventory::from_atoms(&map(&[("A", 1.0e15)])).unwrap();
        let dt = 1.0e5;
        let cref = inv.decay(&sys, dt, TimeUnit::Seconds).unwrap();
        for method in [Method::Bateman, Method::BatemanHp] {
            let got = inv
                .decay_with_method(&sys, dt, TimeUnit::Seconds, method)
                .unwrap();
            for name in ["A", "B", "C"] {
                let rel = (got.atoms[name] - cref.atoms[name]).abs() / cref.atoms[name].max(1e-30);
                assert!(
                    rel < 1e-8,
                    "{name}: {} vs {}",
                    got.atoms[name],
                    cref.atoms[name]
                );
            }
        }
    }

    #[test]
    fn cumulative_matches_quadrature() {
        // Single-nuclide chain: diagonal integral is exact.
        let solo = ChainNuclide {
            name: "Solo".into(),
            half_life: Some(hl(L1)),
            ..Default::default()
        };
        let sys = DepletionSystem::build(
            Chain::from_nuclides(vec![solo]).unwrap(),
            &ReactionRates::new(),
        )
        .unwrap();
        let n0 = vec![1.0e12];
        let dt = 1.0e5;
        let got = cumulative_decays(&sys, &n0, dt).unwrap();
        // Closed form.
        let want = n0[0] * (1.0 - (-L1 * dt).exp());
        assert!((got[0] - want).abs() / want < 1e-12);
        // Numerical quadrature (Simpson, N(t) = N0 exp(-λt), decays = λN).
        let steps = 10_000usize;
        let h = dt / steps as f64;
        let mut sum = 0.0;
        for k in 0..=steps {
            let w = if k == 0 || k == steps {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += w * L1 * n0[0] * (-L1 * h * k as f64).exp();
        }
        let quad = sum * h / 3.0;
        assert!((got[0] - quad).abs() / quad < 1e-9);

        // Stable nuclide reports 0.0.
        let stable = ChainNuclide {
            name: "S".into(),
            ..Default::default()
        };
        let sys_s = DepletionSystem::build(
            Chain::from_nuclides(vec![stable]).unwrap(),
            &ReactionRates::new(),
        )
        .unwrap();
        assert_eq!(cumulative_decays(&sys_s, &[5.0], 100.0).unwrap(), vec![0.0]);

        // Chain case: parent (no in-growth) matches quadrature of CRAM N(t).
        let sys_abc = DepletionSystem::build(abc_chain(), &ReactionRates::new()).unwrap();
        let n0abc = vec![1.0e15, 0.0, 0.0];
        let got_abc = cumulative_decays(&sys_abc, &n0abc, dt).unwrap();
        assert!((got_abc[0] - n0abc[0] * (1.0 - (-L1 * dt).exp())).abs() / got_abc[0] < 1e-12);
        assert_eq!(got_abc[2], 0.0);
    }

    #[test]
    fn k40_style_branching_oracle() {
        // K-40-like topology: ~10.72% EC to Ar, ~89.28% beta- to Ca.
        let parent = ChainNuclide {
            name: "P".into(),
            half_life: Some(1.0e9),
            decay_modes: vec![
                DecayMode {
                    kind: "ec".into(),
                    target: "D1".into(),
                    branching_ratio: 0.1072,
                },
                DecayMode {
                    kind: "beta".into(),
                    target: "D2".into(),
                    branching_ratio: 0.8928,
                },
            ],
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![
            parent,
            ChainNuclide {
                name: "D1".into(),
                ..Default::default()
            },
            ChainNuclide {
                name: "D2".into(),
                ..Default::default()
            },
        ])
        .unwrap();
        let prog = progeny(&chain, "P");
        assert_eq!(prog.len(), 2);
        assert_eq!(prog[0], ("D1".to_string(), 0.1072, "ec".to_string()));
        assert_eq!(prog[1], ("D2".to_string(), 0.8928, "beta".to_string()));
        assert_eq!(branching_fraction(&chain, "P", "D1"), Some(0.1072));
        assert_eq!(branching_fraction(&chain, "P", "D2"), Some(0.8928));
        assert_eq!(branching_fraction(&chain, "P", "Nope"), None);
        assert_eq!(branching_fraction(&chain, "Nope", "D1"), None);
        assert_eq!(decay_mode(&chain, "P", "D1").as_deref(), Some("ec"));
        assert_eq!(decay_mode(&chain, "P", "D2").as_deref(), Some("beta"));
        assert_eq!(decay_mode(&chain, "P", "Nope"), None);
        assert!(progeny(&chain, "Nope").is_empty());
        assert!(progeny(&chain, "D1").is_empty());
        let edges = chain_edges(&chain);
        assert_eq!(
            edges,
            vec![
                ("P".to_string(), "D1".to_string(), 0.1072, "ec".to_string()),
                (
                    "P".to_string(),
                    "D2".to_string(),
                    0.8928,
                    "beta".to_string()
                ),
            ]
        );
    }

    #[test]
    fn csv_round_trip() {
        let inv = DecayInventory::from_atoms(&map(&[("Co60", 1.0e12), ("Ni60", 2.5)])).unwrap();
        let csv = inv.to_csv();
        let back = DecayInventory::from_csv(&csv).unwrap();
        assert_eq!(back, inv);
        // Headerless, unitless rows also parse.
        let bare = DecayInventory::from_csv("Co60,100\nNi60,200,atoms\n").unwrap();
        assert_eq!(bare.atoms["Co60"], 100.0);
        assert_eq!(bare.atoms["Ni60"], 200.0);
    }

    #[test]
    fn fractions_and_half_lives() {
        let sys = co60_sys();
        let inv = DecayInventory::from_atoms(&map(&[("Co60", 1.0e12), ("Ni60", 1.0e12)])).unwrap();
        // Co60 carries all activity (Ni60 stable).
        let af = inv.activity_fractions(&sys).unwrap();
        assert!((af["Co60"] - 1.0).abs() < 1e-12);
        assert_eq!(af["Ni60"], 0.0);
        // Mass/mole fractions sum to one.
        let mf = inv.mass_fractions().unwrap();
        assert!((mf.values().sum::<f64>() - 1.0).abs() < 1e-12);
        let mof = inv.mole_fractions();
        assert!((mof.values().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((mof["Co60"] - 0.5).abs() < 1e-12);
        // Readable half-lives: real isotope vs stable vs unknown.
        let hl_map =
            DecayInventory::from_atoms(&map(&[("Co60", 1.0), ("Ni60", 1.0), ("NotANuclide", 1.0)]))
                .unwrap()
                .half_lives_readable();
        assert!(hl_map["Co60"].ends_with(" y"), "{}", hl_map["Co60"]);
        assert_eq!(hl_map["Ni60"], "stable");
        assert_eq!(hl_map["NotANuclide"], "unknown");
        // All-stable activity fractions are zeros, not NaN.
        let stable_only = DecayInventory::from_atoms(&map(&[("Ni60", 3.0)])).unwrap();
        let af0 = stable_only.activity_fractions(&sys).unwrap();
        assert_eq!(af0["Ni60"], 0.0);
    }

    #[test]
    fn arithmetic() {
        let a = DecayInventory::from_atoms(&map(&[("X", 3.0), ("Y", 1.0)])).unwrap();
        let b = DecayInventory::from_atoms(&map(&[("X", 1.0), ("Z", 4.0)])).unwrap();
        let sum = a.add(&b);
        assert_eq!(sum.atoms["X"], 4.0);
        assert_eq!(sum.atoms["Y"], 1.0);
        assert_eq!(sum.atoms["Z"], 4.0);
        let diff = a.sub(&b);
        assert_eq!(diff.atoms["X"], 2.0);
        assert_eq!(diff.atoms["Y"], 1.0);
        // Clamped at zero: b has 4.0 of Z, a has none.
        assert_eq!(diff.atoms["Z"], 0.0);
        assert_eq!(a.mul(2.0).atoms["X"], 6.0);
        assert_eq!(a.div(2.0).atoms["X"], 1.5);
    }

    #[test]
    fn error_paths() {
        let sys = co60_sys();
        // Unknown nuclide in activity units.
        assert!(
            DecayInventory::from_units(&map(&[("Xx999", 1.0)]), QuantityUnit::Bq, &sys).is_err()
        );
        // Missing mass (synthetic name with no mass entry).
        let synth = DecayInventory::from_atoms(&map(&[("A", 1.0)])).unwrap();
        assert!(synth.masses(QuantityUnit::Gram).is_err());
        assert!(DecayInventory::from_units(&map(&[("A", 1.0)]), QuantityUnit::Gram, &sys).is_err());
        // Stable nuclide at non-zero activity.
        assert!(
            DecayInventory::from_units(&map(&[("Ni60", 1.0)]), QuantityUnit::Bq, &sys).is_err()
        );
        // Stable nuclide at zero activity is fine.
        assert!(DecayInventory::from_units(&map(&[("Ni60", 0.0)]), QuantityUnit::Bq, &sys).is_ok());
        // Bad unit strings.
        assert!(time_unit_from_str("eon").is_err());
        assert!(QuantityUnit::from_str("barn").is_err());
        // Wrong-kind units rejected.
        let inv = DecayInventory::from_atoms(&map(&[("Co60", 1.0)])).unwrap();
        assert!(inv.activities(&sys, QuantityUnit::Gram).is_err());
        assert!(inv.masses(QuantityUnit::Bq).is_err());
        assert!(inv.moles(QuantityUnit::Gram).is_err());
        // Invalid dt.
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(inv.decay(&sys, bad, TimeUnit::Seconds).is_err());
            assert!(cumulative_decays(&sys, &[1.0, 0.0], bad).is_err());
        }
        // Dimension mismatch.
        assert!(cumulative_decays(&sys, &[1.0], 10.0).is_err());
        // Inventory name outside the chain cannot decay.
        assert!(synth.decay(&sys, 1.0, TimeUnit::Seconds).is_err());
        // CSV errors: bad quantity, non-atoms unit, duplicates, negative.
        assert!(DecayInventory::from_csv("Co60,abc\n").is_err());
        assert!(DecayInventory::from_csv("Co60,1,Bq\n").is_err());
        assert!(DecayInventory::from_csv("Co60,1\nCo60,2\n").is_err());
        assert!(DecayInventory::from_csv("Co60,-1\n").is_err());
        assert!(DecayInventory::from_csv("Co60,NaN\n").is_err());
        assert!(DecayInventory::from_csv("onlyonefield\n").is_err());
        // Negative / non-finite construction rejected.
        assert!(DecayInventory::from_atoms(&map(&[("Co60", -1.0)])).is_err());
    }
}
