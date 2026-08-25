//! Naming-dialect conversions: MCNP ZAID, zzllaaam, Serpent, FLUKA, NIST,
//! Cinder, ALARA, sza.
//!
//! Design rule: free functions (or a local extension trait) defined here —
//! do NOT add methods to `NuclideId` in `lib.rs`; do NOT edit `lib.rs`.
//!
//! # Dialect conventions
//!
//! - MCNP ZAID: `Z*1000 + A`, metastable states add `300 + 100*S`
//!   (U-236m → 92636). Am-242/Am-242m are special-cased with swapped
//!   meanings: Am-242m → 95242, Am-242 → 95642. `from_zaid` applies the
//!   standard heuristic distributing `A - 400` excess into successive
//!   metastable states while `A/Z > 3.0` (95942 → Am-242 state 4).
//! - Isomer designator letters use the legacy sequence
//!   `"mnopqrstuvxyz"` (note: no `w`); state `S` maps to letter `S-1`.
//! - zzllaaam: `"ZZ-LL-AAAM"` + lowercase isomer letter
//!   (`"94-Pu-239"`, `"95-Am-242m"`, `"73-Ta-182n"`).
//! - Serpent: `"Ll-AAAM"` + lowercase isomer letter
//!   (`"Pu-239"`, `"Am-242m"`, `"U-236m"`).
//! - FLUKA: 8-character element/isotope names from the FLUKA translation
//!   table (e.g. `"HYDROG-1"`, `"235-U"`).
//!   Upstream quirks preserved: `TRITIUM` maps to H-4 (10040000) and the
//!   `LITHIUM` C literal is octal `030000000`; here it is read as decimal
//!   30_000_000 (Li), which is unreachable either way.
//! - NIST: mass number before symbol, no metastable flag (`"239Pu"`,
//!   `"242Am"`); parsing always yields the ground state.
//! - Cinder: `AAA*10_000 + Z*10 + S` (`aaazzzm`, U-235 → 2350920).
//! - ALARA: lowercase `"ll:AAA"`, no metastable flag (`"pu:239"`).
//! - SZA: `S*1_000_000 + Z*1_000 + A` (SSSZZZAAA; Am-242m → 1095242).
//!
//! # Intentionally skipped corners
//!
//! - Natural elemental nuclides (A = 0, e.g. `"U"` → 920000000) cannot be
//!   represented by [`crate::NuclideId`], which requires `A >= Z >= 1`;
//!   parsers report [`DialectError::NaturalElement`] there and formatters
//!   never emit the `"-nat"` / bare-symbol forms.
//! - Elemental group sets (LAN/ACT/TRU/MA/FP), `abun` tables, `zzzaaa`,
//!   GND, ENSDF, and the ENSDF state-id maps are out of scope.

use std::fmt;

use crate::{element_symbol, element_z, Error, NuclideId};

/// Errors raised by the dialect converters in this module.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DialectError {
    /// Numeric input admits no interpretation in the source dialect.
    NotANuclide(u32),
    /// Input is empty or carries no mass number where one is required.
    MissingMassNumber(String),
    /// Element symbol portion is not recognized.
    UnknownElement(String),
    /// Trailing metastable designator is not one of `mnopqrstuvxyz`.
    BadIsomerLetter(char),
    /// Input denotes a natural element, unrepresentable as a [`NuclideId`].
    NaturalElement(String),
    /// Name is absent from the vendored FLUKA table.
    UnknownFlukaName(String),
    /// Leading `ZZ` block disagrees with the element symbol.
    ZzSymbolMismatch { zz: u32, symbol: String },
    /// Component values were rejected by canonical validation.
    BadComponents(Error),
}

impl fmt::Display for DialectError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotANuclide(v) => {
                write!(f, "value {v} is not interpretable in the source dialect")
            }
            Self::MissingMassNumber(s) => write!(f, "no mass number in `{s}`"),
            Self::UnknownElement(s) => write!(f, "unknown element symbol `{s}`"),
            Self::BadIsomerLetter(c) => write!(f, "invalid metastable designator `{c}`"),
            Self::NaturalElement(s) => {
                write!(f, "natural element `{s}` unrepresentable as a NuclideId")
            }
            Self::UnknownFlukaName(s) => write!(f, "unknown FLUKA name `{s}`"),
            Self::ZzSymbolMismatch { zz, symbol } => {
                write!(
                    f,
                    "leading zz {zz} disagrees with element symbol `{symbol}`"
                )
            }
            Self::BadComponents(e) => write!(f, "invalid nuclide components: {e}"),
        }
    }
}

impl std::error::Error for DialectError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::BadComponents(e) => Some(e),
            _ => None,
        }
    }
}

impl From<Error> for DialectError {
    fn from(e: Error) -> Self {
        Self::BadComponents(e)
    }
}

/// Metastable-state designators (`"mnopqrstuvxyz"`).
/// State `S` (1-based) maps to the letter at index `S - 1`; `w` is absent
/// upstream as well.
const ISOMER_LETTERS: &str = "mnopqrstuvxyz";

/// FLUKA name table (names truncated to FLUKA's 8-character limit; values
/// are canonical nucids).
/// See the module docs for the two preserved upstream quirks.
const FLUKA_NAMES: &[(&str, u32)] = &[
    ("BERYLLIU", 40_000_000),
    ("BARIUM", 560_000_000),
    ("BOHRIUM", 1_070_000_000),
    ("BISMUTH", 830_000_000),
    ("BERKELIU", 970_000_000),
    ("BROMINE", 350_000_000),
    ("RUTHENIU", 440_000_000),
    ("RHENIUM", 750_000_000),
    ("RUTHERFO", 1_040_000_000),
    ("ROENTGEN", 1_110_000_000),
    ("RADIUM", 880_000_000),
    ("RUBIDIUM", 370_000_000),
    ("RADON", 860_000_000),
    ("RHODIUM", 450_000_000),
    ("THULIUM", 690_000_000),
    ("HYDROGEN", 10_000_000),
    ("PHOSPHO", 150_000_000),
    ("GERMANIU", 320_000_000),
    ("GADOLINI", 640_000_000),
    ("GALLIUM", 310_000_000),
    ("OSMIUM", 760_000_000),
    ("HASSIUM", 1_080_000_000),
    ("ZINC", 300_000_000),
    ("HOLMIUM", 670_000_000),
    ("HAFNIUM", 720_000_000),
    ("MERCURY", 800_000_000),
    ("HELIUM", 20_000_000),
    ("PRASEODY", 590_000_000),
    ("PLATINUM", 780_000_000),
    ("239-PU", 940_000_000),
    ("LEAD", 820_000_000),
    ("PROTACTI", 910_000_000),
    ("PALLADIU", 460_000_000),
    ("POLONIUM", 840_000_000),
    ("PROMETHI", 610_000_000),
    ("CARBON", 60_000_000),
    ("POTASSIU", 190_000_000),
    ("OXYGEN", 80_000_000),
    ("SULFUR", 160_000_000),
    ("TUNGSTEN", 740_000_000),
    ("EUROPIUM", 630_000_000),
    ("EINSTEIN", 990_000_000),
    ("ERBIUM", 680_000_000),
    ("MENDELEV", 1_010_000_000),
    ("MAGNESIU", 120_000_000),
    ("MOLYBDEN", 420_000_000),
    ("MANGANES", 250_000_000),
    ("MEITNERI", 1_090_000_000),
    ("URANIUM", 920_000_000),
    ("FRANCIUM", 870_000_000),
    ("IRON", 260_000_000),
    ("FERMIUM", 1_000_000_000),
    ("NICKEL", 280_000_000),
    ("NITROGEN", 70_000_000),
    ("NOBELIUM", 1_020_000_000),
    ("SODIUM", 110_000_000),
    ("NIOBIUM", 410_000_000),
    ("NEODYMIU", 600_000_000),
    ("NEON", 100_000_000),
    ("ZIRCONIU", 400_000_000),
    ("NEPTUNIU", 930_000_000),
    ("BORON", 50_000_000),
    ("COBALT", 270_000_000),
    ("CURIUM", 960_000_000),
    ("FLUORINE", 90_000_000),
    ("CALCIUM", 200_000_000),
    ("CALIFORN", 980_000_000),
    ("CERIUM", 580_000_000),
    ("CADMIUM", 480_000_000),
    ("VANADIUM", 230_000_000),
    ("CESIUM", 550_000_000),
    ("CHROMIUM", 240_000_000),
    ("COPPER", 290_000_000),
    ("STRONTIU", 380_000_000),
    ("KRYPTON", 360_000_000),
    ("SILICON", 140_000_000),
    ("TIN", 500_000_000),
    ("SAMARIUM", 620_000_000),
    ("SCANDIUM", 210_000_000),
    ("ANTIMONY", 510_000_000),
    ("SEABORGI", 1_060_000_000),
    ("SELENIUM", 340_000_000),
    ("YTTERBIU", 700_000_000),
    ("DUBNIUM", 1_050_000_000),
    ("DYSPROSI", 660_000_000),
    ("DARMSTAD", 1_100_000_000),
    ("LANTHANU", 570_000_000),
    ("CHLORINE", 170_000_000),
    ("LITHIUM", 30_000_000),
    ("THALLIUM", 810_000_000),
    ("LUTETIUM", 710_000_000),
    ("LAWRENCI", 1_030_000_000),
    ("THORIUM", 900_000_000),
    ("TITANIUM", 220_000_000),
    ("TELLURIU", 520_000_000),
    ("TERBIUM", 650_000_000),
    ("99-TC", 430_000_000),
    ("TANTALUM", 730_000_000),
    ("ACTINIUM", 890_000_000),
    ("SILVER", 470_000_000),
    ("IODINE", 530_000_000),
    ("IRIDIUM", 770_000_000),
    ("241-AM", 950_000_000),
    ("ALUMINUM", 130_000_000),
    ("ARSENIC", 330_000_000),
    ("ARGON", 180_000_000),
    ("GOLD", 790_000_000),
    ("ASTATINE", 850_000_000),
    ("INDIUM", 490_000_000),
    ("YTTRIUM", 390_000_000),
    ("XENON", 540_000_000),
    ("COPERNIC", 1_120_000_000),
    ("UNUNQUAD", 1_140_000_000),
    ("UNUNHEXI", 1_160_000_000),
    ("HYDROG-1", 10_010_000),
    ("DEUTERIU", 10_020_000),
    ("TRITIUM", 10_040_000),
    ("HELIUM-3", 20_030_000),
    ("HELIUM-4", 20_040_000),
    ("LITHIU-6", 30_060_000),
    ("LITHIU-7", 30_070_000),
    ("BORON-10", 50_100_000),
    ("BORON-11", 50_110_000),
    ("90-SR", 380_900_000),
    ("129-I", 531_290_000),
    ("124-XE", 541_240_000),
    ("126-XE", 541_260_000),
    ("128-XE", 541_280_000),
    ("130-XE", 541_300_000),
    ("131-XE", 541_310_000),
    ("132-XE", 541_320_000),
    ("134-XE", 541_340_000),
    ("135-XE", 541_350_000),
    ("136-XE", 541_360_000),
    ("135-CS", 551_350_000),
    ("137-CS", 551_370_000),
    ("230-TH", 902_300_000),
    ("232-TH", 902_320_000),
    ("233-U", 922_330_000),
    ("234-U", 922_340_000),
    ("235-U", 922_350_000),
    ("238-U", 922_380_000),
];

/// Lowercase isomer designator letter for state `s` (`1 → 'm'`), or `None`.
fn isomer_letter(state: u32) -> Option<char> {
    let idx = state.checked_sub(1)?;
    ISOMER_LETTERS.chars().nth(idx as usize)
}

/// State index for an isomer designator letter (`'m' → 1`), or `None`.
fn isomer_state(letter: char) -> Option<u32> {
    ISOMER_LETTERS
        .to_ascii_lowercase()
        .chars()
        .position(|c| c == letter.to_ascii_lowercase())
        .map(|p| p as u32 + 1)
}

/// Element symbol for a validated atomic number.
fn symbol_of(z: u32) -> &'static str {
    element_symbol(z).expect("NuclideId carries a validated atomic number")
}

/// Canonical (first letter upper, remainder lower) element lookup.
fn z_of_canonical_symbol(sym: &str) -> Option<u32> {
    if sym.is_empty() {
        return None;
    }
    let mut chars = sym.chars();
    let mut canon = String::with_capacity(sym.len());
    if let Some(first) = chars.next() {
        canon.extend(first.to_uppercase());
    }
    canon.push_str(&chars.as_str().to_lowercase());
    element_z(&canon)
}

/// Split `s` into its digit and alphabetic runs.
fn digit_letter_runs(s: &str) -> (String, String) {
    let digits: String = s.chars().filter(|c| c.is_ascii_digit()).collect();
    let letters: String = s.chars().filter(|c| c.is_ascii_alphabetic()).collect();
    (digits, letters)
}

/// Shared parser for dialects written as digits adjacent to an element
/// symbol with no metastable information (NIST, ALARA). Expects separators
/// already removed and `s` uppercased by the caller.
fn id_from_mass_symbol(raw: &str, s: &str) -> Result<NuclideId, DialectError> {
    let (digits, letters) = digit_letter_runs(s);
    if digits.is_empty() {
        return if z_of_canonical_symbol(&letters).is_some() {
            Err(DialectError::NaturalElement(raw.to_string()))
        } else {
            Err(DialectError::UnknownElement(raw.to_string()))
        };
    }
    let z = z_of_canonical_symbol(&letters)
        .ok_or_else(|| DialectError::UnknownElement(raw.to_string()))?;
    let a = digits
        .parse::<u32>()
        .map_err(|_| Error::BadNumber(digits.clone()))?;
    NuclideId::new(z, a, 0).map_err(DialectError::from)
}

/// Convert to the MCNP ZAID form (`to_zaid(U235) == 922350`).
///
/// Metastable states add `300 + 100*S`; Am-242m → 95242 and Am-242 → 95642
/// per the MCNP special case.
pub fn to_zaid(nuc: NuclideId) -> u32 {
    let mut state = nuc.state();
    let mut zaid = nuc.z() * 1_000 + nuc.a();
    if zaid == 95_242 && state < 2 {
        state = (state + 1) % 2;
    }
    if state != 0 {
        zaid += 300 + state * 100;
    }
    zaid
}

/// Interpret a MCNP ZAID (`from_zaid(95642)` → Am-242 ground state).
///
/// Parses MCNP ZAIDs, including the Am-242/242m swap and the
/// `A/Z > 3` metastable redistribution heuristic. Natural elements
/// (`AAA == 0`) are rejected as unrepresentable.
pub fn from_zaid(zaid: u32) -> Result<NuclideId, DialectError> {
    let z = zaid / 1_000;
    let a = zaid % 1_000;
    if z == 0 {
        return Err(DialectError::NotANuclide(zaid));
    }
    if z <= a {
        if a < 400 {
            return if zaid == 95_242 {
                NuclideId::new(95, 242, 1).map_err(DialectError::from)
            } else {
                NuclideId::new(z, a, 0).map_err(DialectError::from)
            };
        }
        if zaid == 95_642 {
            return NuclideId::new(95, 242, 0).map_err(DialectError::from);
        }
        let mut n = ((zaid - 400) * 10_000) + 1;
        loop {
            let aaa = (n / 10_000) % 1_000;
            let zzz = n / 10_000_000;
            if (aaa as f32) / (zzz as f32) <= 3.0 {
                break;
            }
            n -= 999_999;
        }
        return NuclideId::new(n / 10_000_000, (n / 10_000) % 1_000, n % 10)
            .map_err(DialectError::from);
    }
    if a == 0 {
        return Err(DialectError::NaturalElement(symbol_of(z).to_string()));
    }
    Err(DialectError::NotANuclide(zaid))
}

/// zzllaaam form: `"ZZ-LL-AAAM"` plus a lowercase isomer letter
/// (`zzllaaam(Am242m) == "95-Am-242m"`).
pub fn zzllaaam(nuc: NuclideId) -> String {
    let mut out = format!("{}-{}-{}", nuc.z(), symbol_of(nuc.z()), nuc.a());
    if let Some(c) = isomer_letter(nuc.state()) {
        out.push(c);
    }
    out
}

/// Parse `"ZZ-LL-AAAM"` (+ optional isomer letter, case-insensitive)
/// produced by [`zzllaaam`].
pub fn from_zzllaaam(name: &str) -> Result<NuclideId, DialectError> {
    let trimmed = name.trim();
    let parts: Vec<&str> = trimmed.split('-').collect();
    if parts.len() != 3 || trimmed.is_empty() {
        return Err(DialectError::MissingMassNumber(trimmed.to_string()));
    }
    let zz: u32 = parts[0]
        .parse()
        .map_err(|_| Error::BadNumber(parts[0].to_string()))?;
    let body = parts[2];
    if parts[1].eq_ignore_ascii_case("NAT") || body.eq_ignore_ascii_case("NAT") || body.is_empty() {
        return Err(DialectError::NaturalElement(trimmed.to_string()));
    }
    let (state, head) = split_isomer_suffix(body)?;
    if head.is_empty() {
        return Err(DialectError::MissingMassNumber(trimmed.to_string()));
    }
    let a = head
        .parse::<u32>()
        .map_err(|_| Error::BadNumber(head.to_string()))?;
    let z = z_of_canonical_symbol(parts[1])
        .ok_or_else(|| DialectError::UnknownElement(trimmed.to_string()))?;
    if z != zz {
        return Err(DialectError::ZzSymbolMismatch {
            zz,
            symbol: symbol_of(z).to_string(),
        });
    }
    NuclideId::new(z, a, state).map_err(DialectError::from)
}

/// Serpent form: `"Ll-AAAM"` plus a lowercase isomer letter
/// (`serpent(Am242m) == "Am-242m"`).
pub fn serpent(nuc: NuclideId) -> String {
    let mut out = format!("{}-{}", symbol_of(nuc.z()), nuc.a());
    if let Some(c) = isomer_letter(nuc.state()) {
        out.push(c);
    }
    out
}

/// Best-effort parse of a Serpent-style name (`"Am-242m"`, `"He-4"`).
///
/// Parses Serpent names: dashes are ignored, a trailing
/// isomer letter sets the state, and natural-element names are rejected.
pub fn from_serpent(name: &str) -> Result<NuclideId, DialectError> {
    let trimmed = name.trim();
    let s: String = trimmed
        .chars()
        .filter(|c| *c != '-')
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if s.is_empty() {
        return Err(DialectError::MissingMassNumber(trimmed.to_string()));
    }
    let (digits, letters) = digit_letter_runs(&s);
    if digits.is_empty() {
        let base = letters.strip_suffix("NAT").unwrap_or(&letters);
        return if z_of_canonical_symbol(base).is_some() {
            Err(DialectError::NaturalElement(trimmed.to_string()))
        } else {
            Err(DialectError::UnknownElement(trimmed.to_string()))
        };
    }
    let (state, head) = split_isomer_suffix(&s)?;
    let (a_digits, sym) = digit_letter_runs(head);
    if a_digits.is_empty() {
        return Err(DialectError::MissingMassNumber(trimmed.to_string()));
    }
    let z = z_of_canonical_symbol(&sym)
        .ok_or_else(|| DialectError::UnknownElement(trimmed.to_string()))?;
    let a = a_digits
        .parse::<u32>()
        .map_err(|_| Error::BadNumber(a_digits.clone()))?;
    NuclideId::new(z, a, state).map_err(DialectError::from)
}

/// FLUKA material name for `nuc` (e.g. `id_to_fluka(U235) == "235-U"`).
///
/// Only nuclides with an explicit entry in the vendored table resolve;
/// the natural-element rows can never match a [`NuclideId`].
pub fn id_to_fluka(nuc: NuclideId) -> Result<&'static str, DialectError> {
    FLUKA_NAMES
        .iter()
        .find(|(_, nucid)| *nucid == nuc.nucid())
        .map(|(name, _)| *name)
        .ok_or_else(|| DialectError::UnknownFlukaName(nuc.to_name()))
}

/// Resolve a FLUKA name from the vendored table to a [`NuclideId`]
/// (`fluka_to_id("LITHIU-7")` → Li-7).
///
/// Rows denoting natural elements (A = 0) yield
/// [`DialectError::NaturalElement`] rather than an id.
pub fn fluka_to_id(name: &str) -> Result<NuclideId, DialectError> {
    let &(_, nucid) = FLUKA_NAMES
        .iter()
        .find(|(known, _)| *known == name)
        .ok_or_else(|| DialectError::UnknownFlukaName(name.to_string()))?;
    let z = nucid / 10_000_000;
    let a = (nucid / 10_000) % 1_000;
    if a == 0 {
        return Err(DialectError::NaturalElement(name.to_string()));
    }
    NuclideId::new(z, a, nucid % 10).map_err(DialectError::from)
}

/// NIST form: mass number followed by the element symbol, metastable state
/// dropped (`nist(Am242m) == "242Am"`).
pub fn nist(nuc: NuclideId) -> String {
    format!("{}{}", nuc.a(), symbol_of(nuc.z()))
}

/// Parse a NIST-style name (`"239Pu"`, `"4He"`); the result is always a
/// ground state because the dialect carries no state information.
pub fn nist_to_id(name: &str) -> Result<NuclideId, DialectError> {
    let trimmed = name.trim();
    let upper = trimmed.to_ascii_uppercase();
    id_from_mass_symbol(trimmed, &upper)
}

/// Cinder `AAAZZZM` form: `A*10_000 + Z*10 + S` (`to_cinder(U235) == 2350920`).
pub fn to_cinder(nuc: NuclideId) -> u32 {
    nuc.a() * 10_000 + nuc.z() * 10 + nuc.state()
}

/// Interpret a Cinder integer (`from_cinder(2420951)` → Am-242m).
pub fn from_cinder(value: u32) -> Result<NuclideId, DialectError> {
    let state = value % 10;
    let aaazzz = value / 10;
    let z = aaazzz % 1_000;
    let a = aaazzz / 1_000;
    NuclideId::new(z, a, state).map_err(DialectError::from)
}

/// ALARA form: lowercase `"ll:AAA"` with no metastable flag
/// (`alara(Pu239) == "pu:239"`).
pub fn alara(nuc: NuclideId) -> String {
    format!("{}:{}", symbol_of(nuc.z()).to_ascii_lowercase(), nuc.a())
}

/// Parse an ALARA-style name (`"pu:239"`, `"he:4"`); the result is always a
/// ground state because the dialect carries no state information.
pub fn alara_to_id(name: &str) -> Result<NuclideId, DialectError> {
    let trimmed = name.trim();
    let cleaned: String = trimmed
        .chars()
        .filter(|c| *c != ':')
        .map(|c| c.to_ascii_uppercase())
        .collect();
    id_from_mass_symbol(trimmed, &cleaned)
}

/// SZA form: `S*1_000_000 + Z*1_000 + A` (`to_sza(Am242m) == 1095242`).
pub fn to_sza(nuc: NuclideId) -> u32 {
    nuc.state() * 1_000_000 + nuc.z() * 1_000 + nuc.a()
}

/// Interpret an SZA integer (`from_sza(1095242)` → Am-242m).
pub fn from_sza(value: u32) -> Result<NuclideId, DialectError> {
    let state = value / 1_000_000;
    let zzzaaa = value % 1_000_000;
    let z = zzzaaa / 1_000;
    let a = zzzaaa % 1_000;
    NuclideId::new(z, a, state).map_err(DialectError::from)
}

/// Split a trailing isomer designator off `body`, returning the state and
/// the remaining prefix. A trailing digit means the ground state.
fn split_isomer_suffix(body: &str) -> Result<(u32, &str), DialectError> {
    let last = body
        .chars()
        .next_back()
        .ok_or_else(|| DialectError::MissingMassNumber(body.to_string()))?;
    if last.is_ascii_digit() {
        return Ok((0, body));
    }
    match isomer_state(last) {
        Some(state) => {
            let head = &body[..body.len() - last.len_utf8()];
            Ok((state, head))
        }
        None => Err(DialectError::BadIsomerLetter(last)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nid(z: u32, a: u32, s: u32) -> NuclideId {
        NuclideId::new(z, a, s).unwrap()
    }

    #[test]
    fn zaid_matches_mcnp_convention() {
        assert_eq!(to_zaid(nid(1, 1, 0)), 1001);
        assert_eq!(to_zaid(nid(92, 235, 0)), 92_235);
        assert_eq!(to_zaid(nid(92, 236, 1)), 92_636);
        assert_eq!(to_zaid(nid(95, 242, 1)), 95_242);
        assert_eq!(to_zaid(nid(95, 242, 0)), 95_642);
        assert_eq!(to_zaid(nid(2, 4, 0)), 2004);
    }

    #[test]
    fn zaid_edge_cases() {
        assert_eq!(from_zaid(2004).unwrap(), nid(2, 4, 0));
        assert_eq!(from_zaid(1001).unwrap(), nid(1, 1, 0));
        assert_eq!(from_zaid(95_242).unwrap(), nid(95, 242, 1));
        assert_eq!(from_zaid(95_642).unwrap(), nid(95, 242, 0));
        assert_eq!(from_zaid(92_636).unwrap(), nid(92, 236, 1));
        assert_eq!(from_zaid(95_942).unwrap(), nid(95, 242, 4));
        assert_eq!(from_zaid(96_644).unwrap(), nid(96, 244, 1));
    }

    #[test]
    fn zaid_round_trips() {
        for nuc in [
            nid(1, 1, 0),
            nid(92, 235, 0),
            nid(92, 236, 1),
            nid(95, 242, 0),
            nid(95, 242, 1),
            nid(95, 242, 4),
            nid(56, 137, 1),
            nid(73, 182, 2),
        ] {
            assert_eq!(from_zaid(to_zaid(nuc)).unwrap(), nuc);
        }
    }

    #[test]
    fn zaid_error_paths() {
        assert_eq!(from_zaid(92), Err(DialectError::NotANuclide(92)));
        assert_eq!(
            from_zaid(92_000),
            Err(DialectError::NaturalElement("U".to_string()))
        );
        assert_eq!(from_zaid(50_003), Err(DialectError::NotANuclide(50_003)));
    }

    #[test]
    fn zzllaaam_round_trip() {
        assert_eq!(zzllaaam(nid(94, 239, 0)), "94-Pu-239");
        assert_eq!(zzllaaam(nid(95, 242, 1)), "95-Am-242m");
        assert_eq!(zzllaaam(nid(95, 242, 0)), "95-Am-242");
        assert_eq!(zzllaaam(nid(92, 236, 1)), "92-U-236m");
        assert_eq!(zzllaaam(nid(73, 182, 2)), "73-Ta-182n");
    }

    #[test]
    fn zzllaaam_parses() {
        assert_eq!(from_zzllaaam("94-Pu-239").unwrap(), nid(94, 239, 0));
        assert_eq!(from_zzllaaam("95-Am-242m").unwrap(), nid(95, 242, 1));
        assert_eq!(from_zzllaaam("73-Ta-182n").unwrap(), nid(73, 182, 2));
        assert_eq!(from_zzllaaam("95-am-242m").unwrap(), nid(95, 242, 1));
    }

    #[test]
    fn zzllaaam_error_paths() {
        assert_eq!(
            from_zzllaaam("Ta-182b"),
            Err(DialectError::MissingMassNumber("Ta-182b".to_string()))
        );
        assert_eq!(
            from_zzllaaam("95-Pu-239"),
            Err(DialectError::ZzSymbolMismatch {
                zz: 95,
                symbol: "Pu".to_string()
            })
        );
        assert_eq!(
            from_zzllaaam("92-U-nat"),
            Err(DialectError::NaturalElement("92-U-nat".to_string()))
        );
        assert!(matches!(
            from_zzllaaam("94-Xx-239"),
            Err(DialectError::UnknownElement(_))
        ));
    }

    #[test]
    fn serpent_dialect_round_trip() {
        assert_eq!(serpent(nid(94, 239, 0)), "Pu-239");
        assert_eq!(serpent(nid(95, 242, 1)), "Am-242m");
        assert_eq!(serpent(nid(95, 242, 0)), "Am-242");
        assert_eq!(serpent(nid(92, 236, 1)), "U-236m");
        assert_eq!(serpent(nid(73, 182, 2)), "Ta-182n");
    }

    #[test]
    fn serpent_parses() {
        assert_eq!(from_serpent("Pu-239").unwrap(), nid(94, 239, 0));
        assert_eq!(from_serpent("Am-242m").unwrap(), nid(95, 242, 1));
        assert_eq!(from_serpent("He-4").unwrap(), nid(2, 4, 0));
        assert_eq!(from_serpent("U-236m").unwrap(), nid(92, 236, 1));
        assert_eq!(from_serpent("Cm-244m").unwrap(), nid(96, 244, 1));
        assert_eq!(from_serpent("Ta-182n").unwrap(), nid(73, 182, 2));
    }

    #[test]
    fn serpent_error_paths() {
        assert_eq!(
            from_serpent("U-nat"),
            Err(DialectError::NaturalElement("U-nat".to_string()))
        );
        assert_eq!(
            from_serpent("Am-242j"),
            Err(DialectError::BadIsomerLetter('J'))
        );
        assert!(matches!(
            from_serpent("Xx-12"),
            Err(DialectError::UnknownElement(_))
        ));
    }

    #[test]
    fn fluka_to_id_isotopes() {
        assert_eq!(fluka_to_id("LITHIU-7").unwrap(), nid(3, 7, 0));
        assert_eq!(fluka_to_id("HYDROG-1").unwrap(), nid(1, 1, 0));
        assert_eq!(fluka_to_id("235-U").unwrap(), nid(92, 235, 0));
        assert_eq!(fluka_to_id("BORON-10").unwrap(), nid(5, 10, 0));
        assert_eq!(fluka_to_id("HELIUM-4").unwrap(), nid(2, 4, 0));
    }

    #[test]
    fn fluka_to_id_error_paths() {
        assert_eq!(
            fluka_to_id("NOPE"),
            Err(DialectError::UnknownFlukaName("NOPE".to_string()))
        );
        assert_eq!(
            fluka_to_id("URANIUM"),
            Err(DialectError::NaturalElement("URANIUM".to_string()))
        );
    }

    #[test]
    fn id_to_fluka_names() {
        assert_eq!(id_to_fluka(nid(3, 7, 0)), Ok("LITHIU-7"));
        assert_eq!(id_to_fluka(nid(92, 235, 0)), Ok("235-U"));
        assert_eq!(id_to_fluka(nid(1, 1, 0)), Ok("HYDROG-1"));
        assert_eq!(id_to_fluka(nid(2, 4, 0)), Ok("HELIUM-4"));
        assert_eq!(
            id_to_fluka(nid(26, 56, 0)),
            Err(DialectError::UnknownFlukaName("Fe56".to_string()))
        );
    }

    #[test]
    fn nist_dialect_drops_state() {
        assert_eq!(nist(nid(94, 239, 0)), "239Pu");
        assert_eq!(nist(nid(95, 242, 1)), "242Am");
        assert_eq!(nist(nid(2, 4, 0)), "4He");
    }

    #[test]
    fn nist_parses_ground_states() {
        assert_eq!(nist_to_id("4He").unwrap(), nid(2, 4, 0));
        assert_eq!(nist_to_id("244Cm").unwrap(), nid(96, 244, 0));
        assert_eq!(nist_to_id("239Pu").unwrap(), nid(94, 239, 0));
        assert_eq!(nist_to_id("242Am").unwrap(), nid(95, 242, 0));
        assert_eq!(
            nist_to_id("U"),
            Err(DialectError::NaturalElement("U".to_string()))
        );
        assert!(matches!(
            nist_to_id("242Xx"),
            Err(DialectError::UnknownElement(_))
        ));
    }

    #[test]
    fn cinder_dialect_round_trip() {
        assert_eq!(to_cinder(nid(1, 2, 0)), 20_010);
        assert_eq!(to_cinder(nid(95, 242, 1)), 2_420_951);
        assert_eq!(to_cinder(nid(92, 236, 1)), 2_360_921);
        assert_eq!(from_cinder(2_420_951).unwrap(), nid(95, 242, 1));
        assert_eq!(from_cinder(2_360_921).unwrap(), nid(92, 236, 1));
        assert_eq!(from_cinder(2_440_961).unwrap(), nid(96, 244, 1));
        assert!(matches!(
            from_cinder(20),
            Err(DialectError::BadComponents(Error::BadA { .. }))
        ));
    }

    #[test]
    fn alara_dialect_round_trip() {
        assert_eq!(alara(nid(94, 239, 0)), "pu:239");
        assert_eq!(alara(nid(95, 242, 1)), "am:242");
        assert_eq!(alara(nid(2, 4, 0)), "he:4");
        assert_eq!(alara(nid(92, 236, 1)), "u:236");
    }

    #[test]
    fn alara_parses_ground_states() {
        assert_eq!(alara_to_id("pu:239").unwrap(), nid(94, 239, 0));
        assert_eq!(alara_to_id("cm:244").unwrap(), nid(96, 244, 0));
        assert_eq!(alara_to_id("he:4").unwrap(), nid(2, 4, 0));
        assert_eq!(
            alara_to_id("u"),
            Err(DialectError::NaturalElement("u".to_string()))
        );
        assert!(matches!(
            alara_to_id("zz:10"),
            Err(DialectError::UnknownElement(_))
        ));
    }

    #[test]
    fn sza_dialect_round_trip() {
        assert_eq!(to_sza(nid(2, 4, 0)), 2004);
        assert_eq!(to_sza(nid(95, 242, 1)), 1_095_242);
        assert_eq!(to_sza(nid(92, 236, 1)), 1_092_236);
        assert_eq!(to_sza(nid(95, 242, 4)), 4_095_242);
        assert_eq!(from_sza(1_095_242).unwrap(), nid(95, 242, 1));
        assert_eq!(from_sza(2004).unwrap(), nid(2, 4, 0));
        assert_eq!(from_sza(1_096_244).unwrap(), nid(96, 244, 1));
        assert!(matches!(
            from_sza(20),
            Err(DialectError::BadComponents(Error::BadZ(0)))
        ));
    }

    #[test]
    fn dialect_error_display_smoke() {
        let cases = [
            DialectError::NotANuclide(7),
            DialectError::MissingMassNumber("U".into()),
            DialectError::UnknownElement("Xx".into()),
            DialectError::BadIsomerLetter('b'),
            DialectError::NaturalElement("U".into()),
            DialectError::UnknownFlukaName("NOPE".into()),
            DialectError::ZzSymbolMismatch {
                zz: 95,
                symbol: "Pu".into(),
            },
            DialectError::BadComponents(Error::BadZ(0)),
        ];
        for e in cases {
            assert!(!e.to_string().is_empty());
        }
    }
}
