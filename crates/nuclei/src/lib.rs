//! Nuclide identification and naming conventions.
//!
//! Canonical representation is the `nucid`: a single
//! `u32` of the form `(Z*1000 + A) * 10_000 + state`, i.e. the zero-padded
//! six-digit ZZAAAM block followed by a four-digit tail holding the
//! metastable state (e.g. U-235 → 922350000, Am-242m → 952420001).
//! Chosen for compactness, hashing, and direct compatibility with the
//! integer ids used across legacy codes.
//!
//! Scope:
//! - id ↔ name ("U235", "Am242_m1") conversions
//! - id ↔ zzaaam (922350) conversions
//! - element symbol/number tables
//! - naming dialects (MCNP ZAID, Serpent, FLUKA, NIST, Cinder, ALARA), reaction names

use std::fmt;

pub mod armi;
pub mod data;
pub mod dialects;
pub mod particles;
pub mod rxname;

/// Element symbols indexed by atomic number (`ELEMENTS[z]`); index 0 is unused.
pub const ELEMENTS: [&str; 119] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Errors from nuclide parsing/validation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// Atomic number outside 1..=118.
    BadZ(u32),
    /// Mass number smaller than the atomic number.
    BadA {
        /// Atomic number.
        z: u32,
        /// Mass number.
        a: u32,
    },
    /// Mass number above the 3-digit AAA limit (> 999).
    MassNumberTooLarge(u32),
    /// Metastable state index above the supported range (> 9).
    BadState(u32),
    /// Name contained no digits (no mass number).
    MissingMassNumber(String),
    /// Mass number or state component failed to parse as an integer.
    BadNumber(String),
    /// Element symbol not recognized.
    UnknownElement(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::BadZ(z) => write!(f, "atomic number {z} out of range 1..=118"),
            Error::BadA { z, a } => write!(f, "mass number {a} < atomic number {z}"),
            Error::MassNumberTooLarge(a) => write!(f, "mass number {a} > 999 unsupported"),
            Error::BadState(s) => write!(f, "metastable state {s} > 9 unsupported"),
            Error::MissingMassNumber(s) => write!(f, "no mass number in name `{s}`"),
            Error::BadNumber(s) => write!(f, "invalid numeric component `{s}`"),
            Error::UnknownElement(s) => write!(f, "unknown element symbol `{s}`"),
        }
    }
}

impl std::error::Error for Error {}

/// A canonical nuclide identifier.
///
/// Layout (`nucid = (Z*1000 + A) * 10_000 + state`):
/// - H-1   → 10010000
/// - U-235 → 922350000
/// - Am-242m → 952420001
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NuclideId(u32);

impl NuclideId {
    /// Construct and validate a [`NuclideId`] from components.
    ///
    /// Enforces `1 <= Z <= 118`, `Z <= A <= 999`, and `S <= 9`. The
    /// `A <= 999` bound keeps the packed value within `u32` and preserves
    /// the 3-digit AAA invariant used by the zzaaam/zzllaaam dialects.
    pub const fn new(z: u32, a: u32, state: u32) -> Result<Self, Error> {
        if z == 0 || z > 118 {
            return Err(Error::BadZ(z));
        }
        if a < z {
            return Err(Error::BadA { z, a });
        }
        if a > 999 {
            return Err(Error::MassNumberTooLarge(a));
        }
        if state > 9 {
            return Err(Error::BadState(state));
        }
        Ok(Self((z * 1000 + a) * 10_000 + state))
    }

    /// Reconstruct from an existing nucid integer without validation.
    ///
    /// This is a raw bit-cast: invalid bit patterns yield meaningless
    /// components from [`z`](Self::z), [`a`](Self::a), and [`state`](Self::state).
    /// Use [`new`](Self::new) for validated construction.
    ///
    /// The id stays usable: [`to_name`](Self::to_name) and
    /// [`armi::nucid_to_armi_label`] render
    /// such ids with a diagnostic `Z{z}A{a}[m{s}]` fallback instead of
    /// panicking. Use [`try_from_nucid`](Self::try_from_nucid) (or
    /// [`is_valid`](Self::is_valid)) when the integer comes from untrusted
    /// input.
    pub const fn from_nucid(nucid: u32) -> Self {
        Self(nucid)
    }

    /// Reconstruct from a nucid integer with validation.
    ///
    /// Decomposes the integer into `(Z, A, state)` and applies the same
    /// `1 <= Z <= 118`, `Z <= A <= 999`, `S <= 9` checks as
    /// [`new`](Self::new) (kept in sync by inspection; the checks are
    /// inlined because `const fn` cannot match on the `Result`); integers
    /// with a non-canonical tail (the four low digits above 9, so no
    /// single-digit state can explain them) fail with [`Error::BadState`].
    /// Out-of-domain integers fail with the matching [`Error`] instead of
    /// producing an id whose name rendering falls back to the diagnostic form.
    pub const fn try_from_nucid(nucid: u32) -> Result<Self, Error> {
        let tail = nucid % 10_000;
        if tail > 9 {
            return Err(Error::BadState(tail));
        }
        let z = nucid / 10_000_000;
        let a = (nucid % 10_000_000) / 10_000;
        let state = nucid % 10;
        if z == 0 || z > 118 {
            return Err(Error::BadZ(z));
        }
        if a < z {
            return Err(Error::BadA { z, a });
        }
        if a > 999 {
            return Err(Error::MassNumberTooLarge(a));
        }
        if state > 9 {
            return Err(Error::BadState(state));
        }
        Ok(Self(nucid))
    }

    /// Whether this id decomposes into validated `(Z, A, state)` components
    /// (`1 <= Z <= 118`, `Z <= A <= 999`, `S <= 9`, canonical tail).
    ///
    /// Raw ids built by [`from_nucid`](Self::from_nucid) may fail this; every
    /// other constructor guarantees it.
    pub const fn is_valid(&self) -> bool {
        if self.0 % 10_000 > 9 {
            return false;
        }
        let z = self.z();
        let a = self.a();
        let state = self.state();
        z != 0 && z <= 118 && a >= z && a <= 999 && state <= 9
    }

    /// Raw nucid integer (`(Z*1000 + A)*10_000 + state`).
    pub const fn nucid(&self) -> u32 {
        self.0
    }

    /// Atomic number.
    pub const fn z(&self) -> u32 {
        self.0 / 10_000_000
    }

    /// Mass number.
    pub const fn a(&self) -> u32 {
        (self.0 % 10_000_000) / 10_000
    }

    /// Metastable state index (0 = ground).
    pub const fn state(&self) -> u32 {
        self.0 % 10
    }

    /// Six-digit ZZAAAM form (U-235 → 922350, Ba-137m → 561371).
    pub const fn zzaaam(&self) -> u32 {
        self.z() * 10_000 + self.a() * 10 + self.state()
    }

    /// Build from a six-digit ZZAAAM integer.
    pub fn from_zzaaam(v: u32) -> Result<Self, Error> {
        let state = v % 10;
        let rest = v / 10;
        let a = rest % 1_000;
        let z = rest / 1_000;
        Self::new(z, a, state)
    }

    /// Parse a name such as `"U235"`, `"U-235"`, `"u235"`, `"Am242_m1"`,
    /// `"Am-242m"`, `"Am242M"`, or `"Ba137m"`.
    ///
    /// Dashes are ignored and metastable markers are case-insensitive, so
    /// this matches PyNE's `name_to_id` normalization for the common forms.
    pub fn from_name(name: &str) -> Result<Self, Error> {
        let trimmed = name.trim();
        let cleaned: String = trimmed.chars().filter(|&c| c != '-').collect();
        let upper = cleaned.to_ascii_uppercase();
        let digit_start = upper
            .find(|c: char| c.is_ascii_digit())
            .ok_or_else(|| Error::MissingMassNumber(trimmed.to_string()))?;
        let sym_upper = &upper[..digit_start];
        let rest = &upper[digit_start..];

        let sym = canonicalize_symbol(sym_upper);
        let z = element_z(&sym).ok_or_else(|| Error::UnknownElement(sym_upper.to_string()))?;

        // Split mass number from an optional state suffix:
        // "235" | "242_M1" | "242M" | "137M"
        let (a_str, state_str) = if let Some((head, tail)) = rest.split_once('_') {
            // underscore form; tail may start with 'M'
            let tail = tail.strip_prefix('M').unwrap_or(tail);
            (head, Some(tail))
        } else if let Some((head, tail)) = rest.split_once('M') {
            // bare trailing-M form ("137M"); tail may hold the state index
            (head, Some(tail))
        } else {
            (rest, None)
        };

        let a: u32 = a_str
            .parse()
            .map_err(|_| Error::BadNumber(a_str.to_string()))?;
        let state = match state_str {
            None => 0,
            Some("") => 1,
            Some(n) => n.parse().map_err(|_| Error::BadNumber(format!("M{n}")))?,
        };

        Self::new(z, a, state)
    }

    /// GNDS-style name: `"U235"`, `"Am242_m1"`.
    ///
    /// Total over every raw id: validated ids render the canonical name
    /// (unchanged historical spelling, re-parseable by [`from_name`](Self::from_name));
    /// raw ids outside the validated `(Z, A, state)` domain render the
    /// diagnostic fallback `"Z{z}A{a}[m{s}]"`, which `from_name` does not
    /// parse. The fallback exists so display paths over unchecked integers
    /// (decay-table progeny, FFI) can never index `ELEMENTS` out of bounds.
    pub fn to_name(&self) -> String {
        if self.is_valid() {
            // `is_valid` pins `1 <= Z <= 118`, so this index is in bounds.
            let sym = ELEMENTS[self.z() as usize];
            match self.state() {
                0 => format!("{}{}", sym, self.a()),
                s => format!("{}{}_m{}", sym, self.a(), s),
            }
        } else {
            format!("Z{}A{}[m{}]", self.z(), self.a(), self.state())
        }
    }
}

impl fmt::Display for NuclideId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_name())
    }
}

impl std::str::FromStr for NuclideId {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        NuclideId::from_name(s)
    }
}

/// Element symbol for atomic number `z`, or `None`.
pub fn element_symbol(z: u32) -> Option<&'static str> {
    ELEMENTS
        .get(z as usize)
        .and_then(|s| if s.is_empty() { None } else { Some(*s) })
}

/// Atomic number for an element symbol (case-sensitive), or `None`.
pub fn element_z(symbol: &str) -> Option<u32> {
    ELEMENTS.iter().position(|s| *s == symbol).map(|z| z as u32)
}

/// Convert a free-form element symbol to canonical case for lookup.
fn canonicalize_symbol(sym: &str) -> String {
    let mut chars = sym.chars();
    let mut out = String::with_capacity(sym.len());
    if let Some(first) = chars.next() {
        out.extend(first.to_uppercase());
    }
    out.push_str(&chars.as_str().to_lowercase());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_ground_states() {
        assert_eq!(NuclideId::from_name("U235").unwrap().nucid(), 922_350_000);
        assert_eq!(NuclideId::from_name("H1").unwrap().nucid(), 10_010_000);
        assert_eq!(
            NuclideId::from_name("Og294").unwrap().nucid(),
            1_182_940_000
        );
    }

    #[test]
    fn parse_metastables() {
        let am = NuclideId::from_name("Am242_m1").unwrap();
        assert_eq!((am.z(), am.a(), am.state()), (95, 242, 1));
        assert_eq!(am.nucid(), 952_420_001);

        let ba = NuclideId::from_name("Ba137m").unwrap();
        assert_eq!((ba.z(), ba.a(), ba.state()), (56, 137, 1));
        assert_eq!(ba.zzaaam(), 561_371);
    }

    #[test]
    fn round_trip_display() {
        for name in ["U235", "H1", "Am242_m1", "Pu239"] {
            assert_eq!(NuclideId::from_name(name).unwrap().to_name(), name);
        }
    }

    #[test]
    fn zzaaam_round_trip() {
        let u5 = NuclideId::from_name("U235").unwrap();
        assert_eq!(u5.zzaaam(), 922_350);
        assert_eq!(
            NuclideId::from_zzaaam(922_350).map(|n| n.to_name()),
            Ok("U235".to_string())
        );
    }

    #[test]
    fn rejects_bad_input() {
        assert!(matches!(
            NuclideId::from_name("Xx999"),
            Err(Error::UnknownElement(_))
        ));
        assert!(matches!(
            NuclideId::from_name("U"),
            Err(Error::MissingMassNumber(_))
        ));
        assert!(matches!(NuclideId::new(0, 1, 0), Err(Error::BadZ(0))));
        assert!(matches!(NuclideId::new(6, 3, 0), Err(Error::BadA { .. })));
    }

    #[test]
    fn elements_table_sanity() {
        assert_eq!(element_z("U"), Some(92));
        assert_eq!(element_symbol(92), Some("U"));
        assert_eq!(element_z("Xx"), None);
    }

    #[test]
    fn rejects_mass_number_overflow() {
        assert!(matches!(
            NuclideId::new(92, 999_999, 0),
            Err(Error::MassNumberTooLarge(999_999))
        ));
        assert!(matches!(
            NuclideId::from_name("U999999"),
            Err(Error::MassNumberTooLarge(999_999))
        ));
    }

    #[test]
    fn parses_pyne_normalized_forms() {
        assert_eq!(NuclideId::from_name("U-235").unwrap().nucid(), 922_350_000);
        assert_eq!(NuclideId::from_name("u235").unwrap().nucid(), 922_350_000);
        assert_eq!(NuclideId::from_name("Am242M").unwrap().nucid(), 952_420_001);
        assert_eq!(
            NuclideId::from_name("Am-242M").unwrap().nucid(),
            952_420_001
        );
    }

    #[test]
    fn try_from_nucid_validates_raw_integers() {
        // Valid integers pass through untouched.
        assert_eq!(
            NuclideId::try_from_nucid(922_350_000).unwrap(),
            NuclideId::from_nucid(922_350_000)
        );
        assert!(NuclideId::try_from_nucid(10_010_000).unwrap().is_valid());
        assert!(NuclideId::from_nucid(922_350_000).is_valid());
        // Z out of ELEMENTS range (the `to_name` OOB family).
        assert!(!NuclideId::from_nucid(0).is_valid());
        assert!(matches!(NuclideId::try_from_nucid(0), Err(Error::BadZ(0))));
        assert!(matches!(
            NuclideId::try_from_nucid(u32::MAX),
            Err(Error::BadState(7295))
        ));
        assert!(matches!(
            NuclideId::try_from_nucid(1_190_000_000),
            Err(Error::BadZ(119))
        ));
        // A below Z, non-canonical tail (no single-digit state), Z=0 tail.
        assert!(matches!(
            NuclideId::try_from_nucid(920_050_000),
            Err(Error::BadA { z: 92, a: 5 })
        ));
        assert!(matches!(
            NuclideId::try_from_nucid(922_350_010),
            Err(Error::BadState(10))
        ));
        assert!(!NuclideId::from_nucid(920_050_000).is_valid());
        assert!(!NuclideId::from_nucid(922_350_010).is_valid());
    }

    #[test]
    fn to_name_falls_back_for_invalid_raw_ids() {
        // Formerly `ELEMENTS[z]` out-of-bounds panics; now diagnostics.
        assert_eq!(NuclideId::from_nucid(0).to_name(), "Z0A0[m0]");
        assert_eq!(
            NuclideId::from_nucid(u32::MAX).to_name(),
            format!(
                "Z{}A{}[m{}]",
                NuclideId::from_nucid(u32::MAX).z(),
                NuclideId::from_nucid(u32::MAX).a(),
                NuclideId::from_nucid(u32::MAX).state()
            )
        );
        assert_eq!(NuclideId::from_nucid(920_050_000).to_name(), "Z92A5[m0]");
        assert_eq!(NuclideId::from_nucid(922_350_010).to_name(), "Z92A235[m0]");
        // The fallback is diagnostic-only: `from_name` rejects it.
        for raw in [0, u32::MAX, 920_050_000, 1_190_000_000] {
            let name = NuclideId::from_nucid(raw).to_name();
            assert!(NuclideId::from_name(&name).is_err(), "{name}");
            assert_eq!(
                NuclideId::from_nucid(raw).to_string(),
                name,
                "Display follows to_name"
            );
        }
    }

    #[test]
    fn every_validated_id_round_trips_through_name() {
        // Canonical construction paths (incl. the checked raw-integer path).
        let mut ids = vec![
            NuclideId::new(1, 1, 0).unwrap(),
            NuclideId::new(92, 235, 0).unwrap(),
            NuclideId::new(95, 242, 9).unwrap(),
            NuclideId::new(118, 294, 0).unwrap(),
            NuclideId::from_name("Am242_m1").unwrap(),
            NuclideId::from_zzaaam(922_350).unwrap(),
            NuclideId::try_from_nucid(922_350_000).unwrap(),
        ];
        for z in [1, 2, 26, 92, 95, 118] {
            for a in [z, z + 1, 999] {
                for s in [0, 1, 9] {
                    if let Ok(id) = NuclideId::new(z, a.min(999), s) {
                        ids.push(id);
                    }
                }
            }
        }
        for id in ids {
            assert!(id.is_valid());
            assert_eq!(NuclideId::from_name(&id.to_name()).unwrap(), id);
        }
    }

    #[test]
    fn error_arms_construct_and_display() {
        assert!(matches!(
            NuclideId::new(92, 235, 10),
            Err(Error::BadState(10))
        ));
        assert!(matches!(
            NuclideId::from_name("U235_mX"),
            Err(Error::BadNumber(_))
        ));
        assert!(NuclideId::from_name("U235_mX")
            .unwrap_err()
            .to_string()
            .contains("MX"));
        assert!(matches!(
            NuclideId::from_name("U23X5"),
            Err(Error::BadNumber(_))
        ));
        assert!(!Error::BadState(10).to_string().is_empty());
        assert!(!Error::BadNumber("MX".to_string()).to_string().is_empty());
    }
}
