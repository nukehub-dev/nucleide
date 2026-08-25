//! Chemical-formula parsing and element ↔ nuclide expansion.
//!
//! This module bridges the nuclide-level [`Material`] model and the
//! chemistry-level world of formulas such as `"H2O"` or `"Fe2(SO4)3"`, and
//! the MCNP-style natural-element placeholder ids (`z * 10_000_000`,
//! i.e. zaid `z*1000` with `AAA == 0`) produced by `mcnp-io`.
//!
//! Like [`crate::MassProvider`], isotope data is injected through a provider
//! trait ([`AbundanceProvider`]) so nothing here hard-depends on data
//! availability: [`NoAbundances`] fails explicitly with
//! [`FormulaError::NoAbundanceData`], while [`NaturalAbundances`] serves
//! the tabulated natural abundances from `nuclei::data`.
//!
//! ```
//! use material::{Ame2020, Material, NaturalAbundances};
//!
//! let water = Material::from_formula(
//!     "H2O",
//!     &Ame2020,
//!     &NaturalAbundances,
//!     Some(1.0),
//! )
//! .unwrap();
//! let atoms = water.atom_fractions(&Ame2020).unwrap();
//! let h = atoms.iter().filter(|(id, _)| id.z() == 1).map(|(_, f)| f).sum::<f64>();
//! assert!((h - 2.0 / 3.0).abs() < 1e-12);
//! ```

use std::collections::BTreeMap;
use std::sync::OnceLock;

use nuclei::{element_z, NuclideId};
use thiserror::Error;

use crate::Material;
/// Result alias for formula parsing and expansion.
pub type FormulaResult<T> = std::result::Result<T, FormulaError>;

/// Errors from formula parsing and element expansion.
///
/// Distinct from the crate-level [`crate::Error`] so that the parser can
/// report byte positions and unknown symbols without polluting the shared
/// composition-error set; mass-table failures are wrapped in [`Self::Core`].
#[derive(Debug, Error)]
pub enum FormulaError {
    /// The formula is not valid under the supported grammar (see
    /// [`parse_formula`]).
    #[error("formula syntax error at byte {pos}: {message}")]
    ParseError {
        /// Byte offset of the offending character.
        pos: usize,
        /// Human-readable explanation.
        message: String,
    },
    /// An element symbol was not recognized (case-sensitive).
    #[error("unknown element symbol `{0}`")]
    UnknownElement(String),
    /// An element has no tabulated natural isotopes to expand into.
    #[error("no natural abundance data for element Z={0}")]
    NoAbundanceData(u32),
    /// A composition-level failure (e.g. missing atomic mass).
    #[error(transparent)]
    Core(#[from] crate::Error),
}

// ---------------------------------------------------------------------------
// Formula grammar
//
//   formula := segment ("." | "·") segment*        hydrate parts
//   segment := integer? unit*                      leading int only for
//                                                  non-first segments
//   unit     := element integer? | "(" segment ")" integer?
//   element := [A-Z][a-z]?
//
// Examples: H2O, C6H12O6, Ca(OH)2, Fe2(SO4)3, CH3(CH2)6CH3,
// CuSO4·5H2O == CuSO4.5H2O.
// ---------------------------------------------------------------------------

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn parse_error(&self, pos: usize, message: impl Into<String>) -> FormulaError {
        FormulaError::ParseError {
            pos,
            message: message.into(),
        }
    }

    fn error_here(&self, message: impl Into<String>) -> FormulaError {
        self.parse_error(self.pos, message)
    }

    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    /// True when the cursor sits on a hydrate separator (`.` ASCII or `·`
    /// U+00B7, encoded as `C2 B7` in UTF-8).
    fn at_hyphen_dot(&self) -> bool {
        match self.peek() {
            Some(b'.') => true,
            Some(0xC2) => self.src.get(self.pos + 1) == Some(&0xB7),
            _ => false,
        }
    }

    fn advance_over_separator(&mut self) {
        self.pos += if self.src[self.pos] == b'.' { 1 } else { 2 };
    }

    /// Consume a run of ASCII digits as a count, or `None` if absent.
    fn take_count(&mut self) -> FormulaResult<Option<f64>> {
        let start = self.pos;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.pos += 1;
        }
        if start == self.pos {
            return Ok(None);
        }
        let text = std::str::from_utf8(&self.src[start..self.pos]).unwrap_or_default();
        text.parse::<u64>()
            .map(|n| Some(n as f64))
            .map_err(|_| self.parse_error(start, "count too large"))
    }

    /// Consume an element symbol starting at an uppercase byte and its
    /// optional trailing count; returns `(Z, count)`.
    fn take_element(&mut self) -> FormulaResult<(u32, f64)> {
        let start = self.pos;
        self.pos += 1;
        // Prefer a two-letter symbol when a lowercase letter follows.
        let two_letter = self
            .src
            .get(start..start + 2)
            .filter(|bytes| bytes[1].is_ascii_lowercase())
            .and_then(|bytes| std::str::from_utf8(bytes).ok());
        let one_letter = std::str::from_utf8(&self.src[start..start + 1]).ok();
        let (z, matched_len) = match two_letter.and_then(element_z) {
            Some(z) => (Some(z), 2),
            None => (one_letter.and_then(element_z), 1),
        };
        let z = z.ok_or_else(|| {
            let candidate = two_letter.unwrap_or_else(|| one_letter.unwrap_or("?"));
            FormulaError::UnknownElement(candidate.to_string())
        })?;
        self.pos = start + matched_len;
        let count = self.take_count()?.unwrap_or(1.0);
        Ok((z, count))
    }

    /// Consume units until end-of-input, a closing parenthesis, or a hydrate
    /// separator, accumulating counts scaled by `scale` into `out`.
    fn take_units(&mut self, out: &mut BTreeMap<u32, f64>, scale: f64) -> FormulaResult<()> {
        while let Some(c) = self.peek() {
            match c {
                b')' | b'.' => break,
                0xC2 if self.at_hyphen_dot() => break,
                b'(' => {
                    self.pos += 1;
                    let mut inner = BTreeMap::new();
                    self.take_units(&mut inner, 1.0)?;
                    if self.peek() != Some(b')') {
                        return Err(self.error_here("unbalanced parenthesis: expected `)`"));
                    }
                    self.pos += 1;
                    let mult = self.take_count()?.unwrap_or(1.0);
                    for (z, count) in inner {
                        *out.entry(z).or_insert(0.0) += count * mult * scale;
                    }
                }
                b'A'..=b'Z' => {
                    let (z, count) = self.take_element()?;
                    if count > 0.0 {
                        *out.entry(z).or_insert(0.0) += count * scale;
                    }
                }
                b'0'..=b'9' => {
                    return Err(self.error_here("count without a preceding element"));
                }
                _ => {
                    let ch = std::str::from_utf8(&self.src[self.pos..])
                        .unwrap_or("?")
                        .chars()
                        .next()
                        .unwrap_or('?');
                    return Err(self.error_here(format!("unexpected character `{ch}`")));
                }
            }
        }
        Ok(())
    }
}

/// Parse a chemical formula into per-element atom counts.
///
/// Returns `(Z, count)` pairs sorted by atomic number with duplicate
/// elements aggregated (so `"CH3(CH2)6CH3"` yields `[(6, 4.0), (1, 18.0)]`).
///
/// Supported grammar: nested parenthesized groups with multi-digit group
/// counts, one/two-letter case-sensitive element symbols, multi-digit atom
/// counts, and trailing hydrate parts joined by `·` (or plain `.`), each
/// optionally prefixed by an integer multiplier (`CuSO4·5H2O`). Surrounding
/// whitespace is tolerated; embedded whitespace, a leading digit (bare
/// stoichiometric coefficients are not formulas), stray parentheses, and
/// unrecognized symbols are errors.
pub fn parse_formula(formula: &str) -> FormulaResult<Vec<(u32, f64)>> {
    let trimmed = formula.trim();
    let mut p = Parser {
        src: trimmed.as_bytes(),
        pos: 0,
    };
    let mut acc = BTreeMap::new();
    let mut first = true;
    while p.pos < trimmed.len() {
        if !first {
            if !p.at_hyphen_dot() {
                if p.peek() == Some(b')') {
                    return Err(p.error_here("unbalanced parenthesis: unexpected `)`"));
                }
                return Err(p.error_here("expected `.` or `·` hydrate separator"));
            }
            p.advance_over_separator();
        }
        if p.peek().is_some_and(|c| c.is_ascii_digit()) {
            if first {
                return Err(p.parse_error(p.pos, "formula must not begin with a digit"));
            }
            let mult = p.take_count()?.unwrap_or(1.0);
            p.take_units(&mut acc, mult)?;
        } else {
            p.take_units(&mut acc, 1.0)?;
        }
        first = false;
    }
    if first {
        return Err(FormulaError::ParseError {
            pos: 0,
            message: "empty formula".to_string(),
        });
    }
    Ok(acc.into_iter().collect())
}

// ---------------------------------------------------------------------------
// Abundance providers
// ---------------------------------------------------------------------------

/// Source of natural-isotope composition data, element by element.
///
/// Mirrors [`crate::MassProvider`]: expansion code is generic over this so
/// it never depends on the nuclear-data tables being present.
pub trait AbundanceProvider {
    /// Naturally occurring isotopes of element `z` with their abundance
    /// fractions, or `None` when nothing is tabulated for that element.
    fn natural_isotopes(&self, z: u32) -> Option<Vec<(NuclideId, f64)>>;
}

/// An [`AbundanceProvider`] that knows no isotopes.
///
/// Every lookup returns `None`, so expansions fail explicitly with
/// [`FormulaError::NoAbundanceData`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoAbundances;

impl AbundanceProvider for NoAbundances {
    fn natural_isotopes(&self, _z: u32) -> Option<Vec<(NuclideId, f64)>> {
        None
    }
}

/// [`AbundanceProvider`] backed by the natural-abundance table in
/// `nuclei::data`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NaturalAbundances;

impl NaturalAbundances {
    /// Abundance table grouped by atomic number (built once, lazily).
    fn groups() -> &'static BTreeMap<u32, Vec<(NuclideId, f64)>> {
        static GROUPS: OnceLock<BTreeMap<u32, Vec<(NuclideId, f64)>>> = OnceLock::new();
        GROUPS.get_or_init(|| {
            let mut groups: BTreeMap<u32, Vec<(NuclideId, f64)>> = BTreeMap::new();
            for (&nucid, &frac) in nuclei::data::abundance_table() {
                if frac > 0.0 {
                    let id = NuclideId::from_nucid(nucid);
                    groups.entry(id.z()).or_default().push((id, frac));
                }
            }
            groups
        })
    }
}

impl AbundanceProvider for NaturalAbundances {
    fn natural_isotopes(&self, z: u32) -> Option<Vec<(NuclideId, f64)>> {
        Self::groups().get(&z).filter(|v| !v.is_empty()).cloned()
    }
}

// ---------------------------------------------------------------------------
// Expansion helpers
// ---------------------------------------------------------------------------

/// True for natural-element placeholder ids: `z * 10_000_000`, i.e. the
/// nucid of zaid `z*1000` (`AAA == 0`) used throughout `mcnp-io` for bare
/// elemental zaids.
fn is_elemental(id: NuclideId) -> bool {
    id.a() == 0 && id.state() == 0
}

/// Atomic mass of `id`, falling back to the ground state for metastable ids
/// (the AME2020 table stores one mass per (Z, A); ground state and isomers
/// share it).
fn ground_mass(masses: &impl crate::MassProvider, id: NuclideId) -> Option<f64> {
    masses
        .mass(id.nucid())
        .or_else(|| masses.mass(id.nucid() - id.state()))
}

impl Material {
    /// Build a material from a chemical formula, expanding each element into
    /// its naturally occurring isotopes.
    ///
    /// Mirrors [`Material::from_atom_frac`]: atom counts come from the
    /// parsed formula weighted by natural-abundance fractions, stored masses
    /// are `n_i * M_i` using `masses`, and `density` is attached unchanged.
    /// Fails with the parse/abundance variants of [`FormulaError`] for bad
    /// input and with [`FormulaError::Core`] wrapping
    /// [`crate::Error::MissingMass`] when an isotope's mass is unknown.
    pub fn from_formula(
        formula: &str,
        masses: &impl crate::MassProvider,
        abundances: &impl AbundanceProvider,
        density: Option<f64>,
    ) -> FormulaResult<Self> {
        let elements = parse_formula(formula)?;
        let mut atoms = Vec::new();
        for &(z, count) in &elements {
            let isotopes = abundances
                .natural_isotopes(z)
                .ok_or(FormulaError::NoAbundanceData(z))?;
            let total: f64 = isotopes.iter().map(|(_, x)| x).sum();
            if total <= 0.0 {
                return Err(FormulaError::NoAbundanceData(z));
            }
            for (id, frac) in isotopes {
                atoms.push((id, count * frac / total));
            }
        }
        Ok(Material::from_atom_frac(&atoms, masses, density)?)
    }

    /// Replace natural-element placeholder entries with their isotopic
    /// breakdown, preserving each entry's stored mass.
    ///
    /// Placeholders follow the `mcnp-io` inp convention: a bare elemental
    /// zaid (`z*1000`, `AAA == 0`) becomes the nucid `z * 10_000_000`
    /// ([`is_elemental`]). Each placeholder of element `z` holding `g`
    /// grams is replaced by isotope masses `g * x_i * M_i / M̄`, where `x_i`
    /// are the (normalized) natural-abundance fractions and `M̄` the
    /// abundance-weighted mean atomic mass — i.e. the same number of atoms
    /// of each isotope as the elemental entry implied. Explicitly named
    /// nuclides are left untouched, so mixed elemental + isotopic
    /// compositions are supported.
    ///
    /// Fails with [`FormulaError::NoAbundanceData`] when the provider has no
    /// isotopes for an element, or [`FormulaError::Core`] wrapping
    /// [`crate::Error::MissingMass`] when an isotope mass is unknown.
    pub fn expand_elements(
        &mut self,
        masses: &impl crate::MassProvider,
        abundances: &impl AbundanceProvider,
    ) -> FormulaResult<()> {
        let mut expanded = BTreeMap::new();
        for (&id, &grams) in &self.comp {
            if !is_elemental(id) {
                expanded.insert(id, grams);
                continue;
            }
            let z = id.z();
            let isotopes = abundances
                .natural_isotopes(z)
                .ok_or(FormulaError::NoAbundanceData(z))?;
            let total: f64 = isotopes.iter().map(|(_, x)| x).sum();
            if total <= 0.0 {
                return Err(FormulaError::NoAbundanceData(z));
            }
            // Mean atomic mass of the natural element.
            let mean_mass = isotopes
                .iter()
                .map(|&(iso, x)| {
                    ground_mass(masses, iso)
                        .ok_or(crate::Error::MissingMass(iso))
                        .map(|m| x / total * m)
                })
                .sum::<crate::Result<f64>>()
                .map_err(FormulaError::from)?;
            for (iso, x) in isotopes {
                let m = ground_mass(masses, iso).expect("checked by mean_mass loop above");
                expanded.insert(iso, grams * (x / total) * m / mean_mass);
            }
        }
        self.comp = expanded;
        Ok(())
    }

    /// Inverse grouping of [`Material::expand_elements`]: fold every nuclide
    /// into its element's placeholder row keyed by the natural-element id
    /// `z * 10_000_000` (zaid `z*1000`). Placeholder entries already carry
    /// that key and simply accumulate alongside collapsed nuclides. Density
    /// and metadata are preserved; masses are summed exactly.
    pub fn collapse_elements(&self) -> Self {
        let mut comp = BTreeMap::new();
        for (&id, &grams) in &self.comp {
            let key = if is_elemental(id) {
                id
            } else {
                NuclideId::from_nucid(id.z() * 10_000_000)
            };
            *comp.entry(key).or_insert(0.0) += grams;
        }
        let mut out = Material::new();
        out.comp = comp;
        out.set_density(self.density());
        out.set_metadata(self.metadata().cloned());
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Ame2020;

    fn counts(formula: &str) -> Vec<(u32, f64)> {
        parse_formula(formula).unwrap()
    }

    fn assert_counts(formula: &str, expected: &[(u32, f64)]) {
        assert_eq!(counts(formula), expected.to_vec(), "for `{formula}`");
    }

    #[test]
    fn parses_water_and_glucose() {
        assert_counts("H2O", &[(1, 2.0), (8, 1.0)]);
        assert_counts("C6H12O6", &[(1, 12.0), (6, 6.0), (8, 6.0)]);
    }

    #[test]
    fn parses_grouped_and_nested_formulas() {
        // Ca(OH)2 → Ca1 O2 H2
        assert_counts("Ca(OH)2", &[(1, 2.0), (8, 2.0), (20, 1.0)]);
        // Fe2(SO4)3 → Fe2 S3 O12
        assert_counts("Fe2(SO4)3", &[(8, 12.0), (16, 3.0), (26, 2.0)]);
        // Nested groups.
        assert_counts("Mg(NO2)2", &[(7, 2.0), (8, 4.0), (12, 1.0)]);
        assert_counts("U((C)3)2", &[(6, 6.0), (92, 1.0)]);
    }

    #[test]
    fn parses_chained_groups_as_single_elements() {
        // CH3(CH2)6CH3 → C8H18 (n-octane written chain-style)
        assert_counts("CH3(CH2)6CH3", &[(1, 18.0), (6, 8.0)]);
    }

    #[test]
    fn parses_multi_digit_counts_and_hydrates() {
        assert_counts("C12H22O11", &[(1, 22.0), (6, 12.0), (8, 11.0)]);
        let dot = counts("CuSO4·5H2O");
        let ascii = counts("CuSO4.5H2O");
        assert_eq!(dot, ascii);
        assert_eq!(dot, vec![(1, 10.0), (8, 9.0), (16, 1.0), (29, 1.0)]);
        // Multiplier without parentheses after the dot.
        assert_counts("H2O.2H2O", &[(1, 6.0), (8, 3.0)]);
    }

    #[test]
    fn tolerates_surrounding_whitespace_only() {
        assert_eq!(counts("  H2O "), counts("H2O"));
    }

    #[test]
    fn rejects_unbalanced_parentheses() {
        let err = parse_formula("(H2O").unwrap_err();
        assert!(
            matches!(err, FormulaError::ParseError { pos: 4, .. }),
            "{err}"
        );
        let err = parse_formula("H2O)").unwrap_err();
        assert!(matches!(err, FormulaError::ParseError { .. }), "{err}");
    }

    #[test]
    fn rejects_unknown_symbol() {
        match parse_formula("Xx2O").unwrap_err() {
            FormulaError::UnknownElement(s) => assert_eq!(s, "Xx"),
            other => panic!("{other:?}"),
        }
        match parse_formula("Q").unwrap_err() {
            FormulaError::UnknownElement(s) => assert_eq!(s, "Q"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn rejects_leading_digit_empty_and_stray_chars() {
        let err = parse_formula("2H2O").unwrap_err();
        assert!(
            matches!(err, FormulaError::ParseError { pos: 0, .. }),
            "{err}"
        );
        assert!(matches!(
            parse_formula("").unwrap_err(),
            FormulaError::ParseError { .. }
        ));
        let err = parse_formula("H2 O").unwrap_err();
        assert!(
            matches!(err, FormulaError::ParseError { pos: 2, .. }),
            "{err}"
        );
        let err = parse_formula("H1O-1").unwrap_err();
        assert!(matches!(err, FormulaError::ParseError { .. }));
    }

    #[test]
    fn from_formula_water_has_natural_isotopes_and_two_thirds_hydrogen() {
        let water = Material::from_formula("H2O", &Ame2020, &NaturalAbundances, Some(1.0)).unwrap();
        // Natural hydrogen is pure H-1; oxygen carries all three isotopes.
        assert!(water
            .comp
            .contains_key(&NuclideId::from_name("H1").unwrap()));
        for o in ["O16", "O17", "O18"] {
            assert!(
                water.comp.contains_key(&NuclideId::from_name(o).unwrap()),
                "missing {o}"
            );
        }
        // Atom fractions: hydrogen contributes exactly 2 mol per mol water.
        let af = water.atom_fractions(&Ame2020).unwrap();
        let h: f64 = af
            .iter()
            .filter(|(id, _)| id.z() == 1)
            .map(|(_, f)| f)
            .sum();
        assert!((h - 2.0 / 3.0).abs() < 1e-12, "{h}");
        let o: f64 = af
            .iter()
            .filter(|(id, _)| id.z() == 8)
            .map(|(_, f)| f)
            .sum();
        assert!((o - 1.0 / 3.0).abs() < 1e-12);
        // Density passes through untouched.
        assert_eq!(water.density(), Some(1.0));
    }

    #[test]
    fn from_formula_rejects_bad_input_and_missing_data() {
        match Material::from_formula("Xx", &Ame2020, &NaturalAbundances, None).unwrap_err() {
            FormulaError::UnknownElement(s) => assert_eq!(s, "Xx"),
            other => panic!("{other:?}"),
        }
        assert!(matches!(
            Material::from_formula("U", &Ame2020, &NoAbundances, None).unwrap_err(),
            FormulaError::NoAbundanceData(92)
        ));
    }

    #[test]
    fn expand_then_collapse_round_trips_an_elemental_material() {
        // Elemental rows via mcnp-io placeholder ids.
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_nucid(10_000_000), 2.0); // H
        mat.add_nuclide(NuclideId::from_nucid(80_000_000), 16.0); // O
        let original = mat.clone();

        mat.expand_elements(&Ame2020, &NaturalAbundances).unwrap();
        // Isotopes appeared and placeholders vanished.
        assert!(!mat.comp.contains_key(&NuclideId::from_nucid(80_000_000)));
        assert!(mat.comp.contains_key(&NuclideId::from_name("H1").unwrap()));
        assert!(mat.comp.contains_key(&NuclideId::from_name("H2").unwrap()));
        assert!(mat.comp.contains_key(&NuclideId::from_name("O18").unwrap()));

        let back = mat.collapse_elements();
        assert_eq!(
            back.comp.keys().copied().collect::<Vec<_>>(),
            original.comp.keys().copied().collect::<Vec<_>>()
        );
        for (id, m0) in &original.comp {
            let m1 = back.comp[id];
            assert!((m0 - m1).abs() < 1e-9 * m0.abs(), "{id}: {m0} vs {m1}");
        }
    }

    #[test]
    fn expand_preserves_entry_masses_and_leaves_named_nuclides_alone() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_nucid(10_000_000), 18.0); // H, 18 g
        mat.add_nuclide(NuclideId::from_name("Fe56").unwrap(), 5.0);

        mat.expand_elements(&Ame2020, &NaturalAbundances).unwrap();
        close(mat.mass(), 23.0);
        close(
            mat.remove_nuclide(NuclideId::from_name("Fe56").unwrap())
                .unwrap(),
            5.0,
        );
        let h: f64 = mat.comp.values().sum();
        close(h, 18.0);
    }

    #[test]
    fn expand_without_abundances_errors_with_z() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_nucid(920_000_000), 1.0);
        match mat.expand_elements(&Ame2020, &NoAbundances).unwrap_err() {
            FormulaError::NoAbundanceData(z) => assert_eq!(z, 92),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn collapse_folds_named_nuclides_into_placeholder_keys() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 3.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 1.0);
        mat.set_density(Some(19.1));
        let collapsed = mat.collapse_elements();

        let key = NuclideId::from_nucid(920_000_000);
        assert_eq!(collapsed.comp.len(), 1);
        close(collapsed.comp[&key], 4.0);
        // The key really is the z*10_000_000 placeholder form (zaid 92000).
        assert_eq!(key.nucid(), nuclei::element_z("U").unwrap() * 10_000_000);
        assert_eq!(collapsed.density(), Some(19.1));
    }

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }
}
