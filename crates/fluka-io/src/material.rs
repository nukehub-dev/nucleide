//! FLUKA input-deck material cards: built-in material table plus
//! `MATERIAL` / `COMPOUND` card-string generators.
//!
//! Reuses the FLUKA element/isotope name table vendored in
//! [`nuclei::dialects`]; golden card strings are pinned by the vendored
//! fixtures.
//!
//! # Formatting contract (C++ iostream parity)
//!
//! Cards are sequences of 10-character fields. Numeric fields use two
//! upstream styles, reproduced exactly by [`format_field`]:
//!
//! - integral values print as fixed precision-0 with `showpoint`, i.e.
//!   trailing dot (`999.`), right-justified;
//! - fractional values print in C++ default general format at precision 6
//!   significant digits (trailing zeros stripped, scientific below 1e-4
//!   or from 1e6 up), right-justified — e.g. `235.0439…` → `235.044`.
//!
//! Compound fractions always use scientific notation with three decimals
//! (`4.000e-02`). Mass fractions carry a leading `-` sign, atom fractions
//! none (`FracType`). Upstream quirks preserved: density is passed through
//! `sqrt(d²)` (negative densities print as positive), and a compound's
//! own MATERIAL line always uses `z = 1`, `mass = 1`.
//!
//! These are free functions over explicit arguments; comment lines from
//! material metadata are therefore not emitted.

use std::collections::BTreeSet;
use std::fmt;

use nuclei::data as nuc_data;
use nuclei::{dialects, NuclideId};

/// Errors raised while generating FLUKA cards.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Nuclide → FLUKA-name resolution failed (see [`nuclei::dialects`]).
    Naming(dialects::DialectError),
    /// Atomic number has no entry in the FLUKA element-name table.
    UnknownElementZ(u32),
    /// Atomic mass unavailable for this nuclide or natural element.
    NoAtomicMass(String),
    /// A compound card was requested with no components.
    EmptyCompound,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Naming(e) => write!(f, "fluka name lookup failed: {e}"),
            Error::UnknownElementZ(z) => {
                write!(f, "no FLUKA element name for atomic number {z}")
            }
            Error::NoAtomicMass(name) => write!(f, "atomic mass unavailable for `{name}`"),
            Error::EmptyCompound => write!(f, "compound card requested with no components"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Naming(e) => Some(e),
            _ => None,
        }
    }
}

/// The 37 predefined FLUKA materials; anything else needs a MATERIAL card.
///
/// Built-in FLUKA material names (names are
/// truncated to FLUKA's eight-character limit).
pub const BUILTIN_MATERIALS: [&str; 37] = [
    "BLCKHOLE", "VACUUM", "HYDROGEN", "HELIUM", "BERYLLIU", "CARBON", "NITROGEN", "OXYGEN",
    "MAGNESIU", "ALUMINUM", "IRON", "COPPER", "SILVER", "SILICON", "GOLD", "MERCURY", "LEAD",
    "TANTALUM", "SODIUM", "ARGON", "CALCIUM", "TIN", "TUNGSTEN", "TITANIUM", "NICKEL", "WATER",
    "POLYSTYR", "PLASCINT", "PMMA", "BONECOMP", "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE",
    "KAPTON", "POLYETHY", "AIR",
];

/// Whether `name` is one of FLUKA's built-in materials (no card needed).
///
/// Whether the name denotes a built-in FLUKA material (inverse predicate).
pub fn is_fluka_builtin(name: &str) -> bool {
    BUILTIN_MATERIALS.contains(&name)
}

/// FLUKA element names for natural elements, vendored verbatim from the
/// `A == 0` rows of the FLUKA translation table (including rows
/// upstream flags "// No fluka"; they resolve there too). Sorted by Z.
const ELEMENT_NAMES: &[(u32, &str)] = &[
    (1, "HYDROGEN"),
    (2, "HELIUM"),
    (3, "LITHIUM"),
    (4, "BERYLLIU"),
    (5, "BORON"),
    (6, "CARBON"),
    (7, "NITROGEN"),
    (8, "OXYGEN"),
    (9, "FLUORINE"),
    (10, "NEON"),
    (11, "SODIUM"),
    (12, "MAGNESIU"),
    (13, "ALUMINUM"),
    (14, "SILICON"),
    (15, "PHOSPHO"),
    (16, "SULFUR"),
    (17, "CHLORINE"),
    (18, "ARGON"),
    (19, "POTASSIU"),
    (20, "CALCIUM"),
    (21, "SCANDIUM"),
    (22, "TITANIUM"),
    (23, "VANADIUM"),
    (24, "CHROMIUM"),
    (25, "MANGANES"),
    (26, "IRON"),
    (27, "COBALT"),
    (28, "NICKEL"),
    (29, "COPPER"),
    (30, "ZINC"),
    (31, "GALLIUM"),
    (32, "GERMANIU"),
    (33, "ARSENIC"),
    (34, "SELENIUM"),
    (35, "BROMINE"),
    (36, "KRYPTON"),
    (37, "RUBIDIUM"),
    (38, "STRONTIU"),
    (39, "YTTRIUM"),
    (40, "ZIRCONIU"),
    (41, "NIOBIUM"),
    (42, "MOLYBDEN"),
    (43, "99-TC"),
    (44, "RUTHENIU"),
    (45, "RHODIUM"),
    (46, "PALLADIU"),
    (47, "SILVER"),
    (48, "CADMIUM"),
    (49, "INDIUM"),
    (50, "TIN"),
    (51, "ANTIMONY"),
    (52, "TELLURIU"),
    (53, "IODINE"),
    (54, "XENON"),
    (55, "CESIUM"),
    (56, "BARIUM"),
    (57, "LANTHANU"),
    (58, "CERIUM"),
    (59, "PRASEODY"),
    (60, "NEODYMIU"),
    (61, "PROMETHI"),
    (62, "SAMARIUM"),
    (63, "EUROPIUM"),
    (64, "GADOLINI"),
    (65, "TERBIUM"),
    (66, "DYSPROSI"),
    (67, "HOLMIUM"),
    (68, "ERBIUM"),
    (69, "THULIUM"),
    (70, "YTTERBIU"),
    (71, "LUTETIUM"),
    (72, "HAFNIUM"),
    (73, "TANTALUM"),
    (74, "TUNGSTEN"),
    (75, "RHENIUM"),
    (76, "OSMIUM"),
    (77, "IRIDIUM"),
    (78, "PLATINUM"),
    (79, "GOLD"),
    (80, "MERCURY"),
    (81, "THALLIUM"),
    (82, "LEAD"),
    (83, "BISMUTH"),
    (84, "POLONIUM"),
    (85, "ASTATINE"),
    (86, "RADON"),
    (87, "FRANCIUM"),
    (88, "RADIUM"),
    (89, "ACTINIUM"),
    (90, "THORIUM"),
    (91, "PROTACTI"),
    (92, "URANIUM"),
    (93, "NEPTUNIU"),
    (94, "239-PU"),
    (95, "241-AM"),
    (96, "CURIUM"),
    (97, "BERKELIU"),
    (98, "CALIFORN"),
    (99, "EINSTEIN"),
    (100, "FERMIUM"),
    (101, "MENDELEV"),
    (102, "NOBELIUM"),
    (103, "LAWRENCI"),
    (104, "RUTHERFO"),
    (105, "DUBNIUM"),
    (106, "SEABORGI"),
    (107, "BOHRIUM"),
    (108, "HASSIUM"),
    (109, "MEITNERI"),
    (110, "DARMSTAD"),
    (111, "ROENTGEN"),
    (112, "COPERNIC"),
    (114, "UNUNQUAD"),
    (116, "UNUNHEXI"),
];

/// FLUKA name for the natural element with atomic number `z`.
pub fn fluka_element_name(z: u32) -> Option<&'static str> {
    ELEMENT_NAMES
        .iter()
        .find(|(known_z, _)| *known_z == z)
        .map(|(_, name)| *name)
}

/// Which nucleus a card field references.
///
/// FLUKA compounds routinely cite natural elements (`CARBON`, `OXYGEN`),
/// which have no isotope mass and thus cannot be represented by
/// [`NuclideId`]; those use the [`FlukaNuc::Element`] variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlukaNuc {
    /// A concrete nuclide (isotope or metastable).
    Nuclide(NuclideId),
    /// A natural element by atomic number.
    Element(u32),
}

impl From<NuclideId> for FlukaNuc {
    fn from(nuc: NuclideId) -> Self {
        FlukaNuc::Nuclide(nuc)
    }
}

impl From<u32> for FlukaNuc {
    fn from(z: u32) -> Self {
        FlukaNuc::Element(z)
    }
}

impl FlukaNuc {
    /// Resolve to the FLUKA table name for this nucleus.
    pub fn fluka_name(&self) -> Result<&'static str, Error> {
        match self {
            FlukaNuc::Nuclide(nuc) => dialects::id_to_fluka(*nuc).map_err(Error::Naming),
            FlukaNuc::Element(z) => fluka_element_name(*z).ok_or(Error::UnknownElementZ(*z)),
        }
    }

    /// Atomic number of the referenced nucleus.
    pub fn z(&self) -> u32 {
        match self {
            FlukaNuc::Nuclide(nuc) => nuc.z(),
            FlukaNuc::Element(z) => *z,
        }
    }

    /// Sort key: ascending nucid per component
    /// ordering: natural elements sort as their bare `zzz0000000` id,
    /// which is strictly below every isotope of the same element.
    fn sort_key(&self) -> u32 {
        match self {
            FlukaNuc::Nuclide(nuc) => nuc.nucid(),
            FlukaNuc::Element(z) => z * 10_000_000,
        }
    }

    /// Atomic mass in u: AME2020 value for nuclides (ground-state fallback
    /// for metastable ids), abundance-weighted mean for natural elements.
    pub fn atomic_mass(&self) -> Option<f64> {
        match self {
            FlukaNuc::Nuclide(nuc) => nuc_data::atomic_mass(nuc.nucid())
                .or_else(|| nuc_data::atomic_mass(nuc.nucid() - nuc.state())),
            FlukaNuc::Element(z) => natural_mean_atomic_mass(*z),
        }
    }
}

/// Abundance-weighted mean atomic mass of the natural element `z` in u.
pub fn natural_mean_atomic_mass(z: u32) -> Option<f64> {
    let mut weighted = 0.0;
    let mut total_abundance = 0.0;
    for (&nucid, &abundance) in nuc_data::abundance_table() {
        if NuclideId::from_nucid(nucid).z() == z {
            weighted += abundance * nuc_data::atomic_mass(nucid)?;
            total_abundance += abundance;
        }
    }
    (total_abundance > 0.0).then(|| weighted / total_abundance)
}

/// One component line-entry of a COMPOUND card.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Component {
    /// Referenced nucleus.
    pub nuc: FlukaNuc,
    /// Fraction value (sign convention applied at render time).
    pub fraction: f64,
}

impl Component {
    /// Build a component from any [`FlukaNuc`] conversion.
    pub fn new(nuc: impl Into<FlukaNuc>, fraction: f64) -> Self {
        Component {
            nuc: nuc.into(),
            fraction,
        }
    }
}

/// Fraction convention for COMPOUND cards.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FracType {
    /// Mass fractions render with a `-` sign prefix (FLUKA default).
    Mass,
    /// Atom fractions render without a sign prefix.
    Atom,
}

/// Format one numeric field following the FLUKA card convention:
/// right-justified in ten characters, trailing-dot style for integers,
/// otherwise six significant digits with trailing zeros stripped.
pub fn format_field(field: f64) -> String {
    let body = if field == field.trunc() {
        // C++ precision(0) + fixed + showpoint: integer with trailing '.'.
        format!("{field}.")
    } else {
        general_precision6(field)
    };
    format!("{body:>10}")
}

/// C++ default floatfield at precision 6 (i.e. `%g` semantics).
fn general_precision6(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let exponent = value.abs().log10().floor();
    if !(-4.0..6.0).contains(&exponent) {
        return scientific(value, 5);
    }
    let decimals = (5.0 - exponent).max(0.0) as usize;
    let mut text = format!("{value:.decimals$}");
    while text.ends_with('0') || text.ends_with('.') {
        text.pop();
    }
    text
}

/// Scientific notation `d.ddde±XX` with `digits` mantissa decimals and a
/// zero-padded two-digit exponent, as C++ iostreams emit.
fn scientific(value: f64, digits: usize) -> String {
    let rendered = format!("{value:.digits$e}");
    let (mantissa, exponent) = rendered.split_once('e').expect("rust sci format");
    let (sign, magnitude) = match exponent.strip_prefix('-') {
        Some(mag) => ("-", mag),
        None => ("+", exponent),
    };
    // Explicit fill/align (`0>`); the bare `0` flag does not zero-fill
    // string operands in Rust.
    format!("{mantissa}e{sign}{magnitude:0>width$}", width = 2)
}

fn left_10(text: &str) -> String {
    format!("{text:<10}")
}

fn right_10(text: &str) -> String {
    format!("{text:>10}")
}

/// Render a full `MATERIAL` card line (plus newline).
///
/// Build a `MATERIAL` card line without metadata comments:
/// fields are `MATERIAL`, `znum`, `atomic_mass`, `density`, `fid`, two
/// blanks, then the FLUKA name (left-justified). Density passes through
/// the upstream `sqrt(d*d)` idiom, so sign is irrelevant.
pub fn material_line(
    znum: u32,
    atomic_mass: f64,
    density: f64,
    fid: u32,
    fluka_name: &str,
) -> String {
    let mut card = String::with_capacity(96);
    card.push_str(&left_10("MATERIAL"));
    card.push_str(&right_10(&format!("{znum}.")));
    card.push_str(&format_field(atomic_mass));
    card.push_str(&format_field((density * density).sqrt()));
    card.push_str(&right_10(&format!("{fid}.")));
    card.push_str(&right_10(""));
    card.push_str(&right_10(""));
    card.push_str(&left_10(fluka_name));
    card.push('\n');
    card
}

/// Render the MATERIAL record for an elemental nuclide, or an empty
/// string when it is already a built-in FLUKA material.
///
/// Build the `MATERIAL` card string for an element or isotope: the z number comes from
/// the nuclide/element, the atomic mass from [`FlukaNuc::atomic_mass`].
pub fn material_str(fid: u32, nuc: impl Into<FlukaNuc>, density: f64) -> Result<String, Error> {
    let nuc = nuc.into();
    let fluka_name = nuc.fluka_name()?;
    if is_fluka_builtin(fluka_name) {
        // Built-ins need no card; upstream returns empty here too.
        return Ok(String::new());
    }
    let atomic_mass = nuc
        .atomic_mass()
        .ok_or_else(|| Error::NoAtomicMass(fluka_name.to_string()))?;
    Ok(material_line(
        nuc.z(),
        atomic_mass,
        density,
        fid,
        fluka_name,
    ))
}

/// Render the MATERIAL + COMPOUND records describing a compound.
///
/// Build the two-line `COMPOUND` card string with the compound's name
/// and density passed explicitly (upstream reads them from metadata /
/// object state). The compound's own MATERIAL line always uses
/// `z = 1, mass = 1`. Components are emitted sorted by nucid-equivalent
/// key ascending, three per line, remainder on a blank-padded final line;
/// each line is exactly eighty columns.
pub fn compound_str(
    fid: u32,
    compound_name: &str,
    density: f64,
    frac_type: FracType,
    components: &[Component],
) -> Result<String, Error> {
    if components.is_empty() {
        return Err(Error::EmptyCompound);
    }

    // Resolve every name up front so errors never leave partial output.
    let mut sorted: Vec<&Component> = components.iter().collect();
    sorted.sort_by_key(|c| c.nuc.sort_key());
    let names: Vec<&'static str> = sorted
        .iter()
        .map(|c| c.nuc.fluka_name())
        .collect::<Result<_, _>>()?;

    let frac_sign = match frac_type {
        FracType::Mass => "-",
        FracType::Atom => "",
    };
    let frac_field = |fraction: f64| -> String {
        // Upstream streams `sign << value` through `scientific << setprecision(3)`.
        let signed = format!("{}{}", frac_sign, scientific(fraction, 3));
        right_10(&signed)
    };

    let mut card = material_line(1, 1.0, density, fid, compound_name);

    let mut triples = sorted.chunks_exact(3);
    for trio in &mut triples {
        card.push_str(&left_10("COMPOUND"));
        for (component, name) in trio.iter().zip(&names) {
            card.push_str(&frac_field(component.fraction));
            card.push_str(&right_10(name));
        }
        card.push_str(&left_10(compound_name));
        card.push('\n');
    }

    let remainder = triples.remainder();
    if !remainder.is_empty() {
        card.push_str(&left_10("COMPOUND"));
        let skip = sorted.len() - remainder.len();
        for (component, name) in remainder.iter().zip(names[skip..].iter()) {
            card.push_str(&frac_field(component.fraction));
            card.push_str(&right_10(name));
        }
        // Pad the unused component slots, then emit the two trailing
        // blanks upstream always writes before the compound name, keeping
        // every COMPOUND line at eighty columns.
        for _ in 0..(6 - 2 * remainder.len()) {
            card.push_str(&right_10(""));
        }
        card.push_str(&left_10(compound_name));
        card.push('\n');
    }

    Ok(card)
}

/// Set of built-in names, provided for callers that filter repeatedly
/// (mirrors the static `std::set` built in upstream C++).
pub fn builtin_set() -> BTreeSet<&'static str> {
    BUILTIN_MATERIALS.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const U235_NUCID: u32 = 922_350_000;

    #[test]
    fn builtin_table_matches_upstream() {
        assert_eq!(BUILTIN_MATERIALS.len(), 37);
        for name in [
            "BLCKHOLE", "VACUUM", "WATER", "AIR", "LEAD", "PMMA", "POLYETHY",
        ] {
            assert!(is_fluka_builtin(name), "{name} should be builtin");
        }
        for name in ["HYDROG-1", "235-U", "URANIUM", "ORGPOLYM", "", "water"] {
            assert!(!is_fluka_builtin(name), "`{name}` should not be builtin");
        }
        assert_eq!(builtin_set().len(), 37);
    }

    #[test]
    fn element_names_cover_and_resolve() {
        assert_eq!(fluka_element_name(92), Some("URANIUM"));
        assert_eq!(fluka_element_name(8), Some("OXYGEN"));
        assert_eq!(fluka_element_name(26), Some("IRON"));
        assert_eq!(fluka_element_name(94), Some("239-PU"));
        assert_eq!(fluka_element_name(113), None);
        assert_eq!(fluka_element_name(119), None);
    }

    #[test]
    fn fluka_nuc_names_wire_through_dialects() {
        let u235 = NuclideId::from_nucid(U235_NUCID);
        assert_eq!(
            FlukaNuc::from(u235).fluka_name().unwrap(),
            nuclei::dialects::id_to_fluka(u235).unwrap()
        );
        assert_eq!(FlukaNuc::from(u235).fluka_name().unwrap(), "235-U");
        let h1 = NuclideId::from_name("H1").unwrap();
        assert_eq!(FlukaNuc::from(h1).fluka_name().unwrap(), "HYDROG-1");
        // Isotopes absent from the FLUKA table error out.
        let o16 = NuclideId::from_name("O16").unwrap();
        assert!(matches!(
            FlukaNuc::from(o16).fluka_name(),
            Err(Error::Naming(_))
        ));
        assert_eq!(FlukaNuc::from(82_u32).fluka_name().unwrap(), "LEAD");
        assert!(matches!(
            FlukaNuc::from(200_u32).fluka_name(),
            Err(Error::UnknownElementZ(200))
        ));
    }

    #[test]
    fn format_field_integer_branch_trailing_dot() {
        assert_eq!(format_field(999.0), "      999.");
        // Oracle geometry: "1." occupies the last two columns of the
        // ten-character field ("MATERIAL         92.   ...").
        assert_eq!(format_field(1.0), "        1.");
        assert_eq!(format_field(25.0), "       25.");
        assert_eq!(format_field(92.0), "       92.");
        assert_eq!(format_field(-5.0), "       -5.");
    }

    #[test]
    fn format_field_general_six_digits_stripped() {
        // Docstring examples from Material::fluka_format_field (content is
        // right-justified in a ten-character field):
        assert_eq!(format_field(999.12), "    999.12");
        assert_eq!(format_field(999.123), "   999.123");
        assert_eq!(format_field(999.1234), "   999.123");
        // Oracle-derived values:
        assert_eq!(format_field(235.043_928_117), "   235.044");
        assert_eq!(format_field(238.050_788_26), "   238.051");
        assert_eq!(format_field(19.1), "      19.1");
        assert_eq!(format_field(1.007_825_031_898), "   1.00783");
        assert_eq!(format_field(11.35), "     11.35");
    }

    #[test]
    fn scientific_fields_match_oracle_fractions() {
        let atom = |v: f64| right_10(&scientific(v, 3));
        assert_eq!(atom(0.04), " 4.000e-02");
        assert_eq!(atom(0.96), " 9.600e-01");
        assert_eq!(atom(0.1), " 1.000e-01");
        assert_eq!(atom(0.19), " 1.900e-01");
        // Mass type prepends '-' outside the scientific rendering.
        assert_eq!(format!("-{}", scientific(0.04, 3)), "-4.000e-02");
    }

    #[test]
    fn material_line_golden_from_oracle() {
        // Oracle Part I element line (density defaults to ±1 → sqrt trick).
        assert_eq!(
            material_line(92, 235.043_928_117, 1.0, 25, "235-U"),
            "MATERIAL         92.   235.044        1.       25.                    235-U     \n"
        );
        // Oracle collapsed-material line with real density.
        assert_eq!(
            material_line(92, 238.028_91, 19.1, 25, "URANIUM"),
            "MATERIAL         92.   238.029      19.1       25.                    URANIUM   \n"
        );
    }

    #[test]
    fn material_str_empty_for_builtins() {
        assert_eq!(material_str(35, FlukaNuc::from(82_u32), 11.35).unwrap(), "");
        assert_eq!(material_str(1, FlukaNuc::from(79_u32), 19.32).unwrap(), "");
        assert_eq!(material_str(1, FlukaNuc::from(8_u32), 1.0).unwrap(), "");
    }

    #[test]
    fn material_str_isotope_card_golden() {
        let h1 = NuclideId::from_name("H1").unwrap();
        assert_eq!(
            material_str(25, h1, 1.0).unwrap(),
            "MATERIAL          1.   1.00783        1.       25.                    HYDROG-1  \n"
        );
        // Negative density prints positive via the sqrt(d*d) idiom.
        assert_eq!(
            material_str(25, h1, -1.0).unwrap(),
            "MATERIAL          1.   1.00783        1.       25.                    HYDROG-1  \n"
        );
    }

    #[test]
    fn material_str_natural_element_mean_mass() {
        let mean = natural_mean_atomic_mass(92).unwrap();
        assert!(
            (mean - 238.028_91).abs() < 1e-3,
            "natural U mean mass {mean} should be ~238.029"
        );
        // Oracle collapsed-material MATERIAL line, byte for byte.
        assert_eq!(
            material_str(25, FlukaNuc::from(92_u32), 19.1).unwrap(),
            "MATERIAL         92.   238.029      19.1       25.                    URANIUM   \n"
        );
        // Synthetic Z without abundances has no mean mass.
        assert_eq!(natural_mean_atomic_mass(118), None);
    }

    #[test]
    fn compound_str_leu_atom_fraction_oracle() {
        let u235 = NuclideId::from_nucid(U235_NUCID);
        let u238 = NuclideId::from_name("U238").unwrap();
        let card = compound_str(
            27,
            "URANIUM",
            19.1,
            FracType::Atom,
            &[Component::new(u235, 0.04), Component::new(u238, 0.96)],
        )
        .unwrap();
        let expected = concat!(
            "MATERIAL          1.        1.      19.1       27.                    URANIUM   \n",
            "COMPOUND   4.000e-02     235-U 9.600e-01     238-U                    URANIUM   \n",
        );
        assert_eq!(card, expected);
    }

    #[test]
    fn compound_str_leu_mass_fraction_oracle() {
        let u235 = NuclideId::from_nucid(U235_NUCID);
        let u238 = NuclideId::from_name("U238").unwrap();
        let card = compound_str(
            27,
            "URANIUM",
            19.1,
            FracType::Mass,
            &[Component::new(u235, 0.04), Component::new(u238, 0.96)],
        )
        .unwrap();
        let expected = concat!(
            "MATERIAL          1.        1.      19.1       27.                    URANIUM   \n",
            "COMPOUND  -4.000e-02     235-U-9.600e-01     238-U                    URANIUM   \n",
        );
        assert_eq!(card, expected);
    }

    #[test]
    fn compound_str_three_components_full_line_oracle() {
        let card = compound_str(
            25,
            "ORGPOLYM",
            1.0,
            FracType::Mass,
            &[
                // Oracle uses natural H (→ "HYDROGEN"), not an isotope.
                Component::new(FlukaNuc::Element(1), 0.1),
                Component::new(FlukaNuc::Element(8), 0.8),
                Component::new(FlukaNuc::Element(6), 0.1),
            ],
        )
        .unwrap();
        let expected = concat!(
            "MATERIAL          1.        1.        1.       25.                    ORGPOLYM  \n",
            "COMPOUND  -1.000e-01  HYDROGEN-1.000e-01    CARBON-8.000e-01    OXYGENORGPOLYM  \n",
        );
        // Components arrive unsorted but must emit in ascending-nucid order.
        assert_eq!(card, expected);
    }

    #[test]
    fn compound_str_remainder_line_padding() {
        let card = compound_str(
            30,
            "MIXTURE",
            2.7,
            FracType::Atom,
            &[
                Component::new(FlukaNuc::Element(1), 0.1),
                Component::new(FlukaNuc::Element(13), 0.2),
                Component::new(FlukaNuc::Element(26), 0.3),
                Component::new(NuclideId::from_name("U238").unwrap(), 0.4),
            ],
        )
        .unwrap();
        let lines: Vec<&str> = card.lines().collect();
        assert_eq!(lines.len(), 3);
        // Full triple first…
        assert_eq!(
            lines[1],
            "COMPOUND   1.000e-01  HYDROGEN 2.000e-01  ALUMINUM 3.000e-01      IRONMIXTURE   "
        );
        // …then the single-component remainder padded to 80 columns.
        assert_eq!(lines[2].len(), 80);
        assert_eq!(
            lines[2],
            "COMPOUND   4.000e-01     238-U                                        MIXTURE   "
        );
    }

    #[test]
    fn compound_str_errors() {
        let u235 = NuclideId::from_nucid(U235_NUCID);
        assert_eq!(
            compound_str(1, "X", 1.0, FracType::Mass, &[]),
            Err(Error::EmptyCompound)
        );
        // O16 has no FLUKA table entry → naming error before any output.
        let o16 = NuclideId::from_name("O16").unwrap();
        assert!(matches!(
            compound_str(
                1,
                "X",
                1.0,
                FracType::Mass,
                &[Component::new(o16, 1.0), Component::new(u235, 1.0)]
            ),
            Err(Error::Naming(_))
        ));
        // Metastable ids fall back to the ground-state mass.
        let am242m = NuclideId::from_name("Am242_m1").unwrap();
        let mass = FlukaNuc::from(am242m).atomic_mass();
        assert_eq!(mass, nuc_data::atomic_mass(952_420_000));
    }
}
