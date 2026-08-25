//! MCNP input-deck material parsing.
//!
//! Validated against the vendored fixtures
//! `fixtures/mcnp/inp/{mcnp_inp,mcnp_inp_comments}.txt`.
//!
//! # Parsing conventions
//!
//! - **Cell densities**: a line is a cell card when it has > 3 whitespace
//!   tokens whose first two are digits, the third does not start with a
//!   letter (excludes surface cards like `99 7 PX 180`), the line does not
//!   begin with five spaces (continuation), and the second token is not
//!   `"0"` (void). Token 2 is the density; densities referencing the same
//!   material are kept when their relative difference is >= 1e-4.
//! - **Material cards**: first token matching `[mM]<digits>`; the remainder
//!   of the token is the material number. Commented-out cards (`C m4 ...`)
//!   never match because the first *token*, not character, is inspected.
//! - **Card body**: the `$`-truncated card line plus subsequent lines until
//!   a blank line, EOF, or a line that neither begins with five spaces nor
//!   with `c`/`C`. Among those continuations, lines whose first token starts
//!   with `c`/`C` are comments and contribute no data, so commented-out
//!   nuclide lines mid-card are supported.
//! - **Fractions**: `zaid[.suffix] <fraction>` pairs; library keywords
//!   (`NLIB`/`PLIB`/`HLIB`/`PNLIB`/`ELIB`, case-insensitive) consume one
//!   token and carry no composition. Repeated nuclides accumulate
//!   (`+=`, as upstream). The fraction type is the sign of the *first
//!   non-zero* fraction (positive → [`FracKind::Atom`], negative →
//!   [`FracKind::Mass`]); mixed-sign decks follow upstream's
//!   first-nonzero-wins rule instead of failing.
//! - **Comments metadata**: the contiguous block of `c`/`C` lines directly
//!   above the material card, walking upwards until a bare separator line
//!   (`c`/`C` alone), a blank line, or any other card. Stored top-down with
//!   the marker stripped; [`McnpMaterial::name`],
//!   [`McnpMaterial::source`] and [`McnpMaterial::comments_text`] recover
//!   the `name:`/`source:`/`comments:` keys upstream parses out of it.
//! - **Densities**: a single distinct cell density is exposed converted to
//!   g/cm³ (negative → absolute value = mass density; positive → atom
//!   density × 1e24 atoms/cm³ converted through the effective molar mass).
//!   Materials referenced with several distinct densities are upstream
//!   `MultiMaterial`s, which this model cannot represent: [`McnpMaterial::density`]
//!   is `None` there, mirroring the absent-density case.
//!
//! # Documented deviations from the legacy reader
//!
//! - Natural-element zaids (`AAA == 0`, e.g. `1000` H-nat, `8000` O-nat)
//!   cannot be held by a validated [`NuclideId`]; they are carried as
//!   unvalidated placeholders via `NuclideId::from_nucid(zaid * 10_000)`
//!   (`a() == 0`, `state() == 0`) — the same integer as the legacy
//!   natural-element id convention (`zz * 10**7`). Check `a() == 0` before
//!   treating an entry as a specific nuclide.
//! - `table_ids` (the `.15c` suffixes) and library assignments are parsed
//!   and discarded; the typed model has no metadata bag.
//! - Results are a `Vec` in file order rather than a dict keyed by material
//!   number (duplicate material numbers are preserved instead of silently
//!   overwritten).
//! - Upstream crash paths become clean behaviour: a blank line inside a
//!   comment block ends the block (upstream raises `IndexError`), a zero
//!   density compares equal only to itself (upstream divides by zero), and
//!   an all-zero fraction card is a [`Error::BadCard`] rather than an
//!   infinite loop.
//! - `foo=bar` tokens with unrecognized keys raise [`Error::UnknownKeyword`]
//!   where upstream silently ignores them; `comments_text` keeps colons
//!   intact (upstream's `split(":")` join mangles e.g. `http://` URLs).

use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;

use nuclei::data::{abundance_table, atomic_mass, mass_table};
use nuclei::{dialects, NuclideId};

/// Avogadro constant (2019 SI exact value), atoms per mole.
const AVOGADRO: f64 = 6.022_140_76e23;

/// Library-assignment keywords recognized on M cards.
const LIB_KEYWORDS: [&str; 5] = ["NLIB", "PLIB", "HLIB", "PNLIB", "ELIB"];

/// Metadata keys recognized in the comment block above a material card
/// (material densities).
const METADATA_KEYS: [&str; 3] = ["source", "comments", "name"];

/// Errors raised while parsing materials from an MCNP input deck.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Filesystem read failed.
    Io(String),
    /// A material card (or its neighbourhood) is malformed.
    BadCard {
        /// 1-based line number of the material card.
        line: usize,
        /// What was wrong.
        message: String,
    },
    /// A numeric field failed to parse.
    BadNumber {
        /// Which field was being parsed.
        context: &'static str,
        /// The offending text.
        text: String,
    },
    /// A `key=value` token carried an unrecognized key.
    UnknownKeyword {
        /// 1-based line number of the material card.
        line: usize,
        /// The rejected key.
        keyword: String,
    },
    /// A zaid admits no nuclide interpretation (wraps the dialect rules).
    BadZaid {
        /// 1-based line number of the material card.
        line: usize,
        /// The rejected zaid.
        zaid: u32,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadCard { line, message } => {
                write!(f, "bad material card on line {line}: {message}")
            }
            Error::BadNumber { context, text } => {
                write!(f, "cannot parse {context} from `{text}`")
            }
            Error::UnknownKeyword { line, keyword } => {
                write!(
                    f,
                    "unknown keyword `{keyword}` on material card line {line}"
                )
            }
            Error::BadZaid { line, zaid } => {
                write!(
                    f,
                    "uninterpretable zaid {zaid} on material card line {line}"
                )
            }
        }
    }
}

impl std::error::Error for Error {}

/// Whether an M card's fractions are atom or mass fractions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FracKind {
    /// Fractions are atom fractions (first non-zero fraction positive).
    Atom,
    /// Fractions are mass fractions (first non-zero fraction negative,
    /// MCNP prints these negated).
    Mass,
}

/// One MCNP material definition parsed from an input deck.
///
/// See the module docs for the conventions and deviations relative to
/// The reference material-line grammar.
#[derive(Debug, Clone, PartialEq)]
pub struct McnpMaterial {
    /// Material number from the `m<N>` card identifier.
    pub number: u32,
    /// Composition in card order; repeated zaids accumulate as upstream.
    /// Natural-element zaids (`AAA == 0`) appear as placeholder ids with
    /// [`NuclideId::a`] `== 0` (see module docs).
    pub fractions: Vec<(NuclideId, f64)>,
    /// Atom vs mass interpretation, from the sign of the first non-zero
    /// fraction.
    pub fraction_type: FracKind,
    /// Density in g/cm³ when the material is referenced by exactly one
    /// distinct cell density (`None` for unreferenced materials and for
    /// multi-density — upstream `MultiMaterial` — references).
    pub density: Option<f64>,
    /// Comment block directly above the card, top-down, markers stripped.
    pub comments: Vec<String>,
}

impl McnpMaterial {
    /// Value of a `key:` line in the comment block, trimmed.
    fn metadata_value(&self, key: &str) -> Option<&str> {
        self.comments.iter().find_map(|line| {
            let (k, v) = line.split_once(':')?;
            (k.trim().eq_ignore_ascii_case(key)).then_some(v.trim())
        })
    }

    /// `name:` metadata from the comment block (`Some("leu")` upstream).
    pub fn name(&self) -> Option<&str> {
        self.metadata_value("name")
    }

    /// `source:` metadata from the comment block.
    pub fn source(&self) -> Option<&str> {
        self.metadata_value("source")
    }

    /// Free-text `comments:` metadata: the text following the `comments:`
    /// key plus subsequent non-key comment lines, whitespace-normalized
    /// (`None` when the block carries no `comments:` key).
    pub fn comments_text(&self) -> Option<String> {
        let (idx, line) = self.comments.iter().enumerate().find(|(_, l)| {
            l.split(':')
                .next()
                .and_then(|k| k.split_whitespace().next())
                .is_some_and(|k| k.eq_ignore_ascii_case("comments"))
        })?;
        let mut parts = vec![line
            .split_once(':')
            .map(|(_, v)| v.trim())
            .unwrap_or("")
            .to_string()];
        for tail in &self.comments[idx + 1..] {
            let key = tail
                .split(':')
                .next()
                .and_then(|k| k.split_whitespace().next())
                .map(|k| k.to_ascii_lowercase());
            match key {
                Some(k) if METADATA_KEYS.contains(&k.as_str()) => break,
                _ => parts.push((*tail).to_string()),
            }
        }
        Some(parts.join(" ").trim().to_string())
    }
}

/// Parse every material card from MCNP input text.
///
/// Geometry, tallies and all other data cards are ignored; only the cell
/// cards (for densities) and `m<N>` cards are inspected.
pub fn materials_from_inp(text: &str) -> Result<Vec<McnpMaterial>, Error> {
    let lines: Vec<&str> = text.lines().collect();

    let mut cell_densities: BTreeMap<u32, Vec<f64>> = BTreeMap::new();
    let mut cards: Vec<(usize, u32)> = Vec::new();
    for (idx, line) in lines.iter().enumerate() {
        if is_cell_line(line) {
            record_cell_density(line, idx + 1, &mut cell_densities)?;
        }
        if let Some(number) = material_card_number(line, idx + 1)? {
            cards.push((idx, number));
        }
    }

    let mut materials = Vec::with_capacity(cards.len());
    for &(idx, number) in &cards {
        let data_lines = card_data_lines(&lines, idx);
        let fractions = parse_fractions(&data_lines, idx + 1)?;
        let fraction_type = detect_fraction_type(&fractions, idx + 1)?;
        let density = single_density(cell_densities.get(&number), &fractions, fraction_type);
        let comments = comment_block(&lines, idx);
        materials.push(McnpMaterial {
            number,
            fractions,
            fraction_type,
            density,
            comments,
        });
    }
    Ok(materials)
}

/// Read an MCNP input file and parse every material card.
pub fn materials_from_file(path: impl AsRef<Path>) -> Result<Vec<McnpMaterial>, Error> {
    let text = std::fs::read_to_string(path.as_ref()).map_err(|e| Error::Io(e.to_string()))?;
    materials_from_inp(&text)
}

/// Cell card carrying a density assignment.
fn is_cell_line(line: &str) -> bool {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    tokens.len() > 3
        && is_all_digits(tokens[0])
        && is_all_digits(tokens[1])
        && !tokens[2]
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic())
        && !line.starts_with("     ")
        && tokens[1] != "0"
}

fn is_all_digits(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| c.is_ascii_digit())
}

/// Accumulate a cell line's `(material, density)` assignment, dropping
/// duplicates within the upstream 1e-4 relative tolerance.
fn record_cell_density(
    line: &str,
    lineno: usize,
    cell_densities: &mut BTreeMap<u32, Vec<f64>>,
) -> Result<(), Error> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let number: u32 = tokens[1].parse().map_err(|_| Error::BadCard {
        line: lineno,
        message: format!("invalid cell material number `{}`", tokens[1]),
    })?;
    let density: f64 = tokens[2].parse().map_err(|_| Error::BadCard {
        line: lineno,
        message: format!("invalid cell density `{}`", tokens[2]),
    })?;
    let known = cell_densities.entry(number).or_default();
    let duplicate = known
        .iter()
        .any(|d| *d == density || ((density - d) / density).abs() < 1e-4);
    if !duplicate {
        known.push(density);
    }
    Ok(())
}

/// Material number if `line` opens a material card (`[mM]<digits>` token),
/// `Ok(None)` otherwise.
fn material_card_number(line: &str, lineno: usize) -> Result<Option<u32>, Error> {
    let token = match line.split_whitespace().next() {
        Some(t) => t,
        None => return Ok(None),
    };
    let mut chars = token.chars();
    match chars.next() {
        Some('m') | Some('M') => {}
        _ => return Ok(None),
    }
    match chars.next() {
        Some(c) if c.is_ascii_digit() => {}
        _ => return Ok(None),
    }
    token[1..]
        .parse::<u32>()
        .map(Some)
        .map_err(|_| Error::BadCard {
            line: lineno,
            message: format!("invalid material number `{token}`"),
        })
}

/// Raw lines forming the card body: the card line plus continuations up to
/// a blank line, EOF, or a non-continuation line. Comment continuation
/// lines (first token starting `c`/`C`) are dropped from the data, as in
/// The reference material-line grammar.
fn card_data_lines<'a>(lines: &[&'a str], card_idx: usize) -> Vec<&'a str> {
    let mut data = vec![lines[card_idx]];
    let mut i = card_idx + 1;
    while i < lines.len() {
        let line = lines[i];
        if line.split_whitespace().next().is_none() {
            break;
        }
        if !(line.starts_with("     ") || line.starts_with(['c', 'C'])) {
            break;
        }
        let first = line.split_whitespace().next().unwrap_or("").as_bytes()[0];
        if first != b'c' && first != b'C' {
            data.push(line);
        }
        i += 1;
    }
    data
}

/// Parse `zaid[.suffix] fraction` pairs (and library keywords) from the
/// card body; repeated nuclides accumulate.
fn parse_fractions(data_lines: &[&str], lineno: usize) -> Result<Vec<(NuclideId, f64)>, Error> {
    let tokens: Vec<&str> = data_lines
        .iter()
        .map(|l| l.split('$').next().unwrap_or(""))
        .flat_map(|l| l.split_whitespace())
        .skip(1)
        .collect();

    let mut pairs: Vec<(NuclideId, f64)> = Vec::new();
    let mut iter = tokens.into_iter();
    while let Some(token) = iter.next() {
        if token.contains('=') {
            let (key, _) = split_keyword(token).ok_or_else(|| Error::BadCard {
                line: lineno,
                message: format!("malformed keyword token `{token}`"),
            })?;
            if LIB_KEYWORDS.contains(&key.to_ascii_uppercase().as_str()) {
                continue;
            }
            return Err(Error::UnknownKeyword {
                line: lineno,
                keyword: key.to_string(),
            });
        }
        let fraction_token = iter.next().ok_or_else(|| Error::BadCard {
            line: lineno,
            message: format!("missing fraction after `{token}`"),
        })?;
        let zaid_text = token.split('.').next().unwrap_or("");
        let zaid = zaid_text.parse::<u32>().map_err(|_| Error::BadNumber {
            context: "zaid",
            text: token.to_string(),
        })?;
        let nuclide = nuclide_from_zaid(zaid, lineno)?;
        let fraction = fraction_token
            .parse::<f64>()
            .map_err(|_| Error::BadNumber {
                context: "fraction",
                text: fraction_token.to_string(),
            })?;
        match pairs.iter_mut().find(|(n, _)| *n == nuclide) {
            Some(slot) => slot.1 += fraction,
            None => pairs.push((nuclide, fraction)),
        }
    }
    Ok(pairs)
}

/// Split `key=value` requiring exactly one `=` and non-empty halves.
fn split_keyword(token: &str) -> Option<(&str, &str)> {
    let (key, value) = token.split_once('=')?;
    (!key.is_empty() && !value.is_empty()).then_some((key, value))
}

/// ZAID → [`NuclideId`]: natural elements (`AAA == 0`) become placeholder
/// ids (module docs), everything else goes through the shared dialect
/// converter (Am-242 swap, metastable heuristics included).
fn nuclide_from_zaid(zaid: u32, lineno: usize) -> Result<NuclideId, Error> {
    if zaid % 1_000 == 0 && zaid > 0 {
        return Ok(NuclideId::from_nucid(zaid * 10_000));
    }
    dialects::from_zaid(zaid).map_err(|_| Error::BadZaid { line: lineno, zaid })
}

/// Fraction type from the sign of the first non-zero fraction.
fn detect_fraction_type(pairs: &[(NuclideId, f64)], lineno: usize) -> Result<FracKind, Error> {
    match pairs.iter().map(|&(_, f)| f).find(|f| *f != 0.0) {
        Some(f) if f < 0.0 => Ok(FracKind::Mass),
        Some(_) => Ok(FracKind::Atom),
        None => Err(Error::BadCard {
            line: lineno,
            message: "all fractions are zero".to_string(),
        }),
    }
}

/// Converted g/cm³ density when the material has exactly one distinct cell
/// density; `None` for unreferenced and multi-density (MultiMaterial) cases.
fn single_density(
    densities: Option<&Vec<f64>>,
    pairs: &[(NuclideId, f64)],
    kind: FracKind,
) -> Option<f64> {
    let list = densities?;
    if list.len() != 1 {
        return None;
    }
    convert_density(list[0], pairs, kind)
}

/// Convert a cell density to g/cm³: negative densities are mass densities
/// (absolute value); positive ones are total atom densities in atoms/b-cm,
/// converted through the effective molar mass. Yields `None` when a needed
/// atomic mass is missing from the data tables.
fn convert_density(density: f64, pairs: &[(NuclideId, f64)], kind: FracKind) -> Option<f64> {
    if density <= 0.0 {
        return Some(-density);
    }
    let total: f64 = pairs.iter().map(|&(_, f)| f.abs()).sum();
    if total <= 0.0 {
        return None;
    }
    let molar_mass = match kind {
        FracKind::Atom => pairs
            .iter()
            .map(|&(n, f)| component_molar_mass(n).map(|m| m * (f.abs() / total)))
            .sum::<Option<f64>>()?,
        FracKind::Mass => {
            let inverse: Option<f64> = pairs
                .iter()
                .map(|&(n, f)| component_molar_mass(n).map(|m| f.abs() / total / m))
                .sum();
            1.0 / inverse?
        }
    };
    Some(density * 1e24 * molar_mass / AVOGADRO)
}

/// Atomic mass of a composition component; natural-element placeholders
/// resolve to the abundance-weighted elemental average.
fn component_molar_mass(nuclide: NuclideId) -> Option<f64> {
    if nuclide.a() == 0 {
        natural_element_mass(nuclide.z())
    } else {
        atomic_mass(nuclide.nucid())
    }
}

/// Abundance-weighted mean atomic mass of element `z` from the data tables.
fn natural_element_mass(z: u32) -> Option<f64> {
    let masses = mass_table();
    let abundances = abundance_table();
    let mut weighted = 0.0;
    let mut total = 0.0;
    for (&nucid, &mass) in masses {
        if nucid / 10_000_000 == z {
            if let Some(abundance) = abundances.get(&nucid) {
                if *abundance > 0.0 {
                    weighted += abundance * mass;
                    total += abundance;
                }
            }
        }
    }
    (total > 0.0).then(|| weighted / total)
}

/// Comment block directly above the card: `c`/`C`-led lines walked upwards,
/// stopping at a bare `c`/`C` separator, a blank line, or any other card.
/// Returned top-down with the marker token stripped.
fn comment_block(lines: &[&str], card_idx: usize) -> Vec<String> {
    let mut walked = Vec::new();
    let mut i = card_idx;
    while i > 0 {
        i -= 1;
        let line = lines[i];
        let first_token = match line.split_whitespace().next() {
            Some(t) => t,
            None => break,
        };
        if first_token != "c" && first_token != "C" {
            break;
        }
        let stripped = line.trim();
        if stripped == "c" || stripped == "C" {
            break;
        }
        walked.push(stripped[1..].trim().to_string());
    }
    walked.reverse();
    walked
}

#[cfg(test)]
mod tests {
    use super::*;

    fn inp_fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/inp/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn nucid(zaid: u32) -> u32 {
        nuclide_from_zaid(zaid, 1).unwrap().nucid()
    }

    #[test]
    fn mcnp_inp_material_count_numbers_and_leu() {
        let mats = materials_from_file(inp_fixture("mcnp_inp.txt")).unwrap();
        assert_eq!(mats.len(), 3);
        assert_eq!(
            mats.iter().map(|m| m.number).collect::<Vec<_>>(),
            vec![1, 2, 3]
        );

        let leu = &mats[0];
        assert_eq!(
            leu.fractions,
            vec![
                (NuclideId::from_nucid(922_350_000), -0.04),
                (NuclideId::from_nucid(922_380_000), -0.96)
            ]
        );
        assert_eq!(leu.fraction_type, FracKind::Mass);
        assert_eq!(leu.density, Some(19.1));
        assert_eq!(leu.name(), Some("leu"));
        assert_eq!(leu.source(), Some("Some http://URL.com"));
        assert_eq!(
            leu.comments,
            vec![
                "name: leu",
                "source: Some http://URL.com",
                "comments: first line of comments",
                "second comments",
                "third line of comments",
                "forth line of comments",
            ]
        );
    }

    #[test]
    fn mcnp_inp_water_and_default_lib_oracle() {
        let mats =
            materials_from_inp(&std::fs::read_to_string(inp_fixture("mcnp_inp.txt")).unwrap())
                .unwrap();

        // m2: natural H/O, atom fractions 2:1, three distinct cell densities
        // (-0.9, +0.1005 atom, -1.1) -> upstream MultiMaterial -> None here.
        let water = &mats[1];
        assert_eq!(water.number, 2);
        assert_eq!(
            water.fractions,
            vec![
                (NuclideId::from_nucid(10_000_000), 2.0),
                (NuclideId::from_nucid(80_000_000), 1.0)
            ]
        );
        assert_eq!(water.fraction_type, FracKind::Atom);
        assert_eq!(water.density, None);
        assert_eq!(water.name(), Some("water"));
        assert_eq!(water.source(), Some("internet"));
        assert_eq!(
            water.comments_text().as_deref(),
            Some("Here are comments: the comments continue here are more even more")
        );

        // m3: default-library keywords (hlib/nlib/plib) parse and drop.
        let with_lib = &mats[2];
        assert_eq!(with_lib.number, 3);
        assert_eq!(
            with_lib.fractions,
            vec![
                (NuclideId::from_nucid(10_000_000), 2.0),
                (NuclideId::from_nucid(80_000_000), 1.0),
                (NuclideId::from_nucid(60_000_000), 3.0),
            ]
        );
        assert_eq!(with_lib.fraction_type, FracKind::Atom);
        assert_eq!(with_lib.density, Some(1.1));
        assert_eq!(with_lib.comments, Vec::<String>::new());
        assert_eq!(with_lib.name(), None);
        assert_eq!(with_lib.comments_text(), None);
    }

    #[test]
    fn mcnp_inp_comments_repeated_nuclides_accumulate() {
        let mats = materials_from_file(inp_fixture("mcnp_inp_comments.txt")).unwrap();
        assert_eq!(mats.len(), 1);
        let leu = &mats[0];

        // U-238 accumulates -0.94 - 0.01 - 0.01 = -0.96 across comment-split
        // continuations; the commented-out 92233 line contributes nothing.
        assert_eq!(
            leu.fractions,
            vec![
                (NuclideId::from_nucid(922_350_000), -0.04),
                (NuclideId::from_nucid(922_380_000), -0.96)
            ]
        );
        assert_eq!(leu.number, 1);
        assert_eq!(leu.fraction_type, FracKind::Mass);
        assert_eq!(leu.density, Some(19.1));
        assert_eq!(
            leu.comments_text().as_deref(),
            Some("first line of comments second comments third line of comments forth line of comments")
        );
    }

    #[test]
    fn natural_element_zaids_use_placeholder_ids() {
        assert_eq!(nucid(1000), 10_000_000);
        assert_eq!(nucid(6000), 60_000_000);
        assert_eq!(nucid(8000), 80_000_000);
        let h_nat = nuclide_from_zaid(1000, 7).unwrap();
        assert_eq!(h_nat.z(), 1);
        assert_eq!(h_nat.a(), 0);
        assert_eq!(h_nat.state(), 0);
        assert_eq!(dialects::to_zaid(h_nat), 1000);
        assert!(matches!(
            nuclide_from_zaid(50_003, 1),
            Err(Error::BadZaid {
                line: 1,
                zaid: 50_003
            })
        ));
        assert!(matches!(
            nuclide_from_zaid(0, 1),
            Err(Error::BadZaid { .. })
        ));
    }

    #[test]
    fn fraction_type_detection_by_first_nonzero_sign() {
        let atom = materials_from_inp("\nsynth deck\nm1 1001 0.5 2004 0.5\n").unwrap();
        assert_eq!(atom[0].fraction_type, FracKind::Atom);
        assert_eq!(atom[0].density, None);

        let mass = materials_from_inp("\nsynth deck\nm1 1001 -0.5 2004 -0.5\n").unwrap();
        assert_eq!(mass[0].fraction_type, FracKind::Mass);

        // mixed signs follow upstream: the first non-zero decides
        let mixed = materials_from_inp("\nsynth deck\nm1 1001 -1 2004 2\n").unwrap();
        assert_eq!(mixed[0].fraction_type, FracKind::Mass);
    }

    #[test]
    fn repeated_zaid_fractions_accumulate() {
        let deck = "\nsynth deck\nm1\n     1001 0.1\n     1001 0.05\n";
        let mats = materials_from_inp(deck).unwrap();
        assert_eq!(
            mats[0].fractions,
            vec![(NuclideId::from_nucid(10_010_000), 0.150_000_000_000_000_02)]
        );
    }

    #[test]
    fn atom_density_converts_through_effective_molar_mass() {
        let deck = "\nsynth deck\n1 1 0.25 100\nm1 1001 2 8016 1\n";
        let mats = materials_from_inp(deck).unwrap();
        let expected =
            0.25 * ((2.0 * 1.007_825_031_898 + 15.994_914_619_26) / 3.0) * 1e24 / AVOGADRO;
        let got = mats[0].density.unwrap();
        assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
    }

    #[test]
    fn mass_fraction_material_atom_density_converts_harmonically() {
        let deck = "\nsynth deck\n1 1 0.25 100\nm1 1001 -2 8016 -1\n";
        let mats = materials_from_inp(deck).unwrap();
        let mh = 1.007_825_031_898;
        let mo = 15.994_914_619_26;
        let molar_inv = (2.0 / 3.0) / mh + (1.0 / 3.0) / mo;
        let expected = 0.25 * (1.0 / molar_inv) * 1e24 / AVOGADRO;
        let got = mats[0].density.unwrap();
        assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
    }

    #[test]
    fn cell_density_dedup_and_multi_density_none() {
        // two cells within the 1e-4 relative tolerance collapse to one
        let deck = "\nsynth deck\n1 1 -19.1 1\n2 1 -19.100001 1\nm1 92235 1\n";
        let mats = materials_from_inp(deck).unwrap();
        assert_eq!(mats[0].density, Some(19.1));

        // a genuinely different density makes it a MultiMaterial upstream
        let deck = "\nsynth deck\n1 1 -19.1 1\n2 1 -18.0 1\nm1 92235 1\n";
        let mats = materials_from_inp(deck).unwrap();
        assert_eq!(mats[0].density, None);
    }

    #[test]
    fn geometry_only_ptrac_deck_yields_no_materials() {
        let mats = materials_from_file(inp_fixture("mcnp_ptrac_inp.txt")).unwrap();
        assert!(mats.is_empty());
    }

    #[test]
    fn commented_out_cells_and_cards_are_inert() {
        let deck = "\nsynth deck\nc 1 1 -99.0 1\nC 2 1 -2.7 $ commented cell\nC m2 1001 1\nm1 92235 -1 $ trailing comment\n";
        let mats = materials_from_inp(deck).unwrap();
        assert_eq!(mats.len(), 1);
        assert_eq!(mats[0].number, 1);
        assert_eq!(
            mats[0].fractions,
            vec![(NuclideId::from_nucid(922_350_000), -1.0)]
        );
        assert_eq!(mats[0].density, None);
    }

    #[test]
    fn error_garbage_after_material_number() {
        let err = materials_from_inp("\nm1foo 1001 1.0\n").unwrap_err();
        assert_eq!(
            err,
            Error::BadCard {
                line: 2,
                message: "invalid material number `m1foo`".to_string()
            }
        );
    }

    #[test]
    fn error_truncated_pair() {
        let err = materials_from_inp("\nm1 1001\n").unwrap_err();
        assert!(matches!(err, Error::BadCard { line: 2, .. }));
        assert!(err.to_string().contains("missing fraction"));
    }

    #[test]
    fn error_bad_fraction_and_bad_zaid_tokens() {
        let err = materials_from_inp("\nm1 1001 abc\n").unwrap_err();
        assert_eq!(
            err,
            Error::BadNumber {
                context: "fraction",
                text: "abc".to_string()
            }
        );

        let err = materials_from_inp("\nm1 abc 1.0\n").unwrap_err();
        assert_eq!(
            err,
            Error::BadNumber {
                context: "zaid",
                text: "abc".to_string()
            }
        );
    }

    #[test]
    fn error_unknown_keyword() {
        let err = materials_from_inp("\nm1 1001 1.0 foo=bar\n").unwrap_err();
        assert_eq!(
            err,
            Error::UnknownKeyword {
                line: 2,
                keyword: "foo".to_string()
            }
        );
    }

    #[test]
    fn error_all_zero_fractions() {
        let err = materials_from_inp("\nm1 1001 0.0 1002 0\n").unwrap_err();
        assert!(matches!(err, Error::BadCard { line: 2, .. }));
        assert!(err.to_string().contains("all fractions are zero"));
    }

    #[test]
    fn io_error_on_missing_file() {
        assert!(matches!(
            materials_from_file("/definitely/not/here.inp"),
            Err(Error::Io(_))
        ));
    }

    #[test]
    fn error_display_smoke() {
        let cases = [
            Error::Io("boom".to_string()),
            Error::BadCard {
                line: 3,
                message: "why".to_string(),
            },
            Error::BadNumber {
                context: "zaid",
                text: "xx".to_string(),
            },
            Error::UnknownKeyword {
                line: 3,
                keyword: "foo".to_string(),
            },
            Error::BadZaid {
                line: 3,
                zaid: 50_003,
            },
        ];
        for e in cases {
            assert!(!e.to_string().is_empty());
        }
    }
}
