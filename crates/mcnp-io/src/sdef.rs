//! Legacy `SDEF` fixed-source reader: the typed general-source card plus the
//! discrete `SI`/`SP`/`SB` distributions it references.
//!
//! The deck core ([`crate::problem`]) carries `SDEF` as an untyped
//! [`DataCard`](crate::problem::DataCard); this module reads that card (and
//! the distribution cards it points at) into a typed [`SdefProblem`] that
//! re-emits canonical card text through [`SdefProblem::emit`]. Source
//! sampling (drawing particles from the distributions) is transport-code
//! business and is not attempted here.
//!
//! # Accepted subset (everything else is loud, never silently misread)
//!
//! `SDEF` keywords: `POS`, `CELL`, `SURF`, `VEC`, `DIR`, `AXS`, `RAD`,
//! `EXT`, `PHI`, `ERG`, `NRM`, `PAR`, `WGT`, `TME`. (`VEC` rides along because the
//! decay-source emitter this reader round-trips writes `VEC=... DIR=1` lines;
//! `AXS`/`RAD`/`EXT` were added for the plasma-source ring emitter, which
//! writes `POS=... AXS=0 0 1 RAD=D1` delta-ring cards; `PHI` was added for
//! the parametric plasma toroidal-sector emitter, which writes `PHI=Dn`
//! uniform angle-bin marginals.) Each keyword takes
//! either an inline literal (`POS`/`VEC`/`AXS` take exactly three numbers,
//! the rest take exactly one) or a `Dn` distribution reference (`ERG=D1`).
//! Any other keyword — `ARA`, `X`, `Y`, `Z`, `CCC`, abbreviations such as
//! `CEL`, and the rest — is kept verbatim in [`SdefCard::ignored`] as a
//! drift note and never parsed.
//!
//! Distributions: `SIn L ...` (discrete values), `SPn D ...` (discrete
//! probabilities, finite and non-negative), `SBn D ...` (discrete bias
//! values, finite). Any other option letter (`H`, `A`, `C`, `S`, `G`, `V`,
//! `W`, ...) is an [`SdefError::UnsupportedForm`] error, as is a missing
//! option letter (an omitted option does not mean discrete).
//!
//! # Validation ([`parse_sdef_cards`] and [`parse_sdef_text`])
//!
//! - At most one `SDEF` card per card set; a second one is an error.
//! - With no `SDEF` card, [`parse_sdef_cards`] returns `Ok(None)` and any
//!   stray `SI`/`SP`/`SB` cards are left unclaimed (they ride along as
//!   generic data cards, e.g. for out-of-scope `DSn` dependents).
//! - Duplicate keywords on the `SDEF` card, wrong value counts, unparseable
//!   numbers (including non-finite `inf`/`NaN` spellings), `D0` references,
//!   and empty tables are errors.
//! - Every `Dn` reference needs a matching `SIn` card
//!   ([`SdefError::DanglingDistribution`]); `SPn`/`SBn` without `SIn` is
//!   [`SdefError::OrphanDistribution`]; `SPn`/`SBn` entry counts must match
//!   the `SIn` entry count ([`SdefError::LengthMismatch`]).
//!
//! # Canonical emission ([`SdefProblem::emit`])
//!
//! Fields emit in fixed order `POS CELL SURF VEC DIR AXS RAD EXT PHI ERG NRM
//! WGT PAR TME`, each on its own five-space continuation line, except that a
//! literal `VEC` together with a literal `DIR` shares one `VEC=... DIR=...`
//! line — the shape the decay-source emitter writes. Distribution lists wrap
//! at 80 columns with five-space continuations. Floats render in C++
//! default-float precision-6 form (matching the decay-source emitter), so
//! that emitter's output parses and re-emits byte-identically; anything not
//! preserved (letter case of `Dn` references and keywords, `$` comments,
//! blank lines) is a documented normalization, applied loudly by
//! reconstructing the text rather than by editing it in place.
//!
//! # Out of scope
//!
//! `KCODE`/`KSRC` (different cards), `SSR`/`SSW` surface sources (covered
//! by [`crate::surfsrc`]), `FMESH` (mesh tallies), dependent `DSn`
//! distributions, and source sampling. [`DeckProblem`](crate::problem::DeckProblem)
//! exposes the typed view through `DeckProblem::sdef` and enforces this
//! module's validation in `validate`; standalone card text parses through
//! [`parse_sdef_text`].

use std::collections::BTreeMap;
use std::fmt;

use crate::inp::Error as DeckError;
use crate::problem::DataCard;

/// Maximum card-text column before distribution lists wrap (matches the
/// decay-source emitter's fixed-format card width).
const CARD_WIDTH: usize = 80;

/// Continuation indent for wrapped card lines.
const CONTINUATION_INDENT: &str = "     ";

/// Errors raised while reading `SDEF`/`SI`/`SP`/`SB` cards.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SdefError {
    /// Two cards claim the same number in one number space (`SDEF`, `SI1`,
    /// `SP1`, `SB1`, or a repeated `SDEF` keyword).
    DuplicateCard {
        /// Which card or keyword collided (`"SDEF"`, `"SI1"`, `"SDEF POS"`).
        card: String,
        /// 1-based source line of the second claim.
        line: usize,
    },
    /// A known keyword or distribution entry has the wrong shape or value.
    BadValue {
        /// Which card or keyword carries the bad value.
        card: String,
        /// 1-based source line.
        line: usize,
        /// What was wrong.
        message: String,
    },
    /// An `SDEF` keyword references distribution `n` with no `SIn` card.
    DanglingDistribution {
        /// The referenced distribution number.
        number: u32,
        /// 1-based source line of the `SDEF` card.
        line: usize,
    },
    /// An `SPn`/`SBn` card has no matching `SIn` card.
    OrphanDistribution {
        /// Which card is orphaned (`"SP3"`, `"SB3"`).
        card: String,
        /// 1-based source line.
        line: usize,
    },
    /// An `SPn`/`SBn` entry count differs from its `SIn` entry count.
    LengthMismatch {
        /// The distribution number.
        number: u32,
        /// 1-based source line.
        line: usize,
        /// What was wrong.
        message: String,
    },
    /// An `SIn` card carries no entries.
    EmptyDistribution {
        /// The distribution number.
        number: u32,
        /// 1-based source line.
        line: usize,
    },
    /// An `SI`/`SP`/`SB` card uses an option letter outside the discrete
    /// subset (`SIn L`, `SPn D`, `SBn D`).
    UnsupportedForm {
        /// Which card carries the unsupported option.
        card: String,
        /// 1-based source line.
        line: usize,
        /// What was wrong (names the accepted form).
        message: String,
    },
    /// [`parse_sdef_text`] found no `SDEF` card in the text.
    MissingCard {
        /// Which card is missing (always `"SDEF"`).
        card: String,
    },
}

impl fmt::Display for SdefError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SdefError::DuplicateCard { card, line } => {
                write!(f, "duplicate {card} card on line {line}")
            }
            SdefError::BadValue {
                card,
                line,
                message,
            } => {
                write!(f, "bad value on {card} card line {line}: {message}")
            }
            SdefError::DanglingDistribution { number, line } => {
                write!(
                    f,
                    "SDEF card line {line} references missing SI{number} distribution"
                )
            }
            SdefError::OrphanDistribution { card, line } => {
                write!(f, "{card} card line {line} has no matching SI card")
            }
            SdefError::LengthMismatch {
                number,
                line,
                message,
            } => {
                write!(f, "distribution {number} on line {line}: {message}")
            }
            SdefError::EmptyDistribution { number, line } => {
                write!(f, "SI{number} card line {line} has no entries")
            }
            SdefError::UnsupportedForm {
                card,
                line,
                message,
            } => {
                write!(f, "unsupported {card} form on line {line}: {message}")
            }
            SdefError::MissingCard { card } => {
                write!(f, "no {card} card in the text")
            }
        }
    }
}

impl std::error::Error for SdefError {}

/// Line number carried by [`SdefError::MissingCard`] when converting to the
/// deck-level error (there is no source line without input cards).
const MISSING_CARD_LINE: usize = 0;

impl From<SdefError> for DeckError {
    /// Every SDEF failure becomes a deck-level [`DeckError::BadGeometry`]
    /// carrying the original line and message, so `DeckProblem::sdef` and
    /// deck validation fit the existing `Result<_, Error>` patterns.
    fn from(error: SdefError) -> Self {
        let (line, message) = match &error {
            SdefError::DuplicateCard { line, .. }
            | SdefError::BadValue { line, .. }
            | SdefError::DanglingDistribution { line, .. }
            | SdefError::OrphanDistribution { line, .. }
            | SdefError::LengthMismatch { line, .. }
            | SdefError::EmptyDistribution { line, .. }
            | SdefError::UnsupportedForm { line, .. } => (*line, error.to_string()),
            SdefError::MissingCard { .. } => (MISSING_CARD_LINE, error.to_string()),
        };
        DeckError::BadGeometry { line, message }
    }
}

/// An inline literal value or a `Dn` distribution reference.
#[derive(Debug, Clone, PartialEq)]
pub enum SdefRef<T> {
    /// Inline literal value(s) written directly on the `SDEF` card.
    Literal(T),
    /// `Dn`: values come from distribution `n` (`SIn`, plus `SPn`/`SBn`).
    Dist(u32),
}

impl SdefRef<f64> {
    /// Canonical rendering: `D<n>` for references, C++-defaultfloat
    /// precision-6 for literals.
    pub fn render(&self) -> String {
        match self {
            SdefRef::Literal(value) => fmt_g6(*value),
            SdefRef::Dist(number) => format!("D{number}"),
        }
    }

    /// True when this field is a `Dn` distribution reference.
    pub fn is_dist(&self) -> bool {
        matches!(self, SdefRef::Dist(_))
    }
}

impl SdefRef<u32> {
    /// Canonical rendering: `D<n>` for references, plain digits otherwise.
    pub fn render(&self) -> String {
        match self {
            SdefRef::Literal(value) => format!("{value}"),
            SdefRef::Dist(number) => format!("D{number}"),
        }
    }

    /// True when this field is a `Dn` distribution reference.
    pub fn is_dist(&self) -> bool {
        matches!(self, SdefRef::Dist(_))
    }
}

impl SdefRef<[f64; 3]> {
    /// Canonical rendering: `D<n>` for references, three space-separated
    /// precision-6 numbers otherwise.
    pub fn render(&self) -> String {
        match self {
            SdefRef::Literal(value) => format!(
                "{} {} {}",
                fmt_g6(value[0]),
                fmt_g6(value[1]),
                fmt_g6(value[2])
            ),
            SdefRef::Dist(number) => format!("D{number}"),
        }
    }

    /// True when this field is a `Dn` distribution reference.
    pub fn is_dist(&self) -> bool {
        matches!(self, SdefRef::Dist(_))
    }
}

impl SdefRef<String> {
    /// Canonical rendering: `D<n>` for references, the verbatim token
    /// (case preserved) otherwise.
    pub fn render(&self) -> String {
        match self {
            SdefRef::Literal(value) => value.clone(),
            SdefRef::Dist(number) => format!("D{number}"),
        }
    }

    /// True when this field is a `Dn` distribution reference.
    pub fn is_dist(&self) -> bool {
        matches!(self, SdefRef::Dist(_))
    }
}

/// Typed `SDEF` card: one accepted keyword per field, each an inline literal
/// or a `Dn` distribution reference.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SdefCard {
    /// `POS=x y z` source position.
    pub pos: Option<SdefRef<[f64; 3]>>,
    /// `CELL=n` source cell.
    pub cell: Option<SdefRef<u32>>,
    /// `SURF=n` source surface.
    pub surf: Option<SdefRef<u32>>,
    /// `VEC=ux uy uz` reference direction vector.
    pub vec: Option<SdefRef<[f64; 3]>>,
    /// `DIR=d` direction cosine (or `1` with `VEC` for a beam).
    pub dir: Option<SdefRef<f64>>,
    /// `AXS=ax ay az` axis for `RAD`/`EXT` (the plasma-source ring shape).
    pub axs: Option<SdefRef<[f64; 3]>>,
    /// `RAD=r` radial distance from `POS` in the plane perpendicular to
    /// `AXS` (literal or `Dn`).
    pub rad: Option<SdefRef<f64>>,
    /// `EXT=x` axial extent along `AXS` (literal or `Dn`).
    pub ext: Option<SdefRef<f64>>,
    /// `PHI=p` azimuthal angle about `AXS` in radians (literal or `Dn`;
    /// the plasma-source toroidal-sector shape).
    pub phi: Option<SdefRef<f64>>,
    /// `ERG=e` energy or `ERG=Dn` energy distribution.
    pub erg: Option<SdefRef<f64>>,
    /// `NRM=n` direction cosine relative to the surface normal.
    pub nrm: Option<SdefRef<f64>>,
    /// `PAR=p` particle designator (case preserved) or `PAR=Dn`.
    pub par: Option<SdefRef<String>>,
    /// `WGT=w` particle weight.
    pub wgt: Option<SdefRef<f64>>,
    /// `TME=t` time or `TME=Dn` time distribution.
    pub tme: Option<SdefRef<f64>>,
    /// Raw tokens of keywords outside the accepted subset, in card order
    /// (loud drift notes: carried verbatim, never parsed).
    pub ignored: Vec<String>,
    /// 1-based source line of the `SDEF` card.
    pub line: usize,
}

/// One discrete source distribution: `SIn L ...` values with optional
/// `SPn D ...` probabilities and `SBn D ...` bias values.
#[derive(Debug, Clone, PartialEq)]
pub struct SdefDist {
    /// Distribution number (`SIn`/`SPn`/`SBn` share it).
    pub number: u32,
    /// `SIn` discrete values (never empty after validation).
    pub si: Vec<f64>,
    /// `SPn` discrete probabilities (`None` when no `SPn` card).
    pub sp: Option<Vec<f64>>,
    /// `SBn` discrete bias values (`None` when no `SBn` card).
    pub sb: Option<Vec<f64>>,
    /// 1-based source line of the `SIn` card.
    pub line: usize,
}

impl SdefDist {
    /// Space-separated canonical `SIn` values.
    pub fn si_text(&self) -> String {
        self.si
            .iter()
            .map(|v| fmt_g6(*v))
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// Space-separated canonical `SPn` values (`""` when no `SPn` card).
    pub fn sp_text(&self) -> String {
        self.sp
            .as_ref()
            .map(|values| {
                values
                    .iter()
                    .map(|v| fmt_g6(*v))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .unwrap_or_default()
    }

    /// Space-separated canonical `SBn` values (`""` when no `SBn` card).
    pub fn sb_text(&self) -> String {
        self.sb
            .as_ref()
            .map(|values| {
                values
                    .iter()
                    .map(|v| fmt_g6(*v))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .unwrap_or_default()
    }
}

/// Typed `SDEF` source: the card plus its discrete distributions in number
/// order.
#[derive(Debug, Clone, PartialEq)]
pub struct SdefProblem {
    /// The typed `SDEF` card.
    pub card: SdefCard,
    /// Discrete distributions in number order (referenced or not).
    pub dists: Vec<SdefDist>,
}

impl SdefProblem {
    /// Canonical card text: `SDEF` plus `SI`/`SP`/`SB` cards (see the
    /// module docs for field order, wrapping, and number formatting).
    pub fn emit(&self) -> String {
        let card = &self.card;
        let mut out = String::from("SDEF");
        if let Some(pos) = &card.pos {
            out.push_str(&format!(" POS={}", pos.render()));
        }
        if let Some(cell) = &card.cell {
            out.push_str(&format!("\n{CONTINUATION_INDENT}CELL={}", cell.render()));
        }
        if let Some(surf) = &card.surf {
            out.push_str(&format!("\n{CONTINUATION_INDENT}SURF={}", surf.render()));
        }
        match (&card.vec, &card.dir) {
            // The decay-source emitter's joint line shape.
            (Some(vec), Some(dir)) if !vec.is_dist() && !dir.is_dist() => {
                out.push_str(&format!(
                    "\n{CONTINUATION_INDENT}VEC={} DIR={}",
                    vec.render(),
                    dir.render()
                ));
            }
            _ => {
                if let Some(vec) = &card.vec {
                    out.push_str(&format!("\n{CONTINUATION_INDENT}VEC={}", vec.render()));
                }
                if let Some(dir) = &card.dir {
                    out.push_str(&format!("\n{CONTINUATION_INDENT}DIR={}", dir.render()));
                }
            }
        }
        if let Some(axs) = &card.axs {
            out.push_str(&format!("\n{CONTINUATION_INDENT}AXS={}", axs.render()));
        }
        if let Some(rad) = &card.rad {
            out.push_str(&format!("\n{CONTINUATION_INDENT}RAD={}", rad.render()));
        }
        if let Some(ext) = &card.ext {
            out.push_str(&format!("\n{CONTINUATION_INDENT}EXT={}", ext.render()));
        }
        if let Some(phi) = &card.phi {
            out.push_str(&format!("\n{CONTINUATION_INDENT}PHI={}", phi.render()));
        }
        if let Some(erg) = &card.erg {
            out.push_str(&format!("\n{CONTINUATION_INDENT}ERG={}", erg.render()));
        }
        if let Some(nrm) = &card.nrm {
            out.push_str(&format!("\n{CONTINUATION_INDENT}NRM={}", nrm.render()));
        }
        if let Some(wgt) = &card.wgt {
            out.push_str(&format!("\n{CONTINUATION_INDENT}WGT={}", wgt.render()));
        }
        if let Some(par) = &card.par {
            out.push_str(&format!("\n{CONTINUATION_INDENT}PAR={}", par.render()));
        }
        if let Some(tme) = &card.tme {
            out.push_str(&format!("\n{CONTINUATION_INDENT}TME={}", tme.render()));
        }
        for dist in &self.dists {
            push_list(
                &mut out,
                &format!("SI{} L", dist.number),
                &dist.si.iter().map(|v| fmt_g6(*v)).collect::<Vec<_>>(),
            );
            if let Some(sp) = &dist.sp {
                push_list(
                    &mut out,
                    &format!("SP{} D", dist.number),
                    &sp.iter().map(|v| fmt_g6(*v)).collect::<Vec<_>>(),
                );
            }
            if let Some(sb) = &dist.sb {
                push_list(
                    &mut out,
                    &format!("SB{} D", dist.number),
                    &sb.iter().map(|v| fmt_g6(*v)).collect::<Vec<_>>(),
                );
            }
        }
        out
    }
}

/// Read the `SDEF` source (plus its `SI`/`SP`/`SB` distributions) from deck
/// data cards: `Ok(None)` when the deck has no `SDEF` card (stray
/// distributions stay unclaimed generic cards), otherwise the typed source
/// with every validation rule from the module docs enforced.
pub fn parse_sdef_cards(cards: &[DataCard]) -> Result<Option<SdefProblem>, SdefError> {
    let mut sdef: Option<&DataCard> = None;
    for card in cards {
        if card.name.eq_ignore_ascii_case("SDEF") {
            if sdef.is_some() {
                return Err(SdefError::DuplicateCard {
                    card: "SDEF".to_string(),
                    line: card.line,
                });
            }
            sdef = Some(card);
        }
    }
    let Some(card) = sdef else {
        return Ok(None);
    };
    let model = parse_sdef_card(card)?;
    let mut si: BTreeMap<u32, (Vec<f64>, usize)> = BTreeMap::new();
    let mut sp: BTreeMap<u32, (Vec<f64>, usize)> = BTreeMap::new();
    let mut sb: BTreeMap<u32, (Vec<f64>, usize)> = BTreeMap::new();
    for card in cards {
        let Some((kind, number)) = split_dist_name(&card.name) else {
            continue;
        };
        let (target, tag) = match kind {
            DistKind::Si => (&mut si, "SI"),
            DistKind::Sp => (&mut sp, "SP"),
            DistKind::Sb => (&mut sb, "SB"),
        };
        if target.contains_key(&number) {
            return Err(SdefError::DuplicateCard {
                card: format!("{tag}{number}"),
                line: card.line,
            });
        }
        target.insert(number, (parse_dist_values(card, kind)?, card.line));
    }
    let mut dists = Vec::with_capacity(si.len());
    for (number, (values, line)) in &si {
        if values.is_empty() {
            return Err(SdefError::EmptyDistribution {
                number: *number,
                line: *line,
            });
        }
        let probabilities = match sp.remove(number) {
            Some((values, line)) => {
                if values.len() != si[number].0.len() {
                    return Err(SdefError::LengthMismatch {
                        number: *number,
                        line,
                        message: format!(
                            "SP{number} has {} entries for {} SI{number} values",
                            values.len(),
                            si[number].0.len()
                        ),
                    });
                }
                Some(values)
            }
            None => None,
        };
        let biases = match sb.remove(number) {
            Some((values, line)) => {
                if values.len() != si[number].0.len() {
                    return Err(SdefError::LengthMismatch {
                        number: *number,
                        line,
                        message: format!(
                            "SB{number} has {} entries for {} SI{number} values",
                            values.len(),
                            si[number].0.len()
                        ),
                    });
                }
                Some(values)
            }
            None => None,
        };
        dists.push(SdefDist {
            number: *number,
            si: values.clone(),
            sp: probabilities,
            sb: biases,
            line: *line,
        });
    }
    for (orphans, tag) in [(&sp, "SP"), (&sb, "SB")] {
        if let Some((number, (_, line))) = orphans.iter().next() {
            return Err(SdefError::OrphanDistribution {
                card: format!("{tag}{number}"),
                line: *line,
            });
        }
    }
    let known: BTreeMap<u32, &SdefDist> = dists.iter().map(|d| (d.number, d)).collect();
    for number in referenced_dists(&model) {
        if !known.contains_key(&number) {
            return Err(SdefError::DanglingDistribution {
                number,
                line: model.line,
            });
        }
    }
    Ok(Some(SdefProblem { card: model, dists }))
}

/// Parse standalone `SDEF` card text (an `SDEF` card plus `SI`/`SP`/`SB`
/// cards, e.g. the decay-source emitter's output) into the typed source.
///
/// Lines starting with five spaces continue the current card; anything else
/// starts a new card whose first token is its name; blank lines are
/// skipped. [`SdefError::MissingCard`] when no `SDEF` card is present.
pub fn parse_sdef_text(text: &str) -> Result<SdefProblem, SdefError> {
    let mut cards = Vec::new();
    let mut current: Option<(usize, Vec<String>)> = None;
    for (index, line) in text.lines().enumerate() {
        let lineno = index + 1;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with(CONTINUATION_INDENT) && current.is_some() {
            if let Some((_, raw)) = current.as_mut() {
                raw.push(line.to_string());
            }
            continue;
        }
        if let Some((line, raw)) = current.take() {
            cards.push(make_text_card(line, raw));
        }
        current = Some((lineno, vec![line.to_string()]));
    }
    if let Some((line, raw)) = current.take() {
        cards.push(make_text_card(line, raw));
    }
    match parse_sdef_cards(&cards)? {
        Some(problem) => Ok(problem),
        None => Err(SdefError::MissingCard {
            card: "SDEF".to_string(),
        }),
    }
}

/// One text card from its raw lines: the first token (uppercased) is the
/// name, matching [`crate::problem`]'s data-card classification.
fn make_text_card(line: usize, raw: Vec<String>) -> DataCard {
    let first = raw.first().map(String::as_str).unwrap_or("");
    let code = first.split('$').next().unwrap_or("");
    let mut tokens = code.split_whitespace();
    DataCard {
        name: tokens.next().unwrap_or("").to_ascii_uppercase(),
        args: tokens.map(str::to_string).collect(),
        line,
        raw_lines: raw,
    }
}

/// Which distribution family a card name belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DistKind {
    /// `SIn` values.
    Si,
    /// `SPn` probabilities.
    Sp,
    /// `SBn` bias values.
    Sb,
}

/// Split a card name into its distribution family and number (`SI1` to
/// `(Si, 1)`); `None` for anything else (including numberless `SI`).
fn split_dist_name(name: &str) -> Option<(DistKind, u32)> {
    let upper = name.to_ascii_uppercase();
    let (kind, rest) = if let Some(rest) = upper.strip_prefix("SI") {
        (DistKind::Si, rest)
    } else if let Some(rest) = upper.strip_prefix("SP") {
        (DistKind::Sp, rest)
    } else {
        (DistKind::Sb, upper.strip_prefix("SB")?)
    };
    if rest.is_empty() || !rest.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    rest.parse::<u32>().ok().map(|number| (kind, number))
}

/// Parse one `SI`/`SP`/`SB` card body: the leading option letter must be the
/// discrete form (`L` for `SI`, `D` for `SP`/`SB`), followed by finite
/// numbers (`SPn` entries must also be non-negative probabilities).
fn parse_dist_values(card: &DataCard, kind: DistKind) -> Result<Vec<f64>, SdefError> {
    let tag = match kind {
        DistKind::Si => "SI",
        DistKind::Sp => "SP",
        DistKind::Sb => "SB",
    };
    let number_text = card
        .name
        .to_ascii_uppercase()
        .strip_prefix(tag)
        .unwrap_or("")
        .to_string();
    let label = format!("{tag}{number_text}");
    let number: u32 = number_text.parse().map_err(|_| SdefError::BadValue {
        card: label.clone(),
        line: card.line,
        message: format!("bad distribution number `{number_text}` (numbers start at 1)"),
    })?;
    if number < 1 {
        return Err(SdefError::BadValue {
            card: label.clone(),
            line: card.line,
            message: format!("bad distribution number `{number_text}` (numbers start at 1)"),
        });
    }
    let wanted = match kind {
        DistKind::Si => 'L',
        DistKind::Sp | DistKind::Sb => 'D',
    };
    let tokens = crate::semantic::card_tokens(card);
    let Some((option, values)) = tokens.split_first() else {
        return Err(SdefError::UnsupportedForm {
            card: label.clone(),
            line: card.line,
            message: format!("{label} needs the discrete `{wanted}` option with values"),
        });
    };
    if !option.eq_ignore_ascii_case(&wanted.to_string()) {
        let accepted = match kind {
            DistKind::Si => "only the discrete `L` option is supported",
            DistKind::Sp | DistKind::Sb => "only the discrete `D` option is supported",
        };
        return Err(SdefError::UnsupportedForm {
            card: label,
            line: card.line,
            message: format!("{tag} option `{option}` is not supported ({accepted})"),
        });
    }
    let mut out = Vec::with_capacity(values.len());
    for token in values {
        let value: f64 = token.parse().map_err(|_| SdefError::BadValue {
            card: label.clone(),
            line: card.line,
            message: format!("cannot parse `{token}` as a number"),
        })?;
        if !value.is_finite() {
            return Err(SdefError::BadValue {
                card: label.clone(),
                line: card.line,
                message: format!("non-finite entry `{token}`"),
            });
        }
        if kind == DistKind::Sp && value < 0.0 {
            return Err(SdefError::BadValue {
                card: label.clone(),
                line: card.line,
                message: format!("negative SP probability `{token}`"),
            });
        }
        out.push(value);
    }
    Ok(out)
}

/// Parse one `SDEF` data card into the typed card (keyword spanning,
/// duplicate detection, and value shapes).
fn parse_sdef_card(card: &DataCard) -> Result<SdefCard, SdefError> {
    let tokens = crate::semantic::card_tokens(card);
    // Group `KEY=first` heads with their following bare value tokens.
    let mut spans: Vec<(String, String, Vec<String>)> = Vec::new();
    for token in &tokens {
        match token.split_once('=') {
            Some((key, first)) => {
                if key.is_empty() {
                    return Err(bad_value(card, "SDEF", format!("bad keyword `{token}`")));
                }
                let mut values = Vec::new();
                if !first.is_empty() {
                    values.push(first.to_string());
                }
                spans.push((key.to_string(), token.clone(), values));
            }
            None => {
                let Some((_, _, values)) = spans.last_mut() else {
                    return Err(bad_value(
                        card,
                        "SDEF",
                        format!("expected KEY=value, found `{token}`"),
                    ));
                };
                values.push(token.clone());
            }
        }
    }
    let mut model = SdefCard {
        line: card.line,
        ..SdefCard::default()
    };
    for (key, head, values) in &spans {
        let upper = key.to_ascii_uppercase();
        let slot: Option<&str> = match upper.as_str() {
            "POS" | "CELL" | "SURF" | "VEC" | "DIR" | "AXS" | "RAD" | "EXT" | "PHI" | "ERG"
            | "NRM" | "PAR" | "WGT" | "TME" => Some(upper.as_str()),
            _ => None,
        };
        let Some(slot) = slot else {
            // Loud drift note: carried verbatim, never parsed.
            let mut raw = vec![head.clone()];
            raw.extend(values.clone());
            model.ignored.push(raw.join(" "));
            continue;
        };
        if field_is_set(&model, slot) {
            return Err(SdefError::DuplicateCard {
                card: format!("SDEF {slot}"),
                line: card.line,
            });
        }
        match slot {
            "POS" => model.pos = Some(parse_triplet(card, "SDEF POS", values)?),
            "VEC" => model.vec = Some(parse_triplet(card, "SDEF VEC", values)?),
            "AXS" => model.axs = Some(parse_triplet(card, "SDEF AXS", values)?),
            "CELL" => model.cell = Some(parse_uint(card, "SDEF CELL", values)?),
            "SURF" => model.surf = Some(parse_uint(card, "SDEF SURF", values)?),
            "ERG" => model.erg = Some(parse_float(card, "SDEF ERG", values)?),
            "DIR" => model.dir = Some(parse_float(card, "SDEF DIR", values)?),
            "RAD" => model.rad = Some(parse_float(card, "SDEF RAD", values)?),
            "EXT" => model.ext = Some(parse_float(card, "SDEF EXT", values)?),
            "PHI" => model.phi = Some(parse_float(card, "SDEF PHI", values)?),
            "NRM" => model.nrm = Some(parse_float(card, "SDEF NRM", values)?),
            "WGT" => model.wgt = Some(parse_float(card, "SDEF WGT", values)?),
            "TME" => model.tme = Some(parse_float(card, "SDEF TME", values)?),
            "PAR" => {
                if values.len() != 1 {
                    return Err(bad_value(
                        card,
                        "SDEF PAR",
                        format!("needs exactly one value, found {}", values.len()),
                    ));
                }
                model.par = Some(parse_particle(&values[0]));
            }
            _ => unreachable!("slot is one of the fourteen accepted keywords"),
        }
    }
    Ok(model)
}

/// True when the named `SDEF` keyword slot is already set (duplicate
/// detection).
fn field_is_set(model: &SdefCard, slot: &str) -> bool {
    match slot {
        "POS" => model.pos.is_some(),
        "CELL" => model.cell.is_some(),
        "SURF" => model.surf.is_some(),
        "VEC" => model.vec.is_some(),
        "DIR" => model.dir.is_some(),
        "AXS" => model.axs.is_some(),
        "RAD" => model.rad.is_some(),
        "EXT" => model.ext.is_some(),
        "PHI" => model.phi.is_some(),
        "ERG" => model.erg.is_some(),
        "NRM" => model.nrm.is_some(),
        "PAR" => model.par.is_some(),
        "WGT" => model.wgt.is_some(),
        "TME" => model.tme.is_some(),
        _ => false,
    }
}

/// Shorthand for keyword value errors.
fn bad_value(card: &DataCard, what: &str, message: String) -> SdefError {
    SdefError::BadValue {
        card: what.to_string(),
        line: card.line,
        message,
    }
}

/// `Dn` (case-insensitive) with `n >= 1`; `None` for anything else.
fn parse_dist_token(token: &str) -> Option<u32> {
    let rest = token.strip_prefix(['D', 'd'])?;
    if rest.is_empty() || !rest.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    rest.parse::<u32>().ok().filter(|n| *n >= 1)
}

/// One finite float or a `Dn` reference (exactly one token).
fn parse_float(card: &DataCard, what: &str, values: &[String]) -> Result<SdefRef<f64>, SdefError> {
    if values.len() != 1 {
        return Err(bad_value(
            card,
            what,
            format!("needs exactly one value, found {}", values.len()),
        ));
    }
    let token = &values[0];
    if let Some(number) = parse_dist_token(token) {
        return Ok(SdefRef::Dist(number));
    }
    let value: f64 = token
        .parse()
        .map_err(|_| bad_value(card, what, format!("cannot parse `{token}` as a number")))?;
    if !value.is_finite() {
        return Err(bad_value(card, what, format!("non-finite value `{token}`")));
    }
    Ok(SdefRef::Literal(value))
}

/// One non-negative integer or a `Dn` reference (exactly one token).
fn parse_uint(card: &DataCard, what: &str, values: &[String]) -> Result<SdefRef<u32>, SdefError> {
    if values.len() != 1 {
        return Err(bad_value(
            card,
            what,
            format!("needs exactly one value, found {}", values.len()),
        ));
    }
    let token = &values[0];
    if let Some(number) = parse_dist_token(token) {
        return Ok(SdefRef::Dist(number));
    }
    token.parse::<u32>().map(SdefRef::Literal).map_err(|_| {
        bad_value(
            card,
            what,
            format!("cannot parse `{token}` as an integer ≥ 0"),
        )
    })
}

/// Three finite floats or a single `Dn` reference.
fn parse_triplet(
    card: &DataCard,
    what: &str,
    values: &[String],
) -> Result<SdefRef<[f64; 3]>, SdefError> {
    if values.len() == 1 {
        if let Some(number) = parse_dist_token(&values[0]) {
            return Ok(SdefRef::Dist(number));
        }
    }
    if values.len() != 3 {
        return Err(bad_value(
            card,
            what,
            format!("needs exactly three numbers, found {}", values.len()),
        ));
    }
    let mut out = [0.0; 3];
    for (index, token) in values.iter().enumerate() {
        let value: f64 = token
            .parse()
            .map_err(|_| bad_value(card, what, format!("cannot parse `{token}` as a number")))?;
        if !value.is_finite() {
            return Err(bad_value(card, what, format!("non-finite value `{token}`")));
        }
        out[index] = value;
    }
    Ok(SdefRef::Literal(out))
}

/// One particle token (verbatim, case preserved) or a `Dn` reference.
fn parse_particle(token: &str) -> SdefRef<String> {
    match parse_dist_token(token) {
        Some(number) => SdefRef::Dist(number),
        None => SdefRef::Literal(token.to_string()),
    }
}

/// Every distribution number referenced from the `SDEF` card's `Dn` fields.
fn referenced_dists(card: &SdefCard) -> Vec<u32> {
    let mut numbers = Vec::new();
    let mut push = |field: &Option<SdefRef<f64>>| {
        if let Some(SdefRef::Dist(number)) = field {
            numbers.push(*number);
        }
    };
    push(&card.erg);
    push(&card.dir);
    push(&card.nrm);
    push(&card.wgt);
    push(&card.tme);
    push(&card.rad);
    push(&card.ext);
    push(&card.phi);
    for field in [&card.pos, &card.vec, &card.axs] {
        if let Some(SdefRef::Dist(number)) = field {
            numbers.push(*number);
        }
    }
    for field in [&card.cell, &card.surf] {
        if let Some(SdefRef::Dist(number)) = field {
            numbers.push(*number);
        }
    }
    if let Some(SdefRef::Dist(number)) = &card.par {
        numbers.push(*number);
    }
    numbers.sort_unstable();
    numbers.dedup();
    numbers
}

/// Append a distribution card (`head` plus entries), wrapping onto
/// five-space continuation lines past [`CARD_WIDTH`] — the decay-source
/// emitter's list shape, so its output re-emits identically.
fn push_list(out: &mut String, head: &str, entries: &[String]) {
    out.push('\n');
    out.push_str(head);
    let mut col = head.len();
    for entry in entries {
        let need = 1 + entry.len();
        if col + need > CARD_WIDTH {
            out.push('\n');
            out.push_str(CONTINUATION_INDENT);
            col = CONTINUATION_INDENT.len();
        }
        out.push(' ');
        out.push_str(entry);
        col += need;
    }
}

/// Format like C++ `std::ostream <<` default formatting (`%g` with precision
/// 6, trailing zeros stripped): the decay-source emitter's number shape, so
/// its output parses and re-emits byte-identically. This copy must stay in
/// sync with that emitter's formatting.
fn fmt_g6(x: f64) -> String {
    debug_assert!(x.is_finite());
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    const PRECISION: i32 = 6;
    let sci = format!("{:.*e}", (PRECISION - 1) as usize, x);
    let (mant, exp) = sci.split_once('e').expect("scientific format has e");
    let exp: i32 = exp.parse().expect("scientific exponent is an integer");
    let mant = trim_zeros(mant);
    if !(-4..PRECISION).contains(&exp) {
        format!("{mant}e{exp:+03}")
    } else {
        let decimals = (PRECISION - 1 - exp).max(0) as usize;
        trim_zeros(&format!("{x:.decimals$}"))
    }
}

/// Strip trailing fractional zeros (`"1.40000"` to `"1.4"`, `"14.000"` to
/// `"14"`); integer strings pass through.
fn trim_zeros(s: &str) -> String {
    if s.contains('.') {
        let trimmed = s.trim_end_matches('0');
        trimmed.strip_suffix('.').unwrap_or(trimmed).to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::problem::parse_deck;

    /// Single-line isotropic decay-source card (the emitter's golden shape).
    const SINGLE_LINE: &str = "SDEF POS=0 0 0\n     ERG=0.662\n     WGT=1\n     PAR=n";

    /// Multi-line distribution card (the emitter's golden shape).
    const MULTI_LINE: &str = "SDEF POS=0 0 0\n     ERG=D1\n     WGT=1\n     PAR=n\
         \nSI1 L 0.662 1.17 1.33\nSP1 D 0.5 0.25 0.25";

    /// Beam card with the joint `VEC`/`DIR` line (the emitter's golden shape).
    const BEAM_LINE: &str = "SDEF POS=1 2 3\n     VEC=0 0 1 DIR=1\
         \n     ERG=0.662\n     WGT=0.5\n     PAR=p";

    #[test]
    fn single_line_card_round_trips_byte_identical() {
        let problem = parse_sdef_text(SINGLE_LINE).unwrap();
        assert_eq!(problem.card.pos, Some(SdefRef::Literal([0.0, 0.0, 0.0])));
        assert_eq!(problem.card.erg, Some(SdefRef::Literal(0.662)));
        assert_eq!(problem.card.wgt, Some(SdefRef::Literal(1.0)));
        assert_eq!(problem.card.par, Some(SdefRef::Literal("n".to_string())));
        assert!(problem.dists.is_empty());
        assert!(problem.card.ignored.is_empty());
        assert_eq!(problem.emit(), SINGLE_LINE);
    }

    #[test]
    fn multi_line_distribution_card_round_trips_byte_identical() {
        let problem = parse_sdef_text(MULTI_LINE).unwrap();
        assert_eq!(problem.card.erg, Some(SdefRef::Dist(1)));
        assert_eq!(problem.dists.len(), 1);
        assert_eq!(problem.dists[0].number, 1);
        assert_eq!(problem.dists[0].si, vec![0.662, 1.17, 1.33]);
        assert_eq!(
            problem.dists[0].sp,
            Some(vec![0.5, 0.25, 0.25]),
            "probabilities parse in card order"
        );
        assert_eq!(problem.dists[0].sb, None);
        assert_eq!(problem.emit(), MULTI_LINE);
    }

    #[test]
    fn beam_card_joint_vec_dir_line_round_trips() {
        let problem = parse_sdef_text(BEAM_LINE).unwrap();
        assert_eq!(problem.card.vec, Some(SdefRef::Literal([0.0, 0.0, 1.0])));
        assert_eq!(problem.card.dir, Some(SdefRef::Literal(1.0)));
        assert_eq!(problem.emit(), BEAM_LINE);
    }

    #[test]
    fn deck_cards_wire_through_problem_and_validate() {
        let text = "msg\ntitle\n1 0 -1\n\n1 so 1.0\n\nSDEF POS=0 0 0\n     ERG=D1\n     WGT=1\n     PAR=n\nSI1 L 0.662 1.17\nSP1 D 0.5 0.5\n";
        let deck = parse_deck(text).unwrap();
        let problem = deck.sdef().unwrap().expect("deck carries SDEF");
        assert_eq!(problem.card.erg, Some(SdefRef::Dist(1)));
        assert_eq!(problem.dists.len(), 1);
        deck.validate().unwrap();

        let bare = "msg\ntitle\n1 0 -1\n\n1 so 1.0\n\nmode n\n";
        let deck = parse_deck(bare).unwrap();
        assert!(deck.sdef().unwrap().is_none());

        // A dangling distribution reference fails deck validation too.
        let bad = "msg\ntitle\n1 0 -1\n\n1 so 1.0\n\nSDEF ERG=D2\n";
        let deck = parse_deck(bad).unwrap();
        let err = deck.validate().unwrap_err();
        assert!(err.to_string().contains("missing SI2"), "{err}");
    }

    #[test]
    fn unknown_keywords_are_loud_drift_notes() {
        let text = "SDEF POS=0 0 0 ERG=0.662 ARA=1.0 CCC=5";
        let problem = parse_sdef_text(text).unwrap();
        assert_eq!(
            problem.card.ignored,
            vec!["ARA=1.0 1.0".to_string(), "CCC=5 5".to_string()]
        );
        // Ignored keywords do not disturb the typed fields.
        assert_eq!(problem.card.erg, Some(SdefRef::Literal(0.662)));
    }

    #[test]
    fn ring_shape_with_axs_rad_ext_round_trips_byte_identical() {
        // The plasma-source ring emitter's shape: axis + radial delta ring.
        let text = "SDEF POS=0 0 25\n     AXS=0 0 1\n     RAD=D1\n     ERG=14.021\n     WGT=1\
                    \n     PAR=n\nSI1 L 300 300\nSP1 D 0 1";
        let problem = parse_sdef_text(text).unwrap();
        assert_eq!(problem.card.axs, Some(SdefRef::Literal([0.0, 0.0, 1.0])));
        assert_eq!(problem.card.rad, Some(SdefRef::Dist(1)));
        assert!(problem.card.ext.is_none());
        assert_eq!(problem.emit(), text);
        // Case-insensitive keyword input normalizes to the canonical shape.
        let lower = parse_sdef_text("sdef pos=0 0 25\n     axs=0 0 1\n     rad=d1\n     erg=14.021\n     wgt=1\n     par=n\nsi1 l 300 300\nsp1 d 0 1").unwrap();
        assert_eq!(lower.emit(), text);
    }

    #[test]
    fn ext_takes_literals_and_distribution_references() {
        let text = "SDEF POS=0 0 0\
                    \n     AXS=0 0 1\
                    \n     RAD=D1\
                    \n     EXT=D2\
                    \n     ERG=14.1\
                    \nSI1 L 300 300\
                    \nSP1 D 0 1\
                    \nSI2 L -5 5\
                    \nSP2 D 0 1";
        let problem = parse_sdef_text(text).unwrap();
        assert_eq!(problem.card.ext, Some(SdefRef::Dist(2)));
        assert_eq!(problem.dists.len(), 2);
        assert_eq!(problem.emit(), text);

        let bad = parse_sdef_text("SDEF POS=0 0 0 AXS=0 0 RAD=D1").unwrap_err();
        assert!(bad.to_string().contains("exactly three numbers"), "{bad}");
        let dup = parse_sdef_text("SDEF RAD=1 RAD=2").unwrap_err();
        assert!(dup.to_string().contains("duplicate"), "{dup}");
    }

    #[test]
    fn phi_takes_literals_and_distribution_references() {
        // The plasma-source toroidal-sector shape: PHI between EXT and ERG.
        let text = "SDEF POS=0 0 0\
                    \n     AXS=0 0 1\
                    \n     RAD=D1\
                    \n     EXT=D2\
                    \n     PHI=D3\
                    \n     ERG=14.1\
                    \nSI1 L 300 300\
                    \nSP1 D 0 1\
                    \nSI2 L -5 5\
                    \nSP2 D 0 1\
                    \nSI3 L 0.5 1\
                    \nSP3 D 0.5 0.5";
        let problem = parse_sdef_text(text).unwrap();
        assert_eq!(problem.card.phi, Some(SdefRef::Dist(3)));
        assert_eq!(problem.dists.len(), 3);
        assert_eq!(problem.emit(), text);
        // Lowercase normalizes; duplicates stay loud.
        let lower = parse_sdef_text("sdef pos=0 0 0\n     phi=1.5").unwrap();
        assert_eq!(lower.card.phi, Some(SdefRef::Literal(1.5)));
        let dup = parse_sdef_text("SDEF PHI=1 PHI=2").unwrap_err();
        assert!(dup.to_string().contains("duplicate"), "{dup}");
    }

    #[test]
    fn duplicate_sdef_cards_rejected() {
        let deck = parse_deck("msg\ntitle\n1 0 -1\n\n1 so 1\n\nSDEF ERG=1\nSDEF ERG=2\n").unwrap();
        let err = parse_sdef_cards(&deck.data).unwrap_err();
        assert_eq!(
            err,
            SdefError::DuplicateCard {
                card: "SDEF".to_string(),
                line: 8,
            }
        );
    }

    #[test]
    fn duplicate_keywords_rejected() {
        let err = parse_sdef_text("SDEF ERG=1 ERG=2").unwrap_err();
        assert_eq!(
            err,
            SdefError::DuplicateCard {
                card: "SDEF ERG".to_string(),
                line: 1,
            }
        );
    }

    #[test]
    fn non_discrete_distribution_forms_rejected() {
        for (text, fragment) in [
            ("SDEF ERG=D1\nSI1 H 0.662 1.17", "not supported"),
            ("SDEF ERG=D1\nSI1 L 0.662\nSP1 V 0.5", "not supported"),
            ("SDEF ERG=D1\nSI1 0.662", "only the discrete"),
        ] {
            let err = parse_sdef_text(text).unwrap_err();
            assert!(
                matches!(err, SdefError::UnsupportedForm { .. }),
                "{text}: {err}"
            );
            assert!(err.to_string().contains(fragment), "{text}: {err}");
        }
    }

    #[test]
    fn dangling_orphan_and_length_errors() {
        let err = parse_sdef_text("SDEF ERG=D2\nSI1 L 1.0").unwrap_err();
        assert_eq!(err, SdefError::DanglingDistribution { number: 2, line: 1 });

        let err = parse_sdef_text("SDEF POS=0 0 0\nSP1 D 0.5").unwrap_err();
        assert_eq!(
            err,
            SdefError::OrphanDistribution {
                card: "SP1".to_string(),
                line: 2,
            }
        );

        let err = parse_sdef_text("SDEF ERG=D1\nSI1 L 1.0 2.0\nSP1 D 1.0").unwrap_err();
        assert!(
            matches!(err, SdefError::LengthMismatch { number: 1, .. }),
            "{err}"
        );

        let err = parse_sdef_text("SDEF ERG=D0").unwrap_err();
        assert!(matches!(err, SdefError::BadValue { .. }), "{err}");
        assert!(err.to_string().contains("cannot parse"), "{err}");
    }

    #[test]
    fn bad_values_are_loud() {
        for (text, fragment) in [
            ("SDEF POS=0 0", "exactly three numbers"),
            ("SDEF ERG=X", "cannot parse `X`"),
            ("SDEF ERG=inf", "non-finite"),
            ("SDEF CELL=-1", "integer"),
            ("SDEF WGT", "KEY=value"),
            ("SDEF ERG=D1\nSI1 L", "no entries"),
            ("SDEF ERG=D1\nSI1 L 1.0\nSP1 D -0.5", "negative SP"),
            ("SDEF ERG=D1\nSI1 L X", "cannot parse `X`"),
            ("X", "no SDEF card"),
        ] {
            let err = match parse_sdef_text(text) {
                Ok(problem) => panic!("{text} parsed unexpectedly: {problem:?}"),
                Err(err) => err,
            };
            assert!(err.to_string().contains(fragment), "{text}: {err}");
        }
    }

    #[test]
    fn keyword_matching_is_case_insensitive_with_case_normalized_emit() {
        let problem = parse_sdef_text("sdef pos=0 0 0\n     erg=d1\nsi1 l 0.662\nsp1 d 1").unwrap();
        assert_eq!(problem.card.erg, Some(SdefRef::Dist(1)));
        assert_eq!(
            problem.emit(),
            "SDEF POS=0 0 0\n     ERG=D1\nSI1 L 0.662\nSP1 D 1"
        );
    }

    #[test]
    fn sb_biases_ride_along_and_must_match_si_length() {
        let problem =
            parse_sdef_text("SDEF ERG=D1\nSI1 L 1.0 2.0\nSP1 D 0.5 0.5\nSB1 D 1 2").unwrap();
        assert_eq!(problem.dists[0].sb, Some(vec![1.0, 2.0]));
        assert!(problem.emit().contains("\nSB1 D 1 2"));

        let err = parse_sdef_text("SDEF ERG=D1\nSI1 L 1.0 2.0\nSB1 D 1").unwrap_err();
        assert!(matches!(err, SdefError::LengthMismatch { .. }), "{err}");
    }

    #[test]
    fn long_lists_wrap_at_card_width_like_the_emitter() {
        let energies: Vec<String> = (1..=12).map(|i| format!("0.{}", i)).collect();
        let text = format!("SDEF ERG=D1\nSI1 L {}", energies.join(" "));
        let problem = parse_sdef_text(&text).unwrap();
        let emitted = problem.emit();
        for line in emitted.lines() {
            assert!(line.len() <= CARD_WIDTH, "line over width: {line:?}");
        }
        // Re-parsing the emission is a fixed point on card text (source line
        // numbers legitimately shift when wrapping moves entries).
        assert_eq!(parse_sdef_text(&emitted).unwrap().emit(), emitted);
    }

    #[test]
    fn number_formatting_matches_the_emitter_shape() {
        assert_eq!(fmt_g6(0.662), "0.662");
        assert_eq!(fmt_g6(0.5), "0.5");
        assert_eq!(fmt_g6(12.0), "12");
        assert_eq!(fmt_g6(-1.25), "-1.25");
        assert_eq!(fmt_g6(0.0), "0");
        assert_eq!(fmt_g6(1e-5), "1e-05");
    }

    #[test]
    fn sdef_error_converts_to_deck_error_with_line() {
        let err: DeckError = SdefError::DanglingDistribution { number: 2, line: 7 }.into();
        assert!(err.to_string().contains("missing SI2"));
        assert!(matches!(err, DeckError::BadGeometry { line: 7, .. }));
        let err: DeckError = SdefError::MissingCard {
            card: "SDEF".to_string(),
        }
        .into();
        assert!(matches!(err, DeckError::BadGeometry { line: 0, .. }));
    }
}
