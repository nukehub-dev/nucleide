//! MontePy-L3 semantic objects: `MODE`, `TRn`, `U`/`LAT`/`FILL`, `IMP`/`VOL`,
//! and a minimal typed tally model.
//!
//! The byte-exact round-trip core ([`crate::problem`], [`crate::cell`],
//! [`crate::surf`]) is untouched: this module only *reads* the existing
//! [`CellCard`](crate::cell::CellCard) params,
//! [`SurfCard`](crate::surf::SurfCard) pointers, and
//! [`DataCard`](crate::problem::DataCard) cards into typed views, and
//! [`validate_problem`] cross-checks them. [`DeckProblem`](crate::problem::DeckProblem)
//! re-exports thin delegating accessors so callers never touch this module
//! directly.
//!
//! # Ported vs nucleide-defined behaviour
//!
//! Ported MontePy semantics (behaviour matches the reference reader):
//!
//! - `MODE <particles>` has no number and defaults to `{N}` when absent; a
//!   second `MODE` card is an error. All 37 particle shorthands are accepted
//!   (`N P E | Q U V F H L + - X Y O ! < > G / Z K % ^ B _ ~ C W @ D T S A * ? #`).
//! - `TRn ox oy oz [B1..B9] [1|-1]`; `*TRn` gives rotation in degrees.
//!   Displacement is exactly 3 entries (short cards are errors); rotation is
//!   absent or has 5--9 entries; a 13th entry of `1`/`-1` selects the
//!   main-to-aux direction. Inline `FILL` transforms are hidden (numberless)
//!   pass-through objects.
//! - Universes have no card of their own and are numbered `>= 0`. Cell
//!   `U=n` assigns, `U=-n` additionally marks the cell not truncated by its
//!   parent. Data-block `U` cards hold one entry per cell in cell order with
//!   `J` jumps; missing universes (including 0) are auto-created. A second
//!   `U` card, or cell-level and data-level `U` for the same cell, is an error.
//! - `LAT` is `1` (rectangular/hexahedral) or `2` (hexagonal) only, in cells
//!   and in the data block; anything else is an error.
//! - Cell `FILL n [(m|(coords))]` fills with one universe; `FILL i1:i2 j1:j2
//!   k1:k2 u... [(...)]` fills a lattice with a 3-D universe matrix (`0`
//!   means empty). `(m)` references transform `m`, `(coords)` hides an inline
//!   transform. `*FILL` means degrees. Data-block `FILL` cards hold only a
//!   simple per-cell universe list. Fill references to universes or
//!   transforms that do not exist are link errors (the reference reader leaks
//!   a raw key error there).
//! - `TRCL` has no semantic object (kept as a raw passthrough token).
//! - Surface trailing pointers: a positive int is a transform link, a
//!   negative int names the periodic partner surface.
//! - Tallies have no semantic objects in the reference reader (`F`/`FM` fall
//!   to generic data, `DE`/`SDEF`/`FMESH` are refused); the model here is
//!   greenfield (see [`TallyView`]).
//! - Validation mirrors the reference checks: duplicate numbers conflict;
//!   dangling materials, surfaces, complements, periodic partners, and
//!   transforms are malformed input; cell+data redundant definitions and
//!   duplicate `U`/`LAT`/`FILL`/`VOL` cards are malformed input; write-time
//!   state (density/material pairing, empty geometry, surface numbers and
//!   constants, non-empty materials, 3-entry displacement) is enforced.
//!   Particle/mode mismatches are notes, not errors.
//!
//! Nucleide-defined rules (no reference equivalent; enforced by
//! [`validate_problem`] and documented here):
//!
//! - `LAT` without `FILL`, and a `FILL` matrix without `LAT`, on the same
//!   cell are errors. (Single-universe `FILL` without `LAT` is legal.)
//! - Tally numbers must end in a valid type digit (`1, 2, 4, 5, 6, 7, 8`);
//!   `FMn`/`En` cards without a matching `Fn` card are errors.
//! - Data-block `J`/`nJ` jumps and `nR` repeats / `nM` multiplies are expanded
//!   (`nM` multiplies the previous entry by `n`); `nI` interpolation is
//!   rejected with a message.
//! - A data-block `IMP` card needs a particle classifier (`IMP:N ...`); a
//!   bare `VOL` cell keyword means "MCNP calculates it".

use std::collections::{BTreeMap, BTreeSet};

use crate::cell::{CellCard, GeomExpr};
use crate::inp::{Error, McnpMaterial};
use crate::problem::DataCard;
use crate::surf::SurfCard;

/// MCNP particle shorthands (37 single-character tokens).
pub const PARTICLE_SHORTHANDS: [&str; 37] = [
    "N", "P", "E", "|", "Q", "U", "V", "F", "H", "L", "+", "-", "X", "Y", "O", "!", "<", ">", "G",
    "/", "Z", "K", "%", "^", "B", "_", "~", "C", "W", "@", "D", "T", "S", "A", "*", "?", "#",
];

/// True when `token` names an MCNP particle (case-insensitive).
pub fn is_particle_token(token: &str) -> bool {
    PARTICLE_SHORTHANDS.contains(&token.to_ascii_uppercase().as_str())
}

/// Split data-card name into `(prefix, number, classifier, degrees)`.
///
/// The name is the uppercased first token (`F4:N`, `*TR1`, `IMP:N,P`, ...).
/// The prefix is the leading alphabetic run, the number the trailing digits,
/// the classifier the `:` suffix split on `,`, and degrees the leading `*`.
pub struct CardName {
    /// Leading alphabetic run (`F`, `FM`, `TR`, `MODE`, ...).
    pub prefix: String,
    /// Trailing digits, if any.
    pub number: Option<u32>,
    /// `:` suffix split on `,` (particle classifiers).
    pub classifier: Vec<String>,
    /// Leading `*` (`*TRn`, `*FILL`).
    pub degrees: bool,
}

/// Split a data-card name; see [`CardName`].
pub fn split_card_name(name: &str) -> CardName {
    let upper = name.to_ascii_uppercase();
    let (head, classifier) = match upper.split_once(':') {
        Some((h, c)) => (
            h,
            c.split(',')
                .filter(|s| !s.is_empty())
                .map(str::to_string)
                .collect(),
        ),
        None => (upper.as_str(), Vec::new()),
    };
    let (degrees, head) = match head.strip_prefix('*') {
        Some(rest) => (true, rest),
        None => (false, head),
    };
    let split = head
        .char_indices()
        .find(|(_, c)| !c.is_ascii_alphabetic())
        .map(|(i, _)| i)
        .unwrap_or(head.len());
    let (alpha, digits) = head.split_at(split);
    CardName {
        prefix: alpha.to_string(),
        number: if digits.is_empty() {
            None
        } else {
            digits.parse::<u32>().ok()
        },
        classifier,
        degrees,
    }
}

/// One expanded entry of a data-block value list: a jump (`J`) or a literal.
#[derive(Debug, Clone, PartialEq)]
pub enum Slot {
    /// A `J`/`nJ` jump: keep the default for this position.
    Jump,
    /// A literal (or `nR`/`nM`-expanded) token.
    Text(String),
}

/// Expand `J`/`nJ` jumps, `R`/`nR` repeats, and `M`/`nM` multiplies.
///
/// `nM` multiplies the previous entry by `n`; repeats/multiplies with no
/// previous entry (or right after a jump) are errors, as is `nI`
/// interpolation, which has no semantic representation here.
pub fn expand_shortcuts(tokens: &[String], line: usize, context: &str) -> Result<Vec<Slot>, Error> {
    let mut out: Vec<Slot> = Vec::new();
    let mut last: Option<String> = None;
    for token in tokens {
        let upper = token.to_ascii_uppercase();
        if upper == "J" {
            out.push(Slot::Jump);
            last = None;
            continue;
        }
        if let Some(count) = strip_count_suffix(&upper, 'J') {
            for _ in 0..count {
                out.push(Slot::Jump);
            }
            last = None;
            continue;
        }
        if upper == "R" || is_count_suffix(&upper, 'R') {
            let count = if upper == "R" {
                1
            } else {
                upper[..upper.len() - 1]
                    .parse::<usize>()
                    .map_err(|_| bad(line, format!("cannot parse repeat `{token}` on {context}")))?
            };
            let prev = last.clone().ok_or_else(|| {
                bad(
                    line,
                    format!("repeat `{token}` follows a jump on {context}"),
                )
            })?;
            for _ in 0..count {
                out.push(Slot::Text(prev.clone()));
            }
            last = Some(prev);
            continue;
        }
        if let Some(factor) = strip_float_suffix(&upper, 'M') {
            let prev = last.clone().ok_or_else(|| {
                bad(
                    line,
                    format!("multiply `{token}` follows a jump on {context}"),
                )
            })?;
            out.push(Slot::Text(multiply_text(&prev, factor, line, token)?));
            last = out.last().and_then(|s| match s {
                Slot::Text(t) => Some(t.clone()),
                Slot::Jump => None,
            });
            continue;
        }
        if upper == "M" || is_count_suffix(&upper, 'M') {
            return Err(bad(
                line,
                format!("cannot parse multiply `{token}` on {context}"),
            ));
        }
        if upper.contains('I') && is_interpolate(token) {
            return Err(bad(
                line,
                format!("interpolation `{token}` is not supported on {context}"),
            ));
        }
        out.push(Slot::Text(token.clone()));
        last = Some(token.clone());
    }
    Ok(out)
}

/// `nX` with an integer count, else `None`.
fn strip_count_suffix(upper: &str, letter: char) -> Option<usize> {
    if is_count_suffix(upper, letter) {
        upper[..upper.len() - 1].parse::<usize>().ok()
    } else {
        None
    }
}

fn is_count_suffix(upper: &str, letter: char) -> bool {
    upper.len() > 1
        && upper.ends_with(letter)
        && upper[..upper.len() - 1].chars().all(|c| c.is_ascii_digit())
}

/// `nX` with a float factor (`2M`, `0.5M`), else `None`.
fn strip_float_suffix(upper: &str, letter: char) -> Option<f64> {
    if upper.len() > 1 && upper.ends_with(letter) {
        upper[..upper.len() - 1].parse::<f64>().ok()
    } else {
        None
    }
}

/// `prev * factor`, keeping integer spelling when both sides are integers.
fn multiply_text(prev: &str, factor: f64, line: usize, token: &str) -> Result<String, Error> {
    if let (Ok(a), Some(b)) = (prev.parse::<i64>(), as_integer_factor(factor)) {
        return Ok((a * b).to_string());
    }
    let a: f64 = prev.parse().map_err(|_| {
        bad(
            line,
            format!("cannot multiply non-numeric `{prev}` with `{token}`"),
        )
    })?;
    Ok(trim_float(a * factor))
}

/// Integer factor when the float is integral, else `None`.
fn as_integer_factor(factor: f64) -> Option<i64> {
    (factor == factor.trunc()).then_some(factor as i64)
}

/// Shortest rendering of a float (`2.5`, `10`).
fn trim_float(v: f64) -> String {
    if v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

/// Heuristic for interpolate shortcuts (`2I`, `3ILOG`, `2LOG`).
fn is_interpolate(token: &str) -> bool {
    let upper = token.to_ascii_uppercase();
    for suffix in ["ILOG", "LOG", "I"] {
        if upper.len() > suffix.len() && upper.ends_with(suffix) {
            let head = &upper[..upper.len() - suffix.len()];
            if !head.is_empty()
                && head
                    .chars()
                    .all(|c| c.is_ascii_digit() || c == '.' || c == '+' || c == '-')
            {
                return true;
            }
        }
    }
    false
}

fn bad(line: usize, message: String) -> Error {
    Error::BadGeometry { line, message }
}

/// Full token stream of a data card: every `raw_lines` entry with the `$`
/// comment stripped, minus the card name itself (continuations included —
/// [`DataCard::args`] only covers the first line).
pub fn card_tokens(card: &DataCard) -> Vec<String> {
    let mut tokens = Vec::new();
    for line in &card.raw_lines {
        let code = line.split('$').next().unwrap_or("");
        tokens.extend(code.split_whitespace().map(str::to_string));
    }
    if !tokens.is_empty() {
        tokens.remove(0);
    }
    tokens
}

// ---------------------------------------------------------------------------
// MODE
// ---------------------------------------------------------------------------

/// Typed `MODE` card: the problem particle set (default `{N}`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModeView {
    /// Sorted uppercase particle shorthands.
    pub particles: Vec<String>,
}

/// Parse the `MODE` card (at most one per problem).
pub fn parse_mode(data: &[DataCard]) -> Result<ModeView, Error> {
    let mut found: Option<&DataCard> = None;
    for card in data {
        if card.name.eq_ignore_ascii_case("MODE") {
            if found.is_some() {
                return Err(bad(
                    card.line,
                    "duplicate MODE card (only one MODE card per problem)".to_string(),
                ));
            }
            found = Some(card);
        }
    }
    let Some(card) = found else {
        return Ok(ModeView {
            particles: vec!["N".to_string()],
        });
    };
    let mut set = BTreeSet::new();
    for token in card_tokens(card) {
        let particle = token.to_ascii_uppercase();
        if !is_particle_token(&particle) {
            return Err(bad(
                card.line,
                format!("unknown particle `{token}` on MODE card"),
            ));
        }
        set.insert(particle);
    }
    Ok(ModeView {
        particles: set.into_iter().collect(),
    })
}

// ---------------------------------------------------------------------------
// Transforms
// ---------------------------------------------------------------------------

/// Typed `TRn` card (or hidden inline `FILL` transform).
#[derive(Debug, Clone, PartialEq)]
pub struct TransformView {
    /// Transform number (`0` for hidden inline transforms).
    pub number: u32,
    /// `*TRn`: rotation is in degrees, not cosines.
    pub is_in_degrees: bool,
    /// Displacement vector (exactly 3 entries).
    pub displacement: [f64; 3],
    /// Rotation matrix (empty or 5--9 entries).
    pub rotation: Vec<f64>,
    /// Trailing `1` (true) vs `-1` (false); defaults to true.
    pub is_main_to_aux: bool,
    /// Inline `FILL` transform with no card of its own.
    pub hidden: bool,
    /// 1-based source line (`0` for hidden transforms).
    pub line: usize,
}

/// Parse every `TRn` data card.
pub fn parse_transforms(data: &[DataCard]) -> Result<Vec<TransformView>, Error> {
    let mut out = Vec::new();
    let mut seen: BTreeMap<u32, usize> = BTreeMap::new();
    for card in data {
        let name = split_card_name(&card.name);
        if name.prefix != "TR" {
            continue;
        }
        let Some(number) = name.number else { continue };
        if number < 1 {
            return Err(bad(
                card.line,
                format!("invalid transform number `{}`", card.name),
            ));
        }
        if let Some(_first) = seen.insert(number, card.line) {
            return Err(Error::DuplicateNumber {
                kind: "transform",
                number,
            });
        }
        let tokens = card_tokens(card);
        let mut entries = Vec::with_capacity(tokens.len());
        for token in &tokens {
            entries.push(token.parse::<f64>().map_err(|_| {
                bad(
                    card.line,
                    format!("cannot parse `{token}` as a number on TR{number} card"),
                )
            })?);
        }
        out.push(build_transform(
            number,
            name.degrees,
            &entries,
            card.line,
            false,
            &format!("TR{number} card"),
        )?);
    }
    Ok(out)
}

/// Shared `TRn`/hidden-transform entry validation.
fn build_transform(
    number: u32,
    is_in_degrees: bool,
    entries: &[f64],
    line: usize,
    hidden: bool,
    context: &str,
) -> Result<TransformView, Error> {
    if entries.len() < 3 {
        return Err(bad(
            line,
            format!(
                "{context} needs at least 3 displacement entries, found {}",
                entries.len()
            ),
        ));
    }
    if entries.len() > 13 {
        return Err(bad(
            line,
            format!("{context} has too many entries, found {}", entries.len()),
        ));
    }
    let displacement = [entries[0], entries[1], entries[2]];
    let (rotation, is_main_to_aux) = if entries.len() == 13 {
        let flag = entries[12];
        if flag != 1.0 && flag != -1.0 {
            return Err(bad(
                line,
                format!("{context} trailing entry must be 1 or -1, found `{flag}`"),
            ));
        }
        (entries[3..12].to_vec(), flag > 0.0)
    } else {
        (entries[3..].to_vec(), true)
    };
    if !rotation.is_empty() && (rotation.len() < 5 || rotation.len() > 9) {
        return Err(bad(
            line,
            format!(
                "{context} rotation needs 5-9 entries, found {}",
                rotation.len()
            ),
        ));
    }
    Ok(TransformView {
        number,
        is_in_degrees,
        displacement,
        rotation,
        is_main_to_aux,
        hidden,
        line,
    })
}

// ---------------------------------------------------------------------------
// Cell modifiers (U/LAT/FILL/IMP/VOL/TRCL params + data cards)
// ---------------------------------------------------------------------------

/// What a cell `FILL` holds: one universe or a 3-D lattice matrix.
#[derive(Debug, Clone, PartialEq)]
pub enum FillTarget {
    /// `FILL n`.
    Single(u32),
    /// `FILL i1:i2 j1:j2 k1:k2 u...` (`0` entries are empty).
    Matrix {
        /// Minimum indices `[i, j, k]`.
        min_index: [i32; 3],
        /// Maximum indices `[i, j, k]`.
        max_index: [i32; 3],
        /// Universes in `k, j, i` order with `i` fastest (`None` = `0`).
        universes: Vec<Option<u32>>,
    },
}

/// Transform carried by a cell `FILL`.
#[derive(Debug, Clone, PartialEq)]
pub enum FillTransform {
    /// `FILL n (m)`: reference to `TRm`.
    Reference(u32),
    /// `FILL n (coords...)`: hidden inline transform.
    Hidden(TransformView),
}

/// Parsed cell `FILL` (before universe/transform link checks).
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedFill {
    /// Single universe or lattice matrix.
    pub target: FillTarget,
    /// Optional transform reference or hidden transform.
    pub transform: Option<FillTransform>,
    /// `*FILL`: hidden rotation is in degrees.
    pub in_degrees: bool,
}

/// All claimed cell parameters of one cell card.
#[derive(Debug, Clone, PartialEq)]
pub struct CellParams {
    /// Cell number.
    pub cell: u32,
    /// 1-based source line.
    pub line: usize,
    /// `U=n` with the `U=-n` no-truncate flag.
    pub universe: Option<(u32, bool)>,
    /// `LAT=n`.
    pub lattice: Option<u8>,
    /// `FILL ...`.
    pub fill: Option<ParsedFill>,
    /// `(particle, value)` per `IMP:p=...` entry.
    pub importances: Vec<(String, f64)>,
    /// `VOL=v`, or bare `VOL` (`Some(None)` = MCNP calculates it).
    pub volume: Option<Option<f64>>,
    /// Raw `TRCL=...` value (no semantic object exists for it).
    pub trcl: Option<String>,
}

/// Bare cell keywords that take no value (shared with [`crate::cell`]).
const BARE_KEYWORDS: [&str; 4] = ["vol", "pwt", "ext", "fcl"];

/// Parse every cell card's claimed parameters.
pub fn parse_cell_params(cells: &[CellCard]) -> Result<Vec<CellParams>, Error> {
    cells.iter().map(parse_one_cell_params).collect()
}

fn parse_one_cell_params(cell: &CellCard) -> Result<CellParams, Error> {
    let mut out = CellParams {
        cell: cell.num,
        line: cell.line,
        universe: None,
        lattice: None,
        fill: None,
        importances: Vec::new(),
        volume: None,
        trcl: None,
    };
    let mut i = 0;
    while i < cell.params.len() {
        let token = &cell.params[i];
        match token.split_once('=') {
            Some((key, first)) => {
                let (degrees, key) = match key.strip_prefix('*') {
                    Some(rest) => (true, rest),
                    None => (false, key),
                };
                let (key, classifier) = match key.split_once(':') {
                    Some((k, c)) => (k, Some(c)),
                    None => (key, None),
                };
                match key.to_ascii_lowercase().as_str() {
                    "u" => {
                        reject_classifier(cell, classifier, "U")?;
                        out.universe = Some(parse_universe_value(cell, first)?);
                    }
                    "lat" => {
                        reject_classifier(cell, classifier, "LAT")?;
                        out.lattice = Some(parse_lattice_value(cell, first)?);
                    }
                    "fill" => {
                        reject_classifier(cell, classifier, "FILL")?;
                        let (value_tokens, next) = spanning_tokens(&cell.params, i, first);
                        out.fill = Some(parse_cell_fill(cell, value_tokens, degrees)?);
                        i = next;
                        continue;
                    }
                    "trcl" => {
                        reject_classifier(cell, classifier, "TRCL")?;
                        let (value_tokens, next) = spanning_tokens(&cell.params, i, first);
                        out.trcl = Some(value_tokens.join(" "));
                        i = next;
                        continue;
                    }
                    "imp" => {
                        let particles = parse_imp_classifier(cell, classifier)?;
                        let value: f64 = first.parse().map_err(|_| {
                            bad(
                                cell.line,
                                format!(
                                    "cell {} importance must be a number ≥ 0, found `{first}`",
                                    cell.num
                                ),
                            )
                        })?;
                        if value < 0.0 {
                            return Err(bad(
                                cell.line,
                                format!(
                                    "cell {} importance must be a number ≥ 0, found `{first}`",
                                    cell.num
                                ),
                            ));
                        }
                        for particle in particles {
                            out.importances.push((particle, value));
                        }
                    }
                    "vol" => {
                        reject_classifier(cell, classifier, "VOL")?;
                        let value: f64 = first.parse().map_err(|_| {
                            bad(
                                cell.line,
                                format!(
                                    "cell {} volume must be a number ≥ 0.0, found `{first}`",
                                    cell.num
                                ),
                            )
                        })?;
                        if value < 0.0 {
                            return Err(bad(
                                cell.line,
                                format!(
                                    "cell {} volume must be a number ≥ 0.0, found `{first}`",
                                    cell.num
                                ),
                            ));
                        }
                        out.volume = Some(Some(value));
                    }
                    _ => {}
                }
            }
            None => {
                if token.eq_ignore_ascii_case("vol") {
                    out.volume = Some(None);
                }
            }
        }
        i += 1;
    }
    Ok(out)
}

/// Tokens belonging to a multi-word value (`FILL`/`TRCL`): the `=` suffix
/// plus following tokens until the next `key=...` or bare keyword.
fn spanning_tokens(params: &[String], index: usize, first: &str) -> (Vec<String>, usize) {
    let mut tokens = vec![first.to_string()];
    let mut next = index + 1;
    while next < params.len() {
        let token = &params[next];
        if token.contains('=') || BARE_KEYWORDS.iter().any(|k| token.eq_ignore_ascii_case(k)) {
            break;
        }
        tokens.push(token.clone());
        next += 1;
    }
    (tokens, next)
}

fn reject_classifier(cell: &CellCard, classifier: Option<&str>, what: &str) -> Result<(), Error> {
    if classifier.is_some() {
        return Err(bad(
            cell.line,
            format!("cell {} {what} takes no particle classifier", cell.num),
        ));
    }
    Ok(())
}

/// `U=n` / `U=-n` (negative = not truncated by the parent).
fn parse_universe_value(cell: &CellCard, text: &str) -> Result<(u32, bool), Error> {
    let value: i32 = text.parse().map_err(|_| {
        bad(
            cell.line,
            format!(
                "cell {} universe must be an integer ≥ 0, found `{text}`",
                cell.num
            ),
        )
    })?;
    Ok((value.unsigned_abs(), value < 0))
}

/// `LAT=1|2`.
fn parse_lattice_value(cell: &CellCard, text: &str) -> Result<u8, Error> {
    let value: i32 = text.parse().map_err(|_| {
        bad(
            cell.line,
            format!("cell {} LAT must be 1 or 2, found `{text}`", cell.num),
        )
    })?;
    if value != 1 && value != 2 {
        return Err(bad(
            cell.line,
            format!("cell {} LAT must be 1 or 2, found `{text}`", cell.num),
        ));
    }
    Ok(value as u8)
}

/// Classifier of `IMP:p[,q]=...` (required; entries must be particles).
fn parse_imp_classifier(cell: &CellCard, classifier: Option<&str>) -> Result<Vec<String>, Error> {
    let Some(classifier) = classifier else {
        return Err(bad(
            cell.line,
            format!(
                "cell {} IMP needs a particle classifier (imp:n=...)",
                cell.num
            ),
        ));
    };
    let mut particles = Vec::new();
    for part in classifier.split(',') {
        let particle = part.to_ascii_uppercase();
        if !is_particle_token(&particle) {
            return Err(bad(
                cell.line,
                format!("unknown particle `{part}` on cell {} IMP", cell.num),
            ));
        }
        particles.push(particle);
    }
    Ok(particles)
}

/// Parse one cell `FILL` value (single universe or index-range matrix with an
/// optional parenthesized transform).
fn parse_cell_fill(
    cell: &CellCard,
    tokens: Vec<String>,
    in_degrees: bool,
) -> Result<ParsedFill, Error> {
    let joined = tokens.join(" ");
    let error = |message: &str| bad(cell.line, format!("cell {} FILL {message}", cell.num));
    // Split off the single parenthesized transform group, if any.
    let (before, paren) = match joined.find('(') {
        Some(open) => {
            let close = joined
                .rfind(')')
                .ok_or_else(|| error("has unbalanced parentheses"))?;
            if close < open {
                return Err(error("has unbalanced parentheses"));
            }
            let rest = joined[close + 1..].trim();
            if !rest.is_empty() {
                return Err(error(&format!("has trailing text `{rest}`")));
            }
            if joined[open + 1..close].contains('(') {
                return Err(error("has nested parentheses"));
            }
            (
                joined[..open].trim().to_string(),
                Some(joined[open + 1..close].trim().to_string()),
            )
        }
        None => {
            if joined.contains(')') {
                return Err(error("has unbalanced parentheses"));
            }
            (joined.trim().to_string(), None)
        }
    };
    let words: Vec<&str> = before.split_whitespace().collect();
    let target = if words.iter().any(|w| w.contains(':')) {
        parse_fill_matrix(cell, &words)?
    } else {
        if words.len() != 1 || words[0].is_empty() {
            return Err(error("needs exactly one universe outside ranges"));
        }
        let universe: i32 = words[0].parse().map_err(|_| {
            error(&format!(
                "universe must be an integer ≥ 0, found `{}`",
                words[0]
            ))
        })?;
        if universe < 0 {
            return Err(error(&format!(
                "universe must be an integer ≥ 0, found `{}`",
                words[0]
            )));
        }
        FillTarget::Single(universe as u32)
    };
    let transform = match paren {
        Some(group) if !group.is_empty() => Some(parse_fill_transform(cell, &group, in_degrees)?),
        _ => None,
    };
    Ok(ParsedFill {
        target,
        transform,
        in_degrees,
    })
}

/// Maximum total cells in one cell-`FILL` universe matrix.
///
/// A matrix carries one explicit universe entry per cell, so legitimate
/// lattices stay far below this; anything above is a hostile `i:j` range
/// (e.g. `-2147483648:2147483647`, whose width overflows `i32`) and fails
/// with a cap error instead of panicking or allocating gigabytes. Lattices
/// larger than the cap must be built programmatically, not spelled per cell.
const MAX_FILL_MATRIX_CELLS: usize = 1_000_000;

/// `FILL i1:i2 j1:j2 k1:k2 u...` matrix.
fn parse_fill_matrix(cell: &CellCard, words: &[&str]) -> Result<FillTarget, Error> {
    let error = |message: String| bad(cell.line, format!("cell {} FILL {message}", cell.num));
    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut universes: Vec<Option<u32>> = Vec::new();
    for word in words {
        match word.split_once(':') {
            Some((lo, hi)) => {
                let lo: i32 = lo
                    .parse()
                    .map_err(|_| error(format!("bad range `{word}`")))?;
                let hi: i32 = hi
                    .parse()
                    .map_err(|_| error(format!("bad range `{word}`")))?;
                if lo > hi {
                    return Err(error(format!(
                        "minimum {lo} exceeds maximum {hi} in range `{word}`"
                    )));
                }
                ranges.push((lo, hi));
            }
            None => {
                let universe: i32 = word
                    .parse()
                    .map_err(|_| error(format!("universe must be ≥ 0, found `{word}`")))?;
                if universe < 0 {
                    return Err(error(format!("universe must be ≥ 0, found `{word}`")));
                }
                universes.push(if universe == 0 {
                    None
                } else {
                    Some(universe as u32)
                });
            }
        }
    }
    if ranges.len() != 3 {
        return Err(error(format!(
            "matrix needs three i:j ranges, found {}",
            ranges.len()
        )));
    }
    // Widths are computed in i64: `hi - lo` overflows i32 for hostile ranges
    // such as `-2147483648:2147483647`. The total cell count uses checked
    // multiplication and is capped: a matrix is an explicit per-cell universe
    // list, so anything above the cap is a hostile range, not a lattice.
    let mut need: usize = 1;
    for (lo, hi) in &ranges {
        let width = (*hi as i64) - (*lo as i64) + 1;
        debug_assert!(width >= 1, "lo <= hi checked above");
        need = need
            .checked_mul(width as usize)
            .filter(|&n| n <= MAX_FILL_MATRIX_CELLS)
            .ok_or_else(|| {
                error(format!(
                    "matrix {w0} x {w1} x {w2} exceeds the {MAX_FILL_MATRIX_CELLS}-cell cap",
                    w0 = ranges[0].1 as i64 - ranges[0].0 as i64 + 1,
                    w1 = ranges[1].1 as i64 - ranges[1].0 as i64 + 1,
                    w2 = ranges[2].1 as i64 - ranges[2].0 as i64 + 1,
                ))
            })?;
    }
    if universes.len() != need {
        return Err(error(format!(
            "matrix needs {need} universes, found {}",
            universes.len()
        )));
    }
    Ok(FillTarget::Matrix {
        min_index: [ranges[0].0, ranges[1].0, ranges[2].0],
        max_index: [ranges[0].1, ranges[1].1, ranges[2].1],
        universes,
    })
}

/// Parenthesized `FILL` transform: `(m)` references `TRm`, `(coords...)`
/// hides an inline transform (degrees come from the `*FILL` key).
fn parse_fill_transform(
    cell: &CellCard,
    group: &str,
    in_degrees: bool,
) -> Result<FillTransform, Error> {
    let words: Vec<&str> = group.split_whitespace().collect();
    if words.len() == 1 {
        if let Ok(number) = words[0].parse::<u32>() {
            if number >= 1 {
                return Ok(FillTransform::Reference(number));
            }
        }
        return Err(bad(
            cell.line,
            format!(
                "cell {} FILL transform must be a positive integer, found `{}`",
                cell.num, words[0]
            ),
        ));
    }
    let mut entries = Vec::with_capacity(words.len());
    for word in &words {
        entries.push(word.parse::<f64>().map_err(|_| {
            bad(
                cell.line,
                format!("cell {} FILL cannot parse `{word}` as a number", cell.num),
            )
        })?);
    }
    build_transform(
        0,
        in_degrees,
        &entries,
        cell.line,
        true,
        &format!("cell {} FILL hidden transform", cell.num),
    )
    .map(FillTransform::Hidden)
}

// ---------------------------------------------------------------------------
// Universe / lattice / fill / importance / volume views
// ---------------------------------------------------------------------------

/// One auto-created universe with its member cells.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UniverseView {
    /// Universe number (`>= 0`).
    pub number: u32,
    /// Cells assigned to this universe, in file order.
    pub cells: Vec<u32>,
    /// Cells marked `U=-n` (not truncated by the parent).
    pub not_truncated: Vec<u32>,
}

/// One cell `LAT` assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LatticeView {
    /// Cell number.
    pub cell: u32,
    /// `1` = rectangular/hexahedral, `2` = hexagonal.
    pub lattice: u8,
    /// 1-based source line.
    pub line: usize,
}

/// One cell `FILL` assignment with linked transform shape.
#[derive(Debug, Clone, PartialEq)]
pub struct FillView {
    /// Cell number.
    pub cell: u32,
    /// 1-based source line.
    pub line: usize,
    /// Single universe or lattice matrix.
    pub target: FillTarget,
    /// Transform reference or hidden transform.
    pub transform: Option<FillTransform>,
    /// `*FILL`: hidden rotation is in degrees.
    pub in_degrees: bool,
}

/// One cell importance entry.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportanceView {
    /// Cell number.
    pub cell: u32,
    /// Uppercase particle shorthand.
    pub particle: String,
    /// Importance value (`>= 0`).
    pub value: f64,
    /// 1-based source line.
    pub line: usize,
}

/// One manual cell volume.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeView {
    /// Cell number.
    pub cell: u32,
    /// Manual volume (`>= 0`).
    pub volume: f64,
    /// 1-based source line.
    pub line: usize,
}

/// Data-block per-cell entries after shortcut expansion.
struct DataColumn {
    /// Source card.
    line: usize,
    /// Expanded entries (`None` = jump).
    entries: Vec<Option<String>>,
}

/// Collect one data-block per-cell column (`U`, `LAT`, `FILL`, `VOL`).
/// At most one card per problem; duplicates are errors.
fn data_column(
    data: &[DataCard],
    prefix: &str,
    merge_message: &str,
) -> Result<Option<DataColumn>, Error> {
    let mut found: Option<&DataCard> = None;
    for card in data {
        if card.name.eq_ignore_ascii_case(prefix) {
            if found.is_some() {
                return Err(bad(card.line, merge_message.to_string()));
            }
            found = Some(card);
        }
    }
    let Some(card) = found else {
        return Ok(None);
    };
    let slots = expand_shortcuts(&card_tokens(card), card.line, &format!("{prefix} card"))?;
    Ok(Some(DataColumn {
        line: card.line,
        entries: slots
            .into_iter()
            .map(|s| match s {
                Slot::Jump => None,
                Slot::Text(t) => Some(t),
            })
            .collect(),
    }))
}

/// Parse one integer column entry.
fn int_entry(text: &str, line: usize, message: &str) -> Result<i32, Error> {
    text.parse::<i32>()
        .map_err(|_| bad(line, format!("{message}, found `{text}`")))
}

/// Universes with member cells (auto-created from `U` usage, including 0).
pub fn parse_universes(cells: &[CellCard], data: &[DataCard]) -> Result<Vec<UniverseView>, Error> {
    let mods = parse_cell_params(cells)?;
    let column = data_column(data, "U", "Cannot have two universe inputs for the problem")?;
    if let Some(column) = &column {
        if column.entries.len() > cells.len() {
            return Err(bad(
                column.line,
                format!(
                    "U card has more entries ({}) than cells ({})",
                    column.entries.len(),
                    cells.len()
                ),
            ));
        }
        for (index, entry) in column.entries.iter().enumerate() {
            if entry.is_some() && mods[index].universe.is_some() {
                return Err(bad(
                    mods[index].line,
                    format!(
                        "cell {} provided U data when those data were in the data block",
                        mods[index].cell
                    ),
                ));
            }
        }
    }
    // Assignment per cell: (universe, not_truncated).
    let mut assigned: Vec<(u32, bool)> = Vec::with_capacity(cells.len());
    for (index, m) in mods.iter().enumerate() {
        if let Some((universe, notrunc)) = m.universe {
            assigned.push((universe, notrunc));
            continue;
        }
        match column
            .as_ref()
            .and_then(|c| c.entries.get(index))
            .cloned()
            .flatten()
        {
            Some(text) => {
                let value = int_entry(
                    &text,
                    column.as_ref().map(|c| c.line).unwrap_or(0),
                    "U card entries must be integers",
                )?;
                assigned.push((value.unsigned_abs(), value < 0));
            }
            None => assigned.push((0, false)),
        }
    }
    let mut members: BTreeMap<u32, (Vec<u32>, Vec<u32>)> = BTreeMap::new();
    members.entry(0).or_default();
    for (m, (universe, notrunc)) in mods.iter().zip(assigned.iter()) {
        let entry = members.entry(*universe).or_default();
        entry.0.push(m.cell);
        if *notrunc {
            entry.1.push(m.cell);
        }
    }
    Ok(members
        .into_iter()
        .map(|(number, (cells, not_truncated))| UniverseView {
            number,
            cells,
            not_truncated,
        })
        .collect())
}

/// Cell `LAT` assignments (cell level and data block).
pub fn parse_lattices(cells: &[CellCard], data: &[DataCard]) -> Result<Vec<LatticeView>, Error> {
    let mods = parse_cell_params(cells)?;
    let column = data_column(
        data,
        "LAT",
        "Cannot have two lattice inputs for the problem",
    )?;
    if let Some(column) = &column {
        if column.entries.len() > cells.len() {
            return Err(bad(
                column.line,
                format!(
                    "LAT card has more entries ({}) than cells ({})",
                    column.entries.len(),
                    cells.len()
                ),
            ));
        }
        for (index, entry) in column.entries.iter().enumerate() {
            if entry.is_some() && mods[index].lattice.is_some() {
                return Err(bad(
                    mods[index].line,
                    format!(
                        "cell {} provided LAT data when those data were in the data block",
                        mods[index].cell
                    ),
                ));
            }
        }
    }
    let mut out = Vec::new();
    for (index, m) in mods.iter().enumerate() {
        if let Some(lattice) = m.lattice {
            out.push(LatticeView {
                cell: m.cell,
                lattice,
                line: m.line,
            });
            continue;
        }
        let text = column
            .as_ref()
            .and_then(|c| c.entries.get(index))
            .cloned()
            .flatten();
        if let Some(text) = text {
            let line = column.as_ref().map(|c| c.line).unwrap_or(0);
            let value = int_entry(&text, line, "LAT card entries must be 1 or 2")?;
            if value != 1 && value != 2 {
                return Err(bad(
                    line,
                    format!("LAT card entries must be 1 or 2, found `{text}`"),
                ));
            }
            out.push(LatticeView {
                cell: m.cell,
                lattice: value as u8,
                line,
            });
        }
    }
    Ok(out)
}

/// Cell `FILL` assignments (cell level and simple data-block lists).
pub fn parse_fills(cells: &[CellCard], data: &[DataCard]) -> Result<Vec<FillView>, Error> {
    let mods = parse_cell_params(cells)?;
    let column = data_column(data, "FILL", "Cannot have two fill inputs for the problem")?;
    if let Some(column) = &column {
        if column.entries.len() > cells.len() {
            return Err(bad(
                column.line,
                format!(
                    "FILL card has more entries ({}) than cells ({})",
                    column.entries.len(),
                    cells.len()
                ),
            ));
        }
        for (index, entry) in column.entries.iter().enumerate() {
            if entry.is_some() && mods[index].fill.is_some() {
                return Err(bad(
                    mods[index].line,
                    format!(
                        "cell {} provided FILL data when those data were in the data block",
                        mods[index].cell
                    ),
                ));
            }
        }
    }
    let mut out = Vec::new();
    for (index, m) in mods.iter().enumerate() {
        if let Some(fill) = &m.fill {
            out.push(FillView {
                cell: m.cell,
                line: m.line,
                target: fill.target.clone(),
                transform: fill.transform.clone(),
                in_degrees: fill.in_degrees,
            });
            continue;
        }
        let text = column
            .as_ref()
            .and_then(|c| c.entries.get(index))
            .cloned()
            .flatten();
        if let Some(text) = text {
            let line = column.as_ref().map(|c| c.line).unwrap_or(0);
            if text.contains(':') || text.contains('(') || text.contains(')') {
                return Err(bad(
                    line,
                    format!(
                        "data-block FILL takes a simple per-cell universe list, found `{text}`"
                    ),
                ));
            }
            let value = int_entry(
                &text,
                line,
                "FILL card entries must be valid universes (integers ≥ 0)",
            )?;
            if value < 0 {
                return Err(bad(
                    line,
                    format!(
                        "FILL card entries must be valid universes (integers ≥ 0), found `{text}`"
                    ),
                ));
            }
            out.push(FillView {
                cell: m.cell,
                line,
                target: FillTarget::Single(value as u32),
                transform: None,
                in_degrees: false,
            });
        }
    }
    Ok(out)
}

/// Cell importance entries (cell level and data-block `IMP:p` cards).
pub fn parse_importances(
    cells: &[CellCard],
    data: &[DataCard],
) -> Result<Vec<ImportanceView>, Error> {
    let mods = parse_cell_params(cells)?;
    // Data IMP cards keyed by classifier text; the same particle twice is an error.
    let mut data_cards: Vec<(&DataCard, Vec<String>)> = Vec::new();
    for card in data {
        let name = split_card_name(&card.name);
        if name.prefix != "IMP" {
            continue;
        }
        if name.classifier.is_empty() {
            return Err(bad(
                card.line,
                "IMP data card needs a particle classifier (IMP:N ...)".to_string(),
            ));
        }
        for particle in &name.classifier {
            if !is_particle_token(particle) {
                return Err(bad(
                    card.line,
                    format!("unknown particle `{particle}` on IMP card"),
                ));
            }
        }
        data_cards.push((card, name.classifier));
    }
    let mut seen_particles: BTreeSet<String> = BTreeSet::new();
    for (card, classifier) in &data_cards {
        for particle in classifier {
            if !seen_particles.insert(particle.clone()) {
                return Err(bad(
                    card.line,
                    "Cannot have two importance inputs for the same particle type".to_string(),
                ));
            }
        }
    }
    // Single pass: the old `any()` + `find().unwrap()` pair scanned `mods`
    // twice and panicked when the two scans disagreed.
    if let Some(first) = mods.iter().find(|m| !m.importances.is_empty()) {
        if !data_cards.is_empty() {
            return Err(bad(
                first.line,
                format!(
                    "cell {} provided IMP data when those data were in the data block",
                    first.cell
                ),
            ));
        }
    }
    let mut out = Vec::new();
    for m in &mods {
        for (particle, value) in &m.importances {
            out.push(ImportanceView {
                cell: m.cell,
                particle: particle.clone(),
                value: *value,
                line: m.line,
            });
        }
    }
    for (card, classifier) in &data_cards {
        let slots = expand_shortcuts(
            &card_tokens(card),
            card.line,
            &format!("{} card", card.name),
        )?;
        if slots.len() > cells.len() {
            return Err(bad(
                card.line,
                format!(
                    "{} card has more entries ({}) than cells ({})",
                    card.name,
                    slots.len(),
                    cells.len()
                ),
            ));
        }
        for (index, slot) in slots.iter().enumerate() {
            let Slot::Text(text) = slot else { continue };
            let value: f64 = text.parse().map_err(|_| {
                bad(
                    card.line,
                    format!("importances must be ≥ 0, found `{text}`"),
                )
            })?;
            if value < 0.0 {
                return Err(bad(
                    card.line,
                    format!("importances must be ≥ 0, found `{text}`"),
                ));
            }
            for particle in classifier {
                out.push(ImportanceView {
                    cell: mods[index].cell,
                    particle: particle.clone(),
                    value,
                    line: card.line,
                });
            }
        }
    }
    Ok(out)
}

/// Manual cell volumes (cell level and data-block `VOL` cards).
pub fn parse_volumes(cells: &[CellCard], data: &[DataCard]) -> Result<Vec<VolumeView>, Error> {
    let mods = parse_cell_params(cells)?;
    let column = data_column(data, "VOL", "Cannot have two volume inputs for the problem")?;
    // A leading NO disables MCNP calculation; the rest are per-cell values.
    let (column, _no_calc) = match column {
        Some(mut column) => {
            let mut no_calc = false;
            if column
                .entries
                .first()
                .is_some_and(|e| e.as_deref().is_some_and(|t| t.eq_ignore_ascii_case("no")))
            {
                no_calc = true;
                column.entries.remove(0);
            }
            (Some(column), no_calc)
        }
        None => (None, false),
    };
    if let Some(column) = &column {
        if column.entries.len() > cells.len() {
            return Err(bad(
                column.line,
                format!(
                    "VOL card has more entries ({}) than cells ({})",
                    column.entries.len(),
                    cells.len()
                ),
            ));
        }
        for (index, entry) in column.entries.iter().enumerate() {
            if entry.is_some() && mods[index].volume.is_some() {
                return Err(bad(
                    mods[index].line,
                    format!(
                        "cell {} provided VOL data when those data were in the data block",
                        mods[index].cell
                    ),
                ));
            }
        }
    }
    let mut out = Vec::new();
    for (index, m) in mods.iter().enumerate() {
        if let Some(volume) = m.volume {
            if let Some(volume) = volume {
                out.push(VolumeView {
                    cell: m.cell,
                    volume,
                    line: m.line,
                });
            }
            continue;
        }
        let text = column
            .as_ref()
            .and_then(|c| c.entries.get(index))
            .cloned()
            .flatten();
        if let Some(text) = text {
            let line = column.as_ref().map(|c| c.line).unwrap_or(0);
            let value: f64 = text.parse().map_err(|_| {
                bad(
                    line,
                    format!("cell volumes must be numbers ≥ 0.0, found `{text}`"),
                )
            })?;
            if value < 0.0 {
                return Err(bad(
                    line,
                    format!("cell volumes must be numbers ≥ 0.0, found `{text}`"),
                ));
            }
            out.push(VolumeView {
                cell: m.cell,
                volume: value,
                line,
            });
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Tallies (greenfield minimal model)
// ---------------------------------------------------------------------------

/// Valid `F`-tally type digits (there is no `F3` tally).
const TALLY_TYPES: [u32; 7] = [1, 2, 4, 5, 6, 7, 8];

/// Minimal typed tally: an `Fn[:p]` card with its `FMn` multiplier and `En`
/// bins grouped by tally number. `DE`/`DF`/`SDEF`/`FMESH` cards have no
/// semantic objects and stay generic data cards.
#[derive(Debug, Clone, PartialEq)]
pub struct TallyView {
    /// Tally number (`Fn`).
    pub number: u32,
    /// Tally type (units digit: 1, 2, 4, 5, 6, 7, 8).
    pub tally_type: u8,
    /// Particle classifier (`:N`, `:N,P`; empty = problem default).
    pub particles: Vec<String>,
    /// Raw `Fn` entries (cell/surface lists, detector positions, ...).
    pub entries: Vec<String>,
    /// Raw `FMn` multiplier tokens, if any.
    pub fm: Option<Vec<String>>,
    /// Raw `En` energy-bin tokens, if any.
    pub e_bins: Option<Vec<String>>,
    /// 1-based source line of the `Fn` card.
    pub line: usize,
}

/// Group `F`/`FM`/`E` data cards by tally number.
pub fn parse_tallies(data: &[DataCard]) -> Result<Vec<TallyView>, Error> {
    let mut tallies: BTreeMap<u32, TallyView> = BTreeMap::new();
    let mut fm_cards: Vec<&DataCard> = Vec::new();
    let mut e_cards: Vec<&DataCard> = Vec::new();
    for card in data {
        let name = split_card_name(&card.name);
        match name.prefix.as_str() {
            "F" => {
                let Some(number) = name.number else { continue };
                if tallies.contains_key(&number) {
                    return Err(bad(card.line, format!("duplicate F{number} card")));
                }
                let tally_type = number % 10;
                if !TALLY_TYPES.contains(&tally_type) {
                    return Err(bad(
                        card.line,
                        format!(
                            "tally {number} has an invalid type digit (F tallies end in 1, 2, 4, 5, 6, 7, or 8)"
                        ),
                    ));
                }
                for particle in &name.classifier {
                    if !is_particle_token(particle) {
                        return Err(bad(
                            card.line,
                            format!("unknown particle `{particle}` on F{number} card"),
                        ));
                    }
                }
                tallies.insert(
                    number,
                    TallyView {
                        number,
                        tally_type: tally_type as u8,
                        particles: name.classifier.clone(),
                        entries: card_tokens(card),
                        fm: None,
                        e_bins: None,
                        line: card.line,
                    },
                );
            }
            "FM" if name.number.is_some() => fm_cards.push(card),
            "E" if name.number.is_some() => e_cards.push(card),
            "FM" | "E" => {}
            _ => {}
        }
    }
    for card in fm_cards {
        let number = split_card_name(&card.name).number.unwrap_or(0);
        let tally = tallies.get_mut(&number).ok_or_else(|| {
            bad(
                card.line,
                format!("FM{number} has no matching F{number} card"),
            )
        })?;
        if tally.fm.is_some() {
            return Err(bad(card.line, format!("duplicate FM{number} card")));
        }
        tally.fm = Some(card_tokens(card));
    }
    for card in e_cards {
        let number = split_card_name(&card.name).number.unwrap_or(0);
        let tally = tallies.get_mut(&number).ok_or_else(|| {
            bad(
                card.line,
                format!("E{number} has no matching F{number} card"),
            )
        })?;
        if tally.e_bins.is_some() {
            return Err(bad(card.line, format!("duplicate E{number} card")));
        }
        let tokens = card_tokens(card);
        for token in &tokens {
            let upper = token.to_ascii_uppercase();
            let shortcut = upper == "J"
                || strip_count_suffix(&upper, 'J').is_some()
                || upper == "R"
                || is_count_suffix(&upper, 'R')
                || strip_float_suffix(&upper, 'M').is_some()
                || upper == "M";
            if !shortcut && token.parse::<f64>().is_err() {
                return Err(bad(
                    card.line,
                    format!("E{number} bins must be numbers, found `{token}`"),
                ));
            }
        }
        tally.e_bins = Some(tokens);
    }
    Ok(tallies.into_values().collect())
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

/// Collect surface numbers referenced by a cell geometry: `(surfaces, cells)`.
/// Complements (`#n`) reference cells; everything else references surfaces.
fn geom_refs(geom: &GeomExpr, surfs: &mut Vec<i32>, cells: &mut Vec<i32>) {
    match geom {
        GeomExpr::HalfSpace(h) => surfs.push(h.surf),
        GeomExpr::Intersect(parts) => {
            for part in parts {
                geom_refs(part, surfs, cells);
            }
        }
        GeomExpr::Union(a, b) => {
            geom_refs(a, surfs, cells);
            geom_refs(b, surfs, cells);
        }
        GeomExpr::Complement(inner) => {
            let mut inner_surfs = Vec::new();
            geom_refs(inner, &mut inner_surfs, cells);
            // Surface slots inside a complement are cell references (`#n`).
            cells.extend(inner_surfs);
        }
    }
}

/// Validate every L3 semantic rule; see the module docs for the ported vs
/// nucleide-defined split. Particle/mode mismatches are not errors — see
/// [`validation_notes_for`].
pub fn validate_problem(
    cells: &[CellCard],
    surfs: &[SurfCard],
    materials: &[McnpMaterial],
    data: &[DataCard],
) -> Result<(), Error> {
    // Duplicate numbers conflict.
    let mut seen: BTreeSet<u32> = BTreeSet::new();
    for cell in cells {
        if !seen.insert(cell.num) {
            return Err(Error::DuplicateNumber {
                kind: "cell",
                number: cell.num,
            });
        }
    }
    seen.clear();
    for surf in surfs {
        if !seen.insert(surf.num) {
            return Err(Error::DuplicateNumber {
                kind: "surface",
                number: surf.num,
            });
        }
    }
    seen.clear();
    for material in materials {
        if !seen.insert(material.number) {
            return Err(Error::DuplicateNumber {
                kind: "material",
                number: material.number,
            });
        }
    }
    // Parsing also reports duplicate MODE/TRn/cards and malformed values.
    parse_mode(data)?;
    let transforms = parse_transforms(data)?;
    let universes = parse_universes(cells, data)?;
    let lattices = parse_lattices(cells, data)?;
    let fills = parse_fills(cells, data)?;
    parse_importances(cells, data)?;
    parse_volumes(cells, data)?;
    parse_tallies(data)?;

    let universe_numbers: BTreeSet<u32> = universes.iter().map(|u| u.number).collect();
    let transform_numbers: BTreeSet<u32> = transforms.iter().map(|t| t.number).collect();
    let cell_numbers: BTreeSet<u32> = cells.iter().map(|c| c.num).collect();
    let surf_numbers: BTreeSet<u32> = surfs.iter().map(|s| s.num).collect();
    let material_numbers: BTreeSet<u32> = materials.iter().map(|m| m.number).collect();

    // Dangling cell references: materials, surfaces, complements.
    for cell in cells {
        if cell.mat != 0 && !material_numbers.contains(&cell.mat) {
            return Err(bad(
                cell.line,
                format!("cell {} references missing material {}", cell.num, cell.mat),
            ));
        }
        let mut surf_refs = Vec::new();
        let mut cell_refs = Vec::new();
        geom_refs(&cell.geom, &mut surf_refs, &mut cell_refs);
        for surf in surf_refs {
            let number = surf.unsigned_abs();
            if !surf_numbers.contains(&number) {
                return Err(bad(
                    cell.line,
                    format!("cell {} references missing surface {number}", cell.num),
                ));
            }
        }
        for complement in cell_refs {
            let number = complement.unsigned_abs();
            if !cell_numbers.contains(&number) {
                return Err(bad(
                    cell.line,
                    format!("cell {} references missing cell {number}", cell.num),
                ));
            }
        }
    }
    // Dangling surface references: transforms and periodic partners.
    for surf in surfs {
        if let Some(number) = surf.transform {
            if !transform_numbers.contains(&number) {
                return Err(bad(
                    surf.line,
                    format!("surface {} references missing transform {number}", surf.num),
                ));
            }
        }
        if let Some(number) = surf.periodic {
            if !surf_numbers.contains(&number) {
                return Err(bad(
                    surf.line,
                    format!(
                        "surface {} references missing periodic surface {number}",
                        surf.num
                    ),
                ));
            }
        }
    }
    // Dangling fill references (universes auto-created above must cover them).
    for fill in &fills {
        let missing_universe = match &fill.target {
            FillTarget::Single(universe) => {
                (!universe_numbers.contains(universe)).then_some(*universe)
            }
            FillTarget::Matrix { universes, .. } => universes
                .iter()
                .flatten()
                .find(|u| !universe_numbers.contains(u))
                .copied(),
        };
        if let Some(universe) = missing_universe {
            return Err(bad(
                fill.line,
                format!(
                    "cell {} fill references missing universe {universe}",
                    fill.cell
                ),
            ));
        }
        if let Some(FillTransform::Reference(number)) = &fill.transform {
            if !transform_numbers.contains(number) {
                return Err(bad(
                    fill.line,
                    format!(
                        "cell {} fill references missing transform {number}",
                        fill.cell
                    ),
                ));
            }
        }
    }
    // Write-time state checks.
    for cell in cells {
        if cell.mat == 0 && cell.dens.is_some() {
            return Err(bad(
                cell.line,
                format!("Cell {} has a density set but no material", cell.num),
            ));
        }
        if cell.mat != 0 && cell.dens.is_none() {
            return Err(bad(
                cell.line,
                format!("Cell {} has a non-void material but no density", cell.num),
            ));
        }
        if matches!(cell.geom, GeomExpr::Intersect(ref parts) if parts.is_empty()) {
            return Err(bad(
                cell.line,
                format!("Cell {} has no geometry defined", cell.num),
            ));
        }
    }
    for surf in surfs {
        if surf.num < 1 {
            return Err(bad(
                surf.line,
                format!("Surface: {} does not have a valid number set", surf.num),
            ));
        }
        if let Some(arity) = surf.kind.arity() {
            if surf.coeffs.len() != arity {
                return Err(bad(
                    surf.line,
                    format!(
                        "Surface: {} does not have all required constants set",
                        surf.num
                    ),
                ));
            }
        }
    }
    for material in materials {
        if material.fractions.is_empty() && material.number != 0 {
            return Err(bad(
                0,
                format!(
                    "Material: {} does not have any components defined",
                    material.number
                ),
            ));
        }
    }
    for transform in &transforms {
        if transform.displacement.len() != 3 {
            return Err(bad(
                transform.line,
                format!(
                    "Transform: {} does not have a valid displacement vector",
                    transform.number
                ),
            ));
        }
    }
    // Nucleide-defined lattice/fill cross-checks.
    let lattice_cells: BTreeMap<u32, u8> = lattices.iter().map(|l| (l.cell, l.lattice)).collect();
    let fill_cells: BTreeMap<u32, &FillView> = fills.iter().map(|f| (f.cell, f)).collect();
    for lattice in &lattices {
        if lattice.lattice != 1 && lattice.lattice != 2 {
            return Err(bad(
                lattice.line,
                format!(
                    "cell {} LAT must be 1 or 2 (nucleide-defined check)",
                    lattice.cell
                ),
            ));
        }
        if !fill_cells.contains_key(&lattice.cell) {
            return Err(bad(
                lattice.line,
                format!(
                    "cell {} has LAT but no FILL (nucleide-defined check)",
                    lattice.cell
                ),
            ));
        }
    }
    for fill in &fills {
        if matches!(fill.target, FillTarget::Matrix { .. })
            && !lattice_cells.contains_key(&fill.cell)
        {
            return Err(bad(
                fill.line,
                format!(
                    "cell {} has a FILL matrix but no LAT (nucleide-defined check)",
                    fill.cell
                ),
            ));
        }
    }
    Ok(())
}

/// Non-fatal validation notes (particle/mode mismatches). Never fails:
/// parse errors surface through [`validate_problem`] instead.
pub fn validation_notes_for(cells: &[CellCard], data: &[DataCard]) -> Vec<String> {
    let mut notes = BTreeSet::new();
    let (mode, importances, tallies) = match (|| -> Result<_, Error> {
        Ok((
            parse_mode(data)?,
            parse_importances(cells, data)?,
            parse_tallies(data)?,
        ))
    })() {
        Ok(parts) => parts,
        Err(_) => return Vec::new(),
    };
    for importance in &importances {
        if !mode.particles.iter().any(|p| p == &importance.particle) {
            notes.insert(format!(
                "IMP:{} on cell {} is not in MODE ({})",
                importance.particle,
                importance.cell,
                mode.particles.join(" ")
            ));
        }
    }
    for tally in &tallies {
        for particle in &tally.particles {
            if !mode.particles.iter().any(|p| p == particle) {
                notes.insert(format!(
                    "F{}:{} is not in MODE ({})",
                    tally.number,
                    particle,
                    mode.particles.join(" ")
                ));
            }
        }
    }
    notes.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::problem::parse_deck;

    fn data_cards(text: &str) -> Vec<DataCard> {
        parse_deck(text).unwrap().data
    }

    fn deck_l3() -> String {
        std::fs::read_to_string(format!(
            "{}/../../fixtures/mcnp/inp/deck_l3.txt",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap()
    }

    #[test]
    fn particle_table_has_37_shorthands() {
        assert_eq!(PARTICLE_SHORTHANDS.len(), 37);
        for must in ["N", "P", "E", "|", "H", "D", "T", "S", "A", "*", "?", "#"] {
            assert!(is_particle_token(must), "{must}");
            assert!(is_particle_token(&must.to_ascii_lowercase()), "{must}");
        }
        assert!(!is_particle_token("Q2"));
        assert!(!is_particle_token(""));
    }

    #[test]
    fn mode_defaults_and_rejects_duplicates() {
        let deck = parse_deck(&deck_l3()).unwrap();
        assert_eq!(
            parse_mode(&deck.data).unwrap().particles,
            vec!["N".to_string(), "P".to_string()]
        );
        let minimal = "msg\ntitle\n1 0 -1\n\n1 so 1\n\n";
        let deck = parse_deck(minimal).unwrap();
        assert_eq!(
            parse_mode(&deck.data).unwrap().particles,
            vec!["N".to_string()]
        );
        let dup = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nmode n\nmode p\n";
        let err = parse_mode(&data_cards(dup)).unwrap_err();
        assert!(err.to_string().contains("duplicate MODE card"));
        let bad_particle = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nmode n q2\n";
        let err = parse_mode(&data_cards(bad_particle)).unwrap_err();
        assert!(err.to_string().contains("unknown particle `q2`"));
    }

    #[test]
    fn transforms_degrees_flag_and_short_rejection() {
        let deck = parse_deck(&deck_l3()).unwrap();
        let transforms = parse_transforms(&deck.data).unwrap();
        assert_eq!(transforms.len(), 2);
        assert_eq!(transforms[0].number, 1);
        assert!(!transforms[0].is_in_degrees);
        assert_eq!(transforms[0].displacement, [0.0, 0.0, 5.0]);
        assert!(transforms[0].rotation.is_empty());
        assert!(transforms[0].is_main_to_aux);
        assert_eq!(transforms[1].number, 2);
        assert!(transforms[1].is_in_degrees);
        assert_eq!(transforms[1].rotation.len(), 9);
        assert!(transforms[1].is_main_to_aux);

        let short = "msg\ntitle\n1 0 -1\n\n1 so 1\n\ntr1 1.0 2.0\n";
        let err = parse_transforms(&data_cards(short)).unwrap_err();
        assert!(err.to_string().contains("at least 3 displacement entries"));
        let partial = "msg\ntitle\n1 0 -1\n\n1 so 1\n\ntr1 0 0 0 1 0\n";
        let err = parse_transforms(&data_cards(partial)).unwrap_err();
        assert!(err.to_string().contains("rotation needs 5-9 entries"));
        let dup = "msg\ntitle\n1 0 -1\n\n1 so 1\n\ntr1 0 0 0\ntr1 1 1 1\n";
        let err = parse_transforms(&data_cards(dup)).unwrap_err();
        assert_eq!(
            err,
            Error::DuplicateNumber {
                kind: "transform",
                number: 1
            }
        );
        let flag = "msg\ntitle\n1 0 -1\n\n1 so 1\n\ntr3 0 0 0 0 0 1 0 1 0 1 0 0 -1\n";
        let transforms = parse_transforms(&data_cards(flag)).unwrap();
        assert!(!transforms[0].is_main_to_aux);
    }

    #[test]
    fn universes_autocreate_and_data_block_lists() {
        let deck = parse_deck(&deck_l3()).unwrap();
        let universes = parse_universes(&deck.cells, &deck.data).unwrap();
        assert_eq!(
            universes.iter().map(|u| u.number).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        let one = universes.iter().find(|u| u.number == 1).unwrap();
        assert_eq!(one.cells, vec![1, 2]);

        // Data-block U with jumps; -U flag marks no-truncate.
        let text = "msg\ntitle\n1 0 -1 u=-1\n2 0 -2\n\n1 so 1\n2 so 2\n\nU J 2\n";
        let deck = parse_deck(text).unwrap();
        let universes = parse_universes(&deck.cells, &deck.data).unwrap();
        assert_eq!(
            universes.iter().map(|u| u.number).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        let two = universes.iter().find(|u| u.number == 2).unwrap();
        assert_eq!(two.cells, vec![2]);
        let one = universes.iter().find(|u| u.number == 1).unwrap();
        assert_eq!(one.not_truncated, vec![1]);
    }

    #[test]
    fn duplicate_and_redundant_universe_cards_rejected() {
        let dup = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nU 1\nU 2\n";
        let deck = parse_deck(dup).unwrap();
        let err = parse_universes(&deck.cells, &deck.data).unwrap_err();
        assert!(err.to_string().contains("Cannot have two universe inputs"));
        let redundant = "msg\ntitle\n1 0 -1 u=1\n\n1 so 1\n\nU 1\n";
        let deck = parse_deck(redundant).unwrap();
        let err = parse_universes(&deck.cells, &deck.data).unwrap_err();
        assert!(err.to_string().contains("in the data block"));
        // nM multiplies the previous entry.
        let text = "msg\ntitle\n1 0 -1\n2 0 -2\n3 0 -3\n\n1 so 1\n2 so 2\n3 so 3\n\nU J 1 2M\n";
        let deck = parse_deck(text).unwrap();
        let universes = parse_universes(&deck.cells, &deck.data).unwrap();
        assert_eq!(
            universes.iter().map(|u| u.number).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
    }

    #[test]
    fn lattices_accept_only_1_and_2() {
        let deck = parse_deck(&deck_l3()).unwrap();
        let lattices = parse_lattices(&deck.cells, &deck.data).unwrap();
        assert_eq!(lattices.len(), 1);
        assert_eq!(lattices[0].cell, 3);
        assert_eq!(lattices[0].lattice, 1);
        let bad = "msg\ntitle\n1 0 -1 lat=3\n\n1 so 1\n\n";
        let deck = parse_deck(bad).unwrap();
        let err = parse_lattices(&deck.cells, &deck.data).unwrap_err();
        assert!(err.to_string().contains("LAT must be 1 or 2"));
    }

    #[test]
    fn fills_single_matrix_and_hidden_transform() {
        let deck = parse_deck(&deck_l3()).unwrap();
        let fills = parse_fills(&deck.cells, &deck.data).unwrap();
        assert_eq!(fills.len(), 1);
        match &fills[0].target {
            FillTarget::Matrix {
                min_index,
                max_index,
                universes,
            } => {
                assert_eq!(*min_index, [0, 0, 0]);
                assert_eq!(*max_index, [1, 0, 0]);
                assert_eq!(universes, &vec![Some(1), Some(1)]);
            }
            other => panic!("expected matrix, got {other:?}"),
        }
        // Hidden inline transform with *FILL degrees.
        let text = "msg\ntitle\n1 0 -1 *fill=1 (1.5 0.0 0.0)\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let fills = parse_fills(&deck.cells, &deck.data).unwrap();
        assert!(fills[0].in_degrees);
        assert!(matches!(fills[0].transform, Some(FillTransform::Hidden(_))));
        // Transform reference.
        let text = "msg\ntitle\n1 0 -1 fill=1 (3)\n2 0 -2 u=1\n\n1 so 1\n2 so 2\n\nTR3 0 0 0\n";
        let deck = parse_deck(text).unwrap();
        let fills = parse_fills(&deck.cells, &deck.data).unwrap();
        assert_eq!(fills[0].transform, Some(FillTransform::Reference(3)));
        // Complex data-block FILL is rejected.
        let text = "msg\ntitle\n1 0 -1\n2 0 -2\n\n1 so 1\n2 so 2\n\nFILL 0:1 1\n";
        let deck = parse_deck(text).unwrap();
        let err = parse_fills(&deck.cells, &deck.data).unwrap_err();
        assert!(err.to_string().contains("simple per-cell universe list"));
    }

    #[test]
    fn hostile_fill_matrix_ranges_fail_with_cap_error() {
        // `hi - lo` overflows i32 here; old code panicked in debug and
        // wrapped in release. Checked i64 arithmetic fails identically in
        // both profiles with the cap error, before any allocation.
        let text = "msg\ntitle\n1 0 -1 fill=-2147483648:2147483647 0:0 0:0 1 1\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err = parse_fills(&deck.cells, &deck.data).unwrap_err();
        assert!(
            err.to_string().contains("exceeds the 1000000-cell cap"),
            "{err}"
        );
        // Moderate ranges whose checked product still exceeds the cap.
        let text = "msg\ntitle\n1 0 -1 fill=0:999 0:999 0:999 1 1\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err = parse_fills(&deck.cells, &deck.data).unwrap_err();
        assert!(
            err.to_string().contains("exceeds the 1000000-cell cap"),
            "{err}"
        );
        // A matrix at exactly the cap still parses (universes then mismatch).
        let text = "msg\ntitle\n1 0 -1 fill=0:999 0:999 0:0 1 1\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err = parse_fills(&deck.cells, &deck.data).unwrap_err();
        assert!(err.to_string().contains("needs 1000000 universes"), "{err}");
    }

    #[test]
    fn tallies_group_fm_and_e_bins() {
        let deck = parse_deck(&deck_l3()).unwrap();
        let tallies = parse_tallies(&deck.data).unwrap();
        assert_eq!(tallies.len(), 1);
        assert_eq!(tallies[0].number, 4);
        assert_eq!(tallies[0].tally_type, 4);
        assert_eq!(tallies[0].particles, vec!["N".to_string()]);
        assert_eq!(tallies[0].entries, vec!["1".to_string(), "2".to_string()]);
        assert!(tallies[0].fm.is_some());
        assert!(tallies[0].e_bins.is_some());
        let orphan = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nF4:N 1\nFM5 1 1 1\n";
        let err = parse_tallies(&data_cards(orphan)).unwrap_err();
        assert!(err.to_string().contains("FM5 has no matching F5 card"));
        let bad_type = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nF3:N 1\n";
        let err = parse_tallies(&data_cards(bad_type)).unwrap_err();
        assert!(err.to_string().contains("invalid type digit"));
    }

    #[test]
    fn shortcut_expansion_matches_reference_rules() {
        // J jumps, nJ multi-jumps, nR repeats, nM multiplies.
        let slots = expand_shortcuts(
            &[
                "J".to_string(),
                "1".to_string(),
                "2R".to_string(),
                "5M".to_string(),
            ],
            1,
            "U card",
        )
        .unwrap();
        let texts: Vec<Option<&str>> = slots
            .iter()
            .map(|s| match s {
                Slot::Jump => None,
                Slot::Text(t) => Some(t.as_str()),
            })
            .collect();
        assert_eq!(
            texts,
            vec![None, Some("1"), Some("1"), Some("1"), Some("5")]
        );
        assert!(expand_shortcuts(&["2R".to_string()], 1, "U card").is_err());
        assert!(expand_shortcuts(&["2I".to_string()], 1, "E4 card").is_err());
    }

    #[test]
    fn l3_fixture_validates_clean() {
        let deck = parse_deck(&deck_l3()).unwrap();
        validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap();
        assert!(validation_notes_for(&deck.cells, &deck.data).is_empty());
    }

    #[test]
    fn validate_catches_dangling_and_duplicates() {
        // Dangling material.
        let text = "msg\ntitle\n1 5 -1.0 -1\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err =
            validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap_err();
        assert!(err.to_string().contains("missing material 5"));
        // Dangling surface.
        let text = "msg\ntitle\n1 0 -9\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err =
            validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap_err();
        assert!(err.to_string().contains("missing surface 9"));
        // Duplicate cells.
        let text = "msg\ntitle\n1 0 -1\n1 0 -2\n\n1 so 1\n2 so 2\n\n";
        let deck = parse_deck(text).unwrap();
        let err =
            validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap_err();
        assert_eq!(
            err,
            Error::DuplicateNumber {
                kind: "cell",
                number: 1
            }
        );
        // LAT without FILL (nucleide-defined).
        let text = "msg\ntitle\n1 0 -1 lat=1\n\n1 so 1\n\n";
        let deck = parse_deck(text).unwrap();
        let err =
            validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap_err();
        assert!(err.to_string().contains("LAT but no FILL"));
        // Particle/mode mismatch is a note, not an error.
        let text = "msg\ntitle\n1 1 -1.0 -1 imp:p=1\n\n1 so 1\n\nmode n\nm1 92235 1.0\n";
        let deck = parse_deck(text).unwrap();
        validate_problem(&deck.cells, &deck.surfs, &deck.materials, &deck.data).unwrap();
        let notes = validation_notes_for(&deck.cells, &deck.data);
        assert_eq!(notes.len(), 1);
        assert!(notes[0].contains("IMP:P"));
    }
}
