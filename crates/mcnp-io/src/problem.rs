//! MCNP full-deck model: parse, edit, and write back input decks.
//!
//! [`parse_deck`] splits a deck into message/title/cell/surface/data
//! blocks, parses cells ([`crate::cell`]), surfaces ([`crate::surf`]) and
//! materials ([`crate::inp`]), and keeps every card's source lines. Data
//! cards that carry no typed model (mode, kcode, tallies, sdef, `read`
//! includes) ride along as [`DataCard`] passthroughs — `read` includes are
//! never followed.
//!
//! [`write_deck`] re-emits cards verbatim from their source lines, so an
//! unmodified deck round-trips byte-identical (blank-line separators and
//! the trailing newline are canonicalized — author fixtures accordingly).
//! Only cards touched through the setters are re-rendered canonically
//! (single spaces, wrapped at 128 columns with five-space continuations).

use std::collections::BTreeMap;
use std::path::Path;

use crate::cell::{parse_cell_line, CellCard, GeomExpr};
use crate::inp::{self, Error, McnpMaterial};
use crate::surf::{parse_surf_line, SurfCard};

/// Maximum MCNP input line length (MCNP 6.2+).
pub const LINE_LENGTH: usize = 128;

/// Continuation indent used when canonical rendering wraps a long card.
pub const CONTINUATION_INDENT: &str = "     ";

/// An untyped data-block card: everything that is not a material card.
/// Carried verbatim for write-back (`mode`, `kcode`, tallies, `sdef`,
/// `read` includes, ...).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataCard {
    /// First whitespace token, uppercased (`MODE`, `M1`, `F4:N`, ...).
    pub name: String,
    /// Remaining first-line tokens with the `$` comment stripped.
    pub args: Vec<String>,
    /// Source lines forming this card (continuations included).
    pub raw_lines: Vec<String>,
}

/// A parsed MCNP input deck with format-preserving write-back.
#[derive(Debug, Clone, PartialEq)]
pub struct DeckProblem {
    /// Message line (first deck line; may be empty).
    pub message: String,
    /// Title card (second deck line).
    pub title: String,
    /// Cell cards in file order.
    pub cells: Vec<CellCard>,
    /// Surface cards in file order.
    pub surfs: Vec<SurfCard>,
    /// Trailing comment/blank lines of the cell block (after the last
    /// cell card), verbatim.
    pub cell_trailer: Vec<String>,
    /// Trailing comment/blank lines of the surface block, verbatim.
    pub surf_trailer: Vec<String>,
    /// Material cards parsed by [`crate::inp`] over the full deck text.
    pub materials: Vec<McnpMaterial>,
    /// Data-block cards in file order (material cards included, verbatim).
    pub data: Vec<DataCard>,
}

impl DeckProblem {
    /// Find a cell by number.
    pub fn cell(&self, num: u32) -> Option<&CellCard> {
        self.cells.iter().find(|c| c.num == num)
    }

    /// Find a cell by number, mutably.
    fn cell_mut(&mut self, num: u32) -> Result<&mut CellCard, Error> {
        self.cells
            .iter_mut()
            .find(|c| c.num == num)
            .ok_or_else(|| Error::BadCard {
                line: 0,
                message: format!("unknown cell {num}"),
            })
    }

    /// Set a cell's density, re-rendering that card canonically.
    pub fn set_cell_density(&mut self, cell: u32, dens: f64) -> Result<(), Error> {
        if !dens.is_finite() {
            return Err(Error::BadCard {
                line: 0,
                message: format!("non-finite density {dens} for cell {cell}"),
            });
        }
        let target = self.cell_mut(cell)?;
        if target.mat == 0 {
            return Err(Error::BadCard {
                line: 0,
                message: format!("void cell {cell} takes no density"),
            });
        }
        target.dens = Some(dens);
        target.raw_lines = vec![render_cell(target)];
        Ok(())
    }

    /// Set a cell's material number, re-rendering that card canonically.
    pub fn set_cell_material(&mut self, cell: u32, mat: u32) -> Result<(), Error> {
        let target = self.cell_mut(cell)?;
        target.mat = mat;
        if mat == 0 {
            target.dens = None;
        }
        target.raw_lines = vec![render_cell(target)];
        Ok(())
    }

    /// Serialize the deck: verbatim card text with canonical blank-line
    /// separators and a single trailing newline.
    pub fn dumps(&self) -> String {
        let mut out = String::new();
        out.push_str(&self.message);
        out.push('\n');
        out.push_str(&self.title);
        out.push('\n');
        // No unconditional separator here: MCNP starts cells immediately
        // after the title; a blank there (when present) rides along as the
        // first card's prefix lines.
        for cell in &self.cells {
            for line in cell.prefix_lines.iter().chain(cell.raw_lines.iter()) {
                out.push_str(line);
                out.push('\n');
            }
        }
        for line in &self.cell_trailer {
            out.push_str(line);
            out.push('\n');
        }
        out.push('\n');
        for surf in &self.surfs {
            for line in surf.prefix_lines.iter().chain(surf.raw_lines.iter()) {
                out.push_str(line);
                out.push('\n');
            }
        }
        for line in &self.surf_trailer {
            out.push_str(line);
            out.push('\n');
        }
        out.push('\n');
        for card in &self.data {
            for line in &card.raw_lines {
                out.push_str(line);
                out.push('\n');
            }
        }
        out
    }
}

/// Render one cell card canonically, wrapping past [`LINE_LENGTH`] with
/// five-space continuations.
fn render_cell(card: &CellCard) -> String {
    let mut first = format!("{}", card.num);
    first.push_str(&format!(" {}", card.mat));
    if let Some(dens) = card.dens {
        first.push_str(&format!(" {dens}"));
    }
    let geom = card.geom.render();
    if !geom.is_empty() {
        first.push_str(&format!(" {geom}"));
    }
    let mut text = first;
    for param in &card.params {
        if text.len() + 1 + param.len() > LINE_LENGTH {
            text.push('\n');
            text.push_str(CONTINUATION_INDENT);
        } else {
            text.push(' ');
        }
        text.push_str(param);
    }
    // Guard against degenerate empty geometry.
    if matches!(card.geom, GeomExpr::Intersect(ref parts) if parts.is_empty()) {
        return text;
    }
    text
}

/// Parse a full MCNP input deck.
pub fn parse_deck(text: &str) -> Result<DeckProblem, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 2 {
        return Err(Error::BadCard {
            line: lines.len() + 1,
            message: "deck needs at least a message and a title line".to_string(),
        });
    }
    let message = lines[0].to_string();
    let title = lines[1].to_string();
    let (cell_lines, surf_lines, data_lines) = split_blocks(&lines)?;

    let mut cells = Vec::new();
    let mut pending: Vec<String> = Vec::new();
    for (lineno, logical, raw) in logical_cards(&cell_lines) {
        if logical.is_empty() {
            // Comments/blanks attach to the following card (or the block
            // trailer) so write-back preserves deck order.
            pending.extend(raw);
            continue;
        }
        let mut card = parse_cell_line(&logical, lineno)?;
        card.raw_lines = raw;
        card.prefix_lines = std::mem::take(&mut pending);
        cells.push(card);
    }
    let cell_trailer = std::mem::take(&mut pending);
    let mut surfs = Vec::new();
    for (lineno, logical, raw) in logical_cards(&surf_lines) {
        if logical.is_empty() {
            pending.extend(raw);
            continue;
        }
        let mut card = parse_surf_line(&logical, lineno)?;
        card.raw_lines = raw;
        card.prefix_lines = std::mem::take(&mut pending);
        surfs.push(card);
    }
    let surf_trailer = std::mem::take(&mut pending);
    let mut data = Vec::new();
    for (_lineno, logical, raw) in logical_cards(&data_lines) {
        data.push(parse_data_card(&logical, raw));
    }
    let materials = inp::materials_from_inp(text)?;
    Ok(DeckProblem {
        message,
        title,
        cells,
        surfs,
        cell_trailer,
        surf_trailer,
        materials,
        data,
    })
}

/// Read an MCNP input file and parse it into a [`DeckProblem`].
pub fn parse_deck_file(path: impl AsRef<Path>) -> Result<DeckProblem, Error> {
    let text = std::fs::read_to_string(path.as_ref()).map_err(|e| Error::Io(e.to_string()))?;
    parse_deck(&text)
}

/// Serialize a [`DeckProblem`] back to MCNP input text.
pub fn write_deck(problem: &DeckProblem) -> String {
    problem.dumps()
}

/// Numbered source lines of one deck block: `(1-based line, text)`.
type Block<'a> = Vec<(usize, &'a str)>;

/// Split deck lines (after message/title) into cell, surface, and data
/// line ranges at blank-line separators. Extra blanks attach to the
/// following block so write-back stays byte-identical.
fn split_blocks<'a>(lines: &'a [&'a str]) -> Result<(Block<'a>, Block<'a>, Block<'a>), Error> {
    // Collect separator positions: the first blank ends cells, the next
    // blank ends surfaces, everything after is data.
    let mut blanks = Vec::new();
    for (idx, line) in lines.iter().enumerate().skip(2) {
        if line.trim().is_empty() {
            blanks.push(idx);
            if blanks.len() == 2 {
                break;
            }
        }
    }
    if blanks.len() < 2 {
        return Err(Error::BadCard {
            line: lines.len() + 1,
            message: "deck needs blank-line separators between cell, surface, and data blocks"
                .to_string(),
        });
    }
    let numbered =
        |range: std::ops::Range<usize>| range.map(|i| (i + 1, lines[i])).collect::<Vec<_>>();
    Ok((
        numbered(2..blanks[0]),
        numbered(blanks[0] + 1..blanks[1]),
        numbered(blanks[1] + 1..lines.len()),
    ))
}

/// Group numbered block lines into logical cards, keeping each card's raw
/// source lines. Comment and blank lines ride along as their own entries
/// so write-back preserves them; continuations (five-plus leading spaces)
/// join the current card.
fn logical_cards(block: &[(usize, &str)]) -> Vec<(usize, String, Vec<String>)> {
    let mut out = Vec::new();
    let mut current: Option<(usize, String, Vec<String>)> = None;
    let flush = |current: &mut Option<(usize, String, Vec<String>)>,
                 out: &mut Vec<(usize, String, Vec<String>)>| {
        if let Some(done) = current.take() {
            out.push(done);
        }
    };
    for &(lineno, line) in block {
        let no_comment = line.split('$').next().unwrap_or("");
        if line.trim().is_empty() || crate::cell::is_comment(no_comment) {
            flush(&mut current, &mut out);
            out.push((lineno, String::new(), vec![line.to_string()]));
            continue;
        }
        if line.starts_with("     ") {
            if let Some((_, ref mut text, ref mut raw)) = current {
                text.push(' ');
                text.push_str(no_comment.trim_start());
                raw.push(line.to_string());
                continue;
            }
        }
        flush(&mut current, &mut out);
        current = Some((
            lineno,
            no_comment.trim().to_string(),
            vec![line.to_string()],
        ));
    }
    flush(&mut current, &mut out);
    out
}

/// Classify one logical data-block entry: comments/blanks pass through
/// with an empty name; anything else takes its first token as the name.
fn parse_data_card(logical: &str, raw: Vec<String>) -> DataCard {
    if logical.is_empty() {
        return DataCard {
            name: String::new(),
            args: Vec::new(),
            raw_lines: raw,
        };
    }
    let mut tokens = logical.split_whitespace();
    let name = tokens.next().unwrap_or("").to_ascii_uppercase();
    DataCard {
        name,
        args: tokens.map(str::to_string).collect(),
        raw_lines: raw,
    }
}

/// Format-preserving write-back lives here so [`DeckProblem`] stays small.
pub fn write_problem(problem: &DeckProblem) -> String {
    problem.dumps()
}

/// Deck-level material lookup by number.
pub fn deck_materials(problem: &DeckProblem) -> &Vec<McnpMaterial> {
    &problem.materials
}

/// Cell inventory of a deck: `(cell, material)` pairs in file order.
pub fn cell_inventory(problem: &DeckProblem) -> BTreeMap<u32, u32> {
    problem.cells.iter().map(|c| (c.num, c.mat)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        std::fs::read_to_string(format!(
            "{}/../../fixtures/mcnp/inp/{name}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap()
    }

    #[test]
    fn minimal_deck_round_trips_byte_identical() {
        let text = fixture("deck_minimal.txt");
        let problem = parse_deck(&text).unwrap();
        assert_eq!(problem.title, "Minimal pin-cell deck");
        assert_eq!(problem.cells.len(), 3);
        assert_eq!(problem.surfs.len(), 3);
        assert_eq!(problem.materials.len(), 2);
        assert_eq!(write_deck(&problem), text);
    }

    #[test]
    fn macrobody_deck_round_trips() {
        let text = fixture("deck_macro.txt");
        let problem = parse_deck(&text).unwrap();
        assert!(problem
            .surfs
            .iter()
            .any(|s| s.kind == crate::surf::SurfKind::Rpp));
        assert!(problem.surfs.iter().any(|s| s.reflecting));
        assert_eq!(write_deck(&problem), text);
    }

    #[test]
    fn complement_and_params_deck_round_trips() {
        let text = fixture("deck_params.txt");
        let problem = parse_deck(&text).unwrap();
        assert!(problem
            .cells
            .iter()
            .any(|c| matches!(c.geom, GeomExpr::Intersect(_))));
        assert_eq!(write_deck(&problem), text);
    }

    #[test]
    fn setters_rewrite_one_card() {
        let text = fixture("deck_minimal.txt");
        let mut problem = parse_deck(&text).unwrap();
        problem.set_cell_density(1, -7.0).unwrap();
        assert_eq!(problem.cell(1).unwrap().dens, Some(-7.0));
        // Only the edited card differs from the original text.
        let rewritten = write_deck(&problem);
        assert!(rewritten.contains("1 1 -7 -1 imp:n=1"));
        assert_ne!(rewritten, text);

        problem.set_cell_material(2, 1).unwrap();
        assert_eq!(problem.cell(2).unwrap().mat, 1);

        assert!(problem.set_cell_density(99, -1.0).is_err());
        assert!(problem.set_cell_density(3, f64::NAN).is_err());
        // Void cells take no density.
        assert!(problem.set_cell_density(3, -1.0).is_err());
    }

    #[test]
    fn void_cell_to_material_gains_density_slot() {
        let text = fixture("deck_minimal.txt");
        let mut problem = parse_deck(&text).unwrap();
        problem.set_cell_material(3, 2).unwrap();
        let card = problem.cell(3).unwrap();
        assert_eq!(card.mat, 2);
        // No density was ever assigned; the slot stays empty until set.
        assert_eq!(card.dens, None);
    }

    #[test]
    fn structural_errors_carry_lines() {
        assert!(parse_deck("only a message").is_err());
        assert!(parse_deck("msg\ntitle\n1 1 -1.0 -1\n").is_err());
        assert!(parse_deck("msg\ntitle\n\n\nmode n\n")
            .unwrap()
            .cells
            .is_empty());
    }
}
