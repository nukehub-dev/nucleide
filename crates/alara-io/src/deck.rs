//! ALARA input-deck parser: block splitting, typed model, and writer.
//!
//! The deck is parsed line by line (see the ALARA users' guide section on
//! input file syntax). Blocks may appear in any order; blank lines and
//! `#`-comments are skipped; `#include <path>` inlines another file between
//! blocks; variable-length blocks are terminated by `end`. Fixed-size blocks
//! ([`Geometry`], [`FluxDef`], [`Truncation`], library paths, …) carry their
//! arguments on the header line with optional continuation lines, while
//! `end`-terminated blocks ([`Dimension`], [`Volumes`], [`MatLoading`],
//! [`Mixture`], schedules, …) collect body lines up to `end`.
//!
//! [`AlaraDeck::parse`] and [`AlaraDeck::from_file`] reject structural
//! problems (unknown blocks, missing `end`, `dimension` together with
//! `volume`/`volumes`). Dangling cross-references (schedules naming unknown
//! fluxes, `mat_loading` naming unknown mixtures, …) are reported by the
//! separate [`AlaraDeck::validate`] step so callers can inspect partial
//! decks first.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use nucleide_nuclei::NuclideId;

/// Block keywords recognized by the deck parser (all lowercase).
///
/// Both `volume` (used by the shipped ALARA sample inputs) and `volumes`
/// (used by the users' guide) are accepted and normalize to [`Volumes`].
pub const KNOWN_BLOCKS: &[&str] = &[
    "convert_lib",
    "cooling",
    "data_library",
    "dimension",
    "dump_file",
    "element_lib",
    "flux",
    "geometry",
    "ignore",
    "impurity",
    "major_radius",
    "material_lib",
    "mat_loading",
    "minor_radius",
    "mixture",
    "output",
    "pulsehistory",
    "ref_flux_type",
    "schedule",
    "skip_zones",
    "solve_zones",
    "spatial_norm",
    "truncation",
    "volume",
    "volumes",
];

/// Variable-length blocks terminated by the keyword `end`.
const END_BLOCKS: &[&str] = &[
    "cooling",
    "dimension",
    "mat_loading",
    "mixture",
    "output",
    "pulsehistory",
    "schedule",
    "skip_zones",
    "solve_zones",
    "spatial_norm",
    "volume",
    "volumes",
];

/// Maximum `#include` nesting depth (cycle guard).
const MAX_INCLUDE_DEPTH: usize = 32;

/// One raw input block in file order.
///
/// `kind` is the lowercase block keyword (`volume` is normalized to
/// `volumes`); `line` is the 1-based line number of the block header;
/// `body` holds the comment-stripped, trimmed argument lines (for fixed-size
/// blocks the header remainder comes first, then continuation lines; for
/// `end`-terminated blocks the lines between the header and `end`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RawBlock {
    /// Lowercase block keyword.
    pub kind: String,
    /// 1-based line number of the block header.
    pub line: usize,
    /// Comment-stripped argument/body lines.
    pub body: Vec<String>,
}

/// `geometry <option>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Geometry {
    /// Lowercase geometry kind (`point`, `rectangular`, `cylindrical`,
    /// `spherical`, or `torus`).
    pub kind: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One `dimension` zone: `<intervals> <upper boundary>`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DimensionZone {
    /// Number of fine-mesh intervals in this zone.
    pub intervals: u32,
    /// Zone upper boundary in cm.
    pub upper: f64,
}

/// `dimension <axis> [<lower>]` block plus `<intervals> <upper>` lines.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Dimension {
    /// Lowercase axis (`x`, `y`, `z`, `r`, `theta`, or `phi`).
    pub axis: String,
    /// First zone lower boundary in cm (inline or first body line).
    pub lower: f64,
    /// Zones in deck order.
    pub zones: Vec<DimensionZone>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One `volume(s)` entry: `<volume> <zone>`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VolumeEntry {
    /// Interval volume in cm³.
    pub volume: f64,
    /// Zone name (matches a `mat_loading` zone).
    pub zone: String,
}

/// `volume`/`volumes` block: interval volumes and zone membership.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Volumes {
    /// Entries in deck order.
    pub entries: Vec<VolumeEntry>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One `mat_loading` entry: `<zone> <mixture>`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MatLoadingEntry {
    /// Zone name.
    pub zone: String,
    /// Mixture name, or `void` for an empty zone.
    pub mixture: String,
}

/// `mat_loading` block: mixture contained in each zone.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MatLoading {
    /// Entries in deck order.
    pub entries: Vec<MatLoadingEntry>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One `mixture` constituent entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum MixtureEntry {
    /// `material <name> <rel_density> <vol_fraction>`: material-library entry.
    Material {
        /// Material name in the material library.
        name: String,
        /// Density relative to the library density.
        rel_density: f64,
        /// Volume fraction of this material in the mixture.
        vol_fraction: f64,
    },
    /// `element <symbol> <rel_density> <vol_fraction>`: element-library
    /// entry; `symbol` may be a plain symbol (`fe`), an isotope (`U235`),
    /// or a modified symbol (`mn:56`).
    Element {
        /// Element or modified chemical symbol.
        symbol: String,
        /// Density relative to the library density.
        rel_density: f64,
        /// Volume fraction of this element in the mixture.
        vol_fraction: f64,
    },
    /// `like <mixture> <rel_density>`: reuse another mixture definition.
    Like {
        /// Referenced mixture name.
        mixture: String,
        /// Relative density scaling.
        rel_density: f64,
    },
    /// `target <element|isotope> <name>`: reverse-calculation target.
    Target {
        /// Either `element` or `isotope`.
        target_kind: String,
        /// Element symbol or `ZZ-AAA` isotope name.
        name: String,
    },
}

/// `mixture <name>` block: composition of one mixture.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mixture {
    /// Symbolic mixture name.
    pub name: String,
    /// Constituents in deck order.
    pub entries: Vec<MixtureEntry>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `flux <name> <file> <scale> <skip> <format>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FluxDef {
    /// Symbolic flux name referenced by schedule items.
    pub name: String,
    /// Flux file path.
    pub file: String,
    /// Uniform flux scaling factor.
    pub scale: f64,
    /// Number of N-group entries to skip before reading.
    pub skip: usize,
    /// Flux file format (currently `default`).
    pub format: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One raw `schedule` item.
///
/// Six tokens (`<optime> <opunit> <flux> <pulse> <delay> <delayunit>`) denote
/// a pulse on a flux; four tokens (`<sub> <pulse> <delay> <delayunit>`)
/// denote a sub-schedule. Tokens are stored verbatim — full schedule
/// expansion (time-unit conversion, hierarchy flattening) is owned by the
/// schedule module, not the deck parser.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScheduleItemRef {
    /// Whitespace-separated tokens of the item line.
    pub tokens: Vec<String>,
    /// 1-based line number of the item.
    pub line: usize,
}

/// `schedule <name>` block: raw irradiation-history items.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScheduleDef {
    /// Symbolic schedule name.
    pub name: String,
    /// Items in deck order.
    pub items: Vec<ScheduleItemRef>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// One `pulsehistory` level: `<pulses> <delay> <unit>`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PulseLevel {
    /// Number of pulses at this level.
    pub pulses: u64,
    /// Delay between pulses in seconds.
    pub delay_s: f64,
}

/// `pulsehistory <name>` block: multi-level pulsing history.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PulseHistory {
    /// Symbolic pulsing-definition name referenced by schedule items.
    pub name: String,
    /// Levels in deck order.
    pub levels: Vec<PulseLevel>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `cooling` block: after-shutdown cooling times.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Cooling {
    /// Cooling times in seconds, in deck order.
    pub times_s: Vec<f64>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `output <resolution>` block: one output definition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OutputDef {
    /// Lowercase resolution (`interval`, `zone`, or `mixture`).
    pub resolution: String,
    /// Modifier lines in deck order (e.g. `number_density`,
    /// `units Ci m3`, `wdr data/NRCA`), comment-stripped verbatim.
    pub entries: Vec<String>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `truncation <tolerance>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Truncation {
    /// Activation-tree truncation tolerance.
    pub tolerance: f64,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `material_lib <path>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaterialLib {
    /// Library file path.
    pub path: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `element_lib <path>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElementLib {
    /// Library file path.
    pub path: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `data_library <type> <file> [<file2>]` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DataLibrary {
    /// Lowercase library type (`alaralib`, `adjlib`, or `eaflib`).
    pub kind: String,
    /// Library file paths (one for `alaralib`/`adjlib`, two for `eaflib`).
    pub files: Vec<String>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `dump_file <path>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DumpFile {
    /// Dump file path.
    pub path: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `major_radius <value>` fixed block (toroidal geometries).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MajorRadius {
    /// Major radius in cm.
    pub value: f64,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `minor_radius <value>` fixed block (toroidal geometries).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MinorRadius {
    /// Minor radius in cm.
    pub value: f64,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `ignore <tolerance>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ignore {
    /// Relative ignore tolerance.
    pub tolerance: f64,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `ref_flux_type <option>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RefFluxType {
    /// Lowercase option (`max` or `volume_avg`).
    pub kind: String,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `impurity <threshold> <tolerance>` fixed block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Impurity {
    /// Relative-concentration impurity threshold.
    pub threshold: f64,
    /// Truncation tolerance used for impurities.
    pub tolerance: f64,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `solve_zones` block: zones solved in this calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolveZones {
    /// Zone names in deck order.
    pub zones: Vec<String>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `skip_zones` block: zones excluded from this calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SkipZones {
    /// Zone names in deck order.
    pub zones: Vec<String>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `spatial_norm` block: per-interval scalar flux normalizations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpatialNorm {
    /// Normalization values in deck order.
    pub values: Vec<f64>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// `convert_lib <from> <to> <files...>` fixed block (library conversion).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConvertLib {
    /// Lowercase source library format.
    pub from: String,
    /// Lowercase target library format.
    pub to: String,
    /// Library files / output basename arguments.
    pub files: Vec<String>,
    /// 1-based line number of the block header.
    pub line: usize,
}

/// A parsed ALARA input deck.
///
/// `blocks` records every block header in file order; the typed views below
/// hold the parsed blocks in matching order. Later single-occurrence blocks
/// overwrite earlier ones in the typed views while `blocks` keeps the full
/// history.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct AlaraDeck {
    /// Raw blocks in file order.
    pub blocks: Vec<RawBlock>,
    /// `geometry` block, if present.
    pub geometry: Option<Geometry>,
    /// `dimension` blocks in file order.
    pub dimensions: Vec<Dimension>,
    /// `volume`/`volumes` block, if present.
    pub volumes: Option<Volumes>,
    /// `mat_loading` block, if present.
    pub mat_loading: Option<MatLoading>,
    /// `mixture` blocks in file order.
    pub mixtures: Vec<Mixture>,
    /// `flux` blocks in file order.
    pub fluxes: Vec<FluxDef>,
    /// `schedule` blocks in file order.
    pub schedules: Vec<ScheduleDef>,
    /// `pulsehistory` blocks in file order.
    pub pulse_histories: Vec<PulseHistory>,
    /// `cooling` block, if present.
    pub cooling: Option<Cooling>,
    /// `output` blocks in file order.
    pub outputs: Vec<OutputDef>,
    /// `truncation` block, if present.
    pub truncation: Option<Truncation>,
    /// `material_lib` block, if present.
    pub material_lib: Option<MaterialLib>,
    /// `element_lib` block, if present.
    pub element_lib: Option<ElementLib>,
    /// `data_library` block, if present.
    pub data_library: Option<DataLibrary>,
    /// `dump_file` block, if present.
    pub dump_file: Option<DumpFile>,
    /// `major_radius` block, if present.
    pub major_radius: Option<MajorRadius>,
    /// `minor_radius` block, if present.
    pub minor_radius: Option<MinorRadius>,
    /// `ignore` block, if present.
    pub ignore: Option<Ignore>,
    /// `ref_flux_type` block, if present.
    pub ref_flux_type: Option<RefFluxType>,
    /// `impurity` block, if present.
    pub impurity: Option<Impurity>,
    /// `solve_zones` block, if present.
    pub solve_zones: Option<SolveZones>,
    /// `skip_zones` block, if present.
    pub skip_zones: Option<SkipZones>,
    /// `spatial_norm` block, if present.
    pub spatial_norm: Option<SpatialNorm>,
    /// `convert_lib` block, if present.
    pub convert_lib: Option<ConvertLib>,
}

/// Convert a time value to seconds.
///
/// Accepts the ALARA single-character units case-insensitively: `s`
/// (second), `m` (minute, 60 s), `h` (hour), `d` (day), `w` (week, 7 days),
/// `y` (year, defined as 365.25 days to match the schedule-module
/// convention — note the users' guide nominally says 52 weeks), and `c`
/// (century, 100 years). Rejects non-finite or negative values and unknown
/// units with [`Error::BadUnits`].
pub fn parse_time_to_seconds(value: f64, unit: char) -> Result<f64> {
    if !value.is_finite() || value < 0.0 {
        return Err(Error::BadUnits(format!(
            "negative or non-finite time `{value}`"
        )));
    }
    const DAY: f64 = 86_400.0;
    const YEAR: f64 = 365.25 * DAY;
    let factor = match unit.to_ascii_lowercase() {
        's' => 1.0,
        'm' => 60.0,
        'h' => 3_600.0,
        'd' => DAY,
        'w' => 7.0 * DAY,
        'y' => YEAR,
        'c' => 100.0 * YEAR,
        _ => return Err(Error::BadUnits(format!("unknown time unit `{unit}`"))),
    };
    Ok(value * factor)
}

// ---------------------------------------------------------------------------
// Internal line handling
// ---------------------------------------------------------------------------

/// One physical input line with its origin for diagnostics and includes.
struct SourcedLine {
    text: String,
    /// 1-based line number within the file this line came from.
    line: usize,
    /// Directory of the file this line came from; `None` for [`AlaraDeck::parse`].
    dir: Option<PathBuf>,
    /// `#include` nesting depth of this line.
    depth: usize,
}

/// One comment-stripped, non-blank logical line.
struct CodeLine {
    text: String,
    line: usize,
}

/// Streaming line source with `#include` expansion.
struct Parser {
    queue: VecDeque<SourcedLine>,
    /// Line number of the most recently popped line (EOF diagnostics).
    last_line: usize,
}

impl Parser {
    fn new(queue: VecDeque<SourcedLine>) -> Self {
        Self {
            queue,
            last_line: 0,
        }
    }

    /// Pop the next meaningful line, skipping blanks and comments.
    ///
    /// When `allow_include` is set, `#include <path>` directives are expanded
    /// inline (resolved relative to the including file); otherwise they are
    /// rejected, since inclusion must not occur within an input block.
    fn pop_code(&mut self, allow_include: bool) -> Result<Option<CodeLine>> {
        while let Some(current) = self.queue.pop_front() {
            self.last_line = current.line;
            let trimmed = current.text.trim();
            if let Some(rel) = match_include(trimmed) {
                if !allow_include {
                    return Err(Error::Parse {
                        line: current.line,
                        msg: "#include must not occur within an input block".to_string(),
                    });
                }
                let dir = current.dir.clone().ok_or_else(|| Error::Parse {
                    line: current.line,
                    msg: "#include requires from_file (no base directory for in-memory text)"
                        .to_string(),
                })?;
                if current.depth >= MAX_INCLUDE_DEPTH {
                    return Err(Error::Parse {
                        line: current.line,
                        msg: format!("#include nesting exceeds {MAX_INCLUDE_DEPTH} (cycle?)"),
                    });
                }
                let target = dir.join(&rel);
                let included = std::fs::read_to_string(&target)?;
                let inc_dir = target.parent().map(Path::to_path_buf).unwrap_or_default();
                let mut prepend: VecDeque<SourcedLine> = included
                    .lines()
                    .enumerate()
                    .map(|(index, text)| SourcedLine {
                        text: text.to_string(),
                        line: index + 1,
                        dir: Some(inc_dir.clone()),
                        depth: current.depth + 1,
                    })
                    .collect();
                prepend.append(&mut self.queue);
                self.queue = prepend;
                continue;
            }
            let code = strip_comment(&current.text).trim().to_string();
            if code.is_empty() {
                continue;
            }
            return Ok(Some(CodeLine {
                text: code,
                line: current.line,
            }));
        }
        Ok(None)
    }

    /// Push a consumed code line back to the front of the queue.
    fn push_front(&mut self, code: CodeLine) {
        self.queue.push_front(SourcedLine {
            text: code.text,
            line: code.line,
            dir: None,
            depth: 0,
        });
    }
}

/// Strip `#`-comments: a `#` whose preceding text is blank (full-line
/// comment) or ends in whitespace (trailing ` #...` after a single-word
/// value) starts a comment running to end of line. A `#` glued to a token
/// (no preceding whitespace) is kept as ordinary text.
fn strip_comment(line: &str) -> &str {
    for (index, ch) in line.char_indices() {
        if ch == '#' {
            let before = &line[..index];
            let prev_is_space = before
                .chars()
                .next_back()
                .is_none_or(|c| c == ' ' || c == '\t');
            if prev_is_space {
                return before;
            }
        }
    }
    line
}

/// Match an `#include <path>` directive on an already-trimmed line.
///
/// Returns the referenced path with surrounding `<…>`, `"…"`, or `'…'`
/// delimiters removed, or `None` when the line is not a directive (including
/// a bare `#include` with no path, which the caller treats as a comment).
fn match_include(trimmed: &str) -> Option<String> {
    let rest = trimmed.strip_prefix("#include")?;
    if rest.is_empty() || !rest.starts_with([' ', '\t']) {
        return None;
    }
    let path = strip_comment(rest).trim();
    let mut path = path;
    for (open, close) in [("<", ">"), ("\"", "\""), ("'", "'")] {
        if path.len() >= 2 && path.starts_with(open) && path.ends_with(close) {
            path = &path[open.len()..path.len() - close.len()];
            break;
        }
    }
    let path = path.trim();
    if path.is_empty() {
        return None;
    }
    Some(path.to_string())
}

/// Canonical block kind: the users' guide spells it `volumes` while shipped
/// sample inputs use `volume`; both map to the same typed block.
fn canonical_kind(keyword: &str) -> &str {
    if keyword == "volume" {
        "volumes"
    } else {
        keyword
    }
}

/// `(min, max)` argument counts for fixed-size blocks (header remainder plus
/// continuation lines).
fn fixed_arity(kind: &str) -> (usize, usize) {
    match kind {
        "geometry" | "truncation" | "material_lib" | "element_lib" | "dump_file"
        | "major_radius" | "minor_radius" | "ignore" | "ref_flux_type" => (1, 1),
        "flux" => (5, 5),
        "data_library" => (2, 3),
        "impurity" => (2, 2),
        // `convert_lib <from> <to> <files...>`: at least the two format
        // names plus one file argument.
        "convert_lib" => (3, usize::MAX),
        _ => (0, 0),
    }
}

/// Split a code line into `(token, line)` pairs sharing the line number.
fn tokens_with_lines(text: &str, line: usize) -> Vec<(String, usize)> {
    text.split_whitespace()
        .map(|token| (token.to_string(), line))
        .collect()
}

fn parse_float(text: &str, line: usize, what: &str) -> Result<f64> {
    let value: f64 = text.parse().map_err(|_| Error::Parse {
        line,
        msg: format!("expected {what}, found `{text}`"),
    })?;
    if !value.is_finite() {
        return Err(Error::Parse {
            line,
            msg: format!("expected finite {what}, found `{text}`"),
        });
    }
    Ok(value)
}

/// Parse `<value> <unit>` with the deck [`parse_time_to_seconds`] helper;
/// `unit` must be a single character.
fn parse_time_pair(value_text: &str, unit_text: &str, line: usize) -> Result<f64> {
    let value = parse_float(value_text, line, "time value")?;
    let mut chars = unit_text.chars();
    let (unit, extra) = (chars.next(), chars.next());
    match (unit, extra) {
        (Some(unit), None) => parse_time_to_seconds(value, unit).map_err(|error| Error::Parse {
            line,
            msg: error.to_string(),
        }),
        _ => Err(Error::Parse {
            line,
            msg: format!("expected single-character time unit, found `{unit_text}`"),
        }),
    }
}

/// Validate an `element` entry symbol or `target element` name: a plain
/// symbol must be a known element (case-insensitive), an isotope-like symbol
/// (containing digits) must parse via [`NuclideId`], and a modified symbol
/// (`ZZ:XXXX`) needs a known-element base with a non-empty suffix.
fn validate_element_symbol(symbol: &str, line: usize) -> Result<()> {
    if symbol.is_empty() {
        return Err(Error::Parse {
            line,
            msg: "expected element symbol, found empty value".to_string(),
        });
    }
    let (base, modified) = match symbol.split_once(':') {
        Some((base, suffix)) => {
            if suffix.is_empty() {
                return Err(Error::Parse {
                    line,
                    msg: format!("malformed modified element symbol `{symbol}`"),
                });
            }
            (base, true)
        }
        None => (symbol, false),
    };
    if base.is_empty() {
        return Err(Error::Parse {
            line,
            msg: format!("malformed element symbol `{symbol}`"),
        });
    }
    if base.chars().any(|c| c.is_ascii_digit()) {
        if modified {
            return Err(Error::Parse {
                line,
                msg: format!("malformed modified element symbol `{symbol}`"),
            });
        }
        NuclideId::from_name(base).map_err(|_| Error::Parse {
            line,
            msg: format!("unknown isotope `{symbol}`"),
        })?;
        return Ok(());
    }
    if !nucleide_nuclei::ELEMENTS
        .iter()
        .any(|element| element.eq_ignore_ascii_case(base))
    {
        return Err(Error::Parse {
            line,
            msg: format!("unknown element `{symbol}`"),
        });
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// AlaraDeck implementation
// ---------------------------------------------------------------------------

impl AlaraDeck {
    /// Parse a deck from in-memory text.
    ///
    /// `#include` directives are rejected here (no base directory); use
    /// [`AlaraDeck::from_file`] for files that use inclusion.
    pub fn parse(text: &str) -> Result<Self> {
        let queue = text
            .lines()
            .enumerate()
            .map(|(index, line)| SourcedLine {
                text: line.to_string(),
                line: index + 1,
                dir: None,
                depth: 0,
            })
            .collect();
        Self::parse_queue(queue)
    }

    /// Read and parse a deck from disk, resolving `#include` paths relative
    /// to the including file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)?;
        let dir = path.parent().map(Path::to_path_buf);
        let queue = text
            .lines()
            .enumerate()
            .map(|(index, line)| SourcedLine {
                text: line.to_string(),
                line: index + 1,
                dir: dir.clone(),
                depth: 0,
            })
            .collect();
        Self::parse_queue(queue)
    }

    fn parse_queue(queue: VecDeque<SourcedLine>) -> Result<Self> {
        let mut parser = Parser::new(queue);
        let mut deck = Self::default();
        let mut saw_block = false;
        while let Some(header) = parser.pop_code(true)? {
            let mut words = header.text.split_whitespace();
            let keyword_raw = words.next().unwrap_or_default();
            let keyword = keyword_raw.to_ascii_lowercase();
            if keyword == "end" {
                return Err(Error::Parse {
                    line: header.line,
                    msg: "unexpected `end` outside a block".to_string(),
                });
            }
            if !KNOWN_BLOCKS.contains(&keyword.as_str()) {
                return Err(Error::UnknownBlock {
                    block: keyword_raw.to_string(),
                    line: header.line,
                });
            }
            saw_block = true;
            let kind = canonical_kind(&keyword).to_string();
            let inline: Vec<(String, usize)> = words
                .map(|token| (token.to_string(), header.line))
                .collect();
            if END_BLOCKS.contains(&keyword.as_str()) {
                let body = Self::collect_body(&mut parser, &keyword, header.line)?;
                deck.push_end_block(&kind, header.line, inline, body)?;
            } else {
                let args = Self::collect_fixed_args(&mut parser, &keyword, header.line, inline)?;
                deck.push_fixed_block(&kind, header.line, args)?;
            }
        }
        if !saw_block {
            return Err(Error::Parse {
                line: 1,
                msg: "empty deck".to_string(),
            });
        }
        if !deck.dimensions.is_empty() && deck.volumes.is_some() {
            return Err(Error::Parse {
                line: parser.last_line,
                msg: "input defines both `dimension` and `volume(s)` blocks; \
                    only one geometry method is allowed"
                    .to_string(),
            });
        }
        Ok(deck)
    }

    /// Collect body lines up to the terminating `end` (exclusive).
    fn collect_body(
        parser: &mut Parser,
        keyword: &str,
        header_line: usize,
    ) -> Result<Vec<CodeLine>> {
        let mut body = Vec::new();
        loop {
            match parser.pop_code(false)? {
                None => {
                    return Err(Error::Parse {
                        line: parser.last_line,
                        msg: format!(
                            "missing `end` for `{keyword}` block opened at line {header_line}"
                        ),
                    });
                }
                Some(line) => {
                    let first = line.text.split_whitespace().next().unwrap_or_default();
                    if first.eq_ignore_ascii_case("end") {
                        return Ok(body);
                    }
                    body.push(line);
                }
            }
        }
    }

    /// Gather fixed-block arguments: the header remainder plus continuation
    /// lines while more tokens are needed. A following block keyword or `end`
    /// stops gathering (missing arguments are an error); anything else is
    /// consumed as continuation tokens.
    fn collect_fixed_args(
        parser: &mut Parser,
        kind: &str,
        header_line: usize,
        mut args: Vec<(String, usize)>,
    ) -> Result<Vec<(String, usize)>> {
        let (need_min, need_max) = fixed_arity(kind);
        while args.len() < need_min {
            match parser.pop_code(false)? {
                None => {
                    return Err(Error::Parse {
                        line: parser.last_line,
                        msg: format!(
                            "`{kind}` block opened at line {header_line} needs at least \
                            {need_min} argument(s), found {}",
                            args.len()
                        ),
                    });
                }
                Some(line) => {
                    let first = line.text.split_whitespace().next().unwrap_or_default();
                    let first_lower = first.to_ascii_lowercase();
                    if first.eq_ignore_ascii_case("end")
                        || KNOWN_BLOCKS.contains(&first_lower.as_str())
                    {
                        if args.len() < need_min {
                            return Err(Error::Parse {
                                line: line.line,
                                msg: format!(
                                    "`{kind}` block opened at line {header_line} needs at least \
                                    {need_min} argument(s), found {}",
                                    args.len()
                                ),
                            });
                        }
                        parser.push_front(line);
                        break;
                    }
                    args.extend(tokens_with_lines(&line.text, line.line));
                }
            }
        }
        if args.len() > need_max {
            return Err(Error::Parse {
                line: header_line,
                msg: format!(
                    "`{kind}` block takes at most {need_max} argument(s), found {}",
                    args.len()
                ),
            });
        }
        Ok(args)
    }

    /// Build a fixed-size block from gathered argument tokens.
    fn push_fixed_block(
        &mut self,
        kind: &str,
        line: usize,
        args: Vec<(String, usize)>,
    ) -> Result<()> {
        let values: Vec<&str> = args.iter().map(|(token, _)| token.as_str()).collect();
        let lines: Vec<usize> = args.iter().map(|(_, num)| *num).collect();
        let arg_line = |index: usize| lines.get(index).copied().unwrap_or(line);
        let raw_body: Vec<String> = if values.is_empty() {
            Vec::new()
        } else {
            vec![values.join(" ")]
        };
        match kind {
            "geometry" => {
                let name = values[0].to_ascii_lowercase();
                if !["point", "rectangular", "cylindrical", "spherical", "torus"]
                    .contains(&name.as_str())
                {
                    return Err(Error::Parse {
                        line,
                        msg: format!("unknown geometry `{}`", values[0]),
                    });
                }
                self.geometry = Some(Geometry { kind: name, line });
            }
            "truncation" => {
                let tolerance = parse_float(values[0], arg_line(0), "truncation tolerance")?;
                self.truncation = Some(Truncation { tolerance, line });
            }
            "flux" => {
                let scale = parse_float(values[2], arg_line(2), "flux scale")?;
                let skip: usize = values[3].parse().map_err(|_| Error::Parse {
                    line: arg_line(3),
                    msg: format!("expected flux skip value, found `{}`", values[3]),
                })?;
                self.fluxes.push(FluxDef {
                    name: values[0].to_string(),
                    file: values[1].to_string(),
                    scale,
                    skip,
                    format: values[4].to_string(),
                    line,
                });
            }
            "material_lib" => {
                self.material_lib = Some(MaterialLib {
                    path: values[0].to_string(),
                    line,
                });
            }
            "element_lib" => {
                self.element_lib = Some(ElementLib {
                    path: values[0].to_string(),
                    line,
                });
            }
            "data_library" => {
                let lib_kind = values[0].to_ascii_lowercase();
                let files: Vec<String> = values[1..].iter().map(ToString::to_string).collect();
                match lib_kind.as_str() {
                    "alaralib" | "adjlib" if files.len() == 1 => {}
                    "eaflib" if files.len() == 2 => {}
                    _ => {
                        return Err(Error::Parse {
                            line,
                            msg: format!(
                                "data_library type `{}` needs {} file(s), found {}",
                                values[0],
                                if lib_kind == "eaflib" { 2 } else { 1 },
                                files.len()
                            ),
                        });
                    }
                }
                if !["alaralib", "adjlib", "eaflib"].contains(&lib_kind.as_str()) {
                    return Err(Error::Parse {
                        line,
                        msg: format!("unknown data library type `{}`", values[0]),
                    });
                }
                self.data_library = Some(DataLibrary {
                    kind: lib_kind,
                    files,
                    line,
                });
            }
            "dump_file" => {
                self.dump_file = Some(DumpFile {
                    path: values[0].to_string(),
                    line,
                });
            }
            "major_radius" => {
                let value = parse_float(values[0], arg_line(0), "major radius")?;
                self.major_radius = Some(MajorRadius { value, line });
            }
            "minor_radius" => {
                let value = parse_float(values[0], arg_line(0), "minor radius")?;
                self.minor_radius = Some(MinorRadius { value, line });
            }
            "ignore" => {
                let tolerance = parse_float(values[0], arg_line(0), "ignore tolerance")?;
                self.ignore = Some(Ignore { tolerance, line });
            }
            "ref_flux_type" => {
                let name = values[0].to_ascii_lowercase();
                if !["max", "volume_avg"].contains(&name.as_str()) {
                    return Err(Error::Parse {
                        line,
                        msg: format!("unknown reference flux type `{}`", values[0]),
                    });
                }
                self.ref_flux_type = Some(RefFluxType { kind: name, line });
            }
            "impurity" => {
                let threshold = parse_float(values[0], arg_line(0), "impurity threshold")?;
                let tolerance = parse_float(values[1], arg_line(1), "impurity tolerance")?;
                self.impurity = Some(Impurity {
                    threshold,
                    tolerance,
                    line,
                });
            }
            "convert_lib" => {
                let from = values[0].to_ascii_lowercase();
                let to = values[1].to_ascii_lowercase();
                for (raw, role) in [(&from, "source"), (&to, "target")] {
                    if !["alaralib", "adjlib", "eaflib"].contains(&raw.as_str()) {
                        return Err(Error::Parse {
                            line,
                            msg: format!("unknown convert_lib {role} format"),
                        });
                    }
                }
                self.convert_lib = Some(ConvertLib {
                    from,
                    to,
                    files: values[2..].iter().map(ToString::to_string).collect(),
                    line,
                });
            }
            other => {
                return Err(Error::Parse {
                    line,
                    msg: format!("internal error: unhandled fixed block `{other}`"),
                });
            }
        }
        self.blocks.push(RawBlock {
            kind: kind.to_string(),
            line,
            body: raw_body,
        });
        Ok(())
    }

    /// Build an `end`-terminated block from header args and body lines.
    fn push_end_block(
        &mut self,
        kind: &str,
        line: usize,
        inline: Vec<(String, usize)>,
        body: Vec<CodeLine>,
    ) -> Result<()> {
        let raw_body: Vec<String> = body.iter().map(|entry| entry.text.clone()).collect();
        let no_inline = |kind: &str| -> Result<()> {
            if !inline.is_empty() {
                return Err(Error::Parse {
                    line,
                    msg: format!("`{kind}` takes no header arguments"),
                });
            }
            Ok(())
        };
        let one_name = |kind: &str| -> Result<String> {
            match inline.as_slice() {
                [(name, _)] => Ok(name.clone()),
                _ => Err(Error::Parse {
                    line,
                    msg: format!("`{kind}` needs exactly one name argument"),
                }),
            }
        };
        let require_entries = |kind: &str| -> Result<()> {
            if body.is_empty() {
                return Err(Error::Parse {
                    line,
                    msg: format!("`{kind}` block has no entries"),
                });
            }
            Ok(())
        };
        match kind {
            "dimension" => {
                let (axis, lower, lower_line) = match inline.as_slice() {
                    [(axis, _)] => (axis.clone(), None, line),
                    [(axis, _), (lower, lower_line)] => {
                        (axis.clone(), Some((lower.clone(), *lower_line)), line)
                    }
                    _ => {
                        return Err(Error::Parse {
                            line,
                            msg: "`dimension` needs `<axis> [<lower>]` header arguments"
                                .to_string(),
                        });
                    }
                };
                let axis = axis.to_ascii_lowercase();
                if !["x", "y", "z", "r", "theta", "phi"].contains(&axis.as_str()) {
                    return Err(Error::Parse {
                        line,
                        msg: format!("unknown dimension axis `{axis}`"),
                    });
                }
                let mut lower_value: Option<f64> = None;
                let mut lower_seen = false;
                if let Some((text, num)) = lower {
                    lower_value = Some(parse_float(&text, num, "dimension lower bound")?);
                    lower_seen = true;
                    let _ = lower_line;
                }
                let mut zones = Vec::new();
                for entry in &body {
                    let tokens: Vec<&str> = entry.text.split_whitespace().collect();
                    if !lower_seen {
                        if tokens.len() != 1 {
                            return Err(Error::Parse {
                                line: entry.line,
                                msg: format!(
                                    "expected dimension lower bound, found `{}`",
                                    entry.text
                                ),
                            });
                        }
                        lower_value =
                            Some(parse_float(tokens[0], entry.line, "dimension lower bound")?);
                        lower_seen = true;
                        continue;
                    }
                    if tokens.len() != 2 {
                        return Err(Error::Parse {
                            line: entry.line,
                            msg: format!("expected `<intervals> <upper>`, found `{}`", entry.text),
                        });
                    }
                    let intervals: u32 = tokens[0].parse().map_err(|_| Error::Parse {
                        line: entry.line,
                        msg: format!("expected interval count, found `{}`", tokens[0]),
                    })?;
                    if intervals == 0 {
                        return Err(Error::Parse {
                            line: entry.line,
                            msg: "dimension interval count must be at least 1".to_string(),
                        });
                    }
                    let upper = parse_float(tokens[1], entry.line, "zone upper boundary")?;
                    zones.push(DimensionZone { intervals, upper });
                }
                let Some(lower) = lower_value else {
                    return Err(Error::Parse {
                        line,
                        msg: "`dimension` block is missing its lower bound".to_string(),
                    });
                };
                if zones.is_empty() {
                    return Err(Error::Parse {
                        line,
                        msg: "`dimension` block has no zones".to_string(),
                    });
                }
                self.dimensions.push(Dimension {
                    axis,
                    lower,
                    zones,
                    line,
                });
            }
            "volumes" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let mut entries = Vec::new();
                for entry in &body {
                    let tokens: Vec<&str> = entry.text.split_whitespace().collect();
                    if tokens.len() != 2 {
                        return Err(Error::Parse {
                            line: entry.line,
                            msg: format!("expected `<volume> <zone>`, found `{}`", entry.text),
                        });
                    }
                    entries.push(VolumeEntry {
                        volume: parse_float(tokens[0], entry.line, "interval volume")?,
                        zone: tokens[1].to_string(),
                    });
                }
                self.volumes = Some(Volumes { entries, line });
            }
            "mat_loading" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let mut entries = Vec::new();
                for entry in &body {
                    let tokens: Vec<&str> = entry.text.split_whitespace().collect();
                    if tokens.len() != 2 {
                        return Err(Error::Parse {
                            line: entry.line,
                            msg: format!("expected `<zone> <mixture>`, found `{}`", entry.text),
                        });
                    }
                    entries.push(MatLoadingEntry {
                        zone: tokens[0].to_string(),
                        mixture: tokens[1].to_string(),
                    });
                }
                self.mat_loading = Some(MatLoading { entries, line });
            }
            "mixture" => {
                let name = one_name(kind)?;
                require_entries(kind)?;
                let mut entries = Vec::new();
                for entry in &body {
                    entries.push(Self::parse_mixture_entry(&entry.text, entry.line)?);
                }
                self.mixtures.push(Mixture {
                    name,
                    entries,
                    line,
                });
            }
            "solve_zones" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let zones = body
                    .iter()
                    .flat_map(|entry| entry.text.split_whitespace().map(str::to_string))
                    .collect();
                self.solve_zones = Some(SolveZones { zones, line });
            }
            "skip_zones" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let zones = body
                    .iter()
                    .flat_map(|entry| entry.text.split_whitespace().map(str::to_string))
                    .collect();
                self.skip_zones = Some(SkipZones { zones, line });
            }
            "spatial_norm" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let mut values = Vec::new();
                for entry in &body {
                    for token in entry.text.split_whitespace() {
                        values.push(parse_float(token, entry.line, "spatial normalization")?);
                    }
                }
                self.spatial_norm = Some(SpatialNorm { values, line });
            }
            "schedule" => {
                let name = one_name(kind)?;
                require_entries(kind)?;
                let mut items = Vec::new();
                for entry in &body {
                    let tokens: Vec<String> =
                        entry.text.split_whitespace().map(str::to_string).collect();
                    match tokens.len() {
                        4 => {
                            parse_time_pair(&tokens[2], &tokens[3], entry.line)?;
                        }
                        6 => {
                            parse_time_pair(&tokens[0], &tokens[1], entry.line)?;
                            parse_time_pair(&tokens[4], &tokens[5], entry.line)?;
                        }
                        _ => {
                            return Err(Error::Parse {
                                line: entry.line,
                                msg: format!(
                                    "expected 4-token sub-schedule or 6-token pulse item, \
                                    found `{}`",
                                    entry.text
                                ),
                            });
                        }
                    }
                    items.push(ScheduleItemRef {
                        tokens,
                        line: entry.line,
                    });
                }
                self.schedules.push(ScheduleDef { name, items, line });
            }
            "pulsehistory" => {
                let name = one_name(kind)?;
                require_entries(kind)?;
                let mut flat: Vec<(String, usize)> = Vec::new();
                for entry in &body {
                    flat.extend(tokens_with_lines(&entry.text, entry.line));
                }
                if flat.len() % 3 != 0 {
                    return Err(Error::Parse {
                        line,
                        msg: format!(
                            "`pulsehistory` needs `<pulses> <delay> <unit>` triplets, \
                            found {} value(s)",
                            flat.len()
                        ),
                    });
                }
                let mut levels = Vec::new();
                for triplet in flat.chunks(3) {
                    let pulses: u64 = triplet[0].0.parse().map_err(|_| Error::Parse {
                        line: triplet[0].1,
                        msg: format!("expected pulse count, found `{}`", triplet[0].0),
                    })?;
                    let delay_s = parse_time_pair(&triplet[1].0, &triplet[2].0, triplet[1].1)?;
                    levels.push(PulseLevel { pulses, delay_s });
                }
                self.pulse_histories
                    .push(PulseHistory { name, levels, line });
            }
            "cooling" => {
                no_inline(kind)?;
                require_entries(kind)?;
                let mut flat: Vec<(String, usize)> = Vec::new();
                for entry in &body {
                    flat.extend(tokens_with_lines(&entry.text, entry.line));
                }
                if flat.len() % 2 != 0 {
                    return Err(Error::Parse {
                        line,
                        msg: format!(
                            "`cooling` needs `<value> <unit>` pairs, found {} value(s)",
                            flat.len()
                        ),
                    });
                }
                let mut times_s = Vec::new();
                for pair in flat.chunks(2) {
                    times_s.push(parse_time_pair(&pair[0].0, &pair[1].0, pair[0].1)?);
                }
                self.cooling = Some(Cooling { times_s, line });
            }
            "output" => {
                let resolution = one_name(kind)?.to_ascii_lowercase();
                if !["interval", "zone", "mixture"].contains(&resolution.as_str()) {
                    return Err(Error::Parse {
                        line,
                        msg: format!("unknown output resolution `{resolution}`"),
                    });
                }
                require_entries(kind)?;
                self.outputs.push(OutputDef {
                    resolution,
                    entries: raw_body.clone(),
                    line,
                });
            }
            other => {
                return Err(Error::Parse {
                    line,
                    msg: format!("internal error: unhandled end block `{other}`"),
                });
            }
        }
        self.blocks.push(RawBlock {
            kind: kind.to_string(),
            line,
            body: raw_body,
        });
        Ok(())
    }

    /// Parse one `mixture` body line into a [`MixtureEntry`].
    fn parse_mixture_entry(text: &str, line: usize) -> Result<MixtureEntry> {
        let tokens: Vec<&str> = text.split_whitespace().collect();
        let kind = tokens.first().map(|token| token.to_ascii_lowercase());
        match kind.as_deref() {
            Some("material") => {
                if tokens.len() != 4 {
                    return Err(Error::Parse {
                        line,
                        msg: format!(
                            "expected `material <name> <rel_density> <vol_fraction>`, \
                            found `{text}`"
                        ),
                    });
                }
                Ok(MixtureEntry::Material {
                    name: tokens[1].to_string(),
                    rel_density: parse_float(tokens[2], line, "relative density")?,
                    vol_fraction: parse_float(tokens[3], line, "volume fraction")?,
                })
            }
            Some("element") => {
                if tokens.len() != 4 {
                    return Err(Error::Parse {
                        line,
                        msg: format!(
                            "expected `element <symbol> <rel_density> <vol_fraction>`, \
                            found `{text}`"
                        ),
                    });
                }
                validate_element_symbol(tokens[1], line)?;
                Ok(MixtureEntry::Element {
                    symbol: tokens[1].to_string(),
                    rel_density: parse_float(tokens[2], line, "relative density")?,
                    vol_fraction: parse_float(tokens[3], line, "volume fraction")?,
                })
            }
            Some("like") => {
                if tokens.len() != 3 {
                    return Err(Error::Parse {
                        line,
                        msg: format!("expected `like <mixture> <rel_density>`, found `{text}`"),
                    });
                }
                Ok(MixtureEntry::Like {
                    mixture: tokens[1].to_string(),
                    rel_density: parse_float(tokens[2], line, "relative density")?,
                })
            }
            Some("target") => {
                if tokens.len() != 3 {
                    return Err(Error::Parse {
                        line,
                        msg: format!("expected `target <element|isotope> <name>`, found `{text}`"),
                    });
                }
                let target_kind = tokens[1].to_ascii_lowercase();
                if target_kind == "element" {
                    validate_element_symbol(tokens[2], line)?;
                } else if target_kind == "isotope" {
                    NuclideId::from_name(tokens[2]).map_err(|_| Error::Parse {
                        line,
                        msg: format!("unknown target isotope `{}`", tokens[2]),
                    })?;
                } else {
                    return Err(Error::Parse {
                        line,
                        msg: format!(
                            "expected target of type `element` or `isotope`, found `{}`",
                            tokens[1]
                        ),
                    });
                }
                Ok(MixtureEntry::Target {
                    target_kind,
                    name: tokens[2].to_string(),
                })
            }
            _ => Err(Error::Parse {
                line,
                msg: format!(
                    "expected `material`, `element`, `like`, or `target` entry, found `{text}`"
                ),
            }),
        }
    }

    /// Reject a block keyword outside [`KNOWN_BLOCKS`]; `line` is 1-based.
    /// Matching is case-insensitive.
    pub fn check_block(block: &str, line: usize) -> Result<()> {
        if KNOWN_BLOCKS.contains(&block.to_ascii_lowercase().as_str()) {
            Ok(())
        } else {
            Err(Error::UnknownBlock {
                block: block.to_string(),
                line,
            })
        }
    }

    /// Check dangling cross-references between blocks.
    ///
    /// Verifies the `dimension`/`volume(s)` exclusion, `mat_loading` mixture
    /// names (or `void`), `like` mixture references, `volume(s)` and
    /// solve/skip zone names against `mat_loading`, and schedule flux /
    /// pulsing / sub-schedule names against the defined `flux`,
    /// `pulsehistory`, and `schedule` blocks. The first dangling reference
    /// is reported as [`Error::CrossRef`].
    pub fn validate(&self) -> Result<()> {
        if !self.dimensions.is_empty() && self.volumes.is_some() {
            return Err(Error::CrossRef(
                "input defines both `dimension` and `volume(s)` blocks; \
                only one geometry method is allowed"
                    .to_string(),
            ));
        }
        let mixtures: BTreeSet<&str> = self.mixtures.iter().map(|mix| mix.name.as_str()).collect();
        if let Some(loading) = &self.mat_loading {
            for entry in &loading.entries {
                if entry.mixture.eq_ignore_ascii_case("void") {
                    continue;
                }
                if !mixtures.contains(entry.mixture.as_str()) {
                    return Err(Error::CrossRef(format!(
                        "mat_loading zone `{}` names undefined mixture `{}`",
                        entry.zone, entry.mixture
                    )));
                }
            }
        }
        for mix in &self.mixtures {
            for entry in &mix.entries {
                if let MixtureEntry::Like { mixture, .. } = entry {
                    if !mixtures.contains(mixture.as_str()) {
                        return Err(Error::CrossRef(format!(
                            "mixture `{}` reuses undefined mixture `{mixture}`",
                            mix.name
                        )));
                    }
                }
            }
        }
        if let (Some(volumes), Some(loading)) = (&self.volumes, &self.mat_loading) {
            let zones: BTreeSet<&str> = loading
                .entries
                .iter()
                .map(|entry| entry.zone.as_str())
                .collect();
            for entry in &volumes.entries {
                if !zones.contains(entry.zone.as_str()) {
                    return Err(Error::CrossRef(format!(
                        "volumes entry names undefined zone `{}`",
                        entry.zone
                    )));
                }
            }
        }
        if let Some(loading) = &self.mat_loading {
            let zones: BTreeSet<&str> = loading
                .entries
                .iter()
                .map(|entry| entry.zone.as_str())
                .collect();
            for (kind, list) in [
                (
                    "solve_zones",
                    self.solve_zones.as_ref().map(|zone| &zone.zones),
                ),
                (
                    "skip_zones",
                    self.skip_zones.as_ref().map(|zone| &zone.zones),
                ),
            ] {
                if let Some(names) = list {
                    for zone in names {
                        if !zones.contains(zone.as_str()) {
                            return Err(Error::CrossRef(format!(
                                "`{kind}` names undefined zone `{zone}`"
                            )));
                        }
                    }
                }
            }
        }
        let fluxes: BTreeSet<&str> = self.fluxes.iter().map(|flux| flux.name.as_str()).collect();
        let pulses: BTreeSet<&str> = self
            .pulse_histories
            .iter()
            .map(|history| history.name.as_str())
            .collect();
        let schedules: BTreeSet<&str> = self
            .schedules
            .iter()
            .map(|schedule| schedule.name.as_str())
            .collect();
        for schedule in &self.schedules {
            for item in &schedule.items {
                match item.tokens.len() {
                    4 => {
                        if !schedules.contains(item.tokens[0].as_str()) {
                            return Err(Error::CrossRef(format!(
                                "schedule `{}` references undefined sub-schedule `{}`",
                                schedule.name, item.tokens[0]
                            )));
                        }
                        if !pulses.contains(item.tokens[1].as_str()) {
                            return Err(Error::CrossRef(format!(
                                "schedule `{}` references undefined pulsing definition `{}`",
                                schedule.name, item.tokens[1]
                            )));
                        }
                    }
                    6 => {
                        if !fluxes.contains(item.tokens[2].as_str()) {
                            return Err(Error::CrossRef(format!(
                                "schedule `{}` references undefined flux `{}`",
                                schedule.name, item.tokens[2]
                            )));
                        }
                        if !pulses.contains(item.tokens[3].as_str()) {
                            return Err(Error::CrossRef(format!(
                                "schedule `{}` references undefined pulsing definition `{}`",
                                schedule.name, item.tokens[3]
                            )));
                        }
                    }
                    _ => {
                        return Err(Error::CrossRef(format!(
                            "schedule `{}` has a malformed item",
                            schedule.name
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// Resolve the directly nuclide-addressable entries of a mixture to
    /// canonical [`NuclideId`]s.
    ///
    /// `element` entries whose symbol parses as a nuclide (e.g. `U235`,
    /// `mn:56`) and `target isotope` entries contribute their volume
    /// fraction (targets contribute `1.0`); `material`, `like`, and
    /// `target element` entries need library data or hierarchy expansion and
    /// surface as [`Error::CrossRef`], as do unknown mixtures and
    /// unparseable names.
    pub fn mixture_ids(&self, mixture: &str) -> Result<BTreeMap<NuclideId, f64>> {
        let found = self
            .mixtures
            .iter()
            .find(|mix| mix.name == mixture)
            .ok_or_else(|| Error::CrossRef(format!("undefined mixture `{mixture}`")))?;
        let mut ids = BTreeMap::new();
        for entry in &found.entries {
            match entry {
                MixtureEntry::Element {
                    symbol,
                    vol_fraction,
                    ..
                } => {
                    let candidate = match symbol.split_once(':') {
                        Some((base, suffix)) if suffix.chars().all(|c| c.is_ascii_digit()) => {
                            format!("{base}{suffix}")
                        }
                        _ => symbol.clone(),
                    };
                    match NuclideId::from_name(&candidate) {
                        Ok(id) => {
                            ids.insert(id, *vol_fraction);
                        }
                        Err(_) => {
                            return Err(Error::CrossRef(format!(
                                "mixture `{mixture}` entry `{symbol}` is not a resolvable \
                                nuclide (needs element-library data)"
                            )));
                        }
                    }
                }
                MixtureEntry::Target {
                    target_kind, name, ..
                } if target_kind == "isotope" => match NuclideId::from_name(name) {
                    Ok(id) => {
                        ids.insert(id, 1.0);
                    }
                    Err(_) => {
                        return Err(Error::CrossRef(format!(
                            "mixture `{mixture}` names unknown nuclide `{name}`"
                        )));
                    }
                },
                _ => {
                    return Err(Error::CrossRef(format!(
                        "mixture `{mixture}` needs library data or hierarchy expansion \
                        to resolve nuclides"
                    )));
                }
            }
        }
        Ok(ids)
    }

    /// Number of recorded blocks.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// True when no blocks were recorded.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Block kinds in file order (lowercase, `volume` normalized to `volumes`).
    pub fn block_kinds(&self) -> Vec<&str> {
        self.blocks
            .iter()
            .map(|block| block.kind.as_str())
            .collect()
    }

    /// Find mixtures by name.
    pub fn find_mixture(&self, name: &str) -> Option<&Mixture> {
        self.mixtures.iter().find(|mix| mix.name == name)
    }

    /// Find a flux definition by name.
    pub fn find_flux(&self, name: &str) -> Option<&FluxDef> {
        self.fluxes.iter().find(|flux| flux.name == name)
    }

    /// Find a schedule definition by name.
    pub fn find_schedule(&self, name: &str) -> Option<&ScheduleDef> {
        self.schedules.iter().find(|schedule| schedule.name == name)
    }

    /// Find a pulsing history by name.
    pub fn find_pulse_history(&self, name: &str) -> Option<&PulseHistory> {
        self.pulse_histories
            .iter()
            .find(|history| history.name == name)
    }
}

// ---------------------------------------------------------------------------
// Canonical writer
// ---------------------------------------------------------------------------

impl fmt::Display for AlaraDeck {
    /// Re-emit a canonical deck.
    ///
    /// Blocks keep their file order; values are normalized (lowercase
    /// keywords, `volume` spelled `volumes`, times in seconds with `s`
    /// units). The output is not byte-identical to the input, but
    /// `parse(write(parse(x)))` is stable: re-parsing the output yields the
    /// same canonical text.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut dimensions = self.dimensions.iter();
        let mut mixtures = self.mixtures.iter();
        let mut fluxes = self.fluxes.iter();
        let mut schedules = self.schedules.iter();
        let mut pulses = self.pulse_histories.iter();
        let mut outputs = self.outputs.iter();
        for block in &self.blocks {
            match block.kind.as_str() {
                "geometry" => {
                    if let Some(geometry) = &self.geometry {
                        writeln!(f, "geometry {}", geometry.kind)?;
                    }
                }
                "dimension" => {
                    if let Some(dimension) = dimensions.next() {
                        writeln!(f, "dimension {} {}", dimension.axis, dimension.lower)?;
                        for zone in &dimension.zones {
                            writeln!(f, "{} {}", zone.intervals, zone.upper)?;
                        }
                        writeln!(f, "end")?;
                    } else {
                        Self::fmt_raw(f, block, true)?;
                    }
                }
                "volumes" => {
                    if let Some(volumes) = &self.volumes {
                        writeln!(f, "volumes")?;
                        for entry in &volumes.entries {
                            writeln!(f, "{} {}", entry.volume, entry.zone)?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "mat_loading" => {
                    if let Some(loading) = &self.mat_loading {
                        writeln!(f, "mat_loading")?;
                        for entry in &loading.entries {
                            writeln!(f, "{} {}", entry.zone, entry.mixture)?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "mixture" => {
                    if let Some(mix) = mixtures.next() {
                        writeln!(f, "mixture {}", mix.name)?;
                        for entry in &mix.entries {
                            match entry {
                                MixtureEntry::Material {
                                    name,
                                    rel_density,
                                    vol_fraction,
                                } => writeln!(f, "material {name} {rel_density} {vol_fraction}")?,
                                MixtureEntry::Element {
                                    symbol,
                                    rel_density,
                                    vol_fraction,
                                } => writeln!(f, "element {symbol} {rel_density} {vol_fraction}")?,
                                MixtureEntry::Like {
                                    mixture,
                                    rel_density,
                                } => writeln!(f, "like {mixture} {rel_density}")?,
                                MixtureEntry::Target { target_kind, name } => {
                                    writeln!(f, "target {target_kind} {name}")?;
                                }
                            }
                        }
                        writeln!(f, "end")?;
                    } else {
                        Self::fmt_raw(f, block, true)?;
                    }
                }
                "solve_zones" => {
                    if let Some(zones) = &self.solve_zones {
                        writeln!(f, "solve_zones")?;
                        for zone in &zones.zones {
                            writeln!(f, "{zone}")?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "skip_zones" => {
                    if let Some(zones) = &self.skip_zones {
                        writeln!(f, "skip_zones")?;
                        for zone in &zones.zones {
                            writeln!(f, "{zone}")?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "flux" => {
                    if let Some(flux) = fluxes.next() {
                        writeln!(
                            f,
                            "flux {} {} {} {} {}",
                            flux.name, flux.file, flux.scale, flux.skip, flux.format
                        )?;
                    } else {
                        Self::fmt_raw(f, block, false)?;
                    }
                }
                "spatial_norm" => {
                    if let Some(norm) = &self.spatial_norm {
                        writeln!(f, "spatial_norm")?;
                        for value in &norm.values {
                            writeln!(f, "{value}")?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "schedule" => {
                    if let Some(schedule) = schedules.next() {
                        writeln!(f, "schedule {}", schedule.name)?;
                        for item in &schedule.items {
                            writeln!(f, "{}", item.tokens.join(" "))?;
                        }
                        writeln!(f, "end")?;
                    } else {
                        Self::fmt_raw(f, block, true)?;
                    }
                }
                "pulsehistory" => {
                    if let Some(history) = pulses.next() {
                        writeln!(f, "pulsehistory {}", history.name)?;
                        for level in &history.levels {
                            writeln!(f, "{} {} s", level.pulses, level.delay_s)?;
                        }
                        writeln!(f, "end")?;
                    } else {
                        Self::fmt_raw(f, block, true)?;
                    }
                }
                "truncation" => {
                    if let Some(truncation) = &self.truncation {
                        writeln!(f, "truncation {}", truncation.tolerance)?;
                    }
                }
                "impurity" => {
                    if let Some(impurity) = &self.impurity {
                        writeln!(f, "impurity {} {}", impurity.threshold, impurity.tolerance)?;
                    }
                }
                "ignore" => {
                    if let Some(ignore) = &self.ignore {
                        writeln!(f, "ignore {}", ignore.tolerance)?;
                    }
                }
                "ref_flux_type" => {
                    if let Some(reference) = &self.ref_flux_type {
                        writeln!(f, "ref_flux_type {}", reference.kind)?;
                    }
                }
                "cooling" => {
                    if let Some(cooling) = &self.cooling {
                        writeln!(f, "cooling")?;
                        for time in &cooling.times_s {
                            writeln!(f, "{time} s")?;
                        }
                        writeln!(f, "end")?;
                    }
                }
                "output" => {
                    if let Some(output) = outputs.next() {
                        writeln!(f, "output {}", output.resolution)?;
                        for entry in &output.entries {
                            writeln!(f, "{entry}")?;
                        }
                        writeln!(f, "end")?;
                    } else {
                        Self::fmt_raw(f, block, true)?;
                    }
                }
                "material_lib" => {
                    if let Some(library) = &self.material_lib {
                        writeln!(f, "material_lib {}", library.path)?;
                    }
                }
                "element_lib" => {
                    if let Some(library) = &self.element_lib {
                        writeln!(f, "element_lib {}", library.path)?;
                    }
                }
                "data_library" => {
                    if let Some(library) = &self.data_library {
                        writeln!(
                            f,
                            "data_library {} {}",
                            library.kind,
                            library.files.join(" ")
                        )?;
                    }
                }
                "dump_file" => {
                    if let Some(dump) = &self.dump_file {
                        writeln!(f, "dump_file {}", dump.path)?;
                    }
                }
                "major_radius" => {
                    if let Some(radius) = &self.major_radius {
                        writeln!(f, "major_radius {}", radius.value)?;
                    }
                }
                "minor_radius" => {
                    if let Some(radius) = &self.minor_radius {
                        writeln!(f, "minor_radius {}", radius.value)?;
                    }
                }
                "convert_lib" => {
                    if let Some(convert) = &self.convert_lib {
                        writeln!(
                            f,
                            "convert_lib {} {} {}",
                            convert.from,
                            convert.to,
                            convert.files.join(" ")
                        )?;
                    }
                }
                _ => {
                    Self::fmt_raw(f, block, false)?;
                }
            }
        }
        Ok(())
    }
}

impl AlaraDeck {
    /// Fallback writer for blocks without a typed view (only reachable for
    /// programmatically built decks whose `blocks` outnumber typed entries).
    fn fmt_raw(f: &mut fmt::Formatter<'_>, block: &RawBlock, end: bool) -> fmt::Result {
        if block.body.is_empty() {
            writeln!(f, "{}", block.kind)?;
        } else {
            writeln!(f, "{} {}", block.kind, block.body.join(" "))?;
        }
        if end {
            writeln!(f, "end")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/alara/decks")
            .join(name)
    }

    fn day() -> f64 {
        86_400.0
    }

    #[test]
    fn parse_sample2_fixture() {
        let path = fixture("sample2");
        let text = std::fs::read_to_string(&path).unwrap();
        let from_text = AlaraDeck::parse(&text).unwrap();
        let from_file = AlaraDeck::from_file(&path).unwrap();
        assert_eq!(from_text.block_kinds(), from_file.block_kinds());

        assert_eq!(
            from_file.block_kinds(),
            [
                "geometry",
                "dimension",
                "mat_loading",
                "material_lib",
                "element_lib",
                "mixture",
                "mixture",
                "flux",
                "schedule",
                "pulsehistory",
                "dump_file",
                "data_library",
                "cooling",
                "output",
                "truncation",
            ]
        );

        let geometry = from_file.geometry.as_ref().unwrap();
        assert_eq!(geometry.kind, "rectangular");

        let dimension = from_file.dimensions.first().unwrap();
        assert_eq!(dimension.axis, "x");
        assert_eq!(dimension.lower, 0.0);
        assert_eq!(dimension.zones.len(), 2);
        assert_eq!(dimension.zones[0].intervals, 1);
        assert_eq!(dimension.zones[0].upper, 5.0);
        assert_eq!(dimension.zones[1].intervals, 1);
        assert_eq!(dimension.zones[1].upper, 6.0);

        let loading = from_file.mat_loading.as_ref().unwrap();
        assert_eq!(loading.entries.len(), 2);
        assert_eq!(loading.entries[0].zone, "inner_zone");
        assert_eq!(loading.entries[0].mixture, "inner_mix");

        assert_eq!(
            from_file
                .mixtures
                .iter()
                .map(|mix| mix.name.as_str())
                .collect::<Vec<_>>(),
            ["inner_mix", "outer_mix"]
        );
        let inner = from_file.find_mixture("inner_mix").unwrap();
        assert_eq!(inner.entries.len(), 3);
        assert!(matches!(
            &inner.entries[0],
            MixtureEntry::Material { name, rel_density, vol_fraction }
                if name == "Li" && *rel_density == 1.0 && *vol_fraction == 0.75
        ));

        let flux = from_file.find_flux("flux_1").unwrap();
        assert_eq!(flux.file, "data/fluxin1");
        assert_eq!(flux.scale, 1.0);
        assert_eq!(flux.skip, 1);
        assert_eq!(flux.format, "default");

        let schedule = from_file.find_schedule("1_year").unwrap();
        assert_eq!(schedule.items.len(), 1);
        assert_eq!(schedule.items[0].tokens.len(), 6);
        assert_eq!(schedule.items[0].tokens[2], "flux_1");

        let history = from_file.find_pulse_history("steady_state").unwrap();
        assert_eq!(history.levels.len(), 1);
        assert_eq!(history.levels[0].pulses, 1);
        assert_eq!(history.levels[0].delay_s, 0.0);

        let cooling = from_file.cooling.as_ref().unwrap();
        assert_eq!(cooling.times_s.len(), 4);
        assert_eq!(cooling.times_s[0], day());
        assert_eq!(cooling.times_s[1], 100.0 * day());
        assert_eq!(cooling.times_s[2], 365.25 * day());

        assert_eq!(from_file.outputs.len(), 1);
        assert_eq!(from_file.outputs[0].resolution, "interval");
        assert_eq!(
            from_file.outputs[0].entries,
            ["number_density", "specific_activity"]
        );

        assert_eq!(from_file.truncation.as_ref().unwrap().tolerance, 1e-7);
        assert_eq!(
            from_file.material_lib.as_ref().unwrap().path,
            "data/sampleMatlib"
        );
        assert_eq!(
            from_file.element_lib.as_ref().unwrap().path,
            "data/myElelib"
        );
        let data_library = from_file.data_library.as_ref().unwrap();
        assert_eq!(data_library.kind, "alaralib");
        assert_eq!(data_library.files, ["data/truncated_fendl2bin"]);
        assert_eq!(
            from_file.dump_file.as_ref().unwrap().path,
            "dump_files/sample2.dump"
        );

        from_file.validate().unwrap();
        assert_eq!(from_file.len(), 15);
        assert!(!from_file.is_empty());
    }

    #[test]
    fn parse_sample3_fixture() {
        let path = fixture("sample3");
        let text = std::fs::read_to_string(&path).unwrap();
        let from_text = AlaraDeck::parse(&text).unwrap();
        let from_file = AlaraDeck::from_file(&path).unwrap();
        assert_eq!(from_text.block_kinds(), from_file.block_kinds());

        assert_eq!(from_file.geometry.as_ref().unwrap().kind, "cylindrical");
        let dimension = from_file.dimensions.first().unwrap();
        assert_eq!(dimension.axis, "r");
        assert_eq!(dimension.lower, 0.0);
        assert_eq!(dimension.zones.len(), 2);
        assert_eq!(dimension.zones[1].intervals, 2);
        assert_eq!(dimension.zones[1].upper, 15.0);

        assert_eq!(
            from_file
                .mixtures
                .iter()
                .map(|mix| mix.name.as_str())
                .collect::<Vec<_>>(),
            ["inner_mix", "outer_mix"]
        );
        let inner = from_file.find_mixture("inner_mix").unwrap();
        assert_eq!(inner.entries.len(), 3);
        assert!(matches!(
            &inner.entries[2],
            MixtureEntry::Element { symbol, rel_density, vol_fraction }
                if symbol == "ni" && *rel_density == 1.0 && *vol_fraction == 0.70
        ));
        let outer = from_file.find_mixture("outer_mix").unwrap();
        assert!(matches!(
            &outer.entries[0],
            MixtureEntry::Element { symbol, rel_density, vol_fraction }
                if symbol == "fe" && *rel_density == 2.0 && *vol_fraction == 1.0
        ));

        let flux = from_file.find_flux("flux_1").unwrap();
        assert_eq!(flux.file, "data/fluxin2");
        assert_eq!(flux.scale, 1e6);
        assert_eq!(flux.skip, 0);
        assert_eq!(flux.format, "default");

        let schedule = from_file.find_schedule("total").unwrap();
        assert_eq!(schedule.items.len(), 1);

        let history = from_file.find_pulse_history("pulsed_spec").unwrap();
        assert_eq!(history.levels.len(), 1);
        assert_eq!(history.levels[0].pulses, 10);
        assert_eq!(history.levels[0].delay_s, 5.0);

        let cooling = from_file.cooling.as_ref().unwrap();
        assert_eq!(cooling.times_s.len(), 4);
        assert_eq!(cooling.times_s[0], 1.0);
        assert_eq!(cooling.times_s[1], day());
        assert_eq!(cooling.times_s[2], 60.0);
        assert_eq!(cooling.times_s[3], 365.25 * day());

        assert_eq!(from_file.outputs.len(), 1);
        assert_eq!(from_file.outputs[0].resolution, "zone");
        assert_eq!(from_file.outputs[0].entries, ["constituent", "total_heat"]);
        assert_eq!(from_file.truncation.as_ref().unwrap().tolerance, 1e-8);

        from_file.validate().unwrap();
    }

    #[test]
    fn writer_round_trip_is_stable() {
        for name in ["sample2", "sample3"] {
            let text = std::fs::read_to_string(fixture(name)).unwrap();
            let first = AlaraDeck::parse(&text).unwrap().to_string();
            let second = AlaraDeck::parse(&first).unwrap().to_string();
            assert_eq!(first, second, "round trip drifted for {name}");
            // Canonical output re-parses to the same typed deck modulo lines.
            let reparsed = AlaraDeck::parse(&second).unwrap();
            reparsed.validate().unwrap();
            assert!(!reparsed.is_empty());
        }
    }

    #[test]
    fn comments_and_blank_lines_are_ignored() {
        let text = "## full-line comment\n   # indented comment\n\ngeometry rectangular # trailing\n\ntruncation 1e-7\n";
        let deck = AlaraDeck::parse(text).unwrap();
        assert_eq!(deck.block_kinds(), ["geometry", "truncation"]);
        assert_eq!(deck.geometry.as_ref().unwrap().kind, "rectangular");
    }

    #[test]
    fn trailing_comment_after_single_word_value() {
        let text =
            "geometry rectangular\nmixture mix_0\nelement mn:56 1.0 1 # to give atoms/cm3\nend\n";
        let deck = AlaraDeck::parse(text).unwrap();
        let mix = deck.find_mixture("mix_0").unwrap();
        assert!(matches!(
            &mix.entries[0],
            MixtureEntry::Element { symbol, .. } if symbol == "mn:56"
        ));
    }

    #[test]
    fn include_between_blocks_is_inlined() {
        let dir = std::env::temp_dir().join(format!("alara_io_include_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("inc"), "truncation 1e-7\n").unwrap();
        std::fs::write(
            dir.join("main"),
            "geometry rectangular\n#include inc\nmaterial_lib data/lib\n",
        )
        .unwrap();
        let deck = AlaraDeck::from_file(dir.join("main")).unwrap();
        assert_eq!(
            deck.block_kinds(),
            ["geometry", "truncation", "material_lib"]
        );
        assert_eq!(deck.truncation.as_ref().unwrap().tolerance, 1e-7);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn include_resolves_relative_to_including_file() {
        let dir = std::env::temp_dir().join(format!("alara_io_rel_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        std::fs::write(dir.join("sub").join("inc"), "truncation 2e-6\n").unwrap();
        std::fs::write(
            dir.join("main"),
            "geometry point\n#include <sub/inc>\nmaterial_lib data/lib\n",
        )
        .unwrap();
        let deck = AlaraDeck::from_file(dir.join("main")).unwrap();
        assert_eq!(deck.truncation.as_ref().unwrap().tolerance, 2e-6);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn include_inside_block_is_an_error() {
        let dir = std::env::temp_dir().join(format!("alara_io_nest_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("inc"), "truncation 1e-7\n").unwrap();
        std::fs::write(
            dir.join("main"),
            "geometry rectangular\nmixture fuel\n#include inc\nend\n",
        )
        .unwrap();
        let error = AlaraDeck::from_file(dir.join("main")).unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "got {error}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn include_in_memory_text_is_an_error() {
        let error = AlaraDeck::parse("geometry rectangular\n#include other\n").unwrap_err();
        assert!(matches!(error, Error::Parse { .. }));
    }

    #[test]
    fn missing_end_is_an_error() {
        let error =
            AlaraDeck::parse("geometry rectangular\nmixture fuel\nmaterial Li 1 1\n").unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "got {error}");
    }

    #[test]
    fn unknown_block_is_an_error_with_line_number() {
        let error = AlaraDeck::parse("geometry rectangular\nbogus 1\n").unwrap_err();
        match error {
            Error::UnknownBlock { block, line } => {
                assert_eq!(block, "bogus");
                assert_eq!(line, 2);
            }
            other => panic!("expected UnknownBlock, got {other}"),
        }
    }

    #[test]
    fn time_units_convert_to_seconds() {
        assert_eq!(parse_time_to_seconds(1.0, 's').unwrap(), 1.0);
        assert_eq!(parse_time_to_seconds(1.0, 'S').unwrap(), 1.0);
        assert_eq!(parse_time_to_seconds(2.0, 'm').unwrap(), 120.0);
        assert_eq!(parse_time_to_seconds(2.0, 'H').unwrap(), 7_200.0);
        assert_eq!(parse_time_to_seconds(1.0, 'd').unwrap(), 86_400.0);
        assert_eq!(parse_time_to_seconds(1.0, 'w').unwrap(), 7.0 * 86_400.0);
        // One year is 365.25 days by scaffold convention.
        assert_eq!(parse_time_to_seconds(1.0, 'y').unwrap(), 365.25 * 86_400.0);
        assert_eq!(
            parse_time_to_seconds(1.0, 'c').unwrap(),
            100.0 * 365.25 * 86_400.0
        );
        assert!(matches!(
            parse_time_to_seconds(1.0, 'x'),
            Err(Error::BadUnits(_))
        ));
        assert!(matches!(
            parse_time_to_seconds(-1.0, 's'),
            Err(Error::BadUnits(_))
        ));
    }

    #[test]
    fn dimension_and_volumes_conflict_is_an_error() {
        let text = "geometry rectangular\ndimension x 0.0\n1 5.0\nend\nvolumes\n0.5 z1\nend\nmat_loading\nz1 void\nend\n";
        assert!(matches!(AlaraDeck::parse(text), Err(Error::Parse { .. })));
    }

    #[test]
    fn volume_singular_is_accepted() {
        let text = "geometry rectangular\nvolume\n10 zone_1\nend\nmat_loading\nzone_1 void\nend\n";
        let deck = AlaraDeck::parse(text).unwrap();
        assert_eq!(deck.block_kinds(), ["geometry", "volumes", "mat_loading"]);
        assert_eq!(deck.volumes.as_ref().unwrap().entries.len(), 1);
        deck.validate().unwrap();
    }

    #[test]
    fn keywords_are_case_insensitive() {
        let text = "Geometry rectangular\nMixture fuel\nMaterial WATER 1 1\nEND\n";
        let deck = AlaraDeck::parse(text).unwrap();
        assert_eq!(deck.block_kinds(), ["geometry", "mixture"]);
        assert!(deck.find_mixture("fuel").is_some());
    }

    #[test]
    fn parse_rejects_empty_deck() {
        assert!(matches!(
            AlaraDeck::parse("   \n"),
            Err(Error::Parse { line: 1, .. })
        ));
    }

    #[test]
    fn check_block_uses_allowlist() {
        assert!(AlaraDeck::check_block("geometry", 1).is_ok());
        assert!(AlaraDeck::check_block("Geometry", 1).is_ok());
        assert!(AlaraDeck::check_block("volumes", 1).is_ok());
        assert!(matches!(
            AlaraDeck::check_block("bogus", 3),
            Err(Error::UnknownBlock { line: 3, .. })
        ));
    }

    #[test]
    fn validate_catches_dangling_references() {
        let base = "geometry rectangular\nmaterial_lib a\nelement_lib b\ndata_library alaralib c\n";
        // Unknown flux in a schedule item.
        let deck = AlaraDeck::parse(&format!(
            "{base}flux f1 data/x 1 0 default\nschedule s\n1 d f_missing ph 0 s\nend\npulsehistory ph\n1 0 s\nend\n"
        ))
        .unwrap();
        assert!(matches!(deck.validate(), Err(Error::CrossRef(_))));

        // Unknown mixture in mat_loading.
        let deck = AlaraDeck::parse(&format!("{base}mat_loading\nz1 no_such_mix\nend\n")).unwrap();
        assert!(matches!(deck.validate(), Err(Error::CrossRef(_))));

        // Unknown pulsing definition.
        let deck = AlaraDeck::parse(&format!(
            "{base}flux f1 data/x 1 0 default\nschedule s\n1 d f1 no_pulse 0 s\nend\n"
        ))
        .unwrap();
        assert!(matches!(deck.validate(), Err(Error::CrossRef(_))));

        // Unknown sub-schedule.
        let deck = AlaraDeck::parse(&format!(
            "{base}schedule top\nsub_missing ph 1 d\nend\npulsehistory ph\n1 0 s\nend\n"
        ))
        .unwrap();
        assert!(matches!(deck.validate(), Err(Error::CrossRef(_))));

        // `like` pointing at an undefined mixture.
        let deck = AlaraDeck::parse(&format!(
            "{base}mixture a\nlike ghost 1.0\nend\nmat_loading\nz a\nend\n"
        ))
        .unwrap();
        assert!(matches!(deck.validate(), Err(Error::CrossRef(_))));

        // `void` zones and consistent references pass.
        let deck = AlaraDeck::parse(&format!(
            "{base}mixture a\nmaterial WATER 1 1\nend\nmat_loading\nz1 a\nz2 void\nend\n\
            flux f1 data/x 1 0 default\nschedule s\n1 d f1 ph 0 s\nend\npulsehistory ph\n1 0 s\nend\n"
        ))
        .unwrap();
        deck.validate().unwrap();
    }

    #[test]
    fn mixture_ids_resolves_nuclides() {
        let deck = AlaraDeck::parse(
            "geometry rectangular\nmixture fuel\nelement U235 1 0.04\nelement U238 1 0.96\nend\n",
        )
        .unwrap();
        let ids = deck.mixture_ids("fuel").unwrap();
        assert_eq!(ids[&NuclideId::from_name("U235").unwrap()], 0.04);
        assert_eq!(ids[&NuclideId::from_name("U238").unwrap()], 0.96);
        assert!(matches!(
            deck.mixture_ids("missing"),
            Err(Error::CrossRef(_))
        ));
        // Library-backed entries cannot resolve to nuclide ids.
        let deck = AlaraDeck::parse("geometry rectangular\nmixture bad\nmaterial WATER 1 1\nend\n")
            .unwrap();
        assert!(matches!(deck.mixture_ids("bad"), Err(Error::CrossRef(_))));
    }
}
