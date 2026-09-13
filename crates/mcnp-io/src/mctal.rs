//! MCNP MCTAL output parsing — header, tally bodies, mesh tallies, and
//! `kcode` data.
//!
//! Reads the header, per-tally bodies for the standard layout plus
//! rectangular/cylindrical/spherical mesh tallies (`detector_type <= -1`),
//! and `kcode` criticality data.
//!
//! Body card order (`f`/`d`/`u`/`s`/`m`/`c`/`e`/`t`, then `vals`, then an
//! optional `tfc` block) follows the MCTAL tally-block layout: one count
//! line per bin card plus that many values, then a `vals` block of
//! `(value, rel_error)` pairs stored in file order, then the tally
//! fluctuation chart (`tfc` jtf line plus 3–4-float data rows) when present.
//! Upstream PyNE's `Mctal` reads header and `kcode` data only (tally numbers
//! are collected, bodies skipped), so headers and cycles stay
//! byte-compatible with that negative oracle; bodies have no upstream oracle
//! and are validated with hand-built synthetic fixtures instead (closed-form
//! val/err pairing).
//!
//! Supported subset: bin cards spelled `<letter>[t|c][<tally>]` (bare `d` or
//! numbered `d4`, total `ut`/`ut4`, cumulative `uc`/`uc4`, same for
//! `f`/`s`/`m`/`c`/`e`/`t`) with a count plus that many values; a third
//! flag token on `c`/`e`/`t` cards (stored, never interpreted); a `vals`
//! block holding exactly `2 * prod(counts, 0 -> 1)` floats paired as
//! `(value, rel_error)` in file order (mesh tallies multiply by the
//! `ni*nj*nk` mesh cells; pairs stored verbatim, no bin-to-pair mapping
//! assumed); an optional `tfc` block after `vals` for standard tallies
//! (mesh tallies carry no `tfc`); mesh `f` lines carrying the 4-int mesh
//! info (`tally unknown ni nj nk`, bare-`f` or `f<tally>` spellings) plus
//! `(ni+1)+(nj+1)+(nk+1)` cora/b/c bounds. Mesh cell ordering matches the
//! writer loop (`i` fastest, `k` slowest among the mesh axes); the
//! rectangular (`-1`) cora/b/c map onto x/y/z with the same expanded-bounds
//! convention as `meshtal.rs`.
//! Everything else is a loud named-open error, not a silent skip:
//! radiograph/point-detector specials (`detector_type >= 3`, tally names
//! ending in 5 with elided objects), negative-`particle_type` particle lists
//! beyond one plain line, and perturbation (`npert`) bodies.
//!
//! No public MCTAL fixture corpus exists, so validation uses hand-built
//! synthetic files in `fixtures/mcnp/mctal/`: the two legacy kcode-only
//! files exercising the 5-value and 19-value cycle record variants, a
//! body-bearing file with closed-form val/err pairing, a mesh-tally file
//! with closed-form pairing, and a tfc/total-variant/flag file.

use std::fmt;
use std::path::Path;

/// Errors raised while parsing MCTAL files.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    BadStructure(String),
    BadNumber { context: &'static str, text: String },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadStructure(m) => write!(f, "malformed MCTAL: {m}"),
            Error::BadNumber { context, text } => {
                write!(f, "cannot parse {context} from `{text}`")
            }
        }
    }
}

impl std::error::Error for Error {}

fn num(context: &'static str, tok: &str) -> Result<f64, Error> {
    tok.parse::<f64>().map_err(|_| Error::BadNumber {
        context,
        text: tok.to_string(),
    })
}

/// Averaged value with its standard deviation.
pub type AvgStdev = (f64, f64);

/// Per-cycle kcode statistics (19-value variant only).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CycleAverages {
    pub avg_k_col: AvgStdev,
    pub avg_k_abs: AvgStdev,
    pub avg_k_path: AvgStdev,
    pub avg_k_combined: AvgStdev,
    pub avg_k_combined_active: AvgStdev,
    pub prompt_life_combined: AvgStdev,
    pub cycle_histories: f64,
    pub fom: f64,
}

/// One bin card: declared count plus its values verbatim.
/// For `f` the values are cell/surface IDs; for the other cards they are bin
/// boundaries. A count of 0 means a single implicit bin with no values.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BinCard {
    /// Declared bin count from the card line.
    pub count: usize,
    /// Values that followed the card line (length must equal `count`,
    /// except `count == 0` which carries no values).
    pub values: Vec<f64>,
    /// Total (`t`) / cumulative (`c`) variant from the card spelling
    /// (`ut`/`uc`/…, `None` for the plain card). The declared `count`
    /// already includes the total/cumulative bin; the flag is stored
    /// verbatim, never interpreted.
    pub variant: Option<char>,
    /// Third-token flag on `c`/`e`/`t` cards (`cosFlag`/`ergFlag`/`timFlag`
    /// in reader parlance); `None` when absent. Stored verbatim, never
    /// interpreted.
    pub flag: Option<i32>,
}

impl BinCard {
    /// Effective bin count (`0` declares a single implicit bin).
    pub fn bins(&self) -> usize {
        if self.count == 0 {
            1
        } else {
            self.count
        }
    }
}

/// One tally-fluctuation-chart data row: history count, tally value,
/// relative error, and optional figure of merit.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TfcEntry {
    /// History (`nps`) count for the row.
    pub nps: i64,
    /// Tally value at this checkpoint.
    pub value: f64,
    /// Relative error at this checkpoint.
    pub rel_err: f64,
    /// Figure of merit, when the row carries a fourth field.
    pub fom: Option<f64>,
}

/// Tally fluctuation chart (`tfc`) block after a standard tally's `vals`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TfcBlock {
    /// Nine-joint `tfc` bin indices (`jtf` line after the `tfc` keyword).
    pub jtf: Vec<i64>,
    /// One entry per data row until the next `tally`/`kcode` block.
    pub rows: Vec<TfcEntry>,
}

/// One parsed standard-tally body (see module docs for the subset).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TallyBody {
    /// Tally number (must appear in the header `tally_nums` list).
    pub number: u32,
    /// Particle-type token from the `tally` line.
    pub particle_type: i32,
    /// Detector-type token from the `tally` line, when present.
    pub detector_type: Option<i32>,
    /// Particle list for negative `particle_type` (one line of integers;
    /// empty otherwise).
    pub particle_list: Vec<i32>,
    /// `FC` comment lines preceding the `f` card (may be empty).
    pub comment: Vec<String>,
    /// Cell/surface/object IDs.
    pub f: BinCard,
    /// Total-vs-direct bins.
    pub d: BinCard,
    /// User bins.
    pub u: BinCard,
    /// Segment bins.
    pub s: BinCard,
    /// Multiplier bins.
    pub m: BinCard,
    /// Cosine bins.
    pub c: BinCard,
    /// Energy bins.
    pub e: BinCard,
    /// Time bins.
    pub t: BinCard,
    /// `(value, rel_error)` pairs in file order (stored verbatim; no
    /// bin-to-pair mapping is assumed).
    pub vals: Vec<(f64, f64)>,
    /// Tally fluctuation chart block after `vals` (`None` when the tally
    /// carries no `tfc` lines).
    pub tfc: Option<TfcBlock>,
}

/// One parsed mesh-tally body (`detector_type <= -1`).
///
/// The `f` line carries the 4-int mesh info instead of cell IDs; cora/b/c
/// hold `(ni+1)+(nj+1)+(nk+1)` bounds in file order. The remaining cards,
/// `vals` pairing, and total/cumulative/flag spellings match [`TallyBody`].
/// Mesh tallies carry no `tfc` block. Cell ordering is the writer loop
/// order: `i` (cora) fastest, then `j`, then `k` slowest among the mesh
/// axes, inside the `f/d/u/s/m/c/e/t` outer bins.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MeshTallyBody {
    /// Tally number (must appear in the header `tally_nums` list).
    pub number: u32,
    /// Particle-type token from the `tally` line.
    pub particle_type: i32,
    /// Mesh-type token from the `tally` line (`-1` rectangular, `-2`
    /// cylindrical, `-3` spherical, `-4..-6` r/c/smesh; stored verbatim).
    pub detector_type: i32,
    /// Particle list for negative `particle_type` (one line of integers;
    /// empty otherwise).
    pub particle_list: Vec<i32>,
    /// `FC` comment lines preceding the `f` card (may be empty).
    pub comment: Vec<String>,
    /// First mesh-info int (unknown/reserved; stored verbatim).
    pub mesh_unknown: i64,
    /// Mesh bin counts along cora/b/c.
    pub ni: usize,
    /// Mesh bin counts along cora/b/c.
    pub nj: usize,
    /// Mesh bin counts along cora/b/c.
    pub nk: usize,
    /// Cora bounds (`ni + 1` values).
    pub cora: Vec<f64>,
    /// Corb bounds (`nj + 1` values).
    pub corb: Vec<f64>,
    /// Corc bounds (`nk + 1` values).
    pub corc: Vec<f64>,
    /// Total-vs-direct bins.
    pub d: BinCard,
    /// User bins.
    pub u: BinCard,
    /// Segment bins.
    pub s: BinCard,
    /// Multiplier bins.
    pub m: BinCard,
    /// Cosine bins.
    pub c: BinCard,
    /// Energy bins.
    pub e: BinCard,
    /// Time bins.
    pub t: BinCard,
    /// `(value, rel_error)` pairs in file order (stored verbatim).
    pub vals: Vec<(f64, f64)>,
}

impl MeshTallyBody {
    /// `[ni, nj, nk]` mesh cell counts.
    pub fn dims(&self) -> [usize; 3] {
        [self.ni, self.nj, self.nk]
    }

    /// Total mesh cells (`ni*nj*nk`, saturating).
    pub fn num_cells(&self) -> usize {
        self.ni.saturating_mul(self.nj).saturating_mul(self.nk)
    }

    /// Flat mesh index for logical cell `(i, j, k)`: `i` fastest
    /// (`(k*nj + j)*ni + i`), matching the writer loop order.
    pub fn mesh_index(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.nj + j) * self.ni + i
    }

    /// Expected pair count: outer-bin product times mesh cells, saturating
    /// at `usize::MAX` instead of wrapping.
    pub fn expected_pairs(&self) -> usize {
        let mut acc = 1usize;
        for bins in [
            self.d.bins(),
            self.u.bins(),
            self.s.bins(),
            self.m.bins(),
            self.c.bins(),
            self.e.bins(),
            self.t.bins(),
            self.num_cells(),
        ] {
            acc = acc.saturating_mul(bins);
            if acc == usize::MAX {
                break;
            }
        }
        acc
    }

    /// Sum of all tally values (errors excluded).
    pub fn total_val(&self) -> f64 {
        self.vals.iter().map(|(v, _)| v).sum()
    }
}

impl TallyBody {
    /// Expected pair count: product of effective bin counts, saturating at
    /// `usize::MAX` instead of wrapping (absurd products surface as an
    /// overflow error at the `vals` gate, never as a small misparse).
    pub fn expected_pairs(&self) -> usize {
        let mut acc = 1usize;
        for bins in [
            self.f.bins(),
            self.d.bins(),
            self.u.bins(),
            self.s.bins(),
            self.m.bins(),
            self.c.bins(),
            self.e.bins(),
            self.t.bins(),
        ] {
            acc = acc.saturating_mul(bins);
            if acc == usize::MAX {
                break;
            }
        }
        acc
    }

    /// Sum of all tally values (errors excluded).
    pub fn total_val(&self) -> f64 {
        self.vals.iter().map(|(v, _)| v).sum()
    }
}

/// Parsed MCTAL file: header, tally bodies, and kcode data.
#[derive(Debug, Clone, PartialEq)]
pub struct Mctal {
    pub code_name: String,
    pub code_version: String,
    pub code_date: String,
    pub code_time: String,
    /// Dump counter token (kept as string upstream too).
    pub n_dump: String,
    pub n_histories: u64,
    pub n_prn: u32,
    /// Input-deck comment card.
    pub comment: String,
    /// Tally count token from the `tally` line (string upstream as well).
    pub n_tallies: String,
    /// Optional third token of the `tally` line (perturbation count when
    /// present; stored verbatim, PyNE parity is to ignore it — perturbation
    /// bodies themselves are named-open).
    pub npert: Option<String>,
    /// Declared tally numbers line (may be empty).
    pub tally_nums: Vec<u32>,
    /// Parsed standard-tally bodies in file order. Legacy kcode-only files
    /// (including the two synthetic kcode fixtures) carry zero bodies even
    /// when they declare tally numbers; body/count completeness is
    /// caller-side.
    pub tallies: Vec<TallyBody>,
    /// Parsed mesh-tally bodies (`detector_type <= -1`) in file order.
    pub mesh_tallies: Vec<MeshTallyBody>,
    pub n_cycles: usize,
    pub n_inactive: usize,
    /// 0/5 = one 5-float line per cycle; 19 = four lines per cycle.
    pub vars_per_cycle: usize,
    /// keff (collision) per cycle.
    pub k_col: Vec<f64>,
    /// keff (absorption) per cycle.
    pub k_abs: Vec<f64>,
    /// keff (track length) per cycle.
    pub k_path: Vec<f64>,
    pub prompt_life_col: Vec<f64>,
    pub prompt_life_path: Vec<f64>,
    /// Running averages block, present when vars_per_cycle >= 19.
    pub averages: Vec<CycleAverages>,
}

impl Mctal {
    /// Read and parse an MCTAL file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Mctal::parse(&text)
    }

    /// Parse MCTAL text in memory.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let raw: Vec<&str> = text.lines().collect();

        let head: Vec<&str> = raw
            .first()
            .ok_or_else(|| Error::BadStructure("empty file".into()))?
            .split_whitespace()
            .collect();
        let mut pos = 1usize;
        if head.len() < 7 {
            return Err(Error::BadStructure(format!(
                "header needs 7 fields, found {}",
                head.len()
            )));
        }
        let code_name = head[0].to_string();
        let code_version = head[1].to_string();
        let code_date = head[2].to_string();
        let code_time = head[3].to_string();
        let n_dump = head[4].to_string();
        let n_histories = num("n_histories", head[5])? as u64;
        let n_prn = num("n_prn", head[6])? as u32;

        let comment = raw.get(pos).unwrap_or(&"").trim().to_string();
        pos += 1;

        let tally_line: Vec<&str> = raw
            .get(pos)
            .ok_or_else(|| Error::BadStructure("missing tally line".into()))?
            .split_whitespace()
            .collect();
        pos += 1;
        if tally_line.first() != Some(&"tally") && !tally_line.is_empty() {
            return Err(Error::BadStructure(format!(
                "expected `tally` line, got `{}`",
                tally_line.join(" ")
            )));
        }
        let n_tallies = tally_line.get(1).copied().unwrap_or("0").to_string();
        let npert = tally_line.get(2).map(|s| s.to_string());

        let tally_nums: Vec<u32> = match raw.get(pos) {
            Some(l) => {
                pos += 1;
                if l.trim().is_empty() {
                    Vec::new()
                } else {
                    l.split_whitespace()
                        .map(|t| num("tally number", t).map(|v| v as u32))
                        .collect::<Result<_, _>>()?
                }
            }
            None => Vec::new(),
        };

        // Per-tally bodies (verifiable subset; see module docs). Bodies are
        // parsed greedily: while the next non-blank line opens a `tally`
        // block, parse it. Legacy kcode-only files declare tally numbers
        // but carry no bodies — zero bodies then is accepted for backward
        // compatibility (PyNE parity path: jump straight to `kcode`).
        let mut tallies = Vec::new();
        let mut mesh_tallies = Vec::new();
        loop {
            skip_blank(&raw, &mut pos);
            let is_tally = raw
                .get(pos)
                .map(|l| l.split_whitespace().next() == Some("tally"))
                .unwrap_or(false);
            if !is_tally {
                break;
            }
            match parse_tally_block(&raw, &mut pos)? {
                TallyBlock::Standard(body) => {
                    if !tally_nums.contains(&body.number) {
                        return Err(Error::BadStructure(format!(
                            "tally body {} is not in the declared tally list {:?}",
                            body.number, tally_nums
                        )));
                    }
                    if tallies.iter().any(|b: &TallyBody| b.number == body.number)
                        || mesh_tallies
                            .iter()
                            .any(|b: &MeshTallyBody| b.number == body.number)
                    {
                        return Err(Error::BadStructure(format!(
                            "duplicate tally body {}",
                            body.number
                        )));
                    }
                    tallies.push(body);
                }
                TallyBlock::Mesh(body) => {
                    if !tally_nums.contains(&body.number) {
                        return Err(Error::BadStructure(format!(
                            "tally body {} is not in the declared tally list {:?}",
                            body.number, tally_nums
                        )));
                    }
                    if tallies.iter().any(|b: &TallyBody| b.number == body.number)
                        || mesh_tallies
                            .iter()
                            .any(|b: &MeshTallyBody| b.number == body.number)
                    {
                        return Err(Error::BadStructure(format!(
                            "duplicate tally body {}",
                            body.number
                        )));
                    }
                    mesh_tallies.push(body);
                }
            }
        }

        // Like upstream, jump to the kcode line (junk-tolerant: blank and
        // non-kcode lines are skipped).
        let mut kcode: Vec<&str> = Vec::new();
        while pos < raw.len() {
            let w: Vec<&str> = raw[pos].split_whitespace().collect();
            pos += 1;
            if w.first() == Some(&"kcode") {
                kcode = w;
                break;
            }
        }
        if kcode.len() < 4 {
            return Err(Error::BadStructure("missing or short kcode line".into()));
        }
        let n_cycles = num("n_cycles", kcode[1])? as usize;
        let n_inactive = num("n_inactive", kcode[2])? as usize;
        let vars_per_cycle = num("vars_per_cycle", kcode[3])? as usize;

        let per_cycle_lines = match vars_per_cycle {
            0 | 5 => 1,
            19 => 4,
            other => {
                return Err(Error::BadStructure(format!(
                    "unsupported vars_per_cycle {other}"
                )))
            }
        };
        let per_cycle_values = if vars_per_cycle == 0 {
            5
        } else {
            vars_per_cycle
        };

        // Cycle vectors grow with the rows actually read: `n_cycles` comes
        // from the file, so no upfront reservation from the declared count.
        let mut k_col = Vec::new();
        let mut k_abs = Vec::new();
        let mut k_path = Vec::new();
        let mut prompt_life_col = Vec::new();
        let mut prompt_life_path = Vec::new();
        let mut averages = Vec::new();

        for _ in 0..n_cycles {
            let mut values: Vec<f64> = Vec::with_capacity(per_cycle_values);
            for _ in 0..per_cycle_lines {
                let line = raw
                    .get(pos)
                    .ok_or_else(|| Error::BadStructure("cycle data truncated".into()))?;
                pos += 1;
                for t in line.split_whitespace() {
                    values.push(num("cycle value", t)?);
                }
            }
            if values.len() < 5 {
                return Err(Error::BadStructure(format!(
                    "cycle row has {} values, need >= 5",
                    values.len()
                )));
            }
            k_col.push(values[0]);
            k_abs.push(values[1]);
            k_path.push(values[2]);
            prompt_life_col.push(values[3]);
            prompt_life_path.push(values[4]);

            if per_cycle_values > 5 {
                if values.len() < 19 {
                    return Err(Error::BadStructure("19-var cycle row incomplete".into()));
                }
                let pair = |i: usize, j: usize| (values[i], values[j]);
                averages.push(CycleAverages {
                    avg_k_col: pair(5, 6),
                    avg_k_abs: pair(7, 8),
                    avg_k_path: pair(9, 10),
                    avg_k_combined: pair(11, 12),
                    avg_k_combined_active: pair(13, 14),
                    prompt_life_combined: pair(15, 16),
                    cycle_histories: values[17],
                    fom: values[18],
                });
            }
        }

        Ok(Mctal {
            code_name,
            code_version,
            code_date,
            code_time,
            n_dump,
            n_histories,
            n_prn,
            comment,
            n_tallies,
            npert,
            tally_nums,
            tallies,
            mesh_tallies,
            n_cycles,
            n_inactive,
            vars_per_cycle,
            k_col,
            k_abs,
            k_path,
            prompt_life_col,
            prompt_life_path,
            averages,
        })
    }
}

/// Skip blank lines at the cursor.
fn skip_blank(raw: &[&str], pos: &mut usize) {
    while raw.get(*pos).map(|l| l.trim().is_empty()).unwrap_or(false) {
        *pos += 1;
    }
}

/// First whitespace token of a line, lowercased.
fn first_token(line: &str) -> String {
    line.split_whitespace()
        .next()
        .unwrap_or("")
        .to_ascii_lowercase()
}

/// True when `line` opens a section with the given keyword: its first token
/// equals the keyword or starts with it (`tfc` data lines excluded by exact
/// match at call sites).
fn opens(line: &str, keyword: &str) -> bool {
    first_token(line) == keyword
}

/// Parsed card head: bin letter, total/cumulative variant, and optional
/// tally number (`d4` → (`d`, None, Some(4)); `ut` → (`u`, Some('t'), None)).
fn split_card_head(head: &str) -> Option<(char, Option<char>, Option<u32>)> {
    let b = head.as_bytes();
    if b.is_empty() {
        return None;
    }
    let letter = (b[0] as char).to_ascii_lowercase();
    if !matches!(letter, 'f' | 'd' | 'u' | 's' | 'm' | 'c' | 'e' | 't') {
        return None;
    }
    let rest = &head[1..];
    if rest.is_empty() {
        return Some((letter, None, None));
    }
    let (variant, digits) = match rest.as_bytes().first() {
        Some(c) if *c == b't' || *c == b'T' || *c == b'c' || *c == b'C' => {
            (Some((*c as char).to_ascii_lowercase()), &rest[1..])
        }
        _ => (None, rest),
    };
    if digits.is_empty() {
        return Some((letter, variant, None));
    }
    if !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let n: u32 = digits.parse().ok()?;
    Some((letter, variant, Some(n)))
}

/// True when the lowercased first token looks like a bin card: a bin letter
/// with an optional `t`/`c` variant and an optional tally number (`d`,
/// `d4`, `ut`, `ut4`, …). Exact keywords (`tfc`, `tally`, `kcode`, `vals`)
/// never classify as bin cards.
fn is_bin_card(head: &str) -> bool {
    if matches!(head, "tfc" | "tally" | "kcode" | "vals") {
        return false;
    }
    split_card_head(head).is_some()
}

/// Parse one bin card plus its values.
///
/// Accepts `<letter>[t|c][<tally>] <count> [flag]`: bare (`d`) or numbered
/// (`d4`) heads, total (`t`) / cumulative (`c`) variants, and a third flag
/// token on `c`/`e`/`t` cards (stored verbatim). Then collects exactly
/// `count` floats from the following lines (blank lines skipped). A `count`
/// of 0 carries no values. Stops collecting when the next section keyword
/// is reached; anything else is a named-open layout.
fn parse_bin_card(
    raw: &[&str],
    pos: &mut usize,
    letter: char,
    tally: u32,
    context: &'static str,
) -> Result<BinCard, Error> {
    skip_blank(raw, pos);
    let line = raw.get(*pos).ok_or_else(|| {
        Error::BadStructure(format!("tally {tally}: truncated before `{letter}` card"))
    })?;
    let toks: Vec<&str> = line.split_whitespace().collect();
    let got = toks.first().copied().unwrap_or("").to_ascii_lowercase();
    let (got_letter, variant, got_num) = split_card_head(&got).ok_or_else(|| {
        Error::BadStructure(format!(
            "tally {tally}: expected `{letter}{tally}` card, got `{line}`"
        ))
    })?;
    if got_letter != letter {
        return Err(Error::BadStructure(format!(
            "tally {tally}: expected `{letter}{tally}` card, got `{line}`"
        )));
    }
    if let Some(n) = got_num {
        if n != tally {
            return Err(Error::BadStructure(format!(
                "tally {tally}: expected `{letter}{tally}` card, got `{line}`"
            )));
        }
    }
    if toks.len() > 3 {
        return Err(Error::BadStructure(format!(
            "named-open: tally {tally}: `{got}` has too many tokens `{line}`"
        )));
    }
    let flag = if toks.len() == 3 {
        if !matches!(letter, 'c' | 'e' | 't') {
            return Err(Error::BadStructure(format!(
                "named-open: tally {tally}: `{got}` flag tokens `{line}` are not parsed"
            )));
        }
        let f = num(context, toks[2])?;
        if !(f.is_finite() && f.fract() == 0.0) {
            return Err(Error::BadStructure(format!(
                "tally {tally}: `{got}` flag `{}` is not an integer",
                toks[2]
            )));
        }
        // Flags ride the file as small ints; clamp loudly instead of
        // wrapping on hostile input.
        if f < i32::MIN as f64 || f > i32::MAX as f64 {
            return Err(Error::BadStructure(format!(
                "tally {tally}: `{got}` flag `{}` out of range",
                toks[2]
            )));
        }
        Some(f as i32)
    } else {
        None
    };
    let count_tok = toks.get(1).copied().unwrap_or("0");
    // Bin counts come from the file: require a non-negative integer (float
    // truncation and wraparound would silently misparse), and never reserve
    // from the count — allocation follows the bytes actually present.
    let count_f = num(context, count_tok)?;
    if !(count_f.is_finite() && count_f >= 0.0 && count_f.fract() == 0.0) {
        return Err(Error::BadStructure(format!(
            "tally {tally}: `{got}` count `{count_tok}` is not a non-negative integer"
        )));
    }
    let count = count_f as usize;
    *pos += 1;
    let mut values = Vec::new();
    while values.len() < count {
        skip_blank(raw, pos);
        let peek = raw.get(*pos).ok_or_else(|| {
            Error::BadStructure(format!("tally {tally}: truncated in `{letter}` values"))
        })?;
        let head = first_token(peek);
        if head == "vals"
            || head == "tfc"
            || head == "tally"
            || head == "kcode"
            || is_bin_card(&head)
        {
            return Err(Error::BadStructure(format!(
                "tally {tally}: `{letter}` card declares {count} values, found {}",
                values.len()
            )));
        }
        for t in peek.split_whitespace() {
            if values.len() == count {
                break;
            }
            values.push(num(context, t)?);
        }
        *pos += 1;
    }
    Ok(BinCard {
        count,
        values,
        variant,
        flag,
    })
}

/// Parse the optional `tfc` block at the cursor.
///
/// The cursor must sit past a tally's `vals`. When the next non-blank line
/// opens `tfc`, consumes the `tfc` jtf line (exactly 9 ints) plus the
/// following 3–4-float data rows until the next `tally`/`kcode` block;
/// otherwise returns `None` without advancing past blanks.
///
/// Interim: not yet wired into `parse_tally_body` (the `tfc` owner wires it
/// when the block stops being named-open).
#[allow(dead_code)]
fn parse_optional_tfc(
    raw: &[&str],
    pos: &mut usize,
    tally: u32,
) -> Result<Option<TfcBlock>, Error> {
    let save = *pos;
    skip_blank(raw, pos);
    let peek = match raw.get(*pos) {
        Some(l) => *l,
        None => {
            *pos = save;
            return Ok(None);
        }
    };
    if first_token(peek) != "tfc" {
        *pos = save;
        return Ok(None);
    }
    let toks: Vec<&str> = peek.split_whitespace().collect();
    if toks.len() != 10 {
        return Err(Error::BadStructure(format!(
            "tally {tally}: `tfc` jtf line needs 9 ints, got `{peek}`"
        )));
    }
    let mut jtf = Vec::with_capacity(9);
    for tok in &toks[1..] {
        let v: i64 = tok.parse().map_err(|_| Error::BadNumber {
            context: "tfc jtf",
            text: (*tok).to_string(),
        })?;
        jtf.push(v);
    }
    *pos += 1;
    let mut rows = Vec::new();
    loop {
        skip_blank(raw, pos);
        let line = match raw.get(*pos) {
            Some(l) => *l,
            None => break,
        };
        let head = first_token(line);
        if head == "tally" || head == "kcode" {
            break;
        }
        if head == "tfc" || head == "vals" || is_bin_card(&head) {
            return Err(Error::BadStructure(format!(
                "tally {tally}: unexpected `{head}` inside `tfc` data"
            )));
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 3 && toks.len() != 4 {
            return Err(Error::BadStructure(format!(
                "tally {tally}: `tfc` data rows need 3-4 floats, got `{line}`"
            )));
        }
        let nps: i64 = toks[0].parse().map_err(|_| Error::BadNumber {
            context: "tfc nps",
            text: toks[0].to_string(),
        })?;
        let value = num("tfc value", toks[1])?;
        let rel_err = num("tfc rel error", toks[2])?;
        let fom = if toks.len() == 4 {
            Some(num("tfc fom", toks[3])?)
        } else {
            None
        };
        rows.push(TfcEntry {
            nps,
            value,
            rel_err,
            fom,
        });
        *pos += 1;
    }
    Ok(Some(TfcBlock { jtf, rows }))
}

/// One parsed tally block: standard or mesh.
enum TallyBlock {
    Standard(TallyBody),
    Mesh(MeshTallyBody),
}

/// Read the `tally` head line at the cursor (cursor must sit on it) without
/// advancing past anything else.
fn read_tally_head(raw: &[&str], pos: usize) -> Result<(u32, i32, Option<i32>), Error> {
    let line = raw
        .get(pos)
        .ok_or_else(|| Error::BadStructure("truncated tally body".into()))?;
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.first() != Some(&"tally") {
        return Err(Error::BadStructure(format!(
            "expected `tally` line, got `{line}`"
        )));
    }
    if toks.len() < 3 {
        return Err(Error::BadStructure(format!(
            "tally line needs number + particle type, got `{line}`"
        )));
    }
    let number = num("tally number", toks[1])? as u32;
    let particle_type = num("tally particle type", toks[2])? as i32;
    let detector_type = toks
        .get(3)
        .map(|t| num("tally detector type", t).map(|v| v as i32))
        .transpose()?;
    Ok((number, particle_type, detector_type))
}

/// Parse one integer field verbatim from the file (no float truncation,
/// no wraparound: finite, integral, and inside `i64`).
fn int_field(context: &'static str, tok: &str) -> Result<i64, Error> {
    let v = num(context, tok)?;
    if !(v.is_finite() && v.fract() == 0.0) {
        return Err(Error::BadStructure(format!(
            "{context} `{tok}` is not an integer"
        )));
    }
    if v < i64::MIN as f64 || v > i64::MAX as f64 {
        return Err(Error::BadStructure(format!(
            "{context} `{tok}` is out of range"
        )));
    }
    Ok(v as i64)
}

/// Parse one mesh-tally body at the cursor (cursor must sit on its `tally`
/// line with `detector_type <= -1`; blank lines already skipped by the
/// caller).
///
/// The `f` line carries the 4-int mesh info (`f[<tally>] <unknown> <ni>
/// <nj> <nk>`, bare-`f` accepted), followed by `(ni+1)+(nj+1)+(nk+1)`
/// cora/b/c bound floats across as many lines as needed (blank lines
/// skipped; a section keyword ends the bounds run and trips a truncation
/// error). Bounds are partitioned in file order: the first `ni+1` values
/// are `cora`, the next `nj+1` are `corb`, the rest are `corc`. The
/// remaining cards, `vals` pairing, and total/cumulative/flag spellings
/// match the standard body; mesh tallies carry no `tfc` block (a trailing
/// `tfc` is a loud named-open error, not a silent skip).
fn parse_mesh_body(raw: &[&str], pos: &mut usize) -> Result<MeshTallyBody, Error> {
    let (number, particle_type, detector_type) = read_tally_head(raw, *pos)?;
    let mesh_kind = detector_type.unwrap_or(-1);
    *pos += 1;

    let particle_list = if particle_type < 0 {
        read_particle_list(raw, pos, number)?
    } else {
        Vec::new()
    };

    let comment = read_comments(raw, pos, number)?;

    // Mesh `f` line: head plus exactly 4 ints.
    skip_blank(raw, pos);
    let fline = raw.get(*pos).ok_or_else(|| {
        Error::BadStructure(format!("tally {number}: truncated before mesh `f` line"))
    })?;
    let ftoks: Vec<&str> = fline.split_whitespace().collect();
    let fhead = ftoks.first().copied().unwrap_or("").to_ascii_lowercase();
    let (fletter, fvariant, fnum) = split_card_head(&fhead).ok_or_else(|| {
        Error::BadStructure(format!(
            "tally {number}: expected mesh `f{number}` line, got `{fline}`"
        ))
    })?;
    if fletter != 'f' {
        return Err(Error::BadStructure(format!(
            "tally {number}: expected mesh `f{number}` line, got `{fline}`"
        )));
    }
    if let Some(n) = fnum {
        if n != number {
            return Err(Error::BadStructure(format!(
                "tally {number}: expected mesh `f{number}` line, got `{fline}`"
            )));
        }
    }
    if fvariant.is_some() {
        return Err(Error::BadStructure(format!(
            "named-open: tally {number}: mesh total/cumulative `f` variants are not parsed"
        )));
    }
    if ftoks.len() != 5 {
        return Err(Error::BadStructure(format!(
            "tally {number}: mesh `f` line needs 4 ints (`unknown ni nj nk`), got `{fline}`"
        )));
    }
    let mesh_unknown = int_field("mesh unknown", ftoks[1])?;
    let mut dims = [0usize; 3];
    for (i, tok) in ftoks[2..].iter().enumerate() {
        let v = int_field("mesh bin count", tok)?;
        if v < 1 {
            return Err(Error::BadStructure(format!(
                "tally {number}: mesh axis {i} needs >= 1 bins, got `{tok}`"
            )));
        }
        dims[i] = v as usize;
    }
    let (ni, nj, nk) = (dims[0], dims[1], dims[2]);
    *pos += 1;

    // Bounds: allocation follows file bytes, never the declared counts.
    let need_bounds = (ni + 1)
        .checked_add(nj + 1)
        .and_then(|s| s.checked_add(nk + 1))
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: mesh bounds overflow")))?;
    let mut flat: Vec<f64> = Vec::new();
    while flat.len() < need_bounds {
        skip_blank(raw, pos);
        let peek = raw.get(*pos).ok_or_else(|| {
            Error::BadStructure(format!("tally {number}: truncated in mesh bounds"))
        })?;
        let head = first_token(peek);
        if head == "tally" || head == "kcode" || head == "vals" || head == "tfc" {
            return Err(Error::BadStructure(format!(
                "tally {number}: mesh bounds need {need_bounds} floats, found {}",
                flat.len()
            )));
        }
        if is_bin_card(&head) {
            return Err(Error::BadStructure(format!(
                "tally {number}: mesh bounds need {need_bounds} floats, found {} before `{head}`",
                flat.len()
            )));
        }
        for tok in peek.split_whitespace() {
            if flat.len() == need_bounds {
                break;
            }
            flat.push(num("mesh bound", tok)?);
        }
        *pos += 1;
    }
    let corc_split = flat.len() - (nk + 1);
    let corb_split = corc_split - (nj + 1);
    let corc = flat[corc_split..].to_vec();
    let corb = flat[corb_split..corc_split].to_vec();
    let cora = flat[..corb_split].to_vec();

    let d = parse_bin_card(raw, pos, 'd', number, "tally direct bins")?;
    let u = parse_bin_card(raw, pos, 'u', number, "tally user bins")?;
    let s = parse_bin_card(raw, pos, 's', number, "tally segment bins")?;
    let m = parse_bin_card(raw, pos, 'm', number, "tally multiplier bins")?;
    let c = parse_bin_card(raw, pos, 'c', number, "tally cosine bins")?;
    let e = parse_bin_card(raw, pos, 'e', number, "tally energy bins")?;
    let t = parse_bin_card(raw, pos, 't', number, "tally time bins")?;

    let probe = MeshTallyBody {
        number,
        particle_type,
        detector_type: mesh_kind,
        particle_list: particle_list.clone(),
        comment: Vec::new(),
        mesh_unknown,
        ni,
        nj,
        nk,
        cora: Vec::new(),
        corb: Vec::new(),
        corc: Vec::new(),
        d: d.clone(),
        u: u.clone(),
        s: s.clone(),
        m: m.clone(),
        c: c.clone(),
        e: e.clone(),
        t: t.clone(),
        vals: Vec::new(),
    };
    let vals = read_vals_block(raw, pos, number, probe.expected_pairs())?;

    // Mesh tallies carry no `tfc`: a trailing block is named-open, loud.
    skip_blank(raw, pos);
    if let Some(peek) = raw.get(*pos) {
        if opens(peek, "tfc") {
            return Err(Error::BadStructure(format!(
                "named-open: tally {number}: mesh-tally `tfc` blocks are not parsed"
            )));
        }
    }

    Ok(MeshTallyBody {
        number,
        particle_type,
        detector_type: mesh_kind,
        particle_list,
        comment,
        mesh_unknown,
        ni,
        nj,
        nk,
        cora,
        corb,
        corc,
        d,
        u,
        s,
        m,
        c,
        e,
        t,
        vals,
    })
}

/// Parse one tally block at the cursor (cursor must sit on its `tally`
/// line; blank lines already skipped by the caller), dispatching to the
/// standard or mesh body parser on `detector_type`.
fn parse_tally_block(raw: &[&str], pos: &mut usize) -> Result<TallyBlock, Error> {
    let (number, _particle_type, detector_type) = read_tally_head(raw, *pos)?;
    if let Some(d) = detector_type {
        if d <= -1 {
            return parse_mesh_body(raw, pos).map(TallyBlock::Mesh);
        }
        if d >= 3 {
            return Err(Error::BadStructure(format!(
                "named-open: tally {number}: radiograph tallies (detector_type {d}) are not parsed"
            )));
        }
    }
    if number % 10 == 5 {
        return Err(Error::BadStructure(format!(
            "named-open: tally {number}: point-detector tallies are not parsed"
        )));
    }
    parse_tally_body(raw, pos).map(TallyBlock::Standard)
}

/// Read one particle-list line for a negative `particle_type`.
fn read_particle_list(raw: &[&str], pos: &mut usize, number: u32) -> Result<Vec<i32>, Error> {
    skip_blank(raw, pos);
    let pl = raw
        .get(*pos)
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: truncated particle list")))?;
    let head = first_token(pl);
    if head == "tally" || head == "kcode" || head == "vals" || head == "tfc" || is_bin_card(&head) {
        return Err(Error::BadStructure(format!(
            "tally {number}: missing particle list for negative particle type"
        )));
    }
    let mut out = Vec::new();
    for t in pl.split_whitespace() {
        out.push(num("tally particle list", t)? as i32);
    }
    *pos += 1;
    Ok(out)
}

/// Collect `FC` comment lines at the cursor (anything before the `f` card
/// that is not itself a section keyword).
fn read_comments(raw: &[&str], pos: &mut usize, number: u32) -> Result<Vec<String>, Error> {
    let mut comment = Vec::new();
    loop {
        skip_blank(raw, pos);
        let peek = raw.get(*pos).ok_or_else(|| {
            Error::BadStructure(format!("tally {number}: truncated before `f` card"))
        })?;
        let head = first_token(peek);
        if is_bin_card(&head)
            || head == "vals"
            || head == "tfc"
            || head == "tally"
            || head == "kcode"
        {
            break;
        }
        comment.push(peek.trim().to_string());
        *pos += 1;
    }
    Ok(comment)
}

/// Read a `vals` block of exactly `need_pairs` val/err pairs at the cursor
/// (cursor must sit past the `t` card).
fn read_vals_block(
    raw: &[&str],
    pos: &mut usize,
    number: u32,
    need_pairs: usize,
) -> Result<Vec<(f64, f64)>, Error> {
    skip_blank(raw, pos);
    let vline = raw
        .get(*pos)
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: truncated before `vals`")))?;
    if !opens(vline, "vals") {
        return Err(Error::BadStructure(format!(
            "tally {number}: expected `vals` block, got `{vline}`"
        )));
    }
    *pos += 1;
    let need_floats = need_pairs
        .checked_mul(2)
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: bin product overflows")))?;
    let mut flat = Vec::new();
    while flat.len() < need_floats {
        skip_blank(raw, pos);
        let peek = raw
            .get(*pos)
            .ok_or_else(|| Error::BadStructure(format!("tally {number}: truncated in `vals`")))?;
        let head = first_token(peek);
        if head == "tally" || head == "kcode" {
            return Err(Error::BadStructure(format!(
                "tally {number}: `vals` needs {need_floats} floats, found {}",
                flat.len()
            )));
        }
        if head == "tfc" || is_bin_card(&head) || head == "vals" {
            return Err(Error::BadStructure(format!(
                "tally {number}: `vals` needs {need_floats} floats, found {} before `{head}`",
                flat.len()
            )));
        }
        for tok in peek.split_whitespace() {
            if flat.len() == need_floats {
                break;
            }
            flat.push(num("tally vals value", tok)?);
        }
        *pos += 1;
    }
    // `flat.len()` is exact here, so this reservation is file-bounded.
    let mut vals = Vec::with_capacity(flat.len() / 2);
    for pair in flat.chunks_exact(2) {
        vals.push((pair[0], pair[1]));
    }
    Ok(vals)
}

/// Parse one standard-tally body at the cursor (cursor must sit on its
/// `tally` line; blank lines already skipped by the caller).
fn parse_tally_body(raw: &[&str], pos: &mut usize) -> Result<TallyBody, Error> {
    let (number, particle_type, detector_type) = read_tally_head(raw, *pos)?;
    // The dispatcher keeps mesh/radiograph/point-detector arms out of here;
    // re-check defensively so a direct caller still hears a loud named-open.
    if let Some(d) = detector_type {
        if d >= 3 {
            return Err(Error::BadStructure(format!(
                "named-open: tally {number}: radiograph tallies (detector_type {d}) are not parsed"
            )));
        }
    }
    if number % 10 == 5 {
        return Err(Error::BadStructure(format!(
            "named-open: tally {number}: point-detector tallies are not parsed"
        )));
    }
    *pos += 1;

    let particle_list = if particle_type < 0 {
        read_particle_list(raw, pos, number)?
    } else {
        Vec::new()
    };

    let comment = read_comments(raw, pos, number)?;

    let f = parse_bin_card(raw, pos, 'f', number, "tally object bins")?;
    let d = parse_bin_card(raw, pos, 'd', number, "tally direct bins")?;
    let u = parse_bin_card(raw, pos, 'u', number, "tally user bins")?;
    let s = parse_bin_card(raw, pos, 's', number, "tally segment bins")?;
    let m = parse_bin_card(raw, pos, 'm', number, "tally multiplier bins")?;
    let c = parse_bin_card(raw, pos, 'c', number, "tally cosine bins")?;
    let e = parse_bin_card(raw, pos, 'e', number, "tally energy bins")?;
    let t = parse_bin_card(raw, pos, 't', number, "tally time bins")?;

    for (card, bins) in [
        ('f', &f),
        ('d', &d),
        ('u', &u),
        ('s', &s),
        ('m', &m),
        ('c', &c),
        ('e', &e),
        ('t', &t),
    ] {
        if bins.values.len() != bins.count {
            return Err(Error::BadStructure(format!(
                "tally {number}: `{card}` card declares {} values, found {}",
                bins.count,
                bins.values.len()
            )));
        }
    }

    // `vals` block: exactly 2 * prod(bins, 0 -> 1) floats as val/err pairs.
    skip_blank(raw, pos);
    let vline = raw
        .get(*pos)
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: truncated before `vals`")))?;
    if !opens(vline, "vals") {
        if opens(vline, "tfc") {
            return Err(Error::BadStructure(format!(
                "named-open: tally {number}: `tfc` blocks are not parsed"
            )));
        }
        return Err(Error::BadStructure(format!(
            "tally {number}: expected `vals` block, got `{vline}`"
        )));
    }
    *pos += 1;
    let probe = TallyBody {
        number,
        particle_type,
        detector_type,
        particle_list: particle_list.clone(),
        comment: Vec::new(),
        f: f.clone(),
        d: d.clone(),
        u: u.clone(),
        s: s.clone(),
        m: m.clone(),
        c: c.clone(),
        e: e.clone(),
        t: t.clone(),
        vals: Vec::new(),
        tfc: None,
    };
    let need_pairs = probe.expected_pairs();
    let need_floats = need_pairs
        .checked_mul(2)
        .ok_or_else(|| Error::BadStructure(format!("tally {number}: bin product overflows")))?;
    let mut flat = Vec::new();
    while flat.len() < need_floats {
        skip_blank(raw, pos);
        let peek = raw
            .get(*pos)
            .ok_or_else(|| Error::BadStructure(format!("tally {number}: truncated in `vals`")))?;
        let head = first_token(peek);
        if head == "tally" || head == "kcode" {
            return Err(Error::BadStructure(format!(
                "tally {number}: `vals` needs {need_floats} floats, found {}",
                flat.len()
            )));
        }
        if head == "tfc" || is_bin_card(&head) || head == "vals" {
            return Err(Error::BadStructure(format!(
                "tally {number}: `vals` needs {need_floats} floats, found {} before `{head}`",
                flat.len()
            )));
        }
        for tok in peek.split_whitespace() {
            if flat.len() == need_floats {
                break;
            }
            flat.push(num("tally vals value", tok)?);
        }
        *pos += 1;
    }
    // `flat.len()` is exact here, so this reservation is file-bounded.
    let mut vals = Vec::with_capacity(flat.len() / 2);
    for pair in flat.chunks_exact(2) {
        vals.push((pair[0], pair[1]));
    }

    // Optional `tfc` block after `vals` (standard tallies only; mesh
    // tallies reject it in `parse_mesh_body`).
    let tfc = parse_optional_tfc(raw, pos, number)?;

    Ok(TallyBody {
        number,
        particle_type,
        detector_type,
        particle_list,
        comment,
        f,
        d,
        u,
        s,
        m,
        c,
        e,
        t,
        vals,
        tfc,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/mctal/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn kcode5_header_and_cycles() {
        let m = Mctal::from_file(fixture("synthetic_kcode5.mctal")).unwrap();
        assert_eq!(m.code_name, "mcnp");
        assert_eq!(m.code_version, "6.2.0");
        assert_eq!(m.code_date, "05/08/13");
        assert_eq!(m.n_dump, "1");
        assert_eq!(m.n_histories, 100000);
        assert_eq!(m.n_prn, 5);
        assert_eq!(m.comment, "Simple MCNP Example that uses SSW and kcode");
        assert_eq!(m.n_tallies, "0");
        assert_eq!(m.tally_nums, Vec::<u32>::new());
        assert_eq!(m.n_cycles, 8);
        assert_eq!(m.n_inactive, 4);
        assert_eq!(m.vars_per_cycle, 5);
        assert_eq!(m.k_col.len(), 8);
        // Deterministic synthetic values: kc = 0.98 + 0.005*c
        assert_eq!(m.k_col[0], 0.985);
        assert_eq!(m.k_col[7], 1.02);
        assert!((m.prompt_life_col[3] - 5.0e-4 * 4.0).abs() < 1e-12);
        assert!(m.averages.is_empty());
        // Legacy kcode-only file: no bodies parsed.
        assert!(m.tallies.is_empty());
        assert_eq!(m.npert, None);
    }

    #[test]
    fn kcode19_running_averages() {
        let m = Mctal::from_file(fixture("synthetic_kcode19.mctal")).unwrap();
        assert_eq!(m.n_tallies, "2");
        assert_eq!(m.tally_nums, vec![4, 14]);
        assert_eq!(m.n_cycles, 6);
        assert_eq!(m.vars_per_cycle, 19);
        assert_eq!(m.averages.len(), 6);
        let a0 = &m.averages[0];
        assert_eq!(a0.avg_k_col.0, 1.0021);
        assert_eq!(a0.avg_k_col.1, 2e-5);
        assert_eq!(a0.cycle_histories, 5001.0);
        assert_eq!(a0.fom, 42.0);
        assert_eq!(m.averages[5].avg_k_combined_active.0, 1.0125);
        // Declares tallies but carries no bodies (legacy kcode-only path).
        assert!(m.tallies.is_empty());
        assert_eq!(m.npert, None);
    }

    #[test]
    fn missing_kcode_errors() {
        let text = "mcnp 6.2.0 d t 1 10 5\ncomment\ntally 0\n\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("kcode")
        ));
    }

    #[test]
    fn short_header_errors() {
        assert!(matches!(
            Mctal::parse("only one line"),
            Err(Error::BadStructure(_))
        ));
    }

    #[test]
    fn truncated_cycle_data_errors() {
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n\nkcode 3 1 0 5\n 1 1 1 1 1\n 2 2 2 2 2\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("truncated")
        ));
    }

    #[test]
    fn bad_number_reports_context() {
        let text = "mcnp v d t 1 xx 5\nc\ntally 0\n\nkcode 1 0 0 5\n1 1 1 1 1\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadNumber {
                context: "n_histories",
                ..
            })
        ));
    }

    #[test]
    fn error_display_messages() {
        use std::error::Error as _;

        let io = Error::Io("disk".into());
        assert_eq!(io.to_string(), "io error: disk");
        assert!(io.source().is_none());
        assert_eq!(
            Error::BadStructure("no kcode".into()).to_string(),
            "malformed MCTAL: no kcode"
        );
        assert_eq!(
            Error::BadNumber {
                context: "n_cycles",
                text: "xx".into(),
            }
            .to_string(),
            "cannot parse n_cycles from `xx`"
        );
    }

    #[test]
    fn bad_tally_line_errors() {
        let text = "mcnp v d t 1 10 5\nc\nnot_a_tally 0\n\nkcode 1 0 0 5\n1 1 1 1 1\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("expected `tally` line")
        ));
    }

    #[test]
    fn eof_after_tally_line_still_searches_kcode() {
        // No tally-numbers line at all: tally_nums defaults to empty, then
        // the kcode search hits EOF.
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("kcode")
        ));
    }

    #[test]
    fn junk_lines_before_kcode_are_skipped() {
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n4 14\nversion junk\nmore junk\nkcode 2 1 0 5\n 1 1 1 1 1\n 2 2 2 2 2\n";
        let m = Mctal::parse(text).unwrap();
        assert_eq!(m.tally_nums, vec![4, 14]);
        assert_eq!(m.n_tallies, "0");
        assert_eq!(m.k_col, vec![1.0, 2.0]);
    }

    #[test]
    fn unsupported_vars_per_cycle_errors() {
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n\nkcode 1 0 7 5\n1 1 1 1 1\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("unsupported vars_per_cycle 7")
        ));
    }

    #[test]
    fn short_cycle_row_errors() {
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n\nkcode 1 0 5 5\n1 2 3\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("cycle row has 3 values, need >= 5")
        ));
    }

    #[test]
    fn incomplete_19var_row_errors() {
        // Four physical lines per cycle, but only 10 values total.
        let text = "mcnp v d t 1 10 5\nc\ntally 0\n\nkcode 1 0 19\n1 2 3 4 5 6 7 8 9 10\n\n\n\n";
        assert!(matches!(
            Mctal::parse(text),
            Err(Error::BadStructure(m)) if m.contains("19-var cycle row incomplete")
        ));
    }

    #[test]
    fn body_fixture_closed_form_pairs() {
        let m = Mctal::from_file(fixture("synthetic_tally_bodies.mctal")).unwrap();
        // Headers/cycles keep PyNE parity (negative oracle surface).
        assert_eq!(m.code_name, "mcnp");
        assert_eq!(m.n_tallies, "2");
        assert_eq!(m.npert, None);
        assert_eq!(m.tally_nums, vec![4, 14]);
        assert_eq!(m.n_cycles, 2);
        assert_eq!(m.k_col, vec![0.99, 1.01]);

        assert_eq!(m.tallies.len(), 2);
        let t4 = &m.tallies[0];
        assert_eq!(t4.number, 4);
        assert_eq!(t4.particle_type, 1);
        assert_eq!(t4.detector_type, Some(0));
        assert!(t4.particle_list.is_empty());
        assert!(t4.comment.is_empty());
        assert_eq!(t4.f.count, 2);
        assert_eq!(t4.f.values, vec![1.0, 2.0]);
        assert_eq!(t4.d.count, 0);
        assert!(t4.d.values.is_empty());
        assert_eq!(t4.e.count, 2);
        assert_eq!(t4.e.values, vec![0.5, 2.0]);
        assert_eq!(t4.expected_pairs(), 4);
        // Closed-form pairing: val = 10*obj + epos, err = val/8, f outer.
        assert_eq!(
            t4.vals,
            vec![(11.0, 1.375), (12.0, 1.5), (21.0, 2.625), (22.0, 2.75)]
        );
        assert_eq!(t4.total_val(), 66.0);

        let t14 = &m.tallies[1];
        assert_eq!(t14.number, 14);
        assert_eq!(t14.f.values, vec![5.0]);
        assert_eq!(t14.c.count, 1);
        assert_eq!(t14.c.values, vec![0.5]);
        assert_eq!(t14.expected_pairs(), 1);
        assert_eq!(t14.vals, vec![(7.0, 0.875)]);
        assert_eq!(t14.total_val(), 7.0);
    }

    /// Minimal one-tally wrapper used by the body error tests: header plus
    /// `middle` (tally body text) plus a one-cycle kcode tail.
    fn body_probe(middle: &str) -> String {
        format!("mcnp v d t 1 10 5\nc\ntally 1\n4\n{middle}kcode 1 0 5\n1 1 1 1 1\n")
    }

    const MIN_BODY: &str =
        "tally 4 1 0\nf4 1\n 1\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\ne4 0\nt4 0\nvals\n 3.0 0.375\n";

    #[test]
    fn body_minimal_single_bin() {
        let m = Mctal::parse(&body_probe(MIN_BODY)).unwrap();
        assert_eq!(m.tallies.len(), 1);
        assert_eq!(m.tallies[0].expected_pairs(), 1);
        assert_eq!(m.tallies[0].vals, vec![(3.0, 0.375)]);
    }

    #[test]
    fn body_negative_particle_type_reads_list_line() {
        let middle = "tally 24 -1 0\n 1 0\nf24 1\n 1\nd24 0\nu24 0\ns24 0\nm24 0\nc24 0\ne24 0\nt24 0\nvals\n 3.0 0.375\n";
        let text = format!("mcnp v d t 1 10 5\nc\ntally 1\n24\n{middle}kcode 1 0 5\n1 1 1 1 1\n");
        let m = Mctal::parse(&text).unwrap();
        assert_eq!(m.tallies[0].particle_list, vec![1, 0]);
    }

    #[test]
    fn body_fc_comments_collected() {
        let middle =
            "tally 4 1 0\nthis is an FC comment\nf4 1\n 1\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\ne4 0\nt4 0\nvals\n 3.0 0.375\n";
        let m = Mctal::parse(&body_probe(middle)).unwrap();
        assert_eq!(
            m.tallies[0].comment,
            vec!["this is an FC comment".to_string()]
        );
    }

    #[test]
    fn body_tfc_parsed() {
        // `tfc` jtf line (exactly 9 ints) plus 3-float and 4-float rows.
        let middle =
            format!("{MIN_BODY}tfc 1 2 3 4 5 6 7 8 9\n 1000 11.5 0.2\n 2000 11.75 0.15 3.5\n");
        let m = Mctal::parse(&body_probe(&middle)).unwrap();
        let tfc = m.tallies[0].tfc.as_ref().expect("tfc block parsed");
        assert_eq!(tfc.jtf, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(tfc.rows.len(), 2);
        assert_eq!(tfc.rows[0].nps, 1000);
        assert_eq!((tfc.rows[0].value, tfc.rows[0].rel_err), (11.5, 0.2));
        assert_eq!(tfc.rows[0].fom, None);
        assert_eq!(tfc.rows[1].fom, Some(3.5));
        // No `tfc` lines means no block (cursor left before `kcode`).
        let plain = Mctal::parse(&body_probe(MIN_BODY)).unwrap();
        assert!(plain.tallies[0].tfc.is_none());
    }

    #[test]
    fn body_tfc_bad_jtf_errors() {
        let middle = format!("{MIN_BODY}tfc 1 2 3\n");
        assert!(matches!(
            Mctal::parse(&body_probe(&middle)),
            Err(Error::BadStructure(m)) if m.contains("jtf line needs 9 ints")
        ));
    }

    /// Minimal rectangular mesh-tally body (`detector_type -1`, 2x1x1
    /// cells, closed-form val/err pairing like the mesh fixture).
    const MIN_MESH_BODY: &str = "tally 4 1 -1\nf4 0 2 1 1\n 0.0 5.0 10.0\n 0.0 10.0\n 0.0 10.0\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\ne4 0\nt4 0\nvals\n 11.0 1.375 12.0 1.5\n";

    #[test]
    fn mesh_body_parsed() {
        let m = Mctal::parse(&body_probe(MIN_MESH_BODY)).unwrap();
        assert!(m.tallies.is_empty());
        assert_eq!(m.mesh_tallies.len(), 1);
        let b = &m.mesh_tallies[0];
        assert_eq!(b.number, 4);
        assert_eq!(b.detector_type, -1);
        assert_eq!(b.mesh_unknown, 0);
        assert_eq!(b.dims(), [2, 1, 1]);
        assert_eq!(b.num_cells(), 2);
        assert_eq!(b.cora, vec![0.0, 5.0, 10.0]);
        assert_eq!(b.corb, vec![0.0, 10.0]);
        assert_eq!(b.corc, vec![0.0, 10.0]);
        assert_eq!(b.expected_pairs(), 2);
        assert_eq!(b.vals, vec![(11.0, 1.375), (12.0, 1.5)]);
        assert_eq!(b.total_val(), 23.0);
        assert_eq!(b.mesh_index(1, 0, 0), 1);
    }

    #[test]
    fn mesh_body_errors_are_loud() {
        // Short bounds run: 7 floats needed, section keyword ends the run.
        let short = MIN_MESH_BODY.replace(" 0.0 10.0\n 0.0 10.0\nd4", "d4");
        assert!(matches!(
            Mctal::parse(&body_probe(&short)),
            Err(Error::BadStructure(m)) if m.contains("mesh bounds need 7 floats")
        ));
        // Non-integral mesh dimension.
        let bad_dim = MIN_MESH_BODY.replace("f4 0 2 1 1", "f4 0 2.5 1 1");
        assert!(matches!(
            Mctal::parse(&body_probe(&bad_dim)),
            Err(Error::BadStructure(m)) if m.contains("is not an integer")
        ));
        // Mesh tallies carry no `tfc`.
        let with_tfc = format!("{MIN_MESH_BODY}tfc 1 2 3 4 5 6 7 8 9\n 1000 1.0 0.1\n");
        assert!(matches!(
            Mctal::parse(&body_probe(&with_tfc)),
            Err(Error::BadStructure(m)) if m.contains("named-open") && m.contains("mesh-tally `tfc`")
        ));
    }

    #[test]
    fn body_naked_multi_token_f_is_named_open() {
        // A multi-token `f` line on a standard (non-mesh) tally is mesh
        // info without a mesh detector type: loud, never parsed as flags.
        let middle =
            "tally 4 1 0\nf4 1 2 3 4 5\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\ne4 0\nt4 0\nvals\n 3.0 0.375\n";
        assert!(matches!(
            Mctal::parse(&body_probe(middle)),
            Err(Error::BadStructure(m)) if m.contains("named-open") && m.contains("too many tokens")
        ));
    }

    #[test]
    fn body_point_detector_is_named_open() {
        let middle =
            "tally 5 1 2\nf5 1\n 1\nd5 0\nu5 0\ns5 0\nm5 0\nc5 0\ne5 0\nt5 0\nvals\n 3.0 0.375\n";
        let text = format!("mcnp v d t 1 10 5\nc\ntally 1\n5\n{middle}kcode 1 0 5\n1 1 1 1 1\n");
        assert!(matches!(
            Mctal::parse(&text),
            Err(Error::BadStructure(m)) if m.contains("named-open") && m.contains("point-detector")
        ));
    }

    #[test]
    fn body_total_variant_stored_verbatim() {
        // Total-variant cards (`et4`) now parse with the variant stored
        // verbatim (see `BinCard::variant`); they are no longer named-open.
        let middle =
            "tally 4 1 0\nf4 1\n 1\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\net4 2\n 0.5 2.0\nt4 0\nvals\n 3.0 0.375\n 4.0 0.5\n";
        let m = Mctal::parse(&body_probe(middle)).unwrap();
        assert_eq!(m.tallies[0].e.variant, Some('t'));
    }

    #[test]
    fn mesh_fixture_closed_form_pairs() {
        let m = Mctal::from_file(fixture("synthetic_mesh.mctal")).unwrap();
        assert_eq!(m.tally_nums, vec![4]);
        assert!(m.tallies.is_empty());
        assert_eq!(m.mesh_tallies.len(), 1);
        let b = &m.mesh_tallies[0];
        assert_eq!((b.number, b.detector_type, b.mesh_unknown), (4, -1, 0));
        assert_eq!(b.dims(), [2, 1, 1]);
        assert_eq!(b.cora, vec![0.0, 5.0, 10.0]);
        assert_eq!(b.corb, vec![0.0, 10.0]);
        assert_eq!(b.corc, vec![0.0, 10.0]);
        // Closed-form pairing: val = 11 + cell, err = val/8.
        assert_eq!(b.vals, vec![(11.0, 1.375), (12.0, 1.5)]);
        assert_eq!(b.total_val(), 23.0);
        assert_eq!(m.k_col, vec![0.99, 1.01]);
    }

    #[test]
    fn tfc_variant_fixture() {
        let m = Mctal::from_file(fixture("synthetic_tfc_variants.mctal")).unwrap();
        assert_eq!(m.tally_nums, vec![6]);
        assert_eq!(m.tallies.len(), 1);
        let t = &m.tallies[0];
        // Total-variant energy card and flag-carrying time card stored
        // verbatim; pair count already includes the total bin.
        assert_eq!(t.e.variant, Some('t'));
        assert_eq!(t.e.values, vec![0.5, 2.0]);
        assert_eq!(t.t.flag, Some(0));
        assert_eq!(t.vals, vec![(11.0, 1.375), (12.0, 1.5)]);
        let tfc = t.tfc.as_ref().expect("tfc block parsed");
        assert_eq!(tfc.jtf, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(tfc.rows.len(), 2);
        assert_eq!(tfc.rows[1].fom, Some(3.5));
        assert!(m.mesh_tallies.is_empty());
    }

    #[test]
    fn body_vals_truncation_errors() {
        // e4 declares 2 energies (4 pairs) but only 1 pair is present.
        let middle =
            "tally 4 1 0\nf4 2\n 1 2\nd4 0\nu4 0\ns4 0\nm4 0\nc4 0\ne4 2\n 0.5 2.0\nt4 0\nvals\n 11.0 1.375\nkcode 1 0 5\n";
        assert!(matches!(
            Mctal::parse(&body_probe(middle)),
            Err(Error::BadStructure(m)) if m.contains("`vals` needs 8 floats")
        ));
    }

    #[test]
    fn body_hostile_counts_fail_without_huge_allocs() {
        // A hostile `f` count far beyond the bytes present must fail fast
        // (allocation follows file bytes, never the count): the next card
        // boundary trips the shortfall guard.
        let huge = MIN_BODY.replace("f4 1", "f4 1000000000000");
        assert!(matches!(
            Mctal::parse(&body_probe(&huge)),
            Err(Error::BadStructure(m)) if m.contains("declares 1000000000000 values, found 1")
        ));
        // Fractional and negative counts are malformed, not truncated.
        for bad in ["f4 2.5", "f4 -3"] {
            let text = MIN_BODY.replace("f4 1", bad);
            assert!(matches!(
                Mctal::parse(&body_probe(&text)),
                Err(Error::BadStructure(m)) if m.contains("not a non-negative integer")
            ));
        }
        // A hostile kcode cycle count likewise fails on truncation.
        let kc = "mcnp 6.2.0 05/08/13 17:50:49 3 50000 5\ncomment\ntally 0\n\nkcode 1000000000 0 5\n 0.99 0.99 0.99 5.0e-4 5.0e-4\n";
        assert!(matches!(
            Mctal::parse(kc),
            Err(Error::BadStructure(m)) if m.contains("cycle data truncated")
        ));
    }

    #[test]
    fn expected_pairs_saturates_instead_of_wrapping() {
        let big = BinCard {
            count: usize::MAX / 2,
            values: Vec::new(),
            variant: None,
            flag: None,
        };
        let body = TallyBody {
            number: 4,
            particle_type: 1,
            detector_type: None,
            particle_list: Vec::new(),
            comment: Vec::new(),
            f: big.clone(),
            d: big.clone(),
            u: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            s: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            m: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            c: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            e: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            t: BinCard {
                count: 0,
                values: Vec::new(),
                variant: None,
                flag: None,
            },
            vals: Vec::new(),
            tfc: None,
        };
        assert_eq!(body.expected_pairs(), usize::MAX);
    }

    #[test]
    fn body_unknown_tally_number_errors() {
        // Body for tally 6 while only 4 is declared (6 avoids the
        // point-detector named-open path so the membership check fires).
        let body6 = MIN_BODY
            .replace("tally 4", "tally 6")
            .replace("f4", "f6")
            .replace("d4", "d6")
            .replace("u4", "u6")
            .replace("s4", "s6")
            .replace("m4", "m6")
            .replace("c4", "c6")
            .replace("e4", "e6")
            .replace("t4", "t6");
        let text = format!("mcnp v d t 1 10 5\nc\ntally 1\n4\n{body6}kcode 1 0 5\n1 1 1 1 1\n");
        assert!(matches!(
            Mctal::parse(&text),
            Err(Error::BadStructure(m)) if m.contains("not in the declared tally list")
        ));
    }

    #[test]
    fn body_duplicate_number_errors() {
        let text = format!(
            "mcnp v d t 1 10 5\nc\ntally 1\n4\n{MIN_BODY}{MIN_BODY}kcode 1 0 5\n1 1 1 1 1\n"
        );
        assert!(matches!(
            Mctal::parse(&text),
            Err(Error::BadStructure(m)) if m.contains("duplicate tally body")
        ));
    }

    #[test]
    fn body_npert_third_token_stored_verbatim() {
        // Perturbation-style third token: stored, standard bodies still parse.
        let text =
            format!("mcnp v d t 1 10 5\nc\ntally 1 1\n4\n{MIN_BODY}kcode 1 0 5\n1 1 1 1 1\n");
        let m = Mctal::parse(&text).unwrap();
        assert_eq!(m.npert, Some("1".to_string()));
        assert_eq!(m.tallies.len(), 1);
    }
}
