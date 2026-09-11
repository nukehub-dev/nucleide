//! ENDL table reader for the EEDL/EPDL scope (`pyne.endl.Library`).
//!
//! Contract verified against `pyne/endl.py` (EEDL/EPDL scope only —
//! `endl.py:1-15`; the spec PDFs it cites are UNVERIFIED and not used here):
//! - Table framing: a table ends at the first line matching
//!   `END_OF_TABLE_RE = re.compile(" {71}1")`, i.e. 71 spaces then `1`
//!   anchored at the line start ([`is_end_of_table`]).
//! - Header line 1 slices: `[0:6]` zzzaaa nuclide, `[7:9]` incident particle
//!   `yi`, `[10:12]` outgoing particle `yo`, `[13:24]` atomic weight via
//!   [`endftod`], `[25:31]` date, `[31]` iflag.
//! - Header line 2 slices: `[0:2]` reaction descriptor, `[2:5]` reaction
//!   property, `[5:8]` modifier, `[21:32]` subshell indicator `x1` via
//!   [`endftod`]; `x1` is stored as `None` when `rmod == 0`.
//! - Body widths: `NFIELDS_RPROP = {0: 2, 10: 2, 11: 2, 21: 3, 22: 3}`;
//!   every field is a fixed 11-character ENDL number (`fromendl_tok` in
//!   `pyne/_utils.pyx:176-215`).
//! - Lookup: `get_rx(nuc, p_in, rdesc, rprop, x1=None, p_out=None)` returns
//!   the first table matching all given selectors (`endl.py:235-259`).
//!
//! Deviations from upstream (documented, not silent):
//! - The library parses eagerly into [`Table`]s; upstream stores byte limits
//!   and re-reads lazily. Results are identical for well-formed files.
//! - `endftod` mirrors the C++ implementation (`src/utils.cpp:85`: strip
//!   whitespace, empty → 0.0, `D` exponent → `E`, implicit `1.23-45` exponent
//!   via the last interior sign, else plain parse). C `strtod` also accepts a
//!   bare prefix of a malformed field (`"1.23E"` → 1.23); Rust's total parse
//!   cannot, so over-long/exponent-broken fields fall back to the longest
//!   parseable prefix. Synthetic tests pin the documented cases; exotic
//!   malformed fields carry deviation risk.
//! - Nuclide keys are the integer ids (`zzzaaa * 10000`, verified against
//!   `src/nucname.cpp:816`). Name resolution is caller-side: pass an integer
//!   id (e.g. `820000000` for natural Pb) or a fully-specified isotope name.
//! - Body lines are joined with `\n` before applying the `fromendl_tok`
//!   byte arithmetic, so files must use `\n` line endings.

use std::fmt;
use std::path::Path;

/// Reaction property → number of 11-char fields per body line
/// (`NFIELDS_RPROP` in `endl.py:38`).
pub fn fields_for_rprop(rprop: i32) -> Option<usize> {
    match rprop {
        0 | 10 | 11 => Some(2),
        21 | 22 => Some(3),
        _ => None,
    }
}

/// Errors raised while reading ENDL tables.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// I/O failure while reading the file.
    Io(String),
    /// A header line was missing or its slices failed to parse.
    BadHeader {
        /// 1-based line number of the second header line.
        line_no: usize,
        /// Why the header slices failed to parse.
        reason: String,
    },
    /// Reaction property outside `NFIELDS_RPROP`.
    UnknownRprop(i32),
    /// No table stored for this nucleus id.
    UnknownNucleus(i64),
    /// Nucleus known, but no table matched the selectors.
    NoMatch {
        /// Requested nucleus id (`zzzaaa * 10000`).
        nuc: i64,
        /// Requested incident-particle designator.
        p_in: i32,
        /// Requested reaction descriptor.
        rdesc: i32,
        /// Requested reaction property.
        rprop: i32,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadHeader { line_no, reason } => {
                write!(f, "malformed ENDL header near line {line_no}: {reason}")
            }
            Error::UnknownRprop(r) => write!(f, "unknown ENDL reaction property {r}"),
            Error::UnknownNucleus(n) => write!(f, "nucleus {n} does not exist"),
            Error::NoMatch {
                nuc,
                p_in,
                rdesc,
                rprop,
            } => write!(
                f,
                "no table for nucleus {nuc}, incident {p_in}, descriptor {rdesc}, property {rprop}"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// Convert one ENDL number field to f64 (`pyne::endftod`, `src/utils.cpp:85`).
///
/// The caller passes the exact field slice (upstream callers always hand over
/// exactly the 11-character entry). Whitespace is stripped; an empty field is
/// 0.0; an explicit `E/e/D/d` exponent is honored (`D` → `E`); otherwise the
/// last interior `+`/`-` introduces an implicit exponent (`1.23-45` → 1.23e-45).
pub fn endftod(field: &str) -> f64 {
    let compact: String = field.chars().filter(|c| !c.is_whitespace()).collect();
    if compact.is_empty() {
        return 0.0;
    }
    if let Some(pos) = compact.find(['e', 'E', 'd', 'D']) {
        let mut s = compact.clone();
        if s.as_bytes()[pos] == b'd' || s.as_bytes()[pos] == b'D' {
            s.replace_range(pos..=pos, "E");
        }
        if let Ok(v) = s.parse::<f64>() {
            return v;
        }
        return parse_prefix_float(&s);
    }
    // Implicit exponent: last sign that is not the leading character.
    let mut sign_pos: Option<usize> = None;
    for (i, c) in compact.char_indices() {
        if i == 0 {
            continue;
        }
        if c == '+' || c == '-' {
            sign_pos = Some(i);
        }
    }
    if let Some(pos) = sign_pos {
        let mut s = compact.clone();
        s.insert(pos, 'E');
        if let Ok(v) = s.parse::<f64>() {
            return v;
        }
        return parse_prefix_float(&s);
    }
    compact.parse::<f64>().unwrap_or(0.0)
}

/// Longest-parseable-prefix float, mimicking C `strtod` stopping at the first
/// character that cannot extend the number (`"1.23E"` → 1.23).
fn parse_prefix_float(s: &str) -> f64 {
    let mut end = s.len();
    while end > 0 {
        // Step back over one char (fields are ASCII; char-boundary safe).
        end = s[..end]
            .char_indices()
            .next_back()
            .map(|(i, _)| i)
            .unwrap_or(0);
        if let Ok(v) = s[..end].parse::<f64>() {
            return v;
        }
    }
    0.0
}

/// True when `line` ends an ENDL table: `END_OF_TABLE_RE = " {71}1"`
/// matched at the line start (`endl.py:34`).
pub fn is_end_of_table(line: &str) -> bool {
    let b = line.as_bytes();
    b.len() >= 72 && b[..71].iter().all(|c| *c == b' ') && b[71] == b'1'
}

/// One parsed ENDL table: header selectors plus its data rows.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// Nucleus id (`zzzaaa * 10000`).
    pub nuc: i64,
    /// Incident-particle designator (`yi`).
    pub p_in: i32,
    /// Outgoing-particle designator (`yo`).
    pub yo: i32,
    /// Atomic weight from header line 1.
    pub aw: f64,
    /// Evaluation date string.
    pub date: String,
    /// Header flag.
    pub iflag: i32,
    /// Reaction descriptor.
    pub rdesc: i32,
    /// Reaction property.
    pub rprop: i32,
    /// Table modifier (0 ⇒ `x1` is `None`).
    pub rmod: i32,
    /// Atomic subshell indicator (`None` when `rmod == 0`).
    pub x1: Option<i32>,
    /// Data rows, each `fields_for_rprop(rprop)` wide.
    pub data: Vec<Vec<f64>>,
}

/// An ENDL evaluation file: every table parsed eagerly.
#[derive(Debug, Clone, Default)]
pub struct Library {
    tables: Vec<Table>,
}

impl Library {
    /// Read and parse an ENDL file from disk.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Self::parse(&text)
    }

    /// Parse ENDL text in memory (`Library(fh)` header pass, `endl.py:67-129`).
    pub fn parse(text: &str) -> Result<Self, Error> {
        // Work on raw lines (no `\r` stripping: `\n` endings are the contract).
        let lines: Vec<&str> = text.lines().collect();
        let mut tables = Vec::new();
        let mut i = 0;
        while i < lines.len() {
            let line1 = lines[i];
            i += 1;
            if i >= lines.len() {
                // Upstream breaks when the second header line is missing.
                break;
            }
            let line2 = lines[i];
            i += 1;
            let line_no = i - 1;

            let table_no = tables.len();
            let header = parse_header(line1, line2, line_no)
                .map_err(|reason| Error::BadHeader { line_no, reason })?;

            // Skip to the end of the table, collecting body lines.
            let mut body: Vec<&str> = Vec::new();
            while i < lines.len() && !is_end_of_table(lines[i]) {
                body.push(lines[i]);
                i += 1;
            }
            // Consume the terminator line itself when present (EOF also ends
            // a table, mirroring `read_eot = len(line) == 0 or ...`).
            if i < lines.len() {
                i += 1;
            }

            let num_fields =
                fields_for_rprop(header.rprop).ok_or(Error::UnknownRprop(header.rprop))?;
            let data = parse_body(&body, num_fields);
            tables.push(Table {
                nuc: header.nuc,
                p_in: header.p_in,
                yo: header.yo,
                aw: header.aw,
                date: header.date,
                iflag: header.iflag,
                rdesc: header.rdesc,
                rprop: header.rprop,
                rmod: header.rmod,
                x1: if header.rmod != 0 {
                    Some(header.x1)
                } else {
                    None
                },
                data,
            });
            let _ = table_no;
        }
        Ok(Library { tables })
    }

    /// Distinct nucleus ids in file order.
    pub fn nuclides(&self) -> Vec<i64> {
        let mut out = Vec::new();
        for t in &self.tables {
            if !out.contains(&t.nuc) {
                out.push(t.nuc);
            }
        }
        out
    }

    /// Distinct incident particles for one nucleus.
    pub fn pin(&self, nuc: i64) -> Vec<i32> {
        let mut out = Vec::new();
        for t in self.tables.iter().filter(|t| t.nuc == nuc) {
            if !out.contains(&t.p_in) {
                out.push(t.p_in);
            }
        }
        out
    }

    /// All parsed tables (header + data).
    pub fn tables(&self) -> &[Table] {
        &self.tables
    }

    /// Reaction data for one selector set (`get_rx`, `endl.py:235-259`).
    ///
    /// The first table matching nucleus, incident particle, descriptor,
    /// property — plus `x1`/`p_out` when given — wins, exactly like the
    /// upstream `next(...)` over `data_tuples`.
    pub fn get_rx(
        &self,
        nuc: i64,
        p_in: i32,
        rdesc: i32,
        rprop: i32,
        x1: Option<i32>,
        p_out: Option<i32>,
    ) -> Result<&[Vec<f64>], Error> {
        if !self.tables.iter().any(|t| t.nuc == nuc) {
            return Err(Error::UnknownNucleus(nuc));
        }
        self.tables
            .iter()
            .filter(|t| t.nuc == nuc)
            .find(|t| {
                t.p_in == p_in
                    && t.rdesc == rdesc
                    && t.rprop == rprop
                    && (x1.is_none() || t.x1 == x1)
                    && (p_out.is_none() || Some(t.yo) == p_out)
            })
            .map(|t| t.data.as_slice())
            .ok_or(Error::NoMatch {
                nuc,
                p_in,
                rdesc,
                rprop,
            })
    }
}

struct RawHeader {
    nuc: i64,
    p_in: i32,
    yo: i32,
    aw: f64,
    date: String,
    iflag: i32,
    rdesc: i32,
    rprop: i32,
    rmod: i32,
    x1: i32,
}

/// Byte-index a line slice, tolerating short lines (missing ⇒ empty).
fn slice(line: &str, from: usize, to: usize) -> &str {
    let b = line.as_bytes();
    if from >= b.len() {
        return "";
    }
    let end = to.min(b.len());
    if end <= from {
        return "";
    }
    // Header lines are ASCII; a non-boundary cut falls back to empty rather
    // than panicking.
    std::str::from_utf8(&b[from..end]).unwrap_or("")
}

fn parse_int_cell(cell: &str, fallback: i32) -> i32 {
    let t = cell.trim();
    if t.is_empty() {
        return fallback;
    }
    t.parse::<i32>().unwrap_or(fallback)
}

fn parse_header(line1: &str, line2: &str, line_no: usize) -> Result<RawHeader, String> {
    let zzzaaa = parse_int_cell(slice(line1, 0, 6), -1);
    if zzzaaa < 0 {
        return Err(format!(
            "line {line_no}: bad zzzaaa `{}`",
            slice(line1, 0, 6)
        ));
    }
    let aw_slice = slice(line1, 13, 24);
    let aw = if aw_slice.is_empty() {
        -1.0
    } else {
        endftod(aw_slice)
    };
    let x1_slice = slice(line2, 21, 32);
    let x1 = if x1_slice.trim().is_empty() {
        -1
    } else {
        endftod(x1_slice) as i32
    };
    Ok(RawHeader {
        // Verified against `src/nucname.cpp:816`: `zzzaaa_to_id(n) = n*10000`.
        nuc: i64::from(zzzaaa) * 10_000,
        p_in: parse_int_cell(slice(line1, 7, 9), -1),
        yo: parse_int_cell(slice(line1, 10, 12), -1),
        aw,
        date: slice(line1, 25, 31).trim().to_string(),
        iflag: parse_int_cell(slice(line1, 31, 32), 0),
        rdesc: parse_int_cell(slice(line2, 0, 2), -1),
        rprop: parse_int_cell(slice(line2, 2, 5), -1),
        rmod: parse_int_cell(slice(line2, 5, 8), -1),
        x1,
    })
}

/// Parse body lines via the `fromendl_tok` byte arithmetic
/// (`_utils.pyx:176-215`): `line_length = num_fields * 11 + 1`, field `j` of
/// line `i` starts at `j * 11 + i * line_length`; trailing partial bytes are
/// ignored, exactly like `num_lines = len // line_length`.
fn parse_body(lines: &[&str], num_fields: usize) -> Vec<Vec<f64>> {
    if num_fields == 0 {
        return Vec::new();
    }
    let joined = lines.join("\n");
    // fromendl_tok reads the raw file span, whose lines each carry their
    // newline; re-add it so the byte arithmetic lines up.
    let raw = if joined.is_empty() {
        String::new()
    } else {
        joined + "\n"
    };
    let b = raw.as_bytes();
    let line_length = num_fields * 11 + 1;
    let num_lines = b.len() / line_length;
    let mut out = Vec::with_capacity(num_lines);
    for i in 0..num_lines {
        let mut row = Vec::with_capacity(num_fields);
        for j in 0..num_fields {
            let start = j * 11 + i * line_length;
            let field = std::str::from_utf8(&b[start..start + 11]).unwrap_or("");
            row.push(endftod(field));
        }
        out.push(row);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(v: f64) -> String {
        // Render one 11-char implicit-exponent ENDL field (sign of exponent
        // doubles as the separator, exactly like the evaluated files).
        if v == 0.0 {
            return " 0.00000+00".to_string();
        }
        let exp = v.abs().log10().floor() as i32;
        let mant = v / 10f64.powi(exp);
        format!(" {:7.5}{:+03}", mant, exp)
    }

    fn table(
        nuc: &str,
        yi: &str,
        yo: &str,
        rdesc: &str,
        rprop: &str,
        rmod: &str,
        x1: f64,
    ) -> String {
        let line1 = format!("{nuc:<6} {yi:>2} {yo:>2} {} {}0", cell(207.2), "090101");
        let mut line2 = format!("{rdesc:>2}{rprop:>3}{rmod:>3}");
        line2.push_str(&" ".repeat(21 - line2.len()));
        line2.push_str(&cell(x1));
        format!("{line1}\n{line2}\n")
    }

    fn eot() -> String {
        format!("{}1\n", " ".repeat(71))
    }

    #[test]
    fn endftod_vectors_match_cpp_semantics() {
        assert_eq!(endftod(""), 0.0);
        assert_eq!(endftod("           "), 0.0);
        assert_eq!(endftod(" 2.0720+02 "), 207.2);
        assert_eq!(endftod(" 1.00000-05"), 1e-5);
        assert!((endftod(" 3.63530+05") - 363530.0).abs() < 1e-6);
        // Explicit E and D exponents.
        assert_eq!(endftod(" 1.00000E+05"), 1e5);
        assert_eq!(endftod(" 1.00000D+05"), 1e5);
        assert_eq!(endftod(" 1.00000e-05"), 1e-5);
        // Plain decimals and leading signs.
        assert_eq!(endftod(" 207.20000 "), 207.2);
        assert_eq!(endftod("-1.50000+00"), -1.5);
        // strtod-prefix behavior on a broken exponent.
        assert_eq!(endftod(" 1.23000E  "), 1.23);
    }

    #[test]
    fn end_of_table_needs_71_spaces_then_1() {
        assert!(is_end_of_table(&format!("{}1", " ".repeat(71))));
        assert!(is_end_of_table(&format!(
            "{}1 trailing junk",
            " ".repeat(71)
        )));
        assert!(!is_end_of_table(&format!("{}1", " ".repeat(70))));
        assert!(!is_end_of_table(&format!("{}1", " ".repeat(72))));
        assert!(!is_end_of_table(" 1.00000-05 3.63530+05"));
        assert!(!is_end_of_table(""));
    }

    #[test]
    fn fields_for_rprop_mirrors_nfields_rprop() {
        assert_eq!(fields_for_rprop(0), Some(2));
        assert_eq!(fields_for_rprop(10), Some(2));
        assert_eq!(fields_for_rprop(11), Some(2));
        assert_eq!(fields_for_rprop(21), Some(3));
        assert_eq!(fields_for_rprop(22), Some(3));
        assert_eq!(fields_for_rprop(7), None);
    }

    #[test]
    fn header_slices_and_rmod_zero_x1_none() {
        let text = table(" 82000", " 9", " 0", "10", "  0", "  0", 5.0)
            + &format!("{}\n", cell(1e-5) + &cell(363530.0))
            + &eot();
        let lib = Library::parse(&text).unwrap();
        assert_eq!(lib.nuclides(), vec![820000000]);
        let t = &lib.tables()[0];
        assert_eq!(t.p_in, 9);
        assert_eq!(t.yo, 0);
        assert!((t.aw - 207.2).abs() < 1e-9);
        assert_eq!(t.date, "090101");
        assert_eq!(t.iflag, 0);
        assert_eq!(t.rdesc, 10);
        assert_eq!(t.rprop, 0);
        assert_eq!(t.rmod, 0);
        // rmod == 0 forces x1 to None even though the field held 5.0.
        assert_eq!(t.x1, None);
        assert_eq!(t.data.len(), 1);
        assert!((t.data[0][0] - 1e-5).abs() / 1e-5 < 1e-9);
        assert!((t.data[0][1] - 363530.0).abs() / 363530.0 < 1e-9);
    }

    #[test]
    fn get_rx_selects_first_match_and_filters() {
        let mut text = table(" 82000", " 9", " 9", "81", "  0", "  1", 1.0)
            + &format!("{}\n", cell(8.829e-2) + &cell(0.566158))
            + &eot();
        text += &table(" 82000", " 9", "19", "81", "  0", "  1", 2.0);
        text += &format!("{}\n", cell(1e-3) + &cell(2.5));
        text += &eot();
        let lib = Library::parse(&text).unwrap();

        // No selector: first (x1=1) table wins.
        let d = lib.get_rx(820000000, 9, 81, 0, None, None).unwrap();
        assert!((d[0][0] - 8.829e-2).abs() / 8.829e-2 < 1e-9);
        // Subshell + outgoing-particle selectors discriminate.
        let d = lib.get_rx(820000000, 9, 81, 0, Some(2), None).unwrap();
        assert!((d[0][0] - 1e-3).abs() / 1e-3 < 1e-9);
        let d = lib.get_rx(820000000, 9, 81, 0, Some(1), Some(9)).unwrap();
        assert!((d[0][1] - 0.566158).abs() / 0.566158 < 1e-9);
        // Misses are typed.
        assert_eq!(
            lib.get_rx(820000000, 9, 81, 0, Some(99), None),
            Err(Error::NoMatch {
                nuc: 820000000,
                p_in: 9,
                rdesc: 81,
                rprop: 0,
            })
        );
        assert_eq!(
            lib.get_rx(820000000, 8, 81, 0, None, None),
            Err(Error::NoMatch {
                nuc: 820000000,
                p_in: 8,
                rdesc: 81,
                rprop: 0,
            })
        );
        assert_eq!(
            lib.get_rx(920000000, 9, 81, 0, None, None),
            Err(Error::UnknownNucleus(920000000))
        );
    }

    #[test]
    fn three_column_body_parses() {
        let text = table(" 82000", " 9", " 0", "82", " 21", "  0", 0.0)
            + &format!("{}\n", cell(1e-5) + &cell(1e-7) + &cell(4.6132e-8))
            + &format!("{}\n", cell(1e5) + &cell(1e5) + &cell(9.15055e6))
            + &eot();
        let lib = Library::parse(&text).unwrap();
        let d = lib.get_rx(820000000, 9, 82, 21, None, None).unwrap();
        assert_eq!(d.len(), 2);
        assert_eq!(d[0].len(), 3);
        assert!((d[1][2] - 9.15055e6).abs() / 9.15055e6 < 1e-9);
    }

    #[test]
    fn unknown_rprop_errors() {
        let text = table(" 82000", " 9", " 0", "10", " 99", "  0", 0.0) + &eot();
        assert_eq!(Library::parse(&text).unwrap_err(), Error::UnknownRprop(99));
    }

    #[test]
    fn missing_second_header_line_stops() {
        let lib = Library::parse(" 82000  9   0  2.0720+02  0901010\n").unwrap();
        assert!(lib.tables().is_empty());
    }

    #[test]
    fn prefix_float_and_open_success() {
        // Broken implicit exponents fall back to the longest prefix.
        assert_eq!(endftod(" 1.20000+03X"), 1200.0);
        assert_eq!(endftod(" E+3X"), 0.0);
        // Short-line slicing tolerates missing cells.
        assert_eq!(slice("ab", 5, 9), "");
        assert_eq!(slice("ab", 1, 1), "");
        assert_eq!(slice("abcdef", 1, 4), "bcd");
        assert!(parse_body(&[], 0).is_empty());
        // A real file opens through the same parse path.
        let lib = Library::open(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/endl/synthetic_eedl.txt"
        ))
        .unwrap();
        assert!(!lib.tables().is_empty());
    }

    #[test]
    fn short_and_bad_headers() {
        // Truncated line 1: empty atomic-weight and subshell cells.
        let text = " 82000\n        \n".to_string() + &eot();
        assert_eq!(Library::parse(&text).unwrap_err(), Error::UnknownRprop(-1));
        // Non-numeric nucleus id.
        let text = "XXXXXX\n        \n".to_string() + &eot();
        assert!(matches!(
            Library::parse(&text).unwrap_err(),
            Error::BadHeader { .. }
        ));
    }

    #[test]
    fn errors_display_and_helpers() {
        assert!(format!("{}", Error::Io("x".into())).contains("io error"));
        assert!(format!(
            "{}",
            Error::BadHeader {
                line_no: 3,
                reason: "r".into()
            }
        )
        .contains("line 3"));
        assert!(format!("{}", Error::UnknownRprop(7)).contains('7'));
        assert!(format!("{}", Error::UnknownNucleus(1)).contains('1'));
        assert!(format!(
            "{}",
            Error::NoMatch {
                nuc: 1,
                p_in: 2,
                rdesc: 3,
                rprop: 4,
            }
        )
        .contains('1'));
        assert!(matches!(
            Library::open("/nonexistent-endl-file"),
            Err(Error::Io(_))
        ));
        // Empty body parses to zero rows.
        let text = table(" 82000", " 9", " 0", "10", "  0", "  0", 0.0) + &eot();
        let lib = Library::parse(&text).unwrap();
        assert_eq!(lib.pin(820000000), vec![9]);
        assert!(lib.pin(1).is_empty());
        assert!(lib
            .get_rx(820000000, 9, 10, 0, None, None)
            .unwrap()
            .is_empty());
    }
}
