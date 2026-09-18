//! `&INDATA` namelist reader (VMEC-input text grammar).
//!
//! Parses the Fortran namelist block (`&INDATA ... /`) that carries a
//! VMEC equilibrium input: scalar switches (`NFP`, `NCURR`, `MPOL`,
//! `NTOR`, `LFREEB`, ...), power-series profile coefficients (`AM`,
//! `AI`, `AC`, ...), axis and boundary Fourier tables (`RAXIS`,
//! `ZAXIS`, `RBC`/`RBS`/`ZBC`/`ZBS`), and the profile-type / multistage
//! switches (`PMASS_TYPE`, `PCURR_TYPE`, `PIOTA_TYPE`, `NS_ARRAY`,
//! `FTOL_ARRAY`, `NITER_ARRAY`, ...). Names are case-insensitive and
//! stored upper-cased; values keep Fortran number spellings (`1.0D-3`),
//! repeat notation (`3*0.0`), and `.TRUE.`/bare-`T` logicals. Indexed
//! assignment takes scalar indices only — sliced forms like
//! `RBC(0:4,2)=...` are a loud [`Error::SlicedIndex`], per the
//! fixed-boundary conversion stance. Reads data, never solves
//! equilibria.

use std::collections::BTreeMap;

use crate::error::{Error, Result};

/// One namelist value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Integer literal.
    Int(i64),
    /// Real literal (Fortran `D`/`Q` exponents read as `E`).
    Float(f64),
    /// `.TRUE.`/`.FALSE.` (or bare `T`/`F` in value position).
    Bool(bool),
    /// Single- or double-quoted string (original case kept).
    Str(String),
}

impl Value {
    /// Coerce to `i64` (floats accepted only when integral).
    pub fn as_int(&self) -> Result<i64> {
        match self {
            Value::Int(v) => Ok(*v),
            Value::Float(v) if v.is_finite() && v.fract() == 0.0 => Ok(*v as i64),
            _ => Err(Error::BadValue(format!("{self:?} is not an integer"))),
        }
    }

    /// Coerce to `f64` (integers widen exactly).
    pub fn as_float(&self) -> Result<f64> {
        match self {
            Value::Int(v) => Ok(*v as f64),
            Value::Float(v) if v.is_finite() => Ok(*v),
            _ => Err(Error::BadValue(format!("{self:?} is not a float"))),
        }
    }

    /// Coerce to `bool`.
    pub fn as_bool(&self) -> Result<bool> {
        match self {
            Value::Bool(v) => Ok(*v),
            _ => Err(Error::BadValue(format!("{self:?} is not a logical"))),
        }
    }

    /// Borrow the string payload.
    pub fn as_str(&self) -> Result<&str> {
        match self {
            Value::Str(s) => Ok(s),
            _ => Err(Error::BadValue(format!("{self:?} is not a string"))),
        }
    }
}

/// Parsed `&INDATA` block: scalar assignments plus indexed assignments,
/// both keyed by upper-cased variable name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Indata {
    /// Scalar assignments (`NFP = 5`); repeats keep the last value.
    pub scalars: BTreeMap<String, Value>,
    /// Indexed assignments (`AM(3) = 1.0`); repeats overwrite per index.
    /// Whole-array assignment (`AM = 1.0 2.0 ...`) fills from index 0.
    pub indexed: BTreeMap<String, BTreeMap<Vec<i64>, Value>>,
}

impl Indata {
    /// Parse `&INDATA ... /` text.
    pub fn parse(text: &str) -> Result<Indata> {
        parse_indata(text)
    }

    /// Borrow a scalar value by name (case-insensitive).
    pub fn get(&self, name: &str) -> Option<&Value> {
        self.scalars.get(&name.to_ascii_uppercase())
    }

    /// Borrow indexed entries by name (case-insensitive), sorted by index.
    pub fn indexed_entries(&self, name: &str) -> Vec<(&Vec<i64>, &Value)> {
        match self.indexed.get(&name.to_ascii_uppercase()) {
            None => Vec::new(),
            Some(map) => map.iter().collect(),
        }
    }

    fn typed(&self, name: &str, expected: &'static str) -> Result<&Value> {
        self.get(name)
            .ok_or_else(|| Error::MissingVariable(name.to_string()))
            .and_then(|v| {
                let ok = match expected {
                    "an integer" => matches!(v, Value::Int(_) | Value::Float(_)),
                    "a float" => matches!(v, Value::Int(_) | Value::Float(_)),
                    "a logical" => matches!(v, Value::Bool(_)),
                    "a string" => matches!(v, Value::Str(_)),
                    _ => false,
                };
                if ok {
                    Ok(v)
                } else {
                    Err(Error::MismatchedType {
                        name: name.to_ascii_uppercase(),
                        expected,
                    })
                }
            })
    }

    /// Integer switch by name (`NFP`, `NCURR`, `MPOL`, `NTOR`, ...).
    pub fn int(&self, name: &str) -> Result<i64> {
        self.typed(name, "an integer")?.as_int()
    }

    /// Real value by name.
    pub fn float(&self, name: &str) -> Result<f64> {
        self.typed(name, "a float")?.as_float()
    }

    /// Logical switch by name (`LFREEB`, `LASYM`, ...).
    pub fn boolean(&self, name: &str) -> Result<bool> {
        self.typed(name, "a logical")?.as_bool()
    }

    /// String value by name (`MGRID_FILE`, `*_TYPE`, ...).
    pub fn string(&self, name: &str) -> Result<&str> {
        self.typed(name, "a string")?.as_str()
    }

    /// One-dimensional float series (`AM`, `AI`, `AC`, `RAXIS`, ...):
    /// `(index, value)` pairs sorted by index.
    pub fn float_series_1d(&self, name: &str) -> Result<Vec<(i64, f64)>> {
        let map = self
            .indexed
            .get(&name.to_ascii_uppercase())
            .ok_or_else(|| Error::MissingVariable(name.to_string()))?;
        map.iter()
            .map(|(idx, v)| {
                if idx.len() != 1 {
                    return Err(Error::BadValue(format!(
                        "{name} index has {} entries, expected 1",
                        idx.len()
                    )));
                }
                Ok((idx[0], v.as_float()?))
            })
            .collect()
    }

    /// Two-dimensional float table (`RBC`, `RBS`, `ZBC`, `ZBS`):
    /// `((m, n), value)` entries sorted by index pair.
    pub fn float_table_2d(&self, name: &str) -> Result<Vec<((i64, i64), f64)>> {
        let map = self
            .indexed
            .get(&name.to_ascii_uppercase())
            .ok_or_else(|| Error::MissingVariable(name.to_string()))?;
        map.iter()
            .map(|(idx, v)| {
                if idx.len() != 2 {
                    return Err(Error::BadValue(format!(
                        "{name} index has {} entries, expected 2",
                        idx.len()
                    )));
                }
                Ok(((idx[0], idx[1]), v.as_float()?))
            })
            .collect()
    }

    /// `true` when `LFREEB` is set (free-boundary input). Absent means
    /// fixed boundary — the only conversion stance this reader serves.
    pub fn is_free_boundary(&self) -> Result<bool> {
        match self.get("LFREEB") {
            None => Ok(false),
            Some(v) => v.as_bool(),
        }
    }

    /// Raw `NCURR` switch, when present.
    pub fn ncurr(&self) -> Result<Option<i64>> {
        match self.get("NCURR") {
            None => Ok(None),
            Some(v) => v.as_int().map(Some),
        }
    }

    /// Whether `NCURR = 1`, i.e. the current profile is the `I'(s)`
    /// derivative profile (the documented conversion assumption).
    /// `None` when `NCURR` is absent.
    pub fn ncurr_is_iprime(&self) -> Result<Option<bool>> {
        Ok(self.ncurr()?.map(|n| n == 1))
    }

    /// Profile-type switch (`PMASS_TYPE`, `PCURR_TYPE`, `PIOTA_TYPE`),
    /// when present. The conversion stance covers power-series profiles.
    pub fn profile_type(&self, which: &str) -> Result<Option<String>> {
        match self.get(which) {
            None => Ok(None),
            Some(v) => v.as_str().map(|s| Some(s.to_string())),
        }
    }
}

/// Parse `&INDATA ... /` text (see [`Indata::parse`]).
pub fn parse_indata(text: &str) -> Result<Indata> {
    let code = strip_comments(text)?;
    let block = indata_block(&code)?;
    let tokens = lex(block)?;
    parse_block(&tokens)
}

/// Strip `!` comments (outside quotes) from the text.
fn strip_comments(text: &str) -> Result<String> {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars().peekable();
    let mut quote: Option<char> = None;
    while let Some(c) = chars.next() {
        if let Some(q) = quote {
            out.push(c);
            if c == q {
                // Fortran `''` escape inside single-quoted strings.
                if q == '\'' && chars.peek() == Some(&'\'') {
                    out.push(chars.next().expect("peeked"));
                } else {
                    quote = None;
                }
            }
            continue;
        }
        match c {
            '\'' | '"' => {
                quote = Some(c);
                out.push(c);
            }
            '!' => {
                for c2 in chars.by_ref() {
                    if c2 == '\n' {
                        out.push('\n');
                        break;
                    }
                }
            }
            _ => out.push(c),
        }
    }
    if quote.is_some() {
        return Err(Error::BadNamelist("unterminated string".into()));
    }
    Ok(out)
}

/// Cut out the `&INDATA ... /` block (outside quotes).
fn indata_block(code: &str) -> Result<&str> {
    let bytes = code.as_bytes();
    let mut i = 0;
    let mut quote: Option<u8> = None;
    let mut start: Option<usize> = None;
    while i < bytes.len() {
        let c = bytes[i];
        if let Some(q) = quote {
            if c == q {
                if q == b'\'' && i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                    i += 2;
                    continue;
                }
                quote = None;
            }
            i += 1;
            continue;
        }
        match c {
            b'\'' | b'"' => {
                quote = Some(c);
                i += 1;
            }
            b'&' => {
                let rest = &code[i + 1..];
                let name: String = rest
                    .chars()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                    .collect();
                if name.eq_ignore_ascii_case("INDATA") {
                    start = Some(i + 1 + name.len());
                    i += 1 + name.len();
                } else if start.is_none() {
                    // Some other namelist block: skip its `&name`.
                    i += 1 + name.len();
                } else {
                    return Err(Error::BadNamelist("unexpected `&` inside &INDATA".into()));
                }
            }
            b'/' => {
                if let Some(s) = start {
                    return Ok(&code[s..i]);
                }
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }
    if start.is_some() {
        return Err(Error::UnterminatedIndata);
    }
    Err(Error::MissingIndata)
}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Name(String),
    Integer(i64),
    Real(f64),
    Logical(bool),
    Quoted(String),
    Eq,
    LParen,
    RParen,
    Colon,
    Star,
    Comma,
}

fn lex(block: &str) -> Result<Vec<Tok>> {
    let bytes = block.as_bytes();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        match c {
            b' ' | b'\t' | b'\r' | b'\n' => {
                i += 1;
            }
            b'=' => {
                toks.push(Tok::Eq);
                i += 1;
            }
            b'(' => {
                toks.push(Tok::LParen);
                i += 1;
            }
            b')' => {
                toks.push(Tok::RParen);
                i += 1;
            }
            b':' => {
                toks.push(Tok::Colon);
                i += 1;
            }
            b'*' => {
                toks.push(Tok::Star);
                i += 1;
            }
            b',' => {
                toks.push(Tok::Comma);
                i += 1;
            }
            b'\'' | b'"' => {
                let q = c;
                i += 1;
                let mut s = String::new();
                loop {
                    if i >= bytes.len() {
                        return Err(Error::BadNamelist("unterminated string".into()));
                    }
                    if bytes[i] == q {
                        if q == b'\'' && i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                            s.push('\'');
                            i += 2;
                            continue;
                        }
                        i += 1;
                        break;
                    }
                    s.push(bytes[i] as char);
                    i += 1;
                }
                toks.push(Tok::Quoted(s));
            }
            b'.' => {
                // Dotted logical (`.TRUE.`) or a leading-dot number (`.5`).
                if i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit() {
                    let (v, j) = lex_number(bytes, i, false)?;
                    toks.push(v);
                    i = j;
                    continue;
                }
                // Read `.WORD.` explicitly.
                let mut j = i + 1;
                let mut word = String::new();
                while j < bytes.len() && bytes[j] != b'.' {
                    word.push(bytes[j] as char);
                    j += 1;
                }
                if j >= bytes.len() {
                    return Err(Error::BadNamelist("unterminated dotted logical".into()));
                }
                j += 1; // closing dot
                match word.to_ascii_uppercase().as_str() {
                    "T" | "TRUE" => toks.push(Tok::Logical(true)),
                    "F" | "FALSE" => toks.push(Tok::Logical(false)),
                    _ => {
                        return Err(Error::BadNamelist(format!(
                            "unknown dotted word `.{word}.`"
                        )));
                    }
                }
                i = j;
            }
            _ if c.is_ascii_alphabetic() || c == b'_' => {
                let mut j = i;
                while j < bytes.len() && (bytes[j].is_ascii_alphanumeric() || bytes[j] == b'_') {
                    j += 1;
                }
                let word = &block[i..j];
                // Bare T/F in value position are logicals; anything longer
                // (or followed by `(`/`=`) is a name. Single letters that
                // are not T/F are names too.
                let upper = word.to_ascii_uppercase();
                if (upper == "T" || upper == "F") && is_value_position(&toks) {
                    toks.push(Tok::Logical(upper == "T"));
                } else {
                    toks.push(Tok::Name(word.to_string()));
                }
                i = j;
            }
            _ if c.is_ascii_digit() || c == b'+' || c == b'-' => {
                let (v, j) = lex_number(bytes, i, true)?;
                toks.push(v);
                i = j;
            }
            _ => {
                return Err(Error::BadNamelist(format!(
                    "unexpected character `{}`",
                    c as char
                )));
            }
        }
    }
    Ok(toks)
}

/// Whether the next bare word sits in value position (right after `=`,
/// `,`, `(`, `*`, or another value).
fn is_value_position(toks: &[Tok]) -> bool {
    matches!(
        toks.last(),
        Some(Tok::Eq | Tok::Comma | Tok::LParen | Tok::Star)
            | Some(Tok::Integer(_) | Tok::Real(_) | Tok::Logical(_) | Tok::Quoted(_))
    )
}

/// Lex a number starting at `i`; `signed` allows a leading sign.
fn lex_number(bytes: &[u8], i: usize, signed: bool) -> Result<(Tok, usize)> {
    let mut j = i;
    if signed && j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
        j += 1;
    }
    let int_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    let has_int = j > int_start;
    let mut has_frac = false;
    if j < bytes.len() && bytes[j] == b'.' {
        // A trailing dot is a fraction only when followed by a digit,
        // an exponent letter, or end-of-number punctuation; `3.` reads
        // as 3.0 while `3./` still terminates cleanly.
        has_frac = true;
        j += 1;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
    }
    let mut has_exp = false;
    if j < bytes.len() && matches!(bytes[j], b'E' | b'e' | b'D' | b'd' | b'Q' | b'q') {
        has_exp = true;
        j += 1;
        if j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j == exp_start {
            return Err(Error::BadNamelist("malformed exponent".into()));
        }
    }
    if !has_int && !has_frac {
        return Err(Error::BadNamelist("malformed number".into()));
    }
    let raw = std::str::from_utf8(&bytes[i..j])
        .map_err(|_| Error::BadNamelist("non-UTF8 number".into()))?;
    if !has_frac && !has_exp {
        if let Ok(v) = raw.parse::<i64>() {
            return Ok((Tok::Integer(v), j));
        }
    }
    let normalized = raw.replace(['D', 'd', 'Q', 'q'], "E");
    normalized
        .parse::<f64>()
        .map(|v| (Tok::Real(v), j))
        .map_err(|_| Error::BadNamelist(format!("malformed number `{raw}`")))
}

struct Parser<'a> {
    toks: &'a [Tok],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a Tok> {
        self.toks.get(self.pos)
    }

    fn peek2(&self) -> Option<&'a Tok> {
        self.toks.get(self.pos + 1)
    }

    fn next(&mut self) -> Option<&'a Tok> {
        let t = self.toks.get(self.pos)?;
        self.pos += 1;
        Some(t)
    }

    fn skip_commas(&mut self) {
        while self.peek() == Some(&Tok::Comma) {
            self.pos += 1;
        }
    }
}

fn parse_block(toks: &[Tok]) -> Result<Indata> {
    let mut out = Indata::default();
    let mut p = Parser { toks, pos: 0 };
    while p.pos < toks.len() {
        p.skip_commas();
        if p.pos >= toks.len() {
            break;
        }
        let name = match p.next() {
            Some(Tok::Name(n)) => n.to_ascii_uppercase(),
            other => {
                return Err(Error::BadNamelist(format!(
                    "expected a variable name, found {other:?}"
                )));
            }
        };
        // Optional index list.
        let mut indices: Option<Vec<i64>> = None;
        if p.peek() == Some(&Tok::LParen) {
            p.pos += 1;
            let mut idx = Vec::new();
            loop {
                match p.next() {
                    Some(Tok::Integer(v)) => idx.push(*v),
                    Some(Tok::Colon) => return Err(Error::SlicedIndex { name }),
                    other => {
                        return Err(Error::BadNamelist(format!(
                            "bad index in `{name}`, found {other:?}"
                        )));
                    }
                }
                match p.next() {
                    Some(Tok::Comma) => continue,
                    Some(Tok::RParen) => break,
                    Some(Tok::Colon) => return Err(Error::SlicedIndex { name }),
                    other => {
                        return Err(Error::BadNamelist(format!(
                            "bad index list in `{name}`, found {other:?}"
                        )));
                    }
                }
            }
            indices = Some(idx);
        }
        if p.next() != Some(&Tok::Eq) {
            return Err(Error::BadNamelist(format!("`{name}` misses `=`")));
        }
        // Value list: one or more values (repeat `n*v` included) until
        // the next `Name =` / `Name(` assignment or end of block.
        let mut values: Vec<Value> = Vec::new();
        loop {
            p.skip_commas();
            match p.peek() {
                None => break,
                Some(Tok::Name(_)) => {
                    let is_next_assignment = matches!(
                        (p.peek(), p.peek2()),
                        (Some(Tok::Name(_)), Some(Tok::Eq | Tok::LParen))
                    );
                    if is_next_assignment {
                        break;
                    }
                    return Err(Error::BadNamelist(format!(
                        "stray name inside `{name}` values"
                    )));
                }
                _ => {}
            }
            // Repeat prefix `n*value`.
            let mut repeat = 1usize;
            if let Some(Tok::Integer(n)) = p.peek() {
                let n = *n;
                if p.peek2() == Some(&Tok::Star) {
                    if n < 0 {
                        return Err(Error::BadNamelist("negative repeat count".into()));
                    }
                    repeat = n as usize;
                    p.pos += 2;
                }
            }
            let value = match p.next() {
                Some(Tok::Integer(v)) => Value::Int(*v),
                Some(Tok::Real(v)) => Value::Float(*v),
                Some(Tok::Logical(v)) => Value::Bool(*v),
                Some(Tok::Quoted(s)) => Value::Str(s.clone()),
                // Bare T/F lexed as names when not in value position
                // (e.g. right after another value): recover them here.
                Some(Tok::Name(w)) if w.eq_ignore_ascii_case("t") => Value::Bool(true),
                Some(Tok::Name(w)) if w.eq_ignore_ascii_case("f") => Value::Bool(false),
                other => {
                    return Err(Error::BadNamelist(format!(
                        "bad value for `{name}`, found {other:?}"
                    )));
                }
            };
            for _ in 0..repeat {
                values.push(value.clone());
            }
        }
        if values.is_empty() {
            return Err(Error::BadNamelist(format!("`{name}` has no values")));
        }
        match indices {
            None if values.len() == 1 => {
                out.scalars.insert(name, values.pop().expect("one"));
            }
            None => {
                // Whole-array assignment fills from index 0.
                let entry = out.indexed.entry(name).or_default();
                for (k, v) in values.into_iter().enumerate() {
                    entry.insert(vec![k as i64], v);
                }
            }
            Some(idx) if values.len() == 1 => {
                out.indexed
                    .entry(name)
                    .or_default()
                    .insert(idx, values.pop().expect("one"));
            }
            Some(_) => {
                return Err(Error::BadNamelist(format!(
                    "`{name}`: several values need no index list"
                )));
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r"
! A stellarator input in the documented VMEC-input text grammar.
&INDATA
  LFREEB = F
  MGRID_FILE = 'none'
  DELT = 0.9
  NFP = 3
  NCURR = 1
  MPOL = 5 NTOR = 4
  PMASS_TYPE = 'power_series'
  PCURR_TYPE = 'power_series'
  AM = 0.0 1.0D0 3*0.0
  AI(2) = 1.5E-3
  RAXIS = 1.0 0.1
  RBC(0,0) = 1.0
  RBC(1,0) = -2.5D-1
  ZBS(1,0) = 2.5e-1
/
";

    #[test]
    fn parses_stellopt_style_block() {
        let nim = Indata::parse(SAMPLE).unwrap();
        assert_eq!(nim.int("nfp"), Ok(3));
        assert_eq!(nim.int("NFP"), Ok(3));
        assert_eq!(nim.int("mpol"), Ok(5));
        assert_eq!(nim.int("ntor"), Ok(4));
        assert!(!nim.boolean("lfreeb").unwrap());
        assert!(!nim.is_free_boundary().unwrap());
        assert_eq!(nim.float("delt"), Ok(0.9));
        assert_eq!(nim.string("mgrid_file"), Ok("none"));
        assert_eq!(nim.ncurr(), Ok(Some(1)));
        assert_eq!(nim.ncurr_is_iprime(), Ok(Some(true)));
        assert_eq!(
            nim.profile_type("pmass_type").unwrap(),
            Some("power_series".to_string())
        );
        // Whole-array fill from index 0, with D-exponent and repeat runs.
        let am = nim.float_series_1d("am").unwrap();
        assert_eq!(am, vec![(0, 0.0), (1, 1.0), (2, 0.0), (3, 0.0), (4, 0.0)]);
        assert_eq!(nim.float_series_1d("ai").unwrap(), vec![(2, 1.5e-3)]);
        assert_eq!(
            nim.float_table_2d("rbc").unwrap(),
            vec![((0, 0), 1.0), ((1, 0), -0.25)]
        );
        assert_eq!(nim.float_table_2d("zbs").unwrap(), vec![((1, 0), 0.25)]);
    }

    #[test]
    fn rejects_sliced_indices_loudly() {
        let text = "&INDATA\nRBC(0:4,2) = 1.0 2.0 3.0 4.0 5.0\n/\n";
        assert_eq!(
            Indata::parse(text),
            Err(Error::SlicedIndex {
                name: "RBC".to_string()
            })
        );
    }

    #[test]
    fn grammar_edge_vectors() {
        // Missing block.
        assert_eq!(Indata::parse("NFP = 3\n"), Err(Error::MissingIndata));
        // Unterminated block.
        assert_eq!(
            Indata::parse("&INDATA\nNFP = 3\n"),
            Err(Error::UnterminatedIndata)
        );
        // Bare T logical and dotted logicals.
        let nim = Indata::parse("&INDATA\nA = T\nB = .TRUE.\nC = .false.\n/\n").unwrap();
        assert_eq!(nim.boolean("a"), Ok(true));
        assert_eq!(nim.boolean("b"), Ok(true));
        assert_eq!(nim.boolean("c"), Ok(false));
        // Type mismatch is loud and named.
        assert_eq!(
            nim.int("a"),
            Err(Error::MismatchedType {
                name: "A".to_string(),
                expected: "an integer"
            })
        );
        // Missing scalar is loud and named.
        assert_eq!(
            nim.int("zzz"),
            Err(Error::MissingVariable("zzz".to_string()))
        );
        // Free boundary reads through.
        let nim = Indata::parse("&INDATA\nLFREEB = .TRUE.\n/\n").unwrap();
        assert!(nim.is_free_boundary().unwrap());
        // NCURR other than 1 is not the I'(s) profile.
        let nim = Indata::parse("&INDATA\nNCURR = 0\n/\n").unwrap();
        assert_eq!(nim.ncurr_is_iprime(), Ok(Some(false)));
        // Comments inside quoted strings survive.
        let nim = Indata::parse("&INDATA\nS = 'a!b'\n/\n").unwrap();
        assert_eq!(nim.string("s"), Ok("a!b"));
    }
}
