//! XSDIR cross-section index parsing (`Xsdir` / `XsdirTable`),
//! validated against the `dummy_xsdir` fixture.
//!
//! Documented deviations from the legacy reader:
//! - `awr` keys/values are typed (`u32` zaid → `f64` ratio) instead of raw strings.
//! - Optional table fields are `Option<_>` rather than `None`-able Python
//!   attributes; `to_serpent` is fallible when temperature is absent.
//! - Serpent output reproduces Python's `{:.11e}` exponent padding
//!   (`6.44688328094e+15`) byte-for-byte.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::path::Path;

use nuclei::dialects;
use nuclei::NuclideId;

/// Errors raised while parsing or converting xsdir data.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    /// A required header word/section was missing or malformed.
    BadHeader(String),
    /// A numeric field failed to parse.
    BadNumber {
        field: &'static str,
        text: String,
    },
    /// Not enough fields on a directory entry line.
    TooFewFields {
        got: usize,
    },
    /// `to_serpent` needs a temperature; this table has none.
    MissingTemperature,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadHeader(m) => write!(f, "malformed xsdir header: {m}"),
            Error::BadNumber { field, text } => {
                write!(f, "cannot parse {field} from `{text}`")
            }
            Error::TooFewFields { got } => {
                write!(f, "directory entry needs >= 7 fields, found {got}")
            }
            Error::MissingTemperature => {
                write!(
                    f,
                    "table has no temperature; cannot convert to Serpent form"
                )
            }
        }
    }
}

impl std::error::Error for Error {}

fn num_f64(field: &'static str, text: &str) -> Result<f64, Error> {
    text.parse::<f64>().map_err(|_| Error::BadNumber {
        field,
        text: text.to_string(),
    })
}

fn num_i64(field: &'static str, text: &str) -> Result<i64, Error> {
    text.parse::<i64>().map_err(|_| Error::BadNumber {
        field,
        text: text.to_string(),
    })
}

/// One directory entry describing a cross-section table
/// (MCNP5 User's Guide Vol. 3, App. K field names).
#[derive(Debug, Clone, PartialEq)]
pub struct XsdirTable {
    /// ZAID + library id, delimited by `.` (e.g. `"1001.44c"`).
    pub name: String,
    /// Atomic mass ratio.
    pub awr: f64,
    /// Relative path of the file holding the table data.
    pub filename: String,
    /// Access route string; typically `"0"`.
    pub access: String,
    /// 1 = formatted (ascii), 2 = unformatted (binary).
    pub filetype: i64,
    /// Line number (filetype 1) or record number (filetype 2).
    pub address: i64,
    /// Length of the second block of the data table.
    pub tablelength: i64,
    /// Binary only: bytes per record times entries per record.
    pub recordlength: Option<i64>,
    /// Binary only: number of entries per record.
    pub entries: Option<i64>,
    /// Temperature in MeV (neutron data).
    pub temperature: Option<f64>,
    /// True for continuous-energy neutron data with unresolved-resonance
    /// probability tables.
    pub ptable: bool,
}

impl XsdirTable {
    /// The ZAID part of [`Self::name`] (text before the first `.`).
    pub fn zaid(&self) -> &str {
        &self.name[..self.name.find('.').unwrap_or(self.name.len())]
    }

    /// Alias property: the full table name.
    pub fn alias(&self) -> &str {
        &self.name
    }

    /// Serpent table type: 1 continuous (`c`), 2 dosimetry (`y`),
    /// 3 thermal (`t`), else `None`.
    pub fn serpent_type(&self) -> Option<u8> {
        match self.name.chars().last()? {
            'c' => Some(1),
            'y' => Some(2),
            't' => Some(3),
            _ => None,
        }
    }

    /// Metastable flag heuristic: special-cases Am-242 zaids,
    /// otherwise `A > 600` counts as metastable. Only meaningful for
    /// continuous-energy (`c`) tables, hence `Option`.
    pub fn metastable(&self) -> Option<bool> {
        if !self.name.ends_with('c') {
            return None;
        }
        match self.zaid() {
            "95242" => Some(true),
            "95642" => Some(false),
            _ => {
                let full: u32 = self.zaid().parse::<u32>().ok()?;
                Some(full % 1000 > 600)
            }
        }
    }

    /// Serpent directory-entry line. Reproduces the reference formatting exactly,
    /// including Python-style `{:.11e}` exponent padding.
    pub fn to_serpent(&self, directory: &str) -> Result<String, Error> {
        let stype = self.serpent_type().ok_or(Error::MissingTemperature)?; // same failure class as upstream None deref
        let temp_k = self.temperature.ok_or(Error::MissingTemperature)? / 8.617_342_3e-11;
        let dir = if directory.is_empty() {
            String::new()
        } else if directory.ends_with('/') {
            directory.to_string()
        } else {
            format!("{directory}/")
        };
        Ok(format!(
            "{name} {name} {stype} {zaid} {meta} {awr} {tk} {ft} {dir}{file}",
            name = self.name,
            stype = stype,
            zaid = self.zaid(),
            meta = i32::from(self.metastable().unwrap_or(false)),
            awr = self.awr,
            tk = py_exp(temp_k, 11),
            ft = self.filetype - 1,
            dir = dir,
            file = self.filename,
        ))
    }
}

/// Format a float like Python's `format(x, '.{prec}e')`: mantissa with
/// `prec` fraction digits, exponent with sign and at least two digits.
fn py_exp(value: f64, prec: usize) -> String {
    let s = format!("{value:.prec$e}");
    // Rust: "6.44688328094e15"; split at 'e'
    let (mant, exp) = s.split_once('e').expect("lowercase e always present");
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    if digits.len() < 2 {
        format!("{mant}e{sign}0{digits}")
    } else {
        format!("{mant}e{sign}{digits}")
    }
}

/// Parsed contents of an MCNP xsdir index file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Xsdir {
    /// `DATAPATH=` value from the first line, if present.
    pub datapath: Option<String>,
    /// Path the file was read from (`from_file` only).
    pub source_path: Option<String>,
    /// Atomic weight ratios keyed by zaid.
    pub awr: BTreeMap<u32, f64>,
    /// Directory entries in file order.
    pub tables: Vec<XsdirTable>,
}

impl Xsdir {
    /// Read and parse an xsdir file from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        let mut xsdir = Xsdir::parse(&text)?;
        xsdir.source_path = Some(path.display().to_string());
        Ok(xsdir)
    }

    /// Parse xsdir text in memory.
    ///
    /// Line discipline mirrors the legacy reader exactly: line 1 may be blank or carry
    /// `DATAPATH=`; line 2 must be `atomic weight ratios`; AWR pairs run
    /// until an odd-count line or the `directory` marker; a blank line
    /// terminates the directory entries.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let mut lines = text.lines();

        // First section: optional DATAPATH=... (may itself be blank).
        let first = lines
            .next()
            .ok_or_else(|| Error::BadHeader("empty file".into()))?;
        let mut datapath = None;
        let words: Vec<&str> = first.split_whitespace().collect();
        if let Some(w0) = words.first() {
            if w0.to_ascii_lowercase().starts_with("datapath") {
                let idx = first
                    .find('=')
                    .ok_or_else(|| Error::BadHeader("DATAPATH line lacks '='".into()))?;
                datapath = Some(first[idx + 1..].trim().to_string());
            }
        }

        // Second section header: "atomic weight ratios"
        let awr_line = lines
            .next()
            .ok_or_else(|| Error::BadHeader("missing AWR section".into()))?;
        let awr_words: Vec<&str> = awr_line.split_whitespace().collect();
        if awr_words.len() != 3
            || !awr_words[0].eq_ignore_ascii_case("atomic")
            || !awr_words[1].eq_ignore_ascii_case("weight")
            || !awr_words[2].eq_ignore_ascii_case("ratios")
        {
            return Err(Error::BadHeader(format!(
                "expected `atomic weight ratios`, got `{awr_line}`"
            )));
        }

        // AWR pairs until an odd-count line or the `directory` marker.
        let mut awr = BTreeMap::new();
        let mut breaker: Vec<&str>;
        loop {
            let line = match lines.next() {
                Some(l) => l,
                None => return Err(Error::BadHeader("no `directory` section found".into())),
            };
            let w: Vec<&str> = line.split_whitespace().collect();
            if w.is_empty() {
                return Err(Error::BadHeader("blank line inside AWR section".into()));
            }
            if w.len() % 2 != 0 || w.first() == Some(&"directory") {
                breaker = w;
                break;
            }
            for pair in w.chunks(2) {
                let zaid = num_i64("awr zaid", pair[0])? as u32;
                let ratio = num_f64("awr", pair[1])?;
                awr.insert(zaid, ratio);
            }
        }

        // Consume any further pre-directory lines: the AWR
        // section may end on an odd-count line before `directory` appears.
        while breaker.first() != Some(&"directory") {
            let line = match lines.next() {
                Some(l) => l,
                None => return Err(Error::BadHeader("no `directory` section found".into())),
            };
            breaker = line.split_whitespace().collect();
        }

        // Directory entries until EOF or a blank line.
        let mut tables = Vec::new();
        while let Some(line) = lines.next() {
            let mut w: Vec<&str> = line.split_whitespace().collect();
            if w.is_empty() {
                break;
            }
            while w.last() == Some(&"+") {
                w.pop();
                let cont = match lines.next() {
                    Some(l) => l,
                    None => return Err(Error::BadHeader("continuation '+' ends file".into())),
                };
                w.extend(cont.split_whitespace());
            }
            if w.len() < 7 {
                return Err(Error::TooFewFields { got: w.len() });
            }
            let mut t = XsdirTable {
                name: w[0].to_string(),
                awr: num_f64("table awr", w[1])?,
                filename: w[2].to_string(),
                access: w[3].to_string(),
                filetype: num_i64("filetype", w[4])?,
                address: num_i64("address", w[5])?,
                tablelength: num_i64("tablelength", w[6])?,
                recordlength: None,
                entries: None,
                temperature: None,
                ptable: false,
            };
            if w.len() > 7 {
                t.recordlength = Some(num_i64("recordlength", w[7])?);
            }
            if w.len() > 8 {
                t.entries = Some(num_i64("entries", w[8])?);
            }
            if w.len() > 9 {
                t.temperature = Some(num_f64("temperature", w[9])?);
            }
            if w.len() > 10 && w[10] == "ptable" {
                t.ptable = true;
            }
            tables.push(t);
        }

        Ok(Xsdir {
            datapath,
            source_path: None,
            awr,
            tables,
        })
    }

    /// All tables whose name contains `name` (substring semantics).
    pub fn find_table(&self, name: &str) -> Vec<&XsdirTable> {
        self.tables
            .iter()
            .filter(|t| t.name.contains(name))
            .collect()
    }

    /// Distinct nuclides referenced by the directory entries.
    pub fn nucs(&self) -> BTreeSet<NuclideId> {
        self.tables
            .iter()
            .filter_map(|t| t.zaid().parse::<u32>().ok())
            .filter_map(|z| dialects::from_zaid(z).ok())
            .collect()
    }

    /// Serpent xsdata lines for all continuous-energy tables.
    pub fn xsdata_lines(&self) -> Result<Vec<String>, Error> {
        self.tables
            .iter()
            .filter(|t| t.serpent_type() == Some(1))
            .map(|t| t.to_serpent(""))
            .collect()
    }

    /// Write a Serpent xsdata file for all continuous-energy tables.
    pub fn write_xsdata(&self, path: impl AsRef<Path>) -> Result<(), Error> {
        let body = self.xsdata_lines()?;
        let mut out = String::new();
        for line in body {
            out.push_str(&line);
            out.push('\n');
        }
        std::fs::write(path, out).map_err(|e| Error::Io(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/xsdir/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn gen_xsdir() -> Xsdir {
        Xsdir::from_file(fixture_path("dummy_xsdir")).unwrap()
    }

    #[test]
    fn parse_awr_table() {
        let x = gen_xsdir();
        let expected: [(u32, f64); 7] = [
            (1, 1.000000),
            (1000, 0.99931697),
            (3000, 6.88131188),
            (3003, 3.11111111),
            (3004, 4.11111111),
            (3005, 5.111111111),
            (3009, 9.11111111),
        ];
        assert_eq!(x.awr.len(), expected.len());
        for (zaid, v) in expected {
            assert_eq!(x.awr[&zaid], v, "awr[{zaid}]");
        }
    }

    #[test]
    fn parse_directory_tables() {
        let x = gen_xsdir();
        assert_eq!(x.tables.len(), 3);

        let t0 = &x.tables[0];
        assert_eq!(t0.name, "1001.44c");
        assert_eq!(t0.awr, 1.111111);
        assert_eq!(t0.filename, "many_xs/1001.555nc");
        assert_eq!(t0.access, "0");
        assert_eq!(t0.filetype, 1);
        assert_eq!(t0.address, 4);
        assert_eq!(t0.tablelength, 55555);
        assert_eq!(t0.recordlength, Some(0));
        assert_eq!(t0.entries, Some(0));
        assert_eq!(t0.temperature, Some(5.5555e5));
        assert!(!t0.ptable);

        assert_eq!(x.tables[1].filename, "such_data/1001.777nc");
        assert!(x.tables[1].ptable);
        // ptable keyword lands on the '+' continuation line here
        assert_eq!(x.tables[2].filename, "more_data/1001.999nc");
        assert!(x.tables[2].ptable);
    }

    #[test]
    fn find_table_substring() {
        let x = gen_xsdir();
        let hits = x.find_table("1001");
        let names: Vec<&str> = hits.iter().map(|t| t.name.as_str()).collect();
        assert_eq!(names, ["1001.44c", "1001.66c", "1001.70c"]);
        assert!(x.find_table("9999").is_empty());
    }

    #[test]
    fn serpent_conversion_exact_bytes() {
        let x = gen_xsdir();
        let line = x.tables[0].to_serpent(".").unwrap();
        assert_eq!(
            line,
            "1001.44c 1001.44c 1 1001 0 1.111111 6.44688328094e+15 0 ./many_xs/1001.555nc"
        );
    }

    #[test]
    fn xsdata_lines_match_oracle() {
        let x = gen_xsdir();
        let lines = x.xsdata_lines().unwrap();
        assert_eq!(
            lines,
            vec![
                "1001.44c 1001.44c 1 1001 0 1.111111 6.44688328094e+15 0 many_xs/1001.555nc",
                "1001.66c 1001.66c 1 1001 0 1.111111 6.44688328094e+15 0 such_data/1001.777nc",
                "1001.70c 1001.70c 1 1001 0 1.111111 6.44688328094e+15 0 more_data/1001.999nc",
            ]
        );
    }

    #[test]
    fn nucs_set() {
        let x = gen_xsdir();
        let nucs = x.nucs();
        assert_eq!(nucs.len(), 1);
        assert_eq!(nucs.iter().next().unwrap().nucid(), 10010000);
    }

    #[test]
    fn zaid_and_metastable_helpers() {
        let x = gen_xsdir();
        assert_eq!(x.tables[0].zaid(), "1001");
        assert_eq!(x.tables[0].metastable(), Some(false));
        assert_eq!(x.tables[0].serpent_type(), Some(1));
    }

    #[test]
    fn datapath_header_parsed() {
        let text = "datapath=/some/data\natomic weight ratios\ndirectory\n";
        let x = Xsdir::parse(text).unwrap();
        assert_eq!(x.datapath.as_deref(), Some("/some/data"));
        assert!(x.tables.is_empty());
    }

    #[test]
    fn missing_directory_section_errors() {
        let text = "\natomic weight ratios\n1000 1.0\n";
        assert!(matches!(
            Xsdir::parse(text),
            Err(Error::BadHeader(msg)) if msg.contains("directory")
        ));
    }

    #[test]
    fn bad_awr_header_errors() {
        let text = "\natomic weight ration\n";
        assert!(matches!(Xsdir::parse(text), Err(Error::BadHeader(_))));
    }

    #[test]
    fn short_entry_errors() {
        let text = "\natomic weight ratios\ndirectory\n1001.44c 1.0 file 0 1\n";
        assert!(matches!(
            Xsdir::parse(text),
            Err(Error::TooFewFields { got: 5 })
        ));
    }

    #[test]
    fn bad_number_errors() {
        let text = "\natomic weight ratios\ndirectory\n1001.44c xx file 0 1 4 55555\n";
        assert!(matches!(
            Xsdir::parse(text),
            Err(Error::BadNumber {
                field: "table awr",
                ..
            })
        ));
    }

    #[test]
    fn py_exp_matches_python_padding() {
        assert_eq!(py_exp(6.44688328_094e15, 11), "6.44688328094e+15");
        assert_eq!(py_exp(-2.5e-3, 3), "-2.500e-03");
        assert_eq!(py_exp(1.0, 1), "1.0e+00");
    }
}
