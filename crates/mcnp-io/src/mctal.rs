//! MCNP MCTAL output parsing — kcode subset.
//!
//! Reads the header and `kcode` criticality data; tally bodies are skipped.
//! No public MCTAL fixture corpus exists, so validation uses hand-built
//! synthetic files in `fixtures/mcnp/mctal/` exercising both the 5-value
//! and 19-value cycle record variants.

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

/// Parsed MCTAL kcode data.
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
    /// Declared tally numbers line (may be empty).
    pub tally_nums: Vec<u32>,
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
        let mut lines = text.lines();

        let head: Vec<&str> = lines
            .next()
            .ok_or_else(|| Error::BadStructure("empty file".into()))?
            .split_whitespace()
            .collect();
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

        let comment = lines.next().unwrap_or_default().trim().to_string();

        let tally_line: Vec<&str> = lines
            .next()
            .ok_or_else(|| Error::BadStructure("missing tally line".into()))?
            .split_whitespace()
            .collect();
        if tally_line.first() != Some(&"tally") && !tally_line.is_empty() {
            return Err(Error::BadStructure(format!(
                "expected `tally` line, got `{}`",
                tally_line.join(" ")
            )));
        }
        let n_tallies = tally_line.get(1).copied().unwrap_or("0").to_string();

        let tally_nums: Vec<u32> = match lines.next() {
            Some(l) => l
                .split_whitespace()
                .map(|t| num("tally number", t).map(|v| v as u32))
                .collect::<Result<_, _>>()?,
            None => Vec::new(),
        };

        // Tally bodies are skipped (upstream parity).
        // Note: full tally-body skipping requires format knowledge we don't
        // use downstream; like upstream we jump straight to the kcode line.
        let mut kcode: Vec<&str> = Vec::new();
        for l in lines.by_ref() {
            let w: Vec<&str> = l.split_whitespace().collect();
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

        let mut k_col = Vec::with_capacity(n_cycles);
        let mut k_abs = Vec::with_capacity(n_cycles);
        let mut k_path = Vec::with_capacity(n_cycles);
        let mut prompt_life_col = Vec::with_capacity(n_cycles);
        let mut prompt_life_path = Vec::with_capacity(n_cycles);
        let mut averages = Vec::new();

        for _ in 0..n_cycles {
            let mut values: Vec<f64> = Vec::with_capacity(per_cycle_values);
            for _ in 0..per_cycle_lines {
                let line = lines
                    .next()
                    .ok_or_else(|| Error::BadStructure("cycle data truncated".into()))?;
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
            tally_nums,
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
}
