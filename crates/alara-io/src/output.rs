//! ALARA activation-output listings parsed into [`ResponseFrame`].
//!
//! This mirrors the reference `alara_output_processing` Python tool shipped
//! with ALARA (`ALARA/tools/alara_output_processing.py`, vendored under
//! `.research/` for study): `*** <param> ***` headers open a response block,
//! `Interval #n ...` / `Zone #n: ...` / `Material #n: ...` lines select the
//! geometric block, and each `isotope ...` table contributes one
//! [`ResponseRow`] per (nuclide, cooling time), closed by its `total` row.
//!
//! Header mapping rules:
//! - `Number Density`, `Specific Activity`, `Total Decay Heat`, `Alpha Heat`,
//!   `Beta Heat`, `Gamma Heat`, `Contact Dose` map to [`ResponseVar`].
//! - `Folded ... Dose` (adjoint/biological dose) is folded into
//!   [`ResponseVar::ContactDose`]; the raw unit text is preserved in
//!   [`ResponseRow::var_unit`].
//! - `Photon Source Distribution ...` tables are skipped: photon spectra live
//!   in separate `.photonSrc` files parsed by [`crate::photon::PhotonSource`].
//! - `WDR` / `WDR/Clearance index ...` tables are skipped (clearance indices,
//!   not activation responses).
//! - Time columns map `pre-irrad` to `-1` s and `shutdown` to `0` s; every
//!   other column is a `<value> <unit>` pair converted with `s=1, m=60,
//!   h=3600, d=86400, w=7d, y=365d, c=100y` (same factors as the reference
//!   tool; note the deck-duration helper
//!   [`crate::schedule::parse_time_to_seconds`] instead defines a year as
//!   365.25 days).
//! - FISPACT-style `-1` block fields are not modeled: every row carries
//!   a real [`BlockKind`] / name / number triple.
//!
//! No pandas/arrow dependencies: rows live in a plain `Vec`.

use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::error::{Error, Result};
use nuclei::NuclideId;

/// Geometric resolution of one ALARA output table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BlockKind {
    /// Per-interval tables (`Interval #n (Zone: ...)`).
    Interval,
    /// Per-zone tables (`Zone #n: ...`).
    Zone,
    /// Per-material tables (`Material #n: ...`).
    Material,
}

impl BlockKind {
    /// Canonical block name, also accepted by [`ResponseFrame::filter`].
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Interval => "Interval",
            Self::Zone => "Zone",
            Self::Material => "Material",
        }
    }
}

impl std::fmt::Display for BlockKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Response variable tabulated by one ALARA output table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ResponseVar {
    /// `Number Density`.
    NumberDensity,
    /// `Specific Activity`.
    SpecificActivity,
    /// `Total Decay Heat`.
    TotalHeat,
    /// `Alpha Heat`.
    AlphaHeat,
    /// `Beta Heat`.
    BetaHeat,
    /// `Gamma Heat`.
    GammaHeat,
    /// `Contact Dose` (also covers `Folded ... Dose`; see module notes).
    ContactDose,
}

impl ResponseVar {
    /// Canonical ALARA variable name, also accepted by
    /// [`ResponseFrame::filter`].
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NumberDensity => "Number Density",
            Self::SpecificActivity => "Specific Activity",
            Self::TotalHeat => "Total Decay Heat",
            Self::AlphaHeat => "Alpha Heat",
            Self::BetaHeat => "Beta Heat",
            Self::GammaHeat => "Gamma Heat",
            Self::ContactDose => "Contact Dose",
        }
    }
}

impl std::fmt::Display for ResponseVar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One (nuclide, cooling time) response value from an ALARA output table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResponseRow {
    /// Cooling time in seconds (`-1` = pre-irradiation, `0` = shutdown).
    pub time_s: f64,
    /// Raw time-column label (`pre-irrad`, `shutdown`, `1 d`, ...).
    pub time_label: String,
    /// Nuclide name in ALARA dialect, kept verbatim (`h-1`, `fe-55`, `total`).
    pub nuclide: String,
    /// Half-life in seconds (`-1` = stable, `0` = `total` rows).
    pub half_life_s: f64,
    /// Run distinguisher supplied to [`ResponseFrame::parse`].
    pub run_lbl: String,
    /// Geometric resolution of the source table.
    pub block: BlockKind,
    /// Block name (`inner_zone`, `zone_0`, ...).
    pub block_name: String,
    /// Block position (1-based, as printed by ALARA).
    pub block_num: i64,
    /// Response variable of the source table.
    pub variable: ResponseVar,
    /// Unit parsed from `[...]` in the table header.
    pub var_unit: String,
    /// Tabulated response value.
    pub value: f64,
}

impl ResponseRow {
    /// Resolve [`Self::nuclide`] to a canonical [`NuclideId`].
    ///
    /// Fails for `total`/`Other` aggregate rows, which name no nuclide.
    pub fn nuclide_id(&self) -> std::result::Result<NuclideId, nuclei::Error> {
        NuclideId::from_name(&self.nuclide)
    }

    /// True for ALARA `total` aggregate rows.
    pub fn is_total(&self) -> bool {
        self.nuclide.eq_ignore_ascii_case("total")
    }
}

/// All response tables parsed from one ALARA output listing.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ResponseFrame {
    /// Rows in file order, one per (table, nuclide, cooling time).
    pub rows: Vec<ResponseRow>,
}

impl ResponseFrame {
    /// Parse an ALARA output listing from in-memory text.
    ///
    /// `run_lbl` tags every row (run distinguisher, e.g. `"sample2"`).
    /// Tables under skipped parameters (photon source, WDR) are ignored;
    /// blank input is an error.
    pub fn parse(text: &str, run_lbl: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty ALARA output".to_string(),
            });
        }
        let mut frame = Self { rows: Vec::new() };
        let mut param: Option<(ResponseVar, String)> = None;
        let mut block: Option<(BlockKind, String, i64)> = None;
        let mut header: Option<Vec<(f64, String)>> = None;
        let mut in_table = false;

        for (index, raw) in text.lines().enumerate() {
            let line_no = index + 1;
            let line = raw.trim();
            if line.is_empty() {
                // A blank line terminates any open table (its `total` row,
                // when present, always precedes the blank line).
                in_table = false;
                header = None;
                continue;
            }
            if line.starts_with("***") {
                param = map_parameter(line);
                in_table = false;
                header = None;
                continue;
            }
            if let Some(parsed) = parse_block_header(line) {
                block = Some(parsed);
                in_table = false;
                header = None;
                continue;
            }
            if is_block_header(line) {
                // Recognized prefix but unparsable number: stop attributing
                // tables until a well-formed block appears.
                block = None;
                in_table = false;
                header = None;
                continue;
            }
            if line.starts_with("isotope") {
                in_table = true;
                header = match (&param, &block) {
                    (Some(_), Some(_)) => Some(parse_time_columns(line, line_no)?),
                    _ => None,
                };
                continue;
            }
            if !in_table {
                continue;
            }
            if line.starts_with('=') {
                continue;
            }
            let Some(times) = header.as_ref() else {
                // Table under a skipped/unknown parameter: still consume
                // lines until its `total` row so parsing stays in sync.
                if line.starts_with("total") {
                    in_table = false;
                }
                continue;
            };
            let is_total_row = line.starts_with("total");
            let (nuclide, half_life_s, values) = parse_data_row(line, times.len(), line_no)?;
            let (Some((variable, var_unit)), Some((kind, block_name, block_num))) =
                (param.clone(), block.clone())
            else {
                continue;
            };
            for ((time_s, time_label), value) in times.iter().zip(values.iter()) {
                frame.rows.push(ResponseRow {
                    time_s: *time_s,
                    time_label: time_label.clone(),
                    nuclide: nuclide.clone(),
                    half_life_s,
                    run_lbl: run_lbl.to_string(),
                    block: kind,
                    block_name: block_name.clone(),
                    block_num,
                    variable,
                    var_unit: var_unit.clone(),
                    value: *value,
                });
            }
            if is_total_row {
                in_table = false;
                header = None;
            }
        }
        Ok(frame)
    }

    /// Read and parse an ALARA output listing from disk.
    pub fn from_file(path: impl AsRef<Path>, run_lbl: &str) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(&text, run_lbl)
    }

    /// Keep rows whose `field` equals `value`.
    ///
    /// Supported fields: `nuclide` (exact), `variable` (canonical ALARA name,
    /// case-insensitive), `block` (`Interval`/`Zone`/`Material`,
    /// case-insensitive), `block_name` (exact), `block_num` (integer),
    /// `run_lbl`, `time_label` (exact), `time` (seconds, float). Unknown
    /// fields match nothing and yield an empty frame.
    pub fn filter(&self, field: &str, value: &str) -> Self {
        let rows = self
            .rows
            .iter()
            .filter(|row| match field.to_ascii_lowercase().as_str() {
                "nuclide" => row.nuclide == value,
                "variable" => row.variable.as_str().eq_ignore_ascii_case(value),
                "block" => row.block.as_str().eq_ignore_ascii_case(value),
                "block_name" => row.block_name == value,
                "block_num" => value.parse::<i64>().is_ok_and(|num| row.block_num == num),
                "run_lbl" | "run" => row.run_lbl == value,
                "time_label" => row.time_label == value,
                "time" => value.parse::<f64>().is_ok_and(|time| row.time_s == time),
                _ => false,
            })
            .cloned()
            .collect();
        Self { rows }
    }

    /// Keep only the `total` aggregate rows of every table.
    pub fn totals(&self) -> Self {
        Self {
            rows: self
                .rows
                .iter()
                .filter(|row| row.is_total())
                .cloned()
                .collect(),
        }
    }

    /// Sum of `value` over [`ResponseVar::SpecificActivity`] rows.
    pub fn total_activity(&self) -> f64 {
        self.rows
            .iter()
            .filter(|row| row.variable == ResponseVar::SpecificActivity)
            .map(|row| row.value)
            .sum()
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// True when no rows were parsed.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}

/// Historical scaffold name for [`ResponseRow`].
///
/// Kept so the crate root's `ActivationRecord` re-export keeps resolving.
pub type ActivationRecord = ResponseRow;

/// Historical scaffold name for [`ResponseFrame`].
///
/// Kept so the crate root's `ActivationOutput` re-export keeps resolving.
pub type ActivationOutput = ResponseFrame;

/// Conversion factor to seconds for an ALARA cooling-time unit.
///
/// Single letters `s m h d w y c` plus common long names (case-insensitive),
/// with `y = 365 d` and `c = 100 y`, matching the reference Python tool.
pub(crate) fn time_unit_to_seconds(unit: &str) -> Option<f64> {
    match unit.trim().to_ascii_lowercase().as_str() {
        "s" | "sec" | "secs" | "second" | "seconds" => Some(1.0),
        "m" | "min" | "mins" | "minute" | "minutes" => Some(60.0),
        "h" | "hr" | "hrs" | "hour" | "hours" => Some(3_600.0),
        "d" | "day" | "days" => Some(86_400.0),
        "w" | "week" | "weeks" => Some(604_800.0),
        "y" | "yr" | "yrs" | "year" | "years" => Some(31_536_000.0),
        "c" | "century" | "centuries" => Some(3_153_600_000.0),
        _ => None,
    }
}

/// Map a `*** ... ***` header line to its response variable and unit.
///
/// Returns `None` for skipped parameters (photon source, WDR/clearance) and
/// for unrecognized headers such as schedule-hierarchy or target markers.
fn map_parameter(line: &str) -> Option<(ResponseVar, String)> {
    let content = line.trim_matches(|c| c == '*' || c == ' ').trim();
    let lower = content.to_ascii_lowercase();
    // Photon spectra live in separate `.photonSrc` files; WDR tables hold
    // clearance indices rather than activation responses.
    if lower.contains("photon source") || lower.contains("wdr") || lower.contains("clearance") {
        return None;
    }
    let (name, unit) = match content.split_once('[') {
        Some((name, rest)) => {
            let unit = rest.split(']').next().unwrap_or("").trim().to_string();
            (name.trim(), unit)
        }
        None => (content.trim(), String::new()),
    };
    let variable = match name.to_ascii_lowercase().as_str() {
        "number density" => ResponseVar::NumberDensity,
        "specific activity" => ResponseVar::SpecificActivity,
        "total decay heat" => ResponseVar::TotalHeat,
        "alpha heat" => ResponseVar::AlphaHeat,
        "beta heat" => ResponseVar::BetaHeat,
        "gamma heat" => ResponseVar::GammaHeat,
        "contact dose" => ResponseVar::ContactDose,
        _ if name.to_ascii_lowercase().starts_with("folded") => {
            // Folded (adjoint/biological) dose shares the contact-dose
            // semantics for downstream filtering; the raw unit is kept.
            ResponseVar::ContactDose
        }
        _ => return None,
    };
    Some((variable, unit))
}

/// True when `line` opens an interval/zone/material block (well-formed or not).
fn is_block_header(line: &str) -> bool {
    line.starts_with("Interval #") || line.starts_with("Zone #") || line.starts_with("Material #")
}

/// Parse `Interval #n (Zone: x)` / `Zone #n: x` / `Material #n: x` lines.
///
/// Returns `(kind, block_name, block_num)`; the trailing `:` is optional and
/// parenthesized `(Zone: x)` wrappers are unwrapped to `x`.
fn parse_block_header(line: &str) -> Option<(BlockKind, String, i64)> {
    let (kind, rest) = if let Some(rest) = line.strip_prefix("Interval #") {
        (BlockKind::Interval, rest)
    } else if let Some(rest) = line.strip_prefix("Zone #") {
        (BlockKind::Zone, rest)
    } else {
        (BlockKind::Material, line.strip_prefix("Material #")?)
    };
    let rest = rest.trim().trim_end_matches(':').trim();
    let num_end = rest
        .find(|c: char| c.is_whitespace() || c == ':')
        .unwrap_or(rest.len());
    let block_num: i64 = rest[..num_end].parse().ok()?;
    let mut name = rest[num_end..]
        .trim()
        .trim_start_matches(':')
        .trim()
        .to_string();
    if name.starts_with('(') && name.ends_with(')') && name.len() >= 2 {
        name = name[1..name.len() - 1].trim().to_string();
    }
    if let Some((_, tail)) = name.rsplit_once(':') {
        name = tail.trim().to_string();
    }
    Some((kind, name, block_num))
}

/// Parse the time columns of an `isotope t_1/2(s) ...` header line.
fn parse_time_columns(line: &str, line_no: usize) -> Result<Vec<(f64, String)>> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("expected `isotope t_1/2(s) ...`, found `{line}`"),
        });
    }
    let mut rest = tokens[2..].iter().peekable();
    let mut times = Vec::new();
    if rest
        .peek()
        .is_some_and(|token| token.eq_ignore_ascii_case("pre-irrad"))
    {
        times.push((-1.0, "pre-irrad".to_string()));
        rest.next();
    }
    if rest
        .peek()
        .is_some_and(|token| token.eq_ignore_ascii_case("shutdown"))
    {
        times.push((0.0, "shutdown".to_string()));
        rest.next();
    }
    let rest: Vec<&&str> = rest.collect();
    if rest.len() % 2 != 0 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("odd `<value> <unit>` time tokens in `{line}`"),
        });
    }
    for pair in rest.chunks_exact(2) {
        let (value_text, unit) = (*pair[0], *pair[1]);
        let value: f64 = value_text.parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected cooling time, found `{value_text}`"),
        })?;
        let factor = time_unit_to_seconds(unit).ok_or_else(|| Error::Parse {
            line: line_no,
            msg: format!("unknown cooling-time unit `{unit}`"),
        })?;
        times.push((value * factor, format!("{value_text} {unit}")));
    }
    if times.is_empty() {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("no time columns in `{line}`"),
        });
    }
    Ok(times)
}

/// Parse one table data (or `total`) row into nuclide, half-life, and values.
fn parse_data_row(line: &str, n_times: usize, line_no: usize) -> Result<(String, f64, Vec<f64>)> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() != n_times + 2 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!(
                "expected {} entries (nuclide, half-life, {} values), found `{line}`",
                n_times + 2,
                n_times
            ),
        });
    }
    let nuclide = tokens[0].to_string();
    let half_life_s = match tokens[1].to_ascii_lowercase().as_str() {
        "stable" => -1.0,
        "none" => 0.0,
        _ => tokens[1].parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected half-life, found `{}`", tokens[1]),
        })?,
    };
    // `total`/`Other` aggregate rows name no nuclide; everything else must
    // resolve via the shared nuclide table (names kept verbatim).
    if !nuclide.eq_ignore_ascii_case("total") && !nuclide.eq_ignore_ascii_case("other") {
        NuclideId::from_name(&nuclide).map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("unknown nuclide `{nuclide}`"),
        })?;
    }
    let mut values = Vec::with_capacity(n_times);
    for token in &tokens[2..] {
        let value: f64 = token.parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected response value, found `{token}`"),
        })?;
        values.push(value);
    }
    Ok((nuclide, half_life_s, values))
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../fixtures/alara/output/sample2.out"
    );

    fn approx(actual: f64, expected: f64) -> bool {
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0)
    }

    #[test]
    fn parse_sample2_fixture_covers_both_params_and_intervals() {
        let frame = ResponseFrame::from_file(FIXTURE, "sample2").unwrap();
        // 4 tables x 6 cooling times: (194+1) + (100+1) + (121+1) + (54+1).
        assert_eq!(frame.len(), 2838);
        assert!(!frame.is_empty());

        let variables: std::collections::BTreeSet<String> = frame
            .rows
            .iter()
            .map(|row| row.variable.to_string())
            .collect();
        assert_eq!(
            variables,
            [
                "Number Density".to_string(),
                "Specific Activity".to_string()
            ]
            .into_iter()
            .collect()
        );

        assert!(frame
            .rows
            .iter()
            .all(|row| row.block == BlockKind::Interval));
        let blocks: std::collections::BTreeSet<(i64, String)> = frame
            .rows
            .iter()
            .map(|row| (row.block_num, row.block_name.clone()))
            .collect();
        assert_eq!(
            blocks,
            [(1, "inner_zone".to_string()), (2, "outer_zone".to_string())]
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
        );

        let mut times: Vec<f64> = frame.rows.iter().map(|row| row.time_s).collect();
        times.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        times.dedup();
        // pre-irrad, shutdown, 1 d, 100 d, 1 y, 100 y (y = 365 d).
        assert_eq!(
            times,
            [
                -1.0,
                0.0,
                86_400.0,
                8_640_000.0,
                31_536_000.0,
                3_153_600_000.0
            ]
        );

        assert!(frame.rows.iter().all(|row| row.run_lbl == "sample2"));

        // Spot-check a pre-irradiation number density and a total row.
        let h1 = frame
            .filter("nuclide", "h-1")
            .filter("variable", "Number Density")
            .filter("block_num", "1")
            .filter("time", "-1");
        assert_eq!(h1.len(), 1);
        assert!(approx(h1.rows[0].value, 6.6354e21));
        assert_eq!(h1.rows[0].half_life_s, -1.0);
        assert_eq!(h1.rows[0].time_label, "pre-irrad");
        assert_eq!(h1.rows[0].var_unit, "atoms/cm3");

        let totals = frame.totals();
        // One `total` row per table x 6 cooling times.
        assert_eq!(totals.len(), 24);
        assert!(totals.rows.iter().all(|row| row.is_total()));
        assert!(totals.rows.iter().all(|row| row.half_life_s == 0.0));

        // Every table contributes pre-irrad/shutdown/cooling columns.
        for row in frame.filter("nuclide", "total").rows {
            assert!(row.time_s == -1.0 || row.time_s >= 0.0);
        }

        assert!(frame.total_activity() > 0.0);
        assert_eq!(
            frame.filter("block", "Zone").len(),
            0,
            "sample2 holds interval output only"
        );
    }

    #[test]
    fn folded_dose_maps_to_contact_dose() {
        let text = "*** Folded (Adjoint/Biological) Dose [units defined by adjoint flux] : data/adjfile2.gam ***\n\
                      \n\
                      Interval #1 (Zone: inner_zone) :\n\
                      Folded (Adjoint/Biological) Dose [units defined by adjoint flux] : data/adjfile2.gam\n\
                      isotope  t_1/2(s)   pre-irrad   shutdown         1 d\n\
                      ======================================================================================\n\
                      mn-56 \t9.2844e+03  0.0000e+00  1.0000e+05  2.0000e+04\n\
                      ======================================================================================\n\
                      total   0           0.0000e+00  1.0000e+05  2.0000e+04\n";
        let frame = ResponseFrame::parse(text, "fold").unwrap();
        assert_eq!(frame.len(), 6);
        assert!(frame
            .rows
            .iter()
            .all(|row| row.variable == ResponseVar::ContactDose));
        assert_eq!(frame.rows[0].var_unit, "units defined by adjoint flux");
        assert_eq!(frame.totals().len(), 3);
    }

    #[test]
    fn photon_source_and_wdr_tables_are_skipped() {
        let text = "*** Photon Source Distribution [gammas/s/cm] : output/sample9.photonSrc\n\
                     	    with Specific Activity [Bq/cm] ***\n\
                     \n\
                     Zone #1: zone_0\n\
                     Photon Source Distribution [gammas/s/cm] : output/sample9.photonSrc\n\
                     	    with Specific Activity [Bq/cm]\n\
                     isotope  t_1/2(s)   pre-irrad   shutdown\n\
                     ======================================================================================\n\
                     mn-56 \t9.2844e+03  8.6090e+22  6.4272e+18\n\
                     ======================================================================================\n\
                     total   0           8.6090e+22  6.4272e+18\n\
                     \n\
                     *** WDR/Clearance index: data/NRCA ***\n\
                     \n\
                     Zone #1: zone_0\n\
                     WDR/Clearance index: data/NRCA\n\
                     isotope  t_1/2(s)   pre-irrad   shutdown\n\
                     ======================================================================================\n\
                     co-60 \t1.6636e+08  0.0000e+00  8.4814e+00\n\
                     ======================================================================================\n\
                     total   0           0.0000e+00  8.4814e+00\n";
        let frame = ResponseFrame::parse(text, "skip").unwrap();
        assert!(frame.is_empty());
        assert_eq!(frame.total_activity(), 0.0);
    }

    #[test]
    fn zone_and_material_blocks_parse() {
        let text = "*** Specific Activity [Bq/cm3] ***\n\
                      Zone #2: outer_zone\n\
                      isotope  t_1/2(s)   shutdown      5000 s\n\
                      =====\n\
                      mn-56 \t9.2844e+03  1.0000e+05  2.0000e+04\n\
                      =====\n\
                      total   0           1.0000e+05  2.0000e+04\n\
                      \n\
                      *** Number Density [atoms/cm3] ***\n\
                      Material #3: wall_mat\n\
                      isotope  t_1/2(s)   pre-irrad   shutdown\n\
                      =====\n\
                      fe-55 \t8.6314e+07  0.0000e+00  7.4256e+14\n\
                      =====\n\
                      total   0           0.0000e+00  7.4256e+14\n";
        let frame = ResponseFrame::parse(text, "blocks").unwrap();
        assert_eq!(frame.len(), 8);
        let zone = frame.filter("block", "Zone");
        assert_eq!(zone.len(), 4);
        assert_eq!(zone.rows[0].block_name, "outer_zone");
        assert_eq!(zone.rows[0].block_num, 2);
        assert_eq!(zone.rows[0].time_s, 0.0);
        assert_eq!(zone.rows[0].time_label, "shutdown");
        assert_eq!(frame.filter("block", "Material").len(), 4);
        // Nuclide names stay verbatim; metastable-free names resolve.
        assert_eq!(
            frame.filter("nuclide", "fe-55").rows[0]
                .nuclide_id()
                .unwrap(),
            NuclideId::from_name("fe-55").unwrap()
        );
    }

    #[test]
    fn parse_rejects_bad_nuclides_times_and_shapes() {
        let good_header = "*** Number Density [atoms/cm3] ***\nInterval #1 (Zone: inner_zone) :\n";
        assert!(matches!(
            ResponseFrame::parse(
                &format!("{good_header}isotope  t_1/2(s)   shutdown\n=====\nXx999 \t1.0  2.0\n=====\ntotal   0  2.0\n"),
                "r"
            ),
            Err(Error::Parse { .. })
        ));
        // `total` and `Other` aggregates are exempt from nuclide validation.
        let frame = ResponseFrame::parse(
            &format!("{good_header}isotope  t_1/2(s)   shutdown\n=====\nOther   0  1.0\n=====\ntotal   0  1.0\n"),
            "r",
        )
        .unwrap();
        assert_eq!(frame.len(), 2);
        // Odd time tokens, unknown units, and short rows are errors.
        for header in [
            "isotope  t_1/2(s)   shutdown   1",
            "isotope  t_1/2(s)   shutdown   1 fortnight",
            "isotope  t_1/2(s)",
        ] {
            assert!(
                matches!(
                    ResponseFrame::parse(&format!("{good_header}{header}\n"), "r"),
                    Err(Error::Parse { .. })
                ),
                "header `{header}` should fail"
            );
        }
        assert!(matches!(
            ResponseFrame::parse(
                &format!("{good_header}isotope  t_1/2(s)   shutdown   1 d\n=====\nh-1 \t-1  1.0\n=====\ntotal   0  1.0  2.0\n"),
                "r"
            ),
            Err(Error::Parse { .. })
        ));
        // Empty input is an error; unknown filter fields match nothing.
        assert!(matches!(
            ResponseFrame::parse("  \n", "r"),
            Err(Error::Parse { .. })
        ));
        let frame = ResponseFrame::parse(
            &format!("{good_header}isotope  t_1/2(s)   shutdown\n=====\nh-1 \t-1  1.0\n=====\ntotal   0  1.0\n"),
            "r",
        )
        .unwrap();
        assert!(frame.filter("bogus", "x").is_empty());
        assert_eq!(frame.filter("time", "0").len(), 2);
    }
}
