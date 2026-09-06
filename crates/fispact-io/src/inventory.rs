//! FISPACT-II inventory-table parser producing ALARA-compatible response frames.
//!
//! Accepted synthetic layout (see `fixtures/fispact/inventory.fis`):
//! - Cooling-step header: any line containing `time` (case-insensitive) that
//!   carries a `<value> <unit>` pair, e.g. `TIME = 1.0 DAYS`. The characters
//!   `= : , ;` are treated as whitespace before tokenising, and the *last*
//!   `<value> <unit>` pair on the line wins.
//! - Column-title lines (leading `nuclide`/`isotope`, or a line containing
//!   both `atoms` and `activity`) and separator lines (only the characters
//!   `= - * _ . + | / \` plus whitespace) are skipped.
//! - Data rows: `<nuclide> <atoms> <activity_Bq> <heat_W> [...]` free-format
//!   whitespace-separated floats; extra trailing columns (e.g. dose) are
//!   ignored.
//!
//! Mapping to [`nucleide_alara_io::output::ResponseRow`] (11 fields):
//! - `time_s` / `time_label`: converted seconds / `<value> <unit>` label that
//!   preserves the file spelling of both tokens (e.g. `1.0 DAYS`).
//! - `nuclide`: kept verbatim (`h-3`, `fe-55`, `total`); every name except
//!   `total` (case-insensitive) is validated with
//!   [`nucleide_nuclei::NuclideId::from_name`], so lowercase-dash ALARA-style names
//!   resolve while typos fail.
//! - `half_life_s`: always `-1.0`. The inventory table carries no half-life
//!   column, so every row (including `total`) is marked unknown. This differs
//!   from ALARA, where `-1.0` means stable and `0.0` marks `total` rows.
//! - `run_lbl`: caller-supplied tag copied to every row.
//! - `block` / `block_name` / `block_num`: FISPACT has no zone/interval
//!   blocks, so every row uses [`nucleide_alara_io::output::BlockKind::Material`] with
//!   name `"inventory"` and number `-1`. This mirrors the upstream
//!   `alara_output_processing` FISPACT passthrough convention of `-1` block
//!   fields, but keeps the typed enum instead of a raw `-1`.
//! - `variable` / `var_unit` / `value`: one row per variable —
//!   `NumberDensity` in `atoms`, `SpecificActivity` in `Bq`, `TotalHeat` in `W`.
//!
//! Time units (case-insensitive; singular/plural plus common abbreviations):
//! `s = 1`, `m = 60`, `h = 3_600`, `d = 86_400`, `w = 604_800 (7 d)`,
//! `y = 31_536_000 (365 d)`. The 365-day year matches the ALARA reference
//! tool (not 365.25); months and centuries are not accepted here.
//!
//! Skipped: zero-length/blank lines, `#...` and Fortran-style `C ...` comment
//! lines, separator and column-title lines, and non-data prose before the
//! first cooling step. Anything else that looks like a data row outside a
//! cooling step, or any malformed header/row inside a step, is a 1-based
//! [`crate::error::Error::Parse`] error. Whitespace-only input errors at
//! line 1.

use std::path::Path;

use nucleide_alara_io::output::{BlockKind, ResponseFrame, ResponseRow, ResponseVar};
use nucleide_nuclei::NuclideId;

use crate::error::{Error, Result};

/// Block name used for every FISPACT inventory row (no zone/interval blocks).
const BLOCK_NAME: &str = "inventory";
/// Block number used for every FISPACT inventory row (passthrough `-1`).
const BLOCK_NUM: i64 = -1;
/// Half-life marker: the inventory table carries no half-life column.
const UNKNOWN_HALF_LIFE_S: f64 = -1.0;

/// Parse FISPACT-II inventory output text into a response frame.
///
/// `run_lbl` tags every row (run distinguisher, e.g. `"sample"`). See the
/// module docs for the accepted layout and the field-mapping rules.
pub fn parse_to_frame(text: &str, run_lbl: &str) -> Result<ResponseFrame> {
    if text.trim().is_empty() {
        return Err(Error::Parse {
            line: 1,
            msg: "empty FISPACT-II output".to_string(),
        });
    }
    let mut rows: Vec<ResponseRow> = Vec::new();
    let mut current: Option<(f64, String)> = None;

    for (index, raw) in text.lines().enumerate() {
        let line_no = index + 1;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("C ") {
            continue;
        }
        if is_separator(line) || is_column_header(line) {
            continue;
        }
        if line.to_ascii_lowercase().contains("time") {
            let normalized = line.replace(['=', ':', ',', ';'], " ");
            let tokens: Vec<&str> = normalized.split_whitespace().collect();
            if !tokens.iter().any(|t| t.parse::<f64>().is_ok()) {
                // Prose/title line mentioning time but carrying no numeric
                // value (e.g. `FISPACT-II RUN TIME SUMMARY`).
                continue;
            }
            let (time_s, time_label) = parse_cooling_header(&tokens, line_no)?;
            current = Some((time_s, time_label));
            continue;
        }
        match current.clone() {
            None => {
                if looks_like_data_row(line) {
                    return Err(Error::Parse {
                        line: line_no,
                        msg: format!("nuclide row outside cooling step: `{line}`"),
                    });
                }
                // Preamble prose/title before the first cooling step.
            }
            Some((time_s, time_label)) => {
                let (nuclide, atoms, activity, heat) = parse_data_row(line, line_no)?;
                for (variable, var_unit, value) in [
                    (ResponseVar::NumberDensity, "atoms", atoms),
                    (ResponseVar::SpecificActivity, "Bq", activity),
                    (ResponseVar::TotalHeat, "W", heat),
                ] {
                    rows.push(ResponseRow {
                        time_s,
                        time_label: time_label.clone(),
                        nuclide: nuclide.clone(),
                        half_life_s: UNKNOWN_HALF_LIFE_S,
                        run_lbl: run_lbl.to_string(),
                        block: BlockKind::Material,
                        block_name: BLOCK_NAME.to_string(),
                        block_num: BLOCK_NUM,
                        variable,
                        var_unit: var_unit.to_string(),
                        value,
                    });
                }
            }
        }
    }
    Ok(ResponseFrame { rows })
}

/// Read and parse a FISPACT-II inventory output file from disk.
pub fn from_file(path: impl AsRef<Path>, run_lbl: &str) -> Result<ResponseFrame> {
    let text = std::fs::read_to_string(path.as_ref())?;
    parse_to_frame(&text, run_lbl)
}

/// Identify a FISPACT-II output by its `.fis` suffix convention.
pub fn is_fispact_output(path: &str) -> bool {
    path.ends_with(".fis")
}

/// Conversion factor to seconds for a FISPACT cooling-time unit.
///
/// Single letters `s m h d w y` plus common long names (case-insensitive),
/// with `w = 7 d` and `y = 365 d`.
fn time_unit_to_seconds(unit: &str) -> Option<f64> {
    match unit.trim().to_ascii_lowercase().as_str() {
        "s" | "sec" | "secs" | "second" | "seconds" => Some(1.0),
        "m" | "min" | "mins" | "minute" | "minutes" => Some(60.0),
        "h" | "hr" | "hrs" | "hour" | "hours" => Some(3_600.0),
        "d" | "day" | "days" => Some(86_400.0),
        "w" | "week" | "weeks" => Some(604_800.0),
        "y" | "yr" | "yrs" | "year" | "years" => Some(31_536_000.0),
        _ => None,
    }
}

/// Parse the tokenised cooling-step header into seconds and a label.
///
/// Uses the *last* `<float> <unit>` pair on the line. A float followed by an
/// unknown alphabetic token errors as an unknown unit; a trailing float with
/// no following token errors as a missing unit.
fn parse_cooling_header(tokens: &[&str], line_no: usize) -> Result<(f64, String)> {
    let mut last: Option<(f64, String, String)> = None;
    for pair in tokens.windows(2) {
        let Ok(value) = pair[0].parse::<f64>() else {
            continue;
        };
        if let Some(factor) = time_unit_to_seconds(pair[1]) {
            last = Some((value * factor, pair[0].to_string(), pair[1].to_string()));
        }
    }
    if let Some((time_s, value_text, unit_raw)) = last {
        return Ok((time_s, format!("{value_text} {unit_raw}")));
    }
    // No valid pair: report the most specific problem for a 1-based error.
    for pair in tokens.windows(2) {
        if pair[0].parse::<f64>().is_ok()
            && pair[1].chars().any(|c| c.is_ascii_alphabetic())
            && time_unit_to_seconds(pair[1]).is_none()
        {
            return Err(Error::Parse {
                line: line_no,
                msg: format!("unknown cooling-time unit `{}`", pair[1]),
            });
        }
    }
    Err(Error::Parse {
        line: line_no,
        msg: "expected cooling time `<value> <unit>`".to_string(),
    })
}

/// Parse one `<nuclide> <atoms> <activity> <heat> [...]` data row.
///
/// Extra trailing columns are ignored. Every name except `total`
/// (case-insensitive) must resolve via [`NuclideId::from_name`].
fn parse_data_row(line: &str, line_no: usize) -> Result<(String, f64, f64, f64)> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 4 {
        return Err(Error::Parse {
            line: line_no,
            msg: format!("expected `nuclide atoms activity heat`, found `{line}`"),
        });
    }
    let nuclide = tokens[0].to_string();
    if !nuclide.eq_ignore_ascii_case("total") {
        NuclideId::from_name(&nuclide).map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("unknown nuclide `{nuclide}`"),
        })?;
    }
    let mut values = Vec::with_capacity(3);
    for (token, what) in [
        (tokens[1], "atoms"),
        (tokens[2], "activity"),
        (tokens[3], "heat"),
    ] {
        let value: f64 = token.parse().map_err(|_| Error::Parse {
            line: line_no,
            msg: format!("expected {what} value, found `{token}`"),
        })?;
        values.push(value);
    }
    Ok((nuclide, values[0], values[1], values[2]))
}

/// True for separator lines holding no alphanumeric characters.
fn is_separator(line: &str) -> bool {
    !line.is_empty()
        && line.chars().all(|c| {
            c.is_whitespace() || matches!(c, '=' | '-' | '*' | '_' | '.' | '+' | '|' | '/' | '\\')
        })
}

/// True for inventory column-title lines (skipped, never data).
fn is_column_header(line: &str) -> bool {
    let lower = line.to_ascii_lowercase();
    lower.starts_with("nuclide")
        || lower.starts_with("isotope")
        || (lower.contains("atoms") && lower.contains("activity"))
}

/// True when a line outside any cooling step looks like a data-row attempt.
fn looks_like_data_row(line: &str) -> bool {
    match line.split_whitespace().next() {
        Some(first) => {
            first.eq_ignore_ascii_case("total") || first.chars().any(|c| c.is_ascii_digit())
        }
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../fixtures/fispact/inventory.fis"
    );

    fn approx(actual: f64, expected: f64) -> bool {
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0)
    }

    #[test]
    fn rejects_empty() {
        assert!(parse_to_frame("  \n", "run1").is_err());
    }

    #[test]
    fn suffix_check() {
        assert!(is_fispact_output("a.fis"));
        assert!(!is_fispact_output("a.out"));
    }

    #[test]
    fn parses_fixture_into_material_inventory_rows() {
        let frame = from_file(FIXTURE, "synth").unwrap();
        // 3 cooling steps x 4 nuclide rows (incl. `total`) x 3 variables.
        assert_eq!(frame.len(), 36);
        assert!(!frame.is_empty());

        assert!(frame.rows.iter().all(|row| row.run_lbl == "synth"));
        assert!(frame
            .rows
            .iter()
            .all(|row| row.block == BlockKind::Material));
        assert!(frame.rows.iter().all(|row| row.block_name == "inventory"));
        assert!(frame.rows.iter().all(|row| row.block_num == -1));
        // No half-life column in the inventory table: everything is unknown.
        assert!(frame.rows.iter().all(|row| row.half_life_s == -1.0));

        let mut times: Vec<f64> = frame.rows.iter().map(|row| row.time_s).collect();
        times.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        times.dedup();
        assert_eq!(times, [0.0, 86_400.0, 31_536_000.0]);

        let variables: std::collections::BTreeSet<String> = frame
            .rows
            .iter()
            .map(|row| row.variable.to_string())
            .collect();
        assert_eq!(
            variables,
            [
                "Number Density".to_string(),
                "Specific Activity".to_string(),
                "Total Decay Heat".to_string(),
            ]
            .into_iter()
            .collect()
        );

        // Per-variable units are fixed by the mapping rules.
        for row in &frame.rows {
            match row.variable {
                ResponseVar::NumberDensity => assert_eq!(row.var_unit, "atoms"),
                ResponseVar::SpecificActivity => assert_eq!(row.var_unit, "Bq"),
                ResponseVar::TotalHeat => assert_eq!(row.var_unit, "W"),
                _ => panic!("unexpected variable {}", row.variable),
            }
        }

        // Spot-check one cooling step: atoms/activity/heat share the row.
        let h3_atoms = frame
            .rows
            .iter()
            .find(|row| {
                row.nuclide == "h-3"
                    && row.variable == ResponseVar::NumberDensity
                    && row.time_s == 86_400.0
            })
            .unwrap();
        assert!(approx(h3_atoms.value, 9.99e19));
        let h3_act = frame
            .rows
            .iter()
            .find(|row| {
                row.nuclide == "h-3"
                    && row.variable == ResponseVar::SpecificActivity
                    && row.time_s == 86_400.0
            })
            .unwrap();
        assert!(approx(h3_act.value, 1.199e9));
        let h3_heat = frame
            .rows
            .iter()
            .find(|row| {
                row.nuclide == "h-3"
                    && row.variable == ResponseVar::TotalHeat
                    && row.time_s == 86_400.0
            })
            .unwrap();
        assert!(approx(h3_heat.value, 3.49e-3));

        // Lowercase-dash names validate; `total` rows do not name a nuclide.
        assert!(frame
            .rows
            .iter()
            .find(|row| row.nuclide == "fe-55")
            .unwrap()
            .nuclide_id()
            .is_ok());
        let totals = frame.totals();
        assert_eq!(totals.len(), 9);
        assert!(totals.rows.iter().all(|row| row.is_total()));
        assert!(totals.rows.iter().all(|row| row.half_life_s == -1.0));
        assert!(frame
            .rows
            .iter()
            .find(|row| row.nuclide == "total")
            .unwrap()
            .nuclide_id()
            .is_err());
    }

    #[test]
    fn bad_float_reports_one_based_line_number() {
        let text = "TIME = 1.0 DAYS\n\
                    h-3 1.0000E+20 1.2000E+09 3.5000E-03\n\
                    fe-55 2.5000E+19 not_a_float 1.1000E-02\n";
        match parse_to_frame(text, "r") {
            Err(Error::Parse { line, .. }) => assert_eq!(line, 3),
            other => panic!("expected parse error, got {other:?}"),
        }
    }

    #[test]
    fn unknown_time_unit_errors_at_header_line() {
        let text = "TIME = 1.0 fortnights\n\
                    h-3 1.0000E+20 1.2000E+09 3.5000E-03\n";
        match parse_to_frame(text, "r") {
            Err(Error::Parse { line, msg }) => {
                assert_eq!(line, 1);
                assert!(msg.contains("fortnights"), "msg was `{msg}`");
            }
            other => panic!("expected parse error, got {other:?}"),
        }
    }

    #[test]
    fn total_rows_skip_nuclide_validation_but_bad_nuclides_fail() {
        let text = "TIME = 0.0 SECS\n\
                    total 1.3000E+20 5.5000E+09 7.6600E-02\n";
        let frame = parse_to_frame(text, "r").unwrap();
        assert_eq!(frame.len(), 3);
        assert!(frame.rows.iter().all(|row| row.is_total()));

        let bad = "TIME = 0.0 SECS\nXx999 1.0 2.0 3.0\n";
        match parse_to_frame(bad, "r") {
            Err(Error::Parse { line, .. }) => assert_eq!(line, 2),
            other => panic!("expected parse error, got {other:?}"),
        }
    }

    #[test]
    fn skips_comments_blanks_separators_and_column_titles() {
        let text = "# leading comment\n\
                    C Fortran-style comment\n\
                    \n\
                    SYNTHETIC INVENTORY TITLE\n\
                    TIME = 12.0 HOURS\n\
                    NUCLIDE ATOMS ACTIVITY HEAT\n\
                    =====\n\
                    h-3 1.0000E+20 1.2000E+09 3.5000E-03\n\
                    \n\
                    C trailing comment\n";
        let frame = parse_to_frame(text, "r").unwrap();
        assert_eq!(frame.len(), 3);
        assert!(frame.rows.iter().all(|row| row.time_s == 43_200.0));
    }

    #[test]
    fn time_units_convert_with_documented_factors() {
        for (header, expected) in [
            ("TIME = 2.0 SECS", 2.0),
            ("TIME = 2.0 MINS", 120.0),
            ("TIME = 2.0 HOURS", 7_200.0),
            ("TIME = 2.0 DAYS", 172_800.0),
            ("TIME = 2.0 WEEKS", 1_209_600.0),
            ("TIME = 2.0 YEARS", 63_072_000.0),
        ] {
            let text = format!("{header}\nh-3 1.0 2.0 3.0\n");
            let frame = parse_to_frame(&text, "r").unwrap();
            assert_eq!(frame.len(), 3, "header `{header}`");
            assert!(
                frame.rows.iter().all(|row| row.time_s == expected),
                "header `{header}`"
            );
        }
    }

    #[test]
    fn data_row_outside_cooling_step_is_an_error() {
        let text = "h-3 1.0 2.0 3.0\n";
        match parse_to_frame(text, "r") {
            Err(Error::Parse { line, .. }) => assert_eq!(line, 1),
            other => panic!("expected parse error, got {other:?}"),
        }
    }
}
