//! ALARA irradiation/cooling schedules: flat steps plus the hierarchical
//! schedule / pulse-history form.
//!
//! The flat [`Schedule`] parser covers simple `<value> <unit>` duration lists
//! for downstream depletion stepping. The hierarchical form mirrors the ALARA
//! input deck (`schedule` / `pulsehistory` blocks, see the ALARA users guide
//! `inputtext.rst` schedule/pulsehistory/cooling sections):
//!
//! ```text
//! schedule top_sched
//!     sch_1   once   1 d
//!     30 w    my_flux ph_2 0 s
//! end
//!
//! pulsehistory ph_2
//!     10  5.0  m
//! end
//! ```
//!
//! Each [`ScheduleDef`] item is either a [`SchedItem::Pulse`] (operating time,
//! flux, pulsing history, post-item delay) or a [`SchedItem::SubSchedule`]
//! (sub-schedule name, pulsing history, post-item delay). [`expand`] flattens
//! the hierarchy from its top schedule into [`FlatStep`]s following ALARA's
//! solution semantics: a pulse irradiates, then its pulsing history repeats
//! the whole item (`count` copies per level with the level delay between
//! copies), then the post-item delay elapses; a sub-schedule item expands the
//! referenced schedule first and applies its own history/delay around the
//! whole sub-result. Delay steps carry an empty [`FlatStep::flux`] (flux off).
//!
//! Missing schedule/history references and recursive hierarchies surface as
//! [`Error::CrossRef`].

use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::error::{Error, Result};

/// One irradiation or cooling interval.
///
/// Durations are stored in seconds so downstream consumers (notably the
/// `depletion` crate's steppers, which take `dt: f64` seconds) can use them
/// without further conversion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScheduleItem {
    /// Interval duration in seconds.
    pub duration_s: f64,
    /// Flux scaling factor for the interval (`1.0` = full flux, `0.0` = cooling).
    pub flux_scale: f64,
    /// Optional interval label from the deck.
    pub label: String,
}

/// An ALARA irradiation/cooling schedule.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Schedule {
    /// Intervals in deck order.
    pub items: Vec<ScheduleItem>,
}

impl Schedule {
    /// Parse `<value> <unit> [flux_scale] [label ...]` lines from in-memory text.
    ///
    /// Blank lines and `#` comments are skipped. Empty input is an error.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty schedule".to_string(),
            });
        }
        let mut items = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let line_no = index + 1;
            let words: Vec<&str> = line.split_whitespace().collect();
            let (value_text, unit) = match words.as_slice() {
                [value, unit, ..] => (*value, *unit),
                _ => {
                    return Err(Error::Parse {
                        line: line_no,
                        msg: format!("expected `<value> <unit>`, found `{line}`"),
                    });
                }
            };
            let value: f64 = value_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected duration, found `{value_text}`"),
            })?;
            let duration_s = parse_time_to_seconds(value, unit)
                .map_err(|error| Error::BadUnits(format!("line {line_no}: {error}")))?;
            let flux_scale = match words.get(2) {
                None => 1.0,
                Some(text) => text.parse().map_err(|_| Error::Parse {
                    line: line_no,
                    msg: format!("expected flux scale, found `{text}`"),
                })?,
            };
            let label = words.iter().skip(3).copied().collect::<Vec<_>>().join(" ");
            items.push(ScheduleItem {
                duration_s,
                flux_scale,
                label,
            });
        }
        Ok(Self { items })
    }

    /// Read and parse a schedule from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(&text)
    }

    /// Total schedule time in seconds; feeds `depletion` step `dt` values.
    pub fn total_time(&self) -> f64 {
        self.items.iter().map(|item| item.duration_s).sum()
    }

    /// Number of intervals.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// True when no intervals were parsed.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

/// One item of a hierarchical ALARA schedule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SchedItem {
    /// A simple pulse: operate for `op_time_s` at `flux`, repeat per
    /// `history`, then wait `delay_s` (flux off).
    Pulse {
        /// Operating time in seconds.
        op_time_s: f64,
        /// Symbolic flux name referenced by this pulse.
        flux: String,
        /// Symbolic pulsing-history name governing repetition.
        history: String,
        /// Post-item delay in seconds.
        delay_s: f64,
    },
    /// A sub-schedule reference: expand `name`, repeat per `history`, then
    /// wait `delay_s` (flux off).
    SubSchedule {
        /// Symbolic sub-schedule name.
        name: String,
        /// Symbolic pulsing-history name governing repetition.
        history: String,
        /// Post-item delay in seconds.
        delay_s: f64,
    },
}

/// One hierarchical ALARA schedule definition.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ScheduleDef {
    /// Symbolic schedule name.
    pub name: String,
    /// Items in deck order.
    pub items: Vec<SchedItem>,
}

impl ScheduleDef {
    /// Parse one `schedule <name> ... end` block from in-memory text.
    ///
    /// Six-token items (`<op> <unit> <flux> <history> <delay> <unit>`) parse
    /// as pulses; four-token items (`<schedule> <history> <delay> <unit>`)
    /// parse as sub-schedule references. Blank lines and `#` comments are
    /// skipped; the `end` terminator is required.
    pub fn parse(text: &str) -> Result<Self> {
        let lines = meaningful_lines(text);
        let mut lines = lines.into_iter();
        let (line_no, header) = lines.next().ok_or_else(|| Error::Parse {
            line: 1,
            msg: "empty schedule definition".to_string(),
        })?;
        let words: Vec<&str> = header.split_whitespace().collect();
        let name = match words.as_slice() {
            [keyword, name] if keyword.eq_ignore_ascii_case("schedule") => (*name).to_string(),
            _ => {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `schedule <name>`, found `{header}`"),
                });
            }
        };
        let mut items = Vec::new();
        let mut terminated = false;
        for (line_no, line) in lines {
            if line.eq_ignore_ascii_case("end") {
                terminated = true;
                break;
            }
            items.push(parse_sched_item(line, line_no)?);
        }
        if !terminated {
            return Err(Error::Parse {
                line: line_no_of_last(text),
                msg: format!("schedule `{name}` is missing its `end` terminator"),
            });
        }
        Ok(Self { name, items })
    }
}

/// One pulsing level: repeat the governed block `count` times with `delay_s`
/// (flux off) between repetitions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PulseLevel {
    /// Number of repetitions at this level (at least 1).
    pub count: u64,
    /// Delay between repetitions in seconds.
    pub delay_s: f64,
}

/// One ALARA multi-level pulsing history.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PulseHistory {
    /// Symbolic history name.
    pub name: String,
    /// Pulsing levels in deck order (level 0 applies innermost).
    pub levels: Vec<PulseLevel>,
}

impl PulseHistory {
    /// Parse one `pulsehistory <name> ... end` block from in-memory text.
    ///
    /// Each item is a `<count> <delay> <unit>` triplet. Blank lines and `#`
    /// comments are skipped; the `end` terminator is required.
    pub fn parse(text: &str) -> Result<Self> {
        let lines = meaningful_lines(text);
        let mut lines = lines.into_iter();
        let (line_no, header) = lines.next().ok_or_else(|| Error::Parse {
            line: 1,
            msg: "empty pulse history".to_string(),
        })?;
        let words: Vec<&str> = header.split_whitespace().collect();
        let name = match words.as_slice() {
            [keyword, name] if keyword.eq_ignore_ascii_case("pulsehistory") => (*name).to_string(),
            _ => {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `pulsehistory <name>`, found `{header}`"),
                });
            }
        };
        let mut levels = Vec::new();
        let mut terminated = false;
        for (line_no, line) in lines {
            if line.eq_ignore_ascii_case("end") {
                terminated = true;
                break;
            }
            let words: Vec<&str> = line.split_whitespace().collect();
            let [count_text, delay_text, unit] = words.as_slice() else {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `<count> <delay> <unit>`, found `{line}`"),
                });
            };
            let count: u64 = count_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected pulse count, found `{count_text}`"),
            })?;
            if count == 0 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: "pulse count must be at least 1".to_string(),
                });
            }
            let delay: f64 = delay_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected delay, found `{delay_text}`"),
            })?;
            let delay_s = parse_time_to_seconds(delay, unit)
                .map_err(|error| Error::BadUnits(format!("line {line_no}: {error}")))?;
            levels.push(PulseLevel { count, delay_s });
        }
        if !terminated {
            return Err(Error::Parse {
                line: line_no_of_last(text),
                msg: format!("pulse history `{name}` is missing its `end` terminator"),
            });
        }
        Ok(Self { name, levels })
    }
}

/// One flattened irradiation/cooling step.
///
/// Irradiation steps carry their symbolic flux name; delay/cooling steps carry
/// an empty `flux` (see [`FlatStep::is_cooling`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FlatStep {
    /// Step duration in seconds.
    pub duration_s: f64,
    /// Symbolic flux name, or `""` for flux-off delay steps.
    pub flux: String,
}

impl FlatStep {
    /// True for flux-off delay/cooling steps.
    pub fn is_cooling(&self) -> bool {
        self.flux.is_empty()
    }
}

/// Total time in seconds over flattened [`FlatStep`]s.
pub fn total_time(steps: &[FlatStep]) -> f64 {
    steps.iter().map(|step| step.duration_s).sum()
}

/// Expand the top-level schedule of a hierarchy into flat steps.
///
/// The top schedule is the one never referenced as a sub-schedule; an empty
/// set, zero candidates (everything referenced: possible recursion), or
/// several candidates are all errors. See [`expand_from`] for an explicit top.
pub fn expand(schedules: &[ScheduleDef], histories: &[PulseHistory]) -> Result<Vec<FlatStep>> {
    if schedules.is_empty() {
        return Err(Error::Parse {
            line: 1,
            msg: "no schedules to expand".to_string(),
        });
    }
    let referenced: std::collections::HashSet<&str> = schedules
        .iter()
        .flat_map(|sched| sched.items.iter())
        .filter_map(|item| match item {
            SchedItem::SubSchedule { name, .. } => Some(name.as_str()),
            SchedItem::Pulse { .. } => None,
        })
        .collect();
    let tops: Vec<&str> = schedules
        .iter()
        .map(|sched| sched.name.as_str())
        .filter(|name| !referenced.contains(name))
        .collect();
    match tops.as_slice() {
        [top] => expand_from(top, schedules, histories),
        [] => Err(Error::CrossRef(
            "no top-level schedule: every schedule is referenced (possible recursion)".to_string(),
        )),
        _ => Err(Error::CrossRef(format!(
            "ambiguous top-level schedule: {}",
            tops.join(", ")
        ))),
    }
}

/// Expand the hierarchy rooted at `top` into flat steps.
///
/// Unknown schedule/history references and recursive hierarchies surface as
/// [`Error::CrossRef`].
pub fn expand_from(
    top: &str,
    schedules: &[ScheduleDef],
    histories: &[PulseHistory],
) -> Result<Vec<FlatStep>> {
    let mut sched_map: std::collections::HashMap<&str, &ScheduleDef> =
        std::collections::HashMap::new();
    for sched in schedules {
        sched_map.entry(sched.name.as_str()).or_insert(sched);
    }
    let mut hist_map: std::collections::HashMap<&str, &PulseHistory> =
        std::collections::HashMap::new();
    for hist in histories {
        hist_map.entry(hist.name.as_str()).or_insert(hist);
    }
    if !sched_map.contains_key(top) {
        return Err(Error::CrossRef(format!("unknown schedule `{top}`")));
    }
    let mut stack: Vec<String> = Vec::new();
    expand_schedule(top, &sched_map, &hist_map, &mut stack)
}

/// Convert a time value to seconds.
///
/// Accepts case-insensitive `s`/`sec`/`second`(+`s`), `m`/`min`/`minute`(+`s`),
/// `h`/`hr`/`hour`(+`s`), `d`/`day`(+`s`), `w`/`week`(+`s`),
/// `y`/`yr`/`year`(+`s`) with one year defined as 365.25 days, and
/// `c`/`century`(+`ies`) with one century defined as 100 years. Rejects
/// non-finite or negative values and unknown units with [`Error::BadUnits`].
pub fn parse_time_to_seconds(value: f64, unit: &str) -> Result<f64> {
    if !value.is_finite() || value < 0.0 {
        return Err(Error::BadUnits(format!(
            "negative or non-finite time `{value}`"
        )));
    }
    let year = 365.25 * 86_400.0;
    let factor = match unit.trim().to_ascii_lowercase().as_str() {
        "s" | "sec" | "secs" | "second" | "seconds" => 1.0,
        "m" | "min" | "mins" | "minute" | "minutes" => 60.0,
        "h" | "hr" | "hrs" | "hour" | "hours" => 3_600.0,
        "d" | "day" | "days" => 86_400.0,
        "w" | "week" | "weeks" => 604_800.0,
        "y" | "yr" | "yrs" | "year" | "years" => year,
        "c" | "century" | "centuries" => 100.0 * year,
        _ => return Err(Error::BadUnits(format!("unknown time unit `{unit}`"))),
    };
    Ok(value * factor)
}

/// Non-blank, non-`#`-comment lines with 1-based line numbers.
fn meaningful_lines(text: &str) -> Vec<(usize, &str)> {
    text.lines()
        .enumerate()
        .map(|(index, raw)| (index + 1, raw.trim()))
        .filter(|(_, line)| !line.is_empty() && !line.starts_with('#'))
        .collect()
}

/// 1-based line number of the last line (for unterminated-block errors).
fn line_no_of_last(text: &str) -> usize {
    text.lines().count().max(1)
}

/// Parse one schedule item: six tokens = pulse, four tokens = sub-schedule.
fn parse_sched_item(line: &str, line_no: usize) -> Result<SchedItem> {
    let words: Vec<&str> = line.split_whitespace().collect();
    match words.as_slice() {
        [op_text, op_unit, flux, history, delay_text, delay_unit]
            if op_text.parse::<f64>().is_ok() =>
        {
            let op: f64 = op_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected operating time, found `{op_text}`"),
            })?;
            let op_time_s = parse_time_to_seconds(op, op_unit)
                .map_err(|error| Error::BadUnits(format!("line {line_no}: {error}")))?;
            let delay: f64 = delay_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected delay, found `{delay_text}`"),
            })?;
            let delay_s = parse_time_to_seconds(delay, delay_unit)
                .map_err(|error| Error::BadUnits(format!("line {line_no}: {error}")))?;
            Ok(SchedItem::Pulse {
                op_time_s,
                flux: (*flux).to_string(),
                history: (*history).to_string(),
                delay_s,
            })
        }
        [name, history, delay_text, delay_unit] if name.parse::<f64>().is_err() => {
            let delay: f64 = delay_text.parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected delay, found `{delay_text}`"),
            })?;
            let delay_s = parse_time_to_seconds(delay, delay_unit)
                .map_err(|error| Error::BadUnits(format!("line {line_no}: {error}")))?;
            Ok(SchedItem::SubSchedule {
                name: (*name).to_string(),
                history: (*history).to_string(),
                delay_s,
            })
        }
        _ => Err(Error::Parse {
            line: line_no,
            msg: format!(
                "expected `<op> <unit> <flux> <history> <delay> <unit>` or \
                 `<schedule> <history> <delay> <unit>`, found `{line}`"
            ),
        }),
    }
}

/// Recursively expand one schedule; `stack` holds the active reference chain
/// for cycle detection.
fn expand_schedule(
    name: &str,
    sched_map: &std::collections::HashMap<&str, &ScheduleDef>,
    hist_map: &std::collections::HashMap<&str, &PulseHistory>,
    stack: &mut Vec<String>,
) -> Result<Vec<FlatStep>> {
    let Some(def) = sched_map.get(name).copied() else {
        let parent = stack.last().map(String::as_str).unwrap_or("<root>");
        return Err(Error::CrossRef(format!(
            "schedule `{parent}` references unknown sub-schedule `{name}`"
        )));
    };
    if stack.iter().any(|entry| entry == name) {
        return Err(Error::CrossRef(format!(
            "recursive schedule hierarchy involving `{name}`"
        )));
    }
    stack.push(name.to_string());
    let mut steps = Vec::new();
    for item in &def.items {
        match item {
            SchedItem::Pulse {
                op_time_s,
                flux,
                history,
                delay_s,
            } => {
                if *op_time_s < 0.0 {
                    return Err(Error::Parse {
                        line: 1,
                        msg: format!("schedule `{name}` has negative operating time"),
                    });
                }
                steps.push(FlatStep {
                    duration_s: *op_time_s,
                    flux: flux.clone(),
                });
                apply_history(
                    &mut steps,
                    history,
                    &format!("pulse in schedule `{name}`"),
                    hist_map,
                )?;
                push_delay(&mut steps, *delay_s)?;
            }
            SchedItem::SubSchedule {
                name: sub,
                history,
                delay_s,
            } => {
                let mut sub_steps = expand_schedule(sub, sched_map, hist_map, stack)?;
                apply_history(
                    &mut sub_steps,
                    history,
                    &format!("sub-schedule `{sub}` in schedule `{name}`"),
                    hist_map,
                )?;
                steps.append(&mut sub_steps);
                push_delay(&mut steps, *delay_s)?;
            }
        }
    }
    stack.pop();
    Ok(steps)
}

/// Repeat `steps` per the named pulsing history (level 0 innermost).
fn apply_history(
    steps: &mut Vec<FlatStep>,
    history: &str,
    context: &str,
    hist_map: &std::collections::HashMap<&str, &PulseHistory>,
) -> Result<()> {
    let Some(hist) = hist_map.get(history).copied() else {
        return Err(Error::CrossRef(format!(
            "{context} references unknown pulse history `{history}`"
        )));
    };
    for level in &hist.levels {
        if level.delay_s < 0.0 {
            return Err(Error::Parse {
                line: 1,
                msg: format!("pulse history `{history}` has negative delay"),
            });
        }
        match level.count {
            0 => {
                return Err(Error::Parse {
                    line: 1,
                    msg: format!("pulse history `{history}` has a level with zero pulses"),
                });
            }
            1 => {}
            count => {
                let base = steps.clone();
                for _ in 1..count {
                    if level.delay_s > 0.0 {
                        steps.push(FlatStep {
                            duration_s: level.delay_s,
                            flux: String::new(),
                        });
                    }
                    steps.extend(base.iter().cloned());
                }
            }
        }
    }
    Ok(())
}

/// Append a flux-off delay step unless the delay is zero.
fn push_delay(steps: &mut Vec<FlatStep>, delay_s: f64) -> Result<()> {
    if delay_s < 0.0 {
        return Err(Error::Parse {
            line: 1,
            msg: format!("negative schedule delay `{delay_s}`"),
        });
    }
    if delay_s > 0.0 {
        steps.push(FlatStep {
            duration_s: delay_s,
            flux: String::new(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_units_convert_to_seconds() {
        assert_eq!(parse_time_to_seconds(1.0, "s").unwrap(), 1.0);
        assert_eq!(parse_time_to_seconds(2.0, "min").unwrap(), 120.0);
        assert_eq!(parse_time_to_seconds(2.0, "h").unwrap(), 7_200.0);
        assert_eq!(parse_time_to_seconds(1.0, "D").unwrap(), 86_400.0);
        assert_eq!(
            parse_time_to_seconds(1.0, "years").unwrap(),
            365.25 * 86_400.0
        );
    }

    #[test]
    fn week_and_century_units_convert_to_seconds() {
        assert_eq!(parse_time_to_seconds(1.0, "w").unwrap(), 604_800.0);
        assert_eq!(
            parse_time_to_seconds(1.0, "c").unwrap(),
            100.0 * 365.25 * 86_400.0
        );
    }

    #[test]
    fn bad_time_units_error() {
        assert!(matches!(
            parse_time_to_seconds(1.0, "fortnight"),
            Err(Error::BadUnits(_))
        ));
        assert!(matches!(
            parse_time_to_seconds(-1.0, "s"),
            Err(Error::BadUnits(_))
        ));
        assert!(matches!(
            parse_time_to_seconds(f64::NAN, "s"),
            Err(Error::BadUnits(_))
        ));
    }

    #[test]
    fn parse_builds_intervals_in_order() {
        let schedule = Schedule::parse("1.0 h 1.0 burn\n30.0 min 0.0 cool down\n").unwrap();
        assert_eq!(schedule.len(), 2);
        assert!(!schedule.is_empty());
        assert_eq!(schedule.items[0].duration_s, 3_600.0);
        assert_eq!(schedule.items[1].flux_scale, 0.0);
        assert_eq!(schedule.items[1].label, "cool down");
        assert_eq!(schedule.total_time(), 5_400.0);
    }

    #[test]
    fn parse_defaults_flux_scale_and_rejects_bad_lines() {
        let schedule = Schedule::parse("10 s\n").unwrap();
        assert_eq!(schedule.items[0].flux_scale, 1.0);
        assert_eq!(schedule.items[0].label, "");
        assert!(matches!(Schedule::parse("   \n"), Err(Error::Parse { .. })));
        assert!(matches!(
            Schedule::parse("10 fortnights\n"),
            Err(Error::BadUnits(_))
        ));
        assert!(matches!(
            Schedule::parse("soon h\n"),
            Err(Error::Parse { .. })
        ));
    }

    fn histories(names: &[(&str, Vec<(u64, f64)>)]) -> Vec<PulseHistory> {
        names
            .iter()
            .map(|(name, levels)| PulseHistory {
                name: (*name).to_string(),
                levels: levels
                    .iter()
                    .map(|(count, delay_s)| PulseLevel {
                        count: *count,
                        delay_s: *delay_s,
                    })
                    .collect(),
            })
            .collect()
    }

    #[test]
    fn parse_schedule_def_with_pulse_item() {
        let def = ScheduleDef::parse(
            "schedule 1_year\n\
             \t1 y  flux_1  steady_state  0 s\n\
             end\n",
        )
        .unwrap();
        assert_eq!(def.name, "1_year");
        assert_eq!(
            def.items,
            [SchedItem::Pulse {
                op_time_s: 365.25 * 86_400.0,
                flux: "flux_1".to_string(),
                history: "steady_state".to_string(),
                delay_s: 0.0,
            }]
        );
        assert!(matches!(
            ScheduleDef::parse("schedule 1_year\n\t1 y flux_1 steady_state 0 s\n"),
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            ScheduleDef::parse("bogus 1_year\nend\n"),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn parse_schedule_def_with_sub_schedule_items() {
        let def = ScheduleDef::parse(
            "schedule top_sched\n\
                 sch_1   once   1 d\n\
                 sch_2   once   25 d\n\
             end\n",
        )
        .unwrap();
        assert_eq!(def.name, "top_sched");
        assert_eq!(
            def.items,
            [
                SchedItem::SubSchedule {
                    name: "sch_1".to_string(),
                    history: "once".to_string(),
                    delay_s: 86_400.0,
                },
                SchedItem::SubSchedule {
                    name: "sch_2".to_string(),
                    history: "once".to_string(),
                    delay_s: 25.0 * 86_400.0,
                },
            ]
        );
    }

    #[test]
    fn parse_pulse_history_levels() {
        let single = PulseHistory::parse("pulsehistory steady_state\n\t1\t0 s\nend\n").unwrap();
        assert_eq!(single.name, "steady_state");
        assert_eq!(
            single.levels,
            [PulseLevel {
                count: 1,
                delay_s: 0.0
            }]
        );
        let multi = PulseHistory::parse(
            "pulsehistory ph_1\n\
                 1   0.0    s\n\
                 5   5.0    m\n\
                 10  1.0    d\n\
             end\n",
        )
        .unwrap();
        assert_eq!(multi.levels.len(), 3);
        assert_eq!(multi.levels[1].count, 5);
        assert_eq!(multi.levels[1].delay_s, 300.0);
        assert_eq!(multi.levels[2].delay_s, 86_400.0);
        assert!(matches!(
            PulseHistory::parse("pulsehistory ph\n\t0 1 s\nend\n"),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn expand_repeats_pulse_per_history_then_delay() {
        let schedules = [ScheduleDef {
            name: "top".to_string(),
            items: vec![SchedItem::Pulse {
                op_time_s: 100.0,
                flux: "flux_1".to_string(),
                history: "h".to_string(),
                delay_s: 50.0,
            }],
        }];
        let histories = histories(&[("h", vec![(3, 10.0)])]);
        let steps = expand_from("top", &schedules, &histories).unwrap();
        let durations: Vec<f64> = steps.iter().map(|step| step.duration_s).collect();
        assert_eq!(durations, [100.0, 10.0, 100.0, 10.0, 100.0, 50.0]);
        assert_eq!(steps[0].flux, "flux_1");
        assert!(steps[1].is_cooling());
        assert!(!steps[0].is_cooling());
        assert_eq!(total_time(&steps), 370.0);
    }

    #[test]
    fn expand_nests_sub_schedules_with_outer_history() {
        let schedules = [
            ScheduleDef {
                name: "top".to_string(),
                items: vec![SchedItem::SubSchedule {
                    name: "inner".to_string(),
                    history: "twice".to_string(),
                    delay_s: 7.0,
                }],
            },
            ScheduleDef {
                name: "inner".to_string(),
                items: vec![SchedItem::Pulse {
                    op_time_s: 20.0,
                    flux: "f".to_string(),
                    history: "once".to_string(),
                    delay_s: 0.0,
                }],
            },
        ];
        let histories = histories(&[("once", vec![(1, 0.0)]), ("twice", vec![(2, 5.0)])]);
        // `expand` finds the single unreferenced top on its own.
        let steps = expand(&schedules, &histories).unwrap();
        let durations: Vec<f64> = steps.iter().map(|step| step.duration_s).collect();
        assert_eq!(durations, [20.0, 5.0, 20.0, 7.0]);
        assert_eq!(total_time(&steps), 52.0);
    }

    #[test]
    fn expand_detects_cycles() {
        let cyclic = [
            ScheduleDef {
                name: "a".to_string(),
                items: vec![SchedItem::SubSchedule {
                    name: "b".to_string(),
                    history: "once".to_string(),
                    delay_s: 0.0,
                }],
            },
            ScheduleDef {
                name: "b".to_string(),
                items: vec![SchedItem::SubSchedule {
                    name: "a".to_string(),
                    history: "once".to_string(),
                    delay_s: 0.0,
                }],
            },
        ];
        let histories = histories(&[("once", vec![(1, 0.0)])]);
        assert!(matches!(
            expand_from("a", &cyclic, &histories),
            Err(Error::CrossRef(_))
        ));
        // Fully referenced hierarchies have no top candidate either.
        assert!(matches!(
            expand(&cyclic, &histories),
            Err(Error::CrossRef(_))
        ));

        let self_cycle = [ScheduleDef {
            name: "a".to_string(),
            items: vec![SchedItem::SubSchedule {
                name: "a".to_string(),
                history: "once".to_string(),
                delay_s: 0.0,
            }],
        }];
        assert!(matches!(
            expand_from("a", &self_cycle, &histories),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn expand_reports_missing_references() {
        let schedules = [ScheduleDef {
            name: "top".to_string(),
            items: vec![SchedItem::SubSchedule {
                name: "ghost".to_string(),
                history: "once".to_string(),
                delay_s: 0.0,
            }],
        }];
        let histories = histories(&[("once", vec![(1, 0.0)])]);
        assert!(matches!(
            expand_from("top", &schedules, &histories),
            Err(Error::CrossRef(_))
        ));
        assert!(matches!(
            expand_from("ghost", &schedules, &histories),
            Err(Error::CrossRef(_))
        ));

        let pulsed = [ScheduleDef {
            name: "top".to_string(),
            items: vec![SchedItem::Pulse {
                op_time_s: 1.0,
                flux: "f".to_string(),
                history: "ghost_hist".to_string(),
                delay_s: 0.0,
            }],
        }];
        assert!(matches!(
            expand_from("top", &pulsed, &histories),
            Err(Error::CrossRef(_))
        ));
        assert!(matches!(expand(&[], &histories), Err(Error::Parse { .. })));
    }
}
