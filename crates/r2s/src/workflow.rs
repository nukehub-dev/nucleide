//! R2S workflow description: zone-to-flux steps derived from an ALARA deck.
//!
//! [`R2sWorkflow::from_deck`] is the primary constructor. The derivation rules
//! are:
//!
//! - **Steps/zones**: one [`R2sStep`] per `mat_loading` entry whose mixture is
//!   not `void` (case-insensitive), in deck order. Void zones carry no
//!   activation source and are skipped.
//! - **Fluxes**: when the deck defines a single `flux` block it is broadcast
//!   to every step; when several exist the step for zone `Z` takes the flux
//!   block named `Z`, and a zone with no name match is a [`Error::CrossRef`].
//! - **Cooling**: copied verbatim from the deck `cooling` block (`times_s`,
//!   seconds); a deck without a `cooling` block yields an empty `cooling_s`,
//!   which [`R2sWorkflow::validate_against`] rejects.
//! - **Top schedule**: the single schedule never referenced as a sub-schedule
//!   (4-token item) by another schedule — mirroring
//!   [`nucleide_alara_io::schedule::expand`]'s top discovery — or the lone schedule
//!   when the deck defines exactly one. Zero or several unreferenced
//!   schedules are a [`Error::CrossRef`].
//!
//! Transport and activation solving stay inside their respective codes; this
//! module only links zones to spectra and expands the irradiation history.

use std::collections::{BTreeSet, HashSet};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// One R2S workflow step linking a mesh zone to an activation case.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct R2sStep {
    /// Zone name from `mat_loading`.
    pub zone: String,
    /// ALARA flux name for the zone.
    pub flux: String,
}

/// Full R2S workflow: per-zone activation steps plus cooling and schedule.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct R2sWorkflow {
    /// Steps in execution order.
    pub steps: Vec<R2sStep>,
    /// After-shutdown cooling times in seconds (from the deck `cooling` block).
    pub cooling_s: Vec<f64>,
    /// Top-level schedule name (see [`R2sWorkflow::from_deck`]).
    pub top_schedule: String,
}

impl R2sWorkflow {
    /// Derive a workflow from an ALARA deck (see module docs for the rules).
    pub fn from_deck(deck: &nucleide_alara_io::deck::AlaraDeck) -> Result<Self> {
        let loading = deck
            .mat_loading
            .as_ref()
            .ok_or_else(|| Error::Invalid("deck defines no `mat_loading` block".to_string()))?;
        let zones: Vec<String> = loading
            .entries
            .iter()
            .filter(|entry| !entry.mixture.eq_ignore_ascii_case("void"))
            .map(|entry| entry.zone.clone())
            .collect();
        if zones.is_empty() {
            return Err(Error::Invalid(
                "deck defines no non-`void` zones in `mat_loading`".to_string(),
            ));
        }
        if deck.fluxes.is_empty() {
            return Err(Error::Invalid("deck defines no `flux` blocks".to_string()));
        }
        let mut steps = Vec::with_capacity(zones.len());
        for zone in &zones {
            steps.push(R2sStep {
                flux: resolve_flux(&deck.fluxes, zone)?,
                zone: zone.clone(),
            });
        }
        let cooling_s = deck
            .cooling
            .as_ref()
            .map(|cooling| cooling.times_s.clone())
            .unwrap_or_default();
        Ok(Self {
            steps,
            cooling_s,
            top_schedule: top_schedule_name(deck)?,
        })
    }

    /// Validate the workflow itself (rejects empty or blank step lists).
    pub fn validate(&self) -> Result<()> {
        if self.steps.is_empty() {
            return Err(Error::Invalid("empty R2S workflow".to_string()));
        }
        for step in &self.steps {
            if step.zone.trim().is_empty() {
                return Err(Error::Invalid("R2S step names an empty zone".to_string()));
            }
            if step.flux.trim().is_empty() {
                return Err(Error::Invalid(format!(
                    "R2S step for zone `{}` names an empty flux",
                    step.zone
                )));
            }
        }
        Ok(())
    }

    /// Cross-check the workflow against a deck: every step zone must appear
    /// in `mat_loading`, every step flux must be a defined `flux` block, and
    /// both the workflow (`cooling_s`) and the deck (`cooling` block) must
    /// carry a non-empty cooling history.
    pub fn validate_against(&self, deck: &nucleide_alara_io::deck::AlaraDeck) -> Result<()> {
        self.validate()?;
        let loading = deck
            .mat_loading
            .as_ref()
            .ok_or_else(|| Error::Invalid("deck defines no `mat_loading` block".to_string()))?;
        let zones: BTreeSet<&str> = loading
            .entries
            .iter()
            .map(|entry| entry.zone.as_str())
            .collect();
        for step in &self.steps {
            if !zones.contains(step.zone.as_str()) {
                return Err(Error::CrossRef(format!(
                    "workflow zone `{}` is not in deck `mat_loading`",
                    step.zone
                )));
            }
        }
        let fluxes: BTreeSet<&str> = deck.fluxes.iter().map(|flux| flux.name.as_str()).collect();
        for step in &self.steps {
            if !fluxes.contains(step.flux.as_str()) {
                return Err(Error::CrossRef(format!(
                    "workflow flux `{}` for zone `{}` is not a defined deck flux",
                    step.flux, step.zone
                )));
            }
        }
        if self.cooling_s.is_empty() {
            return Err(Error::Invalid(
                "workflow defines no cooling times".to_string(),
            ));
        }
        match &deck.cooling {
            Some(cooling) if !cooling.times_s.is_empty() => Ok(()),
            _ => Err(Error::Invalid("deck defines no cooling times".to_string())),
        }
    }

    /// Expand the deck irradiation hierarchy into flat steps.
    ///
    /// Deck `schedule`/`pulsehistory` blocks are converted to the
    /// [`nucleide_alara_io::schedule`] model and combined with `histories` (deck
    /// histories win on name clashes). When [`Self::top_schedule`] is set the
    /// hierarchy rooted there is expanded via
    /// [`nucleide_alara_io::schedule::expand_from`]; otherwise the single unreferenced
    /// top is discovered via [`nucleide_alara_io::schedule::expand`].
    pub fn expand(
        &self,
        deck: &nucleide_alara_io::deck::AlaraDeck,
        histories: &[nucleide_alara_io::schedule::PulseHistory],
    ) -> Result<Vec<nucleide_alara_io::schedule::FlatStep>> {
        let mut schedules = Vec::with_capacity(deck.schedules.len());
        for def in &deck.schedules {
            schedules.push(convert_schedule(def)?);
        }
        let mut combined = Vec::with_capacity(deck.pulse_histories.len() + histories.len());
        for history in &deck.pulse_histories {
            combined.push(convert_history(history));
        }
        let mut known: HashSet<String> = combined
            .iter()
            .map(|history: &nucleide_alara_io::schedule::PulseHistory| history.name.clone())
            .collect();
        for history in histories {
            if known.insert(history.name.clone()) {
                combined.push(history.clone());
            }
        }
        if self.top_schedule.is_empty() {
            nucleide_alara_io::schedule::expand(&schedules, &combined).map_err(map_alara)
        } else {
            nucleide_alara_io::schedule::expand_from(&self.top_schedule, &schedules, &combined)
                .map_err(map_alara)
        }
    }

    /// Emit one deck clone per workflow step with `solve_zones` narrowed to
    /// that step's zone.
    ///
    /// `AlaraDeck` carries both a typed `solve_zones` field and the raw block
    /// list that the canonical writer iterates, so both are updated: the
    /// typed field is replaced and the existing raw `solve_zones` block(s)
    /// get the new body, or a synthetic one (line `0`, marking it as
    /// generated rather than parsed) is appended when the template has none.
    /// Every other block is shared unchanged; per-zone flux or schedule
    /// tailoring stays with the caller.
    pub fn emit_decks(
        &self,
        template: &nucleide_alara_io::deck::AlaraDeck,
    ) -> Result<Vec<nucleide_alara_io::deck::AlaraDeck>> {
        self.validate()?;
        let mut decks = Vec::with_capacity(self.steps.len());
        for step in &self.steps {
            let mut deck = template.clone();
            let line = deck
                .solve_zones
                .as_ref()
                .map(|zones| zones.line)
                .unwrap_or(0);
            deck.solve_zones = Some(nucleide_alara_io::deck::SolveZones {
                zones: vec![step.zone.clone()],
                line,
            });
            let mut found = false;
            for raw in deck
                .blocks
                .iter_mut()
                .filter(|block| block.kind == "solve_zones")
            {
                raw.body = vec![step.zone.clone()];
                found = true;
            }
            if !found {
                deck.blocks.push(nucleide_alara_io::deck::RawBlock {
                    kind: "solve_zones".to_string(),
                    line,
                    body: vec![step.zone.clone()],
                });
            }
            decks.push(deck);
        }
        Ok(decks)
    }
}

/// Pick the flux block for `zone`: broadcast a lone flux, else match by name.
fn resolve_flux(fluxes: &[nucleide_alara_io::deck::FluxDef], zone: &str) -> Result<String> {
    if fluxes.len() == 1 {
        return Ok(fluxes[0].name.clone());
    }
    fluxes
        .iter()
        .find(|flux| flux.name == zone)
        .map(|flux| flux.name.clone())
        .ok_or_else(|| {
            Error::CrossRef(format!(
                "zone `{zone}` matches no flux block ({} flux blocks, no name match)",
                fluxes.len()
            ))
        })
}

/// Top-schedule discovery mirroring [`nucleide_alara_io::schedule::expand`]: the lone
/// schedule wins outright, otherwise the single schedule never referenced as
/// a 4-token sub-schedule item must be unique.
fn top_schedule_name(deck: &nucleide_alara_io::deck::AlaraDeck) -> Result<String> {
    if deck.schedules.is_empty() {
        return Err(Error::Invalid("deck defines no schedules".to_string()));
    }
    if deck.schedules.len() == 1 {
        return Ok(deck.schedules[0].name.clone());
    }
    let referenced: HashSet<&str> = deck
        .schedules
        .iter()
        .flat_map(|schedule| schedule.items.iter())
        .filter(|item| item.tokens.len() == 4)
        .map(|item| item.tokens[0].as_str())
        .collect();
    let tops: Vec<&str> = deck
        .schedules
        .iter()
        .map(|schedule| schedule.name.as_str())
        .filter(|name| !referenced.contains(name))
        .collect();
    match tops.as_slice() {
        [top] => Ok((*top).to_string()),
        [] => Err(Error::CrossRef(
            "no top-level schedule: every schedule is referenced (possible recursion)".to_string(),
        )),
        _ => Err(Error::CrossRef(format!(
            "ambiguous top-level schedule: {}",
            tops.join(", ")
        ))),
    }
}

/// Convert a raw deck schedule (verbatim token items) to the expandable model.
fn convert_schedule(
    def: &nucleide_alara_io::deck::ScheduleDef,
) -> Result<nucleide_alara_io::schedule::ScheduleDef> {
    let mut items = Vec::with_capacity(def.items.len());
    for item in &def.items {
        let tokens = &item.tokens;
        match tokens.as_slice() {
            [op_text, op_unit, flux, history, delay_text, delay_unit] => {
                items.push(nucleide_alara_io::schedule::SchedItem::Pulse {
                    op_time_s: parse_time(op_text, op_unit, item.line)?,
                    flux: flux.clone(),
                    history: history.clone(),
                    delay_s: parse_time(delay_text, delay_unit, item.line)?,
                });
            }
            [name, history, delay_text, delay_unit] => {
                items.push(nucleide_alara_io::schedule::SchedItem::SubSchedule {
                    name: name.clone(),
                    history: history.clone(),
                    delay_s: parse_time(delay_text, delay_unit, item.line)?,
                });
            }
            _ => {
                return Err(Error::Invalid(format!(
                    "schedule `{}` has a malformed item at line {}",
                    def.name, item.line
                )));
            }
        }
    }
    Ok(nucleide_alara_io::schedule::ScheduleDef {
        name: def.name.clone(),
        items,
    })
}

/// Convert a deck pulsing history to the expandable model.
fn convert_history(
    history: &nucleide_alara_io::deck::PulseHistory,
) -> nucleide_alara_io::schedule::PulseHistory {
    nucleide_alara_io::schedule::PulseHistory {
        name: history.name.clone(),
        levels: history
            .levels
            .iter()
            .map(|level| nucleide_alara_io::schedule::PulseLevel {
                count: level.pulses,
                delay_s: level.delay_s,
            })
            .collect(),
    }
}

/// Parse a `<value> <unit>` pair with the schedule time-unit vocabulary.
fn parse_time(value_text: &str, unit: &str, line: usize) -> Result<f64> {
    let value: f64 = value_text.parse().map_err(|_| {
        Error::Invalid(format!(
            "expected operating time at line {line}, found `{value_text}`"
        ))
    })?;
    nucleide_alara_io::schedule::parse_time_to_seconds(value, unit).map_err(map_alara)
}

/// Map ALARA errors onto the R2S error type.
fn map_alara(error: nucleide_alara_io::Error) -> Error {
    match error {
        nucleide_alara_io::Error::Io(inner) => Error::Io(inner),
        nucleide_alara_io::Error::CrossRef(message) => Error::CrossRef(message),
        other => Error::Invalid(other.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal deck: two `mat_loading` zones (one `void`), one flux, one
    /// schedule, one pulsing history, one cooling time.
    const MINIMAL_DECK: &str = "geometry rectangular\n\
         mat_loading\n\
         zone_a mix_a\n\
         zone_void void\n\
         end\n\
         mixture mix_a\n\
         element U235 1.0 1.0\n\
         end\n\
         flux flux_a data/flux 1.0 0 default\n\
         schedule top_sched\n\
         1 d flux_a once 0 s\n\
         end\n\
         pulsehistory once\n\
         1 0 s\n\
         end\n\
         cooling\n\
         1 d\n\
         end\n";

    fn minimal_deck() -> nucleide_alara_io::deck::AlaraDeck {
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(MINIMAL_DECK).unwrap();
        deck.validate().unwrap();
        deck
    }

    #[test]
    fn rejects_empty_workflow() {
        assert!(R2sWorkflow::default().validate().is_err());
    }

    #[test]
    fn from_deck_derives_zones_flux_cooling_and_top() {
        let workflow = R2sWorkflow::from_deck(&minimal_deck()).unwrap();
        assert_eq!(
            workflow.steps,
            [R2sStep {
                zone: "zone_a".to_string(),
                flux: "flux_a".to_string(),
            }]
        );
        assert_eq!(workflow.cooling_s, [86_400.0]);
        assert_eq!(workflow.top_schedule, "top_sched");
        workflow.validate_against(&minimal_deck()).unwrap();
    }

    #[test]
    fn from_deck_skips_void_zones_only() {
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             z1 void\n\
             z2 void\n\
             end\n",
        )
        .unwrap();
        assert!(matches!(
            R2sWorkflow::from_deck(&deck),
            Err(Error::Invalid(_))
        ));
    }

    #[test]
    fn from_deck_matches_multiple_fluxes_by_zone_name() {
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             zone_a mix_a\n\
             zone_b mix_a\n\
             end\n\
             mixture mix_a\n\
             element U235 1.0 1.0\n\
             end\n\
             flux zone_a data/fa 1.0 0 default\n\
             flux zone_b data/fb 1.0 0 default\n\
             schedule s\n\
             1 d zone_a once 0 s\n\
             end\n\
             pulsehistory once\n\
             1 0 s\n\
             end\n\
             cooling\n\
             1 d\n\
             end\n",
        )
        .unwrap();
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        assert_eq!(workflow.steps.len(), 2);
        assert_eq!(workflow.steps[0].flux, "zone_a");
        assert_eq!(workflow.steps[1].flux, "zone_b");
    }

    #[test]
    fn from_deck_rejects_unmatched_flux() {
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             zone_a mix_a\n\
             end\n\
             mixture mix_a\n\
             element U235 1.0 1.0\n\
             end\n\
             flux other_a data/fa 1.0 0 default\n\
             flux other_b data/fb 1.0 0 default\n\
             schedule s\n\
             1 d other_a once 0 s\n\
             end\n\
             pulsehistory once\n\
             1 0 s\n\
             end\n\
             cooling\n\
             1 d\n\
             end\n",
        )
        .unwrap();
        assert!(matches!(
            R2sWorkflow::from_deck(&deck),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn top_schedule_prefers_unreferenced_over_lone_fallback() {
        // Two schedules where `inner` is referenced: `top` wins.
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             zone_a mix_a\n\
             end\n\
             mixture mix_a\n\
             element U235 1.0 1.0\n\
             end\n\
             flux flux_a data/flux 1.0 0 default\n\
             schedule top\n\
             inner once 1 d\n\
             end\n\
             schedule inner\n\
             1 d flux_a once 0 s\n\
             end\n\
             pulsehistory once\n\
             1 0 s\n\
             end\n\
             cooling\n\
             1 d\n\
             end\n",
        )
        .unwrap();
        assert_eq!(R2sWorkflow::from_deck(&deck).unwrap().top_schedule, "top");

        // Two unreferenced schedules: ambiguous.
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             zone_a mix_a\n\
             end\n\
             mixture mix_a\n\
             element U235 1.0 1.0\n\
             end\n\
             flux flux_a data/flux 1.0 0 default\n\
             schedule s1\n\
             1 d flux_a once 0 s\n\
             end\n\
             schedule s2\n\
             1 d flux_a once 0 s\n\
             end\n\
             pulsehistory once\n\
             1 0 s\n\
             end\n\
             cooling\n\
             1 d\n\
             end\n",
        )
        .unwrap();
        assert!(matches!(
            R2sWorkflow::from_deck(&deck),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn validate_against_rejects_unknown_flux() {
        let deck = minimal_deck();
        let workflow = R2sWorkflow {
            steps: vec![R2sStep {
                zone: "zone_a".to_string(),
                flux: "ghost_flux".to_string(),
            }],
            cooling_s: vec![86_400.0],
            top_schedule: "top_sched".to_string(),
        };
        assert!(matches!(
            workflow.validate_against(&deck),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn validate_against_rejects_unknown_zone() {
        let deck = minimal_deck();
        let workflow = R2sWorkflow {
            steps: vec![R2sStep {
                zone: "ghost_zone".to_string(),
                flux: "flux_a".to_string(),
            }],
            cooling_s: vec![86_400.0],
            top_schedule: "top_sched".to_string(),
        };
        assert!(matches!(
            workflow.validate_against(&deck),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn validate_against_rejects_empty_cooling() {
        let deck = minimal_deck();
        // Empty workflow cooling history.
        let workflow = R2sWorkflow {
            steps: vec![R2sStep {
                zone: "zone_a".to_string(),
                flux: "flux_a".to_string(),
            }],
            cooling_s: Vec::new(),
            top_schedule: "top_sched".to_string(),
        };
        assert!(matches!(
            workflow.validate_against(&deck),
            Err(Error::Invalid(_))
        ));

        // Deck without a cooling block.
        let bare = nucleide_alara_io::deck::AlaraDeck::parse(
            "geometry rectangular\n\
             mat_loading\n\
             zone_a mix_a\n\
             end\n\
             mixture mix_a\n\
             element U235 1.0 1.0\n\
             end\n\
             flux flux_a data/flux 1.0 0 default\n\
             schedule s\n\
             1 d flux_a once 0 s\n\
             end\n\
             pulsehistory once\n\
             1 0 s\n\
             end\n",
        )
        .unwrap();
        let workflow = R2sWorkflow {
            steps: vec![R2sStep {
                zone: "zone_a".to_string(),
                flux: "flux_a".to_string(),
            }],
            cooling_s: vec![1.0],
            top_schedule: "s".to_string(),
        };
        assert!(matches!(
            workflow.validate_against(&bare),
            Err(Error::Invalid(_))
        ));
    }

    #[test]
    fn expand_matches_total_time_spot_check() {
        let deck = minimal_deck();
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        let steps = workflow.expand(&deck, &[]).unwrap();
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].duration_s, 86_400.0);
        assert_eq!(steps[0].flux, "flux_a");
        assert_eq!(nucleide_alara_io::schedule::total_time(&steps), 86_400.0);
    }

    #[test]
    fn expand_accepts_extra_histories() {
        let deck = minimal_deck();
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        let extra = nucleide_alara_io::schedule::PulseHistory {
            name: "extra".to_string(),
            levels: vec![nucleide_alara_io::schedule::PulseLevel {
                count: 1,
                delay_s: 0.0,
            }],
        };
        let steps = workflow.expand(&deck, &[extra]).unwrap();
        assert_eq!(nucleide_alara_io::schedule::total_time(&steps), 86_400.0);
    }

    #[test]
    fn emit_decks_sets_solve_zones_per_zone() {
        let deck = minimal_deck();
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        let emitted = workflow.emit_decks(&deck).unwrap();
        assert_eq!(emitted.len(), 1);
        assert_eq!(
            emitted[0].solve_zones.as_ref().unwrap().zones,
            ["zone_a".to_string()]
        );
        // Raw block mirrors the typed field so the canonical writer emits it.
        let bodies: Vec<Vec<String>> = emitted[0]
            .blocks
            .iter()
            .filter(|block| block.kind == "solve_zones")
            .map(|block| block.body.clone())
            .collect();
        assert_eq!(bodies, [vec!["zone_a".to_string()]]);
        assert_eq!(emitted[0].blocks.last().unwrap().kind, "solve_zones");
        // Template untouched.
        assert!(deck.solve_zones.is_none());
        // Round-trips through the canonical writer.
        let text = emitted[0].to_string();
        let reparsed = nucleide_alara_io::deck::AlaraDeck::parse(&text).unwrap();
        assert_eq!(
            reparsed.solve_zones.as_ref().unwrap().zones,
            ["zone_a".to_string()]
        );
    }

    #[test]
    fn emit_decks_replaces_existing_solve_zones() {
        let deck = nucleide_alara_io::deck::AlaraDeck::parse(&format!(
            "{MINIMAL_DECK}solve_zones\nzone_a\nzone_void\nend\n"
        ))
        .unwrap();
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        let emitted = workflow.emit_decks(&deck).unwrap();
        assert_eq!(emitted.len(), 1);
        assert_eq!(
            emitted[0].solve_zones.as_ref().unwrap().zones,
            ["zone_a".to_string()]
        );
        let raw_count = emitted[0]
            .blocks
            .iter()
            .filter(|block| block.kind == "solve_zones")
            .count();
        assert_eq!(raw_count, 1);
    }

    #[test]
    fn emit_decks_rejects_empty_workflow() {
        let deck = minimal_deck();
        assert!(matches!(
            R2sWorkflow::default().emit_decks(&deck),
            Err(Error::Invalid(_))
        ));
    }
}
