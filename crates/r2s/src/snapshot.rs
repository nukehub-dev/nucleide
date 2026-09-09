//! ARMI database-snapshot → R2S adapter (dict-in, versionless).
//!
//! The caller runs their own version-adapted dump script (their `h5py` plus
//! their ARMI checkout) and passes plain versionless dicts across the
//! boundary. This crate never touches HDF5, never mirrors ARMI `Layout`
//! versioning, and never sees `layout/`, `gridIndex`, `locationType`,
//! `databaseVersion`, serial numbers, or `S`-byte dtypes: if a future ARMI
//! release changes location packing, only the caller's dump script changes,
//! never this schema.
//!
//! # Boundary rules
//!
//! - **Zone ids are opaque strings** (default 1:1 from ARMI block names such
//!   as `B0120-005`). Merging blocks into homogenized zones is caller-side;
//!   nothing is ever inferred here.
//! - **Volumes method only.** The generated deck carries a `volumes` block
//!   and no `dimension` block (the two together are a parse error). The
//!   `geometry` string is still required by parser consumers, so a constant
//!   `rectangular` block is emitted.
//! - **MatLoading is 1:1 zone → mixture.** Each non-empty zone gets its own
//!   mixture (`mix_<zone>`) built from the composition via
//!   [`MixtureEntry::Element`](nucleide_alara_io::deck::MixtureEntry) with
//!   `rel_density = 1.0`. Empty compositions map to the `void` mixture (no
//!   mixture block) and are skipped by
//!   [`R2sWorkflow::from_deck`](crate::workflow::R2sWorkflow::from_deck),
//!   exactly like parsed `void` zones.
//! - **Composition units.** Values arrive as number densities in
//!   atoms/barn-cm and are stored verbatim as the element `vol_fraction`.
//!   That matches [`AlaraDeck::mixture_ids`](nucleide_alara_io::deck::AlaraDeck::mixture_ids)
//!   semantics, which returns the `vol_fraction` unchanged, so a
//!   snapshot-built deck round-trips its number densities through
//!   `mixture_ids`. Note this is a snapshot-path convention: a hand-written
//!   ALARA deck would resolve the same entries through the element library.
//! - **Key rules mirror the `nucleide-emit` ARMI-input rule** (elemental keys
//!   rejected with "expand first", bare `AM242` rejected, explicit `AM242G`
//!   ground accepted); the caller passes post-expansion nuclide keys, never
//!   raw elemental fractions.
//! - **Fluxes are caller-supplied.** One flux block broadcasts to every step,
//!   several match by zone name (see `resolve_flux` in
//!   [`crate::workflow`]); a per-zone `flux` override may additionally pin a
//!   step to a defined flux block.
//! - **Cooling and schedule are caller-supplied.** An empty `cooling_s`
//!   builds a deck with no `cooling` block, which
//!   [`validate_against`](crate::workflow::R2sWorkflow::validate_against)
//!   rejects; `schedule_text` may carry `schedule`/`pulsehistory` blocks only.
//! - **`SolveZones`/`SkipZones`/`SpatialNorm` stay unset**;
//!   [`emit_decks`](crate::workflow::R2sWorkflow::emit_decks) synthesizes
//!   per-zone `solve_zones`.
//!
//! # Follow-up
//!
//! The `wasm` bindings re-implement `from_deck` without depending on `r2s`
//! and do not expose this adapter yet; porting it there is a follow-up, not
//! part of this change.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use nucleide_alara_io::deck::{
    AlaraDeck, Cooling, FluxDef, Geometry, MatLoading, MatLoadingEntry, Mixture, MixtureEntry,
    OutputDef, PulseHistory, PulseLevel, RawBlock, ScheduleDef, ScheduleItemRef, VolumeEntry,
    Volumes,
};
use nucleide_nuclei::{armi::armi_name_to_nucid, NuclideId};

use crate::error::{Error, Result};
use crate::workflow::R2sWorkflow;

/// Geometry string emitted for snapshot decks (volumes method still needs a
/// `geometry` block for parser consumers; the value carries no mesh axes).
pub const SNAPSHOT_GEOMETRY: &str = "rectangular";

/// Default schedule/history names synthesized when no `schedule_text` is given.
pub const SNAPSHOT_SCHEDULE: &str = "snap_schedule";
/// Default pulsing-history name synthesized when no `schedule_text` is given.
pub const SNAPSHOT_HISTORY: &str = "snap_once";

/// One snapshot zone: an opaque id plus a homogenized composition.
///
/// `zbottom_cm`/`ztop_cm`, `material`, `xs_type`, and `temperature_c` are
/// informational metadata carried for the caller's round-trip; the ALARA deck
/// has no per-zone metadata fields, so they are not emitted. `flux_name`, when
/// set, pins the workflow step to that defined flux block (it must exist).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SnapshotZone {
    /// Opaque zone id, 1:1 from the ARMI block name (e.g. `B0120-005`).
    pub zone: String,
    /// Homogenized zone volume in cm³.
    pub volume_cm3: f64,
    /// Axial bottom in cm, if dumped.
    pub zbottom_cm: Option<f64>,
    /// Axial top in cm, if dumped.
    pub ztop_cm: Option<f64>,
    /// ARMI material name, if dumped (informational only).
    pub material: Option<String>,
    /// Cross-section type label, if dumped (informational only).
    pub xs_type: Option<String>,
    /// Temperature in °C, if dumped (informational only).
    pub temperature_c: Option<f64>,
    /// ARMI bare-name → number density in atoms/barn-cm. Empty means void.
    pub composition: Vec<(String, f64)>,
    /// Optional pin to a defined flux block name.
    pub flux_name: Option<String>,
}

/// One caller-supplied flux definition: `<name> <file> <scale> 0 default`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SnapshotFluxDef {
    /// Symbolic flux name referenced by schedule items.
    pub name: String,
    /// Flux file path.
    pub file: String,
    /// Uniform flux scaling factor.
    pub scale: f64,
}

/// Versionless snapshot input: plain dicts from the caller's dump script.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SnapshotInput {
    /// Zones in deck order.
    pub zones: Vec<SnapshotZone>,
    /// Flux definitions (≥1 required).
    pub flux_defs: Vec<SnapshotFluxDef>,
    /// After-shutdown cooling times in seconds (empty fails workflow validation).
    pub cooling_s: Vec<f64>,
    /// Caller-supplied irradiation-history blocks (`schedule`/`pulsehistory`
    /// only); `None`/blank synthesizes a 1-day-per-flux default.
    pub schedule_text: Option<String>,
    /// Optional `output` resolution (`interval`, `zone`, or `mixture`).
    pub output: Option<String>,
}

impl R2sWorkflow {
    /// Derive a workflow from a snapshot-built deck (delegates to
    /// [`from_deck`](R2sWorkflow::from_deck)).
    pub fn from_snapshot(deck: &AlaraDeck) -> Result<Self> {
        Self::from_deck(deck)
    }
}

/// Build a validated [`AlaraDeck`] template from a versionless snapshot.
///
/// Runs [`AlaraDeck::validate`] before returning, so dangling mixture,
/// volume, or schedule references surface here as [`Error::CrossRef`].
/// A missing `cooling` block (empty `cooling_s`) is *not* a deck error; it
/// fails later in `validate_against`, mirroring parsed decks.
pub fn deck_from_snapshot(input: &SnapshotInput) -> Result<AlaraDeck> {
    if input.zones.is_empty() {
        return Err(Error::Invalid("snapshot defines no zones".to_string()));
    }
    if input.flux_defs.is_empty() {
        return Err(Error::Invalid(
            "snapshot defines no flux definitions".to_string(),
        ));
    }
    let mut seen_zones = BTreeSet::new();
    for zone in &input.zones {
        check_zone_id(&zone.zone)?;
        if !seen_zones.insert(zone.zone.as_str()) {
            return Err(Error::Invalid(format!(
                "snapshot zone `{}` appears more than once",
                zone.zone
            )));
        }
        check_volume(&zone.zone, zone.volume_cm3)?;
        check_extent(zone)?;
    }
    let mut seen_fluxes = BTreeSet::new();
    for flux in &input.flux_defs {
        check_flux_name(&flux.name)?;
        if !seen_fluxes.insert(flux.name.as_str()) {
            return Err(Error::Invalid(format!(
                "snapshot flux `{}` appears more than once",
                flux.name
            )));
        }
        if flux.file.trim().is_empty() {
            return Err(Error::Invalid(format!(
                "snapshot flux `{}` names an empty file",
                flux.name
            )));
        }
        if !flux.scale.is_finite() {
            return Err(Error::Invalid(format!(
                "snapshot flux `{}` has non-finite scale {}",
                flux.name, flux.scale
            )));
        }
    }
    for time in &input.cooling_s {
        if !time.is_finite() || *time < 0.0 {
            return Err(Error::Invalid(format!(
                "snapshot cooling time `{time}` is not a non-negative finite value"
            )));
        }
    }

    let mut deck = AlaraDeck::default();
    push_block(&mut deck, "geometry", vec![SNAPSHOT_GEOMETRY.to_string()]);
    deck.geometry = Some(Geometry {
        kind: SNAPSHOT_GEOMETRY.to_string(),
        line: 0,
    });

    let mut volume_entries = Vec::with_capacity(input.zones.len());
    let mut loading_entries = Vec::with_capacity(input.zones.len());
    for zone in &input.zones {
        volume_entries.push(VolumeEntry {
            volume: zone.volume_cm3,
            zone: zone.zone.clone(),
        });
        let resolved = resolve_composition(&zone.zone, &zone.composition)?;
        if resolved.is_empty() {
            loading_entries.push(MatLoadingEntry {
                zone: zone.zone.clone(),
                mixture: "void".to_string(),
            });
        } else {
            let mixture = mixture_name(&zone.zone);
            let entries = resolved
                .into_iter()
                .map(|(id, ndens)| MixtureEntry::Element {
                    symbol: id.to_name(),
                    rel_density: 1.0,
                    vol_fraction: ndens,
                })
                .collect::<Vec<_>>();
            push_block(
                &mut deck,
                "mixture",
                entries
                    .iter()
                    .map(|entry| match entry {
                        MixtureEntry::Element {
                            symbol,
                            rel_density,
                            vol_fraction,
                        } => format!("element {symbol} {rel_density} {vol_fraction}"),
                        _ => unreachable!("snapshot mixtures only hold element entries"),
                    })
                    .collect(),
            );
            deck.mixtures.push(Mixture {
                name: mixture.clone(),
                entries,
                line: 0,
            });
            loading_entries.push(MatLoadingEntry {
                zone: zone.zone.clone(),
                mixture,
            });
        }
    }
    push_block(
        &mut deck,
        "volumes",
        volume_entries
            .iter()
            .map(|entry| format!("{} {}", entry.volume, entry.zone))
            .collect(),
    );
    deck.volumes = Some(Volumes {
        entries: volume_entries,
        line: 0,
    });
    push_block(
        &mut deck,
        "mat_loading",
        loading_entries
            .iter()
            .map(|entry| format!("{} {}", entry.zone, entry.mixture))
            .collect(),
    );
    deck.mat_loading = Some(MatLoading {
        entries: loading_entries,
        line: 0,
    });

    for flux in &input.flux_defs {
        push_block(
            &mut deck,
            "flux",
            vec![format!(
                "{} {} {} 0 default",
                flux.name, flux.file, flux.scale
            )],
        );
        deck.fluxes.push(FluxDef {
            name: flux.name.clone(),
            file: flux.file.clone(),
            scale: flux.scale,
            skip: 0,
            format: "default".to_string(),
            line: 0,
        });
    }

    match input.schedule_text.as_deref() {
        Some(text) if !text.trim().is_empty() => merge_schedule_text(&mut deck, text)?,
        _ => synthesize_schedule(&mut deck, input)?,
    }

    if !input.cooling_s.is_empty() {
        push_block(
            &mut deck,
            "cooling",
            input
                .cooling_s
                .iter()
                .map(|time| format!("{time} s"))
                .collect(),
        );
        deck.cooling = Some(Cooling {
            times_s: input.cooling_s.clone(),
            line: 0,
        });
    }

    if let Some(resolution) = input.output.as_deref() {
        let kind = resolution.trim().to_ascii_lowercase();
        if !["interval", "zone", "mixture"].contains(&kind.as_str()) {
            return Err(Error::Invalid(format!(
                "snapshot output `{resolution}` is not an ALARA resolution \
                 (interval, zone, or mixture)"
            )));
        }
        push_block(&mut deck, "output", vec!["number_density".to_string()]);
        deck.outputs.push(OutputDef {
            resolution: kind,
            entries: vec!["number_density".to_string()],
            line: 0,
        });
    }

    deck.validate().map_err(map_alara)?;
    Ok(deck)
}

/// Build the template deck, derive the workflow (honoring per-zone flux
/// pins), validate it against the template, and emit one deck per step.
///
/// This is `deck_from_snapshot` → `from_deck` → per-zone flux pins →
/// `validate_against` → `emit_decks`, reusing the workflow machinery instead
/// of reimplementing it. Photon assembly stays out: it consumes post-ALARA
/// outputs via [`assemble`](crate::photon::assemble) and
/// [`from_photon_file`](crate::photon::from_photon_file).
pub fn snapshot_workflow(
    input: &SnapshotInput,
) -> Result<(R2sWorkflow, AlaraDeck, Vec<AlaraDeck>)> {
    let deck = deck_from_snapshot(input)?;
    let mut workflow = R2sWorkflow::from_deck(&deck)?;
    let pins: BTreeMap<&str, &str> = input
        .zones
        .iter()
        .filter_map(|zone| {
            zone.flux_name
                .as_deref()
                .map(|flux| (zone.zone.as_str(), flux))
        })
        .collect();
    if !pins.is_empty() {
        let defined: BTreeSet<&str> = deck.fluxes.iter().map(|flux| flux.name.as_str()).collect();
        for step in &mut workflow.steps {
            if let Some(pinned) = pins.get(step.zone.as_str()) {
                if !defined.contains(*pinned) {
                    return Err(Error::CrossRef(format!(
                        "snapshot zone `{}` pins undefined flux `{pinned}`",
                        step.zone
                    )));
                }
                step.flux = (*pinned).to_string();
            }
        }
    }
    workflow.validate_against(&deck)?;
    let decks = workflow.emit_decks(&deck)?;
    Ok((workflow, deck, decks))
}

/// Mixture name for a zone (`mix_<zone>`; the prefix keeps zones named
/// `void` distinct from the `void` mixture).
fn mixture_name(zone: &str) -> String {
    format!("mix_{zone}")
}

/// Record a raw block so the canonical writer replays the typed views in deck
/// order (synthesized blocks use line 0, like `emit_decks` does).
fn push_block(deck: &mut AlaraDeck, kind: &str, body: Vec<String>) {
    deck.blocks.push(RawBlock {
        kind: kind.to_string(),
        line: 0,
        body,
    });
}

/// Zone ids are opaque DECK tokens: non-empty with no whitespace.
fn check_zone_id(zone: &str) -> Result<()> {
    if zone.trim().is_empty() {
        return Err(Error::Invalid("snapshot zone id is empty".to_string()));
    }
    if zone.chars().any(char::is_whitespace) {
        return Err(Error::Invalid(format!(
            "snapshot zone `{zone}` contains whitespace (deck tokens must not)"
        )));
    }
    Ok(())
}

/// Volumes must be finite and positive (ALARA interval volumes in cm³).
fn check_volume(zone: &str, volume: f64) -> Result<()> {
    if !volume.is_finite() || volume <= 0.0 {
        return Err(Error::Invalid(format!(
            "snapshot zone `{zone}` has non-positive non-finite volume {volume} (cm³)"
        )));
    }
    Ok(())
}

/// Axial extents are informational but must be sane when both are dumped.
fn check_extent(zone: &SnapshotZone) -> Result<()> {
    for (label, value) in [
        ("zbottom_cm", zone.zbottom_cm),
        ("ztop_cm", zone.ztop_cm),
        ("temperature_C", zone.temperature_c),
    ] {
        if let Some(v) = value {
            if !v.is_finite() {
                return Err(Error::Invalid(format!(
                    "snapshot zone `{}` has non-finite {label} {v}",
                    zone.zone
                )));
            }
        }
    }
    if let (Some(bottom), Some(top)) = (zone.zbottom_cm, zone.ztop_cm) {
        if top < bottom {
            return Err(Error::Invalid(format!(
                "snapshot zone `{}` has ztop_cm {top} below zbottom_cm {bottom}",
                zone.zone
            )));
        }
    }
    Ok(())
}

/// Flux names are deck tokens referencing schedule items: non-empty, no
/// whitespace.
fn check_flux_name(name: &str) -> Result<()> {
    if name.trim().is_empty() {
        return Err(Error::Invalid("snapshot flux name is empty".to_string()));
    }
    if name.chars().any(char::is_whitespace) {
        return Err(Error::Invalid(format!(
            "snapshot flux `{name}` contains whitespace (deck tokens must not)"
        )));
    }
    Ok(())
}

/// Resolve one zone's composition to canonical ids, summing duplicate keys
/// (e.g. `U235` plus `nU235`) the way a volume homogenization would.
///
/// Key rules mirror the `nucleide-emit` ARMI-input rule: elemental keys fail
/// with "expand first", bare `AM242` fails (pass `AM242M`/`AM242G`
/// explicitly), and anything the nuclei ARMI bridge rejects fails naming the
/// key. Densities must be finite and non-negative; empty input means void.
fn resolve_composition(zone: &str, pairs: &[(String, f64)]) -> Result<BTreeMap<NuclideId, f64>> {
    let mut ids = BTreeMap::new();
    for (key, ndens) in pairs {
        let id = snapshot_key_to_nucid(key)?;
        if !ndens.is_finite() || *ndens < 0.0 {
            return Err(Error::Invalid(format!(
                "snapshot zone `{zone}` key `{key}` has number density {ndens} \
                 (expected finite atoms/barn-cm >= 0)"
            )));
        }
        *ids.entry(id).or_insert(0.0) += *ndens;
    }
    Ok(ids)
}

/// Resolve one ARMI-side key: v1 rejections first, then the nuclei bridge.
/// Mirrors the `nucleide-emit` ARMI-input rule; see the module docs.
fn snapshot_key_to_nucid(key: &str) -> Result<NuclideId> {
    let t = key.trim();
    if let Some(sym) = elemental_symbol(t).or_else(|| strip_db_prefix(t).and_then(elemental_symbol))
    {
        return Err(Error::Invalid(format!(
            "snapshot key `{t}` is elemental (`{sym}`): pass post-expansion nuclide \
             number densities; the adapter never reimplements \
             expandElementalMassFracsToNuclides — expand first"
        )));
    }
    match bare_core(t).as_deref() {
        Some("AM242") => {
            return Err(Error::Invalid(format!(
                "snapshot key `{t}` is ambiguous bare `AM242` (ARMI means the m-state): \
                 pass `AM242M` for Am-242m or `AM242G` for ground explicitly"
            )));
        }
        // Explicit ground alias (the bridge has no `G` suffix form).
        Some("AM242G") => {
            return NuclideId::new(95, 242, 0)
                .map_err(|e| Error::Invalid(format!("snapshot key `{t}`: {e}")));
        }
        _ => {}
    }
    armi_name_to_nucid(t).map_err(|e| Error::Invalid(format!("snapshot key `{t}`: {e}")))
}

/// Canonical element symbol when `t` is a bare elemental key (`ZR` → `Zr`),
/// else `None`. Only pure-letter input qualifies, so no nuclide key with
/// mass digits, tags, or separators can collide.
fn elemental_symbol(t: &str) -> Option<String> {
    if t.is_empty() || !t.bytes().all(|b| b.is_ascii_alphabetic()) {
        return None;
    }
    let mut chars = t.chars();
    let mut canon = String::with_capacity(t.len());
    if let Some(first) = chars.next() {
        canon.extend(first.to_uppercase());
    }
    canon.extend(chars.flat_map(|c| c.to_lowercase()));
    if nucleide_nuclei::element_z(&canon).is_some() {
        Some(canon)
    } else {
        None
    }
}

/// Strip one leading ARMI database `n`/`N` when followed by a letter
/// (`nZr` → `Zr`), mirroring the bridge so prefixed elementals are caught.
fn strip_db_prefix(t: &str) -> Option<&str> {
    let mut chars = t.chars();
    let first = chars.next()?;
    let second = chars.next()?;
    if (first == 'n' || first == 'N') && second.is_ascii_alphabetic() {
        Some(&t[1..])
    } else {
        None
    }
}

/// Uppercased key with the database prefix stripped (`nAm242` → `AM242`),
/// for the bare-`AM242` disambiguation check.
fn bare_core(t: &str) -> Option<String> {
    if t.is_empty() {
        return None;
    }
    let u = t.to_ascii_uppercase();
    let core = match u.strip_prefix('N') {
        Some(rest) if rest.starts_with(|c: char| c.is_ascii_alphabetic()) => rest,
        _ => u.as_str(),
    };
    Some(core.to_string())
}

/// Merge caller-supplied irradiation-history blocks into the deck.
///
/// Only `schedule` and `pulsehistory` blocks may cross this boundary; any
/// other block kind (geometry, fluxes, cooling, …) is a caller error naming
/// the kind, keeping the dict-in schema versionless.
fn merge_schedule_text(deck: &mut AlaraDeck, text: &str) -> Result<()> {
    let fragment =
        AlaraDeck::parse(&format!("geometry {SNAPSHOT_GEOMETRY}\n{text}")).map_err(map_alara)?;
    // Reject the wrapper plus anything outside the irradiation history.
    let mut bad: BTreeSet<&str> = BTreeSet::new();
    let mut geometries = 0;
    for block in &fragment.blocks {
        match block.kind.as_str() {
            "geometry" => geometries += 1,
            "schedule" | "pulsehistory" => {}
            other => {
                bad.insert(other);
            }
        }
    }
    if geometries > 1 {
        bad.insert("geometry");
    }
    if !bad.is_empty() {
        let kinds: Vec<&str> = bad.into_iter().collect();
        return Err(Error::Invalid(format!(
            "snapshot schedule_text must hold schedule/pulsehistory blocks only, \
             found {}",
            kinds.join(", ")
        )));
    }
    for block in fragment.blocks {
        if block.kind == "schedule" || block.kind == "pulsehistory" {
            deck.blocks.push(RawBlock { line: 0, ..block });
        }
    }
    deck.schedules = fragment.schedules;
    deck.pulse_histories = fragment.pulse_histories;
    Ok(())
}

/// Synthesize one 1-day-per-flux schedule plus a single-pulse history, so a
/// snapshot without `schedule_text` still derives a workflow and expands.
fn synthesize_schedule(deck: &mut AlaraDeck, input: &SnapshotInput) -> Result<()> {
    let items = input
        .flux_defs
        .iter()
        .map(|flux| ScheduleItemRef {
            tokens: vec![
                "1".to_string(),
                "d".to_string(),
                flux.name.clone(),
                SNAPSHOT_HISTORY.to_string(),
                "0".to_string(),
                "s".to_string(),
            ],
            line: 0,
        })
        .collect::<Vec<_>>();
    push_block(
        deck,
        "schedule",
        items.iter().map(|item| item.tokens.join(" ")).collect(),
    );
    deck.schedules.push(ScheduleDef {
        name: SNAPSHOT_SCHEDULE.to_string(),
        items,
        line: 0,
    });
    push_block(deck, "pulsehistory", vec!["1 0 s".to_string()]);
    deck.pulse_histories.push(PulseHistory {
        name: SNAPSHOT_HISTORY.to_string(),
        levels: vec![PulseLevel {
            pulses: 1,
            delay_s: 0.0,
        }],
        line: 0,
    });
    Ok(())
}

/// Map ALARA errors onto the R2S error type (mirrors `workflow::map_alara`).
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

    fn zone(zone: &str, volume: f64, composition: &[(&str, f64)]) -> SnapshotZone {
        SnapshotZone {
            zone: zone.to_string(),
            volume_cm3: volume,
            zbottom_cm: None,
            ztop_cm: None,
            material: None,
            xs_type: None,
            temperature_c: None,
            composition: composition
                .iter()
                .map(|(key, ndens)| ((*key).to_string(), *ndens))
                .collect(),
            flux_name: None,
        }
    }

    fn input() -> SnapshotInput {
        SnapshotInput {
            zones: vec![
                zone("B0120-005", 1200.0, &[("U235", 1.0e-3), ("nU238", 2.0e-2)]),
                zone("B0120-006", 800.0, &[("PU239", 5.0e-4)]),
            ],
            flux_defs: vec![SnapshotFluxDef {
                name: "snap_flux".to_string(),
                file: "data/flux".to_string(),
                scale: 1.0,
            }],
            cooling_s: vec![86_400.0],
            schedule_text: None,
            output: None,
        }
    }

    #[test]
    fn broadcast_builds_validated_deck_and_workflow() {
        let (workflow, deck, decks) = snapshot_workflow(&input()).unwrap();
        assert_eq!(workflow.steps.len(), 2);
        assert!(workflow.steps.iter().all(|step| step.flux == "snap_flux"));
        assert_eq!(workflow.cooling_s, [86_400.0]);
        assert_eq!(workflow.top_schedule, SNAPSHOT_SCHEDULE);
        workflow.validate_against(&deck).unwrap();
        assert_eq!(decks.len(), 2);
        assert_eq!(decks[0].solve_zones.as_ref().unwrap().zones, ["B0120-005"]);
        // Number densities round-trip through `mixture_ids` verbatim.
        let ids = deck.mixture_ids("mix_B0120-005").unwrap();
        assert_eq!(ids.len(), 2);
        assert!((ids[&NuclideId::from_name("U235").unwrap()] - 1.0e-3).abs() < 1e-15);
        // Canonical deck text re-parses and re-validates.
        let reparsed = AlaraDeck::parse(&deck.to_string()).unwrap();
        reparsed.validate().unwrap();
        assert_eq!(reparsed.to_string(), deck.to_string());
    }

    #[test]
    fn from_snapshot_delegates_to_from_deck() {
        let deck = deck_from_snapshot(&input()).unwrap();
        assert_eq!(
            R2sWorkflow::from_snapshot(&deck).unwrap(),
            R2sWorkflow::from_deck(&deck).unwrap()
        );
    }

    #[test]
    fn empty_composition_maps_to_void_and_skips() {
        let mut good = input();
        good.zones[1].composition.clear();
        let (workflow, deck, _) = snapshot_workflow(&good).unwrap();
        assert_eq!(workflow.steps.len(), 1);
        assert_eq!(workflow.steps[0].zone, "B0120-005");
        let loading = deck.mat_loading.as_ref().unwrap();
        assert_eq!(loading.entries[1].mixture, "void");
        assert!(deck.find_mixture("mix_B0120-006").is_none());
    }

    #[test]
    fn multiple_fluxes_match_by_zone_name() {
        let mut named = input();
        named.flux_defs = vec![
            SnapshotFluxDef {
                name: "B0120-005".to_string(),
                file: "data/f1".to_string(),
                scale: 1.0,
            },
            SnapshotFluxDef {
                name: "B0120-006".to_string(),
                file: "data/f2".to_string(),
                scale: 2.0,
            },
        ];
        let (workflow, _, _) = snapshot_workflow(&named).unwrap();
        assert_eq!(workflow.steps[0].flux, "B0120-005");
        assert_eq!(workflow.steps[1].flux, "B0120-006");
    }

    #[test]
    fn zone_flux_pin_overrides_broadcast() {
        let mut pinned = input();
        pinned.flux_defs.push(SnapshotFluxDef {
            name: "pin_flux".to_string(),
            file: "data/fpin".to_string(),
            scale: 1.0,
        });
        pinned.zones[0].flux_name = Some("pin_flux".to_string());
        // Two fluxes with no zone-name match: broadcast no longer applies, so
        // the pin is what makes zone B0120-005 resolvable — but B0120-006
        // still matches nothing.
        assert!(matches!(
            snapshot_workflow(&pinned),
            Err(Error::CrossRef(_))
        ));
        // Name every flux after its zone and the pin wins for its zone.
        pinned.flux_defs = vec![
            SnapshotFluxDef {
                name: "B0120-005".to_string(),
                file: "data/f1".to_string(),
                scale: 1.0,
            },
            SnapshotFluxDef {
                name: "B0120-006".to_string(),
                file: "data/f2".to_string(),
                scale: 1.0,
            },
            SnapshotFluxDef {
                name: "pin_flux".to_string(),
                file: "data/fpin".to_string(),
                scale: 1.0,
            },
        ];
        // Three fluxes, no broadcast: B0120-005 name-matches, then the pin
        // overrides it.
        let (workflow, _, _) = snapshot_workflow(&pinned).unwrap();
        assert_eq!(workflow.steps[0].flux, "pin_flux");
        assert_eq!(workflow.steps[1].flux, "B0120-006");
    }

    #[test]
    fn zone_flux_pin_to_unknown_flux_fails() {
        let mut bad = input();
        bad.zones[0].flux_name = Some("ghost".to_string());
        assert!(matches!(snapshot_workflow(&bad), Err(Error::CrossRef(_))));
    }

    #[test]
    fn duplicate_nuclide_keys_sum() {
        let mut dup = input();
        dup.zones[0].composition.push(("nU235".to_string(), 2.0e-3));
        let deck = deck_from_snapshot(&dup).unwrap();
        let ids = deck.mixture_ids("mix_B0120-005").unwrap();
        assert!((ids[&NuclideId::from_name("U235").unwrap()] - 3.0e-3).abs() < 1e-15);
    }

    #[test]
    fn empty_cooling_rejects_in_validate_against() {
        let mut cold = input();
        cold.cooling_s.clear();
        let deck = deck_from_snapshot(&cold).unwrap();
        assert!(deck.cooling.is_none());
        let workflow = R2sWorkflow::from_deck(&deck).unwrap();
        assert!(matches!(
            workflow.validate_against(&deck),
            Err(Error::Invalid(_))
        ));
        assert!(matches!(snapshot_workflow(&cold), Err(Error::Invalid(_))));
    }

    #[test]
    fn bad_inputs_fail() {
        let mut empty = input();
        empty.zones.clear();
        assert!(matches!(deck_from_snapshot(&empty), Err(Error::Invalid(_))));

        let mut noflux = input();
        noflux.flux_defs.clear();
        assert!(matches!(
            deck_from_snapshot(&noflux),
            Err(Error::Invalid(_))
        ));

        let mut dup = input();
        dup.zones.push(dup.zones[0].clone());
        assert!(matches!(deck_from_snapshot(&dup), Err(Error::Invalid(_))));

        let mut badvol = input();
        badvol.zones[0].volume_cm3 = 0.0;
        assert!(matches!(
            deck_from_snapshot(&badvol),
            Err(Error::Invalid(_))
        ));

        let mut badzone = input();
        badzone.zones[0].zone = "has space".to_string();
        assert!(matches!(
            deck_from_snapshot(&badzone),
            Err(Error::Invalid(_))
        ));

        let mut badextent = input();
        badextent.zones[0].zbottom_cm = Some(10.0);
        badextent.zones[0].ztop_cm = Some(5.0);
        assert!(matches!(
            deck_from_snapshot(&badextent),
            Err(Error::Invalid(_))
        ));
    }

    #[test]
    fn composition_key_rules_mirror_emit() {
        for bad in ["DUMP1", "nXx999", "NOTANUCLIDE", ""] {
            let mut keyed = input();
            keyed.zones[0].composition = vec![(bad.to_string(), 1.0e-3)];
            let err = deck_from_snapshot(&keyed).unwrap_err();
            assert!(matches!(err, Error::Invalid(_)), "{bad}: {err}");
            assert!(err.to_string().contains(bad), "{bad}: {err}");
        }
        for elemental in ["ZR", "zr", "nZr", "O"] {
            let mut keyed = input();
            keyed.zones[0].composition = vec![(elemental.to_string(), 1.0e-3)];
            let err = deck_from_snapshot(&keyed).unwrap_err();
            assert!(
                err.to_string().contains("expand first"),
                "{elemental}: {err}"
            );
        }
        for bare in ["AM242", "nAm242"] {
            let mut keyed = input();
            keyed.zones[0].composition = vec![(bare.to_string(), 1.0e-3)];
            let err = deck_from_snapshot(&keyed).unwrap_err();
            assert!(err.to_string().contains("AM242M"), "{bare}: {err}");
        }
        // Explicit isomer states resolve; negative densities fail.
        let mut states = input();
        states.zones[0].composition = vec![
            ("AM242M".to_string(), 1.0e-4),
            ("AM242G".to_string(), 2.0e-4),
        ];
        let deck = deck_from_snapshot(&states).unwrap();
        let ids = deck.mixture_ids("mix_B0120-005").unwrap();
        assert!(ids.contains_key(&NuclideId::new(95, 242, 1).unwrap()));
        assert!(ids.contains_key(&NuclideId::new(95, 242, 0).unwrap()));

        let mut neg = input();
        neg.zones[0].composition = vec![("U235".to_string(), -1.0)];
        assert!(matches!(deck_from_snapshot(&neg), Err(Error::Invalid(_))));
    }

    #[test]
    fn schedule_text_merges_and_rejects_non_history_blocks() {
        let mut custom = input();
        custom.schedule_text = Some(
            "schedule op\n1 d snap_flux once 0 s\nend\npulsehistory once\n1 0 s\nend\n".to_string(),
        );
        let (workflow, deck, _) = snapshot_workflow(&custom).unwrap();
        assert_eq!(workflow.top_schedule, "op");
        assert_eq!(workflow.expand(&deck, &[]).unwrap().len(), 1);

        let mut bad = input();
        bad.schedule_text = Some("cooling\n1 d\nend\n".to_string());
        assert!(matches!(deck_from_snapshot(&bad), Err(Error::Invalid(_))));

        let mut broken = input();
        broken.schedule_text = Some("schedule op\nbogus\nend\n".to_string());
        assert!(matches!(
            deck_from_snapshot(&broken),
            Err(Error::Invalid(_))
        ));

        // Dangling flux references in custom schedules surface as CrossRef.
        let mut dangling = input();
        dangling.schedule_text = Some(
            "schedule op\n1 d ghost_flux once 0 s\nend\npulsehistory once\n1 0 s\nend\n"
                .to_string(),
        );
        assert!(matches!(
            deck_from_snapshot(&dangling),
            Err(Error::CrossRef(_))
        ));
    }

    #[test]
    fn output_resolution_validated() {
        let mut out = input();
        out.output = Some("zone".to_string());
        let deck = deck_from_snapshot(&out).unwrap();
        assert_eq!(deck.outputs[0].resolution, "zone");
        let reparsed = AlaraDeck::parse(&deck.to_string()).unwrap();
        reparsed.validate().unwrap();

        let mut bad = input();
        bad.output = Some("cell".to_string());
        assert!(matches!(deck_from_snapshot(&bad), Err(Error::Invalid(_))));
    }
}
