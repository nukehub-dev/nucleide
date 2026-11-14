#![warn(missing_docs)]
//! Rigorous two-step (R2S) shutdown-dose-rate orchestration.
//!
//! Combines neutron flux meshes (`mcnp-io`), ALARA activation decks
//! (`alara-io`), decay photon sources, and mesh source sampling (`vr-tools`)
//! into one reproducible workflow description. Transport and activation
//! solving stay inside their respective codes; this crate only builds decks,
//! links zones to spectra, and assembles photon sources.

pub mod error;
pub mod photon;
pub mod snapshot;
pub mod sweep;
pub mod tags;
pub mod workflow;

pub use error::{Error, Result};
pub use photon::ZonePhotonSource;
pub use snapshot::{
    deck_from_snapshot, snapshot_inventory, snapshot_workflow, SnapshotFluxDef, SnapshotInput,
    SnapshotZone,
};
pub use sweep::{assemble_results, expand_sweep, SweepAxis, SweepCase};
pub use tags::{
    photon_groups_at, split_zone_totals, split_zone_totals_by_volume, sum_group_strengths,
    tag_zone_totals, VoxelTags,
};
pub use workflow::{R2sStep, R2sWorkflow};
