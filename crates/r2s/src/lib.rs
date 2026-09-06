//! Rigorous two-step (R2S) shutdown-dose-rate orchestration.
//!
//! Combines neutron flux meshes (`mcnp-io`), ALARA activation decks
//! (`alara-io`), decay photon sources, and mesh source sampling (`vr-tools`)
//! into one reproducible workflow description. Transport and activation
//! solving stay inside their respective codes; this crate only builds decks,
//! links zones to spectra, and assembles photon sources.

pub mod error;
pub mod photon;
pub mod workflow;

pub use error::{Error, Result};
pub use photon::ZonePhotonSource;
pub use workflow::{R2sStep, R2sWorkflow};
