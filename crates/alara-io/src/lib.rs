//! ALARA activation-code interop: input-deck, flux, schedule, and output glue.
//!
//! This crate is **glue, not a solver**. It reads the text files ALARA
//! consumes and produces (deck blocks, group fluxes, irradiation schedules,
//! library references, activation and photon-source summaries) and repackages
//! them for downstream workflows. Transmutation solving itself stays inside
//! ALARA; time stepping here only carries durations in plain `f64` seconds,
//! compatible with the `depletion` crate's step sizes.
//!
//! Nuclide identity reuses [`nuclei::NuclideId`]; every fallible entry point
//! returns [`Error`] (never panics) with file/line context where known.
//!
//! ## Out of scope (explicitly)
//!
//! - ALARA solver core (activation/transmutation mathematics)
//! - `.lib` binary cross-section libraries
//! - jALARA (Java GUI / workflow tooling)
//! - `dant2alara` and other Fortran helpers
//! - ALARA Perl post-processing scripts
//!
//! ```
//! use alara_io::parse_time_to_seconds;
//!
//! assert_eq!(parse_time_to_seconds(2.0, "h").unwrap(), 7200.0);
//! ```

pub mod deck;
pub mod error;
pub mod flux;
pub mod libs;
pub mod output;
pub mod photon;
pub mod schedule;

pub use deck::{AlaraDeck, KNOWN_BLOCKS};
pub use error::{Error, Result};
pub use flux::{FluxSpec, FluxSpectra};
pub use libs::{
    canonical_element_z, split_enriched, EleLib, ElementEntry, IsotopeAbund, LibEntry, LibSpec,
    MatEntry, MatLib, MaterialEntry, WdrLib, WdrLimit,
};
pub use output::{
    ActivationOutput, ActivationRecord, BlockKind, ResponseFrame, ResponseRow, ResponseVar,
};
pub use photon::{PhotonGroup, PhotonSource};
pub use schedule::{
    expand, expand_from, parse_time_to_seconds, total_time, FlatStep, PulseHistory, PulseLevel,
    SchedItem, Schedule, ScheduleDef, ScheduleItem,
};
