#![warn(missing_docs)]
//! ALARA activation-code interop: input-deck, flux, schedule, and output glue.
//!
//! This crate is **glue, not a solver**. It reads the text files ALARA
//! consumes and produces (deck blocks, group fluxes, irradiation schedules,
//! library references, activation and photon-source summaries) and repackages
//! them for downstream workflows. Transmutation solving itself stays inside
//! ALARA; time stepping here only carries durations in plain `f64` seconds,
//! compatible with the `depletion` crate's step sizes.
//!
//! Nuclide identity reuses [`nucleide_nuclei::NuclideId`]; every fallible entry point
//! returns [`Error`] (never panics) with file/line context where known.
//!
//! Clearance / waste-classification analytics (clearance index and the
//! sum-of-fractions rule over parsed inventories) live in [`clearance`].
//!
//! Sublet S1+S2+S4+S5 radiological totals (total activity with the IRT
//! α/β/γ split, decay heat per radiation class, and committed
//! ingestion/inhalation hazards over caller inventories) live in
//! [`sublet`]; the S3 gamma dose-rate kernel (slab and point dose over
//! caller gamma groups) lives in [`dose`].
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
//! use nucleide_alara_io::parse_time_to_seconds;
//!
//! assert_eq!(parse_time_to_seconds(2.0, "h").unwrap(), 7200.0);
//! ```

pub mod clearance;
pub mod deck;
pub mod dose;
pub mod error;
pub mod flux;
pub mod libs;
pub mod output;
pub mod photon;
pub mod schedule;
pub mod sublet;

pub use clearance::{
    clearance_index, inventory_from_frame, sum_of_fractions, ClearanceClass, ClearanceTable,
    EsNormMaterial, SumOfFractions,
};
pub use deck::{AlaraDeck, KNOWN_BLOCKS};
pub use dose::{
    dose_point, dose_slab, mixture_mu, DoseGroup, PointDose, SlabDose, DOSE_CONVERSION_C,
    ELEMENTARY_CHARGE, MIN_DISTANCE_M, SLAB_BUILDUP_B,
};
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
pub use sublet::{
    decay_heat, iaea_clearance_index, ingestion_hazard, inhalation_hazard, total_activity,
    transport_ratio, ActivityEntry, ActivityTotal, DecayHeat, DecayHeatEntry, HazardDose,
    HazardEntry, IaeaClearance, IaeaEntry, TransportEntry, TransportRatio, TBQ_TO_BQ,
};
