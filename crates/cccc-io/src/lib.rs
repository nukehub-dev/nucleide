//! CCCC binary-standard readers and PARTISN deck writer.
//!
//! Glue around legacy deterministic-transport files: ISOTXS / RTFLUX / ATFLUX /
//! RZFLUX readers plus a 1D/2D/3D PARTISN input writer with ISOTXS nuclide
//! mapping. No transport solving is performed here.

pub mod error;
pub mod isotxs;
pub mod partisn;
pub mod rtflux;

pub use error::{Error, Result};
pub use isotxs::IsotxsLib;
pub use partisn::PartisnDeck;
pub use rtflux::FluxFile;
