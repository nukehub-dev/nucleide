//! Scoped ORIGEN 2.2 TAPE readers.
//!
//! Covers the decay-data path (`TAPE9`-style decay constants), `TAPE5` input
//! echo, and `TAPE6` output inventories. Full ORIGEN burnup driving stays
//! inside ORIGEN; this crate only reads its text interfaces.

pub mod error;
pub mod tape5;
pub mod tape6;
pub mod tape9;

pub use error::{Error, Result};
pub use tape5::Tape5;
pub use tape6::Tape6;
pub use tape9::Tape9Entry;
