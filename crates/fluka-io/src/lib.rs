//! FLUKA Monte Carlo code interface.
//!
//! - [`usrbin`]: reads USRBIN tally blocks from `.lis` output files
//!   (detector metadata, mesh bounds, track-length data and percentage
//!   errors) into plain vectors — no meshing layer.
//! - [`material`]: the built-in FLUKA material table and generators for
//!   `MATERIAL` / `COMPOUND` input-deck cards with 10-column
//!   Fortran-style fields.

pub mod material;
pub mod usrbin;
