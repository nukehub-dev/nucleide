//! MCNP-family file I/O: xsdir, meshtal, SSW (surfsrc), PTRAC, WWINP,
//! MCTAL, input-deck materials/cells/surfaces/problems, and mesh-to-geometry
//! deck generation.
//!
//! Parsers are validated against the golden-byte fixture set under
//! `fixtures/mcnp/`.

pub mod cell;
pub mod deck;
pub mod endl;
pub mod fortran;
pub mod inp;
pub mod mctal;
pub mod meshtal;
pub mod problem;
pub mod ptrac;
pub mod semantic;
pub mod surf;
pub mod surfsrc;
pub mod wwinp;
pub mod xsdir;
