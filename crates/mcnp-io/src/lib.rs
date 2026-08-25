//! MCNP-family file I/O: xsdir, meshtal, SSW (surfsrc), PTRAC, WWINP,
//! MCTAL, input-deck materials, and mesh-to-geometry deck generation.
//!
//! Parsers are validated against the golden-byte fixture set under
//! `fixtures/mcnp/`.

pub mod deck;
pub mod inp;
pub mod mctal;
pub mod meshtal;
pub mod ptrac;
pub mod surfsrc;
pub mod wwinp;
pub mod xsdir;
