//! Equilibrium data readers: classic-netCDF VMEC `wout` files, the
//! VMEC `&INDATA` input text grammar, and flux-surface Jacobian helpers.
//!
//! Reads equilibrium data, never solves it. The `wout` side is a pure
//! classic-format parser ([`classic`]): CDF-1 (`43 44 46 01`) and CDF-2
//! (`43 44 46 02`) files are accepted — `ncdump -k` must report `classic`
//! or `64-bit offset` — while HDF5-backed netCDF-4 files and CDF-5 files
//! are rejected loudly ([`Error::WrongVariant`]) and stay facade-side
//! conversion only. The text side ([`indata`]) parses the documented
//! `&INDATA` namelist (power-series profiles, `NCURR = 1` as the `I'(s)`
//! profile, fixed-boundary stance, scalar indices only). [`wout`] turns
//! the `gmnc` Jacobian spectrum into point and grid evaluations (J1–J3)
//! for volume weighting; [`wall`] maps a caller birth-rate field onto a
//! caller wall surface in the same flux coordinates (W1–W3) with those
//! Jacobians — closed-form per-cell accumulation, no transport.
//!
//! Synthetic fixtures for the gates below are built from the published
//! variable lists only — never solver output.

#![warn(missing_docs)]

pub mod classic;
pub mod error;
pub mod indata;
pub mod wall;
pub mod wout;

pub use classic::{probe_file as probe_classic_file, ClassicFile, ClassicVariant};
pub use error::{Error, Result, WrongVariant};
pub use indata::{parse_indata, Indata, Value as IndataValue};
pub use wall::{wall_load, WallLoad};
pub use wout::{fourier_jacobian, Matrix, Wout};
