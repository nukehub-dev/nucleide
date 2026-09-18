//! Error type for the `equilib-io` crate.

use thiserror::Error;

/// Result alias for the `equilib-io` crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Wrong on-disk netCDF variant: the classic reader only accepts CDF-1 and
/// CDF-2. Anything else must be converted facade-side (outside Rust —
/// no HDF5 in Rust, by rule) before re-reading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WrongVariant {
    /// HDF5-backed netCDF-4 (8-byte HDF5 signature `89 48 44 46 0D 0A 1A 0A`).
    Hdf5NetCdf4,
    /// 64-bit-data CDF-5 (`43 44 46 05`).
    Cdf5,
}

impl std::fmt::Display for WrongVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WrongVariant::Hdf5NetCdf4 => write!(f, "HDF5-backed netCDF-4"),
            WrongVariant::Cdf5 => write!(f, "CDF-5 64-bit-data"),
        }
    }
}

/// Errors raised while probing a file variant, parsing a classic-netCDF
/// header, extracting `wout` fields, evaluating flux-surface Jacobians, or
/// parsing `&INDATA` namelist text.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum Error {
    /// A file could not be read: `{0}`.
    #[error("equilib-io: file read failed: {0}")]
    Io(String),
    /// Unexpected end of input while reading `{0}`.
    #[error("equilib-io: truncated file while reading {0}")]
    Truncated(&'static str),
    /// The first four bytes are not a classic-netCDF magic
    /// (`43 44 46 01` / `43 44 46 02`); got bytes `{magic:02X?}`.
    #[error("equilib-io: not a classic-netCDF file (magic {magic:02X?})")]
    NotClassic {
        /// The offending first four bytes.
        magic: [u8; 4],
    },
    /// The file is a real netCDF file of the wrong variant (`{0}`).
    /// Convert facade-side (e.g. `ncdump`/`ncap2` to classic or 64-bit
    /// offset) and re-read; the Rust reader never takes HDF5.
    #[error("equilib-io: wrong netCDF variant ({0}); convert to classic/64-bit offset facade-side and re-read")]
    WrongVariant(WrongVariant),
    /// A classic header tag, count, or layout is malformed: `{0}`.
    #[error("equilib-io: malformed classic header: {0}")]
    BadHeader(&'static str),
    /// An unknown netCDF element type code: `{0}` (classic types are 1–6).
    #[error("equilib-io: unknown netCDF type code {0}")]
    UnknownType(u32),
    /// Variable `{0}` sits on the unlimited dimension (a record variable).
    /// Classic `wout` files carry none; convert facade-side.
    #[error("equilib-io: record variable `{0}` unsupported (no unlimited dimension in classic wout files)")]
    RecordVariable(String),
    /// Required dimension `{0}` is absent.
    #[error("equilib-io: missing dimension `{0}`")]
    MissingDimension(String),
    /// Required variable `{0}` is absent.
    #[error("equilib-io: missing variable `{0}`")]
    MissingVariable(String),
    /// An input length is inconsistent: `{what}` holds `{got}` entries,
    /// expected `{expected}`.
    #[error("equilib-io: shape mismatch: {what} holds {got} entries, expected {expected}")]
    BadShape {
        /// Name of the offending input.
        what: &'static str,
        /// Length the input must carry.
        expected: usize,
        /// Length actually supplied.
        got: usize,
    },
    /// Variable `{0}` has an unexpected dimension list (the `wout`
    /// Fourier tables must sit on `(radius, mode)` in that order).
    #[error("equilib-io: variable `{0}` has unexpected dimensions")]
    UnexpectedDims(String),
    /// An index is out of range: `{what}` index `{index}` with length `{len}`.
    #[error("equilib-io: {what} index {index} out of range (length {len})")]
    OutOfRange {
        /// Name of the offending input.
        what: &'static str,
        /// Index actually supplied.
        index: usize,
        /// Length of the input.
        len: usize,
    },
    /// A value is invalid: `{0}`.
    #[error("equilib-io: invalid value: {0}")]
    BadValue(String),
    /// An `&INDATA` value has the wrong kind: `{name}` is not {expected}.
    #[error("equilib-io: `{name}` is not {expected}")]
    MismatchedType {
        /// Upper-cased variable name.
        name: String,
        /// Expected kind (`an integer`, `a float`, ...).
        expected: &'static str,
    },
    /// No `&INDATA` namelist block found in the text.
    #[error("equilib-io: no &INDATA namelist block found")]
    MissingIndata,
    /// The `&INDATA` block never reaches its `/` terminator.
    #[error("equilib-io: unterminated &INDATA block (missing `/`)")]
    UnterminatedIndata,
    /// Sliced index assignment to `{name}` (e.g. `RBC(0:4,2)=...`).
    /// The grammar takes scalar indices only, per the fixed-boundary
    /// conversion stance.
    #[error("equilib-io: sliced index on `{name}` unsupported (scalar indices only)")]
    SlicedIndex {
        /// Upper-cased variable name.
        name: String,
    },
    /// A namelist token, number, or assignment is malformed: `{0}`.
    #[error("equilib-io: malformed namelist: {0}")]
    BadNamelist(String),
}
