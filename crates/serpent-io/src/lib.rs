//! Parsers for [Serpent](https://mon-jeu.vtt.fi/) Monte Carlo MATLAB-style
//! output files (`*_res.m`, `*_dep.m`, `*_det.m`).
//!
//! Serpent writes plain-text "MATLAB" files where every quantity is an
//! assignment such as `NAME = [values];` or, for multi-step result files,
//! repeated `NAME(idx, [1: n]) = [...];` blocks. This crate turns such files
//! into a [`Table`] — a sorted map from variable name to [`Entry`].
//!
//! # Entry points
//!
//! | File kind | Parser | Typical contents |
//! |-----------|--------|------------------|
//! | `*_res.m` | [`parse_res`] | run metadata, k-eigenvalues, cycle statistics |
//! | `*_dep.m` | [`parse_dep`] | depletion inventories (`ZAI`, `NAMES`, densities) |
//! | `*_det.m` | [`parse_det`] | detector bins, values and relative errors |
//!
//! Each parser accepts the file contents as `&str`; use [`from_file`] with a
//! [`Kind`] to read from disk.
//!
//! # Example
//!
//! ```
//! let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/serpent/sample_res.m")).unwrap();
//! let res = serpent_io::parse_res(&text).unwrap();
//! assert_eq!(res.get_f64("IDX").unwrap(), 3.0);
//! let keff = res.get_matrix("IMP_KEFF").unwrap();
//! assert_eq!((keff.row_f64(0).unwrap()[0] * 1e5).round(), 124_207.0);
//! ```

use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::path::Path;

mod dep;
mod det;
pub(crate) mod matlab;
mod res;

pub use dep::parse_dep;
pub use det::parse_det;
pub use res::parse_res;

/// A single MATLAB value: a double or a character string.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Num(f64),
    Str(String),
}

impl Value {
    /// Numeric value, or an error if this is a string.
    pub fn as_f64(&self) -> Result<f64> {
        match self {
            Value::Num(v) => Ok(*v),
            Value::Str(s) => Err(Error::Type {
                context: format!("string value `{s}`"),
                expected: "number",
                found: "string",
            }),
        }
    }

    /// String value, or an error if this is numeric.
    pub fn as_str(&self) -> Result<&str> {
        match self {
            Value::Str(s) => Ok(s),
            Value::Num(v) => Err(Error::Type {
                context: format!("numeric value `{v}`"),
                expected: "string",
                found: "number",
            }),
        }
    }
}

/// Row-major dense matrix of [`Value`]s with recorded dimensions.
///
/// Serpent writes depletion inventories as large rectangular blocks; these
/// are stored as one contiguous vector plus `(rows, cols)` rather than as
/// nested vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<Value>,
}

impl Matrix {
    pub(crate) fn from_parts(rows: usize, cols: usize, data: Vec<Value>) -> Matrix {
        Matrix { rows, cols, data }
    }

    /// Build a matrix from complete rows; errors on ragged or empty input.
    pub fn from_rows(rows: Vec<Vec<Value>>) -> Result<Matrix> {
        let cols = rows.first().map(|r| r.len()).ok_or(Error::Type {
            context: "matrix".into(),
            expected: "at least one row",
            found: "no rows",
        })?;
        if cols == 0 || rows.iter().any(|r| r.len() != cols) {
            return Err(Error::Type {
                context: "matrix".into(),
                expected: "rectangular rows",
                found: "ragged rows",
            });
        }
        let n = rows.len();
        Ok(Matrix {
            rows: n,
            cols,
            data: rows.into_iter().flatten().collect(),
        })
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Flattened row-major backing storage.
    pub fn data(&self) -> &[Value] {
        &self.data
    }

    /// Element access without bounds panics.
    pub fn get(&self, row: usize, col: usize) -> Option<&Value> {
        if row >= self.rows || col >= self.cols {
            return None;
        }
        self.data.get(row * self.cols + col)
    }

    /// Borrowed row slice, bounds-checked.
    pub fn row(&self, row: usize) -> Result<&[Value]> {
        if row >= self.rows {
            return Err(Error::Index {
                index: row,
                length: self.rows,
            });
        }
        Ok(&self.data[row * self.cols..(row + 1) * self.cols])
    }

    /// Numeric copy of one row, bounds-checked.
    pub fn row_f64(&self, row: usize) -> Result<Vec<f64>> {
        self.row(row)?
            .iter()
            .map(Value::as_f64)
            .collect::<Result<Vec<f64>>>()
            .map_err(|e| Error::Type {
                context: "matrix row".into(),
                expected: "numbers",
                found: e.found(),
            })
    }

    /// Numeric value at `(row, col)`, bounds-checked.
    pub fn get_f64(&self, row: usize, col: usize) -> Result<f64> {
        self.row(row)?[col].as_f64()
    }

    /// Full contents as nested vectors of doubles.
    pub fn to_rows_f64(&self) -> Result<Vec<Vec<f64>>> {
        (0..self.rows).map(|r| self.row_f64(r)).collect()
    }

    pub(crate) fn push_row(&mut self, values: Vec<Value>) {
        debug_assert_eq!(values.len(), self.cols);
        self.rows += 1;
        self.data.extend(values);
    }
}

/// One variable's worth of parsed Serpent output.
#[derive(Debug, Clone, PartialEq)]
pub enum Entry {
    Scalar(Value),
    Vector(Vec<Value>),
    Matrix(Matrix),
}

impl Entry {
    /// Kind label used in error messages.
    pub fn type_name(&self) -> &'static str {
        match self {
            Entry::Scalar(_) => "scalar",
            Entry::Vector(_) => "vector",
            Entry::Matrix(_) => "matrix",
        }
    }

    /// Scalar value, or a shape error otherwise.
    pub fn as_f64(&self) -> Result<f64> {
        match self {
            Entry::Scalar(v) => v.as_f64(),
            other => Err(Error::Type {
                context: "entry".into(),
                expected: "scalar",
                found: other.type_name(),
            }),
        }
    }

    /// Scalar string value, or a shape error otherwise.
    pub fn as_str(&self) -> Result<&str> {
        match self {
            Entry::Scalar(Value::Str(s)) => Ok(s),
            other => Err(Error::Type {
                context: "entry".into(),
                expected: "scalar string",
                found: other.type_name(),
            }),
        }
    }

    /// Vector elements, or a shape error otherwise.
    pub fn as_slice(&self) -> Result<&[Value]> {
        match self {
            Entry::Vector(vals) => Ok(vals),
            other => Err(Error::Type {
                context: "entry".into(),
                expected: "vector",
                found: other.type_name(),
            }),
        }
    }

    /// Numeric copies of a vector's elements.
    pub fn as_vec_f64(&self) -> Result<Vec<f64>> {
        self.as_slice()?.iter().map(Value::as_f64).collect()
    }

    /// Matrix view, or a shape error otherwise.
    pub fn as_matrix(&self) -> Result<&Matrix> {
        match self {
            Entry::Matrix(m) => Ok(m),
            other => Err(Error::Type {
                context: "entry".into(),
                expected: "matrix",
                found: other.type_name(),
            }),
        }
    }
}

/// Sorted map from MATLAB variable name to parsed [`Entry`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table(BTreeMap<String, Entry>);

impl Table {
    /// Empty table.
    pub fn new() -> Table {
        Table(BTreeMap::new())
    }

    /// Look up a variable by name.
    pub fn get(&self, key: impl AsRef<str>) -> Option<&Entry> {
        self.0.get(key.as_ref())
    }

    /// Mutable look up by name.
    pub fn get_mut(&mut self, key: impl AsRef<str>) -> Option<&mut Entry> {
        self.0.get_mut(key.as_ref())
    }

    /// Insert a variable, returning the previous one if any.
    pub fn insert(&mut self, key: impl Into<String>, entry: Entry) -> Option<Entry> {
        self.0.insert(key.into(), entry)
    }

    /// True when `key` is present.
    pub fn contains_key(&self, key: impl AsRef<str>) -> bool {
        self.0.contains_key(key.as_ref())
    }

    /// Iterate over `(name, entry)` pairs in sorted order.
    pub fn iter(&self) -> std::collections::btree_map::Iter<'_, String, Entry> {
        self.0.iter()
    }

    /// Iterate over variable names in sorted order.
    pub fn keys(&self) -> std::collections::btree_map::Keys<'_, String, Entry> {
        self.0.keys()
    }

    /// Number of variables.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True when no variables are present.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn require(&self, key: &str) -> Result<&Entry> {
        self.get(key).ok_or_else(|| Error::Missing(key.to_string()))
    }

    /// Scalar numeric value of `key`.
    pub fn get_f64(&self, key: impl AsRef<str>) -> Result<f64> {
        self.require(key.as_ref())?
            .as_f64()
            .map_err(|e| e.with_key(key.as_ref()))
    }

    /// Scalar string value of `key`.
    pub fn get_str(&self, key: impl AsRef<str>) -> Result<&str> {
        self.require(key.as_ref())?
            .as_str()
            .map_err(|e| e.with_key(key.as_ref()))
    }

    /// Numeric copy of the vector stored under `key`.
    pub fn get_vec_f64(&self, key: impl AsRef<str>) -> Result<Vec<f64>> {
        self.require(key.as_ref())?
            .as_vec_f64()
            .map_err(|e| e.with_key(key.as_ref()))
    }

    /// String copies of the vector stored under `key`.
    pub fn get_vec_str(&self, key: impl AsRef<str>) -> Result<Vec<String>> {
        let vals = self
            .require(key.as_ref())?
            .as_slice()
            .map_err(|e| e.with_key(key.as_ref()))?;
        vals.iter()
            .map(|v| v.as_str().map(str::to_string))
            .collect::<Result<Vec<_>>>()
            .map_err(|e| e.with_key(key.as_ref()))
    }

    /// Matrix stored under `key`.
    pub fn get_matrix(&self, key: impl AsRef<str>) -> Result<&Matrix> {
        self.require(key.as_ref())?
            .as_matrix()
            .map_err(|e| e.with_key(key.as_ref()))
    }
}

impl<'a> IntoIterator for &'a Table {
    type Item = (&'a String, &'a Entry);
    type IntoIter = std::collections::btree_map::Iter<'a, String, Entry>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl IntoIterator for Table {
    type Item = (String, Entry);
    type IntoIter = std::collections::btree_map::IntoIter<String, Entry>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl FromIterator<(String, Entry)> for Table {
    fn from_iter<I: IntoIterator<Item = (String, Entry)>>(iter: I) -> Table {
        Table(BTreeMap::from_iter(iter))
    }
}

/// Which kind of Serpent output file to parse in [`from_file`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// `*_res.m` results file.
    Res,
    /// `*_dep.m` depletion file.
    Dep,
    /// `*_det.m` detector file.
    Det,
}

/// Read and parse a Serpent output file of the given kind.
pub fn from_file<P: AsRef<Path>>(path: P, kind: Kind) -> Result<Table> {
    let text = fs::read_to_string(path.as_ref()).map_err(Error::Io)?;
    match kind {
        Kind::Res => parse_res(&text),
        Kind::Dep => parse_dep(&text),
        Kind::Det => parse_det(&text),
    }
}

/// Local error type for parsing and typed lookups.
#[derive(Debug)]
pub enum Error {
    /// Underlying file I/O failure.
    Io(std::io::Error),
    /// Malformed input syntax.
    Syntax {
        /// One-based line number where the problem was noticed.
        line: usize,
        /// Human-readable description.
        message: String,
    },
    /// Requested variable is absent.
    Missing(String),
    /// Variable exists but has the wrong kind/shape.
    Type {
        /// Variable name or value description.
        context: String,
        /// What was requested.
        expected: &'static str,
        /// What was actually present.
        found: &'static str,
    },
    /// Index out of range.
    Index {
        /// Offending index.
        index: usize,
        /// Valid length.
        length: usize,
    },
    /// Right-hand-side expression outside the supported subset.
    UnsupportedExpr(String),
}

impl Error {
    fn with_key(mut self, key: &str) -> Error {
        if let Error::Type {
            ref mut context, ..
        } = self
        {
            if context == "entry" {
                *context = key.to_string();
            }
        }
        self
    }

    fn found(&self) -> &'static str {
        match self {
            Error::Type { found, .. } => found,
            _ => "non-numeric value",
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Syntax { line, message } => write!(f, "syntax error at line {line}: {message}"),
            Error::Missing(key) => write!(f, "variable `{key}` not found"),
            Error::Type {
                context,
                expected,
                found,
            } => {
                write!(f, "`{context}` has {found}, expected {expected}")
            }
            Error::Index { index, length } => {
                write!(f, "index {index} out of range (length {length})")
            }
            Error::UnsupportedExpr(expr) => {
                write!(f, "unsupported expression `{expr}`")
            }
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// Crate-local result alias.
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
pub(crate) mod testing {
    use super::{from_file, Kind, Table};
    use std::path::PathBuf;

    const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/serpent");

    pub(crate) fn load_res(file: &str) -> Table {
        load(Kind::Res, file)
    }

    pub(crate) fn load_dep(file: &str) -> Table {
        load(Kind::Dep, file)
    }

    pub(crate) fn load_det(file: &str) -> Table {
        load(Kind::Det, file)
    }

    fn load(kind: Kind, file: &str) -> Table {
        let path = PathBuf::from(DIR).join(file);
        from_file(path, kind).unwrap_or_else(|e| panic!("failed to parse {file}: {e}"))
    }

    /// Element-wise comparison with a relative tolerance of 1e-9.
    pub(crate) fn assert_close(got: &[f64], expected: &[f64]) {
        assert_eq!(
            got.len(),
            expected.len(),
            "length mismatch: {got:?} vs {expected:?}"
        );
        for (x, y) in got.iter().zip(expected) {
            let tol = 1e-9 * x.abs().max(y.abs()).max(1.0);
            assert!((x - y).abs() <= tol, "value mismatch: {x} vs {y}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Error, Kind};
    use std::path::PathBuf;

    #[test]
    fn accessor_errors_are_typed() {
        let res = crate::testing::load_res("sample_res.m");
        let missing = res.get_f64("no_such_key").unwrap_err();
        assert_eq!(missing.to_string(), "variable `no_such_key` not found");
        for (key, expected) in [
            ("SIX_FF_ETA", "matrix"),
            ("POP", "vector"),
            ("VERSION", "vector"),
        ] {
            match res.get_f64(key) {
                Err(Error::Type { found, .. }) => assert_eq!(found, expected),
                other => panic!("expected Type error for {key}, got {other:?}"),
            }
        }
    }

    #[test]
    fn missing_file_is_an_io_error() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/serpent/definitely_absent.m");
        match crate::from_file(path, Kind::Res) {
            Err(Error::Io(_)) => {}
            other => panic!("expected Io error, got {other:?}"),
        }
    }
}
