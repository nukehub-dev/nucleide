//! Minimal classic-netCDF reader (CDF-1 / CDF-2 only, big-endian).
//!
//! Implements the Unidata classic-format header walk far enough to extract
//! numeric `wout` variables: magic probe, dimension/variable/attribute
//! lists, and contiguous non-record variable data. HDF5-backed netCDF-4
//! and CDF-5 are rejected at the magic probe ([`probe_variant`]); record
//! (unlimited-dimension) variables are rejected at parse time. No
//! compression, no chunking, no groups — none of those exist in the
//! classic model.

use std::collections::BTreeMap;

use crate::error::{Error, Result, WrongVariant};

/// HDF5 file signature (also the first 8 bytes of every netCDF-4 file).
pub const HDF5_SIGNATURE: [u8; 8] = [0x89, 0x48, 0x44, 0x46, 0x0D, 0x0A, 0x1A, 0x0A];

/// Classic-model variant accepted by this reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassicVariant {
    /// CDF-1 classic (`43 44 46 01`); `ncdump -k` reports `classic`.
    Cdf1,
    /// CDF-2 64-bit offset (`43 44 46 02`); `ncdump -k` reports
    /// `64-bit offset`.
    Cdf2,
}

impl ClassicVariant {
    /// The `ncdump -k` spelling for this variant.
    pub fn ncdump_kind(&self) -> &'static str {
        match self {
            ClassicVariant::Cdf1 => "classic",
            ClassicVariant::Cdf2 => "64-bit offset",
        }
    }

    /// Whether 64-bit header fields (lengths, sizes, offsets) apply.
    fn is_64bit_offset(&self) -> bool {
        matches!(self, ClassicVariant::Cdf2)
    }
}

/// Probe the first bytes of a file: accept CDF-1/CDF-2, reject the HDF5
/// signature (netCDF-4) and CDF-5 loudly, and reject anything else.
/// Needs at least 4 bytes; an 8-byte slice distinguishes HDF5 fully.
pub fn probe_variant(bytes: &[u8]) -> Result<ClassicVariant> {
    if bytes.len() < 4 {
        return Err(Error::Truncated("magic"));
    }
    if bytes.len() >= 8 && bytes[..8] == HDF5_SIGNATURE {
        return Err(Error::WrongVariant(WrongVariant::Hdf5NetCdf4));
    }
    if bytes[..3] != [0x43, 0x44, 0x46] {
        let mut magic = [0u8; 4];
        magic.copy_from_slice(&bytes[..4]);
        return Err(Error::NotClassic { magic });
    }
    match bytes[3] {
        0x01 => Ok(ClassicVariant::Cdf1),
        0x02 => Ok(ClassicVariant::Cdf2),
        0x05 => Err(Error::WrongVariant(WrongVariant::Cdf5)),
        _ => {
            let mut magic = [0u8; 4];
            magic.copy_from_slice(&bytes[..4]);
            Err(Error::NotClassic { magic })
        }
    }
}

/// Probe a file on disk by reading its first 8 bytes.
pub fn probe_file(path: &std::path::Path) -> Result<ClassicVariant> {
    let bytes = std::fs::read(path).map_err(|e| Error::Io(e.to_string()))?;
    probe_variant(&bytes)
}

/// Classic element type codes (1–6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NcType {
    /// 8-bit integer.
    Byte,
    /// Character.
    Char,
    /// 16-bit integer.
    Short,
    /// 32-bit integer.
    Int,
    /// 32-bit float.
    Float,
    /// 64-bit float.
    Double,
}

impl NcType {
    fn from_code(code: u32) -> Result<NcType> {
        match code {
            1 => Ok(NcType::Byte),
            2 => Ok(NcType::Char),
            3 => Ok(NcType::Short),
            4 => Ok(NcType::Int),
            5 => Ok(NcType::Float),
            6 => Ok(NcType::Double),
            other => Err(Error::UnknownType(other)),
        }
    }

    /// Storage width of one element in bytes.
    fn width(&self) -> u64 {
        match self {
            NcType::Byte | NcType::Char => 1,
            NcType::Short => 2,
            NcType::Int | NcType::Float => 4,
            NcType::Double => 8,
        }
    }
}

/// Decoded non-record variable payload.
#[derive(Debug, Clone, PartialEq)]
pub enum Decoded {
    /// Raw bytes (`NC_BYTE`).
    Bytes(Vec<u8>),
    /// Integer elements, widened to `i64` (`NC_BYTE`/`NC_SHORT`/`NC_INT`).
    Ints(Vec<i64>),
    /// Floating-point elements, widened to `f64` (`NC_FLOAT`/`NC_DOUBLE`).
    Floats(Vec<f64>),
    /// Character data (`NC_CHAR`), trailing NUL padding stripped.
    Text(String),
}

impl Decoded {
    /// Copy the payload as `f64` (integers widen exactly; text/bytes fail).
    pub fn as_f64_vec(&self) -> Result<Vec<f64>> {
        match self {
            Decoded::Floats(v) => {
                if v.iter().any(|x| !x.is_finite()) {
                    return Err(Error::BadValue("non-finite float payload".into()));
                }
                Ok(v.clone())
            }
            Decoded::Ints(v) => Ok(v.iter().map(|x| *x as f64).collect()),
            Decoded::Bytes(_) => Err(Error::BadValue("byte payload is not numeric".into())),
            Decoded::Text(_) => Err(Error::BadValue("text payload is not numeric".into())),
        }
    }

    /// Copy the payload as `i64` (floats accepted only when integral).
    pub fn as_i64_vec(&self) -> Result<Vec<i64>> {
        match self {
            Decoded::Ints(v) => Ok(v.clone()),
            Decoded::Floats(v) => v
                .iter()
                .map(|x| {
                    if x.is_finite() && x.fract() == 0.0 {
                        Ok(*x as i64)
                    } else {
                        Err(Error::BadValue("non-integral float payload".into()))
                    }
                })
                .collect(),
            Decoded::Bytes(_) => Err(Error::BadValue("byte payload is not numeric".into())),
            Decoded::Text(_) => Err(Error::BadValue("text payload is not numeric".into())),
        }
    }

    /// Read a scalar payload (exactly one element) as `f64`.
    pub fn as_f64_scalar(&self, name: &'static str) -> Result<f64> {
        let v = self.as_f64_vec()?;
        if v.len() != 1 {
            return Err(Error::BadShape {
                what: name,
                expected: 1,
                got: v.len(),
            });
        }
        Ok(v[0])
    }

    /// Read a scalar payload (exactly one element) as `i64`.
    pub fn as_i64_scalar(&self, name: &'static str) -> Result<i64> {
        let v = self.as_i64_vec()?;
        if v.len() != 1 {
            return Err(Error::BadShape {
                what: name,
                expected: 1,
                got: v.len(),
            });
        }
        Ok(v[0])
    }
}

/// A classic dimension: name plus fixed length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Dim {
    /// Dimension name (`radius`, `mn_mode`, `mn_mode_nyq`, ...).
    pub name: String,
    /// Fixed length (a zero length marks the unlimited dimension, which
    /// classic `wout` files never carry).
    pub len: u64,
}

/// A classic variable with its decoded non-record payload.
#[derive(Debug, Clone, PartialEq)]
pub struct Var {
    /// Variable name.
    pub name: String,
    /// Dimension ids in file order (outer slowest, inner fastest).
    pub dim_ids: Vec<u32>,
    /// Element type.
    pub dtype: NcType,
    /// Decoded payload in C (row-major) order.
    pub data: Decoded,
}

/// A parsed classic file: header plus eager non-record payloads.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassicFile {
    /// CDF-1 or CDF-2.
    pub variant: ClassicVariant,
    /// Record count (zero for classic `wout` files, which carry no
    /// unlimited dimension).
    pub numrecs: u64,
    /// Dimensions in file order.
    pub dims: Vec<Dim>,
    /// Non-record variables in file order.
    pub vars: Vec<Var>,
}

impl ClassicFile {
    /// Find a dimension id by name.
    pub fn dim_id(&self, name: &str) -> Result<usize> {
        self.dims
            .iter()
            .position(|d| d.name == name)
            .ok_or_else(|| Error::MissingDimension(name.to_string()))
    }

    /// Find a variable by name.
    pub fn var(&self, name: &str) -> Result<&Var> {
        self.vars
            .iter()
            .find(|v| v.name == name)
            .ok_or_else(|| Error::MissingVariable(name.to_string()))
    }

    /// Dimension lengths keyed by name.
    pub fn dim_lens(&self) -> BTreeMap<&str, u64> {
        self.dims.iter().map(|d| (d.name.as_str(), d.len)).collect()
    }
}

/// Parse a whole classic file from memory.
pub fn parse(bytes: &[u8]) -> Result<ClassicFile> {
    let variant = probe_variant(bytes)?;
    let mut cur = Cursor::new(bytes, 4);
    let numrecs = if variant.is_64bit_offset() {
        cur.u64("numrecs")?
    } else {
        u64::from(cur.u32("numrecs")?)
    };
    let dims = read_dims(&mut cur, variant)?;
    read_attrs(&mut cur, variant)?;
    let raw_vars = read_var_headers(&mut cur, variant)?;
    // A zero-length leading dimension marks the unlimited dimension;
    // classic `wout` files carry none, and record layout is out of scope.
    let unlimited: Option<u32> = dims.iter().position(|d| d.len == 0).map(|i| i as u32);
    let mut vars = Vec::with_capacity(raw_vars.len());
    for raw in raw_vars {
        let is_record = match unlimited {
            Some(uid) => !raw.dim_ids.is_empty() && raw.dim_ids[0] == uid,
            None => false,
        };
        if is_record {
            return Err(Error::RecordVariable(raw.name));
        }
        let data = read_var_data(bytes, &raw, &dims)?;
        vars.push(Var {
            name: raw.name,
            dim_ids: raw.dim_ids,
            dtype: raw.dtype,
            data,
        });
    }
    Ok(ClassicFile {
        variant,
        numrecs,
        dims,
        vars,
    })
}

/// Read a classic file from disk.
pub fn read_file(path: &std::path::Path) -> Result<ClassicFile> {
    let bytes = std::fs::read(path).map_err(|e| Error::Io(e.to_string()))?;
    parse(&bytes)
}

// ---------------------------------------------------------------------------
// Header walk
// ---------------------------------------------------------------------------

const TAG_DIMENSION: u32 = 10;
const TAG_VARIABLE: u32 = 11;
const TAG_ATTRIBUTE: u32 = 12;

struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8], pos: usize) -> Cursor<'a> {
        Cursor { bytes, pos }
    }

    fn take(&mut self, n: usize, what: &'static str) -> Result<&'a [u8]> {
        let end = self.pos.saturating_add(n);
        if end > self.bytes.len() {
            return Err(Error::Truncated(what));
        }
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u32(&mut self, what: &'static str) -> Result<u32> {
        let b = self.take(4, what)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self, what: &'static str) -> Result<u64> {
        let b = self.take(8, what)?;
        Ok(u64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    fn xlen(&mut self, variant: ClassicVariant, what: &'static str) -> Result<u64> {
        if variant.is_64bit_offset() {
            self.u64(what)
        } else {
            Ok(u64::from(self.u32(what)?))
        }
    }

    /// Padded name: u32 length + bytes, advanced to a 4-byte boundary.
    fn name(&mut self) -> Result<String> {
        let len = self.u32("name length")? as usize;
        let raw = self.take(len, "name bytes")?;
        let pad = (4 - (len % 4)) % 4;
        if pad > 0 {
            self.take(pad, "name padding")?;
        }
        String::from_utf8(raw.to_vec()).map_err(|_| Error::BadHeader("non-UTF8 name in header"))
    }
}

fn read_list_tag(cur: &mut Cursor<'_>, expected: u32, what: &'static str) -> Result<Option<u32>> {
    let tag = cur.u32(what)?;
    let count = cur.u32(what)?;
    if tag == 0 && count == 0 {
        return Ok(None);
    }
    if tag != expected {
        return Err(Error::BadHeader(what));
    }
    Ok(Some(count))
}

fn read_dims(cur: &mut Cursor<'_>, variant: ClassicVariant) -> Result<Vec<Dim>> {
    let Some(count) = read_list_tag(cur, TAG_DIMENSION, "dim list")? else {
        return Ok(Vec::new());
    };
    let mut dims = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let name = cur.name()?;
        let len = cur.xlen(variant, "dim length")?;
        dims.push(Dim { name, len });
    }
    Ok(dims)
}

/// Skip one attribute value run (the reader keeps no attributes).
fn skip_attr_values(cur: &mut Cursor<'_>, dtype: NcType, nelems: u64) -> Result<()> {
    let bytes = nelems.saturating_mul(dtype.width());
    let padded = bytes.saturating_add(3) & !3;
    let n: usize = padded
        .try_into()
        .map_err(|_| Error::BadHeader("attribute run too large"))?;
    cur.take(n, "attribute values")?;
    Ok(())
}

fn read_attrs(cur: &mut Cursor<'_>, variant: ClassicVariant) -> Result<()> {
    let Some(count) = read_list_tag(cur, TAG_ATTRIBUTE, "attribute list")? else {
        return Ok(());
    };
    for _ in 0..count {
        let _name = cur.name()?;
        let code = cur.u32("attribute type")?;
        let dtype = NcType::from_code(code)?;
        // Attribute element counts stay 32-bit even in CDF-2.
        let nelems = u64::from(cur.u32("attribute length")?);
        let _ = variant;
        skip_attr_values(cur, dtype, nelems)?;
    }
    Ok(())
}

struct RawVar {
    name: String,
    dim_ids: Vec<u32>,
    dtype: NcType,
    vsize: u64,
    begin: u64,
}

fn read_var_headers(cur: &mut Cursor<'_>, variant: ClassicVariant) -> Result<Vec<RawVar>> {
    let Some(count) = read_list_tag(cur, TAG_VARIABLE, "variable list")? else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let name = cur.name()?;
        let ndims = cur.u32("var ndims")?;
        if ndims > 32 {
            return Err(Error::BadHeader("var ndims out of range"));
        }
        let mut dim_ids = Vec::with_capacity(ndims as usize);
        for _ in 0..ndims {
            dim_ids.push(cur.u32("var dim id")?);
        }
        read_attrs(cur, variant)?;
        let code = cur.u32("var type")?;
        let dtype = NcType::from_code(code)?;
        let vsize = cur.xlen(variant, "var size")?;
        let begin = cur.xlen(variant, "var offset")?;
        out.push(RawVar {
            name,
            dim_ids,
            dtype,
            vsize,
            begin,
        });
    }
    Ok(out)
}

fn read_var_data(bytes: &[u8], raw: &RawVar, dims: &[Dim]) -> Result<Decoded> {
    // Expected payload: product of fixed dimension lengths times width.
    // (Record variables are rejected before this point.)
    let mut nelems: u64 = 1;
    for id in &raw.dim_ids {
        let dim = dims
            .get(*id as usize)
            .ok_or(Error::BadHeader("var dim id out of range"))?;
        nelems = nelems
            .checked_mul(dim.len)
            .ok_or(Error::BadHeader("var shape overflows"))?;
    }
    let want = nelems
        .checked_mul(raw.dtype.width())
        .ok_or(Error::BadHeader("var payload overflows"))?;
    if want != raw.vsize {
        return Err(Error::BadHeader("var size disagrees with shape"));
    }
    let begin: usize = raw
        .begin
        .try_into()
        .map_err(|_| Error::BadHeader("var offset out of range"))?;
    let len: usize = want
        .try_into()
        .map_err(|_| Error::BadHeader("var payload too large"))?;
    let end = begin.saturating_add(len);
    if end > bytes.len() {
        return Err(Error::Truncated("variable data"));
    }
    decode(&bytes[begin..end], raw.dtype, nelems)
}

fn decode(raw: &[u8], dtype: NcType, nelems: u64) -> Result<Decoded> {
    let n: usize = nelems
        .try_into()
        .map_err(|_| Error::BadHeader("var payload too large"))?;
    match dtype {
        NcType::Char => {
            let mut text = String::from_utf8(raw.to_vec())
                .map_err(|_| Error::BadValue("non-UTF8 char payload".into()))?;
            while text.ends_with('\0') {
                text.pop();
            }
            Ok(Decoded::Text(text))
        }
        NcType::Byte => {
            if raw.len() != n {
                return Err(Error::BadHeader("byte payload length"));
            }
            Ok(Decoded::Ints(raw.iter().map(|b| *b as i8 as i64).collect()))
        }
        NcType::Short => {
            if raw.len() != n * 2 {
                return Err(Error::BadHeader("short payload length"));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                out.push(i16::from_be_bytes([raw[2 * i], raw[2 * i + 1]]) as i64);
            }
            Ok(Decoded::Ints(out))
        }
        NcType::Int => {
            if raw.len() != n * 4 {
                return Err(Error::BadHeader("int payload length"));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                out.push(i32::from_be_bytes([
                    raw[4 * i],
                    raw[4 * i + 1],
                    raw[4 * i + 2],
                    raw[4 * i + 3],
                ]) as i64);
            }
            Ok(Decoded::Ints(out))
        }
        NcType::Float => {
            if raw.len() != n * 4 {
                return Err(Error::BadHeader("float payload length"));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                out.push(f32::from_be_bytes([
                    raw[4 * i],
                    raw[4 * i + 1],
                    raw[4 * i + 2],
                    raw[4 * i + 3],
                ]) as f64);
            }
            Ok(Decoded::Floats(out))
        }
        NcType::Double => {
            if raw.len() != n * 8 {
                return Err(Error::BadHeader("double payload length"));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                let mut b = [0u8; 8];
                b.copy_from_slice(&raw[8 * i..8 * i + 8]);
                out.push(f64::from_be_bytes(b));
            }
            Ok(Decoded::Floats(out))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probe_accepts_cdf1_and_cdf2() {
        assert_eq!(
            probe_variant(&[0x43, 0x44, 0x46, 0x01]).unwrap(),
            ClassicVariant::Cdf1
        );
        assert_eq!(
            probe_variant(&[0x43, 0x44, 0x46, 0x02, 0, 0, 0, 0]).unwrap(),
            ClassicVariant::Cdf2
        );
        assert_eq!(ClassicVariant::Cdf1.ncdump_kind(), "classic");
        assert_eq!(ClassicVariant::Cdf2.ncdump_kind(), "64-bit offset");
    }

    #[test]
    fn probe_rejects_hdf5_cdf5_and_garbage() {
        assert_eq!(
            probe_variant(&HDF5_SIGNATURE),
            Err(Error::WrongVariant(WrongVariant::Hdf5NetCdf4))
        );
        assert_eq!(
            probe_variant(&[0x43, 0x44, 0x46, 0x05]),
            Err(Error::WrongVariant(WrongVariant::Cdf5))
        );
        assert!(matches!(
            probe_variant(&[0x42, 0x41, 0x44, 0x21]),
            Err(Error::NotClassic { .. })
        ));
        assert_eq!(probe_variant(&[0x43, 0x44]), Err(Error::Truncated("magic")));
        // Unknown CDF version byte with a CDF magic.
        assert!(matches!(
            probe_variant(&[0x43, 0x44, 0x46, 0x09]),
            Err(Error::NotClassic { .. })
        ));
    }
}
