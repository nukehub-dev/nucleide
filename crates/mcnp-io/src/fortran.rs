//! Shared Fortran unformatted-record framework (`_FortranRecord` /
//! `_BinaryReader.put/get_fortran_record` in `pyne.binaryreader`).
//!
//! Every record on disk is framed as `[i32 len][payload][i32 len]`
//! little-endian; a mismatch between the leading and trailing markers is a
//! typed [`Error::BadRecordMarker`] (the upstream error path itself
//! `AttributeError`s on `num_bytes2` — that bug is deliberately not ported).
//!
//! Scope notes (verified against `pyne/binaryreader.py:175-208`):
//! - Upstream uses the native `struct` `"i"` format with no endianness
//!   parameter; all vendored fixtures are little-endian, so this framework
//!   is little-endian only. There is no big-endian fixture, so a big-endian
//!   real-bytes oracle is UNVERIFIED — big-endian input fails loudly at the
//!   marker check (or via the PTRAC [`crate::ptrac`] `Unsupported` path).
//! - `SSW` ([`crate::surfsrc`]) framing delegates to [`read_record`] /
//!   [`frame_record`]; PTRAC keeps its own width-switching cursor on top of
//!   the same framing contract.
//!
//! Parity quirks documented, not fixed:
//! - `SF_00001` (MCNP6) splits the SSW header over two records; otherwise one
//!   record holds kod/ver/loddat.
//! - A negative stored `np1` (`orignp1 < 0`) signals the extra table-2 record
//!   (cells/particle/macrobody-facet info); the writer emits table 2 only in
//!   that case, unlike the legacy writer which always emits it.
//! - Track record width is `abs(ncrd)` doubles; `nps` shifts in
//!   `combine_files` preserve sign via `copysign`.
//! - PTRAC input-echo repack: 20 f32 halves reinterpreted as 10 f64 flips the
//!   file into 8-byte mode.
//! - Meshtal single-group totals mirror the lone energy group; the column
//!   header normalizes `Rel ` to `Rel_`.
//! - WWINP reproduces Python's `{0:13.5E}` field formatting byte-for-byte.

use std::fmt;
use std::io::{Read, Write};

/// Errors raised while framing or traversing Fortran records.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// I/O failure while reading or writing the stream.
    Io(String),
    /// Leading/trailing record-length markers disagreed.
    BadRecordMarker {
        /// Leading record-length marker.
        lead: i32,
        /// Trailing record-length marker.
        trailer: i32,
    },
    /// Payload exhausted mid-field (`need` bytes wanted, `left` remain).
    ShortRecord {
        /// Bytes wanted.
        need: usize,
        /// Bytes actually present.
        left: usize,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadRecordMarker { lead, trailer } => {
                write!(f, "record markers disagree: {lead} vs {trailer}")
            }
            Error::ShortRecord { need, left } => {
                write!(f, "record too short: need {need} bytes, {left} remain")
            }
        }
    }
}

impl std::error::Error for Error {}

/// One Fortran unformatted record payload with a sequential read cursor
/// (`_FortranRecord.get_*`).
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    bytes: Vec<u8>,
    pos: usize,
}

impl Record {
    /// Wrap an already-framed payload (see [`read_record`]).
    pub fn new(bytes: Vec<u8>) -> Self {
        Record { bytes, pos: 0 }
    }

    /// Raw payload length in bytes.
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// True when the payload holds no bytes.
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Bytes still unread.
    pub fn remaining(&self) -> usize {
        self.bytes.len() - self.pos
    }

    /// Rewind the cursor to the start of the payload.
    pub fn reset(&mut self) {
        self.pos = 0;
    }

    fn take(&mut self, n: usize) -> Result<&[u8], Error> {
        if self.pos + n > self.bytes.len() {
            return Err(Error::ShortRecord {
                need: self.pos + n,
                left: self.bytes.len(),
            });
        }
        let s = &self.bytes[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    /// Fixed-width string field (`get_string`).
    pub fn get_string(&mut self, n: usize) -> String {
        String::from_utf8_lossy(self.take(n).unwrap_or_default()).into_owned()
    }

    /// One little-endian 4-byte integer (`get_int`).
    pub fn get_i32(&mut self) -> Result<i32, Error> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// `n` little-endian 4-byte integers.
    pub fn get_i32_n(&mut self, n: usize) -> Result<Vec<i32>, Error> {
        (0..n).map(|_| self.get_i32()).collect()
    }

    /// One little-endian 8-byte integer (`get_long`).
    pub fn get_i64(&mut self) -> Result<i64, Error> {
        let b = self.take(8)?;
        Ok(i64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    /// One little-endian 4-byte float (`get_float`).
    pub fn get_f32(&mut self) -> Result<f32, Error> {
        let b = self.take(4)?;
        Ok(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// One little-endian 8-byte float (`get_double`).
    pub fn get_f64(&mut self) -> Result<f64, Error> {
        let b = self.take(8)?;
        Ok(f64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    /// `n` little-endian 8-byte floats.
    pub fn get_f64_n(&mut self, n: usize) -> Result<Vec<f64>, Error> {
        (0..n).map(|_| self.get_f64()).collect()
    }

    /// Drain every remaining whole 4-byte integer (trailing extras).
    pub fn drain_i32_extras(&mut self) -> Result<Vec<i32>, Error> {
        let mut v = Vec::new();
        while self.remaining() >= 4 {
            v.push(self.get_i32()?);
        }
        Ok(v)
    }
}

/// Accumulates a single record payload (`_FortranRecord.put_*`).
#[derive(Debug, Clone, Default)]
pub struct RecordSink {
    bytes: Vec<u8>,
}

impl RecordSink {
    /// Empty payload accumulator.
    pub fn new() -> Self {
        RecordSink { bytes: Vec::new() }
    }

    /// Raw string bytes (no length prefix, no padding).
    pub fn put_str(&mut self, s: &str) {
        self.bytes.extend_from_slice(s.as_bytes());
    }

    /// `put_int`: little-endian i32.
    pub fn put_int(&mut self, v: i32) {
        self.bytes.extend_from_slice(&v.to_le_bytes());
    }

    /// `put_long`: little-endian i64.
    pub fn put_long(&mut self, v: i64) {
        self.bytes.extend_from_slice(&v.to_le_bytes());
    }

    /// `put_float`: little-endian f32.
    pub fn put_float(&mut self, v: f32) {
        self.bytes.extend_from_slice(&v.to_le_bytes());
    }

    /// `put_double`: little-endian f64.
    pub fn put_double(&mut self, v: f64) {
        self.bytes.extend_from_slice(&v.to_le_bytes());
    }

    /// Payload length so far.
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// True when nothing has been appended yet.
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Borrow the raw payload bytes.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Frame the payload as `[i32 len][payload][i32 len]`.
    pub fn frame(self) -> Vec<u8> {
        frame_record(&self.bytes)
    }
}

/// Frame `payload` as `[i32 len][payload][i32 len]` little-endian
/// (`put_fortran_record`).
pub fn frame_record(payload: &[u8]) -> Vec<u8> {
    let len = i32::try_from(payload.len())
        .unwrap_or(i32::MAX)
        .to_le_bytes();
    let mut out = Vec::with_capacity(payload.len() + 8);
    out.extend_from_slice(&len);
    out.extend_from_slice(payload);
    out.extend_from_slice(&len);
    out
}

/// Read one `[i32 len][payload][i32 len]` record (`get_fortran_record`).
pub fn read_record<R: Read>(r: &mut R) -> Result<Record, Error> {
    let mut marker = [0u8; 4];
    r.read_exact(&mut marker)
        .map_err(|e| Error::Io(e.to_string()))?;
    let lead = i32::from_le_bytes(marker);
    if lead < 0 {
        return Err(Error::BadRecordMarker {
            lead,
            trailer: lead,
        });
    }
    let mut payload = Vec::with_capacity(lead as usize);
    // `read_to_end` (not `read_exact`) so a truncated payload reports how
    // many bytes were actually present instead of an opaque I/O error.
    let got = r
        .by_ref()
        .take(lead as u64)
        .read_to_end(&mut payload)
        .map_err(|e| Error::Io(e.to_string()))?;
    if got != lead as usize {
        return Err(Error::ShortRecord {
            need: lead as usize,
            left: got,
        });
    }
    r.read_exact(&mut marker)
        .map_err(|e| Error::Io(e.to_string()))?;
    let trailer = i32::from_le_bytes(marker);
    if lead != trailer {
        return Err(Error::BadRecordMarker { lead, trailer });
    }
    Ok(Record::new(payload))
}

/// Write one framed record to a stream.
pub fn write_record<W: Write>(w: &mut W, payload: &[u8]) -> Result<(), Error> {
    w.write_all(&frame_record(payload))
        .map_err(|e| Error::Io(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_round_trip_preserves_payload() {
        let payload = b"Hello World!".to_vec();
        let framed = frame_record(&payload);
        assert_eq!(framed.len(), payload.len() + 8);
        let mut cursor = std::io::Cursor::new(&framed);
        let rec = read_record(&mut cursor).unwrap();
        assert_eq!(rec.bytes, payload);
    }

    #[test]
    fn empty_payload_frames_to_zero_markers() {
        let framed = frame_record(&[]);
        assert_eq!(framed, vec![0, 0, 0, 0, 0, 0, 0, 0]);
        let mut cursor = std::io::Cursor::new(&framed);
        let rec = read_record(&mut cursor).unwrap();
        assert!(rec.is_empty());
    }

    #[test]
    fn marker_mismatch_is_typed_not_attribute_error() {
        // Upstream raises ValueError through a path that itself
        // AttributeErrors (`self.num_bytes2`); here the disagreement is a
        // typed error carrying both markers.
        let mut framed = frame_record(b"abcd");
        let last = framed.len() - 1;
        framed[last] = framed[last].wrapping_add(1);
        let mut cursor = std::io::Cursor::new(&framed);
        // Trailer bytes went from [04 00 00 00] to [04 00 00 01].
        assert_eq!(
            read_record(&mut cursor),
            Err(Error::BadRecordMarker {
                lead: 4,
                trailer: 4 + (1 << 24),
            })
        );
    }

    #[test]
    fn negative_lead_marker_is_rejected() {
        let mut bad = (-3i32).to_le_bytes().to_vec();
        bad.extend_from_slice(&[0u8; 4]);
        let mut cursor = std::io::Cursor::new(&bad);
        assert_eq!(
            read_record(&mut cursor),
            Err(Error::BadRecordMarker {
                lead: -3,
                trailer: -3,
            })
        );
        // Empty input fails on the leading marker read itself.
        let mut cursor = std::io::Cursor::new(Vec::new());
        assert!(matches!(read_record(&mut cursor), Err(Error::Io(_))));
    }

    #[test]
    fn truncated_payload_reports_need_and_left() {
        let mut framed = frame_record(b"abcdefgh");
        framed.truncate(4 + 5); // markers + 5 of 8 payload bytes, no trailer
        let mut cursor = std::io::Cursor::new(&framed);
        assert_eq!(
            read_record(&mut cursor),
            Err(Error::ShortRecord { need: 8, left: 5 })
        );
    }

    #[test]
    // 3.14 is the literal from the PyNE oracle vector (its f32 bytes are
    // pinned below); it must not become `PI`.
    #[allow(clippy::approx_constant)]
    fn mixed_record_matches_pyne_oracle_bytes() {
        // Byte vector from pyne `test_binaryreader.py::test_write_FR_mixed_record`:
        // put_int(8), put_string("Hello World!"), put_double([1.6e-19, 6.02e23]),
        // put_float(3.14) over 4 + 12 + 16 + 4 = 36 payload bytes.
        let expected: Vec<u8> = b"\x08\x00\x00\x00Hello World!#B\x92\x0c\xa1\x9c\x07<a\xd3\
            \xa8\x10\x9f\xde\xdfD\xc3\xf5H@"
            .to_vec();
        let mut sink = RecordSink::new();
        sink.put_int(8);
        sink.put_str("Hello World!");
        sink.put_double(1.6e-19);
        sink.put_double(6.02e23);
        sink.put_float(3.14);
        assert_eq!(sink.len(), 36);
        assert_eq!(sink.bytes(), &expected[..]);

        // And the same bytes decode back through the cursor.
        let mut rec = Record::new(expected);
        assert_eq!(rec.get_i32().unwrap(), 8);
        assert_eq!(rec.get_string(12), "Hello World!");
        let doubles = rec.get_f64_n(2).unwrap();
        assert!((doubles[0] - 1.6e-19).abs() / 1.6e-19 < 1e-12);
        assert!((doubles[1] - 6.02e23).abs() / 6.02e23 < 1e-12);
        assert!((rec.get_f32().unwrap() - 3.14).abs() < 1e-6);
        assert_eq!(rec.remaining(), 0);
    }

    #[test]
    fn read_past_end_is_short_record() {
        let mut rec = Record::new(vec![1, 2, 3]);
        assert_eq!(rec.get_i32(), Err(Error::ShortRecord { need: 4, left: 3 }));
    }

    #[test]
    fn write_record_frames_to_stream() {
        let mut out = Vec::new();
        write_record(&mut out, b"xyz").unwrap();
        assert_eq!(&out[4..7], b"xyz");
        let mut cursor = std::io::Cursor::new(&out);
        assert_eq!(read_record(&mut cursor).unwrap().bytes, b"xyz");
    }

    #[test]
    fn errors_display_and_accessors() {
        assert!(format!("{}", Error::Io("x".into())).contains("io error"));
        assert!(format!(
            "{}",
            Error::BadRecordMarker {
                lead: 1,
                trailer: 2
            }
        )
        .contains("1 vs 2"));
        assert!(format!("{}", Error::ShortRecord { need: 4, left: 1 }).contains("need 4"));
        let mut sink = RecordSink::new();
        assert!(sink.is_empty());
        assert_eq!(sink.len(), 0);
        sink.put_int(-7);
        sink.put_long(-123456789012345);
        assert!(!sink.is_empty());
        assert_eq!(sink.len(), 12);
        assert_eq!(sink.bytes().len(), 12);
        let mut rec = Record::new(sink.bytes().to_vec());
        assert_eq!(rec.len(), 12);
        assert_eq!(rec.remaining(), 12);
        assert_eq!(rec.get_i32_n(1).unwrap(), vec![-7]);
        assert_eq!(rec.get_i64().unwrap(), -123456789012345);
        assert_eq!(rec.remaining(), 0);
        rec.reset();
        assert_eq!(rec.remaining(), 12);
    }
}
