//! MCPL Monte Carlo Particle List read/write.
//!
//! Implements the binary format described in Kittelmann et al.,
//! *Monte Carlo Particle Lists: MCPL*, Computer Physics Communications 218
//! (2017) pp. 17–42 (<https://doi.org/10.1016/j.cpc.2017.04.012>), full text
//! at <https://mctools.github.io/mcpl/mcpl.pdf>, format summary at
//! <https://mctools.github.io/mcpl/format/> (upstream project at
//! <https://github.com/mctools/mcpl>, Apache-2.0): format versions 2 and 3,
//! the logical particle fields (kinetic energy in MeV, polarisation vector,
//! position in cm, unit direction, time in ms, weight, PDG code, user flags),
//! and the per-file option flags (userflags / polarisation /
//! single-vs-double precision / universal PDG / universal weight) with
//! little-endian storage. No upstream code is vendored here.
//!
//! Validation status: byte-exact round-trips on hand-built synthetic records
//! (closed-form packing math), `stat:sum` syntax validated per the upstream
//! rules, plus record-level interop cross-checked both directions against
//! the upstream 2.2.8 implementation (upstream-C-written files read here;
//! files written here read by the upstream Python reader and `mcpltool`).
//! The container harness keeps a record-level cross-check
//! (`validation/mcpl_vs_refs.py`, loud SKIP when the oracle is absent).
//!
//! ## Wire layout (v1 scope, all verified)
//!
//! ```text
//! magic      4 bytes  "MCPL"
//! version    3 bytes  ASCII "002" | "003" (we write "003")
//! endian     1 byte   'L' (little-endian; 'B' files are rejected)
//! nparticles u64 LE   particle count (patched at write time)
//! options    8 x u32 LE: ncomments, nblobs, has_userflags,
//!              has_polarisation, single_prec, universal_pdgcode bits,
//!              particle_size, has_universal_weight
//! uweight    f64 LE, present iff has_universal_weight != 0
//! srcname    u32 LE length + bytes
//! comments   ncomments x (u32 LE length + bytes)
//! blob keys  nblobs x (u32 LE length + bytes)
//! blobs      nblobs x (u32 LE length + raw bytes)
//! particles  nparticles records, each:
//!              [polarisation 3 x fp]? position 3 x fp,
//!              packed dir+ekin 3 x fp, time 1 x fp,
//!              [weight 1 x fp]? [pdgcode i32]? [userflags u32]?
//!              (fp = f32 when single_prec != 0, else f64)
//! ```
//!
//! Direction/energy packing (format 3, "Adaptive Projection Packing"): the
//! unit direction is reduced to two floats plus one sign bit, and the kinetic
//! energy rides in that sign slot via `copysign(ekin, sign)`. Format 2
//! ("Octahedral Packing") is read-only; we always write format 3.
//!
//! ## Named-open items (explicitly out of v1 scope)
//!
//! - Big-endian (`'B'`) files: rejected with [`Error::UnsupportedEndianness`]
//!   (upstream byteswaps; this reader stays little-endian only).
//! - `stat:sum` merge/scale semantics (files merged by third-party tools):
//!   such comments round-trip verbatim and their syntax is validated
//!   ([`statsum_validate`], [`statsum_comment`]) but values are never
//!   interpreted.
//! - SSW↔MCPL conversion beyond the neutron/gamma-only v1 in [`ssw`]:
//!   other particle kinds are named errors, never silent skips.
//! - `.gz` compression levels: gzip transport is transparent, but compressed
//!   bytes are never asserted (encoder settings differ across writers).

pub mod ssw;

use std::io::{Read, Write};
use std::path::Path;

use thiserror::Error;

/// Format version written by this crate (current per the paper/upstream).
pub const FORMAT_VERSION_WRITE: u16 = 3;
/// Lowest format version readable by this crate.
pub const FORMAT_VERSION_MIN: u16 = 2;
/// Magic bytes opening every MCPL file.
pub const MAGIC: [u8; 4] = *b"MCPL";
/// Reader tolerance for unit directions (matches upstream `1.0e-5`).
pub const UNIT_TOL: f64 = 1.0e-5;

/// Errors raised while reading or writing MCPL data.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum Error {
    /// Underlying filesystem, gzip, or stream failure with context.
    #[error("io error: {0}")]
    Io(String),
    /// First four bytes are not the `MCPL` magic.
    #[error("not an MCPL file (bad magic)")]
    BadMagic,
    /// Three-digit format version is not 2 or 3.
    #[error("unsupported MCPL format version {0} (supported: 2, 3)")]
    UnsupportedVersion(u16),
    /// Big-endian files are detected but not supported.
    #[error("big-endian MCPL files are not supported (endian flag '{0}')")]
    UnsupportedEndianness(char),
    /// Byte stream ends inside the named header or particle region.
    #[error("truncated file while reading {0}")]
    Truncated(String),
    /// A length-prefixed string field is not valid UTF-8.
    #[error("string field {0} holds non-UTF8 bytes")]
    InvalidUtf8(String),
    /// A length-prefixed string field contains an embedded NUL byte.
    #[error("string field {0} contains an embedded NUL byte")]
    EmbeddedNul(String),
    /// Stored `particle_size` disagrees with the size implied by the flags.
    #[error("particle record size mismatch: header says {found}, flags imply {expected}")]
    ParticleSizeMismatch {
        /// Record size implied by the option flags.
        expected: usize,
        /// Stored `particle_size` from the file header.
        found: usize,
    },
    /// Bytes remain after the declared particle records.
    #[error("trailing bytes: file has {found} bytes, {particle} records end at {expected}")]
    TrailingBytes {
        /// Offset where the declared records end.
        expected: usize,
        /// Total uncompressed file length.
        found: usize,
        /// Declared particle count.
        particle: u64,
    },
    /// A decoded direction deviates from unit length beyond [`UNIT_TOL`].
    #[error("particle {0} direction is not a unit vector (|d|^2 = {1})")]
    NonUnitDirection(usize, f64),
    /// A particle carries a negative kinetic energy.
    #[error("particle {0} has negative kinetic energy {1}")]
    NegativeEnergy(usize, f64),
    /// A file-wide weight must be positive and finite.
    #[error("universal weight must be positive and finite, got {0}")]
    BadUniversalWeight(f64),
    /// A file-wide PDG code must be non-zero.
    #[error("universal PDG code must be non-zero")]
    BadUniversalPdg,
    /// A `stat:sum:` comment fails the upstream syntax
    /// (`stat:sum:<key>:<24-char value>` with key
    /// `[a-zA-Z][a-zA-Z0-9_]*` and a finite `-1` or `>= 0` value).
    /// Upstream rejects such files on open and refuses to write such
    /// comments; this crate does the same instead of parsing leniently.
    #[error("invalid stat:sum comment {index}: {reason}")]
    BadStatSum {
        /// Position of the offending comment in the header.
        index: usize,
        /// Why the comment is malformed.
        reason: String,
    },
    /// Two header blobs share a key (upstream rejects this on write).
    #[error("duplicate blob key {0:?}")]
    DuplicateBlobKey(String),
    /// Two `stat:sum:` comments carry the same key (upstream errors at
    /// open time).
    #[error("comment {0} duplicates stat:sum key {1:?}")]
    DuplicateStatSum(usize, String),
}

impl Error {
    fn truncated(what: impl Into<String>) -> Self {
        Error::Truncated(what.into())
    }
}

/// One logical particle.
///
/// Units follow the paper and the `mcpl.h` API: kinetic energy in MeV,
/// position in cm, time in ms. `polarisation` is meaningful only when the
/// file header carries it (otherwise zeros); `weight`/`pdgcode`/`userflags`
/// are per-particle unless the header fixes them file-wide (universal).
#[derive(Debug, Clone, PartialEq)]
pub struct Particle {
    /// Kinetic energy in MeV.
    pub ekin: f64,
    /// Polarisation vector (meaningful only when the header carries it).
    pub polarisation: [f64; 3],
    /// Position in cm.
    pub position: [f64; 3],
    /// Unit direction vector.
    pub direction: [f64; 3],
    /// Time in ms.
    pub time: f64,
    /// Statistical weight (per-particle unless the header fixes it file-wide).
    pub weight: f64,
    /// PDG particle code (per-particle unless the header fixes it file-wide).
    pub pdgcode: i32,
    /// User flags (zero unless the header stores them).
    pub userflags: u32,
}

impl Particle {
    /// A neutron at rest at the origin (unit direction +z).
    pub fn neutron() -> Self {
        Particle {
            ekin: 0.0,
            polarisation: [0.0; 3],
            position: [0.0; 3],
            direction: [0.0, 0.0, 1.0],
            time: 0.0,
            weight: 1.0,
            pdgcode: 2112,
            userflags: 0,
        }
    }
}

/// A named binary blob from the file header.
#[derive(Debug, Clone, PartialEq)]
pub struct Blob {
    /// Blob key string.
    pub key: String,
    /// Raw blob bytes.
    pub data: Vec<u8>,
}

/// Parsed MCPL header block.
///
/// `version` is 2 or 3 on read; writers always emit 3. `nparticles` is the
/// stored count on read; on write the particle slice length wins (see
/// [`write_to`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    /// Format version: 2 or 3 on read; writers always emit 3.
    pub version: u16,
    /// Whether per-particle user flags are stored.
    pub has_userflags: bool,
    /// Whether per-particle polarisation vectors are stored.
    pub has_polarisation: bool,
    /// `false` (default) = single precision, `true` = double precision.
    pub double_prec: bool,
    /// `None` = per-particle PDG codes; `Some(v != 0)` = file-wide code.
    pub universal_pdgcode: Option<i32>,
    /// `None` = per-particle weights; `Some(w)` = file-wide weight.
    pub universal_weight: Option<f64>,
    /// Source-name string from the header.
    pub srcname: String,
    /// Header comment strings (round-tripped verbatim, never interpreted).
    pub comments: Vec<String>,
    /// Header binary blobs.
    pub blobs: Vec<Blob>,
    /// Stored particle count on read; particle slice length wins on write.
    pub nparticles: u64,
}

impl Default for Header {
    fn default() -> Self {
        Header {
            version: FORMAT_VERSION_WRITE,
            has_userflags: false,
            has_polarisation: false,
            double_prec: false,
            universal_pdgcode: None,
            universal_weight: None,
            srcname: "unknown".to_string(),
            comments: Vec::new(),
            blobs: Vec::new(),
            nparticles: 0,
        }
    }
}

impl Header {
    /// Byte size of one particle record implied by the option flags.
    pub fn particle_size(&self) -> usize {
        let fp = if self.double_prec { 8 } else { 4 };
        let mut n = 7 * fp;
        if self.has_polarisation {
            n += 3 * fp;
        }
        if self.universal_pdgcode.is_none() {
            n += 4;
        }
        if self.universal_weight.is_none() {
            n += fp;
        }
        if self.has_userflags {
            n += 4;
        }
        n
    }

    fn validate_for_write(&self) -> Result<(), Error> {
        if let Some(w) = self.universal_weight {
            if !(w.is_finite() && w > 0.0) {
                return Err(Error::BadUniversalWeight(w));
            }
        }
        if self.universal_pdgcode == Some(0) {
            return Err(Error::BadUniversalPdg);
        }
        validate_statsum_comments(&self.comments)?;
        let mut seen: Vec<&str> = Vec::new();
        for b in &self.blobs {
            if seen.contains(&b.key.as_str()) {
                return Err(Error::DuplicateBlobKey(b.key.clone()));
            }
            seen.push(b.key.as_str());
        }
        Ok(())
    }

    /// Serialize just the header block (everything up to the first particle).
    pub fn header_block(&self) -> Result<Vec<u8>, Error> {
        self.validate_for_write()?;
        let mut out = Vec::new();
        out.extend_from_slice(&MAGIC);
        let v = FORMAT_VERSION_WRITE;
        out.push(b'0' + (v / 100) as u8);
        out.push(b'0' + ((v / 10) % 10) as u8);
        out.push(b'0' + (v % 10) as u8);
        out.push(b'L');
        out.extend_from_slice(&self.nparticles.to_le_bytes());
        let opts = [
            self.comments.len() as u32,
            self.blobs.len() as u32,
            u32::from(self.has_userflags),
            u32::from(self.has_polarisation),
            u32::from(!self.double_prec),
            self.universal_pdgcode.unwrap_or(0) as u32,
            self.particle_size() as u32,
            u32::from(self.universal_weight.is_some()),
        ];
        for o in opts {
            out.extend_from_slice(&o.to_le_bytes());
        }
        if let Some(w) = self.universal_weight {
            out.extend_from_slice(&w.to_le_bytes());
        }
        put_string(&mut out, "srcname", &self.srcname)?;
        for (i, c) in self.comments.iter().enumerate() {
            put_string(&mut out, &format!("comment {i}"), c)?;
        }
        for b in &self.blobs {
            put_string(&mut out, "blob key", &b.key)?;
        }
        for b in &self.blobs {
            let len = u32::try_from(b.data.len())
                .map_err(|_| Error::truncated("blob larger than 4 GiB"))?;
            out.extend_from_slice(&len.to_le_bytes());
            out.extend_from_slice(&b.data);
        }
        Ok(out)
    }
}

fn put_string(out: &mut Vec<u8>, what: &str, s: &str) -> Result<(), Error> {
    if s.as_bytes().contains(&0) {
        return Err(Error::EmbeddedNul(what.to_string()));
    }
    // Upstream caps in-memory header strings at 65534 bytes.
    if s.len() > 65534 {
        return Err(Error::truncated(format!(
            "{what} exceeds the 65534-byte header string limit"
        )));
    }
    let len = s.len() as u32;
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(s.as_bytes());
    Ok(())
}

/// A parsed MCPL file: header plus the raw uncompressed bytes.
#[derive(Debug, Clone)]
pub struct McplFile {
    /// Source path when opened from disk, `None` for in-memory parses.
    pub path: Option<String>,
    /// Parsed header block.
    pub header: Header,
    data: Vec<u8>,
    header_len: usize,
}

impl McplFile {
    /// Open a file (`.gz` suffix reads through gzip transparently).
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let raw = std::fs::read(path).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
        let data = maybe_gunzip(path, raw)?;
        Self::from_bytes(data).map(|mut s| {
            s.path = Some(path.display().to_string());
            s
        })
    }

    /// Parse uncompressed MCPL bytes.
    pub fn from_bytes(data: Vec<u8>) -> Result<Self, Error> {
        let (header, header_len) = read_header(&data)?;
        Ok(McplFile {
            path: None,
            header,
            data,
            header_len,
        })
    }

    /// Raw uncompressed bytes of the whole file.
    pub fn raw_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Decode all particle records.
    pub fn particles(&self) -> Result<Vec<Particle>, Error> {
        let size = self.header.particle_size();
        let n = self.header.nparticles as usize;
        let want_end =
            self.header_len
                .checked_add(n.checked_mul(size).ok_or_else(|| {
                    Error::truncated("particle section larger than address space")
                })?)
                .ok_or_else(|| Error::truncated("particle section larger than address space"))?;
        if self.data.len() < want_end {
            return Err(Error::truncated(format!(
                "particle {} of {}",
                (self.data.len().saturating_sub(self.header_len)) / size.max(1),
                n
            )));
        }
        if self.data.len() > want_end {
            return Err(Error::TrailingBytes {
                expected: want_end,
                found: self.data.len(),
                particle: n as u64,
            });
        }
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let off = self.header_len + i * size;
            out.push(decode_particle(&self.header, &self.data[off..off + size])?);
        }
        Ok(out)
    }
}

fn maybe_gunzip(path: &Path, raw: Vec<u8>) -> Result<Vec<u8>, Error> {
    if path.extension().is_some_and(|e| e == "gz") {
        let mut dec = flate2::read::GzDecoder::new(&raw[..]);
        let mut out = Vec::new();
        dec.read_to_end(&mut out)
            .map_err(|e| Error::Io(format!("gzip decode {}: {e}", path.display())))?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

// ── Direction packing ─────────────────────────────────────────────
// Adaptive Projection Packing (format 3): pack the unit direction into two
// floats plus one sign bit; the energy rides in that sign slot. Branch choice
// is by largest-magnitude component; the largest component `z` is stored as
// `1/z` (which exceeds 1 in magnitude) so the reader can recover the branch.

/// Pack a unit direction with Adaptive Projection Packing (format 3).
fn pack_adaptproj(dir: [f64; 3]) -> [f64; 3] {
    let (ax, ay, az) = (dir[0].abs(), dir[1].abs(), dir[2].abs());
    if az < ax.max(ay) {
        let invz = if dir[2] != 0.0 {
            1.0 / dir[2]
        } else {
            f64::INFINITY
        };
        if ax >= ay {
            [invz, dir[1], 1.0f64.copysign(dir[0])]
        } else {
            [dir[0], invz, 1.0f64.copysign(dir[1])]
        }
    } else {
        [dir[0], dir[1], 1.0f64.copysign(dir[2])]
    }
}

/// Unpack an Adaptive-Projection triple (third slot carries only a sign).
fn unpack_adaptproj(packed: [f64; 3]) -> [f64; 3] {
    debug_assert!(packed[2] == 1.0 || packed[2] == -1.0);
    if packed[0].abs() > 1.0 {
        let z = 1.0 / packed[0];
        let x = packed[2] * (1.0 - (packed[1] * packed[1] + z * z)).max(0.0).sqrt();
        [x, packed[1], z]
    } else if packed[1].abs() > 1.0 {
        let z = 1.0 / packed[1];
        let y = packed[2] * (1.0 - (packed[0] * packed[0] + z * z)).max(0.0).sqrt();
        [packed[0], y, z]
    } else {
        let z = packed[2]
            * (1.0 - (packed[0] * packed[0] + packed[1] * packed[1]))
                .max(0.0)
                .sqrt();
        [packed[0], packed[1], z]
    }
}

/// Unpack an octahedral-packed pair (format 2, read-only).
fn unpack_oct(packed: [f64; 3]) -> [f64; 3] {
    let z = 1.0 - packed[0].abs() - packed[1].abs();
    // Upstream folds with `(v >= 0 ? +1 : -1)`; Rust `signum` maps
    // `0.0` to `0.0`, so spell the comparison out for exact parity.
    let sign = |v: f64| if v >= 0.0 { 1.0 } else { -1.0 };
    let (x, y) = if z < 0.0 {
        (
            (1.0 - packed[1].abs()) * sign(packed[0]),
            (1.0 - packed[0].abs()) * sign(packed[1]),
        )
    } else {
        (packed[0], packed[1])
    };
    let n = 1.0 / (x * x + y * y + z * z).sqrt();
    [x * n, y * n, z * n]
}

// ── Low-level codec ───────────────────────────────────────────────

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }

    fn take(&mut self, n: usize, what: &str) -> Result<&'a [u8], Error> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or_else(|| Error::truncated(what))?;
        if end > self.buf.len() {
            return Err(Error::truncated(what));
        }
        let s = &self.buf[self.pos..end];
        self.pos = end;
        Ok(s)
    }

    fn u32(&mut self, what: &str) -> Result<u32, Error> {
        Ok(u32::from_le_bytes(self.take(4, what)?.try_into().unwrap()))
    }

    fn u64(&mut self, what: &str) -> Result<u64, Error> {
        Ok(u64::from_le_bytes(self.take(8, what)?.try_into().unwrap()))
    }

    fn f32(&mut self, what: &str) -> Result<f64, Error> {
        Ok(f32::from_le_bytes(self.take(4, what)?.try_into().unwrap()) as f64)
    }

    fn f64(&mut self, what: &str) -> Result<f64, Error> {
        Ok(f64::from_le_bytes(self.take(8, what)?.try_into().unwrap()))
    }

    fn i32(&mut self, what: &str) -> Result<i32, Error> {
        Ok(i32::from_le_bytes(self.take(4, what)?.try_into().unwrap()))
    }

    fn string(&mut self, what: &str) -> Result<String, Error> {
        let n = self.u32(what)? as usize;
        let raw = self.take(n, what)?;
        if raw.contains(&0) {
            return Err(Error::EmbeddedNul(what.to_string()));
        }
        String::from_utf8(raw.to_vec()).map_err(|_| Error::InvalidUtf8(what.to_string()))
    }
}

fn read_header(data: &[u8]) -> Result<(Header, usize), Error> {
    let mut c = Cursor::new(data);
    let magic = c.take(4, "magic")?;
    if magic != MAGIC {
        return Err(Error::BadMagic);
    }
    let vbytes = c.take(3, "version")?;
    if !vbytes.iter().all(u8::is_ascii_digit) {
        return Err(Error::BadMagic);
    }
    let version = u16::from(vbytes[0] - b'0') * 100
        + u16::from(vbytes[1] - b'0') * 10
        + u16::from(vbytes[2] - b'0');
    if version != 2 && version != 3 {
        return Err(Error::UnsupportedVersion(version));
    }
    let endian = c.take(1, "endianness")?[0] as char;
    if endian == 'B' {
        return Err(Error::UnsupportedEndianness('B'));
    }
    if endian != 'L' {
        return Err(Error::BadMagic);
    }
    let nparticles = c.u64("nparticles")?;
    let ncomments = c.u32("ncomments")?;
    let nblobs = c.u32("nblobs")?;
    let has_userflags = c.u32("has_userflags")? != 0;
    let has_polarisation = c.u32("has_polarisation")? != 0;
    let single_prec = c.u32("single_prec")?;
    let universal_pdgcode_raw = c.u32("universal_pdgcode")?;
    let particle_size_stored = c.u32("particle_size")? as usize;
    let has_uweight = c.u32("has_universal_weight")? != 0;

    // Provisional header for the size cross-check.
    let provisional = Header {
        version,
        has_userflags,
        has_polarisation,
        double_prec: single_prec == 0,
        universal_pdgcode: {
            let v = universal_pdgcode_raw as i32;
            if v == 0 {
                None
            } else {
                Some(v)
            }
        },
        universal_weight: None,
        srcname: String::new(),
        comments: Vec::new(),
        blobs: Vec::new(),
        nparticles,
    };
    // `particle_size` is implied by the flags; a mismatch means the file
    // uses an option combination this reader does not understand. The
    // weight slot is present exactly when there is no universal weight,
    // so the provisional header must already reflect the flag.
    let implied = Header {
        universal_weight: has_uweight.then_some(0.0),
        ..provisional
    }
    .particle_size();
    if implied != particle_size_stored {
        return Err(Error::ParticleSizeMismatch {
            expected: implied,
            found: particle_size_stored,
        });
    }

    let universal_weight = if has_uweight {
        Some(c.f64("universal weight")?)
    } else {
        None
    };
    let srcname = c.string("srcname")?;
    let mut comments = Vec::with_capacity(ncomments.min(1 << 20) as usize);
    for i in 0..ncomments {
        comments.push(c.string(&format!("comment {i}"))?);
    }
    validate_statsum_comments(&comments)?;
    let mut keys = Vec::with_capacity(nblobs.min(1 << 20) as usize);
    for i in 0..nblobs {
        keys.push(c.string(&format!("blob key {i}"))?);
    }
    let mut blobs = Vec::with_capacity(nblobs.min(1 << 20) as usize);
    for (i, key) in keys.into_iter().enumerate() {
        let n = c.u32(&format!("blob {i} length"))? as usize;
        let raw = c.take(n, &format!("blob {i} data"))?.to_vec();
        blobs.push(Blob { key, data: raw });
    }
    let header = Header {
        version,
        has_userflags,
        has_polarisation,
        double_prec: single_prec == 0,
        universal_pdgcode: {
            let v = universal_pdgcode_raw as i32;
            if v == 0 {
                None
            } else {
                Some(v)
            }
        },
        universal_weight,
        srcname,
        comments,
        blobs,
        nparticles,
    };
    Ok((header, c.pos))
}

/// Validate every `stat:sum:` comment (upstream syntax, enforced on open
/// and on write) and reject duplicated keys (upstream errors at open time).
fn validate_statsum_comments(comments: &[String]) -> Result<(), Error> {
    let mut seen: Vec<&str> = Vec::new();
    for (i, cmt) in comments.iter().enumerate() {
        if !cmt.starts_with("stat:sum:") {
            continue;
        }
        let key = statsum_validate(cmt).map_err(|reason| Error::BadStatSum { index: i, reason })?;
        if seen.contains(&key) {
            return Err(Error::DuplicateStatSum(i, key.to_string()));
        }
        seen.push(key);
    }
    Ok(())
}

/// Check one `stat:sum:<key>:<24-char value>` comment against the upstream
/// syntax (key `[a-zA-Z][a-zA-Z0-9_]*` of 1–64 chars; value field exactly
/// 24 chars holding a finite `-1` or `>= 0` double over `[0-9.\-+eE ]`).
/// Returns the key. Anything else starting with `stat:` (but not
/// `stat:sum:`) is only reserved prose upstream and passes through here.
pub fn statsum_validate(comment: &str) -> Result<&str, String> {
    let rest = comment
        .strip_prefix("stat:sum:")
        .ok_or_else(|| "missing stat:sum: prefix".to_string())?;
    let (key, field) = rest
        .split_once(':')
        .ok_or_else(|| "did not find colon separating key and value".to_string())?;
    if key.is_empty() {
        return Err("empty key".to_string());
    }
    if key.len() > 64 {
        return Err("key length exceeds 64 characters".to_string());
    }
    let mut bytes = key.bytes();
    match bytes.next() {
        Some(b) if b.is_ascii_alphabetic() => {}
        _ => {
            return Err("key does not adhere to naming [a-zA-Z][a-zA-Z0-9_]*".to_string());
        }
    }
    if !bytes.all(|b| b.is_ascii_alphanumeric() || b == b'_') {
        return Err("key does not adhere to naming [a-zA-Z][a-zA-Z0-9_]*".to_string());
    }
    if field.len() != 24 {
        return Err("value field is not exactly 24 characters wide".to_string());
    }
    let trimmed = field.trim_matches(' ');
    if trimmed.is_empty() {
        return Err("value field missing actual value".to_string());
    }
    if !trimmed
        .bytes()
        .all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'-' | b'+' | b'e' | b'E'))
    {
        return Err("value field holds forbidden characters".to_string());
    }
    let value: f64 = trimmed
        .parse()
        .map_err(|_| "could not decode contents of value field".to_string())?;
    if !value.is_finite() {
        return Err("value field holds forbidden value (NaN or infinity)".to_string());
    }
    if !(value >= 0.0 || value == -1.0) {
        return Err("value field must hold non-negative value or -1".to_string());
    }
    Ok(key)
}

/// Build a well-formed `stat:sum:` comment for `key` and `value`
/// (key `[a-zA-Z][a-zA-Z0-9_]*`, value finite and `>= 0` or `-1`).
/// The 24-char value field is emitted losslessly: plain notation when it
/// fits, scientific notation otherwise, and an error when neither holds the
/// value exactly. The result always passes [`statsum_validate`].
pub fn statsum_comment(key: &str, value: f64) -> Result<String, Error> {
    if key.is_empty() || key.len() > 64 {
        return Err(Error::BadStatSum {
            index: 0,
            reason: "empty key or key length exceeds 64 characters".to_string(),
        });
    }
    let bad_key = || Error::BadStatSum {
        index: 0,
        reason: "key does not adhere to naming [a-zA-Z][a-zA-Z0-9_]*".to_string(),
    };
    let mut bytes = key.bytes();
    match bytes.next() {
        Some(b) if b.is_ascii_alphabetic() => {}
        _ => return Err(bad_key()),
    }
    if !bytes.all(|b| b.is_ascii_alphanumeric() || b == b'_') {
        return Err(bad_key());
    }
    if !(value.is_finite() && (value >= 0.0 || value == -1.0)) {
        return Err(Error::BadStatSum {
            index: 0,
            reason: "value must be finite and non-negative or -1".to_string(),
        });
    }
    // Upstream special-cases zero so no negative-zero sign leaks in.
    if value == 0.0 {
        return Ok(format!("stat:sum:{key}:{}", " ".repeat(23) + "0"));
    }
    for field in [
        format!("{value:24}"),
        format!("{value:24.17e}"),
        format!("{value:24.15e}"),
    ] {
        if field.len() != 24 {
            continue;
        }
        let back: f64 = match field.trim_matches(' ').parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if back == value {
            let comment = format!("stat:sum:{key}:{field}");
            statsum_validate(&comment).map_err(|reason| Error::BadStatSum { index: 0, reason })?;
            return Ok(comment);
        }
    }
    Err(Error::BadStatSum {
        index: 0,
        reason: "value has no lossless 24-char encoding".to_string(),
    })
}

fn decode_particle(h: &Header, rec: &[u8]) -> Result<Particle, Error> {
    let mut c = Cursor::new(rec);
    let single = !h.double_prec;
    let mut get_fp = |what: &str| -> Result<f64, Error> {
        if single {
            c.f32(what)
        } else {
            c.f64(what)
        }
    };
    let mut polarisation = [0.0; 3];
    if h.has_polarisation {
        for v in &mut polarisation {
            *v = get_fp("polarisation")?;
        }
    }
    let mut position = [0.0; 3];
    for v in &mut position {
        *v = get_fp("position")?;
    }
    let mut pack = [0.0; 3];
    for v in &mut pack {
        *v = get_fp("packed direction/energy")?;
    }
    let time = get_fp("time")?;
    let weight = if let Some(w) = h.universal_weight {
        w
    } else {
        get_fp("weight")?
    };
    let pdgcode = if let Some(p) = h.universal_pdgcode {
        p
    } else {
        c.i32("pdgcode")?
    };
    let userflags = if h.has_userflags {
        c.u32("userflags")?
    } else {
        0
    };

    let (ekin, direction) = if h.version >= 3 {
        let ekin = pack[2].abs();
        let sign = if pack[2].is_sign_negative() {
            -1.0
        } else {
            1.0
        };
        (ekin, unpack_adaptproj([pack[0], pack[1], sign]))
    } else {
        let mut direction = unpack_oct(pack);
        let mut ekin = pack[2];
        if pack[2].is_sign_negative() {
            ekin = -ekin;
            direction[2] = 0.0;
        }
        (ekin, direction)
    };
    Ok(Particle {
        ekin,
        polarisation,
        position,
        direction,
        time,
        weight,
        pdgcode,
        userflags,
    })
}

fn encode_particle(h: &Header, p: &Particle, index: usize, out: &mut Vec<u8>) -> Result<(), Error> {
    let dir2 = p.direction[0] * p.direction[0]
        + p.direction[1] * p.direction[1]
        + p.direction[2] * p.direction[2];
    if (dir2 - 1.0).abs() > UNIT_TOL {
        return Err(Error::NonUnitDirection(index, dir2));
    }
    if p.ekin < 0.0 {
        return Err(Error::NegativeEnergy(index, p.ekin));
    }
    let mut pack = pack_adaptproj(p.direction);
    pack[2] = p.ekin.copysign(pack[2]);

    if h.double_prec {
        if h.has_polarisation {
            for v in p.polarisation {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        for v in p.position {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for v in pack {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&p.time.to_le_bytes());
        if h.universal_weight.is_none() {
            out.extend_from_slice(&p.weight.to_le_bytes());
        }
    } else {
        if h.has_polarisation {
            for v in p.polarisation {
                out.extend_from_slice(&(v as f32).to_le_bytes());
            }
        }
        for v in p.position {
            out.extend_from_slice(&(v as f32).to_le_bytes());
        }
        for v in pack {
            out.extend_from_slice(&(v as f32).to_le_bytes());
        }
        out.extend_from_slice(&(p.time as f32).to_le_bytes());
        if h.universal_weight.is_none() {
            out.extend_from_slice(&(p.weight as f32).to_le_bytes());
        }
    }
    if h.universal_pdgcode.is_none() {
        out.extend_from_slice(&p.pdgcode.to_le_bytes());
    }
    if h.has_userflags {
        out.extend_from_slice(&p.userflags.to_le_bytes());
    }
    Ok(())
}

/// Encode a complete file: header block plus particle records.
///
/// The stored `nparticles` is the length of `particles`; any count already
/// in `header.nparticles` is ignored so callers can build the header before
/// knowing the tally.
pub fn encode_file(header: &Header, particles: &[Particle]) -> Result<Vec<u8>, Error> {
    let mut h = header.clone();
    h.nparticles = particles.len() as u64;
    let mut out = h.header_block()?;
    out.reserve(particles.len() * h.particle_size());
    for (i, p) in particles.iter().enumerate() {
        encode_particle(&h, p, i, &mut out)?;
    }
    Ok(out)
}

/// Write a complete MCPL file (header block plus particles) to a stream.
pub fn write_to<W: Write>(w: &mut W, header: &Header, particles: &[Particle]) -> Result<(), Error> {
    let bytes = encode_file(header, particles)?;
    w.write_all(&bytes).map_err(|e| Error::Io(e.to_string()))
}

/// Convenience wrapper around [`write_to`] that writes to a file path.
///
/// A `.gz` suffix compresses through gzip transparently (compressed bytes
/// are encoder-dependent and never asserted); other paths write raw MCPL.
pub fn write_to_path<P: AsRef<Path>>(
    path: P,
    header: &Header,
    particles: &[Particle],
) -> Result<(), Error> {
    let bytes = encode_file(header, particles)?;
    if path.as_ref().extension().is_some_and(|e| e == "gz") {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(&bytes)
            .map_err(|e| Error::Io(e.to_string()))?;
        let compressed = enc.finish().map_err(|e| Error::Io(e.to_string()))?;
        std::fs::write(path, compressed).map_err(|e| Error::Io(e.to_string()))
    } else {
        std::fs::write(path, bytes).map_err(|e| Error::Io(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axis_particles() -> Vec<Particle> {
        vec![
            Particle {
                ekin: 2.5,
                polarisation: [0.0; 3],
                position: [1.0, -2.0, 0.5],
                direction: [0.0, 0.0, 1.0],
                time: 3.0,
                weight: 1.0,
                pdgcode: 2112,
                userflags: 0,
            },
            Particle {
                ekin: 0.662,
                polarisation: [0.0; 3],
                position: [0.0, 0.0, 0.0],
                direction: [1.0, 0.0, 0.0],
                time: 0.0,
                weight: 0.5,
                pdgcode: 22,
                userflags: 7,
            },
        ]
    }

    fn assert_round_trip_bytes(header: &Header, particles: &[Particle]) {
        let bytes = encode_file(header, particles).unwrap();
        let file = McplFile::from_bytes(bytes.clone()).unwrap();
        assert_eq!(file.header.version, FORMAT_VERSION_WRITE);
        assert_eq!(file.header.nparticles, particles.len() as u64);
        // Header identity (ignoring the patched count, already checked).
        let mut want_h = header.clone();
        want_h.nparticles = particles.len() as u64;
        let mut got_h = file.header.clone();
        assert_eq!(got_h.srcname, want_h.srcname);
        assert_eq!(got_h.comments, want_h.comments);
        assert_eq!(got_h.blobs, want_h.blobs);
        got_h.nparticles = want_h.nparticles;
        assert_eq!(got_h, want_h);
        // Byte-exact re-emit from the decoded particles.
        let got_p = file.particles().unwrap();
        let rewritten = encode_file(&file.header, &got_p).unwrap();
        assert_eq!(rewritten.len(), bytes.len(), "length mismatch");
        assert_eq!(rewritten, bytes, "byte mismatch");
    }

    #[test]
    fn round_trip_single_basic() {
        let h = Header::default();
        assert_round_trip_bytes(&h, &axis_particles());
    }

    #[test]
    fn round_trip_double_polar_user() {
        let h = Header {
            double_prec: true,
            has_polarisation: true,
            has_userflags: true,
            srcname: "nucleide-test".to_string(),
            comments: vec!["synthetic double-prec probe".to_string()],
            blobs: vec![Blob {
                key: "probe".to_string(),
                data: vec![1, 2, 3, 4],
            }],
            ..Header::default()
        };
        let mut ps = axis_particles();
        ps[0].polarisation = [0.1, 0.2, 0.3];
        ps[1].polarisation = [-0.5, 0.0, 0.5];
        ps[1].userflags = 0xdead_beef;
        assert_round_trip_bytes(&h, &ps);
    }

    #[test]
    fn round_trip_universal() {
        let h = Header {
            universal_pdgcode: Some(2112),
            universal_weight: Some(1.5),
            srcname: "universal".to_string(),
            ..Header::default()
        };
        // Per-particle pdg/weight ride along in memory but are not stored.
        assert_round_trip_bytes(&h, &axis_particles());
        let bytes = encode_file(&h, &axis_particles()).unwrap();
        let file = McplFile::from_bytes(bytes).unwrap();
        let ps = file.particles().unwrap();
        assert!(ps.iter().all(|p| p.pdgcode == 2112));
        assert!(ps.iter().all(|p| p.weight == 1.5));
    }

    #[test]
    fn header_block_prefix_matches_file() {
        let h = Header {
            srcname: "prefix".to_string(),
            comments: vec!["a".to_string(), statsum_comment("nps", 1.0).unwrap()],
            ..Header::default()
        };
        let mut hh = h.clone();
        hh.nparticles = 2;
        let block = hh.header_block().unwrap();
        let bytes = encode_file(&h, &axis_particles()).unwrap();
        assert_eq!(&bytes[..block.len()], &block[..]);
    }

    #[test]
    fn adaptproj_round_trip_precise() {
        // Axis + z-dominant vectors stay in-branch through f32 quantization.
        for dir in [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.1, 0.2, 0.9746794344808963],
        ] {
            let packed = pack_adaptproj(dir);
            let sign = if packed[2].is_sign_negative() {
                -1.0
            } else {
                1.0
            };
            let unpacked = unpack_adaptproj([packed[0], packed[1], sign]);
            for (a, b) in dir.iter().zip(unpacked.iter()) {
                assert!((a - b).abs() < 1e-12, "{dir:?} -> {unpacked:?}");
            }
            // Through single precision the direction survives to ~1e-7.
            let q = [packed[0] as f32 as f64, packed[1] as f32 as f64, sign];
            let uq = unpack_adaptproj(q);
            for (a, b) in dir.iter().zip(uq.iter()) {
                assert!((a - b).abs() < 1e-6, "{dir:?} -> {uq:?}");
            }
        }
    }

    #[test]
    fn oct_unpack_axis_vectors() {
        for (packed, want) in [
            ([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
            ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
        ] {
            let got = unpack_oct(packed);
            for (a, b) in got.iter().zip(want.iter()) {
                assert!((a - b).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn v2_file_reads_oct_direction() {
        // Hand-framed format-2 file: octahedral (0,0) unpacks to +z.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"MCPL002L");
        bytes.extend_from_slice(&1u64.to_le_bytes());
        let opts = [0u32, 0, 0, 0, 1, 0, 7 * 4 + 4 + 4, 0];
        for o in opts {
            bytes.extend_from_slice(&o.to_le_bytes());
        }
        put_string(&mut bytes, "srcname", "v2probe").unwrap();
        // position(0,0,0) + oct-pack(0,0) with ekin 1.25 in slot 2 + time 0 + weight 1.
        for v in [0.0f32, 0.0, 0.0, 0.0, 0.0, 1.25, 0.0, 1.0] {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        bytes.extend_from_slice(&2112i32.to_le_bytes());
        let file = McplFile::from_bytes(bytes).unwrap();
        assert_eq!(file.header.version, 2);
        let ps = file.particles().unwrap();
        assert_eq!(ps.len(), 1);
        assert!((ps[0].ekin - 1.25).abs() < 1e-9);
        assert!((ps[0].direction[2] - 1.0).abs() < 1e-9);
        assert_eq!(ps[0].pdgcode, 2112);
    }

    #[test]
    fn gzip_round_trip_semantic() {
        let h = Header::default();
        let dir = std::env::temp_dir().join("nucleide_mcpl_rt.mcpl.gz");
        write_to_path(&dir, &h, &axis_particles()).unwrap();
        let file = McplFile::open(&dir).unwrap();
        let ps = file.particles().unwrap();
        assert_eq!(ps.len(), 2);
        assert!((ps[0].ekin - 2.5).abs() < 1e-6);
        assert!((ps[1].direction[0] - 1.0).abs() < 1e-6);
        std::fs::remove_file(&dir).ok();
    }

    #[test]
    fn rejects_bad_magic() {
        let mut bytes = encode_file(&Header::default(), &axis_particles()).unwrap();
        bytes[0] = b'X';
        assert_eq!(McplFile::from_bytes(bytes).unwrap_err(), Error::BadMagic);
    }

    #[test]
    fn rejects_big_endian() {
        let mut bytes = encode_file(&Header::default(), &axis_particles()).unwrap();
        bytes[7] = b'B';
        assert_eq!(
            McplFile::from_bytes(bytes).unwrap_err(),
            Error::UnsupportedEndianness('B')
        );
    }

    #[test]
    fn rejects_unknown_version() {
        let mut bytes = encode_file(&Header::default(), &axis_particles()).unwrap();
        bytes[4] = b'0';
        bytes[5] = b'0';
        bytes[6] = b'9';
        assert_eq!(
            McplFile::from_bytes(bytes).unwrap_err(),
            Error::UnsupportedVersion(9)
        );
    }

    #[test]
    fn rejects_truncation() {
        let bytes = encode_file(&Header::default(), &axis_particles()).unwrap();
        let cut = McplFile::from_bytes(bytes[..bytes.len() / 2].to_vec());
        match cut {
            // A header-only prefix still parses; the missing records fail.
            Ok(file) => assert!(matches!(file.particles().unwrap_err(), Error::Truncated(_))),
            Err(Error::Truncated(_)) => {}
            Err(other) => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn writer_rejects_non_unit_direction() {
        let mut ps = axis_particles();
        ps[0].direction = [1.0, 1.0, 1.0];
        let err = encode_file(&Header::default(), &ps).unwrap_err();
        assert!(matches!(err, Error::NonUnitDirection(0, _)));
    }

    #[test]
    fn writer_rejects_negative_energy() {
        let mut ps = axis_particles();
        ps[0].ekin = -1.0;
        let err = encode_file(&Header::default(), &ps).unwrap_err();
        assert_eq!(err, Error::NegativeEnergy(0, -1.0));
    }

    #[test]
    fn writer_rejects_bad_universal_weight() {
        let h = Header {
            universal_weight: Some(0.0),
            ..Header::default()
        };
        assert_eq!(
            encode_file(&h, &axis_particles()).unwrap_err(),
            Error::BadUniversalWeight(0.0)
        );
    }

    #[test]
    fn rejects_duplicate_statsum() {
        // Write-side rejection first (upstream refuses these on write too).
        let h = Header {
            comments: vec![
                statsum_comment("nps", 1.0).unwrap(),
                statsum_comment("nps", 2.0).unwrap(),
            ],
            ..Header::default()
        };
        assert!(matches!(
            encode_file(&h, &[]).unwrap_err(),
            Error::DuplicateStatSum(1, _)
        ));
        // Read-side rejection: splice a duplicate key into otherwise valid
        // bytes (equal-length keys keep every offset intact).
        let h = Header {
            comments: vec![
                statsum_comment("nps", 1.0).unwrap(),
                statsum_comment("aaa", 2.0).unwrap(),
            ],
            ..Header::default()
        };
        let mut bytes = encode_file(&h, &[]).unwrap();
        let needle = b"stat:sum:aaa:";
        let pos = bytes
            .windows(needle.len())
            .position(|w| w == needle)
            .expect("spliced key must occur");
        bytes[pos + 9..pos + 12].copy_from_slice(b"nps");
        assert!(matches!(
            McplFile::from_bytes(bytes).unwrap_err(),
            Error::DuplicateStatSum(1, _)
        ));
    }

    #[test]
    fn statsum_comment_round_trip_values() {
        // Hand-computed encodings: plain notation when it fits, scientific
        // otherwise; every emission revalidates under the upstream syntax.
        for (key, value, want_field) in [
            ("nps", 1.0, " ".repeat(23) + "1"),
            ("nps", 0.0, " ".repeat(23) + "0"),
            ("nps", -1.0, " ".repeat(22) + "-1"),
        ] {
            let comment = statsum_comment(key, value).unwrap();
            let (head, field) = comment.rsplit_once(':').unwrap();
            assert_eq!(head, "stat:sum:nps");
            assert_eq!(field, want_field);
            assert_eq!(field.len(), 24);
            assert_eq!(statsum_validate(&comment).unwrap(), "nps");
        }
        let big = statsum_comment("nps", 1.0e300).unwrap();
        assert!(big.ends_with("e300"));
        assert_eq!(statsum_validate(&big).unwrap(), "nps");
    }

    #[test]
    fn rejects_malformed_statsum() {
        // Legacy `=` separator, digit-leading key, short value field,
        // forbidden value character, and over-long key are upstream-invalid.
        for bad in [
            "stat:sum:nps=1".to_string(),
            "stat:sum:1nps:                       1".to_string(),
            "stat:sum:nps:1".to_string(),
            format!("stat:sum:nps:{}1x", " ".repeat(22)),
            format!("stat:sum:{}:{}", "k".repeat(65), " ".repeat(23) + "1"),
        ] {
            assert!(
                statsum_validate(&bad).is_err(),
                "accepted malformed comment: {bad:?}"
            );
        }
        assert!(statsum_comment("1nps", 1.0).is_err());
        assert!(statsum_comment("nps", -2.0).is_err());
        assert!(statsum_comment("nps", f64::NAN).is_err());
        assert!(statsum_comment("nps", f64::INFINITY).is_err());
        // Malformed comments are rejected on write, matching upstream.
        let h = Header {
            comments: vec!["stat:sum:nps=1".to_string()],
            ..Header::default()
        };
        assert!(matches!(
            h.header_block().unwrap_err(),
            Error::BadStatSum { .. }
        ));
        // `stat:`-prefixed non-sum prose stays verbatim (upstream warns only).
        let h = Header {
            comments: vec!["stat:future-thing".to_string()],
            ..Header::default()
        };
        let bytes = encode_file(&h, &[]).unwrap();
        let file = McplFile::from_bytes(bytes).unwrap();
        assert_eq!(file.header.comments, vec!["stat:future-thing".to_string()]);
    }

    #[test]
    fn rejects_duplicate_blob_keys() {
        let h = Header {
            blobs: vec![
                Blob {
                    key: "k".to_string(),
                    data: vec![1],
                },
                Blob {
                    key: "k".to_string(),
                    data: vec![2],
                },
            ],
            ..Header::default()
        };
        assert!(matches!(
            encode_file(&h, &[]).unwrap_err(),
            Error::DuplicateBlobKey(_)
        ));
    }
}
