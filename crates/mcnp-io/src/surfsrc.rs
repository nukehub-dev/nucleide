//! MCNP SSW surface-source file reading (`SurfSrc`, a.k.a. ssr/ssw),
//! validated against the vendored fixtures (`mcnp5_surfsrc.w`,
//! `mcnp6_surfsrc.w`, `mcnpx_surfsrc.w`, `mcnp_surfsrc_onetrack.w`).
//!
//! Format notes:
//! - Little-endian Fortran unformatted records `[i32 len][payload][i32 len]`.
//! - Header branching follows the code identifier: `SF_00001` (MCNP6) splits
//!   the header into two records; otherwise one record holds kod/ver/loddat.
//!   Version quirks preserved verbatim: MCNP-2.6.0 stores `np1`/`nrss` as
//!   i32 while plain MCNP5 stores them as i64.
//! - A negative original `np1` signals an extra table-2 record
//!   (cells/particle/macrobody-facet info).
//! - Deviations: `print_header` fixes the mis-templated counts line of the legacy writer;
//!   absent table-2 fields compare as unequal rather than crashing.

use std::cmp::Ordering;
use std::fmt;
use std::io::{Read, Write};
use std::path::Path;

use crate::fortran::{self, read_record, RecordSink};

/// Errors raised while reading or writing SSW data.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    /// Leading/trailing record-length markers disagreed.
    BadRecordMarker {
        lead: i32,
        trailer: i32,
    },
    /// Payload exhausted mid-field.
    ShortRecord {
        need: usize,
        left: usize,
    },
    /// Unsupported code/version signature.
    UnsupportedVersion(String),
    /// Writer: negative `orignp1` promises a table-2 record but its fields
    /// were never parsed (or were stripped) from the header.
    MissingTable2,
    /// Writer: tracklist length disagrees with the header's `nrss`.
    TrackCountMismatch {
        expected: u64,
        found: usize,
    },
    /// Writer: a track's record width disagrees with `abs(ncrd)`.
    TrackRecordWidth {
        index: usize,
        expected: usize,
        found: usize,
    },
    /// `combine_files`: headers disagree on a compared field, carry
    /// unsupported `niwr` extras, or promise an unrepresentable layout.
    Incompatible(String),
}

impl From<fortran::Error> for Error {
    fn from(e: fortran::Error) -> Self {
        match e {
            fortran::Error::Io(m) => Error::Io(m),
            fortran::Error::BadRecordMarker { lead, trailer } => {
                Error::BadRecordMarker { lead, trailer }
            }
            fortran::Error::ShortRecord { need, left } => Error::ShortRecord { need, left },
        }
    }
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
            Error::UnsupportedVersion(v) => write!(f, "MCNP version `{v}` not supported"),
            Error::MissingTable2 => write!(
                f,
                "negative orignp1 requires table-2 fields (niwr/mipts/kjaq)"
            ),
            Error::TrackCountMismatch { expected, found } => {
                write!(f, "tracklist holds {found} tracks, header says {expected}")
            }
            Error::TrackRecordWidth {
                index,
                expected,
                found,
            } => write!(
                f,
                "track {index} holds {found} values, ncrd says {expected}"
            ),
            Error::Incompatible(m) => write!(f, "surface-source files incompatible: {m}"),
        }
    }
}

impl std::error::Error for Error {}

/// One surface entry from the header's per-surface records.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceSurf {
    pub id: i32,
    pub facet_id: i32,
    pub surface_type: i32,
    pub num_params: usize,
    pub surf_params: Vec<f64>,
}

/// Parsed SSW header block.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfSrcHeader {
    pub kod: String,
    pub ver: String,
    pub loddat: String,
    pub idtm: String,
    pub probid: String,
    pub aid: String,
    pub knod: i32,
    /// Histories used to generate the source (absolute value).
    pub np1: i64,
    /// Signed `np1` exactly as stored (negative ⇒ table-2 present).
    pub orignp1: i64,
    pub nrss: i64,
    pub ncrd: i32,
    pub njsw: i32,
    pub niss: i64,
    /// Present only when the file carries table 2.
    pub niwr: Option<i32>,
    pub mipts: Option<i32>,
    pub kjaq: Option<i32>,
    pub table1extra: Vec<i32>,
    pub table2extra: Vec<i32>,
    pub surflist: Vec<SourceSurf>,
    pub summary_table: Vec<i32>,
    pub summary_extra: Vec<i32>,
}

impl SurfSrcHeader {
    /// Comparison semantics: header identity only;
    /// tracklists are never compared.
    pub fn cmp_semantics(&self, other: &Self) -> Ordering {
        self.kod
            .cmp(&other.kod)
            .then_with(|| self.ver.cmp(&other.ver))
            .then_with(|| self.loddat.cmp(&other.loddat))
            .then_with(|| self.ncrd.cmp(&other.ncrd))
            .then_with(|| self.njsw.cmp(&other.njsw))
            .then_with(|| self.np1.cmp(&other.np1))
            .then_with(|| self.nrss.cmp(&other.nrss))
            .then_with(|| self.niss.cmp(&other.niss))
            .then_with(|| self.niwr.cmp(&other.niwr))
            .then_with(|| self.mipts.cmp(&other.mipts))
            .then_with(|| self.kjaq.cmp(&other.kjaq))
            .then_with(|| cmp_surflist(&self.surflist, &other.surflist))
    }

    /// Informative header rendering (counts line corrected from the
    /// mis-templated upstream version).
    pub fn print_header(&self) -> String {
        let mut s = format!(
            "Code: {} (version: {}) [{}]\n",
            self.kod, self.ver, self.loddat
        );
        s += &format!(
            "Problem info: ({}) {}\n{}\n",
            self.idtm, self.probid, self.aid
        );
        s += &format!("Showing dump #{}\n", self.knod);
        s += &format!(
            "{} histories, {} tracks, {} record size, {} surfaces, {} histories\n",
            self.np1, self.nrss, self.ncrd, self.njsw, self.niss
        );
        s += &format!(
            "{} cells, source particle: {}, macrobody facet flag: {}\n",
            self.niwr
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".into()),
            self.mipts
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".into()),
            self.kjaq
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".into()),
        );
        for i in &self.surflist {
            s += &format!(
                "Surface {}: facet {}, type {} with {} parameters: ",
                i.id, i.facet_id, i.surface_type, i.num_params
            );
            s += "(";
            for p in &i.surf_params {
                s += &format!(" {p}");
            }
            s += ")\n";
        }
        s += &format!("Summary Table: {:?}", self.summary_table);
        s
    }
}

/// One track record from the tracklist.
#[derive(Debug, Clone, PartialEq)]
pub struct TrackData {
    pub record: Vec<f64>,
    pub nps: f64,
    pub bitarray: f64,
    pub wgt: f64,
    pub erg: f64,
    pub tme: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub u: f64,
    pub v: f64,
    pub cs: f64,
    pub w: f64,
}

impl TrackData {
    /// Cell decode: `abs(bitarray) // 8 % 100000000`.
    pub fn cell(&self) -> f64 {
        ((self.bitarray.abs() as u64) / 8 % 100_000_000) as f64
    }
}

/// A parsed MCNP surface-source file.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfSrc {
    pub path: Option<String>,
    pub header: SurfSrcHeader,
    data: Vec<u8>,
}

impl SurfSrc {
    /// Open a file and eagerly read its header block.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let data =
            std::fs::read(path).map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Self::from_bytes(data).map(|mut s| {
            s.path = Some(path.display().to_string());
            s
        })
    }

    /// Parse a surface-source file from its raw bytes.
    pub fn from_bytes(data: Vec<u8>) -> Result<Self, Error> {
        let mut cursor = std::io::Cursor::new(&data);
        let header = read_header(&mut cursor)?;
        Ok(SurfSrc {
            path: None,
            header,
            data,
        })
    }

    /// Read all track records following the header.
    pub fn read_tracklist(&self) -> Result<Vec<TrackData>, Error> {
        let mut cursor = std::io::Cursor::new(&self.data);
        skip_header_records(&mut cursor, &self.header)?;
        let ncrd_abs = self.header.ncrd.unsigned_abs() as usize;
        let mut tracks = Vec::with_capacity(self.header.nrss.max(0) as usize);
        for _ in 0..self.header.nrss {
            let mut rec = read_record(&mut cursor)?;
            let record = rec.get_f64_n(ncrd_abs)?;
            tracks.push(build_track(record));
        }
        Ok(tracks)
    }
}

impl TrackData {
    /// Number of doubles in a track record (abs(ncrd) at read time).
    pub const RECORD_WIDTH: usize = 11;

    /// Build a track from its raw record fields.
    pub fn from_record(record: Vec<f64>) -> Self {
        build_track(record)
    }
}

fn build_track(record: Vec<f64>) -> TrackData {
    let g = |i: usize| record.get(i).copied().unwrap_or_default();
    let nps = g(0);
    let bitarray = g(1);
    let u = g(8);
    let v = g(9);
    let cs = g(10);
    let wgt = g(2);
    let erg = g(3);
    let tme = g(4);
    let x = g(5);
    let y = g(6);
    let z = g(7);
    let w = (1.0 - u * u - v * v).max(0.0).sqrt().copysign(bitarray);
    TrackData {
        record,
        nps,
        bitarray,
        wgt,
        erg,
        tme,
        x,
        y,
        z,
        u,
        v,
        cs,
        w,
    }
}

/// Re-read (and discard) the header records so a fresh file cursor lands at
/// the start of the tracklist.
fn skip_header_records<R: Read>(f: &mut R, h: &SurfSrcHeader) -> Result<(), Error> {
    let _ = read_record(f)?; // header (or first of two for SF_00001)
    if h.kod.contains("SF_00001") {
        let _ = read_record(f)?; // second header record
    }
    let _ = read_record(f)?; // table 1
    if h.orignp1 < 0 {
        let _ = read_record(f)?; // table 2
    }
    for _ in 0..h.njsw.max(0) {
        let _ = read_record(f)?; // surface records
    }
    if let Some(n) = h.niwr {
        for _ in 0..n.max(0) {
            let _ = read_record(f)?; // unhandled extras
        }
    }
    let _ = read_record(f)?; // summary table
    Ok(())
}

/// Field-wise surface-list comparison using IEEE total ordering for params
/// (compares id/facet/type/num_params/params in that order).
fn cmp_surflist(a: &[SourceSurf], b: &[SourceSurf]) -> Ordering {
    for (sa, sb) in a.iter().zip(b.iter()) {
        let o = sa
            .id
            .cmp(&sb.id)
            .then_with(|| sa.facet_id.cmp(&sb.facet_id))
            .then_with(|| sa.surface_type.cmp(&sb.surface_type))
            .then_with(|| sa.num_params.cmp(&sb.num_params))
            .then_with(|| {
                sa.surf_params
                    .iter()
                    .zip(sb.surf_params.iter())
                    .map(|(x, y)| x.total_cmp(y))
                    .find(|o| *o != Ordering::Equal)
                    .unwrap_or(Ordering::Equal)
            });
        if o != Ordering::Equal {
            return o;
        }
    }
    a.len().cmp(&b.len())
}

/// Parse the full header block (4+ Fortran records) from a stream.
fn read_header<R: Read>(f: &mut R) -> Result<SurfSrcHeader, Error> {
    let mut header_rec = read_record(f)?;
    let kod = header_rec.get_string(8);

    // Strings common to every layout.
    let (ver, loddat): (String, String);
    let idtm: String;
    let probid: String;
    let aid: String;
    let knod: i32;

    // Table-1 counters.
    let np1_raw: i64;
    let nrss: i64;
    let ncrd: i32;
    let njsw: i32;
    let niss: i64;
    let table1extra: Vec<i32>;

    if kod.contains("SF_00001") {
        // MCNP6: header split over two records.
        let mut h2 = read_record(f)?;
        ver = h2.get_string(12);
        loddat = h2.get_string(9);
        idtm = h2.get_string(19);
        probid = h2.get_string(19);
        aid = h2.get_string(80);
        knod = h2.get_i32()?;

        let mut t1 = read_record(f)?;
        np1_raw = t1.get_i32()? as i64;
        let _notsure0 = t1.get_i32()?;
        nrss = t1.get_i32()? as i64;
        let _notsure1 = t1.get_i32()?;
        ncrd = t1.get_i32()?;
        njsw = t1.get_i32()?;
        niss = t1.get_i32()? as i64;
        table1extra = t1.drain_i32_extras()?;
    } else {
        ver = header_rec.get_string(5);
        if ver.contains("2.6.0") {
            loddat = header_rec.get_string(28);
            idtm = header_rec.get_string(19);
            probid = header_rec.get_string(19);
            aid = header_rec.get_string(80);
            knod = header_rec.get_i32()?;

            let mut t1 = read_record(f)?;
            np1_raw = t1.get_i32()? as i64;
            nrss = t1.get_i32()? as i64;
            ncrd = t1.get_i32()?;
            njsw = t1.get_i32()?;
            niss = t1.get_i32()? as i64;
            table1extra = t1.drain_i32_extras()?;
        } else if ver.contains('5') {
            loddat = header_rec.get_string(8);
            idtm = header_rec.get_string(19);
            probid = header_rec.get_string(19);
            aid = header_rec.get_string(80);
            knod = header_rec.get_i32()?;

            let mut t1 = read_record(f)?;
            // Plain MCNP5 stores the two big counters as i64.
            np1_raw = t1.get_i64()?;
            nrss = t1.get_i64()?;
            ncrd = t1.get_i32()?;
            njsw = t1.get_i32()?;
            niss = t1.get_i32()? as i64;
            table1extra = t1.drain_i32_extras()?;
        } else {
            return Err(Error::UnsupportedVersion(ver.trim().to_string()));
        }
    }

    let mut niwr = None;
    let mut mipts = None;
    let mut kjaq = None;
    let mut table2extra: Vec<i32> = Vec::new();

    if np1_raw < 0 {
        let mut t2 = read_record(f)?;
        niwr = Some(t2.get_i32()?);
        mipts = Some(t2.get_i32()?);
        kjaq = Some(t2.get_i32()?);
        table2extra = t2.drain_i32_extras()?;
    }

    // Per-surface records.
    let kjaq_flag = kjaq.unwrap_or(0);
    let mut surflist = Vec::with_capacity(njsw.max(0) as usize);
    for _ in 0..njsw {
        let mut rec = read_record(f)?;
        let id = rec.get_i32()?;
        let facet_id = if kjaq_flag == 1 { rec.get_i32()? } else { -1 };
        let surface_type = rec.get_i32()?;
        let num_params = rec.get_i32()?.max(0) as usize;
        let surf_params = rec.get_f64_n(num_params)?;
        surflist.push(SourceSurf {
            id,
            facet_id,
            surface_type,
            num_params,
            surf_params,
        });
    }

    // Extra unhandled records between surfaces and summary.
    if let Some(n) = niwr {
        for _ in 0..n.max(0) {
            let _ = read_record(f)?;
        }
    }

    // Summary record: fixed count plus trailing extras.
    let mipts_v = mipts.unwrap_or(0).max(0) as usize;
    let njsw_v = njsw.max(0) as usize;
    let niwr_v = niwr.map(|v| v.max(0)).unwrap_or(0) as usize;
    let summary_count = (2 + 4 * mipts_v) * (njsw_v + niwr_v) + 1;
    let mut summary_rec = read_record(f)?;
    let mut summary_table = Vec::with_capacity(summary_count);
    for _ in 0..summary_count {
        summary_table.push(summary_rec.get_i32()?);
    }
    let summary_extra = summary_rec.drain_i32_extras()?;

    Ok(SurfSrcHeader {
        kod,
        ver,
        loddat,
        idtm,
        probid,
        aid,
        knod,
        np1: np1_raw.abs(),
        orignp1: np1_raw,
        nrss,
        ncrd,
        njsw,
        niss,
        niwr,
        mipts,
        kjaq,
        table1extra,
        table2extra,
        surflist,
        summary_table,
        summary_extra,
    })
}

// Writer for `put_header` / `put_table_1` / `put_table_2` /
// `put_surface_info` / `put_summary` / `write_tracklist`, built on the shared
// [`RecordSink`] accumulator (`_FortranRecord.put_*`); only the
// version-specific counter widths below are SSW-specific. Every record is
// framed as `[i32 len][payload][i32 len]` little-endian and the version-
// specific counter widths are preserved: MCNPX-2.6.0 stores `np1`/`nrss`
// as i32 while plain MCNP5 and MCNP6 (`SF_00001`) store them as i64 (the
// on-disk MCNP6 layout of two i32 halves is exactly the sign-extended i64).
//
// Deviations from upstream:
// - Table 2 is emitted only when `orignp1 < 0`, keeping writer and reader
//   symmetric. The legacy writer always emits table 2, which only works
//   for files that carry one.
// - Header strings are joined as parsed (padding preserved).
// - Records counted by `niwr > 0` between surfaces and the summary are not
//   reproduced: their contents were discarded at parse time (the legacy
//   reader warns and discards too). All vendored fixtures have `niwr == 0`,
//   so round trips
//   stay byte-exact.

fn put_header(out: &mut Vec<u8>, h: &SurfSrcHeader, layout: &Layout) {
    match layout {
        Layout::SplitHeader => {
            // First record holds only the code identifier.
            let mut rec = RecordSink::new();
            rec.put_str(&h.kod);
            *out = rec.frame();

            // Second record carries everything else plus the dump number.
            let mut rec = RecordSink::new();
            rec.put_str(&h.ver);
            rec.put_str(&h.loddat);
            rec.put_str(&h.idtm);
            rec.put_str(&h.probid);
            rec.put_str(&h.aid);
            rec.put_int(h.knod);
            out.extend_from_slice(&rec.frame());
        }
        Layout::Mcnpx260 | Layout::Mcnp5 => {
            let mut rec = RecordSink::new();
            rec.put_str(&h.kod);
            rec.put_str(&h.ver);
            rec.put_str(&h.loddat);
            rec.put_str(&h.idtm);
            rec.put_str(&h.probid);
            rec.put_str(&h.aid);
            rec.put_int(h.knod);
            *out = rec.frame();
        }
    }
}

fn put_table_1(out: &mut Vec<u8>, h: &SurfSrcHeader, layout: &Layout) {
    let mut rec = RecordSink::new();
    match layout {
        // MCNPX 2.6.0 keeps the two big counters narrow.
        Layout::Mcnpx260 => {
            rec.put_int(h.orignp1 as i32);
            rec.put_int(h.nrss as i32);
        }
        // Plain MCNP5 stores i64; MCNP6's split i32 halves are the sign
        // extension of the same i64 value, so writing i64 reproduces both.
        Layout::SplitHeader | Layout::Mcnp5 => {
            rec.put_long(h.orignp1);
            rec.put_long(h.nrss);
        }
    }
    rec.put_int(h.ncrd);
    rec.put_int(h.njsw);
    rec.put_int(h.niss as i32); // MCNP needs 'int', could be 'long'? (upstream)
    for e in &h.table1extra {
        rec.put_int(*e);
    }
    out.extend_from_slice(&rec.frame());
}

/// Writes the optional table-2 record; errors if the header promises one but
/// its fields are missing.
fn put_table_2(out: &mut Vec<u8>, h: &SurfSrcHeader) -> Result<(), Error> {
    let (Some(niwr), Some(mipts), Some(kjaq)) = (h.niwr, h.mipts, h.kjaq) else {
        return Err(Error::MissingTable2);
    };
    let mut rec = RecordSink::new();
    rec.put_int(niwr);
    rec.put_int(mipts);
    rec.put_int(kjaq);
    for e in &h.table2extra {
        rec.put_int(*e);
    }
    out.extend_from_slice(&rec.frame());
    Ok(())
}

fn put_surface_info(out: &mut Vec<u8>, h: &SurfSrcHeader) {
    for s in &h.surflist {
        let mut rec = RecordSink::new();
        rec.put_int(s.id);
        if h.kjaq == Some(1) {
            rec.put_int(s.facet_id); // macrobody facet flag present
        }
        rec.put_int(s.surface_type);
        rec.put_int(s.num_params as i32);
        for p in &s.surf_params {
            rec.put_double(*p);
        }
        out.extend_from_slice(&rec.frame());
    }
}

fn put_summary(out: &mut Vec<u8>, h: &SurfSrcHeader) {
    let mut rec = RecordSink::new();
    for v in h.summary_table.iter().chain(h.summary_extra.iter()) {
        rec.put_int(*v);
    }
    out.extend_from_slice(&rec.frame());
}

/// Which on-disk header layout a file uses; mirrors the reader branching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// MCNP6 `SF_00001`: header split over two records.
    SplitHeader,
    /// MCNPX 2.6.0: single header record, i32 counters.
    Mcnpx260,
    /// Plain MCNP5: single header record, i64 counters.
    Mcnp5,
}

fn detect_layout(h: &SurfSrcHeader) -> Result<Layout, Error> {
    if h.kod.contains("SF_00001") {
        Ok(Layout::SplitHeader)
    } else if h.ver.contains("2.6.0") {
        Ok(Layout::Mcnpx260)
    } else if h.ver.contains('5') {
        Ok(Layout::Mcnp5)
    } else {
        Err(Error::UnsupportedVersion(h.ver.trim().to_string()))
    }
}

impl SurfSrcHeader {
    /// Serialize just the header block (all records up to and including the
    /// summary), validating version support and table-2 availability.
    pub fn header_block(&self) -> Result<Vec<u8>, Error> {
        let layout = detect_layout(self)?;
        if self.orignp1 < 0 && !self.has_table2() {
            return Err(Error::MissingTable2);
        }
        let mut out = Vec::new();
        emit_block(&mut out, self, &layout)?;
        Ok(out)
    }

    /// True when every table-2 field was parsed from the header.
    fn has_table2(&self) -> bool {
        self.niwr.is_some() && self.mipts.is_some() && self.kjaq.is_some()
    }
}

fn validate_for_write(h: &SurfSrcHeader) -> Result<(), Error> {
    detect_layout(h)?;
    if h.orignp1 < 0 && !h.has_table2() {
        return Err(Error::MissingTable2);
    }
    Ok(())
}

fn emit_block(out: &mut Vec<u8>, h: &SurfSrcHeader, layout: &Layout) -> Result<(), Error> {
    put_header(out, h, layout);
    put_table_1(out, h, layout);
    if h.orignp1 < 0 {
        put_table_2(out, h)?;
    }
    put_surface_info(out, h);
    put_summary(out, h);
    Ok(())
}

/// Write a complete SSW file (header block plus tracklist) to a stream.
///
/// Implements the reference `write_header` + `write_tracklist` behavior. The
/// tracklist length must equal `nrss` and every track record width must equal
/// `abs(ncrd)`; otherwise a validation error is returned before any bytes are
/// produced.
pub fn write_to<W: Write>(
    w: &mut W,
    header: &SurfSrcHeader,
    tracks: &[TrackData],
) -> Result<(), Error> {
    let bytes = encode_file(header, tracks)?;
    w.write_all(&bytes).map_err(|e| Error::Io(e.to_string()))
}

/// Convenience wrapper around [`write_to`] that writes to a file path.
pub fn write_to_path<P: AsRef<Path>>(
    path: P,
    header: &SurfSrcHeader,
    tracks: &[TrackData],
) -> Result<(), Error> {
    let bytes = encode_file(header, tracks)?;
    std::fs::write(path, bytes).map_err(|e| Error::Io(e.to_string()))
}

fn encode_file(header: &SurfSrcHeader, tracks: &[TrackData]) -> Result<Vec<u8>, Error> {
    validate_for_write(header)?;
    let ncrd_abs = header.ncrd.unsigned_abs() as usize;
    if tracks.len() as u64 != header.nrss.max(0) as u64 {
        return Err(Error::TrackCountMismatch {
            expected: header.nrss.max(0) as u64,
            found: tracks.len(),
        });
    }
    for (i, t) in tracks.iter().enumerate() {
        if t.record.len() != ncrd_abs {
            return Err(Error::TrackRecordWidth {
                index: i,
                expected: ncrd_abs,
                found: t.record.len(),
            });
        }
    }

    let layout = detect_layout(header)?;
    let mut out = Vec::new();
    emit_block(&mut out, header, &layout)?;

    // Track records: `ncrd` doubles each, in particle order.
    for t in tracks {
        let mut rec = RecordSink::new();
        for v in &t.record {
            rec.put_double(*v);
        }
        out.extend_from_slice(&rec.frame());
    }
    Ok(out)
}

// ── Combine ───────────────────────────────────────────────────────────────
//
// Port of `scripts/ssw_combine.py::combine_multiple_ss_files`: merges several
// SSW files into one with a combined header. Header compatibility follows
// `_compare_compatible` (kod/ver/loddat, mipts, njsw/niwr/kjaq, per-surface
// records); unlike `cmp_semantics` — and like upstream — counters
// (`np1`/`nrss`/`niss`/`ncrd`) are *not* compared.
//
// Combined header rules (upstream `:65-66`):
// - `orignp1` is the *signed* sum of the inputs' `orignp1` (note `orignp1`,
//   not the absolute `np1`); `np1` is its absolute value. `nrss` is the plain
//   sum. Everything else comes from the first file.
// - Track `nps` values of files after the first shift by the cumulative
//   *absolute* `np1` (`:46-49`), preserving each track's sign via `copysign`
//   (`:98-104`).
//
// Deviations from upstream (documented, not silent):
// - Incompatible headers are a typed [`Error::Incompatible`] instead of
//   printing and returning `False`.
// - Table 2 follows this crate's `orignp1 < 0` rule, not the legacy
//   always-emit: combining table-2-carrying files whose signed `orignp1` sum
//   is non-negative is rejected, since the writer could not reproduce the
//   table-2 record the inputs carried.
// - `njsw..njsw+niwr` extras are unsupported on both sides (upstream prints
//   "unsupported entries" and drops them too); combining files with
//   `niwr != 0` is rejected so no combine silently discards records.

/// Combine several SSW files into one (`scripts/ssw_combine.py` port).
///
/// Reads every input header, rejects incompatible sets, then writes `output`
/// with the summed header and concatenated (nps-shifted) tracklists.
pub fn combine_files<P: AsRef<Path>>(inputs: &[P], output: impl AsRef<Path>) -> Result<(), Error> {
    if inputs.is_empty() {
        return Err(Error::Incompatible("need at least one input file".into()));
    }
    let srcs: Vec<SurfSrc> = inputs
        .iter()
        .map(SurfSrc::open)
        .collect::<Result<Vec<_>, _>>()?;
    let first = &srcs[0].header;
    for other in srcs.iter().skip(1) {
        check_compatible(first, &other.header)?;
    }
    if first.niwr.unwrap_or(0) != 0 {
        return Err(Error::Incompatible(format!(
            "niwr = {} extras are unsupported for combine (would be dropped)",
            first.niwr.unwrap_or(0)
        )));
    }

    let mut header = first.clone();
    let orignp1_sum: i64 = srcs.iter().map(|s| s.header.orignp1).sum();
    let carried_table2 = srcs.iter().any(|s| s.header.orignp1 < 0);
    if carried_table2 && orignp1_sum >= 0 {
        return Err(Error::Incompatible(format!(
            "signed orignp1 sum {orignp1_sum} >= 0 cannot carry the inputs' table-2 record"
        )));
    }
    header.orignp1 = orignp1_sum;
    header.np1 = orignp1_sum.abs();
    header.nrss = srcs.iter().map(|s| s.header.nrss).sum();

    // Cumulative absolute-np1 offsets; the first file's tracks are unshifted.
    let mut tracks = Vec::with_capacity(header.nrss.max(0) as usize);
    let mut offset: i64 = 0;
    for (file_idx, src) in srcs.iter().enumerate() {
        for mut t in src.read_tracklist()? {
            if file_idx > 0 {
                let shifted = t.nps + (offset as f64).copysign(t.nps);
                if !t.record.is_empty() {
                    t.record[0] = shifted;
                }
                t.nps = shifted;
            }
            tracks.push(t);
        }
        offset += src.header.np1;
    }

    write_to_path(output, &header, &tracks)
}

/// Header compatibility for combining (`_compare_compatible`).
fn check_compatible(first: &SurfSrcHeader, other: &SurfSrcHeader) -> Result<(), Error> {
    let bad = |what: &str| Error::Incompatible(what.to_string());
    if other.kod != first.kod {
        return Err(bad(&format!(
            "kod {:?} != {:?}",
            other.kod.trim_end(),
            first.kod.trim_end()
        )));
    }
    if other.ver != first.ver {
        return Err(bad(&format!(
            "ver {:?} != {:?}",
            other.ver.trim_end(),
            first.ver.trim_end()
        )));
    }
    if other.loddat != first.loddat {
        return Err(bad("loddat mismatch"));
    }
    if other.mipts != first.mipts {
        return Err(bad("mipts (particle type) mismatch"));
    }
    if other.njsw != first.njsw {
        return Err(bad("njsw (surface count) mismatch"));
    }
    if other.niwr != first.niwr {
        return Err(bad("niwr (cell count) mismatch"));
    }
    if other.kjaq != first.kjaq {
        return Err(bad("kjaq (macrobody facet flag) mismatch"));
    }
    if other.surflist.len() != first.surflist.len() {
        return Err(bad("surface list length mismatch"));
    }
    for (n, (a, b)) in first.surflist.iter().zip(other.surflist.iter()).enumerate() {
        if a.id != b.id
            || a.facet_id != b.facet_id
            || a.surface_type != b.surface_type
            || a.num_params != b.num_params
            || a.surf_params != b.surf_params
        {
            return Err(bad(&format!("surface record {n} mismatch")));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/ssw/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn mcnp5_header_matches_oracle() {
        let s = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap();
        let h = &s.header;
        assert_eq!(h.kod, "mcnp    ");
        assert_eq!(h.ver, "5    ");
        assert_eq!(h.loddat, "01232009");
        assert_eq!(h.idtm, " 10/31/11 13:52:39 ");
        assert_eq!(h.probid, " 10/31/11 13:52:35 ");
        assert_eq!(
            h.aid,
            "c Test deck with H20 cube, point n source, SSW of top surface interactions      "
        );
        assert_eq!(h.knod, 2);
        assert_eq!(h.np1, 1000);
        assert_eq!(h.nrss, 173);
        assert_eq!(h.ncrd, 11);
        assert_eq!(h.njsw, 1);
        assert_eq!(h.niss, 173);
        // Table 2 present (orignp1 < 0 upstream).
        assert_eq!(h.orignp1, -1000);
        assert_eq!(h.niwr, Some(0));
        assert_eq!(h.mipts, Some(3));
        assert_eq!(h.kjaq, Some(0));
        assert_eq!(h.surflist.len(), 1);
    }

    #[test]
    fn mcnp6_header_matches_oracle() {
        let s = SurfSrc::open(fixture("mcnp6_surfsrc.w")).unwrap();
        let h = &s.header;
        assert_eq!(h.kod, "SF_00001");
        assert_eq!(h.ver, "mcnp    6   ");
        assert_eq!(h.loddat, " 05/08/13");
        assert_eq!(h.idtm, " 11/18/13 17:50:49 ");
        assert_eq!(h.probid, " 11/18/13 17:50:43 ");
        assert!(h.aid.starts_with("Simple MCNP Example that uses SSW"));
        assert_eq!(h.knod, 2);
        assert_eq!(h.np1, 10000);
        assert_eq!(h.nrss, 1710);
        assert_eq!(h.ncrd, -11);
        assert_eq!(h.njsw, 1);
        assert_eq!(h.niss, 1701);
        assert_eq!(h.mipts, Some(37));
        assert_eq!(h.kjaq, Some(0));
    }

    #[test]
    fn mcnpx_header_matches_oracle() {
        let s = SurfSrc::open(fixture("mcnpx_surfsrc.w")).unwrap();
        let h = &s.header;
        assert_eq!(h.kod, "mcnpx   ");
        assert_eq!(h.ver, "2.6.0");
        assert_eq!(h.loddat, "Wed Apr 09 08:00:00 MST 2008");
        assert_eq!(h.idtm, "  10/28/13 02:16:22");
        assert_eq!(h.probid, "  10/28/13 02:16:16");
        assert!(h.aid.starts_with("Simple MCNP Example that uses SSW"));
        assert_eq!(h.knod, 2);
        assert_eq!(h.np1, 10000);
        assert_eq!(h.nrss, 1658);
        assert_eq!(h.ncrd, 11);
        assert_eq!(h.njsw, 1);
        assert_eq!(h.niss, 1652);
        assert_eq!(h.mipts, Some(35));
        assert_eq!(h.kjaq, Some(0));
    }

    #[test]
    fn onetrack_tracklist_values_match_oracle() {
        let s = SurfSrc::open(fixture("mcnp_surfsrc_onetrack.w")).unwrap();
        let tracks = s.read_tracklist().unwrap();
        assert_eq!(tracks.len(), 1);
        let t = &tracks[0];
        assert_eq!(t.nps, 1.0);
        assert!(((t.bitarray - 8.000048e6) / 8.000048e6).abs() < 1e-6);
        assert!((t.wgt - 0.99995639).abs() < 1e-7);
        assert!((t.erg - 5.54203947).abs() < 1e-7);
        assert!((t.tme - 0.17144023).abs() < 1e-7);
        assert!((t.x - (-8.05902e-2)).abs() < 1e-7);
        assert!((t.y - 3.122666098).abs() < 1e-7);
        assert!((t.z - 5.0).abs() < 1e-7);
        assert!((t.u - (-0.35133163)).abs() < 1e-7);
        assert!((t.v - 0.48465036).abs() < 1e-7);
        assert!((t.cs - 0.80104937).abs() < 1e-7);
        assert!((t.w - 0.80104937).abs() < 1e-7);
    }

    #[test]
    fn track_count_matches_nrss_for_mcnp5() {
        let s = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap();
        let tracks = s.read_tracklist().unwrap();
        assert_eq!(tracks.len(), 173);
        // First track spot-check following the printed-data conventions.
        assert_eq!(tracks[0].nps as i64, tracks[0].record[0] as i64);
    }

    #[test]
    fn compare_same_file_equal_cross_version_unequal() {
        let a = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap();
        let b = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap();
        assert_eq!(a.header.cmp_semantics(&b.header), Ordering::Equal);

        let c = SurfSrc::open(fixture("mcnp6_surfsrc.w")).unwrap();
        assert_ne!(a.header.cmp_semantics(&c.header), Ordering::Equal);
        let d = SurfSrc::open(fixture("mcnpx_surfsrc.w")).unwrap();
        assert_ne!(a.header.cmp_semantics(&d.header), Ordering::Equal);
    }

    #[test]
    fn print_header_contains_expected_substrings() {
        let s = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap();
        let txt = s.header.print_header();
        assert!(txt.starts_with("Code: mcnp     (version: 5    ) [01232009]"));
        assert!(txt.contains("Showing dump #2"));
        assert!(txt.contains("1000 histories, 173 tracks, 11 record size"));
        // Surface record: id=6, facet=-1 (kjaq=0), type=4, one param (5.0).
        assert!(txt.contains("Surface 6: facet -1, type 4 with 1 parameters: ( 5)"));
    }

    #[test]
    fn truncated_file_errors_not_panics() {
        let full = std::fs::read(fixture("mcnp5_surfsrc.w")).unwrap();
        let dir = std::env::temp_dir();
        let p = dir.join("nucleide_trunc_ssw.w");
        std::fs::write(&p, &full[..full.len() / 3]).unwrap();
        let r = SurfSrc::open(&p);
        if let Ok(s) = r {
            // Header may survive truncation; tracklist must not panic.
            let _ = s.read_tracklist();
        }
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn cell_decode_from_bitarray() {
        let mut rec = vec![0.0; 11];
        rec[1] = 8.000048e6;
        let t = build_track(rec);
        // abs(bitarray)//8 % 1e8 == 1000006
        assert_eq!(t.cell(), 1_000_006.0);
    }

    // ── Writer round trips ──────────────────────────────────────────────

    /// Parse a fixture, re-serialize header + tracklist, and require the
    /// output bytes to equal the original file exactly.
    fn assert_round_trip_bytes(name: &str) {
        let path = fixture(name);
        let original = std::fs::read(&path).unwrap();
        let s = SurfSrc::open(&path).unwrap();
        let tracks = s.read_tracklist().unwrap();

        let dir = std::env::temp_dir();
        let out = dir.join(format!("nucleide_rt_{name}"));
        write_to_path(&out, &s.header, &tracks).unwrap();
        let rewritten = std::fs::read(&out).unwrap();

        // Byte-level check: parsed fields alone must reproduce every byte
        // (all vendored fixtures carry niwr == 0, so no discarded records).
        assert_eq!(rewritten.len(), original.len(), "{name}: length mismatch");
        assert_eq!(rewritten, original, "{name}: byte mismatch");

        // Semantic check: re-parse what we wrote and compare via the
        // comparison semantics plus full track equality.
        let reparsed = SurfSrc::open(&out).unwrap();
        assert_eq!(
            s.header.cmp_semantics(&reparsed.header),
            Ordering::Equal,
            "{name}: cmp_semantics diverged"
        );
        assert_eq!(reparsed.header, s.header, "{name}: header fields diverged");
        let rt_tracks = reparsed.read_tracklist().unwrap();
        assert_eq!(rt_tracks, tracks, "{name}: tracklists diverged");
        std::fs::remove_file(&out).ok();
    }

    #[test]
    fn round_trip_mcnp5_bytes_identical() {
        assert_round_trip_bytes("mcnp5_surfsrc.w");
    }

    #[test]
    fn round_trip_mcnp6_bytes_identical() {
        assert_round_trip_bytes("mcnp6_surfsrc.w");
    }

    #[test]
    fn round_trip_mcnpx_bytes_identical() {
        assert_round_trip_bytes("mcnpx_surfsrc.w");
    }

    #[test]
    fn round_trip_onetrack_bytes_identical() {
        assert_round_trip_bytes("mcnp_surfsrc_onetrack.w");
    }

    #[test]
    fn header_block_alone_round_trips() {
        for name in [
            "mcnp5_surfsrc.w",
            "mcnp6_surfsrc.w",
            "mcnpx_surfsrc.w",
            "mcnp_surfsrc_onetrack.w",
        ] {
            let original = std::fs::read(fixture(name)).unwrap();
            let s = SurfSrc::open(fixture(name)).unwrap();
            let block = s.header.header_block().unwrap();
            let nrss = s.header.nrss as usize;
            let ncrd = s.header.ncrd.unsigned_abs() as usize;
            let expected_len = block.len() + nrss * (8 + 8 * ncrd);
            assert_eq!(
                original.len(),
                expected_len,
                "{name}: header block + nrss*track records must span the file"
            );
            assert_eq!(&original[..block.len()], &block[..], "{name}");
        }
    }

    #[test]
    fn writer_rejects_unsupported_version() {
        let mut h = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap().header;
        h.kod = "mcnp4   ".into();
        h.ver = "4    ".into();
        let err =
            write_to_path(std::env::temp_dir().join("nucleide_bad_ver.w"), &h, &[]).unwrap_err();
        assert_eq!(err, Error::UnsupportedVersion("4".into()));
    }

    #[test]
    fn writer_rejects_missing_table2_fields() {
        let mut h = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap().header;
        assert!(h.orignp1 < 0);
        h.niwr = None; // strip table-2 state behind the writer's back
        let err = h.header_block().unwrap_err();
        assert_eq!(err, Error::MissingTable2);
    }

    #[test]
    fn writer_rejects_track_count_mismatch() {
        let s = SurfSrc::open(fixture("mcnp_surfsrc_onetrack.w")).unwrap();
        let mut tracks = s.read_tracklist().unwrap();
        tracks.clear(); // nrss == 1, no tracks supplied
        let err = write_to(&mut std::io::sink(), &s.header, &tracks).unwrap_err();
        assert_eq!(
            err,
            Error::TrackCountMismatch {
                expected: 1,
                found: 0
            }
        );
    }

    #[test]
    fn writer_rejects_wrong_record_width() {
        let s = SurfSrc::open(fixture("mcnp_surfsrc_onetrack.w")).unwrap();
        let mut tracks = s.read_tracklist().unwrap();
        tracks[0].record.truncate(10); // ncrd == 11
        let err = write_to(&mut std::io::sink(), &s.header, &tracks).unwrap_err();
        assert_eq!(
            err,
            Error::TrackRecordWidth {
                index: 0,
                expected: 11,
                found: 10
            }
        );
    }

    // ── Combine (`ssw_combine.py` port) ───────────────────────────────────

    #[test]
    fn combine_onetrack_with_itself_sums_signed_headers() {
        let dir = std::env::temp_dir();
        let out = dir.join("nucleide_combine_pair.w");
        combine_files(
            &[
                fixture("mcnp_surfsrc_onetrack.w"),
                fixture("mcnp_surfsrc_onetrack.w"),
            ],
            &out,
        )
        .unwrap();

        let first = SurfSrc::open(fixture("mcnp_surfsrc_onetrack.w")).unwrap();
        let merged = SurfSrc::open(&out).unwrap();
        // orignp1 sums signed; np1 is its absolute value; nrss sums plain.
        assert_eq!(
            merged.header.orignp1,
            first.header.orignp1 + first.header.orignp1
        );
        assert_eq!(merged.header.np1, 2 * first.header.np1);
        assert_eq!(merged.header.nrss, 2 * first.header.nrss);

        // First file's tracks verbatim; later tracks shift nps by the
        // cumulative absolute np1, preserving each track's sign.
        let before = first.read_tracklist().unwrap();
        let after = merged.read_tracklist().unwrap();
        assert_eq!(after.len(), 2);
        assert_eq!(after[0], before[0]);
        let want = before[0].nps + (first.header.np1 as f64).copysign(before[0].nps);
        assert_eq!(after[1].nps, want);
        assert_eq!(after[1].record[0], want);

        // The merged file re-parses and round-trips byte-exactly.
        let tracks = merged.read_tracklist().unwrap();
        let dir2 = std::env::temp_dir();
        let rt = dir2.join("nucleide_combine_pair_rt.w");
        write_to_path(&rt, &merged.header, &tracks).unwrap();
        assert_eq!(std::fs::read(&rt).unwrap(), std::fs::read(&out).unwrap());
        std::fs::remove_file(&out).ok();
        std::fs::remove_file(&rt).ok();
    }

    #[test]
    fn combine_single_file_is_identity() {
        let dir = std::env::temp_dir();
        let out = dir.join("nucleide_combine_single.w");
        combine_files(&[fixture("mcnp5_surfsrc.w")], &out).unwrap();
        assert_eq!(
            std::fs::read(&out).unwrap(),
            std::fs::read(fixture("mcnp5_surfsrc.w")).unwrap()
        );
        std::fs::remove_file(&out).ok();
    }

    #[test]
    fn combine_incompatible_versions_rejected() {
        let dir = std::env::temp_dir();
        let out = dir.join("nucleide_combine_bad.w");
        let err = combine_files(
            &[fixture("mcnp5_surfsrc.w"), fixture("mcnpx_surfsrc.w")],
            &out,
        )
        .unwrap_err();
        assert!(matches!(err, Error::Incompatible(_)));
        std::fs::remove_file(&out).ok();
    }

    #[test]
    fn combine_needs_at_least_one_input() {
        let dir = std::env::temp_dir();
        let out = dir.join("nucleide_combine_empty.w");
        let err = combine_files::<String>(&[], &out).unwrap_err();
        assert!(matches!(err, Error::Incompatible(_)));
        std::fs::remove_file(&out).ok();
    }

    // ── Synthetic files (hand-built records, no fixtures) ────────────────

    fn synth_rec_i32(vals: &[i32]) -> Vec<u8> {
        let mut s = RecordSink::new();
        for v in vals {
            s.put_int(*v);
        }
        s.frame()
    }

    /// Parameters for a synthetic plain-MCNP5-layout SSW file, built from
    /// hand-framed little-endian records.
    #[derive(Clone)]
    struct SynthSsw {
        ver: String,
        loddat: String,
        orignp1: i64,
        nrss: i64,
        /// `(niwr, mipts, kjaq)`; `None` writes no table-2 record (requires
        /// `orignp1 >= 0`).
        table2: Option<(i32, i32, i32)>,
        /// Records between surfaces and summary; the count must equal `niwr`.
        niwr_extras: Vec<Vec<i32>>,
        surfaces: Vec<SourceSurf>,
        tracks: Vec<Vec<f64>>,
    }

    impl SynthSsw {
        /// A small compatible pair base: one surface, one track, table 2 with
        /// `niwr == 0`.
        fn base() -> Self {
            SynthSsw {
                ver: "5    ".into(),
                loddat: "01012000".into(),
                orignp1: -1000,
                nrss: 1,
                table2: Some((0, 3, 0)),
                niwr_extras: Vec::new(),
                surfaces: vec![SourceSurf {
                    id: 6,
                    facet_id: -1,
                    surface_type: 4,
                    num_params: 1,
                    surf_params: vec![5.0],
                }],
                tracks: vec![(0..11).map(|i| i as f64).collect()],
            }
        }

        fn build(&self) -> Vec<u8> {
            if let Some((niwr, _, _)) = self.table2 {
                assert!(self.orignp1 < 0, "test bug: table 2 requires orignp1 < 0");
                assert_eq!(
                    niwr.max(0) as usize,
                    self.niwr_extras.len(),
                    "test bug: extras must match niwr"
                );
            } else {
                assert!(self.orignp1 >= 0, "test bug: negative np1 requires table 2");
            }
            let (niwr, mipts, kjaq) = self.table2.unwrap_or((0, 0, 0));
            let njsw = self.surfaces.len() as i32;
            let ncrd = self
                .tracks
                .first()
                .map_or(TrackData::RECORD_WIDTH, Vec::len) as i32;

            let mut out = Vec::new();
            // Header record: kod(8) ver(5) loddat(8) idtm(19) probid(19)
            // aid(80) knod — the plain-MCNP5 single-record layout.
            let mut h = RecordSink::new();
            h.put_str("mcnp    ");
            h.put_str(&self.ver);
            h.put_str(&self.loddat);
            h.put_str(&" ".repeat(19));
            h.put_str(&" ".repeat(19));
            h.put_str(&" ".repeat(80));
            h.put_int(2); // knod
            out.extend_from_slice(&h.frame());
            // Table 1: plain MCNP5 stores np1/nrss as i64.
            let mut t1 = RecordSink::new();
            t1.put_long(self.orignp1);
            t1.put_long(self.nrss);
            t1.put_int(ncrd);
            t1.put_int(njsw);
            t1.put_int(self.nrss as i32); // niss
            out.extend_from_slice(&t1.frame());
            if self.orignp1 < 0 {
                let mut t2 = RecordSink::new();
                t2.put_int(niwr);
                t2.put_int(mipts);
                t2.put_int(kjaq);
                out.extend_from_slice(&t2.frame());
            }
            for s in &self.surfaces {
                let mut r = RecordSink::new();
                r.put_int(s.id);
                if kjaq == 1 {
                    r.put_int(s.facet_id); // macrobody facet flag present
                }
                r.put_int(s.surface_type);
                r.put_int(s.num_params as i32);
                for p in &s.surf_params {
                    r.put_double(*p);
                }
                out.extend_from_slice(&r.frame());
            }
            for extra in &self.niwr_extras {
                out.extend_from_slice(&synth_rec_i32(extra));
            }
            let summary_count = (2 + 4 * mipts.max(0)) * (njsw.max(0) + niwr.max(0)) + 1;
            let mut sum = RecordSink::new();
            for i in 0..summary_count {
                sum.put_int(i);
            }
            out.extend_from_slice(&sum.frame());
            for t in &self.tracks {
                let mut r = RecordSink::new();
                for v in t {
                    r.put_double(*v);
                }
                out.extend_from_slice(&r.frame());
            }
            out
        }
    }

    fn write_temp_ssw(name: &str, bytes: &[u8]) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, bytes).unwrap();
        p.display().to_string()
    }

    #[test]
    fn synth_no_table2_round_trips_byte_exact() {
        let mut spec = SynthSsw::base();
        spec.orignp1 = 1000;
        spec.table2 = None;
        let bytes = spec.build();
        let parsed = SurfSrc::from_bytes(bytes.clone()).unwrap();
        assert_eq!(parsed.header.orignp1, 1000);
        assert_eq!(parsed.header.niwr, None);
        let tracks = parsed.read_tracklist().unwrap();
        let mut out = Vec::new();
        write_to(&mut out, &parsed.header, &tracks).unwrap();
        assert_eq!(out, bytes);
    }

    #[test]
    fn reader_discards_niwr_extra_records() {
        // `niwr > 0` records between the surface list and the summary are
        // discarded at parse (the legacy reader warns and discards too); the
        // writer does not reproduce them — the pinned documented deviation.
        // kjaq = 1 also exercises the macrobody-facet surface-record layout.
        let mut spec = SynthSsw::base();
        spec.table2 = Some((2, 0, 1));
        spec.niwr_extras = vec![vec![777_111], vec![777_222]];
        spec.surfaces[0].facet_id = 2;
        let bytes = spec.build();

        let parsed = SurfSrc::from_bytes(bytes.clone()).unwrap();
        let h = &parsed.header;
        assert_eq!(h.niwr, Some(2));
        assert_eq!(h.mipts, Some(0));
        assert_eq!(h.kjaq, Some(1));
        assert_eq!(h.surflist[0].facet_id, 2);
        // (2 + 4*mipts) * (njsw + niwr) + 1 entries, extras accounted.
        assert_eq!(h.summary_table.len(), 7);

        let tracks = parsed.read_tracklist().unwrap();
        assert_eq!(tracks.len(), 1);
        assert_eq!(
            tracks[0].record,
            (0..11).map(|i| i as f64).collect::<Vec<_>>()
        );

        // The rewrite drops the extra records (the pinned documented
        // deviation): two 12-byte framed records vanish, the header still
        // reports niwr = 2, so the rewritten file deliberately does not
        // re-parse — combine refuses such inputs for exactly this reason.
        let mut rewritten = Vec::new();
        write_to(&mut rewritten, &parsed.header, &tracks).unwrap();
        let extra_framed = 2 * (4 + 4 + 4);
        assert_eq!(rewritten.len(), bytes.len() - extra_framed);
        assert!(SurfSrc::from_bytes(rewritten.clone()).is_err());
        let marker = 777_111i32.to_le_bytes();
        assert!(bytes.windows(4).any(|w| w == marker));
        assert!(!rewritten.windows(4).any(|w| w == marker));
        let marker = 777_222i32.to_le_bytes();
        assert!(bytes.windows(4).any(|w| w == marker));
        assert!(!rewritten.windows(4).any(|w| w == marker));
    }

    #[test]
    fn combine_rejects_niwr_extras() {
        // Compatible inputs carrying niwr extras are rejected rather than
        // silently dropping the unsupported records.
        let mut spec = SynthSsw::base();
        spec.table2 = Some((2, 3, 0));
        spec.niwr_extras = vec![vec![777_111], vec![777_222]];
        let bytes = spec.build();
        let a = write_temp_ssw("nucleide_synth_niwr_a.w", &bytes);
        let b = write_temp_ssw("nucleide_synth_niwr_b.w", &bytes);
        let out = std::env::temp_dir().join("nucleide_synth_niwr_out.w");
        let err = combine_files(&[a, b], &out).unwrap_err();
        assert_eq!(
            err,
            Error::Incompatible(
                "niwr = 2 extras are unsupported for combine (would be dropped)".into()
            )
        );
        assert!(!out.exists());
        for p in ["nucleide_synth_niwr_a.w", "nucleide_synth_niwr_b.w"] {
            std::fs::remove_file(std::env::temp_dir().join(p)).ok();
        }
    }

    #[test]
    fn combine_flags_each_header_mismatch() {
        // Each check_compatible arm, one diverging field at a time.
        let base = SynthSsw::base();
        let mut ver = base.clone();
        ver.ver = "5.1  ".into();
        let mut loddat = base.clone();
        loddat.loddat = "02022002".into();
        let mut mipts = base.clone();
        mipts.table2 = Some((0, 4, 0));
        let mut njsw = base.clone();
        njsw.surfaces.push(SourceSurf {
            id: 7,
            facet_id: -1,
            surface_type: 4,
            num_params: 1,
            surf_params: vec![6.0],
        });
        let mut niwr = base.clone();
        niwr.table2 = Some((1, 3, 0));
        niwr.niwr_extras = vec![vec![777_111]];
        let mut kjaq = base.clone();
        kjaq.table2 = Some((0, 3, 1));
        kjaq.surfaces[0].facet_id = 2;
        let mut params = base.clone();
        params.surfaces[0].surf_params = vec![6.0];
        let cases: Vec<(SynthSsw, &str)> = vec![
            (ver, "ver \"5.1\" != \"5\""),
            (loddat, "loddat mismatch"),
            (mipts, "mipts (particle type) mismatch"),
            (njsw, "njsw (surface count) mismatch"),
            (niwr, "niwr (cell count) mismatch"),
            (kjaq, "kjaq (macrobody facet flag) mismatch"),
            (params, "surface record 0 mismatch"),
        ];
        for (i, (case, msg)) in cases.into_iter().enumerate() {
            let a = write_temp_ssw(&format!("nucleide_synth_mm{i}_a.w"), &base.build());
            let b = write_temp_ssw(&format!("nucleide_synth_mm{i}_b.w"), &case.build());
            let out = std::env::temp_dir().join(format!("nucleide_synth_mm{i}_out.w"));
            let err = combine_files(&[a, b], &out).unwrap_err();
            assert_eq!(err, Error::Incompatible(msg.into()), "case {i}");
            assert!(!out.exists());
            for suffix in ["a.w", "b.w"] {
                std::fs::remove_file(
                    std::env::temp_dir().join(format!("nucleide_synth_mm{i}_{suffix}")),
                )
                .ok();
            }
        }
    }

    #[test]
    fn reader_rejects_unsupported_version() {
        let mut spec = SynthSsw::base();
        spec.ver = "4    ".into();
        let err = SurfSrc::from_bytes(spec.build()).unwrap_err();
        assert_eq!(err, Error::UnsupportedVersion("4".into()));
    }

    #[test]
    fn big_endian_framing_fails_loudly() {
        // A big-endian framed header record: the LE lead read is 8 << 24, so
        // parsing fails loudly (ShortRecord) instead of silently accepting
        // swapped bytes. The framework is little-endian only by construction.
        let mut be = vec![0, 0, 0, 8];
        be.extend_from_slice(b"mcnp    ");
        be.extend_from_slice(&[0, 0, 0, 8]);
        be.extend_from_slice(&[0u8; 64]);
        let err = SurfSrc::from_bytes(be).unwrap_err();
        assert!(matches!(err, Error::ShortRecord { .. }), "{err:?}");
    }

    #[test]
    fn fortran_errors_convert_to_typed_variants() {
        // Empty input: the leading marker read fails as I/O.
        assert!(matches!(SurfSrc::from_bytes(Vec::new()), Err(Error::Io(_))));
        // Declared 16-byte payload with only 8 present.
        let mut short = 16i32.to_le_bytes().to_vec();
        short.extend_from_slice(&[0u8; 8]);
        assert!(matches!(
            SurfSrc::from_bytes(short),
            Err(Error::ShortRecord { .. })
        ));
        // Lead/trailer disagreement on the header record converts verbatim.
        let mut bad = SynthSsw::base().build();
        let header_payload = 8 + 5 + 8 + 19 + 19 + 80 + 4;
        let last_trailer_byte = header_payload + 4 + 4 - 1;
        bad[last_trailer_byte] ^= 0xFF;
        assert!(matches!(
            SurfSrc::from_bytes(bad),
            Err(Error::BadRecordMarker { .. })
        ));
    }

    #[test]
    fn track_data_from_record_matches_builder() {
        let rec: Vec<f64> = (0..TrackData::RECORD_WIDTH).map(|i| i as f64).collect();
        let t = TrackData::from_record(rec.clone());
        assert_eq!(t.nps, rec[0]);
        assert_eq!(t.bitarray, rec[1]);
        assert_eq!(t.record, rec);
        let w = (1.0 - rec[8] * rec[8] - rec[9] * rec[9])
            .max(0.0)
            .sqrt()
            .copysign(rec[1]);
        assert_eq!(t.w, w);
    }

    #[test]
    fn cmp_semantics_diverges_on_surface_fields() {
        let a = SurfSrc::open(fixture("mcnp5_surfsrc.w")).unwrap().header;
        let mut b = a.clone();
        b.surflist[0].surf_params[0] = 99.0;
        assert_ne!(a.cmp_semantics(&b), Ordering::Equal);
        let mut c = a.clone();
        c.surflist[0].id = 7;
        assert_ne!(a.cmp_semantics(&c), Ordering::Equal);
    }

    #[test]
    fn writer_validate_rejects_stripped_table2_fields() {
        // Same strip as writer_rejects_missing_table2_fields, but through
        // write_to so validate_for_write (not header_block) raises it.
        let parsed = SurfSrc::from_bytes(SynthSsw::base().build()).unwrap();
        let mut h = parsed.header.clone();
        assert!(h.orignp1 < 0);
        h.niwr = None;
        let tracks = parsed.read_tracklist().unwrap();
        let err = write_to(&mut std::io::sink(), &h, &tracks).unwrap_err();
        assert_eq!(err, Error::MissingTable2);
    }

    #[test]
    fn write_to_surfaces_io_errors() {
        struct FailingWriter;
        impl std::io::Write for FailingWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                let _ = buf;
                Err(std::io::Error::other("disk full"))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::other("disk full"))
            }
        }
        let parsed = SurfSrc::from_bytes(SynthSsw::base().build()).unwrap();
        let tracks = parsed.read_tracklist().unwrap();
        assert!(FailingWriter.flush().is_err());
        match write_to(&mut FailingWriter, &parsed.header, &tracks) {
            Err(Error::Io(m)) => assert!(m.contains("disk full"), "{m}"),
            other => panic!("expected Io error, got {other:?}"),
        }
    }

    #[test]
    fn error_display_messages() {
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
        assert!(format!("{}", Error::UnsupportedVersion("5".into()))
            .contains("MCNP version `5` not supported"));
        assert!(format!("{}", Error::MissingTable2).contains("negative orignp1 requires table-2"));
        assert!(format!(
            "{}",
            Error::TrackCountMismatch {
                expected: 1,
                found: 0
            }
        )
        .contains("tracklist holds 0 tracks, header says 1"));
        assert!(format!(
            "{}",
            Error::TrackRecordWidth {
                index: 0,
                expected: 11,
                found: 10
            }
        )
        .contains("track 0 holds 10 values, ncrd says 11"));
        assert!(format!("{}", Error::Incompatible("z".into()))
            .contains("surface-source files incompatible: z"));
    }
}
