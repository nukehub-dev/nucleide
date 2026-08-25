//! MCNP PTRAC event-file reading.
//!
//! Validated against the vendored golden fixtures
//! (`mcnp_ptrac_i{4,8}_little.ptrac`, MCNP6 variants).
//!
//! Format notes (verified against fixture hex dumps):
//! - Fortran records `[i32 len][payload][i32 len]`, little-endian. Record
//!   markers stay 4-byte even in 8-byte-number files; only payload widths
//!   change.
//! - File opens with a sentinel record holding -1; endianness is probed via
//!   its leading marker (must be 4).
//! - The PTRAC *input echo* record is 10 floats; reading it as f32 yields 20
//!   halves on 8-byte files, which triggers the `I8` mode (repack pairs of
//!   f32 into f64).
//! - Variable-count headers follow; NPS lines are integers at the active
//!   width, event lines are floats at the active width.

use std::fmt;
use std::path::Path;

/// Number-width layout detected from the input-echo record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// 4-byte integers and floats.
    I4LittleEndian,
    /// 8-byte integers and floats (record markers still 4-byte).
    I8LittleEndian,
}

/// Errors raised while parsing PTRAC files.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    /// Big-endian or otherwise unsupported layout.
    Unsupported(String),
    /// Payload exhausted mid-record.
    Truncated,
    /// Structural inconsistency (marker mismatch, bad count).
    BadStructure(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::Unsupported(m) => write!(f, "unsupported PTRAC format: {m}"),
            Error::Truncated => write!(f, "file truncated mid-record"),
            Error::BadStructure(m) => write!(f, "malformed PTRAC structure: {m}"),
        }
    }
}

impl std::error::Error for Error {}

/// PTRAC variable-id code → field-name table.
const VARIABLE_MAPPINGS: &[(i32, &str)] = &[
    (1, "nps"),
    (3, "ncl"),
    (4, "nsf"),
    (8, "node"),
    (9, "nsr"),
    (10, "nxs"),
    (11, "ntyn"),
    (12, "nsf"),
    (16, "ipt"),
    (17, "ncl"),
    (18, "mat"),
    (19, "ncp"),
    (20, "xxx"),
    (21, "yyy"),
    (22, "zzz"),
    (23, "uuu"),
    (24, "vvv"),
    (25, "www"),
    (26, "erg"),
    (27, "wgt"),
    (28, "tme"),
];

fn mapping(code: i32) -> Option<&'static str> {
    VARIABLE_MAPPINGS
        .iter()
        .find(|(c, _)| *c == code)
        .map(|(_, n)| *n)
}

/// Per-event-type variable counts.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct VariableNums {
    pub nps: usize,
    pub src: usize,
    pub bnk: usize,
    pub sur: usize,
    pub col: usize,
    pub ter: usize,
}

/// Variable-id lists per event type, in file order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct VariableIds {
    pub nps: Vec<i32>,
    pub src: Vec<i32>,
    pub bnk: Vec<i32>,
    pub sur: Vec<i32>,
    pub col: Vec<i32>,
    pub ter: Vec<i32>,
}

/// One decoded particle event.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Event {
    /// Event class (1000 src, 3000 surface, 4000 collision, 5000 termination,
    /// else bank).
    pub event_type: i32,
    fields: Vec<(&'static str, f64)>,
}

impl Event {
    /// Named field lookup (`"xxx"`, `"erg"`, `"wgt"`, ...).
    pub fn get(&self, name: &str) -> Option<f64> {
        self.fields
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| *v)
    }

    /// Field iterator in variable order.
    pub fn iter(&self) -> impl Iterator<Item = (&'static str, f64)> + '_ {
        self.fields.iter().copied()
    }
}

/// A parsed PTRAC file.
#[derive(Debug, Clone)]
pub struct PtracFile {
    /// Version banner string from record 1.
    pub mcnp_version_info: String,
    /// Problem title card (trimmed).
    pub problem_title: String,
    /// Detected number-width layout.
    pub format: Format,
    pub variable_nums: VariableNums,
    pub variable_ids: VariableIds,
    data: Vec<u8>,
    events_start: usize,
}

impl PtracFile {
    /// Read a file and parse its headers.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let data =
            std::fs::read(path).map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        PtracFile::from_bytes(data)
    }

    /// Parse PTRAC bytes in memory.
    pub fn from_bytes(data: Vec<u8>) -> Result<Self, Error> {
        let mut c = Cursor::new(&data, 0);

        // Sentinel probe. The opening record stores -1 either as i32
        // (4-byte builds: 12-byte record) or i64 (8-byte builds: 16-byte
        // record); markers stay 4-byte in both. Distinguish by checking
        // whether the 4 bytes after a 4-byte payload repeat it.
        let format = probe_sentinel(&c)?;
        c.set_width(format);
        c.pos += match format {
            Format::I4LittleEndian => 12,
            Format::I8LittleEndian => 16,
        };

        // Version banner + title.
        let mcnp_version_info = c.read_string()?;
        let mut problem_title = c.read_string()?;
        while problem_title.ends_with(' ') || problem_title.ends_with('\t') {
            problem_title.pop();
        }
        while problem_title.starts_with(' ') || problem_title.starts_with('\t') {
            problem_title.remove(0);
        }

        let echo_payload_start = c.pos + 4; // skip leading marker
        let mut floats = c.read_float_record()?;

        // Width detection: the input echo holds exactly 10 values at
        // the native width; reading 20 f32 halves on an 8-byte build flips
        // the mode and the original byte pairs get reinterpreted as f64.
        let mut format = format;
        if floats.len() != 10 {
            format = Format::I8LittleEndian;
            let payload = c
                .data
                .get(echo_payload_start..echo_payload_start + floats.len() * 4)
                .ok_or(Error::Truncated)?;
            floats = payload
                .chunks_exact(8)
                .map(|b| f64::from_le_bytes(b.try_into().expect("8 bytes")))
                .collect();
        }
        c.set_width(format);

        // Walk the n-values / values variable-spec structure.
        let num_variables = floats[0] as usize;
        if !(1..=64).contains(&num_variables) {
            return Err(Error::BadStructure(format!(
                "input spec declares {num_variables} variables"
            )));
        }
        let mut current_pos = 1;
        let mut current_variable = 1usize;
        while current_variable <= num_variables {
            if current_pos >= floats.len() {
                floats.extend(c.read_fixed_floats(10)?);
            }
            let n = floats[current_pos] as usize;
            if current_variable < num_variables && current_pos + n + 1 >= floats.len() {
                floats.extend(c.read_fixed_floats(10)?);
            }
            current_pos += n + 1;
            current_variable += 1;
        }

        // Variable-count header.
        let counts: [i64; 11];
        match format {
            Format::I8LittleEndian => {
                let ver_token = mcnp_version_info
                    .chars()
                    .skip(8)
                    .take(5)
                    .collect::<String>()
                    .trim()
                    .to_string();
                if ver_token == "6" || ver_token == "6.mpi" {
                    counts = c.read_counts_mixed(true)?;
                } else {
                    counts = c.read_counts_mixed(false)?;
                }
            }
            Format::I4LittleEndian => {
                let raw = c.read_fixed_i32(20)?;
                counts = core::array::from_fn(|i| raw[i] as i64);
            }
        }

        let variable_nums = VariableNums {
            nps: counts[0] as usize,
            src: (counts[1] + counts[2]) as usize,
            bnk: (counts[3] + counts[4]) as usize,
            sur: (counts[5] + counts[6]) as usize,
            col: (counts[7] + counts[8]) as usize,
            ter: (counts[9] + counts[10]) as usize,
        };

        // Variable-id list.
        let all_ids: Vec<i32> = match format {
            Format::I4LittleEndian => {
                let total: usize = counts[..11].iter().map(|v| *v as usize).sum();
                c.read_fixed_i32(total)?.to_vec()
            }
            Format::I8LittleEndian => {
                let n_q = counts[0] as usize;
                let n_i: usize = counts[1..11].iter().map(|v| *v as usize).sum();
                c.read_mixed_ids(n_q, n_i)?
            }
        };
        if all_ids.len() < variable_nums.total() {
            return Err(Error::BadStructure("fewer ids than declared counts".into()));
        }
        let mut ids_iter = all_ids.into_iter();
        let mut drain = |n: usize| -> Vec<i32> { ids_iter.by_ref().take(n).collect() };
        let variable_ids = VariableIds {
            nps: drain(variable_nums.nps),
            src: drain(variable_nums.src),
            bnk: drain(variable_nums.bnk),
            sur: drain(variable_nums.sur),
            col: drain(variable_nums.col),
            ter: drain(variable_nums.ter),
        };

        Ok(PtracFile {
            mcnp_version_info,
            problem_title,
            format,
            variable_nums,
            variable_ids,
            events_start: c.pos,
            data,
        })
    }

    /// Decode every event in the file.
    ///
    /// Loop structure mirrors the reference writer: an NPS integer record primes
    /// each history; float event records then chain via their first value
    /// until it reads 9000 (end-of-history sentinel).
    pub fn events(&self) -> Result<Vec<Event>, Error> {
        let mut c = Cursor::new(&self.data, self.events_start);
        c.set_width(self.format);
        let mut out = Vec::new();

        loop {
            if c.remaining() == 0 {
                break;
            }
            if c.remaining() < 8 {
                return Err(Error::Truncated);
            }
            let nps_line = c.read_int_record()?;
            let mut next = *nps_line
                .get(1)
                .ok_or_else(|| Error::BadStructure("NPS line shorter than 2 entries".into()))?;

            while next != 9000 {
                let evt_line = c.read_float_record()?;
                if evt_line.is_empty() {
                    return Err(Error::BadStructure("empty event line".into()));
                }
                let e = classify(next);
                let ids = match e {
                    "src" => &self.variable_ids.src,
                    "sur" => &self.variable_ids.sur,
                    "col" => &self.variable_ids.col,
                    "ter" => &self.variable_ids.ter,
                    _ => &self.variable_ids.bnk,
                };
                let mut ev = Event {
                    event_type: next as i32,
                    fields: Vec::new(),
                };
                for (i, code) in ids.iter().enumerate().skip(1) {
                    if let Some(name) = mapping(*code) {
                        if let Some(v) = evt_line.get(i) {
                            ev.fields.push((name, *v));
                        }
                    }
                }
                out.push(ev);
                next = evt_line[0] as i64;
            }
        }
        Ok(out)
    }
}

impl VariableNums {
    fn total(&self) -> usize {
        self.nps + self.src + self.bnk + self.sur + self.col + self.ter
    }
}

fn classify(event_type: i64) -> &'static str {
    match event_type {
        1000 => "src",
        3000 => "sur",
        4000 => "col",
        5000 => "ter",
        _ => "bnk",
    }
}

/// Classify the number width from the sentinel record's byte pattern.
fn probe_sentinel(c: &Cursor<'_>) -> Result<Format, Error> {
    let lead = c.peek_i32()?;
    if lead != 4 {
        return Err(Error::Unsupported(format!(
            "expected leading record marker 4, found {lead} (big-endian?)"
        )));
    }
    let d = c.data;
    let p = c.pos;
    let word = |o: usize| -> [u8; 4] {
        d.get(p + o..p + o + 4)
            .and_then(|b| b.try_into().ok())
            .unwrap_or([0; 4])
    };
    let marker = [4u8, 0, 0, 0];
    if word(8) == marker {
        return Ok(Format::I4LittleEndian); // [4][-1 i32][4]
    }
    if d.len() >= p + 16 && word(12) == marker {
        return Ok(Format::I8LittleEndian); // [4][-1 i64][4]
    }
    Err(Error::BadStructure("unrecognised sentinel layout".into()))
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
    int_w: usize,
    flt_w: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8], pos: usize) -> Self {
        Self {
            data,
            pos,
            int_w: 4,
            flt_w: 4,
        }
    }

    fn set_width(&mut self, f: Format) {
        let w = match f {
            Format::I4LittleEndian => 4,
            Format::I8LittleEndian => 8,
        };
        self.int_w = w;
        self.flt_w = w;
    }

    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    fn peek_i32(&self) -> Result<i32, Error> {
        let b = self
            .data
            .get(self.pos..self.pos + 4)
            .ok_or(Error::Truncated)?;
        Ok(i32::from_le_bytes(b.try_into().expect("4 bytes")))
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], Error> {
        let b = self
            .data
            .get(self.pos..self.pos + n)
            .ok_or(Error::Truncated)?;
        self.pos += n;
        Ok(b)
    }

    fn read_marker(&mut self) -> Result<usize, Error> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes(b.try_into().expect("4 bytes")) as usize)
    }

    /// `[len][payload][len]` → payload slice.
    fn record(&mut self) -> Result<&'a [u8], Error> {
        let lead = self.read_marker()?;
        let payload = self.take(lead)?;
        let trailer = self.read_marker()?;
        if trailer != lead {
            return Err(Error::BadStructure(format!(
                "record markers disagree: {lead} vs {trailer}"
            )));
        }
        Ok(payload)
    }

    fn read_string(&mut self) -> Result<String, Error> {
        Ok(String::from_utf8_lossy(self.record()?).into_owned())
    }

    fn read_float_record(&mut self) -> Result<Vec<f64>, Error> {
        let p = self.record()?;
        Ok(p.chunks_exact(self.flt_w)
            .map(|b| match self.flt_w {
                8 => f64::from_le_bytes(b.try_into().expect("8 bytes")),
                _ => f32::from_le_bytes(b.try_into().expect("4 bytes")) as f64,
            })
            .collect())
    }

    fn read_fixed_floats(&mut self, n: usize) -> Result<Vec<f64>, Error> {
        let mut out = Vec::with_capacity(n);
        let payload = self.record()?;
        let want = n * self.flt_w;
        if payload.len() < want {
            return Err(Error::Truncated);
        }
        for chunk in payload[..want].chunks_exact(self.flt_w) {
            out.push(match self.flt_w {
                8 => f64::from_le_bytes(chunk.try_into().expect("8")),
                _ => f32::from_le_bytes(chunk.try_into().expect("4")) as f64,
            });
        }
        Ok(out)
    }

    fn read_int_record(&mut self) -> Result<Vec<i64>, Error> {
        let p = self.record()?;
        Ok(p.chunks_exact(self.int_w)
            .map(|b| match self.int_w {
                8 => i64::from_le_bytes(b.try_into().expect("8 bytes")),
                _ => i32::from_le_bytes(b.try_into().expect("4 bytes")) as i64,
            })
            .collect())
    }

    fn read_fixed_i32(&mut self, n: usize) -> Result<Vec<i32>, Error> {
        let payload = self.record()?;
        if payload.len() < n * 4 {
            return Err(Error::Truncated);
        }
        Ok(payload[..n * 4]
            .chunks_exact(4)
            .map(|b| i32::from_le_bytes(b.try_into().expect("4")))
            .collect())
    }

    /// 8-byte-file counts: MCNP6 puts the first count in 4 bytes then ten
    /// 8-byte counts; others use eleven 8-byte counts; trailing 4-byte
    /// extras fill to the record end either way.
    fn read_counts_mixed(&mut self, mcnp6: bool) -> Result<[i64; 11], Error> {
        let payload = self.record()?;
        let mut out = [0i64; 11];
        if mcnp6 {
            out[0] = i32::from_le_bytes(
                payload
                    .get(0..4)
                    .ok_or(Error::Truncated)?
                    .try_into()
                    .expect("4"),
            ) as i64;
            for (i, slot) in out.iter_mut().skip(1).enumerate() {
                let s = 4 + i * 8;
                *slot = i64::from_le_bytes(
                    payload
                        .get(s..s + 8)
                        .ok_or(Error::Truncated)?
                        .try_into()
                        .expect("8"),
                );
            }
        } else {
            for (i, slot) in out.iter_mut().enumerate() {
                let s = i * 8;
                *slot = i64::from_le_bytes(
                    payload
                        .get(s..s + 8)
                        .ok_or(Error::Truncated)?
                        .try_into()
                        .expect("8"),
                );
            }
        }
        Ok(out)
    }

    /// 8-byte-file id list: first `n_q` ids are 8-byte, remaining `n_i` are
    /// 4-byte.
    fn read_mixed_ids(&mut self, n_q: usize, n_i: usize) -> Result<Vec<i32>, Error> {
        let payload = self.record()?;
        let mut out = Vec::with_capacity(n_q + n_i);
        let q_bytes = n_q * 8;
        if payload.len() < q_bytes + n_i * 4 {
            return Err(Error::Truncated);
        }
        for b in payload[..q_bytes].chunks_exact(8) {
            out.push(i64::from_le_bytes(b.try_into().expect("8")) as i32);
        }
        for b in payload[q_bytes..q_bytes + n_i * 4].chunks_exact(4) {
            out.push(i32::from_le_bytes(b.try_into().expect("4")));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/ptrac/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn i4() -> PtracFile {
        PtracFile::open(fixture("mcnp_ptrac_i4_little.ptrac")).unwrap()
    }

    #[test]
    fn titles_match_oracle_both_widths() {
        // The fixture decks carry a fixed problem-title card; assert its
        // stable prefix/suffix rather than embedding the historical text.
        let expected_prefix = "Generate a well-defined PTRAC file";
        assert_eq!(
            &i4().problem_title[..expected_prefix.len()],
            expected_prefix
        );
        let p8 = PtracFile::open(fixture("mcnp_ptrac_i8_little.ptrac")).unwrap();
        assert_eq!(&p8.problem_title[..expected_prefix.len()], expected_prefix);
        assert!(p8.problem_title.ends_with("test cases"));
    }

    #[test]
    fn formats_detected() {
        assert_eq!(i4().format, Format::I4LittleEndian);
        let p8 = PtracFile::open(fixture("mcnp_ptrac_i8_little.ptrac")).unwrap();
        assert_eq!(p8.format, Format::I8LittleEndian);
    }

    #[test]
    fn variable_counts_sane() {
        let p = i4();
        // From the fixture's count record: nps=2 and the six event types.
        assert!(p.variable_nums.nps >= 2);
        assert!(p.variable_nums.src > 0);
        assert_eq!(p.variable_ids.nps.len(), p.variable_nums.nps);
        let total: usize = p.variable_nums.total();
        let listed = p.variable_ids.nps.len()
            + p.variable_ids.src.len()
            + p.variable_ids.bnk.len()
            + p.variable_ids.sur.len()
            + p.variable_ids.col.len()
            + p.variable_ids.ter.len();
        assert_eq!(total, listed);
    }

    #[test]
    fn first_event_matches_oracle() {
        let p = i4();
        let events = p.events().unwrap();
        assert!(!events.is_empty());
        let e0 = &events[0];
        assert_eq!(e0.event_type, 1000); // src event
        assert_eq!(e0.get("xxx"), Some(0.0));
        assert_eq!(e0.get("yyy"), Some(0.0));
        assert_eq!(e0.get("zzz"), Some(0.0));
    }

    #[test]
    fn i8_twin_decodes_equivalently() {
        let a = i4().events().unwrap();
        let b = PtracFile::open(fixture("mcnp_ptrac_i8_little.ptrac"))
            .unwrap()
            .events()
            .unwrap();
        assert_eq!(a.len(), b.len());
        // Same source deck; the two files were written with different PTRAC
        // cards so variable lists differ — compare the type chain only.
        for (ea, eb) in a.iter().zip(b.iter()) {
            assert_eq!(ea.event_type, eb.event_type);
        }
        // Spot-check equivalent physics on a shared field of the first src.
        assert_eq!(a[0].get("xxx"), b[0].get("xxx"));
    }

    #[test]
    fn mcnp6_serial_parses() {
        // Despite the "i4" in its filename, this fixture carries an
        // 80-byte input echo (10 f64) — content-based detection is what
        // counts, matching the reference reader.
        let p = PtracFile::open(fixture("mcnp6_serial_ptrac_i4_little.ptrac")).unwrap();
        assert_eq!(p.format, Format::I8LittleEndian);
        let events = p.events().unwrap();
        assert!(!events.is_empty());
        assert!(events.iter().all(|e| e.event_type > 0));
    }

    #[test]
    fn truncated_bytes_error_not_panic() {
        let full = std::fs::read(fixture("mcnp_ptrac_i4_little.ptrac")).unwrap();
        for cut in [4usize, 40] {
            let r = PtracFile::from_bytes(full[..full.len() - cut].to_vec());
            match r {
                Ok(p) => {
                    // Header survived; the event stream must not panic.
                    let _ = p.events();
                }
                Err(Error::Truncated) => {}
                Err(other) => panic!("cut={cut} unexpected error: {other:?}"),
            }
        }
    }

    #[test]
    fn big_endian_marker_rejected() {
        let mut data = 4i32.to_be_bytes().to_vec();
        data.extend(std::iter::repeat_n(0u8, 64));
        assert!(matches!(
            PtracFile::from_bytes(data),
            Err(Error::Unsupported(_))
        ));
    }

    #[test]
    fn variable_mapping_table_covers_oracle_codes() {
        assert_eq!(mapping(20), Some("xxx"));
        assert_eq!(mapping(27), Some("wgt"));
        assert_eq!(mapping(28), Some("tme"));
        assert_eq!(mapping(999), None);
    }
}
