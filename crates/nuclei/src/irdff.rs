//! IAEA IRDFF-II dosimetry cross sections in the SAND-II 725-group structure.
//!
//! [IRDFF-II] (International Reactor Dosimetry and Fusion File, January 2020
//! release, data files updated December 9, 2020; the files self-identify as
//! `IRDFF-2.0`) is the IAEA metrology library for foil-activation neutron
//! dosimetry, 0–60 MeV (Trkov et al., Nucl. Data Sheets 163 (2020) 1, open
//! access arXiv:1909.03336). This module parses the pre-grouped
//! `IRDFF-II.g725` member of the official `IRDFF-II_g725.zip` distribution
//! into caller-ready response rows for spectrum unfolding: one cross-section
//! vector per named reaction, each aligned to the 725-group structure, plus
//! the 726 group boundaries in eV.
//!
//! # Distribution contract
//!
//! Nothing from the IAEA files is vendored into this repository (IAEA
//! copyright, permission-based reproduction). The zip is fetched at runtime
//! by the Python layer (`nucleide.data.fetch_irdff`), hash-pinned against
//! the published SHA-256 recorded there and in `python/nucleide/data.py`,
//! and member text is passed here for parsing. Tests use synthetic
//! hand-built sections in the exact layout below.
//!
//! # Release and group structure pins
//!
//! - Release: IRDFF-II, IAEA, January 2020; cross-section files updated
//!   2020-12-09 (zip member dated 2020-12-08). Primary reference Trkov et
//!   al. 2020 (above).
//! - URL: `https://www-nds.iaea.org/IRDFF/IRDFF-II_g725.zip`
//!   (SHA-256 `6ec2b33c0f67bed46d46be062a24ccedaa5ffea9bbba919958da4b1349f48c85`,
//!   pinned in `nucleide.data.IRDFF_SHA256`).
//! - Group structure: the SAND-II 725-group structure
//!   (`1e-5 eV` to `60 MeV`, 726 boundaries), published as
//!   `https://www-nds.iaea.org/IRDFF/SAND_725.egb` (SHA-256
//!   `0ba231a7dbd0a8c72664e1ce7db0bd3fa06502d407c32cbe3b8336cafb189552`,
//!   verified against the energies carried in the data sections during the
//!   0.13.0 cycle; the boundaries are parsed from the data file itself, so
//!   the egb is a reference cross-check, not a runtime fetch).
//! - Offline or hash-mismatch is a loud `RuntimeError` from the fetch layer
//!   naming the failing URL or both digests (mirroring `fetch_fgr15`).
//!
//! # v1 named reaction subset
//!
//! The pack covers the classic foil-activation set of [`V1_REACTIONS`]
//! (eight reactions, all `MF=3` cross-section sections):
//!
//! | name | reaction | MAT | MT | role |
//! |---|---|---|---|---|
//! | `au197_ng` | 197Au(n,g)198Au | 7925 | 102 | thermal/epithermal capture (4.9 eV resonance) |
//! | `in115_ng` | 115In(n,g)116In | 4931 | 102 | thermal/epithermal capture (1.46 eV resonance) |
//! | `u235_nf` | 235U(n,f) | 9228 | 18 | thermal + fast fission |
//! | `u238_nf` | 238U(n,f) | 9237 | 18 | fast fission (incl. subthreshold) |
//! | `fe56_np` | 56Fe(n,p)56Mn | 2631 | 103 | threshold ~2.9 MeV |
//! | `ni58_np` | 58Ni(n,p)58Co | 2825 | 103 | threshold ~0.4 MeV |
//! | `al27_na` | 27Al(n,a)24Na | 1325 | 107 | threshold ~3.2 MeV |
//! | `na23_n2n` | 23Na(n,2n)22Na | 1125 | 16 | threshold ~12.9 MeV |
//!
//! # v2 named extension set
//!
//! [`V2_REACTIONS`] adds the 26 further `MF=3` dosimetry sections below,
//! from the same `IRDFF-II.g725` member of the same pinned zip (same URL +
//! SHA-256 pins, same cache, same loud offline/mismatch errors — the fetch
//! mechanics are unchanged). The list is the classic foil/dosimetry
//! complement of the v1 set: the remaining threshold (n,p)/(n,a)/(n,2n)
//! monitors, the capture foils, the fission standards, and the
//! high-threshold bismuth monitors. Every (MAT, MT) key was verified
//! present as an `MF=3` section of the pinned distribution during the
//! 0.14.0 expansion cycle, and every section parses with the unchanged
//! `MF=3` machinery (full-range capture/fission rows anchor the structure,
//! threshold rows are contiguous runs ending at the top-boundary
//! terminator — no new row shapes, no parser changes).
//!
//! | name | reaction | MAT | MT | role |
//! |---|---|---|---|---|
//! | `f19_n2n` | 19F(n,2n)18F | 925 | 16 | high-threshold monitor |
//! | `b10_na` | 10B(n,a)7Li | 525 | 107 | thermal 1/v alpha monitor |
//! | `mg24_np` | 24Mg(n,p)24Na | 1225 | 103 | threshold monitor |
//! | `al27_np` | 27Al(n,p)27Mg | 1325 | 103 | threshold monitor |
//! | `si28_np` | 28Si(n,p)28Al | 1425 | 103 | threshold monitor |
//! | `p31_np` | 31P(n,p)31Si | 1525 | 103 | threshold monitor |
//! | `s32_np` | 32S(n,p)32P | 1625 | 103 | threshold monitor |
//! | `sc45_ng` | 45Sc(n,g)46Sc | 2125 | 102 | thermal/epithermal capture foil |
//! | `ti46_np` | 46Ti(n,p)46Sc | 2225 | 103 | threshold monitor |
//! | `ti47_np` | 47Ti(n,p)47Sc | 2228 | 103 | threshold monitor |
//! | `ti48_np` | 48Ti(n,p)48Sc | 2231 | 103 | threshold monitor |
//! | `mn55_n2n` | 55Mn(n,2n)54Mn | 2525 | 16 | threshold monitor |
//! | `fe54_np` | 54Fe(n,p)54Mn | 2625 | 103 | threshold monitor |
//! | `co59_ng` | 59Co(n,g)60Co | 2725 | 102 | thermal capture foil |
//! | `co59_np` | 59Co(n,p)59Fe | 2725 | 103 | threshold monitor |
//! | `ni58_n2n` | 58Ni(n,2n)57Ni | 2825 | 16 | high-threshold monitor |
//! | `cu63_na` | 63Cu(n,a)60Co | 2925 | 107 | threshold monitor |
//! | `zn64_np` | 64Zn(n,p)64Cu | 3025 | 103 | threshold monitor |
//! | `in113_ng` | 113In(n,g)114mIn | 4925 | 102 | thermal/epithermal capture foil |
//! | `ta181_ng` | 181Ta(n,g)182Ta | 7328 | 102 | thermal/epithermal capture foil |
//! | `w186_ng` | 186W(n,g)187W | 7443 | 102 | thermal/epithermal capture foil |
//! | `th232_nf` | 232Th(n,f) | 9040 | 18 | fast fission chamber |
//! | `np237_nf` | 237Np(n,f) | 9346 | 18 | fast fission chamber |
//! | `pu239_nf` | 239Pu(n,f) | 9437 | 18 | thermal + fast fission chamber |
//! | `bi209_n2n` | 209Bi(n,2n)208Bi | 8325 | 16 | high-threshold monitor |
//! | `bi209_n3n` | 209Bi(n,3n)207Bi | 8325 | 17 | high-threshold monitor |
//!
//! The full registry is `V1_REACTIONS` ∪ `V2_REACTIONS` (34 reactions).
//! Higher-order bismuth sections (`MT=37/152/153`) and the gas-production,
//! damage, disappearance, and kerma sections (`MT=1/2/101/105/205/207/800/801`)
//! are not dosimetry response rows and stay out of the registry.
//!
//! Isomer-production dosimetry reactions (115In(n,n′)115mIn,
//! 103Rh(n,n′)103mRh) are stored in this distribution as `MF=10` sections
//! in pointwise (not group-aligned) form, so they stay out of the v1 pack;
//! the parser handles `MF=3` sections only and says so loudly.
//!
//! # `MF=3` section layout (as distributed)
//!
//! 80-column ENDF-style records; the trailing control triplet is MAT
//! (cols 67–70), MF (cols 71–72), MT (cols 73–75), sequence (cols 76–80):
//!
//! ```text
//!  79197.0000 195.274000          0          0          0          07925 3  1    1   <- head: ZA, AWR, 0, L2, 0, 0
//!   6512340.00  6512340.00          0          0          1        7267925 3  1    2   <- ctrl1: ..., 0, 0, 1, NG2
//!          726          1                                            7925 3  1    3   <- ctrl2: NG2, 1
//!   1.00000E-5 4894.41045 1.05000E-5 4778.69934 1.10000E-5 4671.979847925 3  1    4   <- (E, sigma) pairs, 3 per line
//!   ...
//!   58000000.0 5.28690E-6 60000000.0        0.0                      7925 3  1  242   <- final pair: (60 MeV, 0.0) terminator
//! ```
//!
//! - `NG2` (ctrl1 field 6) is the exact pair count; ctrl2 repeats it with
//!   a trailing `1`. Fields are eleven columns each; floats may use the
//!   Fortran embedded exponent (`2.589913-5` = `2.589913e-5`).
//! - The head card's third integer field (`L2`) is `0` on most sections
//!   but `99` on three sections of the pinned file (`56Fe(n,p)` MAT 2631,
//!   `197Au(n,2n)` MAT 7925 MT 16, and the natural-boron total MAT 528
//!   MT 1 — only the first is in the registry). The field is unused
//!   downstream; the parser accepts `0` or `99` there and rejects anything
//!   else loudly, like every other head-card deviation.
//! - The pair energies are the SAND-II 725-group **lower boundaries** in
//!   eV. Full-range reactions (capture, fission) list all 726 boundaries;
//!   threshold reactions list only the contiguous run from their first
//!   above-threshold group through 60 MeV (groups below are zero). The
//!   final pair is always `(60 MeV, 0.0)` — a structural terminator, not a
//!   group value.
//! - Parsing is strict and loud: control cards must match the layout,
//!   exactly `NG2` pairs must be present, every value must be finite, and
//!   every energy must lie on the anchor row's group structure
//!   (contiguous, ending at the top boundary). Rows are never silently
//!   truncated or padded.
//!
//! Metrology/screening context: these are evaluated dosimetry cross
//! sections for foil unfolding, not transport-grade data for a specific
//! facility model.
//!
//! [IRDFF-II]: https://www-nds.iaea.org/IRDFF/

use std::collections::{HashMap, HashSet};
use std::fmt;

/// Result alias for IRDFF-II g725 parsing.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors from parsing IRDFF-II `MF=3` group cross-section sections.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// A requested reaction has no `MF=3` section in the file.
    MissingSection {
        /// Registry name of the missing reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
    },
    /// A requested reaction appears twice in the request or in the file.
    DuplicateSection {
        /// Registry name of the duplicated reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
    },
    /// A section head or control card is not the documented six-field
    /// layout, or the two `NG2` copies disagree.
    MalformedControl {
        /// Registry name of the offending reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// 1-based source line number.
        line: usize,
        /// The offending line text.
        text: String,
    },
    /// The section holds fewer or more `(E, sigma)` pairs than `NG2`.
    BadPairCount {
        /// Registry name of the offending reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// Declared pair count.
        expected: usize,
        /// Parsed pair count.
        found: usize,
    },
    /// A pair line has an asymmetric blank field, an unparseable float, a
    /// non-finite value, or a non-positive/non-finite energy.
    MalformedPair {
        /// Registry name of the offending reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// 1-based source line number.
        line: usize,
        /// The offending line text.
        text: String,
    },
    /// A pair energy is not a group boundary of the anchor row's structure,
    /// breaks contiguity, or the row does not end at the top boundary.
    Misaligned {
        /// Registry name of the offending reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// The offending energy in eV.
        energy: f64,
    },
    /// The anchor row (the wanted row with the most pairs) does not define
    /// exactly [`GROUP_COUNT`] + 1 boundaries, so the file is not in the
    /// SAND-II 725-group structure this pack pins.
    BadAnchor {
        /// Registry name of the would-be anchor reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// Boundaries the anchor actually carries.
        boundaries: usize,
    },
    /// The section's final pair is not the `(top boundary, 0.0)`
    /// terminator the distribution documents.
    BadTerminator {
        /// Registry name of the offending reaction.
        name: &'static str,
        /// ENDF MAT number.
        mat: u32,
        /// ENDF MT number.
        mt: u32,
        /// Terminator energy in eV.
        energy: f64,
        /// Terminator value in barns.
        value: f64,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingSection { name, mat, mt } => {
                write!(f, "IRDFF-II g725: no MF=3 section for {name} (MAT {mat}, MT {mt})")
            }
            Self::DuplicateSection { name, mat, mt } => {
                write!(
                    f,
                    "IRDFF-II g725: duplicate MF=3 section for {name} (MAT {mat}, MT {mt})"
                )
            }
            Self::MalformedControl {
                name,
                mat,
                mt,
                line,
                text,
            } => write!(
                f,
                "IRDFF-II g725: malformed control card for {name} (MAT {mat}, MT {mt}) at line {line}: `{text}`"
            ),
            Self::BadPairCount {
                name,
                mat,
                mt,
                expected,
                found,
            } => write!(
                f,
                "IRDFF-II g725: section {name} (MAT {mat}, MT {mt}) holds {found} pairs, expected {expected}"
            ),
            Self::MalformedPair {
                name,
                mat,
                mt,
                line,
                text,
            } => write!(
                f,
                "IRDFF-II g725: malformed pair line for {name} (MAT {mat}, MT {mt}) at line {line}: `{text}`"
            ),
            Self::Misaligned {
                name,
                mat,
                mt,
                energy,
            } => write!(
                f,
                "IRDFF-II g725: energy {energy:e} eV of {name} (MAT {mat}, MT {mt}) is off the SAND-II 725-group structure"
            ),
            Self::BadAnchor {
                name,
                mat,
                mt,
                boundaries,
            } => write!(
                f,
                "IRDFF-II g725: anchor row {name} (MAT {mat}, MT {mt}) carries {boundaries} boundaries, expected {} (SAND-II 725-group structure)",
                GROUP_COUNT + 1
            ),
            Self::BadTerminator {
                name,
                mat,
                mt,
                energy,
                value,
            } => write!(
                f,
                "IRDFF-II g725: terminator of {name} (MAT {mat}, MT {mt}) is ({energy:e} eV, {value} b), expected (top boundary, 0.0)"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// Number of energy groups of the pinned SAND-II group structure.
pub const GROUP_COUNT: usize = 725;

/// One named reaction of the v1 foil pack: the ENDF section key plus a
/// human-readable reaction title.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IrdffReaction {
    /// Canonical snake-case pack key (`au197_ng`, ...).
    pub name: &'static str,
    /// ENDF MAT number of the parent material.
    pub mat: u32,
    /// ENDF MT number of the reaction.
    pub mt: u32,
    /// Human-readable reaction description.
    pub title: &'static str,
}

impl IrdffReaction {
    /// Const registry constructor.
    pub const fn new(name: &'static str, mat: u32, mt: u32, title: &'static str) -> Self {
        Self {
            name,
            mat,
            mt,
            title,
        }
    }
}

/// The v1 named foil-activation subset (see module docs for the table and
/// provenance). All eight are `MF=3` sections of `IRDFF-II.g725`.
pub const V1_REACTIONS: [IrdffReaction; 8] = [
    IrdffReaction::new("au197_ng", 7925, 102, "197Au(n,g)198Au capture foil"),
    IrdffReaction::new("in115_ng", 4931, 102, "115In(n,g)116In capture foil"),
    IrdffReaction::new("u235_nf", 9228, 18, "235U(n,f) fission chamber"),
    IrdffReaction::new("u238_nf", 9237, 18, "238U(n,f) fast fission chamber"),
    IrdffReaction::new("fe56_np", 2631, 103, "56Fe(n,p)56Mn threshold monitor"),
    IrdffReaction::new("ni58_np", 2825, 103, "58Ni(n,p)58Co threshold monitor"),
    IrdffReaction::new("al27_na", 1325, 107, "27Al(n,a)24Na threshold monitor"),
    IrdffReaction::new(
        "na23_n2n",
        1125,
        16,
        "23Na(n,2n)22Na high-threshold monitor",
    ),
];

/// The v2 named extension set (see module docs for the table and
/// provenance): 26 further `MF=3` dosimetry sections of `IRDFF-II.g725`.
/// The full registry is `V1_REACTIONS` plus these entries.
pub const V2_REACTIONS: [IrdffReaction; 26] = [
    IrdffReaction::new("f19_n2n", 925, 16, "19F(n,2n)18F high-threshold monitor"),
    IrdffReaction::new("b10_na", 525, 107, "10B(n,a)7Li thermal monitor"),
    IrdffReaction::new("mg24_np", 1225, 103, "24Mg(n,p)24Na threshold monitor"),
    IrdffReaction::new("al27_np", 1325, 103, "27Al(n,p)27Mg threshold monitor"),
    IrdffReaction::new("si28_np", 1425, 103, "28Si(n,p)28Al threshold monitor"),
    IrdffReaction::new("p31_np", 1525, 103, "31P(n,p)31Si threshold monitor"),
    IrdffReaction::new("s32_np", 1625, 103, "32S(n,p)32P threshold monitor"),
    IrdffReaction::new("sc45_ng", 2125, 102, "45Sc(n,g)46Sc capture foil"),
    IrdffReaction::new("ti46_np", 2225, 103, "46Ti(n,p)46Sc threshold monitor"),
    IrdffReaction::new("ti47_np", 2228, 103, "47Ti(n,p)47Sc threshold monitor"),
    IrdffReaction::new("ti48_np", 2231, 103, "48Ti(n,p)48Sc threshold monitor"),
    IrdffReaction::new("mn55_n2n", 2525, 16, "55Mn(n,2n)54Mn threshold monitor"),
    IrdffReaction::new("fe54_np", 2625, 103, "54Fe(n,p)54Mn threshold monitor"),
    IrdffReaction::new("co59_ng", 2725, 102, "59Co(n,g)60Co capture foil"),
    IrdffReaction::new("co59_np", 2725, 103, "59Co(n,p)59Fe threshold monitor"),
    IrdffReaction::new(
        "ni58_n2n",
        2825,
        16,
        "58Ni(n,2n)57Ni high-threshold monitor",
    ),
    IrdffReaction::new("cu63_na", 2925, 107, "63Cu(n,a)60Co threshold monitor"),
    IrdffReaction::new("zn64_np", 3025, 103, "64Zn(n,p)64Cu threshold monitor"),
    IrdffReaction::new("in113_ng", 4925, 102, "113In(n,g)114mIn capture foil"),
    IrdffReaction::new("ta181_ng", 7328, 102, "181Ta(n,g)182Ta capture foil"),
    IrdffReaction::new("w186_ng", 7443, 102, "186W(n,g)187W capture foil"),
    IrdffReaction::new("th232_nf", 9040, 18, "232Th(n,f) fast fission chamber"),
    IrdffReaction::new("np237_nf", 9346, 18, "237Np(n,f) fast fission chamber"),
    IrdffReaction::new("pu239_nf", 9437, 18, "239Pu(n,f) fission chamber"),
    IrdffReaction::new(
        "bi209_n2n",
        8325,
        16,
        "209Bi(n,2n)208Bi high-threshold monitor",
    ),
    IrdffReaction::new(
        "bi209_n3n",
        8325,
        17,
        "209Bi(n,3n)207Bi high-threshold monitor",
    ),
];

/// One parsed reaction row: the registry entry plus its group cross
/// sections, expanded to [`GROUP_COUNT`] values (zero below a threshold
/// reaction's first listed group).
#[derive(Debug, Clone, PartialEq)]
pub struct IrdffRow {
    reaction: IrdffReaction,
    sigma: Vec<f64>,
}

impl IrdffRow {
    /// The registry entry (name, MAT, MT, title).
    pub fn reaction(&self) -> IrdffReaction {
        self.reaction
    }

    /// Group-averaged cross sections in barns, one per group, aligned to
    /// [`IrdffPack::bounds`] (group `g` spans `[bounds[g], bounds[g+1])`).
    pub fn sigma(&self) -> &[f64] {
        &self.sigma
    }
}

/// A parsed IRDFF-II g725 pack: the group structure plus one row per
/// requested reaction, in request order.
#[derive(Debug, Clone, PartialEq)]
pub struct IrdffPack {
    bounds: Vec<f64>,
    rows: Vec<IrdffRow>,
}

impl IrdffPack {
    /// Group boundaries in eV, ascending, exactly [`GROUP_COUNT`] + 1
    /// entries (SAND-II 725-group structure: `1e-5` eV to `60` MeV).
    pub fn bounds(&self) -> &[f64] {
        &self.bounds
    }

    /// Parsed rows in request order.
    pub fn rows(&self) -> &[IrdffRow] {
        &self.rows
    }

    /// The row for `name`, or `None` when the name is not in the request.
    pub fn row(&self, name: &str) -> Option<&IrdffRow> {
        self.rows.iter().find(|r| r.reaction.name == name)
    }

    /// Caller-ready response matrix: one cross-section row per requested
    /// reaction (request order), each [`GROUP_COUNT`] values — the shape
    /// the SAND-II iterator in `nucleide-unfold` takes unchanged.
    pub fn response(&self) -> Vec<Vec<f64>> {
        self.rows.iter().map(|r| r.sigma.clone()).collect()
    }
}

/// Parse an ENDF-6 eleven-column float field, including the Fortran
/// embedded-exponent spelling without `E` (`2.589913-5`).
fn endf_f64(field: &str) -> Option<f64> {
    let s = field.trim();
    if s.is_empty() {
        return None;
    }
    if s.contains(['e', 'E']) {
        return s.parse().ok();
    }
    let mut split = None;
    for (i, c) in s.char_indices().skip(1) {
        if c == '+' || c == '-' {
            split = Some(i);
        }
    }
    if let Some(i) = split {
        let (mantissa, exponent) = s.split_at(i);
        // Reconstruct a standard float literal so the parse is correctly
        // rounded (a manual `m * 10^e` product carries double rounding).
        format!("{mantissa}e{exponent}").parse().ok()
    } else {
        s.parse().ok()
    }
}

/// Split a record into its six eleven-column fields (short lines pad with
/// empty fields).
fn fields11(line: &str) -> [&str; 6] {
    let mut out = [""; 6];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = line.get(i * 11..(i + 1) * 11).unwrap_or("");
    }
    out
}

/// Record control fields: (MAT, MF, MT, sequence) from columns 67–80.
fn control(line: &str) -> Option<(u32, u32, u32, u32)> {
    if line.len() < 75 {
        return None;
    }
    let mat = line.get(66..70)?.trim().parse().ok()?;
    let mf = line.get(70..72)?.trim().parse().ok()?;
    let mt = line.get(72..75)?.trim().parse().ok()?;
    let seq = line.get(75..80)?.trim().parse().ok()?;
    Some((mat, mf, mt, seq))
}

/// A raw section freshly read off the file, before structure alignment.
struct RawSection {
    reaction: IrdffReaction,
    /// Pair energies in eV as printed (lower boundaries of consecutive
    /// groups); the last entry is the structural terminator energy.
    energies: Vec<f64>,
    /// Group cross sections in barns (index `k` belongs to the group above
    /// `energies[k]`); the last entry is the terminator value.
    values: Vec<f64>,
}

/// Parse the `MF=3` sections of `wanted` reactions from `IRDFF-II.g725`
/// member text into an [`IrdffPack`].
///
/// The group structure is anchored by the first wanted row carrying
/// [`GROUP_COUNT`] + 1 boundaries (a full-range capture/fission row): its
/// energies must be exactly those boundaries, strictly increasing, and every other
/// wanted row must be a contiguous run of those boundaries ending at the
/// top boundary, with the documented `(top, 0.0)` terminator. Anything
/// else is a loud named error; rows are never silently truncated or
/// padded. An empty `wanted` yields an empty pack.
pub fn parse_g725(text: &str, wanted: &[IrdffReaction]) -> Result<IrdffPack> {
    if wanted.is_empty() {
        return Ok(IrdffPack {
            bounds: Vec::new(),
            rows: Vec::new(),
        });
    }
    let mut seen_wanted: HashSet<(u32, u32)> = HashSet::new();
    for reaction in wanted {
        if !seen_wanted.insert((reaction.mat, reaction.mt)) {
            return Err(Error::DuplicateSection {
                name: reaction.name,
                mat: reaction.mat,
                mt: reaction.mt,
            });
        }
    }
    let by_key: HashMap<(u32, u32), IrdffReaction> =
        wanted.iter().map(|r| ((r.mat, r.mt), *r)).collect();

    let mut found: Vec<RawSection> = Vec::new();
    let mut lines = text.lines().enumerate().peekable();
    while let Some((lineno, line)) = lines.next() {
        let Some((mat, mf, mt, seq)) = control(line) else {
            continue;
        };
        if mf != 3 || mt == 0 || seq != 1 {
            continue;
        }
        let Some(reaction) = by_key.get(&(mat, mt)).copied() else {
            continue;
        };
        if found.iter().any(|s| s.reaction == reaction) {
            return Err(Error::DuplicateSection {
                name: reaction.name,
                mat,
                mt,
            });
        }
        // Head card: ZA and AWR floats, three zero integer fields, and
        // the L2 field (accepted as 0 or the distribution's documented 99).
        let head = fields11(line);
        let head_ok = head[2].trim() == "0"
            && (head[3].trim() == "0" || head[3].trim() == "99")
            && head[4..].iter().all(|field| field.trim() == "0")
            && head[0..2]
                .iter()
                .all(|field| field.trim().parse::<f64>().is_ok_and(|v| v.is_finite()));
        if !head_ok {
            return Err(Error::MalformedControl {
                name: reaction.name,
                mat,
                mt,
                line: lineno + 1,
                text: line.trim().to_string(),
            });
        }
        // ctrl1: two floats, two zero integer fields, a constant, then NG2.
        let Some((ctrl1_line, ctrl1)) = lines.next() else {
            return Err(Error::MalformedControl {
                name: reaction.name,
                mat,
                mt,
                line: lineno + 1,
                text: line.trim().to_string(),
            });
        };
        let c1 = fields11(ctrl1);
        let ng2: usize = match (endf_f64(c1[0]), endf_f64(c1[1]), c1[2].trim(), c1[3].trim()) {
            (Some(_), Some(_), "0", "0") => match c1[5].trim().parse() {
                Ok(n) if n > 0 => n,
                _ => {
                    return Err(Error::MalformedControl {
                        name: reaction.name,
                        mat,
                        mt,
                        line: ctrl1_line + 1,
                        text: ctrl1.trim().to_string(),
                    })
                }
            },
            _ => {
                return Err(Error::MalformedControl {
                    name: reaction.name,
                    mat,
                    mt,
                    line: ctrl1_line + 1,
                    text: ctrl1.trim().to_string(),
                })
            }
        };
        // ctrl2: NG2 repeated, then a constant 1.
        let Some((ctrl2_line, ctrl2)) = lines.next() else {
            return Err(Error::MalformedControl {
                name: reaction.name,
                mat,
                mt,
                line: ctrl1_line + 1,
                text: ctrl1.trim().to_string(),
            });
        };
        let c2 = fields11(ctrl2);
        if c2[0].trim() != c1[5].trim() || c2[1].trim() != "1" {
            return Err(Error::MalformedControl {
                name: reaction.name,
                mat,
                mt,
                line: ctrl2_line + 1,
                text: ctrl2.trim().to_string(),
            });
        }
        // Pair region: NG2 (E, sigma) pairs, three per line, until the
        // record triplet changes or the file ends.
        let mut energies = Vec::with_capacity(ng2);
        let mut values = Vec::with_capacity(ng2);
        while let Some(&(n, next)) = lines.peek() {
            let Some((m2, mf2, mt2, _)) = control(next) else {
                break;
            };
            if m2 != mat || mf2 != 3 || mt2 != mt {
                break;
            }
            let _ = lines.next();
            let f = fields11(next);
            for k in (0..6).step_by(2) {
                let (e_f, v_f) = (endf_f64(f[k]), endf_f64(f[k + 1]));
                match (e_f, v_f) {
                    (Some(e), Some(v)) => {
                        if !e.is_finite() || !v.is_finite() || e <= 0.0 || v < 0.0 {
                            return Err(Error::MalformedPair {
                                name: reaction.name,
                                mat,
                                mt,
                                line: n + 1,
                                text: next.trim().to_string(),
                            });
                        }
                        energies.push(e);
                        values.push(v);
                    }
                    (None, None) => {} // short final line
                    _ => {
                        return Err(Error::MalformedPair {
                            name: reaction.name,
                            mat,
                            mt,
                            line: n + 1,
                            text: next.trim().to_string(),
                        })
                    }
                }
            }
        }
        if energies.len() != ng2 {
            return Err(Error::BadPairCount {
                name: reaction.name,
                mat,
                mt,
                expected: ng2,
                found: energies.len(),
            });
        }
        found.push(RawSection {
            reaction,
            energies,
            values,
        });
    }

    for reaction in wanted {
        if !found.iter().any(|s| s.reaction == *reaction) {
            return Err(Error::MissingSection {
                name: reaction.name,
                mat: reaction.mat,
                mt: reaction.mt,
            });
        }
    }

    // Anchor: the first wanted row carrying the most pairs (a full-range
    // capture/fission row when the file is the pinned structure). A tied
    // row that disagrees with the anchor fails the alignment check below,
    // loudly.
    let max_pairs = found
        .iter()
        .map(|s| s.energies.len())
        .max()
        .expect("wanted non-empty and MissingSection checked above");
    let anchor = wanted
        .iter()
        .find_map(|r| {
            found
                .iter()
                .find(|s| s.reaction == *r && s.energies.len() == max_pairs)
        })
        .expect("wanted non-empty and MissingSection checked above");
    if anchor.energies.len() != GROUP_COUNT + 1 {
        return Err(Error::BadAnchor {
            name: anchor.reaction.name,
            mat: anchor.reaction.mat,
            mt: anchor.reaction.mt,
            boundaries: anchor.energies.len(),
        });
    }
    for w in anchor.energies.windows(2) {
        if w[0] >= w[1] {
            return Err(Error::Misaligned {
                name: anchor.reaction.name,
                mat: anchor.reaction.mat,
                mt: anchor.reaction.mt,
                energy: w[1],
            });
        }
    }
    let position: HashMap<u64, usize> = anchor
        .energies
        .iter()
        .enumerate()
        .map(|(i, &e)| (e.to_bits(), i))
        .collect();

    let mut rows = Vec::with_capacity(found.len());
    for reaction in wanted {
        let section = found
            .iter()
            .find(|s| s.reaction == *reaction)
            .expect("MissingSection checked above");
        let top = anchor.energies[GROUP_COUNT];
        let mut sigma = vec![0.0; GROUP_COUNT];
        let mut previous: Option<usize> = None;
        for (k, &e) in section.energies.iter().enumerate() {
            let Some(&g) = position.get(&e.to_bits()) else {
                return Err(Error::Misaligned {
                    name: section.reaction.name,
                    mat: section.reaction.mat,
                    mt: section.reaction.mt,
                    energy: e,
                });
            };
            if previous.is_some_and(|prev| g != prev + 1) {
                return Err(Error::Misaligned {
                    name: section.reaction.name,
                    mat: section.reaction.mat,
                    mt: section.reaction.mt,
                    energy: e,
                });
            }
            previous = Some(g);
            if k + 1 == section.energies.len() {
                // Terminator: top boundary, value 0.0, no group above it.
                if e != top {
                    return Err(Error::Misaligned {
                        name: section.reaction.name,
                        mat: section.reaction.mat,
                        mt: section.reaction.mt,
                        energy: e,
                    });
                }
                let v = section.values[k];
                if v != 0.0 {
                    return Err(Error::BadTerminator {
                        name: section.reaction.name,
                        mat: section.reaction.mat,
                        mt: section.reaction.mt,
                        energy: e,
                        value: v,
                    });
                }
            } else {
                sigma[g] = section.values[k];
            }
        }
        rows.push(IrdffRow {
            reaction: *reaction,
            sigma,
        });
    }

    Ok(IrdffPack {
        bounds: anchor.energies.clone(),
        rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic registry for tests (mirrors the real one in shape only;
    /// the real registry is exercised end-to-end by the live-fetch gate).
    const AU: IrdffReaction = IrdffReaction::new("au197_ng", 7925, 102, "synthetic Au capture");
    const IN: IrdffReaction = IrdffReaction::new("in115_ng", 4931, 102, "synthetic In capture");
    const FE: IrdffReaction = IrdffReaction::new("fe56_np", 2631, 103, "synthetic Fe threshold");
    const NA: IrdffReaction = IrdffReaction::new("na23_n2n", 1125, 16, "synthetic Na");

    /// Synthetic hand-built group structure with the same shape as the
    /// SAND-II 725-group structure (726 strictly increasing boundaries);
    /// the real boundaries are IAEA data and never appear in tests.
    /// Values are canonicalized through the eleven-column printing so the
    /// parsed pack reproduces them bit-for-bit.
    fn synthetic_bounds() -> Vec<f64> {
        (0..=GROUP_COUNT)
            .map(|g| {
                let raw = 10f64.powf(-5.0 + 11.0 * g as f64 / GROUP_COUNT as f64);
                endf_f64(&fmt_field(raw)).expect("synthetic field must parse")
            })
            .collect()
    }

    /// Synthetic values for a full-range row: `base` + group index, with
    /// the structural terminator slot (index GROUP_COUNT) forced to 0.0.
    fn full_values(base: f64) -> Vec<f64> {
        let mut vals: Vec<f64> = (0..GROUP_COUNT + 1).map(|g| base + g as f64).collect();
        vals[GROUP_COUNT] = 0.0;
        vals
    }

    /// Format one field in the distribution's eleven-column style,
    /// including the embedded-exponent spelling for sub-unit magnitudes.
    fn fmt_field(x: f64) -> String {
        let s = if x == 0.0 {
            "0.0".to_string()
        } else if x.abs() < 1.0 {
            let raw = format!("{x:.6e}");
            let (m, e) = raw.split_once('e').unwrap();
            format!(
                "{}{}{}",
                m,
                if e.starts_with('-') { "-" } else { "+" },
                e[1..].trim_start_matches('0')
            )
        } else if x.abs() < 1.0e5 {
            format!("{x:.5}")
        } else {
            format!("{x:.1}")
        };
        format!("{s:>11}")
    }

    /// Build one synthetic MF=3 section in the exact distributed layout:
    /// head card, ctrl1, ctrl2, then `energies.len()` (E, sigma) pairs at
    /// three per line (the final line padded with blanks when short).
    fn section_text(za: f64, mat: u32, mt: u32, energies: &[f64], values: &[f64]) -> String {
        assert_eq!(energies.len(), values.len());
        let zero = format!("{:>11}", 0);
        let mut out = String::new();
        // Head card, then the two control cards, in the distributed layout.
        out.push_str(&format!(
            "{}{}{zero}{zero}{zero}{zero}{mat:>4} 3{mt:>3}    1\n",
            fmt_field(za),
            fmt_field(100.0),
        ));
        let ctrl1 = format!(
            "{}{}{zero}{zero}{:>11}{:>11}",
            fmt_field(0.0),
            fmt_field(0.0),
            1,
            energies.len()
        );
        out.push_str(&format!("{ctrl1}{mat:>4} 3{mt:>3}    2\n"));
        let ctrl2 = format!("{:>11}{:>11}", energies.len(), 1);
        out.push_str(&format!("{ctrl2}{mat:>4} 3{mt:>3}    3\n"));
        for (chunk_idx, chunk) in energies.chunks(3).enumerate() {
            let mut fields = String::new();
            for (i, &e) in chunk.iter().enumerate() {
                fields.push_str(&fmt_field(e));
                fields.push_str(&fmt_field(values[chunk_idx * 3 + i]));
            }
            for _ in chunk.len()..3 {
                fields.push_str(&" ".repeat(22));
            }
            out.push_str(&format!(
                "{fields}{mat:>4} 3{mt:>3}{seq:>5}\n",
                seq = chunk_idx + 4
            ));
        }
        out
    }

    /// Two full-range rows (Au, In) over the synthetic structure plus one
    /// threshold row (Fe) starting at group 500 — the real pack's shape.
    fn synthetic_pack_text() -> String {
        let bounds = synthetic_bounds();
        let fe_energies = &bounds[500..];
        let mut fe_vals: Vec<f64> = (0..fe_energies.len()).map(|k| 300.0 + k as f64).collect();
        fe_vals[fe_energies.len() - 1] = 0.0;
        let mut text = String::new();
        text.push_str(&section_text(
            79197.0,
            AU.mat,
            AU.mt,
            &bounds,
            &full_values(100.0),
        ));
        text.push_str(&section_text(
            49115.0,
            IN.mat,
            IN.mt,
            &bounds,
            &full_values(200.0),
        ));
        text.push_str(&section_text(26056.0, FE.mat, FE.mt, fe_energies, &fe_vals));
        text
    }

    #[test]
    fn parses_full_range_and_threshold_rows() {
        let bounds = synthetic_bounds();
        let pack = parse_g725(&synthetic_pack_text(), &[AU, IN, FE]).unwrap();
        assert_eq!(pack.bounds(), bounds.as_slice());
        assert_eq!(pack.rows().len(), 3);

        // Full-range Au row: terminator dropped, GROUP_COUNT entries.
        let au = pack.row("au197_ng").unwrap();
        assert_eq!(au.sigma().len(), GROUP_COUNT);
        assert_eq!(au.sigma()[0], 100.0);
        assert_eq!(au.sigma()[GROUP_COUNT - 1], 824.0);
        assert_eq!(au.reaction().mat, 7925);
        assert_eq!(au.reaction().mt, 102);

        // Second full-range row aligns to the same structure.
        let in_row = pack.row("in115_ng").unwrap();
        assert_eq!(in_row.sigma()[10], 210.0);

        // Threshold row: zeros below group 500, then contiguous values.
        let fe = pack.row("fe56_np").unwrap();
        assert!(fe.sigma()[..500].iter().all(|&v| v == 0.0));
        assert_eq!(fe.sigma()[500], 300.0);
        assert_eq!(fe.sigma()[GROUP_COUNT - 1], 300.0 + 224.0);

        // Caller-ready response matrix in request order.
        let response = pack.response();
        assert_eq!(response.len(), 3);
        assert!(response.iter().all(|r| r.len() == GROUP_COUNT));
        assert_eq!(response[2][500], 300.0);
    }

    #[test]
    fn request_order_is_preserved() {
        let pack = parse_g725(&synthetic_pack_text(), &[FE, AU]).unwrap();
        assert_eq!(pack.rows()[0].reaction().name, "fe56_np");
        assert_eq!(pack.rows()[1].reaction().name, "au197_ng");
    }

    #[test]
    fn missing_section_is_named() {
        match parse_g725(&synthetic_pack_text(), &[AU, NA]).unwrap_err() {
            Error::MissingSection { name, mat, mt } => {
                assert_eq!(name, "na23_n2n");
                assert_eq!((mat, mt), (1125, 16));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn duplicate_request_and_duplicate_section_error() {
        // Duplicate in the request itself.
        match parse_g725("", &[AU, AU]).unwrap_err() {
            Error::DuplicateSection { name, .. } => assert_eq!(name, "au197_ng"),
            other => panic!("{other:?}"),
        }
        // Duplicate section in the file.
        let dup = format!("{}{}", synthetic_pack_text(), synthetic_pack_text());
        match parse_g725(&dup, &[AU]).unwrap_err() {
            Error::DuplicateSection { name, mat, mt } => {
                assert_eq!((name, mat, mt), ("au197_ng", 7925, 102));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn malformed_control_cards_are_loud() {
        let text = synthetic_pack_text();
        // ctrl2's NG2 repeat disagrees with ctrl1.
        let broken = text.replacen("        726          1", "        725          1", 1);
        match parse_g725(&broken, &[AU]).unwrap_err() {
            Error::MalformedControl { name, line, .. } => {
                assert_eq!(name, "au197_ng");
                assert_eq!(line, 3);
            }
            other => panic!("{other:?}"),
        }
        // Non-zero integer field in the head card (first of the four).
        let zero_field = format!("{:>11}", 0);
        let head_prefix = format!("{}{}{zero_field}", fmt_field(79197.0), fmt_field(100.0));
        let broken_head = head_prefix.replacen(&zero_field, &format!("{:>11}", 7), 1);
        let broken = text.replacen(&head_prefix, &broken_head, 1);
        match parse_g725(&broken, &[AU]).unwrap_err() {
            Error::MalformedControl { name, line, .. } => {
                assert_eq!(name, "au197_ng");
                assert_eq!(line, 1);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn head_l2_quirk_is_accepted_but_still_pinned() {
        // The pinned distribution carries L2=99 (not 0) in three MF=3 head
        // cards, notably v1's 56Fe(n,p): the parser accepts 0 or 99 there
        // and rejects anything else loudly.
        let bounds = synthetic_bounds();
        let mut fe_vals = vec![3.0; bounds.len() - 500];
        *fe_vals.last_mut().unwrap() = 0.0;
        let base = format!(
            "{}{}",
            section_text(79197.0, AU.mat, AU.mt, &bounds, &full_values(100.0)),
            section_text(26056.0, FE.mat, FE.mt, &bounds[500..], &fe_vals),
        );
        let marker = format!("{:>4} 3{:>3}{:>5}", FE.mat, FE.mt, 1);
        let head_idx = base
            .lines()
            .position(|l| l.ends_with(&marker))
            .expect("synthetic Fe head card");
        let with_l2 = |l2: u32| {
            base.lines()
                .enumerate()
                .map(|(i, l)| {
                    if i == head_idx {
                        format!("{}{:>11}{}", &l[..33], l2, &l[44..])
                    } else {
                        l.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join("\n")
                + "\n"
        };
        let pack = parse_g725(&with_l2(99), &[AU, FE]).unwrap();
        let fe = pack.row("fe56_np").unwrap();
        assert!(fe.sigma()[..500].iter().all(|&v| v == 0.0));
        assert_eq!(fe.sigma()[500], 3.0);
        match parse_g725(&with_l2(98), &[AU, FE]).unwrap_err() {
            Error::MalformedControl { name, line, .. } => {
                assert_eq!(name, "fe56_np");
                assert_eq!(line, head_idx + 1);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn pair_count_mismatch_is_loud() {
        let bounds = synthetic_bounds();
        // Rebuild the Au section without its final data line: NG2 still
        // declares 726 pairs but only 723 are present.
        let mut full = String::new();
        for line in section_text(79197.0, AU.mat, AU.mt, &bounds, &full_values(1.0))
            .lines()
            .take(3 + 241)
        {
            full.push_str(line);
            full.push('\n');
        }
        full.push_str(&section_text(
            49115.0,
            IN.mat,
            IN.mt,
            &bounds,
            &full_values(2.0),
        ));
        match parse_g725(&full, &[AU]).unwrap_err() {
            Error::BadPairCount {
                name,
                expected,
                found,
                ..
            } => {
                assert_eq!(name, "au197_ng");
                assert_eq!((expected, found), (726, 723));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn malformed_pair_lines_are_loud() {
        let text = synthetic_pack_text();
        // Unparseable value field (the group-1 cross section, unique to
        // the pair region).
        let bad = text.replacen(
            &fmt_field(101.0),
            &fmt_field(101.0).replacen("101", "not", 1),
            1,
        );
        assert!(matches!(
            parse_g725(&bad, &[AU]).unwrap_err(),
            Error::MalformedPair { .. }
        ));
        // Non-finite value.
        let bad = text.replacen(&fmt_field(101.0), "    1.0+999", 1);
        assert!(matches!(
            parse_g725(&bad, &[AU]).unwrap_err(),
            Error::MalformedPair { .. }
        ));
        // Negative cross section: physics forbids it (σ ≥ 0), so the
        // parser rejects it here rather than downstream in `unfold`.
        let bad = text.replacen(&fmt_field(101.0), &fmt_field(-1.0), 1);
        assert!(matches!(
            parse_g725(&bad, &[AU]).unwrap_err(),
            Error::MalformedPair { .. }
        ));
    }

    #[test]
    fn misaligned_energies_are_loud() {
        let bounds = synthetic_bounds();
        // Shift one interior energy of the *In* row off the Au anchor's
        // structure while keeping the row ordered (bounds are ~3.5% apart,
        // so x1.02 stays in the gap): the parser must reject the In row as
        // off-structure at the shifted value. (A shifted energy inside the
        // anchor row itself would merely redefine the structure, which is
        // why the anchor is always a full-range row from the pinned file.)
        let shifted = endf_f64(&fmt_field(bounds[1] * 1.02)).unwrap();
        let mut in_energies = bounds.clone();
        in_energies[1] = shifted;
        let text = format!(
            "{}{}",
            section_text(79197.0, AU.mat, AU.mt, &bounds, &full_values(100.0)),
            section_text(49115.0, IN.mat, IN.mt, &in_energies, &full_values(200.0))
        );
        match parse_g725(&text, &[AU, IN]).unwrap_err() {
            Error::Misaligned { name, energy, .. } => {
                assert_eq!(name, "in115_ng");
                assert_eq!(energy, shifted);
            }
            other => panic!("{other:?}"),
        }
        // Threshold row not ending at the top boundary.
        let mut full = section_text(79197.0, AU.mat, AU.mt, &bounds, &full_values(1.0));
        full.push_str(&section_text(
            26056.0,
            FE.mat,
            FE.mt,
            &bounds[500..724],
            &vec![3.0; 224],
        ));
        match parse_g725(&full, &[AU, FE]).unwrap_err() {
            Error::Misaligned { name, .. } => assert_eq!(name, "fe56_np"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn bad_anchor_and_bad_terminator_are_loud() {
        // No full-range row: the largest row carries 226 boundaries.
        let bounds = synthetic_bounds();
        let fe_only = section_text(26056.0, FE.mat, FE.mt, &bounds[500..], &{
            let mut v = vec![3.0; 226];
            v[225] = 0.0;
            v
        });
        match parse_g725(&fe_only, &[FE]).unwrap_err() {
            Error::BadAnchor {
                name, boundaries, ..
            } => {
                assert_eq!(name, "fe56_np");
                assert_eq!(boundaries, 226);
            }
            other => panic!("{other:?}"),
        }
        // Terminator value nonzero.
        let mut vals = full_values(1.0);
        vals[GROUP_COUNT] = 7.0;
        let bad_term = section_text(79197.0, AU.mat, AU.mt, &bounds, &vals);
        match parse_g725(&bad_term, &[AU]).unwrap_err() {
            Error::BadTerminator { name, value, .. } => {
                assert_eq!(name, "au197_ng");
                assert_eq!(value, 7.0);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn empty_request_yields_empty_pack() {
        let pack = parse_g725("anything", &[]).unwrap();
        assert!(pack.bounds().is_empty());
        assert!(pack.rows().is_empty());
        assert!(pack.response().is_empty());
    }

    #[test]
    fn embedded_exponent_fields_parse() {
        assert_eq!(endf_f64(" 2.589913-5"), Some(2.589913e-5));
        assert_eq!(endf_f64(" 4.632640-6"), Some(4.632640e-6));
        assert_eq!(endf_f64("-2913100.00"), Some(-2913100.00));
        assert_eq!(endf_f64(" 1.00000E-5"), Some(1.0e-5));
        assert_eq!(endf_f64(" 60000000.0"), Some(6.0e7));
        assert_eq!(endf_f64(" .001020819"), Some(0.001020819));
        assert_eq!(endf_f64("           "), None);
        assert_eq!(endf_f64("  not-a-num"), None);
    }

    #[test]
    fn v1_registry_is_well_formed() {
        let mut keys = HashSet::new();
        for r in V1_REACTIONS {
            assert!(keys.insert(r.name));
            assert!(r.title.contains("(n,"));
            assert!(r.mat > 0 && r.mt > 0);
        }
        assert_eq!(V1_REACTIONS.len(), 8);
    }

    #[test]
    fn v2_registry_is_well_formed_and_disjoint_from_v1() {
        let mut keys = HashSet::new();
        for r in V1_REACTIONS {
            keys.insert((r.name, r.mat, r.mt));
        }
        for r in V2_REACTIONS {
            assert!(r.title.contains("(n,"));
            assert!(r.mat > 0 && r.mt > 0);
            assert!(
                keys.insert((r.name, r.mat, r.mt)),
                "V2 entry duplicates the registry: {}",
                r.name
            );
        }
        assert_eq!(V2_REACTIONS.len(), 26);
    }
}
