//! SSW surface-source ↔ MCPL particle-list conversion (neutron/gamma-only v1).
//!
//! Thin conversion layer over [`crate::Header`]/[`crate::Particle`] and the
//! `mcnp-io` SSW header/track types. All business logic lives here; the
//! Python facades only move dicts across the boundary.
//!
//! ## Field mapping
//!
//! | SSW track field | MCPL particle field | Rule |
//! | --- | --- | --- |
//! | `erg` (MeV) | `ekin` (MeV) | verbatim both directions |
//! | `tme` (shakes) | `time` (ms) | shakes → ms `× 1e-5` ([`SHAKES_TO_MS`]), ms → shakes `× 1e5` (`MS_TO_SHAKES`) |
//! | surface id | `userflags` | SSW → MCPL only when [`Ssw2McplOptions::surf_to_userflags`] is set (default); MCPL → SSW reads `userflags` unless [`Mcpl2SswOptions::surface`] overrides |
//! | particle kind | `pdgcode` | neutron ↔ 2112, gamma ↔ 22; anything else is [`SswError::UnsupportedPdg`] |
//! | `wgt`, `x`/`y`/`z`, `u`/`v`/`cs` | `weight`, `position`, `direction` | verbatim |
//!
//! ## Explicit caller parameters instead of guessed decodes
//!
//! The only per-track packed-field layout derivable from the local `mcnp-io`
//! code is [`TrackData::cell`](nucleide_mcnp_io::surfsrc::TrackData::cell)
//! (`abs(bitarray) // 8 % 100000000`) plus the
//! `w = sqrt(1-u²-v²)·sign(bitarray)` derivation in its track builder. No
//! surface-id (`isurf`) or particle-type (`rawtype`) decode exists locally,
//! so both ride in as explicit caller parameters on every [`SswTrack`]
//! (`surf`, `kind`) instead of being guessed out of `bitarray`.
//!
//! MCPL → SSW direction consequences, all documented rather than inferred:
//!
//! - The 11-double SSW record carries no per-track type slot (particle type
//!   is file-level `mipts` in the header, passed through opaquely from the
//!   reference), so `kind` only gates the conversion: non-neutron/gamma PDG
//!   codes are rejected, never mistyped.
//! - `bitarray` magnitude is 0 (the track's cell and upstream type word are
//!   unknown to the converter); its sign mirrors the direction z so the
//!   reader's `copysign`-derived `w` matches `cs`. Third-party tools that
//!   decode the particle type out of `bitarray` will skip these tracks —
//!   use the nucleide leg for round-trips (named-open, see below).
//! - `nps` is the 1-based particle index (MCPL carries no history number).
//! - No direction cosine is forced: the local SSW writer stores track
//!   records verbatim, so `(u, v, cs)` propagate the particle direction
//!   unchanged. Upstream `mcpl2ssw` instead forces the stored `cs` slot to
//!   1.0 (verified against the 2.2.8 binary); v1 keeps the true cosine and
//!   documents the difference rather than destroying information.
//! - The output header clones the reference (code/version/deck passthrough)
//!   with the count fields patched to the particle count — `nrss`, `np1`,
//!   and `orignp1` with the reference's table-2 sign preserved — exactly as
//!   upstream `mcpl2ssw` reports it ("N particles (nrss) and N histories
//!   (np1)", verified against the 2.2.8 binary). `niss` passes through
//!   (upstream behavior there is unverified).
//!
//! ## Named-open items (explicitly out of v1 scope)
//!
//! - Particle types beyond neutron/gamma ([`SswError::UnsupportedPdg`]).
//! - Polarisation (always zeros; MCPL output never sets `has_polarisation`)
//!   and universal PDG/weight codes (always per-particle).
//! - Transport semantics; this only re-homes already-transported tracks.
//! - Reference headers whose `abs(ncrd)` is not
//!   [`TrackData::RECORD_WIDTH`](nucleide_mcnp_io::surfsrc::TrackData::RECORD_WIDTH)
//!   ([`SswError::UnsupportedRecordWidth`]).
//! - Upstream type-word round-trips: v1 never writes cell/type words into
//!   `bitarray` (writing them would transcribe the upstream rawtype table
//!   into this crate, which the license-clean v1 cut forbids), so
//!   third-party `bitarray`-decoding tools skip v1-written tracks.
//!   Owner/next step: revisit when a licensed type-table source exists.

use nucleide_mcnp_io::surfsrc::{SurfSrcHeader, TrackData};
use thiserror::Error;

use crate::{encode_file, Blob, Header, Particle};

/// Shakes → milliseconds (1 shake = 1e-8 s = 1e-5 ms).
pub const SHAKES_TO_MS: f64 = 1e-5;
/// Milliseconds → shakes.
pub const MS_TO_SHAKES: f64 = 1e5;
/// PDG code for the neutron (only non-photon kind v1 converts).
pub const PDG_NEUTRON: i32 = 2112;
/// PDG code for the gamma (only non-neutron kind v1 converts).
pub const PDG_GAMMA: i32 = 22;
/// Largest deck blob accepted by [`Ssw2McplOptions::deck_blob`] (100 MiB).
pub const MAX_DECK_BLOB_BYTES: usize = 100 * 1024 * 1024;
/// Largest surface id accepted on the MCPL → SSW path.
pub const MAX_SURFACE_ID: u32 = 999_999;

/// SSW particle kind (neutron/gamma-only v1).
///
/// The kind is an explicit caller parameter on [`SswTrack`]: no particle-type
/// decode exists in the local SSW reader, so v1 never guesses it from
/// `bitarray` (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SswParticleKind {
    /// Neutron (PDG 2112).
    Neutron,
    /// Gamma (PDG 22).
    Gamma,
}

impl SswParticleKind {
    /// PDG code for this kind (neutron 2112, gamma 22).
    pub fn pdg(self) -> i32 {
        match self {
            SswParticleKind::Neutron => PDG_NEUTRON,
            SswParticleKind::Gamma => PDG_GAMMA,
        }
    }

    /// Map a PDG code back to a kind; anything else is [`SswError::UnsupportedPdg`].
    pub fn from_pdg(pdg: i32) -> Result<Self, SswError> {
        match pdg {
            PDG_NEUTRON => Ok(SswParticleKind::Neutron),
            PDG_GAMMA => Ok(SswParticleKind::Gamma),
            other => Err(SswError::UnsupportedPdg(other)),
        }
    }

    /// `"neutron"` or `"gamma"` (the only spellings the Python facade accepts).
    pub fn as_str(self) -> &'static str {
        match self {
            SswParticleKind::Neutron => "neutron",
            SswParticleKind::Gamma => "gamma",
        }
    }

    /// Parse `"neutron"`/`"gamma"`; anything else is `None` (the facade turns
    /// this into a `ValueError` naming the accepted values).
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "neutron" => Some(SswParticleKind::Neutron),
            "gamma" => Some(SswParticleKind::Gamma),
            _ => None,
        }
    }
}

/// One SSW-side track with its explicit conversion parameters.
///
/// `ekin` is in MeV, `time_shakes` in shakes, `position` in cm; `direction`
/// must be a unit vector. `surf`/`kind` are caller-supplied (never decoded
/// from `bitarray`; see the module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct SswTrack {
    /// Kinetic energy in MeV (maps verbatim to MCPL `ekin`).
    pub ekin: f64,
    /// Track time in shakes (maps to MCPL `time` in ms via [`SHAKES_TO_MS`]).
    pub time_shakes: f64,
    /// Position in cm (maps verbatim to MCPL `position`).
    pub position: [f64; 3],
    /// Unit direction (maps verbatim to MCPL `direction`).
    pub direction: [f64; 3],
    /// Statistical weight (maps verbatim to MCPL `weight`).
    pub weight: f64,
    /// Surface id (stored as MCPL `userflags` when
    /// [`Ssw2McplOptions::surf_to_userflags`] is set).
    pub surf: u32,
    /// Particle kind (maps to the MCPL `pdgcode`).
    pub kind: SswParticleKind,
}

/// Optional MCPL header blob embedding the SSW deck text.
#[derive(Debug, Clone, PartialEq)]
pub struct DeckBlob {
    /// Blob key string (caller-chosen, e.g. `"ssw_deck"`).
    pub key: String,
    /// Raw deck bytes (must fit [`MAX_DECK_BLOB_BYTES`]).
    pub data: Vec<u8>,
}

/// Options for [`ssw2mcpl`].
#[derive(Debug, Clone, PartialEq)]
pub struct Ssw2McplOptions {
    /// Emit double precision (`false` = single-prec default).
    pub double_prec: bool,
    /// Store each track's `surf` as the particle `userflags` (default `true`;
    /// also drives the output `has_userflags` flag).
    pub surf_to_userflags: bool,
    /// Gzip-compress the [`ssw2mcpl_bytes`] output (pure [`ssw2mcpl`]
    /// particles are unaffected; transport-only).
    pub gzip: bool,
    /// Optional deck blob embedded in the MCPL header.
    pub deck_blob: Option<DeckBlob>,
    /// MCPL `srcname` string.
    pub srcname: String,
    /// MCPL header comments (round-tripped verbatim, never interpreted).
    pub comments: Vec<String>,
}

impl Default for Ssw2McplOptions {
    fn default() -> Self {
        Ssw2McplOptions {
            double_prec: false,
            surf_to_userflags: true,
            gzip: false,
            deck_blob: None,
            srcname: "ssw2mcpl".to_string(),
            comments: Vec::new(),
        }
    }
}

/// Options for [`mcpl2ssw`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mcpl2SswOptions {
    /// Surface id stamped on every output track (`None` = each particle's
    /// `userflags`; must lie in `[1, 999999]` either way).
    pub surface: Option<u32>,
}

/// Errors raised by SSW ↔ MCPL conversion.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum SswError {
    /// MCPL encode failure inside [`ssw2mcpl_bytes`].
    #[error(transparent)]
    Mcpl(#[from] crate::Error),
    /// Filesystem or gzip failure while writing conversion output.
    #[error("io error: {0}")]
    Io(String),
    /// Particle kind outside the neutron/gamma-only v1 scope.
    #[error("unsupported PDG code {0} (v1 converts neutrons/2112 and gammas/22 only)")]
    UnsupportedPdg(i32),
    /// Surface id outside `[1, 999999]`.
    #[error("surface id {0} outside [1, 999999]")]
    SurfaceIdOutOfRange(u32),
    /// `userflags` is 0 while no [`Mcpl2SswOptions::surface`] override is set.
    #[error("particle {index} has userflags 0 with no surface override set")]
    MissingSurfaceId {
        /// Position of the offending particle in the input slice.
        index: usize,
    },
    /// Deck blob larger than [`MAX_DECK_BLOB_BYTES`].
    #[error("deck blob holds {bytes} bytes, above the {limit}-byte cap")]
    DeckBlobTooLarge {
        /// Offered blob length in bytes.
        bytes: usize,
        /// Cap in bytes ([`MAX_DECK_BLOB_BYTES`]).
        limit: usize,
    },
    /// Reference header record width is not the 11-double track layout.
    #[error("reference ncrd width {found} != {expected} (v1 converts 11-double tracks only)")]
    UnsupportedRecordWidth {
        /// Required width ([`TrackData::RECORD_WIDTH`](nucleide_mcnp_io::surfsrc::TrackData::RECORD_WIDTH)).
        expected: usize,
        /// Reference header's `abs(ncrd)`.
        found: usize,
    },
    /// Input direction deviates from unit length beyond [`crate::UNIT_TOL`].
    #[error("track {0} direction is not a unit vector (|d|^2 = {1})")]
    NonUnitDirection(usize, f64),
    /// Input kinetic energy is negative.
    #[error("track {0} has negative kinetic energy {1}")]
    NegativeEnergy(usize, f64),
}

/// Convert SSW-side tracks to an MCPL header plus particle records.
///
/// Energy maps verbatim (MeV), time maps shakes → ms (`× 1e-5`), surface ids
/// map to `userflags` when [`Ssw2McplOptions::surf_to_userflags`] is set.
/// Direction/energy are validated eagerly (unit [`crate::UNIT_TOL`],
/// non-negative) so the pure function fails exactly where the MCPL writer
/// would.
pub fn ssw2mcpl(
    tracks: &[SswTrack],
    options: &Ssw2McplOptions,
) -> Result<(Header, Vec<Particle>), SswError> {
    if let Some(blob) = &options.deck_blob {
        if blob.data.len() > MAX_DECK_BLOB_BYTES {
            return Err(SswError::DeckBlobTooLarge {
                bytes: blob.data.len(),
                limit: MAX_DECK_BLOB_BYTES,
            });
        }
    }
    let mut particles = Vec::with_capacity(tracks.len());
    for (i, t) in tracks.iter().enumerate() {
        let dir2 = t.direction[0] * t.direction[0]
            + t.direction[1] * t.direction[1]
            + t.direction[2] * t.direction[2];
        if (dir2 - 1.0).abs() > crate::UNIT_TOL {
            return Err(SswError::NonUnitDirection(i, dir2));
        }
        if t.ekin < 0.0 {
            return Err(SswError::NegativeEnergy(i, t.ekin));
        }
        particles.push(Particle {
            ekin: t.ekin,
            polarisation: [0.0; 3],
            position: t.position,
            direction: t.direction,
            time: t.time_shakes * SHAKES_TO_MS,
            weight: t.weight,
            pdgcode: t.kind.pdg(),
            userflags: if options.surf_to_userflags { t.surf } else { 0 },
        });
    }
    let header = Header {
        has_userflags: options.surf_to_userflags,
        double_prec: options.double_prec,
        srcname: options.srcname.clone(),
        comments: options.comments.clone(),
        blobs: options
            .deck_blob
            .as_ref()
            .map(|b| Blob {
                key: b.key.clone(),
                data: b.data.clone(),
            })
            .into_iter()
            .collect(),
        nparticles: particles.len() as u64,
        ..Header::default()
    };
    Ok((header, particles))
}

/// Encode [`ssw2mcpl`] output to MCPL file bytes.
///
/// Uses the existing MCPL writer (always format 3); when
/// [`Ssw2McplOptions::gzip`] is set the bytes are gzip-compressed with the
/// same encoder settings as [`crate::write_to_path`] (compressed bytes stay
/// encoder-dependent and are asserted by magic only, never byte-exact).
pub fn ssw2mcpl_bytes(tracks: &[SswTrack], options: &Ssw2McplOptions) -> Result<Vec<u8>, SswError> {
    let (header, particles) = ssw2mcpl(tracks, options)?;
    let bytes = encode_file(&header, &particles)?;
    if !options.gzip {
        return Ok(bytes);
    }
    use std::io::Write as _;
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(&bytes)
        .map_err(|e| SswError::Io(e.to_string()))?;
    enc.finish().map_err(|e| SswError::Io(e.to_string()))
}

/// Resolve the output surface id for one particle: the
/// [`Mcpl2SswOptions::surface`] override when set, else the particle's
/// `userflags`. Both paths require `[1, 999999]`; `userflags == 0` without
/// an override is [`SswError::MissingSurfaceId`].
fn resolve_surface(index: usize, p: &Particle, options: &Mcpl2SswOptions) -> Result<u32, SswError> {
    let surf = options.surface.unwrap_or(p.userflags);
    if options.surface.is_none() && p.userflags == 0 {
        return Err(SswError::MissingSurfaceId { index });
    }
    if !(1..=MAX_SURFACE_ID).contains(&surf) {
        return Err(SswError::SurfaceIdOutOfRange(surf));
    }
    Ok(surf)
}

/// Convert MCPL particles to an SSW header plus track records against a
/// reference SSW header.
///
/// The output header clones `reference` (code/version/deck passthrough) with
/// the count fields patched to the particle count: `nrss` and `np1` become
/// the particle count and `orignp1` keeps the reference's table-2 sign.
/// Surface ids come from
/// `userflags` or [`Mcpl2SswOptions::surface`]; PDG codes outside 2112/22
/// are [`SswError::UnsupportedPdg`]. Energy maps verbatim (MeV), time maps
/// ms → shakes (`× 1e5`), `(u, v, cs)` carry the particle direction
/// unchanged (no cosine is forced; see the module docs).
pub fn mcpl2ssw(
    particles: &[Particle],
    reference: &SurfSrcHeader,
    options: &Mcpl2SswOptions,
) -> Result<(SurfSrcHeader, Vec<TrackData>), SswError> {
    let width = reference.ncrd.unsigned_abs() as usize;
    if width != TrackData::RECORD_WIDTH {
        return Err(SswError::UnsupportedRecordWidth {
            expected: TrackData::RECORD_WIDTH,
            found: width,
        });
    }
    let mut tracks = Vec::with_capacity(particles.len());
    for (i, p) in particles.iter().enumerate() {
        // Gate the kind even though the 11-double record has no per-track
        // type slot: mistyping a proton as a neutron must be loud, not silent.
        let _kind = SswParticleKind::from_pdg(p.pdgcode)?;
        let _surf = resolve_surface(i, p, options)?;
        let dir2 = p.direction[0] * p.direction[0]
            + p.direction[1] * p.direction[1]
            + p.direction[2] * p.direction[2];
        if (dir2 - 1.0).abs() > crate::UNIT_TOL {
            return Err(SswError::NonUnitDirection(i, dir2));
        }
        if p.ekin < 0.0 {
            return Err(SswError::NegativeEnergy(i, p.ekin));
        }
        // Cell unknown to the converter (magnitude 0); the sign mirrors the
        // direction z so the reader's copysign-derived `w` matches `cs`.
        let bitarray = if p.direction[2].is_sign_negative() {
            -0.0
        } else {
            0.0
        };
        let record = vec![
            (i + 1) as f64, // nps: 1-based particle index
            bitarray,
            p.weight,
            p.ekin,
            p.time * MS_TO_SHAKES,
            p.position[0],
            p.position[1],
            p.position[2],
            p.direction[0],
            p.direction[1],
            p.direction[2],
        ];
        tracks.push(TrackData::from_record(record));
    }
    let mut header = reference.clone();
    let n = particles.len() as i64;
    header.nrss = n;
    // Upstream patches the history counts to the converted tally (verified:
    // "N particles (nrss) and N histories (np1)"); the table-2 sign rides
    // along so the cloned deck layout stays writable.
    header.np1 = n;
    header.orignp1 = if reference.orignp1 < 0 { -n } else { n };
    Ok((header, tracks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nucleide_mcnp_io::surfsrc::{SourceSurf, SurfSrc};

    fn synth_tracks() -> Vec<SswTrack> {
        vec![
            SswTrack {
                ekin: 2.5,
                time_shakes: 3.0e5, // == 3.0 ms
                position: [1.0, -2.0, 0.5],
                direction: [0.0, 0.0, 1.0],
                weight: 1.0,
                surf: 100,
                kind: SswParticleKind::Neutron,
            },
            SswTrack {
                ekin: 0.662,
                time_shakes: 0.0,
                position: [0.0, 0.0, 0.0],
                direction: [1.0, 0.0, 0.0],
                weight: 0.5,
                surf: 200,
                kind: SswParticleKind::Gamma,
            },
        ]
    }

    fn synth_reference() -> SurfSrcHeader {
        SurfSrcHeader {
            kod: "mcnp    ".to_string(),
            ver: "5    ".to_string(),
            loddat: "01012026".to_string(),
            idtm: " 09/12/26 00:00:00 ".to_string(),
            probid: " 09/12/26 00:00:00 ".to_string(),
            aid: "synthetic reference".to_string(),
            knod: 2,
            np1: 1000,
            orignp1: -1000,
            nrss: 2,
            ncrd: 11,
            njsw: 1,
            niss: 2,
            niwr: Some(0),
            mipts: Some(3),
            kjaq: Some(0),
            table1extra: Vec::new(),
            table2extra: Vec::new(),
            surflist: vec![SourceSurf {
                id: 100,
                facet_id: -1,
                surface_type: 1,
                num_params: 1,
                surf_params: vec![0.0],
            }],
            summary_table: vec![0; 15],
            summary_extra: Vec::new(),
        }
    }

    #[test]
    fn ssw2mcpl_closed_form_defaults() {
        let (h, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        assert_eq!(ps.len(), 2);
        assert_eq!(h.nparticles, 2);
        assert!(!h.double_prec);
        assert!(h.has_userflags);
        assert_eq!(h.srcname, "ssw2mcpl");
        assert!(h.blobs.is_empty());
        // Neutron: energy verbatim, shakes→ms within one f64 rounding
        // (3.0e5 * 1e-5 == 3.0000000000000004, not 3.0).
        assert_eq!(ps[0].ekin, 2.5);
        assert!((ps[0].time - 3.0).abs() < 1e-12);
        assert_eq!(ps[0].pdgcode, PDG_NEUTRON);
        assert_eq!(ps[0].userflags, 100);
        assert_eq!(ps[0].position, [1.0, -2.0, 0.5]);
        assert_eq!(ps[0].direction, [0.0, 0.0, 1.0]);
        assert_eq!(ps[0].weight, 1.0);
        // Gamma.
        assert_eq!(ps[1].ekin, 0.662);
        assert_eq!(ps[1].time, 0.0);
        assert_eq!(ps[1].pdgcode, PDG_GAMMA);
        assert_eq!(ps[1].userflags, 200);
    }

    #[test]
    fn ssw2mcpl_surf_flag_off_clears_userflags() {
        let opts = Ssw2McplOptions {
            surf_to_userflags: false,
            ..Ssw2McplOptions::default()
        };
        let (h, ps) = ssw2mcpl(&synth_tracks(), &opts).unwrap();
        assert!(!h.has_userflags);
        assert!(ps.iter().all(|p| p.userflags == 0));
        // Energies still verbatim.
        assert_eq!(ps[0].ekin, 2.5);
    }

    #[test]
    fn ssw2mcpl_double_prec_and_deck_blob() {
        let opts = Ssw2McplOptions {
            double_prec: true,
            deck_blob: Some(DeckBlob {
                key: "ssw_deck".to_string(),
                data: b"c synthetic deck".to_vec(),
            }),
            srcname: "probe".to_string(),
            comments: vec!["synthetic".to_string()],
            ..Ssw2McplOptions::default()
        };
        let (h, ps) = ssw2mcpl(&synth_tracks(), &opts).unwrap();
        assert!(h.double_prec);
        assert_eq!(h.blobs.len(), 1);
        assert_eq!(h.blobs[0].key, "ssw_deck");
        assert_eq!(h.blobs[0].data, b"c synthetic deck");
        assert_eq!(ps.len(), 2);
        // Byte-exact re-emit through the existing writer.
        let bytes = encode_file(&h, &ps).unwrap();
        let file = crate::McplFile::from_bytes(bytes.clone()).unwrap();
        assert!(file.header.double_prec);
        assert_eq!(file.header.blobs, h.blobs);
        let rewritten = encode_file(&file.header, &file.particles().unwrap()).unwrap();
        assert_eq!(rewritten, bytes);
    }

    #[test]
    fn ssw2mcpl_rejects_oversize_deck_blob() {
        let opts = Ssw2McplOptions {
            deck_blob: Some(DeckBlob {
                key: "big".to_string(),
                data: vec![0u8; MAX_DECK_BLOB_BYTES + 1],
            }),
            ..Ssw2McplOptions::default()
        };
        let err = ssw2mcpl(&synth_tracks(), &opts).unwrap_err();
        assert_eq!(
            err,
            SswError::DeckBlobTooLarge {
                bytes: MAX_DECK_BLOB_BYTES + 1,
                limit: MAX_DECK_BLOB_BYTES,
            }
        );
    }

    #[test]
    fn ssw2mcpl_bytes_gzip_magic_and_round_trip() {
        let opts = Ssw2McplOptions {
            gzip: true,
            ..Ssw2McplOptions::default()
        };
        let bytes = ssw2mcpl_bytes(&synth_tracks(), &opts).unwrap();
        // Gzip magic only (encoder bytes never asserted); decode manually
        // since McplFile::open gunzips on the `.gz` extension only.
        assert_eq!(&bytes[..2], b"\x1f\x8b");
        use std::io::Read as _;
        let mut dec = flate2::read::GzDecoder::new(&bytes[..]);
        let mut raw = Vec::new();
        dec.read_to_end(&mut raw).unwrap();
        let file = crate::McplFile::from_bytes(raw).unwrap();
        let ps = file.particles().unwrap();
        assert_eq!(ps.len(), 2);
        assert!((ps[0].ekin - 2.5).abs() < 1e-9);
        assert_eq!(ps[0].pdgcode, PDG_NEUTRON);
        assert_eq!(ps[0].userflags, 100);
        // Raw (non-gzip) bytes parse directly.
        let raw_opts = Ssw2McplOptions::default();
        let raw = ssw2mcpl_bytes(&synth_tracks(), &raw_opts).unwrap();
        assert_eq!(&raw[..4], b"MCPL");
    }

    #[test]
    fn ssw2mcpl_rejects_bad_tracks() {
        let mut tracks = synth_tracks();
        tracks[0].direction = [1.0, 1.0, 1.0];
        assert!(matches!(
            ssw2mcpl(&tracks, &Ssw2McplOptions::default()).unwrap_err(),
            SswError::NonUnitDirection(0, _)
        ));
        let mut tracks = synth_tracks();
        tracks[1].ekin = -0.5;
        assert_eq!(
            ssw2mcpl(&tracks, &Ssw2McplOptions::default()).unwrap_err(),
            SswError::NegativeEnergy(1, -0.5)
        );
    }

    #[test]
    fn mcpl2ssw_closed_form_userflags_surface() {
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        let (h, tracks) = mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap();
        // Reference passthrough: every string/surface preserved; the count
        // fields name the converted tally (upstream-agreed: np1 == nrss == n,
        // orignp1 keeps the table-2 sign).
        let reference = synth_reference();
        assert_eq!(h.kod, reference.kod);
        assert_eq!(h.ver, reference.ver);
        assert_eq!(h.surflist, reference.surflist);
        assert_eq!(h.np1, 2);
        assert_eq!(h.orignp1, -2);
        assert_eq!(h.nrss, 2);
        assert_eq!(tracks.len(), 2);
        // Neutron: energy verbatim, ms→shakes within f64 rounding of the
        // shakes→ms step above.
        assert_eq!(tracks[0].erg, 2.5);
        assert!((tracks[0].tme - 3.0e5).abs() / 3.0e5 < 1e-12);
        assert_eq!(tracks[0].wgt, 1.0);
        assert_eq!((tracks[0].x, tracks[0].y, tracks[0].z), (1.0, -2.0, 0.5));
        assert_eq!((tracks[0].u, tracks[0].v, tracks[0].cs), (0.0, 0.0, 1.0));
        assert_eq!(tracks[0].w, 1.0);
        assert_eq!(tracks[0].nps, 1.0);
        // Gamma: direction preserved, no cosine forced.
        assert_eq!(tracks[1].erg, 0.662);
        assert_eq!((tracks[1].u, tracks[1].v, tracks[1].cs), (1.0, 0.0, 0.0));
        assert_eq!(tracks[1].w, 0.0);
        assert_eq!(tracks[1].nps, 2.0);
        // Full SSW write + reparse round-trip through the existing writer.
        let mut bytes = Vec::new();
        nucleide_mcnp_io::surfsrc::write_to(&mut bytes, &h, &tracks).unwrap();
        let reparsed = SurfSrc::from_bytes(bytes).unwrap();
        assert_eq!(reparsed.header.nrss, 2);
        assert_eq!(reparsed.read_tracklist().unwrap(), tracks);
    }

    #[test]
    fn mcpl2ssw_negative_z_sign_rides_bitarray() {
        let (_, mut ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        ps[0].direction = [0.0, 0.0, -1.0];
        let (_, tracks) = mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap();
        assert_eq!(tracks[0].cs, -1.0);
        assert_eq!(tracks[0].w, -1.0);
        assert!(tracks[0].bitarray.is_sign_negative());
        assert_eq!(tracks[0].bitarray.abs(), 0.0);
    }

    #[test]
    fn mcpl2ssw_positive_orignp1_stays_positive() {
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        let mut reference = synth_reference();
        reference.orignp1 = 1000;
        reference.np1 = 1000;
        reference.niwr = None;
        reference.mipts = None;
        reference.kjaq = None;
        let (h, tracks) = mcpl2ssw(&ps, &reference, &Mcpl2SswOptions::default()).unwrap();
        assert_eq!((h.np1, h.orignp1, h.nrss), (2, 2, 2));
        assert_eq!(tracks.len(), 2);
    }

    #[test]
    fn mcpl2ssw_surface_override_and_range() {
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        // Explicit override wins over userflags (which stay [100, 200]).
        let opts = Mcpl2SswOptions { surface: Some(7) };
        let (_, tracks) = mcpl2ssw(&ps, &synth_reference(), &opts).unwrap();
        assert_eq!(tracks.len(), 2);
        // Out-of-range override rejected even with valid userflags present.
        for bad in [0, MAX_SURFACE_ID + 1] {
            let opts = Mcpl2SswOptions { surface: Some(bad) };
            assert_eq!(
                mcpl2ssw(&ps, &synth_reference(), &opts).unwrap_err(),
                SswError::SurfaceIdOutOfRange(bad)
            );
        }
        // userflags 0 without override is a named error, not a silent 0.
        let mut noflag = ps.clone();
        noflag[0].userflags = 0;
        assert_eq!(
            mcpl2ssw(&noflag, &synth_reference(), &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::MissingSurfaceId { index: 0 }
        );
    }

    #[test]
    fn mcpl2ssw_rejects_unsupported_pdg_and_width() {
        let (_, mut ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        ps[0].pdgcode = 2212; // proton: named-open, never mistyped as neutron
        assert_eq!(
            mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::UnsupportedPdg(2212)
        );
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        let mut reference = synth_reference();
        reference.ncrd = 10;
        assert_eq!(
            mcpl2ssw(&ps, &reference, &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::UnsupportedRecordWidth {
                expected: TrackData::RECORD_WIDTH,
                found: 10,
            }
        );
        // Particle-level validation mirrors the writer gates.
        let mut bad = ps.clone();
        bad[1].direction = [0.0, 0.0, 0.5];
        assert!(matches!(
            mcpl2ssw(&bad, &synth_reference(), &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::NonUnitDirection(1, _)
        ));
    }

    #[test]
    fn kind_pdg_round_trip() {
        assert_eq!(SswParticleKind::Neutron.pdg(), 2112);
        assert_eq!(SswParticleKind::Gamma.pdg(), 22);
        assert_eq!(
            SswParticleKind::from_pdg(2112).unwrap(),
            SswParticleKind::Neutron
        );
        assert_eq!(
            SswParticleKind::from_pdg(22).unwrap(),
            SswParticleKind::Gamma
        );
        assert_eq!(
            SswParticleKind::from_pdg(11).unwrap_err(),
            SswError::UnsupportedPdg(11)
        );
        assert_eq!(
            SswParticleKind::parse("neutron"),
            Some(SswParticleKind::Neutron)
        );
        assert_eq!(
            SswParticleKind::parse("gamma"),
            Some(SswParticleKind::Gamma)
        );
        assert_eq!(SswParticleKind::parse("proton"), None);
        assert_eq!(SswParticleKind::parse("p"), None);
    }

    #[test]
    fn fixture_reference_round_trips_both_directions() {
        // Committed synthetic pair (see fixtures/mcpl/ssw_conversion/README.md):
        // reference.w holds the two hand-framed tracks; the (surf, kind)
        // pairing is the explicit caller parameter documented there.
        let dir = format!(
            "{}/../../fixtures/mcpl/ssw_conversion",
            env!("CARGO_MANIFEST_DIR")
        );
        let ssw = SurfSrc::open(format!("{dir}/reference.w")).unwrap();
        let raw = ssw.read_tracklist().unwrap();
        assert_eq!(raw.len(), 2);
        let paired = [
            (100u32, SswParticleKind::Neutron),
            (200u32, SswParticleKind::Gamma),
        ];
        let tracks: Vec<SswTrack> = raw
            .iter()
            .zip(paired)
            .map(|(t, (surf, kind))| SswTrack {
                ekin: t.erg,
                time_shakes: t.tme,
                position: [t.x, t.y, t.z],
                direction: [t.u, t.v, t.cs],
                weight: t.wgt,
                surf,
                kind,
            })
            .collect();
        let opts = Ssw2McplOptions::default();
        let bytes = ssw2mcpl_bytes(&tracks, &opts).unwrap();
        let expected = std::fs::read(format!("{dir}/ssw2mcpl_expected.mcpl")).unwrap();
        assert_eq!(
            bytes, expected,
            "conversion output drifted from golden bytes"
        );
        // And back: MCPL → SSW against the same reference header.
        let file = crate::McplFile::from_bytes(bytes).unwrap();
        let (h, back) = mcpl2ssw(
            &file.particles().unwrap(),
            &ssw.header,
            &Mcpl2SswOptions::default(),
        )
        .unwrap();
        assert_eq!(h.nrss, 2);
        for (rt, want) in back.iter().zip(raw.iter()) {
            // Single-prec MCPL packing quantizes energy/time (abs 1e-6
            // covers the f32 step); geometry rides on exact 0/1 packing.
            assert!((rt.erg - want.erg).abs() < 1e-6);
            assert!((rt.tme - want.tme).abs() < 1e-6);
            assert!((rt.wgt - want.wgt).abs() < 1e-9);
            assert_eq!((rt.x, rt.y, rt.z), (want.x, want.y, want.z));
            assert_eq!((rt.u, rt.v, rt.cs), (want.u, want.v, want.cs));
        }
    }
}
