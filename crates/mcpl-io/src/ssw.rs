//! SSW surface-source ↔ MCPL particle-list conversion (SSW-PDG table v2).
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
//! | particle kind | `pdgcode` | SSW-PDG table ([`SswParticleKind`]): neutron 2112, gamma 22, electron 11, positron −11, proton 2212; anything else is [`SswError::UnsupportedPdg`] |
//! | `wgt`, `x`/`y`/`z`, `u`/`v`/`cs` | `weight`, `position`, `direction` | verbatim, except opt-in [`Mcpl2SswOptions::force_cs_to_one`] |
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
//!   reference), so `kind` only gates the conversion: PDG codes outside the
//!   SSW-PDG table are rejected, never mistyped.
//! - `bitarray` magnitude is 0 (the track's cell and upstream type word are
//!   unknown to the converter); its sign mirrors the direction z so the
//!   reader's `copysign`-derived `w` matches `cs`. Third-party tools that
//!   decode the particle type out of `bitarray` will skip these tracks —
//!   use the nucleide leg for round-trips (named-open, see below).
//! - `nps` is the 1-based particle index (MCPL carries no history number).
//! - No direction cosine is forced by default: the local SSW writer stores
//!   track records verbatim, so `(u, v, cs)` propagate the particle direction
//!   unchanged. Upstream `mcpl2ssw` instead forces the stored `cs` slot to
//!   1.0 (verified against the 2.2.8 binary); opt in with
//!   [`Mcpl2SswOptions::force_cs_to_one`] to reproduce that byte (`cs = 1.0`,
//!   `u`/`v` still verbatim) instead of keeping the true cosine.
//! - The output header clones the reference (code/version/deck passthrough)
//!   with the count fields patched to the particle count — `nrss`, `np1`,
//!   and `orignp1` with the reference's table-2 sign preserved — exactly as
//!   upstream `mcpl2ssw` reports it ("N particles (nrss) and N histories
//!   (np1)", verified against the 2.2.8 binary). `niss` passes through by
//!   default; set [`Mcpl2SswOptions::niss_override`] to stamp an explicit
//!   value instead (upstream 2.2.8 leaves `niss` on the reference header
//!   untouched — verified by running the 2.2.8 `mcpl2ssw` script over the
//!   synthetic pair and diffing the header — so the passthrough default is
//!   the compatible spelling and the override is the tally-convention opt-in).
//! - Polarisation has no SSW slot: MCPL → SSW drops it, loudly by default
//!   (any non-zero input vector is [`SswError::PolarisationPresent`] unless
//!   [`Mcpl2SswOptions::allow_polarisation`] opts into the drop). The
//!   reverse leg carries a uniform vector only when
//!   [`Ssw2McplOptions::polarisation`] is set (otherwise zeros with
//!   `has_polarisation == false`). Universal PDG/weight inputs arrive
//!   already unfolded by the MCPL reader, so their decoded per-particle
//!   values propagate verbatim; the reverse leg emits file-wide codes only
//!   when [`Ssw2McplOptions::universal_pdg`]/[`universal_weight`](Ssw2McplOptions::universal_weight)
//!   opt in (mixed inputs are [`SswError::MixedPdgForUniversal`]/
//!   [`SswError::MixedWeightForUniversal`]).
//!
//! ## Named-open items (explicitly out of scope)
//!
//! - Particle types beyond the SSW-PDG table ([`SswError::UnsupportedPdg`]).
//!   Accepted set decision: neutron/gamma/electron/positron/proton — the
//!   five MCPL PDG codes with identical 11-double SSW geometry (energy,
//!   time, position, direction all verbatim). Heavier ions and mesons stay
//!   loud errors, never silent skips; no transport semantics ride along.
//! - Transport semantics; this only re-homes already-transported tracks.
//! - Reference headers whose `abs(ncrd)` is not
//!   [`TrackData::RECORD_WIDTH`](nucleide_mcnp_io::surfsrc::TrackData::RECORD_WIDTH)
//!   ([`SswError::UnsupportedRecordWidth`]).
//! - Upstream type-word round-trips: this converter never writes cell/type words into
//!   `bitarray` (writing them would transcribe the upstream rawtype table
//!   into this crate, which the license-clean cut forbids), so
//!   third-party `bitarray`-decoding tools skip converter-written tracks.
//!   Owner/next step: revisit when a licensed type-table source exists.

use nucleide_mcnp_io::surfsrc::{SurfSrcHeader, TrackData};
use thiserror::Error;

use crate::{encode_file, Blob, Header, Particle};

/// Shakes → milliseconds (1 shake = 1e-8 s = 1e-5 ms).
pub const SHAKES_TO_MS: f64 = 1e-5;
/// Milliseconds → shakes.
pub const MS_TO_SHAKES: f64 = 1e5;
/// PDG code for the neutron.
pub const PDG_NEUTRON: i32 = 2112;
/// PDG code for the gamma.
pub const PDG_GAMMA: i32 = 22;
/// PDG code for the electron.
pub const PDG_ELECTRON: i32 = 11;
/// PDG code for the positron.
pub const PDG_POSITRON: i32 = -11;
/// PDG code for the proton.
pub const PDG_PROTON: i32 = 2212;
/// Largest deck blob accepted by [`Ssw2McplOptions::deck_blob`] (100 MiB).
pub const MAX_DECK_BLOB_BYTES: usize = 100 * 1024 * 1024;
/// Largest surface id accepted on the MCPL → SSW path.
pub const MAX_SURFACE_ID: u32 = 999_999;

/// SSW particle kind (SSW-PDG table: n/γ/e⁻/e⁺/p).
///
/// The kind is an explicit caller parameter on [`SswTrack`]: no particle-type
/// decode exists in the local SSW reader, so the converter never guesses it
/// from `bitarray` (see the module docs). All five kinds share the identical
/// 11-double SSW geometry (energy, time, position, direction verbatim); the
/// kind only gates the PDG mapping, never the layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SswParticleKind {
    /// Neutron (PDG 2112).
    Neutron,
    /// Gamma (PDG 22).
    Gamma,
    /// Electron (PDG 11).
    Electron,
    /// Positron (PDG −11).
    Positron,
    /// Proton (PDG 2212).
    Proton,
}

impl SswParticleKind {
    /// All accepted kinds in stable facade order.
    pub const ALL: [SswParticleKind; 5] = [
        SswParticleKind::Neutron,
        SswParticleKind::Gamma,
        SswParticleKind::Electron,
        SswParticleKind::Positron,
        SswParticleKind::Proton,
    ];

    /// PDG code for this kind.
    pub fn pdg(self) -> i32 {
        match self {
            SswParticleKind::Neutron => PDG_NEUTRON,
            SswParticleKind::Gamma => PDG_GAMMA,
            SswParticleKind::Electron => PDG_ELECTRON,
            SswParticleKind::Positron => PDG_POSITRON,
            SswParticleKind::Proton => PDG_PROTON,
        }
    }

    /// Map a PDG code back to a kind; anything else is [`SswError::UnsupportedPdg`].
    pub fn from_pdg(pdg: i32) -> Result<Self, SswError> {
        match pdg {
            PDG_NEUTRON => Ok(SswParticleKind::Neutron),
            PDG_GAMMA => Ok(SswParticleKind::Gamma),
            PDG_ELECTRON => Ok(SswParticleKind::Electron),
            PDG_POSITRON => Ok(SswParticleKind::Positron),
            PDG_PROTON => Ok(SswParticleKind::Proton),
            other => Err(SswError::UnsupportedPdg(other)),
        }
    }

    /// Accepted spellings (the only strings the Python facade accepts).
    pub fn as_str(self) -> &'static str {
        match self {
            SswParticleKind::Neutron => "neutron",
            SswParticleKind::Gamma => "gamma",
            SswParticleKind::Electron => "electron",
            SswParticleKind::Positron => "positron",
            SswParticleKind::Proton => "proton",
        }
    }

    /// Parse an accepted spelling; anything else is `None` (the facade turns
    /// this into a `ValueError` naming the accepted values).
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "neutron" => Some(SswParticleKind::Neutron),
            "gamma" => Some(SswParticleKind::Gamma),
            "electron" => Some(SswParticleKind::Electron),
            "positron" => Some(SswParticleKind::Positron),
            "proton" => Some(SswParticleKind::Proton),
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
    /// Uniform polarisation stamped on every output particle (`None` =
    /// zeros with `has_polarisation == false`, the default). `Some(v)`
    /// sets `has_polarisation` and stores `v` per particle; `v` must be
    /// finite (see [`SswError::InvalidPolarisation`]).
    pub polarisation: Option<[f64; 3]>,
    /// Emit a file-wide PDG code (`None` = per-particle codes, the default).
    /// `true` requires every track to share one kind (see
    /// [`SswError::MixedPdgForUniversal`]).
    pub universal_pdg: bool,
    /// Emit a file-wide weight (`true` requires every track weight to be
    /// exactly equal; see [`SswError::MixedWeightForUniversal`]).
    pub universal_weight: bool,
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
            polarisation: None,
            universal_pdg: false,
            universal_weight: false,
        }
    }
}

/// Options for [`mcpl2ssw`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mcpl2SswOptions {
    /// Surface id stamped on every output track (`None` = each particle's
    /// `userflags`; must lie in `[1, 999999]` either way).
    pub surface: Option<u32>,
    /// Reproduce the upstream `mcpl2ssw` 2.2.8 `cs` spelling: force the
    /// stored `cs` slot to `1.0` (`u`/`v` stay verbatim). Default `false`
    /// keeps the true direction cosine.
    pub force_cs_to_one: bool,
    /// Stamp `niss` on the output header (`None` = pass the reference
    /// header's `niss` through untouched, the default and the upstream-2.2.8
    /// spelling). `Some(v)` requires `v >= 0` (see
    /// [`SswError::NissOutOfRange`]).
    pub niss_override: Option<i64>,
    /// Allow polarised MCPL inputs (dropped, since SSW has no slot).
    /// Default `false` rejects any non-zero input vector with
    /// [`SswError::PolarisationPresent`]; `true` drops silently.
    pub allow_polarisation: bool,
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
    /// Particle kind outside the SSW-PDG table.
    #[error("unsupported PDG code {0} (accepted: 2112 neutron, 22 gamma, 11 electron, -11 positron, 2212 proton)")]
    UnsupportedPdg(i32),
    /// `universal_pdg` requested but the tracks mix PDG codes.
    #[error("universal PDG requested but tracks mix PDG codes {0:?} (pass one kind or leave universal_pdg off)")]
    MixedPdgForUniversal(Vec<i32>),
    /// `universal_weight` requested but the track weights differ.
    #[error("universal weight requested but track weights differ (pass equal weights or leave universal_weight off)")]
    MixedWeightForUniversal,
    /// `polarisation` vector holds a non-finite entry.
    #[error("polarisation vector {0:?} holds a non-finite entry")]
    InvalidPolarisation([f64; 3]),
    /// MCPL → SSW input carries polarisation with no SSW slot to hold it.
    #[error("particle {index} carries non-zero polarisation (pass allow_polarisation to drop it)")]
    PolarisationPresent {
        /// Position of the offending particle in the input slice.
        index: usize,
    },
    /// `niss_override` is negative.
    #[error("niss override {0} is negative (pass >= 0 or leave unset for passthrough)")]
    NissOutOfRange(i64),
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
    #[error("reference ncrd width {found} != {expected} (converts 11-double tracks only)")]
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
/// would. Opt-in [`Ssw2McplOptions::polarisation`] stamps one uniform vector
/// on every particle (otherwise zeros); opt-in
/// [`Ssw2McplOptions::universal_pdg`]/[`universal_weight`](Ssw2McplOptions::universal_weight)
/// collapse the file-wide codes (mixed inputs are named errors, never silent
/// per-particle fallbacks).
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
    if let Some(pol) = options.polarisation {
        if !pol.iter().all(|v| v.is_finite()) {
            return Err(SswError::InvalidPolarisation(pol));
        }
    }
    if options.universal_pdg && !tracks.is_empty() {
        let first = tracks[0].kind.pdg();
        if tracks.iter().any(|t| t.kind.pdg() != first) {
            let mut pdgs: Vec<i32> = tracks.iter().map(|t| t.kind.pdg()).collect();
            pdgs.sort_unstable();
            pdgs.dedup();
            return Err(SswError::MixedPdgForUniversal(pdgs));
        }
    }
    if options.universal_weight && !tracks.is_empty() {
        let first = tracks[0].weight;
        if tracks.iter().any(|t| t.weight != first) {
            return Err(SswError::MixedWeightForUniversal);
        }
    }
    let pol = options.polarisation.unwrap_or([0.0; 3]);
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
            polarisation: pol,
            position: t.position,
            direction: t.direction,
            time: t.time_shakes * SHAKES_TO_MS,
            weight: t.weight,
            pdgcode: t.kind.pdg(),
            userflags: if options.surf_to_userflags { t.surf } else { 0 },
        });
    }
    let universal_pdgcode = if options.universal_pdg && !tracks.is_empty() {
        Some(tracks[0].kind.pdg())
    } else {
        None
    };
    let universal_weight = if options.universal_weight && !tracks.is_empty() {
        Some(tracks[0].weight)
    } else {
        None
    };
    let header = Header {
        has_userflags: options.surf_to_userflags,
        has_polarisation: options.polarisation.is_some(),
        double_prec: options.double_prec,
        universal_pdgcode,
        universal_weight,
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
/// `niss` passes through unless [`Mcpl2SswOptions::niss_override`] stamps an
/// explicit value. Surface ids come from
/// `userflags` or [`Mcpl2SswOptions::surface`]; PDG codes outside the
/// SSW-PDG table are [`SswError::UnsupportedPdg`]. Energy maps verbatim
/// (MeV), time maps ms → shakes (`× 1e5`), `(u, v, cs)` carry the particle
/// direction unchanged unless [`Mcpl2SswOptions::force_cs_to_one`] reproduces
/// the upstream `cs = 1.0` spelling. Non-zero polarisation is
/// [`SswError::PolarisationPresent`] unless
/// [`Mcpl2SswOptions::allow_polarisation`] drops it.
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
    if let Some(niss) = options.niss_override {
        if niss < 0 {
            return Err(SswError::NissOutOfRange(niss));
        }
    }
    let mut tracks = Vec::with_capacity(particles.len());
    for (i, p) in particles.iter().enumerate() {
        // Gate the kind even though the 11-double record has no per-track
        // type slot: mistyping an ion as a neutron must be loud, not silent.
        let _kind = SswParticleKind::from_pdg(p.pdgcode)?;
        let _surf = resolve_surface(i, p, options)?;
        if !options.allow_polarisation && p.polarisation != [0.0, 0.0, 0.0] {
            return Err(SswError::PolarisationPresent { index: i });
        }
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
        // direction z so the reader's copysign-derived `w` matches `cs`
        // (under force_cs_to_one the stored cs diverges by upstream design).
        let bitarray = if p.direction[2].is_sign_negative() {
            -0.0
        } else {
            0.0
        };
        let cs = if options.force_cs_to_one {
            1.0
        } else {
            p.direction[2]
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
            cs,
        ];
        tracks.push(TrackData::from_record(record));
    }
    let mut header = reference.clone();
    let n = particles.len() as i64;
    header.nrss = n;
    // Upstream patches the history counts to the converted tally (verified:
    // "N particles (nrss) and N histories (np1)"); the table-2 sign rides
    // along so the cloned deck layout stays writable. niss passes through
    // (upstream 2.2.8 leaves it untouched — verified over the synthetic pair)
    // unless the caller stamps an explicit override.
    header.np1 = n;
    header.orignp1 = if reference.orignp1 < 0 { -n } else { n };
    if let Some(niss) = options.niss_override {
        header.niss = niss;
    }
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
        let opts = Mcpl2SswOptions {
            surface: Some(7),
            ..Default::default()
        };
        let (_, tracks) = mcpl2ssw(&ps, &synth_reference(), &opts).unwrap();
        assert_eq!(tracks.len(), 2);
        // Out-of-range override rejected even with valid userflags present.
        for bad in [0, MAX_SURFACE_ID + 1] {
            let opts = Mcpl2SswOptions {
                surface: Some(bad),
                ..Default::default()
            };
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
        ps[0].pdgcode = 211; // pion: outside the SSW-PDG table, never mistyped
        assert_eq!(
            mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::UnsupportedPdg(211)
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
        assert_eq!(SswParticleKind::Electron.pdg(), 11);
        assert_eq!(SswParticleKind::Positron.pdg(), -11);
        assert_eq!(SswParticleKind::Proton.pdg(), 2212);
        for kind in SswParticleKind::ALL {
            assert_eq!(SswParticleKind::from_pdg(kind.pdg()).unwrap(), kind);
            assert_eq!(SswParticleKind::parse(kind.as_str()), Some(kind));
        }
        assert_eq!(
            SswParticleKind::from_pdg(211).unwrap_err(),
            SswError::UnsupportedPdg(211)
        );
        assert_eq!(
            SswParticleKind::from_pdg(1_000_020_040).unwrap_err(),
            SswError::UnsupportedPdg(1_000_020_040)
        );
        assert_eq!(SswParticleKind::parse("pion"), None);
        assert_eq!(SswParticleKind::parse("p"), None);
    }

    #[test]
    fn extended_kinds_round_trip_both_directions_verbatim() {
        // All five SSW-PDG kinds share the identical 11-double geometry:
        // energy/time/position/direction/weight ride verbatim both ways.
        let tracks: Vec<SswTrack> = SswParticleKind::ALL
            .iter()
            .enumerate()
            .map(|(i, kind)| SswTrack {
                ekin: 1.0 + i as f64 * 0.5,
                time_shakes: 1.0e5 * (i as f64 + 1.0),
                position: [i as f64, -(i as f64), 0.25 * i as f64],
                direction: [0.0, 0.0, 1.0],
                weight: 0.75,
                surf: 100 + i as u32,
                kind: *kind,
            })
            .collect();
        let (h, ps) = ssw2mcpl(&tracks, &Ssw2McplOptions::default()).unwrap();
        assert_eq!(
            ps.iter().map(|p| p.pdgcode).collect::<Vec<_>>(),
            vec![2112, 22, 11, -11, 2212]
        );
        assert!(!h.has_polarisation);
        let (oh, back) = mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap();
        assert_eq!((oh.nrss, oh.np1), (5, 5));
        assert_eq!(back.len(), 5);
        for (rt, want) in back.iter().zip(tracks.iter()) {
            assert_eq!(rt.erg, want.ekin);
            assert!((rt.tme - want.time_shakes).abs() / want.time_shakes < 1e-12);
            assert_eq!(rt.wgt, want.weight);
            assert_eq!((rt.u, rt.v, rt.cs), (0.0, 0.0, 1.0));
        }
        // Full SSW write + reparse round-trip stays stable.
        let mut bytes = Vec::new();
        nucleide_mcnp_io::surfsrc::write_to(&mut bytes, &oh, &back).unwrap();
        let reparsed = SurfSrc::from_bytes(bytes).unwrap();
        assert_eq!(reparsed.read_tracklist().unwrap(), back);
    }

    #[test]
    fn mcpl2ssw_force_cs_to_one_matches_upstream_spelling() {
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        // Default keeps the true cosine (gamma track rides +x, cs == 0).
        let (_, verbatim) = mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap();
        assert_eq!(verbatim[1].cs, 0.0);
        assert_eq!((verbatim[1].u, verbatim[1].v), (1.0, 0.0));
        // Opt-in reproduces the upstream 2.2.8 spelling: cs forced to 1.0,
        // u/v still verbatim.
        let opts = Mcpl2SswOptions {
            force_cs_to_one: true,
            ..Default::default()
        };
        let (_, forced) = mcpl2ssw(&ps, &synth_reference(), &opts).unwrap();
        assert_eq!(forced[0].cs, 1.0);
        assert_eq!(forced[1].cs, 1.0);
        assert_eq!((forced[1].u, forced[1].v), (1.0, 0.0));
        // And back through the writer/reparse leg the forced value persists.
        let (h, _) = mcpl2ssw(&ps, &synth_reference(), &opts).unwrap();
        let mut bytes = Vec::new();
        nucleide_mcnp_io::surfsrc::write_to(&mut bytes, &h, &forced).unwrap();
        let reparsed = SurfSrc::from_bytes(bytes).unwrap();
        assert_eq!(reparsed.read_tracklist().unwrap()[1].cs, 1.0);
    }

    #[test]
    fn mcpl2ssw_niss_passthrough_default_override_opt_in() {
        let (_, ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        // Default: reference niss rides through untouched (upstream-2.2.8
        // spelling, verified over the synthetic pair).
        let mut reference = synth_reference();
        reference.niss = 17;
        let (h, _) = mcpl2ssw(&ps, &reference, &Mcpl2SswOptions::default()).unwrap();
        assert_eq!(h.niss, 17);
        // Opt-in: stamp an explicit non-negative value.
        let opts = Mcpl2SswOptions {
            niss_override: Some(2),
            ..Default::default()
        };
        let (h, back) = mcpl2ssw(&ps, &reference, &opts).unwrap();
        assert_eq!(h.niss, 2);
        assert_eq!(back.len(), 2);
        // Negative override is a named error, never a silent wrap.
        let bad = Mcpl2SswOptions {
            niss_override: Some(-1),
            ..Default::default()
        };
        assert_eq!(
            mcpl2ssw(&ps, &reference, &bad).unwrap_err(),
            SswError::NissOutOfRange(-1)
        );
    }

    #[test]
    fn ssw2mcpl_polarisation_and_universal_opt_in() {
        // Polarisation opt-in stamps one uniform vector + the header flag.
        let opts = Ssw2McplOptions {
            polarisation: Some([0.1, 0.2, 0.3]),
            ..Ssw2McplOptions::default()
        };
        let (h, ps) = ssw2mcpl(&synth_tracks(), &opts).unwrap();
        assert!(h.has_polarisation);
        assert!(ps.iter().all(|p| p.polarisation == [0.1, 0.2, 0.3]));
        // Byte-exact re-emit through the MCPL writer (polarised layout).
        let bytes = encode_file(&h, &ps).unwrap();
        let file = crate::McplFile::from_bytes(bytes.clone()).unwrap();
        assert!(file.header.has_polarisation);
        let rewritten = encode_file(&file.header, &file.particles().unwrap()).unwrap();
        assert_eq!(rewritten, bytes);
        // Non-finite vectors are a named error.
        let bad = Ssw2McplOptions {
            polarisation: Some([f64::NAN, 0.0, 0.0]),
            ..Ssw2McplOptions::default()
        };
        assert!(matches!(
            ssw2mcpl(&synth_tracks(), &bad).unwrap_err(),
            SswError::InvalidPolarisation(_)
        ));
        // Universal PDG opt-in over a single-kind pair.
        let single_kind = vec![
            SswTrack {
                kind: SswParticleKind::Neutron,
                ..synth_tracks()[0].clone()
            },
            SswTrack {
                ekin: 0.5,
                kind: SswParticleKind::Neutron,
                ..synth_tracks()[1].clone()
            },
        ];
        let uopts = Ssw2McplOptions {
            universal_pdg: true,
            ..Ssw2McplOptions::default()
        };
        let (uh, ups) = ssw2mcpl(&single_kind, &uopts).unwrap();
        assert_eq!(uh.universal_pdgcode, Some(PDG_NEUTRON));
        assert_eq!(ups.len(), 2);
        // Mixed kinds under universal_pdg are loud, never silent per-particle.
        assert_eq!(
            ssw2mcpl(&synth_tracks(), &uopts).unwrap_err(),
            SswError::MixedPdgForUniversal(vec![22, 2112])
        );
        // Universal weight opt-in over equal weights; mixed weights are loud.
        let equal_w = vec![
            SswTrack {
                weight: 1.0,
                ..synth_tracks()[0].clone()
            },
            SswTrack {
                weight: 1.0,
                surf: 200,
                kind: SswParticleKind::Gamma,
                ..synth_tracks()[1].clone()
            },
        ];
        let wopts = Ssw2McplOptions {
            universal_weight: true,
            ..Ssw2McplOptions::default()
        };
        let (wh, _) = ssw2mcpl(&equal_w, &wopts).unwrap();
        assert_eq!(wh.universal_weight, Some(1.0));
        assert_eq!(
            ssw2mcpl(&synth_tracks(), &wopts).unwrap_err(),
            SswError::MixedWeightForUniversal
        );
    }

    #[test]
    fn mcpl2ssw_polarisation_gate_and_allow_opt_in() {
        let (h, mut ps) = ssw2mcpl(&synth_tracks(), &Ssw2McplOptions::default()).unwrap();
        assert!(!h.has_polarisation);
        ps[0].polarisation = [0.0, 0.0, 1.0];
        // Default: loud, never a silent drop.
        assert_eq!(
            mcpl2ssw(&ps, &synth_reference(), &Mcpl2SswOptions::default()).unwrap_err(),
            SswError::PolarisationPresent { index: 0 }
        );
        // Opt-in: dropped (SSW has no slot), geometry still verbatim.
        let opts = Mcpl2SswOptions {
            allow_polarisation: true,
            ..Default::default()
        };
        let (_, tracks) = mcpl2ssw(&ps, &synth_reference(), &opts).unwrap();
        assert_eq!(tracks.len(), 2);
        assert_eq!(tracks[0].erg, 2.5);
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
        // Committed goldens for the fidelity-tail opt-ins (same reference.w,
        // documented pairings in fixtures/mcpl/ssw_conversion/README.md).
        for (name, kinds, opts) in [
            (
                "ssw2mcpl_extended_expected.mcpl",
                vec!["electron", "positron"],
                Ssw2McplOptions::default(),
            ),
            (
                "ssw2mcpl_polarised_expected.mcpl",
                vec!["neutron", "gamma"],
                Ssw2McplOptions {
                    polarisation: Some([0.1, 0.2, 0.3]),
                    ..Ssw2McplOptions::default()
                },
            ),
            (
                "ssw2mcpl_universal_pdg_expected.mcpl",
                vec!["neutron", "neutron"],
                Ssw2McplOptions {
                    universal_pdg: true,
                    ..Ssw2McplOptions::default()
                },
            ),
        ] {
            let paired: Vec<(u32, SswParticleKind)> = kinds
                .iter()
                .zip([100u32, 200u32])
                .map(|(k, s)| (s, SswParticleKind::parse(k).unwrap()))
                .collect();
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
            let bytes = ssw2mcpl_bytes(&tracks, &opts).unwrap();
            let expected = std::fs::read(format!("{dir}/{name}")).unwrap();
            assert_eq!(bytes, expected, "{name} drifted from golden bytes");
        }
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
