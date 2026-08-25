//! Particle names and cross-code translation (MCNP/MCNP6/FLUKA/GEANT4).
//!
//! The canonical set is the 32 fundamental particles numbered by the Berkeley
//! Particle Data Centre (PDC) scheme. Heavy ions are valid *specs*
//! (see [`is_heavy_ion`]) but are not representable as a [`ParticleId`] —
//! they carry PDC number 0.
//!
//! Registry notes:
//! - `AntiSigmaMinus` carries PDC `-3112` (duplicating `3112` would be
//!   physically wrong and break id ↔ name injectivity).
//! - FLUKA/Geant4 strings with stray spaces in some legacy tables
//!   (`"Anti Tauon"`, `"KaonZero Short"`, ...) are attached to their intended
//!   variants, so the translations actually resolve.
//! - `AntiKaonZero` has no Geant4 string in legacy tables and keeps none here.

use std::fmt;

use crate::NuclideId;

/// Fundamental particles (Berkeley PDC numbering).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ParticleId {
    Electron,
    Positron,
    ElectronNeutrino,
    ElectronAntiNeutrino,
    Muon,
    AntiMuon,
    MuonNeutrino,
    MuonAntiNeutrino,
    Tauon,
    AntiTauon,
    TauNeutrino,
    TauAntiNeutrino,
    /// Gauge boson; aliases `"Gamma"` / `"X-Ray"` / `"g"`.
    Photon,
    /// Charged pion (identified with the negative pion).
    Pion,
    AntiPion,
    Kaon,
    AntiKaon,
    KaonZeroShort,
    KaonZero,
    AntiKaonZero,
    Neutron,
    AntiNeutron,
    /// Aliases `"Hydrogen"` / `"Protium"` / `"p"`; also any H-1 nuclide spec.
    Proton,
    AntiProton,
    Lambda,
    AntiLambda,
    SigmaMinus,
    AntiSigmaMinus,
    SigmaPlus,
    AntiSigmaPlus,
    SigmaZero,
    AntiSigmaZero,
}

/// Every canonical particle, in declaration order.
pub const ALL: [ParticleId; 32] = [
    ParticleId::Electron,
    ParticleId::Positron,
    ParticleId::ElectronNeutrino,
    ParticleId::ElectronAntiNeutrino,
    ParticleId::Muon,
    ParticleId::AntiMuon,
    ParticleId::MuonNeutrino,
    ParticleId::MuonAntiNeutrino,
    ParticleId::Tauon,
    ParticleId::AntiTauon,
    ParticleId::TauNeutrino,
    ParticleId::TauAntiNeutrino,
    ParticleId::Photon,
    ParticleId::Pion,
    ParticleId::AntiPion,
    ParticleId::Kaon,
    ParticleId::AntiKaon,
    ParticleId::KaonZeroShort,
    ParticleId::KaonZero,
    ParticleId::AntiKaonZero,
    ParticleId::Neutron,
    ParticleId::AntiNeutron,
    ParticleId::Proton,
    ParticleId::AntiProton,
    ParticleId::Lambda,
    ParticleId::AntiLambda,
    ParticleId::SigmaMinus,
    ParticleId::AntiSigmaMinus,
    ParticleId::SigmaPlus,
    ParticleId::AntiSigmaPlus,
    ParticleId::SigmaZero,
    ParticleId::AntiSigmaZero,
];

/// Alternate spellings accepted by [`ParticleId::parse`] beyond the canonical
/// names from the legacy alternate-name table, plus the common one-letter physics symbols.
static ALIASES: [(&str, ParticleId); 13] = [
    ("Hydrogen", ParticleId::Proton),
    ("Protium", ParticleId::Proton),
    ("Beta", ParticleId::Electron),
    ("Beta-", ParticleId::Electron),
    ("Beta+", ParticleId::Positron),
    ("Gamma", ParticleId::Photon),
    ("X-Ray", ParticleId::Photon),
    ("n", ParticleId::Neutron),
    ("p", ParticleId::Proton),
    ("e", ParticleId::Electron),
    ("e-", ParticleId::Electron),
    ("e+", ParticleId::Positron),
    ("g", ParticleId::Photon),
];

impl ParticleId {
    /// Canonical CamelCase name (`"Neutron"`, `"AntiProton"`, ...).
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Electron => "Electron",
            Self::Positron => "Positron",
            Self::ElectronNeutrino => "ElectronNeutrino",
            Self::ElectronAntiNeutrino => "ElectronAntiNeutrino",
            Self::Muon => "Muon",
            Self::AntiMuon => "AntiMuon",
            Self::MuonNeutrino => "MuonNeutrino",
            Self::MuonAntiNeutrino => "MuonAntiNeutrino",
            Self::Tauon => "Tauon",
            Self::AntiTauon => "AntiTauon",
            Self::TauNeutrino => "TauNeutrino",
            Self::TauAntiNeutrino => "TauAntiNeutrino",
            Self::Photon => "Photon",
            Self::Pion => "Pion",
            Self::AntiPion => "AntiPion",
            Self::Kaon => "Kaon",
            Self::AntiKaon => "AntiKaon",
            Self::KaonZeroShort => "KaonZeroShort",
            Self::KaonZero => "KaonZero",
            Self::AntiKaonZero => "AntiKaonZero",
            Self::Neutron => "Neutron",
            Self::AntiNeutron => "AntiNeutron",
            Self::Proton => "Proton",
            Self::AntiProton => "AntiProton",
            Self::Lambda => "Lambda",
            Self::AntiLambda => "AntiLambda",
            Self::SigmaMinus => "Sigma-",
            Self::AntiSigmaMinus => "AntiSigma-",
            Self::SigmaPlus => "Sigma+",
            Self::AntiSigmaPlus => "AntiSigma+",
            Self::SigmaZero => "Sigma",
            Self::AntiSigmaZero => "AntiSigmaZero",
        }
    }

    /// Human-readable description.
    #[must_use]
    pub const fn describe(self) -> &'static str {
        match self {
            Self::Electron => "Electron",
            Self::Positron => "Positron",
            Self::ElectronNeutrino => "Electron Neutrino",
            Self::ElectronAntiNeutrino => "Electron Anti Neutrino",
            Self::Muon => "Muon",
            Self::AntiMuon => "Anti Muon",
            Self::MuonNeutrino => "Muon Neutrino",
            Self::MuonAntiNeutrino => "Muon Anti Neutrino",
            Self::Tauon => "Tauon",
            Self::AntiTauon => "Anti Tauon",
            Self::TauNeutrino => "Tau Neutrino",
            Self::TauAntiNeutrino => "Tau Anti Neutrino",
            Self::Photon => "Photon",
            Self::Pion => "Pion",
            Self::AntiPion => "Anti Pion",
            Self::Kaon => "Kaon",
            Self::AntiKaon => "Anti Kaon",
            Self::KaonZeroShort => "Kaon Zero Short",
            Self::KaonZero => "Kaon Zero",
            Self::AntiKaonZero => "Anti Kaon Zero",
            Self::Neutron => "Neutron",
            Self::AntiNeutron => "Anti Neutron",
            Self::Proton => "Proton",
            Self::AntiProton => "Anti Proton",
            Self::Lambda => "Lambda",
            Self::AntiLambda => "Anti Lambda",
            Self::SigmaMinus => "Sigma-",
            Self::AntiSigmaMinus => "Anti Sigma-",
            Self::SigmaPlus => "Sigma+",
            Self::AntiSigmaPlus => "Anti Sigma+",
            Self::SigmaZero => "Sigma",
            Self::AntiSigmaZero => "Anti Sigma Zero",
        }
    }

    /// Berkeley PDC integer id.
    #[must_use]
    pub const fn pdc(self) -> i32 {
        match self {
            Self::Electron => 11,
            Self::Positron => -11,
            Self::ElectronNeutrino => 12,
            Self::ElectronAntiNeutrino => -12,
            Self::Muon => 13,
            Self::AntiMuon => -13,
            Self::MuonNeutrino => 14,
            Self::MuonAntiNeutrino => -14,
            Self::Tauon => 15,
            Self::AntiTauon => -15,
            Self::TauNeutrino => 16,
            Self::TauAntiNeutrino => -16,
            Self::Photon => 22,
            Self::Pion => 211,
            Self::AntiPion => -211,
            Self::Kaon => 321,
            Self::AntiKaon => -321,
            Self::KaonZeroShort => 310,
            Self::KaonZero => 311,
            Self::AntiKaonZero => -311,
            Self::Neutron => 2112,
            Self::AntiNeutron => -2112,
            Self::Proton => 2212,
            Self::AntiProton => -2212,
            Self::Lambda => 3122,
            Self::AntiLambda => -3122,
            Self::SigmaMinus => 3112,
            Self::AntiSigmaMinus => -3112,
            Self::SigmaPlus => 3222,
            Self::AntiSigmaPlus => -3222,
            Self::SigmaZero => 3212,
            Self::AntiSigmaZero => -3212,
        }
    }

    /// Particle from a PDC number, or `None`.
    #[must_use]
    pub fn from_pdc(n: i32) -> Option<Self> {
        ALL.iter().copied().find(|p| p.pdc() == n)
    }

    /// MCNP5 particle designator; `None` when MCNP5 cannot score it
    /// (only `n`, `p`, `e`; `"?"` otherwise).
    #[must_use]
    pub const fn mcnp(self) -> Option<&'static str> {
        match self {
            Self::Neutron => Some("n"),
            Self::Photon => Some("p"),
            Self::Electron => Some("e"),
            _ => None,
        }
    }

    /// MCNP6 particle designator; adds proton `"h"` over MCNP5.
    #[must_use]
    pub const fn mcnp6(self) -> Option<&'static str> {
        match self {
            Self::Neutron => Some("n"),
            Self::Photon => Some("p"),
            Self::Electron => Some("e"),
            Self::Proton => Some("h"),
            _ => None,
        }
    }

    /// FLUKA particle name; `None` when unlisted upstream.
    #[must_use]
    pub const fn fluka(self) -> Option<&'static str> {
        match self {
            Self::Electron => Some("ELECTRON"),
            Self::Positron => Some("POSITRON"),
            Self::ElectronNeutrino => Some("NEUTRIE"),
            Self::ElectronAntiNeutrino => Some("ANEUTRIE"),
            Self::Muon => Some("MUON+"),
            Self::AntiMuon => Some("MUON-"),
            Self::MuonNeutrino => Some("NEUTRIM"),
            Self::MuonAntiNeutrino => Some("ANEUTRIM"),
            Self::Tauon => Some("TAU+"),
            Self::AntiTauon => Some("TAU-"),
            Self::TauNeutrino => Some("NEUTRIT"),
            Self::TauAntiNeutrino => Some("ANEUTRIT"),
            Self::Photon => Some("PHOTON"),
            Self::Pion => Some("PION-"),
            Self::AntiPion => Some("PION+"),
            Self::Kaon => Some("KAON+"),
            Self::AntiKaon => Some("KAON-"),
            Self::KaonZeroShort => Some("KAONSHRT"),
            Self::KaonZero => Some("KAONZERO"),
            Self::AntiKaonZero => Some("AKAONZER"),
            Self::Neutron => Some("NEUTRON"),
            Self::AntiNeutron => Some("ANEUTRON"),
            Self::Proton => Some("PROTON"),
            Self::AntiProton => Some("APROTON"),
            Self::Lambda => Some("LAMBDA"),
            Self::AntiLambda => Some("ALAMBDA"),
            Self::SigmaMinus => Some("SIGMA-"),
            Self::AntiSigmaMinus => Some("ASIGMA-"),
            Self::SigmaPlus => Some("SIGMA+"),
            Self::AntiSigmaPlus => Some("ASIGMA+"),
            Self::SigmaZero => Some("SIGMAZER"),
            Self::AntiSigmaZero => Some("ASIGMAZE"),
        }
    }

    /// GEANT4 particle name; `None` when unlisted upstream (heavy ions are
    /// handled separately by transport codes as `GenericIon`).
    #[must_use]
    pub const fn geant4(self) -> Option<&'static str> {
        match self {
            Self::Electron => Some("e-"),
            Self::Positron => Some("e+"),
            Self::ElectronNeutrino => Some("nu_e"),
            Self::ElectronAntiNeutrino => Some("anti_nu_e"),
            Self::Muon => Some("mu+"),
            Self::AntiMuon => Some("mu-"),
            Self::MuonNeutrino => Some("nu_mu"),
            Self::MuonAntiNeutrino => Some("anti_nu_mu"),
            Self::Tauon => Some("tau+"),
            Self::AntiTauon => Some("tau-"),
            Self::TauNeutrino => Some("nu_tau"),
            Self::TauAntiNeutrino => Some("anti_nu_tau"),
            Self::Photon => Some("gamma"),
            Self::Pion => Some("pi-"),
            Self::AntiPion => Some("pi+"),
            Self::Kaon => Some("kaon+"),
            Self::AntiKaon => Some("kaon-"),
            Self::KaonZeroShort => Some("kaon0S"),
            Self::KaonZero => Some("kaon0"),
            // No GEANT4 string upstream for AntiKaonZero.
            Self::AntiKaonZero => None,
            Self::Neutron => Some("neutron"),
            Self::AntiNeutron => Some("anti_neutron"),
            Self::Proton => Some("proton"),
            Self::AntiProton => Some("anti_proton"),
            Self::Lambda => Some("lambda"),
            Self::AntiLambda => Some("anti_lambda"),
            Self::SigmaMinus => Some("sigma-"),
            Self::AntiSigmaMinus => Some("anti_sigma-"),
            Self::SigmaPlus => Some("sigma+"),
            Self::AntiSigmaPlus => Some("anti_sigma+"),
            Self::SigmaZero => Some("sigma0"),
            Self::AntiSigmaZero => Some("anti_sigma0"),
        }
    }

    /// Parse a particle spec: canonical name, alternate name, PDC number,
    /// or an H-1 nuclide spec (`"H1"`, `"1H"`, `"10010000"`), which maps to
    /// [`ParticleId::Proton`] exactly. Heavy ions are rejected
    /// here — see [`is_heavy_ion`].
    pub fn parse(spec: &str) -> Result<Self, Error> {
        let s = spec.trim();
        if s.is_empty() {
            return Err(Error::NotAParticle(spec.to_string()));
        }
        match s {
            "Electron" => return Ok(Self::Electron),
            "Positron" => return Ok(Self::Positron),
            "ElectronNeutrino" => return Ok(Self::ElectronNeutrino),
            "ElectronAntiNeutrino" => return Ok(Self::ElectronAntiNeutrino),
            "Muon" => return Ok(Self::Muon),
            "AntiMuon" => return Ok(Self::AntiMuon),
            "MuonNeutrino" => return Ok(Self::MuonNeutrino),
            "MuonAntiNeutrino" => return Ok(Self::MuonAntiNeutrino),
            "Tauon" => return Ok(Self::Tauon),
            "AntiTauon" => return Ok(Self::AntiTauon),
            "TauNeutrino" => return Ok(Self::TauNeutrino),
            "TauAntiNeutrino" => return Ok(Self::TauAntiNeutrino),
            "Photon" => return Ok(Self::Photon),
            "Pion" => return Ok(Self::Pion),
            "AntiPion" => return Ok(Self::AntiPion),
            "Kaon" => return Ok(Self::Kaon),
            "AntiKaon" => return Ok(Self::AntiKaon),
            "KaonZeroShort" => return Ok(Self::KaonZeroShort),
            "KaonZero" => return Ok(Self::KaonZero),
            "AntiKaonZero" => return Ok(Self::AntiKaonZero),
            "Neutron" => return Ok(Self::Neutron),
            "AntiNeutron" => return Ok(Self::AntiNeutron),
            "Proton" => return Ok(Self::Proton),
            "AntiProton" => return Ok(Self::AntiProton),
            "Lambda" => return Ok(Self::Lambda),
            "AntiLambda" => return Ok(Self::AntiLambda),
            "Sigma-" => return Ok(Self::SigmaMinus),
            "AntiSigma-" => return Ok(Self::AntiSigmaMinus),
            "Sigma+" => return Ok(Self::SigmaPlus),
            "AntiSigma+" => return Ok(Self::AntiSigmaPlus),
            "Sigma" => return Ok(Self::SigmaZero),
            "AntiSigmaZero" => return Ok(Self::AntiSigmaZero),
            _ => {}
        }
        if let Some((_, p)) = ALIASES.iter().find(|(a, _)| *a == s) {
            return Ok(*p);
        }
        if let Some(p) = ALL
            .iter()
            .copied()
            .find(|p| p.name().eq_ignore_ascii_case(s))
        {
            return Ok(p);
        }
        if let Some((_, p)) = ALIASES.iter().find(|(a, _)| a.eq_ignore_ascii_case(s)) {
            return Ok(*p);
        }
        let (neg, digits) = match s.strip_prefix('-') {
            Some(d) => (true, d),
            None => (false, s),
        };
        if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
            if let Ok(v) = digits.parse::<i64>() {
                let signed = if neg { -v } else { v };
                if signed >= i32::MIN as i64 && signed <= i32::MAX as i64 {
                    if let Some(p) = Self::from_pdc(signed as i32) {
                        return Ok(p);
                    }
                }
                if !neg {
                    if let Some(nuc) = nucid_from_digits(digits) {
                        if is_ground_hydrogen(nuc) {
                            return Ok(Self::Proton);
                        }
                    }
                }
            }
            return Err(Error::NotAParticle(spec.to_string()));
        }
        if let Some(nuc) = nuclide_from_spec(s) {
            if is_ground_hydrogen(nuc) {
                return Ok(Self::Proton);
            }
            return Err(Error::NotAParticle(spec.to_string()));
        }
        Err(Error::NotAParticle(spec.to_string()))
    }
}

impl fmt::Display for ParticleId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl std::str::FromStr for ParticleId {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

/// Errors from particle-spec parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The spec is not a fundamental particle (heavy ions fall here; they
    /// are recognized by [`is_heavy_ion`] instead).
    NotAParticle(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NotAParticle(s) => write!(f, "not a valid particle name `{s}`"),
        }
    }
}

impl std::error::Error for Error {}

fn is_ground_hydrogen(nuc: NuclideId) -> bool {
    nuc.z() == 1 && nuc.a() == 1 && nuc.state() == 0
}

/// Nuclide behind a text spec, accepting both `"Na22"` (via
/// [`NuclideId::from_name`]) and reversed `"22Na"` forms.
fn nuclide_from_spec(s: &str) -> Option<NuclideId> {
    if let Ok(nuc) = NuclideId::from_name(s) {
        return Some(nuc);
    }
    let split = s.find(|c: char| c.is_ascii_alphabetic())?;
    if split == 0 || !s[..split].bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let (num, sym) = s.split_at(split);
    if sym.len() > 2 || !sym.bytes().all(|b| b.is_ascii_alphabetic()) {
        return None;
    }
    NuclideId::from_name(&format!("{sym}{num}")).ok()
}

/// Nuclide behind a pure-digit spec interpreted as a raw nucid.
fn nucid_from_digits(digits: &str) -> Option<NuclideId> {
    let v: u64 = digits.parse().ok()?;
    if v > u32::MAX as u64 {
        return None;
    }
    let nuc = NuclideId::from_nucid(v as u32);
    let (z, a, state) = (nuc.z(), nuc.a(), nuc.state());
    if (1..=118).contains(&z) && a >= z && state <= 9 {
        Some(nuc)
    } else {
        None
    }
}

fn spec_nuclide(s: &str) -> Option<NuclideId> {
    let nuc = nuclide_from_spec(s);
    if nuc.is_some() {
        return nuc;
    }
    if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) {
        return nucid_from_digits(s);
    }
    None
}

/// Is the spec a fundamental particle name/alias/PDC number, or a nuclide
/// (hydrogen or heavy ion)?
#[must_use]
pub fn is_valid(spec: &str) -> bool {
    ParticleId::parse(spec).is_ok() || is_heavy_ion(spec)
}

/// Is the integer a registered PDC number?
#[must_use]
pub fn is_valid_pdc(n: i32) -> bool {
    ParticleId::from_pdc(n).is_some()
}

/// Is the spec ground-state hydrogen — `"Proton"`, `"Hydrogen"`, `"Protium"`,
/// `"H1"`, `"1H"`, or nucid `10010000`?
#[must_use]
pub fn is_hydrogen(spec: &str) -> bool {
    let s = spec.trim();
    if ParticleId::parse(s) == Ok(ParticleId::Proton) {
        return true;
    }
    spec_nuclide(s).is_some_and(is_ground_hydrogen)
}

/// Is the spec a nuclide heavier than ground-state hydrogen? Heavy ions are
/// outside the PDC scheme (assigned PDC number 0).
#[must_use]
pub fn is_heavy_ion(spec: &str) -> bool {
    let s = spec.trim();
    spec_nuclide(s).is_some_and(|nuc| !is_ground_hydrogen(nuc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_canonical_names() {
        assert_eq!(ParticleId::parse("Neutron"), Ok(ParticleId::Neutron));
        assert_eq!(ParticleId::parse("Proton"), Ok(ParticleId::Proton));
        assert_eq!(ParticleId::parse("AntiProton"), Ok(ParticleId::AntiProton));
        assert_eq!(ParticleId::parse("Photon"), Ok(ParticleId::Photon));
        assert_eq!(ParticleId::parse("Sigma-"), Ok(ParticleId::SigmaMinus));
        assert_eq!(
            ParticleId::parse("ElectronNeutrino"),
            Ok(ParticleId::ElectronNeutrino)
        );
    }

    #[test]
    fn parses_alternate_names() {
        assert_eq!(ParticleId::parse("Hydrogen"), Ok(ParticleId::Proton));
        assert_eq!(ParticleId::parse("Protium"), Ok(ParticleId::Proton));
        assert_eq!(ParticleId::parse("Beta"), Ok(ParticleId::Electron));
        assert_eq!(ParticleId::parse("Beta-"), Ok(ParticleId::Electron));
        assert_eq!(ParticleId::parse("Beta+"), Ok(ParticleId::Positron));
        assert_eq!(ParticleId::parse("Gamma"), Ok(ParticleId::Photon));
        assert_eq!(ParticleId::parse("X-Ray"), Ok(ParticleId::Photon));
    }

    #[test]
    fn parses_case_insensitively() {
        assert_eq!(ParticleId::parse("neutron"), Ok(ParticleId::Neutron));
        assert_eq!(ParticleId::parse("NEUTRON"), Ok(ParticleId::Neutron));
        assert_eq!(ParticleId::parse("gamma"), Ok(ParticleId::Photon));
        assert_eq!(ParticleId::parse("beta-"), Ok(ParticleId::Electron));
        assert_eq!(ParticleId::parse("x-ray"), Ok(ParticleId::Photon));
    }

    #[test]
    fn parses_short_physics_symbols() {
        assert_eq!(ParticleId::parse("n"), Ok(ParticleId::Neutron));
        assert_eq!(ParticleId::parse("p"), Ok(ParticleId::Proton));
        assert_eq!(ParticleId::parse("e"), Ok(ParticleId::Electron));
        assert_eq!(ParticleId::parse("e-"), Ok(ParticleId::Electron));
        assert_eq!(ParticleId::parse("e+"), Ok(ParticleId::Positron));
        assert_eq!(ParticleId::parse("g"), Ok(ParticleId::Photon));
    }

    #[test]
    fn parses_numeric_specs() {
        assert_eq!(ParticleId::parse("2112"), Ok(ParticleId::Neutron));
        assert_eq!(ParticleId::parse("22"), Ok(ParticleId::Photon));
        assert_eq!(ParticleId::parse("-2212"), Ok(ParticleId::AntiProton));
        assert_eq!(ParticleId::parse("10010000"), Ok(ParticleId::Proton));
        assert!(matches!(
            ParticleId::parse("999999"),
            Err(Error::NotAParticle(_))
        ));
        assert!(matches!(
            ParticleId::parse("-999999"),
            Err(Error::NotAParticle(_))
        ));
    }

    #[test]
    fn rejects_non_particles_but_keeps_heavy_ions_detectable() {
        assert!(matches!(
            ParticleId::parse("Waka waka"),
            Err(Error::NotAParticle(_))
        ));
        assert!(matches!(
            ParticleId::parse("22Na"),
            Err(Error::NotAParticle(_))
        ));
        assert!(matches!(ParticleId::parse(""), Err(Error::NotAParticle(_))));
        assert!(is_heavy_ion("22Na"));
    }

    #[test]
    fn every_variant_round_trips_through_name_and_pdc() {
        for p in ALL {
            assert_eq!(ParticleId::parse(p.name()), Ok(p));
            assert_eq!(format!("{p}"), p.name());
            assert_eq!(ParticleId::from_pdc(p.pdc()), Some(p));
            let mut hit = 0;
            for q in ALL {
                if q.pdc() == p.pdc() {
                    hit += 1;
                }
            }
            assert_eq!(hit, 1, "duplicate PDC number on {}", p.name());
        }
        assert_eq!(ParticleId::from_pdc(0), None);
        assert_eq!(ALL.len(), 32);
    }

    #[test]
    fn mcnp_translations() {
        assert_eq!(ParticleId::Neutron.mcnp(), Some("n"));
        assert_eq!(ParticleId::Photon.mcnp(), Some("p"));
        assert_eq!(ParticleId::Electron.mcnp(), Some("e"));
        assert_eq!(ParticleId::Proton.mcnp(), None);
        assert_eq!(ParticleId::from_pdc(2112).unwrap().mcnp(), Some("n"));
    }

    #[test]
    fn mcnp6_translations() {
        assert_eq!(ParticleId::Neutron.mcnp6(), Some("n"));
        assert_eq!(ParticleId::Photon.mcnp6(), Some("p"));
        assert_eq!(ParticleId::Electron.mcnp6(), Some("e"));
        assert_eq!(ParticleId::Proton.mcnp6(), Some("h"));
        assert_eq!(ParticleId::parse("Hydrogen").unwrap().mcnp6(), Some("h"));
    }

    #[test]
    fn fluka_translations() {
        assert_eq!(ParticleId::Neutron.fluka(), Some("NEUTRON"));
        assert_eq!(ParticleId::Photon.fluka(), Some("PHOTON"));
        assert_eq!(ParticleId::Electron.fluka(), Some("ELECTRON"));
        assert_eq!(ParticleId::Proton.fluka(), Some("PROTON"));
        assert_eq!(
            ParticleId::parse("Beta-").unwrap().fluka(),
            Some("ELECTRON")
        );
        assert_eq!(
            ParticleId::parse("Hydrogen").unwrap().fluka(),
            Some("PROTON")
        );
        assert_eq!(ParticleId::Pion.fluka(), Some("PION-"));
        assert_eq!(ParticleId::AntiPion.fluka(), Some("PION+"));
        assert_eq!(ParticleId::AntiTauon.fluka(), Some("TAU-"));
        assert_eq!(ParticleId::SigmaZero.fluka(), Some("SIGMAZER"));
    }

    #[test]
    fn geant4_translations() {
        assert_eq!(ParticleId::Neutron.geant4(), Some("neutron"));
        assert_eq!(ParticleId::Photon.geant4(), Some("gamma"));
        assert_eq!(ParticleId::Electron.geant4(), Some("e-"));
        assert_eq!(ParticleId::Positron.geant4(), Some("e+"));
        assert_eq!(ParticleId::Proton.geant4(), Some("proton"));
        assert_eq!(ParticleId::parse("Beta-").unwrap().geant4(), Some("e-"));
        assert_eq!(
            ParticleId::parse("Hydrogen").unwrap().geant4(),
            Some("proton")
        );
        assert_eq!(ParticleId::KaonZeroShort.geant4(), Some("kaon0S"));
        assert_eq!(ParticleId::AntiSigmaZero.geant4(), Some("anti_sigma0"));
        assert_eq!(ParticleId::AntiKaonZero.geant4(), None);
    }

    #[test]
    fn validity_predicates() {
        for spec in [
            "Proton",
            "Protium",
            "Hydrogen",
            "Neutron",
            "AntiProton",
            "AntiNeutron",
            "H1",
        ] {
            assert!(is_valid(spec), "{spec}");
        }
        assert!(!is_valid("Waka waka"));
        assert!(is_valid_pdc(2212));
        assert!(is_valid_pdc(-2212));
        assert!(!is_valid_pdc(0));
        assert!(!is_valid_pdc(42));
    }

    #[test]
    fn hydrogen_and_heavy_ion_predicates() {
        for spec in ["Proton", "Hydrogen", "Protium", "1H", "H1", "10010000"] {
            assert!(!is_heavy_ion(spec), "{spec}");
        }
        for spec in ["2H", "H2", "3He", "He3", "22Na", "Na22"] {
            assert!(is_heavy_ion(spec), "{spec}");
        }
        for spec in ["Proton", "Hydrogen", "Protium", "1H", "H1", "10010000"] {
            assert!(is_hydrogen(spec), "{spec}");
        }
        assert!(!is_hydrogen("2H"));
        assert!(!is_hydrogen("Neutron"));
        assert!(!is_heavy_ion("Neutron"));
    }

    #[test]
    fn describes_particles() {
        assert_eq!(ParticleId::Photon.describe(), "Photon");
        assert_eq!(ParticleId::Proton.describe(), "Proton");
        assert_eq!(ParticleId::Neutron.describe(), "Neutron");
        assert_eq!(ParticleId::Electron.describe(), "Electron");
        assert_eq!(
            ParticleId::ElectronAntiNeutrino.describe(),
            "Electron Anti Neutrino"
        );
    }

    #[test]
    fn error_display() {
        assert_eq!(
            Error::NotAParticle("waka".into()).to_string(),
            "not a valid particle name `waka`"
        );
    }
}
