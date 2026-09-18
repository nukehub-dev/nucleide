//! Parametric tokamak plasma source: caller-supplied profiles over Miller
//! geometry, reactivity-weighted neutron emission, seeded sampling, and
//! histogram source-card emission.
//!
//! The source model is the public
//! ITER/EU-DEMO parametrization of Fausser et al., Fus. Eng. Des. **87**
//! (2012) 787 — [`MillerGeometry`] flux surfaces ([`crate::miller`]) with
//! L/H/A-mode density/temperature profiles ([`crate::profile`]) — and the
//! neutron emission strength
//!
//! ```text
//! S(r) = f_fuel · n_i(r)² · ⟨σv⟩(T_i(r))          [neutrons/s/cm³, relative]
//! ```
//!
//! with `⟨σv⟩` the Bosch & Hale fits ([`crate::reactivity`]) and
//! `f_fuel = 1/4` for equimolar D-T (`n_D = n_T = n_i/2`) or `1/2` for pure
//! D-D (the ½ avoids double-counting identical reactant pairs, and follows
//! the Fausser/`openmc-plasma-source` strength convention). Birth positions
//! are drawn ∝ `S(r)·R·|J|` — the volume element of [`crate::miller`] —
//! and birth energies from the local-ion-temperature Ballabio Gaussian of
//! [`crate::reaction`].
//!
//! # Fuel mixtures (Eriksson normalization)
//!
//! [`FuelMixture`] generalizes the two single-fuel factors to arbitrary D/T
//! fractions `f_D + f_T = 1` at the shared ion temperature of the profiles
//! (mixture v1). With total ion density `n`,
//!
//! ```text
//! S = n² · [f_D·f_T·⟨σv⟩_DT(T_i) + (f_D²/2)·⟨σv⟩_DD(T_i)]
//! ```
//!
//! — the standard `n_a·n_b·⟨σv⟩` rate with the `1/(1+δ_ab)` same-species
//! guard, the DRESS rate construction (Eriksson et al., Comput. Phys.
//! Commun. **199** (2016) 40, and the mixture form upstream
//! `openmc-plasma-source` implements): the D-T term pairs every D with
//! every T (`f_D·f_T`), while the D-D term halves the `f_D²` self-pair
//! count so identical reactants are not double-counted. Both branches keep
//! the landed neutron-branch convention — the T-T branch and the D(d,p)T
//! proton branch stay out of scope.
//!
//! Hand vectors (flat profile, `T_i = 20` keV, `n = 1e20` m⁻³; the pinned
//! reactivity goldens `⟨σv⟩_DT = 4.330220397689175e-22` and
//! `⟨σv⟩_DD = 2.602582958721524e-24` m³/s):
//!
//! ```text
//! f_D = f_T = 1/2:      S = 1082555099422.2937 + 3253228698.4019046
//!                       = 1085808328120.6956   [neutrons/s/cm³, relative]
//! f_D = 1 (f_T = 0):    S = 13012914793.607618
//! f_D = 0.7, f_T = 0.3: S = 909346283514.7268 + 6376328248.867733
//!                       = 915722611763.5945
//! ```
//!
//! Exact recovery anchors (regression gates, not approximations): the D-T
//! branch term at `f_D = f_T = 1/2` is bit-for-bit the landed equimolar
//! kernel (`0.5·0.5 == 0.25` exactly, same expression tree), and the total
//! at `f_D = 1` is bit-for-bit the landed pure D-D kernel. Configurations
//! with `fuel_mixture: None` keep the landed expression trees untouched.
//!
//! # Fusion branches: T-T neutrons and D(d,p)T protons (v1 stance)
//!
//! Two branches sit behind the landed two-neutron-branch kernel:
//!
//! ## T-T neutrons: normalization pinned, transport deferred
//!
//! The T-T neutron term enters the mixture sum as a third term with neutron
//! multiplicity 2 folded in — the same rule the `openmc-plasma-source`
//! oracle implements (`tokamak_source`: T-T fuel density `n_T²/2`, neutron
//! source density doubled since T(t,2n)⁴He releases two neutrons):
//!
//! ```text
//! S = n² · [f_D·f_T·⟨σv⟩_DT + (f_D²/2)·⟨σv⟩_DD + f_T²·⟨σv⟩_TT]
//!                              [neutrons/s/cm³, relative]
//! ```
//!
//! reacting at `T_T` (both reactants are tritium, so — like the D-D arm at
//! `T_D` — no pair effective temperature applies). The neutron coefficient
//! is [`FuelMixture::tt_neutron_coefficient`] (`f_T²`): exactly `0.0` at
//! `f_T = 0`, `1.0` for pure tritium.
//!
//! v1 does NOT transport T-T neutrons through the sampler or the cards,
//! because no publishable closed form exists for either half of the branch:
//! Bosch & Hale, Nucl. Fusion **32** (1992) 611 covers only D(d,n)³He,
//! D(d,p)T, T(d,n)⁴He, and ³He(d,p)⁴He (no T-T reactivity fit), and Ballabio
//! et al., Nucl. Fusion **38** (1998) 1723 Table III covers the two-body
//! D-D/D-T lines (the T-T three-body breakup continuum has no Gaussian
//! line). The oracle carries T-T as vendored tables instead (NeSST `reac_TT`
//! Hale spline plus Brune/Eriksson/Gatu-Johnson continuum files), and
//! vendored plasma data stays out of scope — a Gaussian T-T line would be a
//! guess, never shipped. The deferred branch is quantified, not ignored:
//! the validation oracle compares the landed two-branch kernel against the
//! three-branch upstream model on a tritium-rich blend (container-only),
//! and a pure-tritium mixture stays loud (zero neutron strength).
//!
//! Zero-tritium recovery anchor (regression gate): a mixture with `f_T = 0`
//! draws exactly `0.0` from the absent branches, so the no-T mixture
//! reproduces the landed two-branch kernel bit-for-bit.
//!
//! ## D(d,p)T protons: accounting output, sampler untouched
//!
//! The D(d,p)T proton branch shares the D-D total rate 50/50 with the
//! neutron branch (Bosch & Hale resolve the two D-D branches to ~4–8% at
//! 10–20 keV; the pinned v1 convention reports them equal —
//! [`ParametricPlasmaConfig::proton_strength_density`] IS the D-D neutron
//! branch density, bit-for-bit). Protons are bookkeeping alongside the
//! neutron source, never particles: the sampler draws no proton energies,
//! the cards carry no proton distributions, and proton transport stays out
//! of scope. [`ParametricPlasmaConfig::total_proton_strength`] integrates
//! the proton density over the same `S·R·|J|` volume (sector-aware), so the
//! `proton_per_neutron` bookkeeping ratio is exactly `1.0` for pure D-D.
//!
//! Hand vectors (flat profile, `T_i = 20` keV, `n = 1e20` m⁻³; the pinned
//! `⟨σv⟩_DD = 2.602582958721524e-24` m³/s):
//!
//! ```text
//! pure D-D:             P = 13012914793.607618 = S_DD  →  P/N = 1.0
//! f_D = 0.7, f_T = 0.3: P = 6376328248.867733 (the landed D-D hand vector)
//! f_D = f_T = 1/2:      P = 3253228698.4019046
//! T-T neutron coeff:    f_T² = 0.25 / 0.0 / 0.09 / 1.0
//!                       (equimolar / no-T / 70-30 / pure-T)
//! ```
//!
//! # Per-species ion temperatures (mixture v2)
//!
//! [`SpeciesIonTemperatures`] generalizes the mixture rule to distinct
//! Maxwellian species temperatures `(T_D, T_T)` in keV — uniform scalars,
//! flat in minor radius. Each pair reacts at its own relative-temperature
//! closed form: two Maxwellians at `(T_D, T_T)` have a Maxwellian
//! relative-velocity distribution at the mass-weighted effective temperature
//! (Eriksson et al., Comput. Phys. Commun. **199** (2016) 40 — the full
//! Eriksson generalization for distinct Maxwellians):
//!
//! ```text
//! S = n² · [f_D·f_T·⟨σv⟩_DT(T_DT) + (f_D²/2)·⟨σv⟩_DD(T_D)]
//! T_DT = (m_D·T_T + m_T·T_D) / (m_D + m_T),  m_D = 2, m_T = 3 (mass numbers)
//!      = T_D + (2/5)·(T_T − T_D)   (this spelling, so T_D = T_T recovers T_D
//!                                   bit-for-bit: the difference is exactly 0)
//! ```
//!
//! The D-D branch is exact (both reactants are deuterium at `T_D`); the D-T
//! branch evaluates the landed Bosch–Hale fit at `T_DT`. Per-branch Ballabio
//! spectra follow the same pair temperatures — the D-T line at `T_DT`, the
//! D-D line at `T_D` — so the sampler's branch roulette, the card energy
//! marginals, and the axis summary all stay consistent with the strength
//! rule. The D-T spectrum arm is an effective-temperature convention:
//! Ballabio Table III was fitted to single-temperature Maxwellians, so away
//! from `T_D = T_T` the Gaussian-at-`T_DT` line captures the leading
//! mass-weighted relative-temperature dependence while the residual against
//! the true distinct-temperature spectrum (the full Eriksson numerical
//! integration over arbitrary reactant distributions) is out of scope to
//! quantify — non-Maxwellian reactants stay a loud boundary.
//!
//! Hand vectors (flat profile, `n = 1e20` m⁻³, `f_D = 0.7`, `f_T = 0.3`,
//! `T_D = 20` keV, `T_T = 30` keV; `T_DT = 24.0` keV exactly, pinned
//! reactivity `⟨σv⟩_DT(24) = 5.414667327922193e-22` m³/s alongside the
//! landed `⟨σv⟩_DD(20) = 2.602582958721524e-24` m³/s):
//!
//! ```text
//! S = 1137080138863.6606 + 6376328248.867733 = 1143456467112.5283
//! ```
//!
//! Exact recovery anchor (regression gate): `T_D = T_T` reproduces the
//! shared-temperature mixture kernel bit-for-bit — strength densities and
//! the axis spectrum summary — because `T_DT` computes to `T_D` exactly.
//! The pair requires a fuel mixture (a single-fuel config with the pair set
//! is a loud [`Error::NotYetSupported`]); the profile ion temperature is
//! then unused for rate and spectrum (still validated), while the density
//! profile keeps shaping `S(r) ∝ n(r)²`.
//!
//! # Deuterium hot tail (mixture v3)
//!
//! [`DeuteriumTail`] pins one parametrized non-Maxwellian reactant shape — a
//! single-tail-temperature deuterium fraction, the auxiliary-heated
//! (NBI/minority-ICRH) tail approximated as a hot Maxwellian sub-population
//! — within the Eriksson et al., Comput. Phys. Commun. **199** (2016) 40
//! arbitrary-distribution framework. The deuterium distribution splits into a
//! bulk `(1 − η)` at `T_D` and a tail `η` at `T_tail` (uniform scalars, flat
//! in minor radius); tritium stays Maxwellian at `T_T`. Each sub-population
//! pair reacts at its own mass-weighted relative-temperature closed form
//! (the mixture-v2 rule applied per pair):
//!
//! ```text
//! S = n² · [f_D·f_T·((1−η)·⟨σv⟩_DT(T_DT) + η·⟨σv⟩_DT(T_DTt))
//!          + (f_D²/2)·((1−η)²·⟨σv⟩_DD(T_D) + 2η(1−η)·⟨σv⟩_DD(T_mix)
//!                       + η²·⟨σv⟩_DD(T_tail))]
//! T_DT   = T_D + (2/5)·(T_T − T_D)         (bulk D–T, the v2 spelling)
//! T_DTt  = T_tail + (2/5)·(T_T − T_tail)   (tail D–T)
//! T_mix  = (T_D + T_tail)/2                (bulk–tail D–D cross pair:
//!                                          equal masses, so the relative
//!                                          temperature is the plain mean)
//! ```
//!
//! Per-branch Ballabio spectra follow the same sub-pair temperatures — the
//! D-T arm carries two lines (at `T_DT` and `T_DTt`), the D-D arm three (at
//! `T_D`, `T_mix`, `T_tail`) — weighted by the sub-rate rule, so the
//! sampler's branch roulette, the card energy marginals, and the axis summary
//! all stay consistent with the strength rule. The tail-line arm is an
//! effective-temperature convention with the same error stance as mixture v2:
//! Ballabio Table III was fitted to single-temperature Maxwellians, so the
//! Gaussian-at-effective-`T` lines capture the leading mass-weighted
//! relative-temperature dependence while the residual against the true
//! distinct-temperature spectrum (the full Eriksson numerical integration
//! over arbitrary reactant distributions) is out of scope to quantify.
//!
//! Hand vectors (flat profile, `n = 1e20` m⁻³, `f_D = 0.7`, `f_T = 0.3`,
//! `T_D = 20` keV, `T_T = 30` keV, `η = 0.05`, `T_tail = 60` keV;
//! `T_DT = 24.0`, `T_DTt = 48.0`, `T_mix = 40.0` keV exactly, pinned
//! reactivities `⟨σv⟩_DT(48) = 8.557148874535018e-22`,
//! `⟨σv⟩_DD(40) = 8.234915801644433e-24`,
//! `⟨σv⟩_DD(60) = 1.4470008265609934e-23` m³/s alongside the landed
//! `⟨σv⟩_DT(24) = 5.414667327922193e-22` and
//! `⟨σv⟩_DD(20) = 2.602582958721524e-24` m³/s):
//!
//! ```text
//! S = 1170076195103.095 + 7759941698.06273 = 1177836136801.1577
//! ```
//!
//! Exact recovery anchor (regression gate): `η = 0` reproduces the
//! no-tail mixture kernel bit-for-bit — strength densities, sampled stream,
//! and emitted cards — because `(1 − 0) = 1` and every `η`-scaled term is
//! exactly `0.0` (the `η == 0` path additionally skips the tail reactivity
//! evaluations, so any in-range `T_tail` spelling recovers). The tail
//! requires a fuel mixture (a single-fuel config with the tail set is a loud
//! [`Error::NotYetSupported`]); it composes with [`SpeciesIonTemperatures`]
//! (`T_D`/`T_T` from the pair when set, else the shared profile temperature)
//! and with [`ToroidalSector`] (the angle draw is untouched — the tail adds
//! no new RNG draws on the mixture path beyond the landed one-per-particle
//! branch roulette).
//!
//! # Toroidal sectors
//!
//! [`ToroidalSector`] restricts birth positions to the toroidal interval
//! `[start_angle, start_angle + rotation_angle)` (radians, measured from the
//! machine `+x` axis toward `+y`, the sampler's azimuth convention). The
//! landed kernel `S(r)·R·|J|` is independent of the toroidal angle, so the
//! toroidal integral contributes `rotation_angle` instead of `2π`:
//!
//! ```text
//! total(sector) = total(full torus) · rotation_angle / 2π
//! ```
//!
//! and birth angles are uniform over the sector. The radial, poloidal, and
//! spectral physics are untouched — a sector only rescales the total and
//! remaps the one toroidal-angle draw — so `(r, θ, E)` moments match the
//! full-torus quadrature and only the angle marginal changes.
//!
//! Hand vectors (flat L-mode profile, `T_i = 20` keV, `n = 1e20` m⁻³; the
//! pinned D-T branch strength `1082555099422.2937` from the mixture vectors
//! above, equimolar D-T config):
//!
//! ```text
//! rotation = π (half torus):     S(r) unchanged; total = full · 0.5 exactly
//! rotation = π/2 (quarter torus): total = full · 0.25 exactly
//! ```
//!
//! Exact recovery anchor (regression gate, not an approximation): a sector
//! with `start_angle = 0` and `rotation_angle = 2π` (exactly) reproduces the
//! landed full-torus kernel bit-for-bit — the same `2.0·π` constant scales
//! the total, `0.0 + 2π·u == 2π·u` maps the same toroidal draw, no new RNG
//! draws exist on any path, and card emission omits the angle marginal, so
//! the sampled stream and the emitted cards are identical. A full-rotation
//! sector with nonzero start is distribution-identical (not bit-for-bit:
//! the angle origin shifts).
//!
//! Card emission carries a partial sector as a fourth uniform angle-bin
//! marginal (`PHI=D4` on SDEF, `phi d4` on Serpent, bin centers over
//! `[start, start + rotation)` radians) plus a `toroidal sector` drift row;
//! a full rotation (or no sector) emits the landed three-marginal card.
//!
//! Sector angles are validated loudly at construction: non-finite angles
//! give [`Error::NonFinite`]; `start_angle` outside `[0, 2π)` or
//! `rotation_angle` outside `(0, 2π]` gives [`Error::InvalidSector`].
//! Angles are never renormalized or wrapped — pass a normalized sector.
//!
//! Profiles are caller inputs; nothing here computes profiles or solves an
//! equilibrium. The documented loud boundary ([`Error::NotYetSupported`]):
//! reactant distributions beyond one deuterium hot-tail fraction
//! ([`DeuteriumTail`] pins the single v1 tail shape — the rest of the full
//! Eriksson generalization stays out),
//! T-T neutron transport (normalization pinned above; no publishable fit or
//! line exists), and proton transport (the D(d,p)T proton *rate* is
//! accounted via [`ParametricPlasmaConfig::proton_strength_density`] /
//! [`ParametricPlasmaConfig::total_proton_strength`]; no proton particles
//! are sampled and no proton distributions reach the cards).
//! Mixture fractions, species temperatures, tail parameters, and sector angles
//! themselves are validated loudly at construction: non-finite, negative,
//! non-summing, or out-of-range values never reach the sampler.
//!
//! # Card emission
//!
//! Transport source cards cannot represent the correlated `(r, z)` joint
//! (nor the position–energy correlation), so [`emission_histograms`] renders
//! the *marginals*: radial and vertical histograms plus the global marginal
//! energy spectrum. The drift report quantifies what the card preserves and
//! what it loses (truncation; the joint-correlation distance).

use std::f64::consts::PI;

use crate::miller::MillerGeometry;
use crate::profile::{self, DensityProfile, ProfileMode, TemperatureProfile};
use crate::sample::{Particle, Rng};
use crate::{Error, FusionReaction, Result};

/// Fine-grid resolution for the radial weight table.
const R_GRID: usize = 256;
/// Poloidal quadrature points per radial cell.
const THETA_GRID: usize = 256;

/// D/T fuel mixture for the parametric source: atom fractions of the total
/// ion density with `f_D + f_T = 1`, reacting at the shared ion temperature
/// of the profiles (mixture v1 — see the module rustdoc for the pinned
/// normalization, hand vectors, and recovery anchors).
///
/// Construction is loud: NaN/infinite fractions give [`Error::NonFinite`],
/// negative or non-summing pairs give [`Error::InvalidFuelMixture`]. The
/// fractions are used exactly as given — nothing is renormalized.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FuelMixture {
    /// Deuterium fraction `f_D` of the total ion density.
    pub f_deuterium: f64,
    /// Tritium fraction `f_T` of the total ion density.
    pub f_tritium: f64,
}

/// Absolute tolerance on the `f_D + f_T = 1` precondition. It admits
/// last-ulp missums (e.g. `1/3 + 2/3` in f64) and never renormalizes — the
/// rate rule plugs the fractions in as given.
const FRACTION_SUM_TOL: f64 = 1e-12;

impl FuelMixture {
    /// New validated mixture; errors are loud and named.
    pub fn new(f_deuterium: f64, f_tritium: f64) -> Result<Self> {
        let mixture = Self {
            f_deuterium,
            f_tritium,
        };
        mixture.validate()?;
        Ok(mixture)
    }

    /// Validate both fractions: finite, non-negative, and summing to 1
    /// within [`FRACTION_SUM_TOL`].
    pub fn validate(&self) -> Result<()> {
        if !self.f_deuterium.is_finite() {
            return Err(Error::NonFinite("deuterium fraction"));
        }
        if !self.f_tritium.is_finite() {
            return Err(Error::NonFinite("tritium fraction"));
        }
        if self.f_deuterium < 0.0 || self.f_tritium < 0.0 {
            return Err(Error::InvalidFuelMixture("fuel fractions must be >= 0"));
        }
        if (self.f_deuterium + self.f_tritium - 1.0).abs() > FRACTION_SUM_TOL {
            return Err(Error::InvalidFuelMixture("fuel fractions must sum to 1"));
        }
        Ok(())
    }

    /// Branch weights of the mixture rate rule (module rustdoc): the D-T
    /// pair weight `f_D·f_T` and the same-species-guarded D-D weight
    /// `f_D²/2`. At the recovery anchors these are exact — `(0.5, 0.5)` →
    /// `(0.25, 0.125)`, `(1, 0)` → `(0.0, 0.5)` — bit-for-bit the landed
    /// single-fuel factors.
    pub fn branch_weights(&self) -> (f64, f64) {
        (
            self.f_deuterium * self.f_tritium,
            self.f_deuterium * self.f_deuterium / 2.0,
        )
    }

    /// T-T neutron coefficient of the pinned three-branch normalization
    /// (module rustdoc): `f_T²` — the `f_T²/2` same-species-guarded reaction
    /// weight with the T(t,2n)⁴He neutron multiplicity 2 folded in, matching
    /// the oracle's `½·n_T²·reac_TT·2` rule. Exactly `0.0` at `f_T = 0` (the
    /// zero-tritium recovery anchor) and `1.0` for pure tritium. v1 carries
    /// the coefficient only — no `⟨σv⟩_TT` fit exists in publishable closed
    /// form (module rustdoc), so T-T neutron transport stays a loud
    /// boundary and this method is the machine-readable pin for the future
    /// branch term.
    pub fn tt_neutron_coefficient(&self) -> f64 {
        self.f_tritium * self.f_tritium
    }
}

/// Toroidal sector for the parametric source: birth positions are uniform
/// in toroidal angle over `[start_angle, start_angle + rotation_angle)`
/// (radians; see the module rustdoc for the pinned normalization, hand
/// vectors, and recovery anchor).
///
/// Construction is loud: NaN/infinite angles give [`Error::NonFinite`],
/// `start_angle` outside `[0, 2π)` or `rotation_angle` outside `(0, 2π]`
/// gives [`Error::InvalidSector`]. Angles are used exactly as given —
/// nothing is renormalized or wrapped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToroidalSector {
    /// Sector start angle `φ₀` \[rad\], in `[0, 2π)`.
    pub start_angle: f64,
    /// Sector extent \[rad\], in `(0, 2π]` (`2π` exactly is the full torus).
    pub rotation_angle: f64,
}

impl ToroidalSector {
    /// New validated sector; errors are loud and named.
    pub fn new(start_angle: f64, rotation_angle: f64) -> Result<Self> {
        let sector = Self {
            start_angle,
            rotation_angle,
        };
        sector.validate()?;
        Ok(sector)
    }

    /// Validate both angles: finite, `start_angle` in `[0, 2π)`,
    /// `rotation_angle` in `(0, 2π]`.
    pub fn validate(&self) -> Result<()> {
        if !self.start_angle.is_finite() {
            return Err(Error::NonFinite("toroidal sector start angle"));
        }
        if !self.rotation_angle.is_finite() {
            return Err(Error::NonFinite("toroidal sector rotation angle"));
        }
        if !(0.0..2.0 * PI).contains(&self.start_angle) {
            return Err(Error::InvalidSector("start_angle must be in [0, 2π)"));
        }
        if self.rotation_angle <= 0.0 || self.rotation_angle > 2.0 * PI {
            return Err(Error::InvalidSector("rotation_angle must be in (0, 2π]"));
        }
        Ok(())
    }

    /// Fraction of the full torus this sector covers
    /// (`rotation_angle / 2π`): the total-strength scale factor.
    pub fn fraction(&self) -> f64 {
        self.rotation_angle / (2.0 * PI)
    }

    /// True for the exact full-rotation spelling (`rotation_angle == 2π`):
    /// card emission omits the angle marginal and the kernel recovers the
    /// landed full-torus form (module rustdoc anchor).
    pub fn is_full_rotation(&self) -> bool {
        self.rotation_angle == 2.0 * PI
    }

    /// Map a unit uniform draw to a birth toroidal angle \[rad\].
    pub fn map_angle(&self, u: f64) -> f64 {
        self.start_angle + self.rotation_angle * u
    }
}

/// Per-species ion temperatures for the D/T fuel mixture: uniform
/// deuterium/tritium temperatures in keV, flat in minor radius (mixture v2 —
/// see the module rustdoc for the pinned normalization, hand vectors, and
/// recovery anchor).
///
/// Construction is loud: NaN/infinite temperatures give [`Error::NonFinite`],
/// negative temperatures give [`Error::NegativeIonTemperature`]. The pair
/// requires [`ParametricPlasmaConfig::fuel_mixture`] (a single-fuel config
/// with the pair set is [`Error::NotYetSupported`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpeciesIonTemperatures {
    /// Deuterium ion temperature `T_D` \[keV\]; non-negative.
    pub deuterium_kev: f64,
    /// Tritium ion temperature `T_T` \[keV\]; non-negative.
    pub tritium_kev: f64,
}

impl SpeciesIonTemperatures {
    /// D–T mass weight `m_D/(m_D + m_T)` with mass numbers `m_D = 2`,
    /// `m_T = 3` (exact integers, documented convention).
    const DEUTERIUM_WEIGHT: f64 = 2.0 / 5.0;

    /// New validated pair; errors are loud and named.
    pub fn new(deuterium_kev: f64, tritium_kev: f64) -> Result<Self> {
        let pair = Self {
            deuterium_kev,
            tritium_kev,
        };
        pair.validate()?;
        Ok(pair)
    }

    /// Validate both temperatures: finite and non-negative.
    pub fn validate(&self) -> Result<()> {
        if !self.deuterium_kev.is_finite() {
            return Err(Error::NonFinite("deuterium ion temperature"));
        }
        if !self.tritium_kev.is_finite() {
            return Err(Error::NonFinite("tritium ion temperature"));
        }
        if self.deuterium_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(self.deuterium_kev));
        }
        if self.tritium_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(self.tritium_kev));
        }
        Ok(())
    }

    /// D–T effective temperature \[keV\] for a species pair: the
    /// mass-weighted relative temperature
    /// `(m_D·T_T + m_T·T_D)/(m_D + m_T)` spelled as
    /// `T_D + (2/5)·(T_T − T_D)`, so equal temperatures recover `T_D`
    /// bit-for-bit (the difference is exactly `0.0`; module rustdoc anchor).
    pub fn dt_effective_kev(deuterium_kev: f64, tritium_kev: f64) -> f64 {
        deuterium_kev + Self::DEUTERIUM_WEIGHT * (tritium_kev - deuterium_kev)
    }

    /// D–T effective temperature \[keV\] of this pair (see
    /// [`Self::dt_effective_kev`]).
    pub fn dt_effective(&self) -> f64 {
        Self::dt_effective_kev(self.deuterium_kev, self.tritium_kev)
    }
}

/// Deuterium hot-tail fraction for the D/T fuel mixture: a uniform tail
/// population `η` at `T_tail` \[keV\], flat in minor radius, on top of the
/// bulk deuterium at `T_D` (mixture v3 — see the module rustdoc for the
/// pinned normalization, hand vectors, and recovery anchor).
///
/// Construction is loud: NaN/infinite parameters give [`Error::NonFinite`], a
/// fraction outside `[0, 1]` gives [`Error::InvalidTail`], a negative tail
/// temperature gives [`Error::NegativeIonTemperature`]. The fraction is used
/// exactly as given — nothing is renormalized. The tail requires
/// [`ParametricPlasmaConfig::fuel_mixture`] (a single-fuel config with the
/// tail set is [`Error::NotYetSupported`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeuteriumTail {
    /// Tail fraction `η` of the deuterium population in `[0, 1]` (`0`
    /// recovers the no-tail kernel bit-for-bit).
    pub fraction: f64,
    /// Tail temperature `T_tail` \[keV\]; non-negative.
    pub temperature_kev: f64,
}

impl DeuteriumTail {
    /// New validated tail; errors are loud and named.
    pub fn new(fraction: f64, temperature_kev: f64) -> Result<Self> {
        let tail = Self {
            fraction,
            temperature_kev,
        };
        tail.validate()?;
        Ok(tail)
    }

    /// Validate both parameters: finite, fraction in `[0, 1]`, non-negative
    /// tail temperature.
    pub fn validate(&self) -> Result<()> {
        if !self.fraction.is_finite() {
            return Err(Error::NonFinite("deuterium tail fraction"));
        }
        if !self.temperature_kev.is_finite() {
            return Err(Error::NonFinite("deuterium tail temperature"));
        }
        if !(0.0..=1.0).contains(&self.fraction) {
            return Err(Error::InvalidTail(
                "deuterium tail fraction must be in [0, 1]",
            ));
        }
        if self.temperature_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(self.temperature_kev));
        }
        Ok(())
    }
}

/// Parametric tokamak plasma source configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParametricPlasmaConfig {
    /// Miller flux-surface geometry.
    pub geometry: MillerGeometry,
    /// Confinement-mode profile family (L/H/A).
    pub mode: ProfileMode,
    /// Ion-density profile parameters (Fausser convention, m⁻³).
    pub ion_density: DensityProfile,
    /// Ion-temperature profile parameters (keV).
    pub ion_temperature: TemperatureProfile,
    /// Pedestal radius `r_ped` \[cm\]; required in `(0, a_minor)` for H/A.
    pub pedestal_radius_cm: f64,
    /// Fuel reaction: D-T (equimolar) or D-D. Used only when
    /// [`Self::fuel_mixture`] is `None`; a mixture activates both branches
    /// through [`FuelMixture::branch_weights`] and this field plays no
    /// role in the strength or spectrum model.
    pub fuel: FusionReaction,
    /// Optional D/T fuel mixture (module rustdoc normalization). `None` —
    /// the default — keeps the landed single-fuel kernels bit-for-bit.
    pub fuel_mixture: Option<FuelMixture>,
    /// Optional per-species ion temperatures (module rustdoc normalization).
    /// `None` — the default — reacts at the shared profile ion temperature.
    /// When set, a fuel mixture is required and the profile ion temperature
    /// is unused for rate and spectrum (still validated).
    pub species_temperatures: Option<SpeciesIonTemperatures>,
    /// Optional deuterium hot-tail fraction (module rustdoc normalization).
    /// `None` — the default — keeps the landed mixture kernels bit-for-bit.
    /// When set, a fuel mixture is required; it composes with
    /// [`Self::species_temperatures`] (`T_D`/`T_T` from the pair when set,
    /// else the shared profile temperature).
    pub tail: Option<DeuteriumTail>,
    /// Optional toroidal sector (module rustdoc normalization). `None` —
    /// the default — keeps the landed full-torus kernel bit-for-bit.
    pub sector: Option<ToroidalSector>,
    /// Particle weight carried by the sampler and the emitted cards.
    pub weight: f64,
}

impl ParametricPlasmaConfig {
    /// Validate the full configuration.
    pub fn validate(&self) -> Result<()> {
        self.geometry.validate()?;
        profile::validate_profiles(
            self.mode,
            &self.ion_density,
            &self.ion_temperature,
            self.geometry.minor_radius_cm,
            self.pedestal_radius_cm,
        )?;
        if !self.pedestal_radius_cm.is_finite() {
            return Err(Error::InvalidProfile("pedestal radius"));
        }
        if !self.weight.is_finite() || self.weight <= 0.0 {
            return Err(Error::NonPositiveWeight(self.weight));
        }
        if let Some(mixture) = &self.fuel_mixture {
            mixture.validate()?;
        }
        if let Some(pair) = &self.species_temperatures {
            pair.validate()?;
            if self.fuel_mixture.is_none() {
                return Err(Error::NotYetSupported(
                    "per-species ion temperatures need a D/T fuel mixture",
                ));
            }
        }
        if let Some(tail) = &self.tail {
            tail.validate()?;
            if self.fuel_mixture.is_none() {
                return Err(Error::NotYetSupported(
                    "a deuterium hot tail needs a D/T fuel mixture",
                ));
            }
        }
        if let Some(sector) = &self.sector {
            sector.validate()?;
        }
        Ok(())
    }

    /// Ion density \[m⁻³\] at minor radius `r` \[cm\].
    pub fn density_m3(&self, r: f64) -> f64 {
        profile::density(
            self.mode,
            &self.ion_density,
            self.geometry.minor_radius_cm,
            self.pedestal_radius_cm,
            r,
        )
    }

    /// Ion temperature \[keV\] at minor radius `r` \[cm\].
    pub fn temperature_kev(&self, r: f64) -> f64 {
        profile::temperature(
            self.mode,
            &self.ion_temperature,
            self.geometry.minor_radius_cm,
            self.pedestal_radius_cm,
            r,
        )
    }

    /// Species ion temperatures `(T_D, T_T)` \[keV\] at minor radius `r`
    /// \[cm\]: the uniform [`SpeciesIonTemperatures`] pair when set, else
    /// the shared profile ion temperature twice.
    fn species_temperatures_at(&self, r: f64) -> (f64, f64) {
        match &self.species_temperatures {
            Some(pair) => (pair.deuterium_kev, pair.tritium_kev),
            None => {
                let ti_kev = self.temperature_kev(r);
                (ti_kev, ti_kev)
            }
        }
    }

    /// (D-T, D-D) branch strength densities \[neutrons/s/cm³, arbitrary
    /// global scale\] at minor radius `r` \[cm\]: each branch `k·n²·⟨σv⟩`
    /// with density in cm⁻³ and reactivity in cm³/s, built in the landed
    /// left-associated expression tree so that single-fuel configurations
    /// and the module-rustdoc recovery anchors reproduce bit-for-bit. A
    /// mixture reacts the D-T branch at the pair effective temperature
    /// [`SpeciesIonTemperatures::dt_effective_kev`] and the D-D branch at
    /// `T_D` (module rustdoc); without a pair both reduce to the shared
    /// profile temperature exactly. A [`DeuteriumTail`] splits the deuterium
    /// population into bulk/tail sub-pairs at their own effective
    /// temperatures (module rustdoc); `η == 0` takes the landed mixture
    /// expressions verbatim (bit-for-bit recovery for any `T_tail`
    /// spelling). Zero
    /// temperature ⇒ zero strength on both branches (cold separatrix makes
    /// no neutrons).
    fn branch_strength_density(&self, r: f64) -> Result<(f64, f64)> {
        let n_m3 = self.density_m3(r);
        let (t_d, t_t) = self.species_temperatures_at(r);
        let n_cm3 = n_m3 * 1e-6;
        match &self.fuel_mixture {
            Some(mixture) => {
                let (k_dt, k_dd) = mixture.branch_weights();
                let t_eff = SpeciesIonTemperatures::dt_effective_kev(t_d, t_t);
                let tail_eta = self.tail.map(|tail| tail.fraction).unwrap_or(0.0);
                if tail_eta == 0.0 {
                    let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(t_eff)? * 1e6;
                    let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(t_d)? * 1e6;
                    return Ok((k_dt * n_cm3 * n_cm3 * sv_dt, k_dd * n_cm3 * n_cm3 * sv_dd));
                }
                let tail = self.tail.expect("tail fraction nonzero implies a tail");
                let eta = tail.fraction;
                let t_tail = tail.temperature_kev;
                let t_eff_tail = SpeciesIonTemperatures::dt_effective_kev(t_tail, t_t);
                let t_mix = 0.5 * (t_d + t_tail);
                let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(t_eff)? * 1e6;
                let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(t_eff_tail)? * 1e6;
                let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(t_d)? * 1e6;
                let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(t_mix)? * 1e6;
                let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(t_tail)? * 1e6;
                let s_dt = k_dt * n_cm3 * n_cm3 * ((1.0 - eta) * sv_dt + eta * sv_dt_tail);
                let s_dd = k_dd
                    * n_cm3
                    * n_cm3
                    * ((1.0 - eta) * (1.0 - eta) * sv_dd
                        + 2.0 * eta * (1.0 - eta) * sv_dd_mix
                        + eta * eta * sv_dd_tail);
                Ok((s_dt, s_dd))
            }
            None => {
                let reactivity = self.fuel.reactivity_m3_per_s(t_d)? * 1e6;
                match self.fuel {
                    FusionReaction::Dt => Ok((0.25 * n_cm3 * n_cm3 * reactivity, 0.0)),
                    FusionReaction::Dd => Ok((0.0, 0.5 * n_cm3 * n_cm3 * reactivity)),
                }
            }
        }
    }

    /// Relative neutron source density \[neutrons/s/cm³, arbitrary global
    /// scale\] at minor radius `r` \[cm\]: the single-fuel
    /// `f_fuel·n²·⟨σv⟩`, or the sum of the two mixture branches
    /// (module rustdoc rate rule) when `fuel_mixture` is set.
    pub fn strength_density(&self, r: f64) -> Result<f64> {
        let (s_dt, s_dd) = self.branch_strength_density(r)?;
        Ok(s_dt + s_dd)
    }

    /// Relative D(d,p)T proton production density \[protons/s/cm³, same
    /// arbitrary global scale as [`Self::strength_density`] at minor radius
    /// `r` \[cm\]: the pinned v1 accounting convention (module rustdoc) —
    /// the proton branch shares the D-D total rate 50/50 with the neutron
    /// branch, so this IS the D-D neutron-branch density, bit-for-bit
    /// (the branch tuple's second element, not a re-expression). Zero for
    /// D-T-only configurations. Accounting output only: the sampler and the
    /// cards never see it.
    pub fn proton_strength_density(&self, r: f64) -> Result<f64> {
        let (_, s_dd) = self.branch_strength_density(r)?;
        Ok(s_dd)
    }

    /// Total relative proton production strength ∭ P·R·|J| da dθ dφ over the
    /// same volume element as [`Self::total_strength`] (full torus, or the
    /// sector's `rotation_angle` when one is set — the landed `2.0·π`
    /// constant is kept verbatim on the `None` path). Same arbitrary global
    /// scale as the neutron total, so ratios (`proton_per_neutron`) are
    /// meaningful; pure D-D gives exactly the neutron total. D-T-only
    /// configurations make no protons: exact `0.0` (fuel logic, not a
    /// table-build error). A mixture that makes D-D neutrons but has zero
    /// total strength elsewhere (cold profiles) still fails loudly through
    /// the table build.
    pub fn total_proton_strength(&self) -> Result<f64> {
        if let Some(sector) = &self.sector {
            sector.validate()?;
        }
        let makes_protons = match (&self.fuel_mixture, self.fuel) {
            (Some(mixture), _) => mixture.f_deuterium > 0.0,
            (None, FusionReaction::Dd) => true,
            (None, FusionReaction::Dt) => false,
        };
        if !makes_protons {
            return Ok(0.0);
        }
        let table = WeightTable::build_with(self, |r| self.proton_strength_density(r))?;
        let toroidal = match &self.sector {
            None => 2.0 * PI,
            Some(sector) => sector.rotation_angle,
        };
        Ok(toroidal * table.masses.iter().sum::<f64>())
    }

    /// Neutron-branch spectra at species temperatures `(t_d, t_t)` \[keV\]:
    /// `(reaction, weight, mean, sigma)` tuples with weights summing to 1.
    /// Single fuel → one unit-weight branch at the shared temperature; a
    /// mixture → D-T and D-D branches weighted by the rate rule, the D-T
    /// arm at the pair effective temperature and the D-D arm at `t_d`
    /// (the module rustdoc normalization). A [`DeuteriumTail`] expands the
    /// mixture to bulk/tail sub-branches — D-T at `T_DT`/`T_DTt`, D-D at
    /// `T_D`/`T_mix`/`T_tail` — weighted by the sub-rate rule (module
    /// rustdoc); `η == 0` returns the landed two-branch vector verbatim
    /// (bit-for-bit recovery). Where the total mixture rate
    /// underflows — both reactivities zero, a cold annulus — the D-T branch
    /// is returned by convention: the mixture spectrum degenerates to the
    /// 14.021 MeV line exactly as the single-fuel kernels degenerate at
    /// `T_i = 0`.
    fn spectrum_branches(
        &self,
        t_d: f64,
        t_t: f64,
    ) -> Result<Vec<(FusionReaction, f64, f64, f64)>> {
        match &self.fuel_mixture {
            None => {
                let (mu, sigma) = self.fuel.moments_mev(t_d)?;
                Ok(vec![(self.fuel, 1.0, mu, sigma)])
            }
            Some(mixture) => {
                let t_eff = SpeciesIonTemperatures::dt_effective_kev(t_d, t_t);
                let (k_dt, k_dd) = mixture.branch_weights();
                let tail_eta = self.tail.map(|tail| tail.fraction).unwrap_or(0.0);
                if tail_eta == 0.0 {
                    let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(t_eff)?;
                    let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(t_d)?;
                    let w_dt = k_dt * sv_dt;
                    let w_dd = k_dd * sv_dd;
                    let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(t_eff)?;
                    if w_dt + w_dd <= 0.0 {
                        return Ok(vec![(FusionReaction::Dt, 1.0, mu_dt, sigma_dt)]);
                    }
                    let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(t_d)?;
                    let p_dt = w_dt / (w_dt + w_dd);
                    return Ok(vec![
                        (FusionReaction::Dt, p_dt, mu_dt, sigma_dt),
                        (FusionReaction::Dd, 1.0 - p_dt, mu_dd, sigma_dd),
                    ]);
                }
                let tail = self.tail.expect("tail fraction nonzero implies a tail");
                let eta = tail.fraction;
                let t_tail = tail.temperature_kev;
                let t_eff_tail = SpeciesIonTemperatures::dt_effective_kev(t_tail, t_t);
                let t_mix = 0.5 * (t_d + t_tail);
                let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(t_eff)?;
                let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(t_eff_tail)?;
                let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(t_d)?;
                let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(t_mix)?;
                let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(t_tail)?;
                let subs = [
                    (FusionReaction::Dt, k_dt * (1.0 - eta) * sv_dt, t_eff),
                    (FusionReaction::Dt, k_dt * eta * sv_dt_tail, t_eff_tail),
                    (
                        FusionReaction::Dd,
                        k_dd * (1.0 - eta) * (1.0 - eta) * sv_dd,
                        t_d,
                    ),
                    (
                        FusionReaction::Dd,
                        k_dd * 2.0 * eta * (1.0 - eta) * sv_dd_mix,
                        t_mix,
                    ),
                    (FusionReaction::Dd, k_dd * eta * eta * sv_dd_tail, t_tail),
                ];
                let total: f64 = subs.iter().map(|s| s.1).sum();
                let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(t_eff)?;
                if total <= 0.0 {
                    return Ok(vec![(FusionReaction::Dt, 1.0, mu_dt, sigma_dt)]);
                }
                let mut branches = Vec::with_capacity(subs.len());
                for (reaction, weight, ti) in subs {
                    let (mu, sigma) = reaction.moments_mev(ti)?;
                    branches.push((reaction, weight / total, mu, sigma));
                }
                Ok(branches)
            }
        }
    }

    /// Spectrum summary at the magnetic axis (`r = 0`): `(nominal, mean,
    /// sigma)` \[MeV\]. Single fuel → the reaction's nominal line and
    /// Ballabio moments; a mixture → the two-branch Gaussian-mixture
    /// moments (between-branch variance included) and the dominant
    /// branch's nominal line (D-T preferred on exact ties).
    pub fn axis_spectrum_summary(&self) -> Result<(f64, f64, f64)> {
        let (t_d, t_t) = self.species_temperatures_at(0.0);
        match &self.fuel_mixture {
            None => {
                let (mean, sigma) = self.fuel.moments_mev(t_d)?;
                Ok((self.fuel.nominal_energy_mev(), mean, sigma))
            }
            Some(_) => {
                let branches = self.spectrum_branches(t_d, t_t)?;
                let mean: f64 = branches.iter().map(|b| b.1 * b.2).sum();
                let second: f64 = branches.iter().map(|b| b.1 * (b.3 * b.3 + b.2 * b.2)).sum();
                let variance = (second - mean * mean).max(0.0);
                // Strictly-greater scan: D-T (listed first) wins exact ties.
                let mut nominal = FusionReaction::Dt.nominal_energy_mev();
                let mut best = 0.0f64;
                for &(reaction, weight, _, _) in &branches {
                    if weight > best {
                        best = weight;
                        nominal = reaction.nominal_energy_mev();
                    }
                }
                Ok((nominal, mean, variance.sqrt()))
            }
        }
    }

    /// Tail-only neutron source density \[neutrons/s/cm³, same arbitrary
    /// global scale as [`Self::strength_density\]] at minor radius `r` \[cm\]:
    /// the `η`-scaled sub-rate of the module-rustdoc rule — the tail D-T
    /// arm plus the cross and tail-tail D-D arms. `0.0` without a tail or at
    /// `η == 0`. Card-drift input only: the sampler and the strength kernel
    /// never call it.
    pub fn tail_strength_density(&self, r: f64) -> Result<f64> {
        let Some(tail) = &self.tail else {
            return Ok(0.0);
        };
        let eta = tail.fraction;
        if eta == 0.0 {
            return Ok(0.0);
        }
        let Some(mixture) = &self.fuel_mixture else {
            return Ok(0.0);
        };
        let (k_dt, k_dd) = mixture.branch_weights();
        let n_cm3 = self.density_m3(r) * 1e-6;
        let (t_d, t_t) = self.species_temperatures_at(r);
        let t_tail = tail.temperature_kev;
        let t_eff_tail = SpeciesIonTemperatures::dt_effective_kev(t_tail, t_t);
        let t_mix = 0.5 * (t_d + t_tail);
        let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(t_eff_tail)? * 1e6;
        let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(t_mix)? * 1e6;
        let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(t_tail)? * 1e6;
        let s_dt = k_dt * n_cm3 * n_cm3 * eta * sv_dt_tail;
        let s_dd =
            k_dd * n_cm3 * n_cm3 * (2.0 * eta * (1.0 - eta) * sv_dd_mix + eta * eta * sv_dd_tail);
        Ok(s_dt + s_dd)
    }

    /// Total relative tail-neutron strength ∭ T·R·|J| da dθ dφ over the same
    /// volume element as [`Self::total_strength`] (sector-aware). Same
    /// arbitrary global scale, so the tail neutron share
    /// (`total_tail_strength / total_strength`) is meaningful; `0.0`
    /// without a tail. Card-drift input only.
    pub fn total_tail_strength(&self) -> Result<f64> {
        if let Some(sector) = &self.sector {
            sector.validate()?;
        }
        if self.tail.is_none_or(|tail| tail.fraction == 0.0) {
            return Ok(0.0);
        }
        let table = WeightTable::build_with(self, |r| self.tail_strength_density(r))?;
        let toroidal = match &self.sector {
            None => 2.0 * PI,
            Some(sector) => sector.rotation_angle,
        };
        Ok(toroidal * table.masses.iter().sum::<f64>())
    }

    /// Total relative source strength ∭ S·R·|J| da dθ dφ. The toroidal
    /// integral contributes `2π` on the full torus or `rotation_angle` on a
    /// sector (module rustdoc) — the landed `2.0·π` constant is kept
    /// verbatim on the `None` path so the full-torus value is untouched,
    /// and an exact-`2π` sector scales by the identical constant
    /// (bit-for-bit recovery). Arbitrary global scale, but ratios between
    /// configurations are meaningful.
    pub fn total_strength(&self) -> Result<f64> {
        if let Some(sector) = &self.sector {
            sector.validate()?;
        }
        let table = WeightTable::build(self)?;
        let toroidal = match &self.sector {
            None => 2.0 * PI,
            Some(sector) => sector.rotation_angle,
        };
        Ok(toroidal * table.masses.iter().sum::<f64>())
    }
}

/// Piecewise-constant inverse-CDF table over `[0, a_minor]`: cell masses and
/// per-cell poloidal CDFs (so θ sampling needs no per-particle quadrature).
#[derive(Debug)]
struct WeightTable {
    /// Cell edges (R_GRID + 1).
    edges: Vec<f64>,
    /// Cell masses `m_i ≈ S(r_i)·∫R|J|dθ·Δr` (unnormalized).
    masses: Vec<f64>,
    /// Cumulative mass (R_GRID + 1, first 0).
    cumulative: Vec<f64>,
    /// Per-cell poloidal CDFs: THETA_GRID + 1 entries each (first 0, last 1).
    theta_cdfs: Vec<Vec<f64>>,
}

impl WeightTable {
    fn build(config: &ParametricPlasmaConfig) -> Result<Self> {
        Self::build_with(config, |r| config.strength_density(r))
    }

    /// Build over an explicit radial density: the landed [`Self::build`]
    /// passes the neutron [`ParametricPlasmaConfig::strength_density`];
    /// proton accounting passes
    /// [`ParametricPlasmaConfig::proton_strength_density`]. Same grid, same
    /// expression tree — the neutron path is bit-identical to the landed
    /// build.
    fn build_with(
        config: &ParametricPlasmaConfig,
        density: impl Fn(f64) -> Result<f64>,
    ) -> Result<Self> {
        let a = config.geometry.minor_radius_cm;
        let dr = a / R_GRID as f64;
        let mut edges = Vec::with_capacity(R_GRID + 1);
        let mut masses = Vec::with_capacity(R_GRID);
        let mut cumulative = Vec::with_capacity(R_GRID + 1);
        let mut theta_cdfs = Vec::with_capacity(R_GRID);
        let mut acc = 0.0;
        cumulative.push(0.0);
        for i in 0..R_GRID {
            let r = (i as f64 + 0.5) * dr;
            edges.push(i as f64 * dr);
            let strength = density(r)?;
            // Poloidal marginal ∫ R|J| dθ at cell midpoint (uniform θ grid,
            // trapezoid on a periodic integrand).
            let dtheta = 2.0 * PI / THETA_GRID as f64;
            let mut weights = Vec::with_capacity(THETA_GRID);
            let mut w_acc = 0.0;
            for j in 0..THETA_GRID {
                let theta = (j as f64 + 0.5) * dtheta;
                let w = config.geometry.volume_element(r, theta);
                weights.push(w);
                w_acc += w;
            }
            // A zero (or non-finite) accumulated weight would make the
            // per-cell theta CDF a 0/0 division (and the sampler's
            // partial_cmp().unwrap() a panic path); reject loudly at build
            // time instead of relying on the global-strength check below.
            if w_acc.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
                return Err(Error::InvalidProfile(
                    "parametric source cell has zero poloidal weight",
                ));
            }
            let mass = strength * w_acc * dtheta * dr;
            masses.push(mass);
            acc += mass;
            cumulative.push(acc);
            let mut cdf = Vec::with_capacity(THETA_GRID + 1);
            let mut c = 0.0;
            cdf.push(0.0);
            for w in weights {
                c += w / w_acc;
                cdf.push(c);
            }
            theta_cdfs.push(cdf);
        }
        edges.push(a);
        if acc.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || !acc.is_finite() {
            return Err(Error::InvalidProfile(
                "total parametric source strength is zero (check profiles and fuel)",
            ));
        }
        Ok(Self {
            edges,
            masses,
            cumulative,
            theta_cdfs,
        })
    }

    /// Sample a minor radius from the cell masses.
    fn sample_r(&self, u: f64) -> f64 {
        let total = self.cumulative[R_GRID];
        let target = u * total;
        // `total_cmp` (not `partial_cmp().unwrap()`): total order, no panic
        // path even if a corrupt table ever held NaN (construction rejects
        // non-finite masses, so this is defense-in-depth).
        let idx = match self.cumulative.binary_search_by(|v| v.total_cmp(&target)) {
            Ok(i) => i.min(R_GRID - 1),
            Err(i) => i.saturating_sub(1).min(R_GRID - 1),
        };
        let lo = self.edges[idx];
        let hi = self.edges[idx + 1];
        let frac = (target - self.cumulative[idx])
            / (self.cumulative[idx + 1] - self.cumulative[idx]).max(1e-300);
        lo + frac.min(1.0) * (hi - lo)
    }

    /// Sample a poloidal angle from the cell's poloidal CDF.
    fn sample_theta(&self, cell: usize, u: f64) -> f64 {
        let cdf = &self.theta_cdfs[cell];
        let idx = match cdf.binary_search_by(|v| v.total_cmp(&u)) {
            Ok(i) => i.min(THETA_GRID - 1),
            Err(i) => i.saturating_sub(1).min(THETA_GRID - 1),
        };
        let dtheta = 2.0 * PI / THETA_GRID as f64;
        let lo = idx as f64 * dtheta;
        let frac = (u - cdf[idx]) / (cdf[idx + 1] - cdf[idx]).max(1e-300);
        lo + frac.min(1.0) * dtheta
    }
}

/// Deterministic sampler for a parametric plasma source.
pub struct ParametricSampler {
    config: ParametricPlasmaConfig,
    table: WeightTable,
    rng: Rng,
}

impl ParametricSampler {
    /// New sampler over `config` pinned to `seed`.
    pub fn new(config: ParametricPlasmaConfig, seed: u64) -> Result<Self> {
        config.validate()?;
        let table = WeightTable::build(&config)?;
        Ok(Self {
            config,
            table,
            rng: Rng::new(seed),
        })
    }

    /// Borrow the underlying configuration.
    pub fn config(&self) -> &ParametricPlasmaConfig {
        &self.config
    }

    /// Sample one particle: `r` from the strength-weighted radial CDF, `θ`
    /// from the cell's poloidal CDF, `φ` uniform over the full torus or the
    /// sector (module rustdoc); energy from the local
    /// ion-temperature Ballabio Gaussian; isotropic direction.
    ///
    /// The sector remaps the same single toroidal draw (`start + rotation·u`;
    /// `0 + 2π·u == 2π·u`, so the exact full-rotation spelling reproduces
    /// the landed stream bit-for-bit) and consumes no new RNG draws on any
    /// path.
    ///
    /// The `Err(_) => 0.0` arm is unreachable-in-practice: the config is
    /// validated at construction and `temperature_kev(r)` evaluates a
    /// caller profile already accepted by `WeightTable::build`, so
    /// `moments_mev` cannot fail here. The zero fallback (not a panic) keeps
    /// the infallible sampling stream total; a fallible caller should use
    /// [`ParametricSampler::try_sample`].
    pub fn sample(&mut self) -> Particle {
        match self.try_sample() {
            Ok(particle) => particle,
            Err(_) => Particle {
                position_cm: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
                energy_mev: 0.0,
                weight: self.config.weight,
            },
        }
    }

    /// Draw the birth spectrum moments `(mean, sigma)` \[MeV\] at species
    /// temperatures `(t_d, t_t)`: the config fuel's Ballabio moments, or —
    /// for a fuel mixture — the branch roulette weighted by the rate rule
    /// (module rustdoc) followed by the chosen branch's moments.
    /// The single-fuel path consumes no RNG, so the landed sampling stream
    /// is preserved bit-for-bit; the mixture path draws one extra uniform
    /// per particle. Both reactivities zero (a cold annulus) degenerates to
    /// the D-T line by the documented convention of
    /// [`ParametricPlasmaConfig::spectrum_branches`].
    fn sample_spectrum_moments(&mut self, t_d: f64, t_t: f64) -> Result<(f64, f64)> {
        let branches = self.config.spectrum_branches(t_d, t_t)?;
        if branches.len() == 1 {
            return Ok((branches[0].2, branches[0].3));
        }
        let u = self.rng.uniform();
        let mut acc = 0.0;
        let mut picked = branches[0];
        for &branch in &branches {
            acc += branch.1;
            if u < acc {
                picked = branch;
                break;
            }
        }
        Ok((picked.2, picked.3))
    }

    /// Fallible single-particle sample: like [`ParametricSampler::sample`]
    /// but surfaces a spectrum-moment failure as a loud [`Error`] instead of
    /// the documented zero-energy fallback.
    pub fn try_sample(&mut self) -> Result<Particle> {
        let r = self.table.sample_r(self.rng.uniform());
        let dr = self.config.geometry.minor_radius_cm / R_GRID as f64;
        let cell = ((r / dr).floor() as usize).min(R_GRID - 1);
        let theta = self.table.sample_theta(cell, self.rng.uniform());
        let phi = match &self.config.sector {
            None => 2.0 * PI * self.rng.uniform(),
            Some(sector) => sector.map_angle(self.rng.uniform()),
        };
        let (big_r, z) = self.config.geometry.map(r, theta);
        let (t_d, t_t) = self.config.species_temperatures_at(r);
        let (mean, sigma) = self.sample_spectrum_moments(t_d, t_t)?;
        let energy_mev = if sigma > 0.0 {
            mean + sigma * self.rng.standard_normal()
        } else {
            mean
        };
        Ok(Particle {
            position_cm: [big_r * phi.cos(), big_r * phi.sin(), z],
            direction: self.rng.isotropic_direction(),
            energy_mev,
            weight: self.config.weight,
        })
    }

    /// Sample `n` particles.
    pub fn sample_n(&mut self, n: usize) -> Vec<Particle> {
        (0..n).map(|_| self.sample()).collect()
    }
}

/// A binned marginal distribution for source-card emission.
#[derive(Debug, Clone, PartialEq)]
pub struct BinnedDistribution {
    /// Bin centers, ascending.
    pub centers: Vec<f64>,
    /// Normalized bin masses (sum 1.0).
    pub masses: Vec<f64>,
}

/// Marginal histograms of the parametric source for card emission:
/// radial birth profile, vertical birth profile, and the global marginal
/// energy spectrum (built on the fine grid; `bins` controls the output
/// binning of each marginal).
#[derive(Debug, Clone, PartialEq)]
pub struct EmissionHistograms {
    /// Radial marginal (birth minor-radius profile).
    pub radial: BinnedDistribution,
    /// Vertical marginal (birth Z profile over `[-z_max, z_max]`).
    pub vertical: BinnedDistribution,
    /// Global marginal energy spectrum \[MeV\].
    pub energy: BinnedDistribution,
    /// Uniform toroidal-angle marginal over
    /// `[start, start + rotation)` \[rad\] (`bins` bins): `Some` for a
    /// partial sector, `None` on the full torus (no sector or an exact
    /// full rotation — card emission then keeps the landed three-marginal
    /// shape bit-for-bit).
    pub toroidal: Option<BinnedDistribution>,
    /// Captured energy probability mass (1 − tail truncation).
    pub energy_coverage: f64,
    /// Half the L1 distance between the true `(r, z)` birth joint and the
    /// product of the two marginals — the correlation information a
    /// product-form source card cannot carry (0 = independent).
    pub joint_correlation: f64,
}

/// Energy tabulation half-width in sigma (same convention as the ring/point
/// spectrum tabulation).
const ENERGY_WIDTH_SIGMA: f64 = 4.0;

/// Build the marginal histograms for one configuration.
pub fn emission_histograms(
    config: &ParametricPlasmaConfig,
    bins: usize,
) -> Result<EmissionHistograms> {
    config.validate()?;
    let bins = bins.max(4);
    let g = &config.geometry;
    let a = g.minor_radius_cm;
    let dr = a / R_GRID as f64;
    let dtheta = 2.0 * PI / THETA_GRID as f64;

    let mut radial_m = vec![0.0f64; bins];
    let mut vertical_m = vec![0.0f64; bins];
    let z_max = g.z_max();
    let mut energy_lo = f64::INFINITY;
    let mut energy_hi = f64::NEG_INFINITY;
    let mut cell_weights: Vec<f64> = Vec::with_capacity(R_GRID * THETA_GRID);
    let mut cell_r: Vec<f64> = Vec::with_capacity(R_GRID * THETA_GRID);
    let mut cell_z: Vec<f64> = Vec::with_capacity(R_GRID * THETA_GRID);
    let mut total = 0.0f64;
    for i in 0..R_GRID {
        let r = (i as f64 + 0.5) * dr;
        let strength = config.strength_density(r)?;
        let (t_d, t_t) = config.species_temperatures_at(r);
        for &(_, _, mu, sigma) in &config.spectrum_branches(t_d, t_t)? {
            energy_lo = energy_lo.min(mu - ENERGY_WIDTH_SIGMA * sigma);
            energy_hi = energy_hi.max(mu + ENERGY_WIDTH_SIGMA * sigma);
        }
        for j in 0..THETA_GRID {
            let theta = (j as f64 + 0.5) * dtheta;
            let (_, z) = g.map(r, theta);
            let w = strength * g.volume_element(r, theta) * dtheta * dr;
            total += w;
            cell_weights.push(w);
            cell_r.push(r);
            cell_z.push(z);
        }
    }
    if total.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return Err(Error::InvalidProfile(
            "total parametric source strength is zero",
        ));
    }
    let energy_lo = if energy_lo.is_finite() {
        energy_lo
    } else {
        0.0
    };
    let energy_hi = if energy_hi.is_finite() {
        energy_hi
    } else {
        1.0
    };

    // Bin the spatial marginals and accumulate the joint for the
    // correlation distance.
    let mut joint = vec![vec![0.0f64; bins]; bins];
    for k in 0..cell_weights.len() {
        let (w, r, z) = (cell_weights[k], cell_r[k], cell_z[k]);
        let ri = ((r / a) * bins as f64).floor() as usize;
        let ri = ri.min(bins - 1);
        let zi = (((z + z_max) / (2.0 * z_max)) * bins as f64).floor() as usize;
        let zi = zi.min(bins - 1);
        radial_m[ri] += w;
        vertical_m[zi] += w;
        joint[ri][zi] += w / total;
    }
    let radial: Vec<f64> = radial_m.iter().map(|m| m / total).collect();
    let vertical: Vec<f64> = vertical_m.iter().map(|m| m / total).collect();
    // joint_correlation = 1/2 Σ |p_true(i,j) − p_r(i) p_z(j)|
    let mut corr = 0.0;
    for i in 0..bins {
        for j in 0..bins {
            corr += (joint[i][j] - radial[i] * vertical[j]).abs();
        }
    }
    corr *= 0.5;

    // Marginal energy spectrum: Σ_cells w_cell · [Φ(hi; μ,σ) − Φ(lo; μ,σ)].
    let de = (energy_hi - energy_lo) / bins as f64;
    let mut energy_m = vec![0.0f64; bins];
    let mut energy_coverage = 0.0f64;
    for i in 0..R_GRID {
        let r = (i as f64 + 0.5) * dr;
        let (t_d, t_t) = config.species_temperatures_at(r);
        let branches = config.spectrum_branches(t_d, t_t)?;
        let cell_mass: f64 = cell_weights[i * THETA_GRID..(i + 1) * THETA_GRID]
            .iter()
            .sum();
        for (b, mass) in energy_m.iter_mut().enumerate() {
            let lo = energy_lo + b as f64 * de;
            let hi = lo + de;
            // Branch-mixture bin probability: Σ w_b · [Φ(hi; μ_b, σ_b) −
            // Φ(lo; μ_b, σ_b)] (the monoenergetic arm per branch).
            let mut p = 0.0;
            for &(_, weight, mu, sigma) in &branches {
                let pb = if sigma > 0.0 {
                    crate::spectrum::normal_cdf((hi - mu) / sigma)
                        - crate::spectrum::normal_cdf((lo - mu) / sigma)
                } else if lo < mu && mu <= hi {
                    1.0
                } else {
                    0.0
                };
                p += weight * pb;
            }
            *mass += cell_mass * p;
        }
        let mut covered = 0.0;
        for &(_, weight, mu, sigma) in &branches {
            let cb = if sigma > 0.0 {
                crate::spectrum::normal_cdf((energy_hi - mu) / sigma)
                    - crate::spectrum::normal_cdf((energy_lo - mu) / sigma)
            } else {
                1.0
            };
            covered += weight * cb;
        }
        energy_coverage += cell_mass * covered;
    }
    energy_coverage /= total;
    let energy_centers: Vec<f64> = (0..bins)
        .map(|b| energy_lo + (b as f64 + 0.5) * de)
        .collect();
    let energy_masses: Vec<f64> = energy_m.iter().map(|m| m / total).collect();

    // Toroidal-angle marginal: the birth angle is uniform over the sector
    // (module rustdoc), so the marginal is exact, not sampled — bin centers
    // over [start, start + rotation) with equal masses. The full torus
    // carries no marginal (the landed card shape is untouched).
    let toroidal = match config.sector {
        Some(sector) if !sector.is_full_rotation() => {
            let width = sector.rotation_angle / bins as f64;
            Some(BinnedDistribution {
                centers: (0..bins)
                    .map(|b| sector.start_angle + (b as f64 + 0.5) * width)
                    .collect(),
                masses: vec![1.0 / bins as f64; bins],
            })
        }
        _ => None,
    };

    Ok(EmissionHistograms {
        radial: BinnedDistribution {
            centers: (0..bins)
                .map(|b| (b as f64 + 0.5) * (a / bins as f64))
                .collect(),
            masses: radial,
        },
        vertical: BinnedDistribution {
            centers: (0..bins)
                .map(|b| -z_max + (b as f64 + 0.5) * (2.0 * z_max / bins as f64))
                .collect(),
            masses: vertical,
        },
        energy: BinnedDistribution {
            centers: energy_centers,
            masses: energy_masses,
        },
        toroidal,
        energy_coverage,
        joint_correlation: corr,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::miller::MillerGeometry;

    /// ITER-ish synthetic H-mode case (hand round numbers).
    pub(crate) fn iter_h_mode() -> ParametricPlasmaConfig {
        ParametricPlasmaConfig {
            geometry: MillerGeometry {
                major_radius_cm: 620.0,
                minor_radius_cm: 200.0,
                elongation: 1.85,
                triangularity: 0.35,
                shafranov_factor_cm: 15.0,
            },
            mode: ProfileMode::H,
            ion_density: DensityProfile {
                centre_m3: 1.2e20,
                peaking_factor: 1.1,
                pedestal_m3: 4.0e19,
                separatrix_m3: 3.0e19,
            },
            ion_temperature: TemperatureProfile {
                centre_kev: 28.0,
                peaking_factor: 2.5,
                beta: 2.0,
                pedestal_kev: 4.0,
                separatrix_kev: 0.1,
            },
            pedestal_radius_cm: 150.0,
            fuel: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            tail: None,
            sector: None,
            weight: 1.0,
        }
    }

    /// Flat-profile configuration (constant n and T): strength ∝ volume
    /// element, so moments have closed forms (see module docs and gates).
    fn flat_config(kappa: f64, delta: f64, esh: f64) -> ParametricPlasmaConfig {
        let mut c = iter_h_mode();
        c.mode = ProfileMode::L;
        c.geometry.elongation = kappa;
        c.geometry.triangularity = delta;
        c.geometry.shafranov_factor_cm = esh;
        c.ion_density = DensityProfile {
            centre_m3: 1.0e20,
            peaking_factor: 0.0,
            pedestal_m3: 1.0,
            separatrix_m3: 1.0,
        };
        c.ion_temperature = TemperatureProfile {
            centre_kev: 20.0,
            peaking_factor: 0.0,
            beta: 1.0,
            pedestal_kev: 20.0,
            separatrix_kev: 20.0,
        };
        c
    }

    // --- G-strength: hand strength cases -------------------------------------

    #[test]
    fn strength_density_hand_cases() {
        let c = iter_h_mode();
        // n(0) = 1.2e20 m-3, T(0) = 28 keV.
        let s0 = c.strength_density(0.0).unwrap();
        let n_cm3 = 1.2e20 * 1e-6;
        let sv = FusionReaction::Dt.reactivity_m3_per_s(28.0).unwrap() * 1e6;
        let want = 0.25 * n_cm3 * n_cm3 * sv;
        assert!((s0 - want).abs() < 1e-12 * want, "{s0} vs {want}");
        // Cold separatrix: T(a) = 0.1 keV, tiny but non-zero strength.
        assert!(c.strength_density(200.0).unwrap() > 0.0);
        // L-mode flat density 1e20, T 20: exact factor 1/4·(1e14)^2·<σv>(20).
        let f = flat_config(1.85, 0.35, 15.0);
        let sf = f.strength_density(123.0).unwrap();
        let sv20 = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap() * 1e6;
        let wantf = 0.25 * 1e28 * sv20;
        assert!((sf - wantf).abs() < 1e-12 * wantf, "{sf} vs {wantf}");
    }

    #[test]
    fn total_strength_is_positive_and_scales_with_reactivity() {
        let hot = flat_config(1.0, 0.0, 0.0);
        let mut cold = hot;
        cold.ion_temperature.centre_kev = 10.0;
        cold.ion_temperature.peaking_factor = 0.0;
        cold.ion_temperature.pedestal_kev = 10.0;
        cold.ion_temperature.separatrix_kev = 10.0;
        let s_hot = hot.total_strength().unwrap();
        let s_cold = cold.total_strength().unwrap();
        let ratio = s_hot / s_cold;
        let want = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap()
            / FusionReaction::Dt.reactivity_m3_per_s(10.0).unwrap();
        assert!((ratio - want).abs() < 1e-3 * want, "{ratio} vs {want}");
    }

    // --- G-mix: mixture normalization, recovery anchors, loud errors -------

    #[test]
    fn mixture_branch_weights_are_exact() {
        // Hand vectors at exact equality: the recovery anchors pinned in
        // the module rustdoc.
        let half = FuelMixture::new(0.5, 0.5).unwrap();
        assert_eq!(half.branch_weights(), (0.25, 0.125));
        let pure_d = FuelMixture::new(1.0, 0.0).unwrap();
        assert_eq!(pure_d.branch_weights(), (0.0, 0.5));
        let pure_t = FuelMixture::new(0.0, 1.0).unwrap();
        assert_eq!(pure_t.branch_weights(), (0.0, 0.0));
    }

    #[test]
    fn mixture_strength_hand_vectors() {
        // Flat L-mode profile: n = 1e20 m-3, T = 20 keV everywhere — the
        // module rustdoc hand vectors. Exact expression equality against
        // the branch decomposition, plus the pinned decimal goldens at
        // 1e-12 (the reactivity-transcription gate precision).
        let flat = flat_config(1.85, 0.0, 0.0);
        let n_cm3 = flat.density_m3(123.0) * 1e-6;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap() * 1e6;
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap() * 1e6;

        let mut half = flat;
        half.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        let s = half.strength_density(123.0).unwrap();
        assert_eq!(
            s,
            0.25 * n_cm3 * n_cm3 * sv_dt + 0.125 * n_cm3 * n_cm3 * sv_dd
        );
        let want = 1085808328120.6956_f64; // module rustdoc hand vector
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");

        let mut pure = flat;
        pure.fuel_mixture = Some(FuelMixture::new(1.0, 0.0).unwrap());
        let s = pure.strength_density(123.0).unwrap();
        assert_eq!(s, 0.5 * n_cm3 * n_cm3 * sv_dd);
        let want = 13012914793.607618_f64;
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");

        let mut blend = flat;
        blend.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        let s = blend.strength_density(123.0).unwrap();
        let (k_dt, k_dd) = FuelMixture::new(0.7, 0.3).unwrap().branch_weights();
        assert_eq!(
            s,
            k_dt * n_cm3 * n_cm3 * sv_dt + k_dd * n_cm3 * n_cm3 * sv_dd
        );
        let want = 915722611763.5945_f64;
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");
    }

    #[test]
    fn equimolar_mixture_dt_branch_recovers_landed_kernel_bit_for_bit() {
        let landed = iter_h_mode(); // fuel: Dt, fuel_mixture: None
        let mut mixture = landed;
        mixture.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            let (s_dt, s_dd) = mixture.branch_strength_density(r).unwrap();
            assert_eq!(
                s_dt,
                landed.strength_density(r).unwrap(),
                "D-T branch must equal the landed equimolar kernel at r={r}"
            );
            assert!(
                s_dd > 0.0,
                "the D-D branch contributes at equimolar (module rustdoc): r={r}"
            );
        }
    }

    #[test]
    fn pure_deuterium_mixture_recovers_landed_kernel_bit_for_bit() {
        let mut landed = iter_h_mode();
        landed.fuel = FusionReaction::Dd;
        let mut mixture = landed;
        mixture.fuel_mixture = Some(FuelMixture::new(1.0, 0.0).unwrap());
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            assert_eq!(
                mixture.strength_density(r).unwrap(),
                landed.strength_density(r).unwrap(),
                "pure-D total must equal the landed D-D kernel at r={r}"
            );
        }
    }

    // --- G-branches: T-T coefficient, D(d,p)T proton accounting, zero-T anchor

    #[test]
    fn tt_neutron_coefficient_hand_vectors_are_exact() {
        // Pinned three-branch normalization (module rustdoc): f_T² with the
        // multiplicity folded in — exact at every anchor.
        assert_eq!(
            FuelMixture::new(0.5, 0.5).unwrap().tt_neutron_coefficient(),
            0.25
        );
        assert_eq!(
            FuelMixture::new(1.0, 0.0).unwrap().tt_neutron_coefficient(),
            0.0
        );
        assert_eq!(
            FuelMixture::new(0.7, 0.3).unwrap().tt_neutron_coefficient(),
            0.09
        );
        assert_eq!(
            FuelMixture::new(0.0, 1.0).unwrap().tt_neutron_coefficient(),
            1.0
        );
    }

    #[test]
    fn proton_strength_is_the_dd_branch_bit_for_bit() {
        // The pinned 50/50 accounting convention (module rustdoc): the
        // proton density IS the D-D neutron-branch density — asserted with
        // assert_eq, not a tolerance, on both the flat hand-vector config
        // and the H-mode config, with and without a fuel mixture.
        let flat = flat_config(1.85, 0.0, 0.0);
        let n_cm3 = flat.density_m3(123.0) * 1e-6;
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap() * 1e6;
        let mut blend = flat;
        blend.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        let p = blend.proton_strength_density(123.0).unwrap();
        let (_, s_dd) = blend.branch_strength_density(123.0).unwrap();
        assert_eq!(p, s_dd);
        assert_eq!(p, 0.7 * 0.7 / 2.0 * n_cm3 * n_cm3 * sv_dd);
        let want = 6376328248.867733_f64; // module rustdoc hand vector
        assert!((p - want).abs() < 1e-12 * want, "{p} vs {want}");

        // Single-fuel D-D: proton density equals the landed kernel.
        let mut dd = flat;
        dd.fuel = FusionReaction::Dd;
        for &r in &[0.0, 50.0, 123.0, 199.9] {
            assert_eq!(
                dd.proton_strength_density(r).unwrap(),
                dd.strength_density(r).unwrap(),
            );
        }
        // D-T-only configurations make no protons: exact zero, no residue.
        let dt = flat; // fuel: Dt, fuel_mixture: None
        for &r in &[0.0, 123.0, 199.9] {
            assert_eq!(dt.proton_strength_density(r).unwrap(), 0.0);
        }
        assert_eq!(dt.total_proton_strength().unwrap(), 0.0);
        let mut equimolar = flat;
        equimolar.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        assert_eq!(
            equimolar.proton_strength_density(123.0).unwrap(),
            3253228698.4019046_f64
        );
        // H-mode mixture: bit-identity holds off the flat profile too.
        let mut h = iter_h_mode();
        h.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            let (_, s_dd) = h.branch_strength_density(r).unwrap();
            assert_eq!(h.proton_strength_density(r).unwrap(), s_dd);
        }
    }

    #[test]
    fn no_tritium_mixture_recovers_two_branch_kernel_bit_for_bit() {
        // Zero-tritium recovery anchor (module rustdoc): with f_T = 0 the
        // absent D-T and T-T branches contribute exactly 0.0, so the no-T
        // mixture reproduces the landed two-branch kernel bit-for-bit —
        // neutron total and proton density alike.
        let mut landed = iter_h_mode();
        landed.fuel = FusionReaction::Dd;
        let mut mixture = landed;
        mixture.fuel_mixture = Some(FuelMixture::new(1.0, 0.0).unwrap());
        assert_eq!(
            mixture.fuel_mixture.unwrap().tt_neutron_coefficient(),
            0.0,
            "the deferred T-T term vanishes exactly at f_T = 0"
        );
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            assert_eq!(
                mixture.strength_density(r).unwrap(),
                landed.strength_density(r).unwrap(),
            );
            assert_eq!(
                mixture.proton_strength_density(r).unwrap(),
                landed.strength_density(r).unwrap(),
            );
        }
        assert_eq!(
            mixture.total_proton_strength().unwrap(),
            landed.total_strength().unwrap(),
            "pure-D proton total must equal the landed D-D neutron total"
        );
    }

    #[test]
    fn total_proton_strength_matches_branch_share() {
        // Flat 70/30 blend: uniform temperature, so the proton share of the
        // neutron total equals the D-D branch weight ratio (grid-quadrature
        // precision, 1e-12); pure D-D gives proton/neutron == 1 within
        // summation order.
        let flat = flat_config(1.85, 0.0, 0.0);
        let mut blend = flat;
        blend.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        let (k_dt, k_dd) = FuelMixture::new(0.7, 0.3).unwrap().branch_weights();
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let p_dd = k_dd * sv_dd / (k_dt * sv_dt + k_dd * sv_dd);
        let n = blend.total_strength().unwrap();
        let p = blend.total_proton_strength().unwrap();
        assert!(
            (p / n - p_dd).abs() < 1e-12,
            "proton share {} vs {p_dd}",
            p / n
        );
        let mut dd = flat;
        dd.fuel = FusionReaction::Dd;
        let rel = (dd.total_proton_strength().unwrap() / dd.total_strength().unwrap() - 1.0).abs();
        assert!(rel < 1e-12, "pure-D proton/neutron {rel}");
    }

    #[test]
    fn mixture_fractions_are_loud() {
        assert_eq!(
            FuelMixture::new(f64::NAN, 0.5),
            Err(Error::NonFinite("deuterium fraction"))
        );
        assert_eq!(
            FuelMixture::new(0.5, f64::INFINITY),
            Err(Error::NonFinite("tritium fraction"))
        );
        assert_eq!(
            FuelMixture::new(-0.1, 1.1),
            Err(Error::InvalidFuelMixture("fuel fractions must be >= 0"))
        );
        assert_eq!(
            FuelMixture::new(0.5, 0.6),
            Err(Error::InvalidFuelMixture("fuel fractions must sum to 1"))
        );
        // Last-ulp missums are admitted (fractions used as given).
        assert!(FuelMixture::new(1.0 / 3.0, 2.0 / 3.0).is_ok());
        // Config-level validation surfaces the same named errors.
        let mut c = iter_h_mode();
        c.fuel_mixture = Some(FuelMixture {
            f_deuterium: -0.5,
            f_tritium: 1.5,
        });
        assert!(matches!(
            c.validate(),
            Err(Error::InvalidFuelMixture("fuel fractions must be >= 0"))
        ));
    }

    // --- G-species: per-species ion temperatures, recovery anchor, loud errors

    #[test]
    fn species_effective_temperature_hand_vectors() {
        use SpeciesIonTemperatures as ST;
        // Equal temperatures recover exactly (bit-for-bit anchor input).
        assert_eq!(ST::dt_effective_kev(20.0, 20.0), 20.0);
        assert_eq!(ST::dt_effective_kev(0.0, 0.0), 0.0);
        // Hand vector: T_D = 20, T_T = 30 → 20 + (2/5)·10 = 24.0 exactly.
        assert_eq!(ST::dt_effective_kev(20.0, 30.0), 24.0);
        assert_eq!(
            ST::dt_effective_kev(20.0, 30.0),
            20.0 + (2.0 / 5.0) * (30.0 - 20.0)
        );
        // Mass weighting favors the heavier species' partner: heating T
        // moves T_DT less than heating D by the same amount.
        assert!(ST::dt_effective_kev(20.0, 30.0) < ST::dt_effective_kev(30.0, 20.0));
        assert_eq!(ST::new(20.0, 30.0).unwrap().dt_effective(), 24.0);
    }

    #[test]
    fn species_strength_hand_vectors() {
        // Flat L-mode profile: n = 1e20 m⁻³ everywhere; f_D = 0.7,
        // f_T = 0.3, T_D = 20 keV, T_T = 30 keV — the module rustdoc hand
        // vector (T_DT = 24 keV). Exact expression equality against the
        // branch decomposition, plus the pinned decimal golden at 1e-12.
        let flat = flat_config(1.85, 0.0, 0.0);
        let mut species = flat;
        species.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        species.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        let n_cm3 = flat.density_m3(123.0) * 1e-6;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap() * 1e6;
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap() * 1e6;
        let s = species.strength_density(123.0).unwrap();
        assert_eq!(
            s,
            0.7 * 0.3 * n_cm3 * n_cm3 * sv_dt + 0.7 * 0.7 / 2.0 * n_cm3 * n_cm3 * sv_dd
        );
        let want = 1143456467112.5283_f64; // module rustdoc hand vector
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");
    }

    #[test]
    fn species_equal_temperatures_recover_shared_mixture_bit_for_bit() {
        // T_D = T_T = profile T recovers the shared-temperature mixture
        // kernel bit-for-bit: T_DT computes to T_D exactly, so strengths
        // and the axis spectrum summary are identical.
        let mut shared = iter_h_mode();
        shared.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        let mut species = shared;
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            let ti = shared.temperature_kev(r);
            species.species_temperatures = Some(SpeciesIonTemperatures::new(ti, ti).unwrap());
            assert_eq!(
                species.strength_density(r).unwrap(),
                shared.strength_density(r).unwrap(),
                "per-species strength must equal the shared kernel at r={r}"
            );
        }
        let ti0 = shared.temperature_kev(0.0);
        species.species_temperatures = Some(SpeciesIonTemperatures::new(ti0, ti0).unwrap());
        assert_eq!(
            species.axis_spectrum_summary().unwrap(),
            shared.axis_spectrum_summary().unwrap(),
        );
        // The equimolar pair-equal total also matches the landed hand
        // vector decimal golden.
        let mut flat = flat_config(1.85, 0.0, 0.0);
        flat.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        flat.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 20.0).unwrap());
        let s = flat.strength_density(123.0).unwrap();
        let want = 1085808328120.6956_f64;
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");
    }

    #[test]
    fn species_branch_spectra_follow_pair_temperatures() {
        // The D-T line sits at T_DT = 24 keV, the D-D line at T_D = 20 keV
        // (module rustdoc convention): per-branch Ballabio moments pinned
        // to the independent goldens.
        let flat = flat_config(1.85, 0.0, 0.0);
        let mut species = flat;
        species.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        species.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        let branches = species.spectrum_branches(20.0, 30.0).unwrap();
        assert_eq!(branches.len(), 2);
        let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(24.0).unwrap();
        let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        assert_eq!((branches[0].2, branches[0].3), (mu_dt, sigma_dt));
        assert_eq!((branches[1].2, branches[1].3), (mu_dd, sigma_dd));
        assert!((mu_dt - 14.077934372230146).abs() < 1e-12, "{mu_dt}");
        assert!((sigma_dt - 0.37003905085517863).abs() < 1e-12, "{sigma_dt}");
    }

    #[test]
    fn species_temperatures_are_loud() {
        assert_eq!(
            SpeciesIonTemperatures::new(f64::NAN, 20.0),
            Err(Error::NonFinite("deuterium ion temperature"))
        );
        assert_eq!(
            SpeciesIonTemperatures::new(20.0, f64::INFINITY),
            Err(Error::NonFinite("tritium ion temperature"))
        );
        assert_eq!(
            SpeciesIonTemperatures::new(-1.0, 20.0),
            Err(Error::NegativeIonTemperature(-1.0))
        );
        assert_eq!(
            SpeciesIonTemperatures::new(20.0, -0.5),
            Err(Error::NegativeIonTemperature(-0.5))
        );
        // Config-level validation surfaces the same named errors.
        let mut c = iter_h_mode();
        c.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        c.species_temperatures = Some(SpeciesIonTemperatures {
            deuterium_kev: 20.0,
            tritium_kev: -1.0,
        });
        assert!(matches!(
            c.validate(),
            Err(Error::NegativeIonTemperature(_))
        ));
        // The pair without a fuel mixture is a loud scope error, never a
        // silent single-fuel fallback.
        let mut single = iter_h_mode();
        single.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        assert!(matches!(single.validate(), Err(Error::NotYetSupported(_))));
        assert!(matches!(
            ParametricSampler::new(single, 7),
            Err(Error::NotYetSupported(_))
        ));
    }

    // --- G-tail: deuterium hot-tail fraction, recovery anchor, loud errors ---

    /// Flat 70/30 blend at T_D = 20 keV, T_T = 30 keV with the pinned
    /// tail (η = 0.05 at T_tail = 60 keV — the module rustdoc hand vector).
    fn tail_config() -> ParametricPlasmaConfig {
        let mut c = flat_config(1.85, 0.0, 0.0);
        c.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        c.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        c.tail = Some(DeuteriumTail::new(0.05, 60.0).unwrap());
        c
    }

    #[test]
    fn tail_strength_hand_vectors() {
        // Flat L-mode profile: n = 1e20 m⁻³ everywhere; the module rustdoc
        // hand vector (T_DT = 24, T_DTt = 48, T_mix = 40 keV). Exact
        // expression equality against the sub-rate decomposition, plus the
        // pinned decimal golden at 1e-12.
        let c = tail_config();
        let n_cm3 = 1e14_f64;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap() * 1e6;
        let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(48.0).unwrap() * 1e6;
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap() * 1e6;
        let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(40.0).unwrap() * 1e6;
        let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(60.0).unwrap() * 1e6;
        let eta = 0.05_f64;
        let (s_dt, s_dd) = c.branch_strength_density(123.0).unwrap();
        assert_eq!(
            s_dt,
            0.7 * 0.3 * n_cm3 * n_cm3 * ((1.0 - eta) * sv_dt + eta * sv_dt_tail)
        );
        assert_eq!(
            s_dd,
            0.7 * 0.7 / 2.0
                * n_cm3
                * n_cm3
                * ((1.0 - eta) * (1.0 - eta) * sv_dd
                    + 2.0 * eta * (1.0 - eta) * sv_dd_mix
                    + eta * eta * sv_dd_tail)
        );
        let want_dt = 1170076195103.095_f64; // sub-rate hand vector (1e-12)
        assert!(
            (s_dt - want_dt).abs() < 1e-12 * want_dt,
            "{s_dt} vs {want_dt}"
        );
        let want = 1177836136801.1577_f64; // module rustdoc hand vector
        let s = c.strength_density(123.0).unwrap();
        assert_eq!(s, s_dt + s_dd);
        assert!((s - want).abs() < 1e-12 * want, "{s} vs {want}");
        // The tail reactivity goldens pin the transcription at the new
        // sub-pair temperatures.
        let want_tail = 8.557148874535017e-16_f64;
        assert!((sv_dt_tail - want_tail).abs() < 1e-15 * want_tail);
        let want_mix = 8.234915801644433e-18_f64;
        assert!((sv_dd_mix - want_mix).abs() < 1e-15 * want_mix);
        let want_ddt = 1.4470008265609935e-17_f64;
        assert!((sv_dd_tail - want_ddt).abs() < 1e-15 * want_ddt);
        // Proton accounting stays bit-identical with the (tail-inclusive)
        // D-D neutron branch.
        assert_eq!(c.proton_strength_density(123.0).unwrap(), s_dd);
        // Tail-only density is the η-scaled sub-rate (the bulk depletes by
        // the complementary share, so total − base is smaller — assert the
        // direct sub-rate expression instead, exact tree).
        assert_eq!(
            c.tail_strength_density(123.0).unwrap(),
            0.7 * 0.3 * n_cm3 * n_cm3 * eta * sv_dt_tail
                + 0.7 * 0.7 / 2.0
                    * n_cm3
                    * n_cm3
                    * (2.0 * eta * (1.0 - eta) * sv_dd_mix + eta * eta * sv_dd_tail)
        );
    }

    #[test]
    fn tail_sub_branch_spectra_follow_sub_pair_temperatures() {
        // D-T carries two lines (24 / 48 keV), D-D three (20 / 40 / 60 keV),
        // weighted by the sub-rate rule; per-branch Ballabio moments pinned
        // to the independent goldens.
        let c = tail_config();
        let branches = c.spectrum_branches(20.0, 30.0).unwrap();
        assert_eq!(branches.len(), 5);
        let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(24.0).unwrap();
        let (mu_dt_tail, sigma_dt_tail) = FusionReaction::Dt.moments_mev(48.0).unwrap();
        let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        let (mu_dd_mix, sigma_dd_mix) = FusionReaction::Dd.moments_mev(40.0).unwrap();
        let (mu_dd_tail, sigma_dd_tail) = FusionReaction::Dd.moments_mev(60.0).unwrap();
        assert_eq!((branches[0].2, branches[0].3), (mu_dt, sigma_dt));
        assert_eq!((branches[1].2, branches[1].3), (mu_dt_tail, sigma_dt_tail));
        assert_eq!((branches[2].2, branches[2].3), (mu_dd, sigma_dd));
        assert_eq!((branches[3].2, branches[3].3), (mu_dd_mix, sigma_dd_mix));
        assert_eq!((branches[4].2, branches[4].3), (mu_dd_tail, sigma_dd_tail));
        assert_eq!(branches[0].0, FusionReaction::Dt);
        assert_eq!(branches[4].0, FusionReaction::Dd);
        assert!(
            (mu_dt_tail - 14.104551866140476).abs() < 1e-12,
            "{mu_dt_tail}"
        );
        assert!((sigma_dt_tail - 0.5241295764702804).abs() < 1e-12);
        assert!((mu_dd_mix - 2.553611770233373).abs() < 1e-12, "{mu_dd_mix}");
        assert!(
            (mu_dd_tail - 2.5984132221035106).abs() < 1e-12,
            "{mu_dd_tail}"
        );
        // Weights sum to 1 and match the sub-rate rule.
        let sum: f64 = branches.iter().map(|b| b.1).sum();
        assert!((sum - 1.0).abs() < 1e-12, "branch weights sum {sum}");
        let eta = 0.05_f64;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap();
        let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(48.0).unwrap();
        let w_dt = 0.7 * 0.3 * ((1.0 - eta) * sv_dt + eta * sv_dt_tail);
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(40.0).unwrap();
        let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(60.0).unwrap();
        let w_dd = 0.7 * 0.7 / 2.0
            * ((1.0 - eta) * (1.0 - eta) * sv_dd
                + 2.0 * eta * (1.0 - eta) * sv_dd_mix
                + eta * eta * sv_dd_tail);
        assert!((branches[0].1 + branches[1].1 - w_dt / (w_dt + w_dd)).abs() < 1e-12);
        assert!(
            (branches[2].1 + branches[3].1 + branches[4].1 - w_dd / (w_dt + w_dd)).abs() < 1e-12
        );
    }

    #[test]
    fn zero_tail_fraction_recovers_no_tail_kernel_bit_for_bit() {
        // Module-rustdoc anchor: η = 0 reproduces the no-tail mixture kernel
        // bit-for-bit — strengths, axis summary, sampled stream, and emitted
        // cards — for any T_tail spelling (even one whose reactivity would be
        // out of domain, since the tail evaluations are skipped).
        let mut base = flat_config(1.85, 0.0, 0.0);
        base.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        base.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        let mut zero = base;
        zero.tail = Some(DeuteriumTail::new(0.0, 1.0e6).unwrap());
        for &r in &[0.0, 37.5, 100.0, 150.0, 199.9, 200.0] {
            assert_eq!(
                zero.strength_density(r).unwrap(),
                base.strength_density(r).unwrap(),
                "tail-free strength must equal the no-tail kernel at r={r}"
            );
            assert_eq!(
                zero.proton_strength_density(r).unwrap(),
                base.proton_strength_density(r).unwrap(),
            );
        }
        assert_eq!(
            zero.axis_spectrum_summary().unwrap(),
            base.axis_spectrum_summary().unwrap(),
        );
        assert_eq!(
            zero.spectrum_branches(20.0, 30.0).unwrap(),
            base.spectrum_branches(20.0, 30.0).unwrap(),
        );
        let a = ParametricSampler::new(base, 17).unwrap().sample_n(512);
        let b = ParametricSampler::new(zero, 17).unwrap().sample_n(512);
        assert_eq!(a, b);
        assert_eq!(
            emission_histograms(&zero, 15).unwrap(),
            emission_histograms(&base, 15).unwrap()
        );
        // No tail drift row at η = 0: the emitted cards (text and drift) are
        // the landed ones.
        let cards_base = crate::emit_sdef_parametric(&base, 5, 15).unwrap();
        let cards_zero = crate::emit_sdef_parametric(&zero, 5, 15).unwrap();
        assert_eq!(cards_zero.text, cards_base.text);
        assert_eq!(cards_zero.drift, cards_base.drift);
        let serp_base = crate::emit_serpent_parametric(&base, 15).unwrap();
        let serp_zero = crate::emit_serpent_parametric(&zero, 15).unwrap();
        assert_eq!(serp_zero.text, serp_base.text);
        assert_eq!(serp_zero.drift, serp_base.drift);
    }

    #[test]
    fn tail_sampler_fires_all_sub_branches() {
        // The pinned tail config: all five sub-lines must fire in the stream
        // at the sub-rate-rule shares, and the mean energy must match the
        // five-branch weighted mean (module rustdoc).
        let config = tail_config();
        let mut sampler = ParametricSampler::new(config, 13).unwrap();
        let particles = sampler.sample_n(N);
        let eta = 0.05_f64;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap();
        let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(48.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(40.0).unwrap();
        let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(60.0).unwrap();
        let w = [
            0.7 * 0.3 * (1.0 - eta) * sv_dt,
            0.7 * 0.3 * eta * sv_dt_tail,
            0.7 * 0.7 / 2.0 * (1.0 - eta) * (1.0 - eta) * sv_dd,
            0.7 * 0.7 / 2.0 * 2.0 * eta * (1.0 - eta) * sv_dd_mix,
            0.7 * 0.7 / 2.0 * eta * eta * sv_dd_tail,
        ];
        let total: f64 = w.iter().sum();
        let p: Vec<f64> = w.iter().map(|x| x / total).collect();
        // D-D share (all three D-D sub-lines sit below 10 MeV).
        let p_dd = p[2] + p[3] + p[4];
        let below = particles.iter().filter(|p| p.energy_mev < 10.0).count() as f64;
        let frac_dd = below / N as f64;
        let se = (p_dd * (1.0 - p_dd) / N as f64).sqrt();
        assert!(
            (frac_dd - p_dd).abs() < 8.0 * se,
            "D-D branch fraction {frac_dd} vs {p_dd} ± {se}"
        );
        // Five-branch weighted mean and variance.
        let temps = [24.0, 48.0, 20.0, 40.0, 60.0];
        let reactions = [
            FusionReaction::Dt,
            FusionReaction::Dt,
            FusionReaction::Dd,
            FusionReaction::Dd,
            FusionReaction::Dd,
        ];
        let mut want = 0.0;
        let mut second = 0.0;
        for (i, (&reaction, &ti)) in reactions.iter().zip(temps.iter()).enumerate() {
            let (mu, sigma) = reaction.moments_mev(ti).unwrap();
            want += p[i] * mu;
            second += p[i] * (sigma * sigma + mu * mu);
        }
        let var = second - want * want;
        let mean: f64 = particles.iter().map(|p| p.energy_mev).sum::<f64>() / N as f64;
        let se_mean = var.sqrt() / (N as f64).sqrt();
        assert!(
            (mean - want).abs() < 6.0 * se_mean,
            "<E> {mean} vs {want} ± {se_mean}"
        );
        // Deterministic per seed on the tail path too.
        let a = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        let b = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        assert_eq!(a, b);
    }

    #[test]
    fn tail_emission_carries_tail_drift_row() {
        let config = tail_config();
        let hist = emission_histograms(&config, 61).unwrap();
        let energy_sum: f64 = hist.energy.masses.iter().sum();
        assert!((energy_sum - hist.energy_coverage).abs() < 1e-12);
        assert!(hist.energy_coverage > 0.999);
        // The tail D-T sub-line (48 keV: mean 14.10 MeV, sigma 0.52 MeV)
        // widens the tabulated window past the bulk line and lifts the
        // marginal mean: both are deterministic histogram properties.
        let mut base = config;
        base.tail = None;
        let hist_base = emission_histograms(&base, 61).unwrap();
        assert!(hist.energy.centers.last().unwrap() > hist_base.energy.centers.last().unwrap());
        // The tail boosts the D-D sub-rate more than the D-T line shifts up
        // (hot deuterons meet deuterons at T_mix/T_tail), so the below-10
        // MeV histogram mass grows while the window top extends past the
        // bulk line — both deterministic card properties.
        let mass_low = |h: &EmissionHistograms| {
            h.energy
                .masses
                .iter()
                .zip(&h.energy.centers)
                .filter(|(_, &c)| c < 10.0)
                .map(|(&m, _)| m)
                .sum::<f64>()
        };
        assert!(mass_low(&hist) > mass_low(&hist_base));
        // Both dialects carry the tail drift row with the integrated share.
        let sdef = crate::emit_sdef_parametric(&config, 5, 15).unwrap();
        sdef.verify_round_trip().unwrap();
        let quantities: Vec<_> = sdef
            .drift
            .rows
            .iter()
            .map(|r| r.quantity.as_str())
            .collect();
        assert_eq!(
            quantities,
            [
                "emission probability",
                "spatial marginals",
                "joint correlation",
                "deuterium tail"
            ]
        );
        let serpent = crate::emit_serpent_parametric(&config, 15).unwrap();
        assert_eq!(
            serpent.drift.rows.last().unwrap().quantity,
            "deuterium tail"
        );
        // The integrated tail share on the flat uniform plasma matches the
        // sub-rate hand ratio (grid-independent: every cell shares one
        // value): the tail sub-branches carry ~7.8% of the births while the
        // net rate increase is smaller (the bulk depletes).
        let share = config.total_tail_strength().unwrap() / config.total_strength().unwrap();
        let eta = 0.05_f64;
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap();
        let sv_dt_tail = FusionReaction::Dt.reactivity_m3_per_s(48.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd_mix = FusionReaction::Dd.reactivity_m3_per_s(40.0).unwrap();
        let sv_dd_tail = FusionReaction::Dd.reactivity_m3_per_s(60.0).unwrap();
        let bulk =
            0.7 * 0.3 * (1.0 - eta) * sv_dt + 0.7 * 0.7 / 2.0 * (1.0 - eta) * (1.0 - eta) * sv_dd;
        let tail = 0.7 * 0.3 * eta * sv_dt_tail
            + 0.7 * 0.7 / 2.0 * (2.0 * eta * (1.0 - eta) * sv_dd_mix + eta * eta * sv_dd_tail);
        let want_share = tail / (bulk + tail);
        assert!((share - want_share).abs() < 1e-9, "{share} vs {want_share}");
        assert!((share - 0.07798654309040301).abs() < 1e-9, "{share}");
    }

    #[test]
    fn tail_parameters_are_loud() {
        assert_eq!(
            DeuteriumTail::new(f64::NAN, 60.0),
            Err(Error::NonFinite("deuterium tail fraction"))
        );
        assert_eq!(
            DeuteriumTail::new(0.05, f64::INFINITY),
            Err(Error::NonFinite("deuterium tail temperature"))
        );
        assert_eq!(
            DeuteriumTail::new(-0.1, 60.0),
            Err(Error::InvalidTail(
                "deuterium tail fraction must be in [0, 1]"
            ))
        );
        assert_eq!(
            DeuteriumTail::new(1.1, 60.0),
            Err(Error::InvalidTail(
                "deuterium tail fraction must be in [0, 1]"
            ))
        );
        assert_eq!(
            DeuteriumTail::new(0.05, -1.0),
            Err(Error::NegativeIonTemperature(-1.0))
        );
        // Boundary fractions are admitted exactly.
        assert!(DeuteriumTail::new(0.0, 60.0).is_ok());
        assert!(DeuteriumTail::new(1.0, 60.0).is_ok());
        // Config-level validation surfaces the same named errors.
        let mut c = iter_h_mode();
        c.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        c.tail = Some(DeuteriumTail {
            fraction: 2.0,
            temperature_kev: 60.0,
        });
        assert!(matches!(c.validate(), Err(Error::InvalidTail(_))));
        let mut c = iter_h_mode();
        c.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        c.tail = Some(DeuteriumTail {
            fraction: f64::NAN,
            temperature_kev: 60.0,
        });
        assert!(matches!(c.validate(), Err(Error::NonFinite(_))));
        // The tail without a fuel mixture is a loud scope error, never a
        // silent single-fuel fallback.
        let mut single = iter_h_mode();
        single.tail = Some(DeuteriumTail::new(0.05, 60.0).unwrap());
        assert!(matches!(single.validate(), Err(Error::NotYetSupported(_))));
        assert!(matches!(
            ParametricSampler::new(single, 7),
            Err(Error::NotYetSupported(_))
        ));
    }

    // --- G-sector: toroidal-sector normalization, recovery, loud angles ---

    #[test]
    fn sector_strength_scales_by_rotation_over_two_pi() {
        // The kernel is φ-independent, so the toroidal integral is the
        // rotation itself: power-of-two fractions are exact in binary FP
        // (scaling by 1/2 commutes with the one rounding), general angles
        // at transcription precision.
        let full = flat_config(1.85, 0.0, 0.0);
        let total_full = full.total_strength().unwrap();
        let mut half = full;
        half.sector = Some(ToroidalSector::new(0.0, PI).unwrap());
        assert_eq!(half.total_strength().unwrap(), total_full * 0.5);
        let mut quarter = full;
        quarter.sector = Some(ToroidalSector::new(1.0, PI / 2.0).unwrap());
        assert_eq!(quarter.total_strength().unwrap(), total_full * 0.25);
        // Start shifts births, never the total.
        let mut shifted = full;
        shifted.sector = Some(ToroidalSector::new(2.0, PI).unwrap());
        assert_eq!(
            shifted.total_strength().unwrap(),
            half.total_strength().unwrap()
        );
        // General angle: ratio at 1e-12.
        let mut general = full;
        general.sector = Some(ToroidalSector::new(0.5, 1.0).unwrap());
        let ratio = general.total_strength().unwrap() / total_full;
        let want = 1.0 / (2.0 * PI);
        assert!((ratio - want).abs() < 1e-12 * want, "{ratio} vs {want}");
        // The local density is sector-independent (only the toroidal
        // integral rescales).
        assert_eq!(
            general.strength_density(123.0).unwrap(),
            full.strength_density(123.0).unwrap()
        );
    }

    #[test]
    fn full_rotation_sector_recovers_full_torus_bit_for_bit() {
        // Module-rustdoc anchor: start 0 + exact 2π reproduces the landed
        // kernel — total, sampled stream, and card marginals.
        let full = iter_h_mode();
        let mut sector = full;
        sector.sector = Some(ToroidalSector::new(0.0, 2.0 * PI).unwrap());
        assert_eq!(
            sector.total_strength().unwrap(),
            full.total_strength().unwrap()
        );
        let a = ParametricSampler::new(full, 17).unwrap().sample_n(256);
        let b = ParametricSampler::new(sector, 17).unwrap().sample_n(256);
        assert_eq!(a, b);
        assert_eq!(
            emission_histograms(&sector, 15).unwrap(),
            emission_histograms(&full, 15).unwrap()
        );
    }

    #[test]
    fn sector_sampler_birth_angles_are_uniform_over_sector() {
        // start 0.5 rad, rotation 1.5 rad: every birth lands in-sector and
        // the mean angle matches the uniform closed form.
        let mut config = flat_config(1.85, 0.0, 0.0);
        let (start, rotation) = (0.5, 1.5);
        config.sector = Some(ToroidalSector::new(start, rotation).unwrap());
        let mut sampler = ParametricSampler::new(config, 23).unwrap();
        let particles = sampler.sample_n(N);
        let mut sum = 0.0;
        for p in &particles {
            let mut phi = p.position_cm[1].atan2(p.position_cm[0]);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }
            assert!(
                phi >= start && phi < start + rotation,
                "birth angle {phi} outside [{start}, {})",
                start + rotation
            );
            sum += phi;
        }
        let mean = sum / N as f64;
        let want = start + rotation / 2.0;
        let se = rotation / (12.0 * N as f64).sqrt();
        assert!(
            (mean - want).abs() < 8.0 * se,
            "mean angle {mean} vs {want} ± {se}"
        );
    }

    #[test]
    fn sector_emission_histograms_carry_uniform_angle_marginal() {
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.sector = Some(ToroidalSector::new(0.5, 1.5).unwrap());
        let hist = emission_histograms(&config, 15).unwrap();
        let angle = hist.toroidal.as_ref().expect("partial sector marginal");
        assert_eq!(angle.masses.len(), 15);
        assert!(angle.masses.iter().all(|&m| m == 1.0 / 15.0));
        assert!(angle.centers.windows(2).all(|w| w[0] < w[1]));
        assert!((angle.centers[0] - (0.5 + 0.05)).abs() < 1e-12);
        assert!((angle.centers[14] - (0.5 + 1.5 - 0.05)).abs() < 1e-12);
        let sum: f64 = angle.masses.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "angle masses sum {sum}");
        // Full torus carries no marginal: the landed card shape is untouched.
        let full = emission_histograms(&iter_h_mode(), 15).unwrap();
        assert!(full.toroidal.is_none());
        let mut rotation = flat_config(1.85, 0.0, 0.0);
        rotation.sector = Some(ToroidalSector::new(1.0, 2.0 * PI).unwrap());
        assert!(emission_histograms(&rotation, 15)
            .unwrap()
            .toroidal
            .is_none());
    }

    #[test]
    fn sector_angles_are_loud() {
        assert_eq!(
            ToroidalSector::new(f64::NAN, 1.0),
            Err(Error::NonFinite("toroidal sector start angle"))
        );
        assert_eq!(
            ToroidalSector::new(0.0, f64::INFINITY),
            Err(Error::NonFinite("toroidal sector rotation angle"))
        );
        assert_eq!(
            ToroidalSector::new(-0.1, 1.0),
            Err(Error::InvalidSector("start_angle must be in [0, 2π)"))
        );
        assert_eq!(
            ToroidalSector::new(2.0 * PI, 1.0),
            Err(Error::InvalidSector("start_angle must be in [0, 2π)"))
        );
        assert_eq!(
            ToroidalSector::new(0.0, 0.0),
            Err(Error::InvalidSector("rotation_angle must be in (0, 2π]"))
        );
        assert_eq!(
            ToroidalSector::new(0.0, -1.0),
            Err(Error::InvalidSector("rotation_angle must be in (0, 2π]"))
        );
        assert_eq!(
            ToroidalSector::new(0.0, 2.0 * PI + 0.1),
            Err(Error::InvalidSector("rotation_angle must be in (0, 2π]"))
        );
        // Boundary values are admitted exactly.
        assert!(ToroidalSector::new(0.0, 2.0 * PI).is_ok());
        assert!(ToroidalSector::new(0.0, 1e-300).is_ok());
        // Config-level validation surfaces the same named errors.
        let mut c = iter_h_mode();
        c.sector = Some(ToroidalSector {
            start_angle: -1.0,
            rotation_angle: 1.0,
        });
        assert!(matches!(c.validate(), Err(Error::InvalidSector(_))));
        let mut c = iter_h_mode();
        c.sector = Some(ToroidalSector {
            start_angle: f64::NAN,
            rotation_angle: 1.0,
        });
        assert!(matches!(c.validate(), Err(Error::NonFinite(_))));
        // A bad sector is loud through total_strength too.
        let mut c = iter_h_mode();
        c.sector = Some(ToroidalSector {
            start_angle: 0.0,
            rotation_angle: f64::NAN,
        });
        assert!(c.total_strength().is_err());
    }

    // --- G-flat: sampled moments vs closed forms at δ = 0 ---------------------

    const N: usize = 200_000;

    #[test]
    fn flat_profile_sampled_moments_match_closed_form() {
        // δ = 0, esh = 0: J = κr, volume weight = R·κ·r.
        // <r²> = a²/2 (closed); <Z> = 0; <R> = R0 (flat strength).
        let config = flat_config(1.85, 0.0, 0.0);
        let mut sampler = ParametricSampler::new(config, 11).unwrap();
        let particles = sampler.sample_n(N);
        let a = config.geometry.minor_radius_cm;
        let r0 = config.geometry.major_radius_cm;
        let mut r2 = 0.0;
        let mut z_sum = 0.0;
        let mut big_r = 0.0;
        let kappa = config.geometry.elongation;
        for p in &particles {
            let x = p.position_cm[0];
            let y = p.position_cm[1];
            // Recover the sampled minor radius: at δ = esh = 0 the map
            // gives (R − R0) = r·cos θ and Z = κ·r·sin θ.
            let r = ((x.hypot(y) - r0).powi(2) + (p.position_cm[2] / kappa).powi(2)).sqrt();
            r2 += r * r;
            z_sum += p.position_cm[2];
            big_r += x.hypot(y);
        }
        let tol = 8.0 * a / (2.0_f64 * N as f64).sqrt();
        assert!(
            (r2 / N as f64 - a * a / 2.0).abs() < tol * a,
            "<r²> {}",
            r2 / N as f64
        );
        assert!(z_sum.abs() / (N as f64) < tol);
        // Birth major radius is volume-element weighted: <R> = ∫∫R²J/∫∫RJ
        // = R₀ + a²/(4R₀) exactly at δ = esh = 0 (outboard volume bias).
        let want_r = r0 + a * a / (4.0 * r0);
        assert!(
            (big_r / N as f64 - want_r).abs() < tol,
            "<R> {} vs {want_r}",
            big_r / N as f64
        );
    }

    #[test]
    fn dd_config_samples_dd_birth_energies() {
        // Regression: the sampler must draw from `config.fuel`, not a
        // hardcoded reaction — DT here would sit ~11.6 MeV too high.
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.fuel = FusionReaction::Dd;
        let mut sampler = ParametricSampler::new(config, 7).unwrap();
        let particles = sampler.sample_n(N);
        let (mu, sigma) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        let mean: f64 = particles.iter().map(|p| p.energy_mev).sum::<f64>() / N as f64;
        let se = sigma / (N as f64).sqrt();
        assert!((mean - mu).abs() < 6.0 * se, "<E> {mean} vs {mu} ± {se}");
    }

    #[test]
    fn mixture_sampler_fires_both_branches() {
        // 70/30 D/T flat plasma at T_i = 20 keV: the D-D branch carries
        // p_dd = w_dd/(w_dt+w_dd) ≈ 6.96e-3 of the birth weight (module
        // rustdoc hand vectors), so both lines must appear in the stream
        // and the mean energy must match the branch-weighted mean.
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        let mut sampler = ParametricSampler::new(config, 13).unwrap();
        let particles = sampler.sample_n(N);
        let (k_dt, k_dd) = FuelMixture::new(0.7, 0.3).unwrap().branch_weights();
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let p_dd = k_dd * sv_dd / (k_dt * sv_dt + k_dd * sv_dd);
        let below = particles.iter().filter(|p| p.energy_mev < 10.0).count() as f64;
        let frac_dd = below / N as f64;
        let se = (p_dd * (1.0 - p_dd) / N as f64).sqrt();
        assert!(
            (frac_dd - p_dd).abs() < 8.0 * se,
            "D-D branch fraction {frac_dd} vs {p_dd} ± {se}"
        );
        // Branch-weighted mean and two-branch standard deviation at the
        // shared T_i = 20 keV.
        let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(20.0).unwrap();
        let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        let want = p_dd * mu_dd + (1.0 - p_dd) * mu_dt;
        let var = p_dd * (sigma_dd * sigma_dd + mu_dd * mu_dd)
            + (1.0 - p_dd) * (sigma_dt * sigma_dt + mu_dt * mu_dt)
            - want * want;
        let mean: f64 = particles.iter().map(|p| p.energy_mev).sum::<f64>() / N as f64;
        let se_mean = var.sqrt() / (N as f64).sqrt();
        assert!(
            (mean - want).abs() < 6.0 * se_mean,
            "<E> {mean} vs {want} ± {se_mean}"
        );
    }

    #[test]
    fn species_sampler_fires_both_branches_at_pair_temperatures() {
        // 70/30 D/T flat plasma at T_D = 20 keV, T_T = 30 keV: the D-T
        // branch reacts at T_DT = 24 keV and the D-D branch at T_D, so both
        // lines must appear in the stream and the mean energy must match
        // the pair-temperature branch-weighted mean (module rustdoc).
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        config.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 30.0).unwrap());
        let mut sampler = ParametricSampler::new(config, 13).unwrap();
        let particles = sampler.sample_n(N);
        let (k_dt, k_dd) = FuelMixture::new(0.7, 0.3).unwrap().branch_weights();
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(24.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let p_dd = k_dd * sv_dd / (k_dt * sv_dt + k_dd * sv_dd);
        let below = particles.iter().filter(|p| p.energy_mev < 10.0).count() as f64;
        let frac_dd = below / N as f64;
        let se = (p_dd * (1.0 - p_dd) / N as f64).sqrt();
        assert!(
            (frac_dd - p_dd).abs() < 8.0 * se,
            "D-D branch fraction {frac_dd} vs {p_dd} ± {se}"
        );
        let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(24.0).unwrap();
        let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        let want = p_dd * mu_dd + (1.0 - p_dd) * mu_dt;
        let var = p_dd * (sigma_dd * sigma_dd + mu_dd * mu_dd)
            + (1.0 - p_dd) * (sigma_dt * sigma_dt + mu_dt * mu_dt)
            - want * want;
        let mean: f64 = particles.iter().map(|p| p.energy_mev).sum::<f64>() / N as f64;
        let se_mean = var.sqrt() / (N as f64).sqrt();
        assert!(
            (mean - want).abs() < 6.0 * se_mean,
            "<E> {mean} vs {want} ± {se_mean}"
        );
        // Deterministic per seed on the species path too.
        let a = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        let b = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        assert_eq!(a, b);
    }

    #[test]
    fn mixture_sampling_is_deterministic_and_loud_when_inert() {
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.fuel_mixture = Some(FuelMixture::new(0.5, 0.5).unwrap());
        let a = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        let b = ParametricSampler::new(config, 31).unwrap().sample_n(256);
        assert_eq!(a, b);
        // Pure tritium models no reactions (no T-T branch): the strength
        // integral is zero, so construction is loud — never a dead stream.
        let mut inert = flat_config(1.85, 0.0, 0.0);
        inert.fuel_mixture = Some(FuelMixture::new(0.0, 1.0).unwrap());
        assert!(matches!(
            ParametricSampler::new(inert, 31),
            Err(Error::InvalidProfile(_))
        ));
    }

    #[test]
    fn mixture_emission_histograms_cover_both_lines() {
        let mut config = flat_config(1.85, 0.0, 0.0);
        config.fuel_mixture = Some(FuelMixture::new(0.7, 0.3).unwrap());
        let hist = emission_histograms(&config, 61).unwrap();
        let energy_sum: f64 = hist.energy.masses.iter().sum();
        assert!(
            (energy_sum - hist.energy_coverage).abs() < 1e-12,
            "energy {energy_sum}"
        );
        // Both lines carry mass: the below-10 MeV share matches the D-D
        // branch weight of the rate rule (histogram discretization only).
        let (k_dt, k_dd) = FuelMixture::new(0.7, 0.3).unwrap().branch_weights();
        let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(20.0).unwrap();
        let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(20.0).unwrap();
        let p_dd = k_dd * sv_dd / (k_dt * sv_dt + k_dd * sv_dd);
        let mass_low: f64 = hist
            .energy
            .masses
            .iter()
            .zip(&hist.energy.centers)
            .filter(|(_, &c)| c < 10.0)
            .map(|(&m, _)| m)
            .sum();
        assert!(
            (mass_low / energy_sum - p_dd).abs() < 0.02,
            "D-D histogram share {} vs {p_dd}",
            mass_low / energy_sum
        );
        // The mixture axis summary carries between-branch variance: sigma
        // far above either single-line Ballabio width, and the mean sits
        // between the two line means (closer to the dominant D-T line).
        let (_, mean, sigma) = config.axis_spectrum_summary().unwrap();
        let (mu_dt, _) = FusionReaction::Dt.moments_mev(20.0).unwrap();
        let (mu_dd, _) = FusionReaction::Dd.moments_mev(20.0).unwrap();
        assert!(mu_dd < mean && mean < mu_dt, "mixture mean {mean}");
        assert!((mean - mu_dt).abs() < (mean - mu_dd).abs());
        assert!(sigma > 0.5, "mixture sigma {sigma}");
        assert!(hist.energy_coverage > 0.999);
    }

    #[test]
    fn sampled_mean_energy_matches_strength_weighted_temperature() {
        // Global mean birth energy ≈ ⟨μ(T_i(r))⟩ weighted by S·R·J —
        // computed here by independent fine quadrature.
        let config = iter_h_mode();
        let mut sampler = ParametricSampler::new(config, 5).unwrap();
        let particles = sampler.sample_n(N);
        let e_mean: f64 = particles.iter().map(|p| p.energy_mev).sum::<f64>() / N as f64;
        let quad = mean_energy_quadrature(&config);
        let sigma_max = 0.4;
        let tol = 8.0 * sigma_max / (N as f64).sqrt();
        assert!(
            (e_mean - quad).abs() < tol,
            "mean E {e_mean} vs quadrature {quad}"
        );
    }

    /// Independent fine-quadrature reference for the strength-weighted mean
    /// Ballabio mean energy.
    fn mean_energy_quadrature(config: &ParametricPlasmaConfig) -> f64 {
        let g = &config.geometry;
        let n_r = 400;
        let n_t = 400;
        let dr = g.minor_radius_cm / n_r as f64;
        let dt = 2.0 * PI / n_t as f64;
        let mut num = 0.0;
        let mut den = 0.0;
        for i in 0..n_r {
            let r = (i as f64 + 0.5) * dr;
            let s = config.strength_density(r).unwrap();
            let ti = config.temperature_kev(r);
            let (mu, _) = config.fuel.moments_mev(ti).unwrap();
            for j in 0..n_t {
                let theta = (j as f64 + 0.5) * dt;
                let w = s * g.volume_element(r, theta) * dr * dt;
                num += w * mu;
                den += w;
            }
        }
        num / den
    }

    #[test]
    fn determinism_and_cell_consistency() {
        let config = iter_h_mode();
        let a = ParametricSampler::new(config, 99).unwrap().sample_n(64);
        let b = ParametricSampler::new(config, 99).unwrap().sample_n(64);
        assert_eq!(a, b);
        let c = ParametricSampler::new(config, 98).unwrap().sample_n(64);
        assert_ne!(a, c);
    }

    // --- G-hist: emission histograms -------------------------------------------

    #[test]
    fn emission_histograms_are_normalized_and_consistent() {
        let config = iter_h_mode();
        let hist = emission_histograms(&config, 21).unwrap();
        for dist in [&hist.radial, &hist.vertical] {
            let sum: f64 = dist.masses.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "masses sum {sum}");
            assert!(dist.centers.windows(2).all(|w| w[0] < w[1]));
            assert!(dist.masses.iter().all(|&m| m >= 0.0));
        }
        // The energy marginal keeps the tabulation truncation: its masses
        // sum to the captured coverage (same convention as the ring/point
        // spectrum tabulation), reported as drift.
        let energy_sum: f64 = hist.energy.masses.iter().sum();
        assert!(
            (energy_sum - hist.energy_coverage).abs() < 1e-12,
            "energy {energy_sum}"
        );
        assert!(hist.energy.centers.windows(2).all(|w| w[0] < w[1]));
        // Radial mass concentrates in the core (peaked H-mode profiles).
        let core_mass: f64 = hist.radial.masses[..7].iter().sum();
        assert!(core_mass > 0.3, "core mass {core_mass}");
        // Vertical span matches the geometry: centers within ±κa.
        let z_max = config.geometry.z_max();
        assert!(hist.vertical.centers.iter().all(|&z| z.abs() <= z_max));
        // Energy spectrum covers the DT line with truncation drift only.
        assert!(hist.energy_coverage > 0.999);
        // Joint correlation distance is a probability metric in [0, 1).
        assert!(hist.joint_correlation >= 0.0 && hist.joint_correlation < 1.0);
    }

    #[test]
    fn invalid_configs_are_loud() {
        let mut c = iter_h_mode();
        c.pedestal_radius_cm = 250.0;
        assert!(matches!(c.validate(), Err(Error::InvalidProfile(_))));
        let mut c = iter_h_mode();
        c.geometry.triangularity = 1.5;
        assert!(matches!(c.validate(), Err(Error::InvalidGeometry(_))));
        let mut c = iter_h_mode();
        c.ion_density.centre_m3 = -1.0;
        assert!(matches!(c.validate(), Err(Error::InvalidProfile(_))));
    }

    #[test]
    fn zero_poloidal_cell_weight_is_loud_at_build() {
        // A zero minor radius collapses every cell's poloidal Jacobian sum
        // (`w_acc == 0`), the 0/0 CDF division behind the sampler's
        // partial_cmp().unwrap() panic path. WeightTable::build is reached
        // without validation (e.g. via total_strength), so the guard must
        // fire on its own.
        let mut c = iter_h_mode();
        c.geometry.minor_radius_cm = 0.0;
        let err = WeightTable::build(&c).unwrap_err();
        assert!(
            matches!(err, Error::InvalidProfile(m) if m.contains("zero poloidal weight")),
            "{err}"
        );
        // Same via the public total_strength entry point.
        assert!(
            matches!(c.total_strength(), Err(Error::InvalidProfile(_))),
            "{:?}",
            c.total_strength()
        );
    }
}
