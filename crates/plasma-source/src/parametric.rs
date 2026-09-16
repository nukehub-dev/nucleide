//! Parametric tokamak plasma source: caller-supplied profiles over Miller
//! geometry, reactivity-weighted neutron emission, seeded sampling, and
//! histogram source-card emission.
//!
//! This is the second Cycle-01 landing. The source model is the public
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
//! Profiles are caller inputs; nothing here computes profiles or solves an
//! equilibrium. The documented loud boundary ([`Error::NotYetSupported`]):
//! reactant distributions beyond the shared-temperature Maxwellian mixture
//! above (per-species ion temperatures, non-Maxwellian tails — the full
//! Eriksson generalization), the T-T and D(d,p)T branches, and toroidal
//! sectors (caller-side rejection of the sampled `φ` remains the spelling).
//! Mixture fractions themselves are validated loudly at construction:
//! non-finite, negative, or non-summing pairs never reach the sampler.
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

    /// (D-T, D-D) branch strength densities \[neutrons/s/cm³, arbitrary
    /// global scale\] at minor radius `r` \[cm\]: each branch `k·n²·⟨σv⟩`
    /// with density in cm⁻³ and reactivity in cm³/s, built in the landed
    /// left-associated expression tree so that single-fuel configurations
    /// and the module-rustdoc recovery anchors reproduce bit-for-bit. Zero
    /// temperature ⇒ zero strength on both branches (cold separatrix makes
    /// no neutrons).
    fn branch_strength_density(&self, r: f64) -> Result<(f64, f64)> {
        let n_m3 = self.density_m3(r);
        let ti_kev = self.temperature_kev(r);
        let n_cm3 = n_m3 * 1e-6;
        match &self.fuel_mixture {
            Some(mixture) => {
                let (k_dt, k_dd) = mixture.branch_weights();
                let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(ti_kev)? * 1e6;
                let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(ti_kev)? * 1e6;
                Ok((k_dt * n_cm3 * n_cm3 * sv_dt, k_dd * n_cm3 * n_cm3 * sv_dd))
            }
            None => {
                let reactivity = self.fuel.reactivity_m3_per_s(ti_kev)? * 1e6;
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

    /// Neutron-branch spectra at ion temperature `ti_kev` \[keV\]:
    /// `(reaction, weight, mean, sigma)` tuples with weights summing to 1.
    /// Single fuel → one unit-weight branch; a mixture → D-T and D-D
    /// branches weighted by the rate rule at `ti_kev` (the module rustdoc
    /// normalization). Where the total mixture rate underflows — both
    /// reactivities zero, a cold annulus — the D-T branch is returned by
    /// convention: the mixture spectrum degenerates to the 14.021 MeV line
    /// exactly as the single-fuel kernels degenerate at `T_i = 0`.
    fn spectrum_branches(&self, ti_kev: f64) -> Result<Vec<(FusionReaction, f64, f64, f64)>> {
        match &self.fuel_mixture {
            None => {
                let (mu, sigma) = self.fuel.moments_mev(ti_kev)?;
                Ok(vec![(self.fuel, 1.0, mu, sigma)])
            }
            Some(mixture) => {
                let (k_dt, k_dd) = mixture.branch_weights();
                let sv_dt = FusionReaction::Dt.reactivity_m3_per_s(ti_kev)?;
                let sv_dd = FusionReaction::Dd.reactivity_m3_per_s(ti_kev)?;
                let w_dt = k_dt * sv_dt;
                let w_dd = k_dd * sv_dd;
                let (mu_dt, sigma_dt) = FusionReaction::Dt.moments_mev(ti_kev)?;
                if w_dt + w_dd <= 0.0 {
                    return Ok(vec![(FusionReaction::Dt, 1.0, mu_dt, sigma_dt)]);
                }
                let (mu_dd, sigma_dd) = FusionReaction::Dd.moments_mev(ti_kev)?;
                let p_dt = w_dt / (w_dt + w_dd);
                Ok(vec![
                    (FusionReaction::Dt, p_dt, mu_dt, sigma_dt),
                    (FusionReaction::Dd, 1.0 - p_dt, mu_dd, sigma_dd),
                ])
            }
        }
    }

    /// Spectrum summary at the magnetic axis (`r = 0`): `(nominal, mean,
    /// sigma)` \[MeV\]. Single fuel → the reaction's nominal line and
    /// Ballabio moments; a mixture → the two-branch Gaussian-mixture
    /// moments (between-branch variance included) and the dominant
    /// branch's nominal line (D-T preferred on exact ties).
    pub fn axis_spectrum_summary(&self) -> Result<(f64, f64, f64)> {
        let ti_kev = self.temperature_kev(0.0);
        match &self.fuel_mixture {
            None => {
                let (mean, sigma) = self.fuel.moments_mev(ti_kev)?;
                Ok((self.fuel.nominal_energy_mev(), mean, sigma))
            }
            Some(_) => {
                let branches = self.spectrum_branches(ti_kev)?;
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

    /// Total relative source strength ∭ S·R·|J| da dθ dφ (the toroidal
    /// integral contributes the factor `2π`). Arbitrary global scale, but
    /// ratios between configurations are meaningful.
    pub fn total_strength(&self) -> Result<f64> {
        let table = WeightTable::build(self)?;
        Ok(2.0 * PI * table.masses.iter().sum::<f64>())
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
            let strength = config.strength_density(r)?;
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
    /// from the cell's poloidal CDF, `φ` uniform; energy from the local
    /// ion-temperature Ballabio Gaussian; isotropic direction.
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

    /// Draw the birth spectrum moments `(mean, sigma)` \[MeV\] at ion
    /// temperature `ti_kev`: the config fuel's Ballabio moments, or — for a
    /// fuel mixture — the branch roulette weighted by the rate rule at
    /// `ti_kev` (module rustdoc) followed by the chosen branch's moments.
    /// The single-fuel path consumes no RNG, so the landed sampling stream
    /// is preserved bit-for-bit; the mixture path draws one extra uniform
    /// per particle. Both reactivities zero (a cold annulus) degenerates to
    /// the D-T line by the documented convention of
    /// [`ParametricPlasmaConfig::spectrum_branches`].
    fn sample_spectrum_moments(&mut self, ti_kev: f64) -> Result<(f64, f64)> {
        let branches = self.config.spectrum_branches(ti_kev)?;
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
        let phi = 2.0 * PI * self.rng.uniform();
        let (big_r, z) = self.config.geometry.map(r, theta);
        let ti_kev = self.config.temperature_kev(r);
        let (mean, sigma) = self.sample_spectrum_moments(ti_kev)?;
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
        let ti = config.temperature_kev(r);
        for &(_, _, mu, sigma) in &config.spectrum_branches(ti)? {
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
        let ti = config.temperature_kev(r);
        let branches = config.spectrum_branches(ti)?;
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
