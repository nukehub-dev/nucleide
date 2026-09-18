//! Arbitrary-3D birth-rate lattice source: caller-supplied point clouds with
//! field-period symmetry reduction, a seeded sampler, and histogram
//! source-card emission.
//!
//! Every tokamak assumption is lifted
//! except the physics kernel itself. The caller supplies the full 3D birth
//! distribution as an explicit point list — positions × relative birth rates ×
//! per-point ion temperatures — and the landed machinery (D/T mixture rates,
//! Ballabio spectra, seeded sampling, card emission with drift reports, the
//! MCPL caller-side rule) applies per point.
//!
//! # Lattice contract
//!
//! A [`LatticePoint`] is one node of the caller cloud:
//!
//! ```text
//! position_cm: [x, y, z]  [cm]   (arbitrary 3D, no axisymmetry assumed)
//! rate:        w >= 0            (relative birth-rate weight, dimensionless;
//!                                arbitrary global scale — only ratios matter)
//! ion_temperature_kev: Ti >= 0   [keV] (drives the node's Ballabio line;
//!                                unused for rate/spectrum when a uniform
//!                                species-temperature pair is set, still
//!                                validated)
//! ```
//!
//! The lattice is a point list, not a regular grid: stellarator birth clouds
//! are rarely Cartesian, and a point list carries the Stellaris precedent
//! (30 radial × 50 poloidal × 100 toroidal point sources, 14.06 MeV
//! monoenergetic — here the monoenergetic limit is `Ti = 0` on the D-T line,
//! nominal 14.021 MeV) without imposing a topology. Per-point spectra are the
//! landed Ballabio Gaussians at the node's ion temperature; a [`FuelMixture`]
//! on the config activates both branches through the same Eriksson/DRESS
//! rate rule as the parametric source (per-point branch roulette), and a
//! uniform [`SpeciesIonTemperatures`] pair reacts every node at
//! `(T_D, T_T)` exactly like the parametric pair rule (equal-T recovery is
//! bit-for-bit, see below).
//!
//! Hand vectors (D-T, `Ti = 0`, nominal 14.021 MeV):
//!
//! ```text
//! A = (300, 0, 25) rate 2, B = (−300, 0, 25) rate 1: total = 3.0 exactly;
//! mean birth position = (100, 0, 25); mean energy = 14.021 MeV exactly.
//! ```
//!
//! # Symmetry reduction
//!
//! [`LatticeSymmetry`] declares field-period symmetry with `field_periods =
//! n ≥ 1` and a `base_angle` (radians, machine `+x` toward `+y`, the sampler
//! azimuth convention). The config points are then the BASE-SECTOR cloud
//! (one field period); toroidal angles identify modulo `period = 2π/n` about
//! `base_angle`, and physical rates accumulate by orbit mean (exact `f64`
//! orbit sum divided by the period count — the per-sector rate, so the
//! symmetric total `n·Σbase` reproduces the full-cloud total). [`expand_lattice`] replicates the base cloud orbit-major
//! (`copy 0` is the untouched base point, bit-for-bit; copy `k` rotates by
//! `k·period` through the shared [`rotate_xy`] used by the sampler), and
//! [`fold_lattice`] inverts it on orbit-major full clouds (chunked
//! consecutive groups of `n`, rates summed by left fold, `Ti` equality
//! enforced exactly — a non-symmetric input is a loud
//! [`Error::InvalidSymmetry`], never a silently merged cloud).
//!
//! Bit-identity gate (not distributionICTRIBUTion identity): sampling the
//! folded base cloud with symmetry set reproduces sampling the expanded full
//! cloud without symmetry bit-for-bit (same seed), because the sampler picks
//! the field-period copy from the within-bin fraction of the single index
//! draw — no extra RNG draws on any path (a lone base node draws the copy
//! fraction directly, which is the same single selection draw the expanded
//! full cloud spends on its index) — and copy-0 rotation is the
//! identity. With unit rates and power-of-two periods the bin-boundary
//! arithmetic is exact (division/multiplication by powers of two); general
//! rates agree up to the usual 1-ulp boundary adjacency, so the pinned gate
//! uses exact-friendly values. `n = 1` is the identity (fold/expand are
//! bit-for-bit no-ops); the sampler total is `n·Σrates` with symmetry,
//! `Σrates` without.
//!
//! # Recovery anchors
//!
//! - Single-point lattice ≡ point source: one node with `Ti` reproduces the
//!   landed [`crate::SourceSampler`] point stream bit-for-bit (the sampler
//!   skips the index draw for a single node without symmetry, so draw
//!   accounting matches the point path exactly; single-fuel only — a mixture
//!   adds the branch roulette draw on both paths identically only through
//!   the lattice sampler, so the anchor is pinned single-fuel).
//! - Axisymmetric limit converges (moments, not streams — different sampler
//!   paths): a uniform ring lattice reproduces the analytic ring moments
//!   (`<R>`, `z`, Ballabio mean) to sampling tolerance, and the total equals
//!   the summed rates exactly.
//! - Equal species pair `(T, T)` reproduces the shared-temperature lattice
//!   kernel bit-for-bit (stream and cards), the parametric spelling.
//!
//! # Card emission
//!
//! A discrete point cloud has no spelling in the typed `SDEF` subset (and no
//! Serpent point-list reader exists in the workspace), so emission follows
//! the parametric precedent: product-form *marginals* — cylindrical-`R`
//! histogram over `[0, R_max]`, vertical histogram over the cloud `z` span,
//! global rate-weighted marginal energy spectrum — with the same drift rows
//! (energy coverage, spatial marginals, joint correlation) plus a `lattice`
//! row (node count, total rate, symmetry spelling). The `(R, z, E)`
//! correlation the card drops is quantified, never hidden.
//!
//! Loud boundary ([`Error::NotYetSupported`]): everything the parametric
//! source defers (T-T neutron transport, proton transport, non-Maxwellian
//! tails beyond the parametric deuterium hot-tail fraction —
//! [`crate::parametric::DeuteriumTail`] is parametric-only, there is no tail
//! spelling on lattice configs — equilibrium solving, CAD/DAGMC) stays
//! deferred here — the per-node
//! proton *rate* is accounted through the D-D branch share like the
//! parametric bookkeeping while sampler and cards stay neutron-only. No HDF5,
//! no vendored data, no equilibrium-file consumption.
//!
//! Units: cm, MeV, keV (card convention); rates dimensionless relative
//! weights; reactivity exposed in m³/s where mixtures need it.

use std::f64::consts::PI;

use crate::parametric::{FuelMixture, SpeciesIonTemperatures};
use crate::sample::{Particle, Rng};
use crate::{BinnedDistribution, Error, FusionReaction, Result};

/// One node of the caller-supplied birth-rate cloud.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticePoint {
    /// Birth position `[x, y, z]` \[cm\]; each finite.
    pub position_cm: [f64; 3],
    /// Relative birth-rate weight (`>= 0`, finite; arbitrary global scale).
    pub rate: f64,
    /// Node ion temperature \[keV\] driving the Ballabio line (`>= 0`,
    /// finite; `0` is the monoenergetic nominal line).
    pub ion_temperature_kev: f64,
}

impl LatticePoint {
    /// Validate one node; errors are loud and named.
    pub fn validate(&self) -> Result<()> {
        for (field, value) in [
            ("lattice point x", self.position_cm[0]),
            ("lattice point y", self.position_cm[1]),
            ("lattice point z", self.position_cm[2]),
        ] {
            if !value.is_finite() {
                return Err(Error::NonFinite(field));
            }
        }
        if !self.rate.is_finite() {
            return Err(Error::NonFinite("lattice point rate"));
        }
        if self.rate < 0.0 {
            return Err(Error::InvalidLattice("lattice point rate must be >= 0"));
        }
        if !self.ion_temperature_kev.is_finite() {
            return Err(Error::NonFinite("lattice point ion temperature"));
        }
        if self.ion_temperature_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(self.ion_temperature_kev));
        }
        Ok(())
    }
}

/// Field-period symmetry declaration for a lattice source: the config points
/// are the base-sector cloud (one field period spanning
/// `[base_angle, base_angle + period)`), physically replicated `field_periods`
/// times around the machine axis (module rustdoc).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticeSymmetry {
    /// Field-period count `n ≥ 1` (`1` is the identity: no replication).
    pub field_periods: u32,
    /// Base-sector start angle \[rad\] (finite; measured from machine `+x`
    /// toward `+y`).
    pub base_angle: f64,
}

impl LatticeSymmetry {
    /// New validated symmetry; errors are loud and named.
    pub fn new(field_periods: u32, base_angle: f64) -> Result<Self> {
        let symmetry = Self {
            field_periods,
            base_angle,
        };
        symmetry.validate()?;
        Ok(symmetry)
    }

    /// Validate: `field_periods ≥ 1`, finite base angle.
    pub fn validate(&self) -> Result<()> {
        if self.field_periods < 1 {
            return Err(Error::InvalidSymmetry("field periods must be >= 1"));
        }
        if !self.base_angle.is_finite() {
            return Err(Error::NonFinite("lattice symmetry base angle"));
        }
        Ok(())
    }

    /// Field period `2π/n` \[rad\].
    pub fn period(&self) -> f64 {
        2.0 * PI / self.field_periods as f64
    }

    /// True for `field_periods == 1` (fold/expand are no-ops, totals equal
    /// the bare rate sum).
    pub fn is_trivial(&self) -> bool {
        self.field_periods == 1
    }

    /// Total-strength scale factor (`field_periods` as `f64`).
    pub fn copies(&self) -> f64 {
        self.field_periods as f64
    }
}

/// Rotate an `(x, y)` pair about the machine `z` axis by `angle` \[rad\].
/// Shared by [`expand_lattice`] and the sampler so expanded clouds and
/// symmetry-replicated draws agree bit-for-bit.
pub fn rotate_xy(x: f64, y: f64, angle: f64) -> (f64, f64) {
    let (sin, cos) = angle.sin_cos();
    (x * cos - y * sin, x * sin + y * cos)
}

/// Arbitrary-3D birth-rate lattice source configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeSourceConfig {
    /// Caller cloud: the full sampling distribution without symmetry, or the
    /// base-sector cloud with symmetry (module rustdoc).
    pub points: Vec<LatticePoint>,
    /// Fuel reaction: D-T or D-D. Used only when [`Self::fuel_mixture`] is
    /// `None`; a mixture activates both branches through
    /// [`FuelMixture::branch_weights`] and this field plays no role.
    pub reaction: FusionReaction,
    /// Optional D/T fuel mixture (the parametric Eriksson/DRESS rule,
    /// evaluated per node at the node temperatures).
    pub fuel_mixture: Option<FuelMixture>,
    /// Optional uniform per-species ion temperatures (the parametric pair
    /// rule — every node reacts at `(T_D, T_T)`; node temperatures are then
    /// unused for rate and spectrum, still validated).
    pub species_temperatures: Option<SpeciesIonTemperatures>,
    /// Optional field-period symmetry (module rustdoc).
    pub symmetry: Option<LatticeSymmetry>,
    /// Particle weight carried by the sampler and the emitted cards.
    pub weight: f64,
}

impl LatticeSourceConfig {
    /// Validate the full configuration (fields only; the positive-total-rate
    /// precondition is enforced at sampler/table build like the parametric
    /// zero-strength gate).
    pub fn validate(&self) -> Result<()> {
        if self.points.is_empty() {
            return Err(Error::InvalidLattice("lattice needs at least one point"));
        }
        for point in &self.points {
            point.validate()?;
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
        if let Some(symmetry) = &self.symmetry {
            symmetry.validate()?;
        }
        if !self.weight.is_finite() || self.weight <= 0.0 {
            return Err(Error::NonPositiveWeight(self.weight));
        }
        Ok(())
    }

    /// Species ion temperatures `(T_D, T_T)` \[keV\] at one node: the uniform
    /// pair when set, else the node temperature twice.
    fn species_temperatures_at(&self, node_ti_kev: f64) -> (f64, f64) {
        match &self.species_temperatures {
            Some(pair) => (pair.deuterium_kev, pair.tritium_kev),
            None => (node_ti_kev, node_ti_kev),
        }
    }

    /// Bare rate sum `Σrates` over the config points (finite, `>= 0`).
    pub fn base_total_rate(&self) -> Result<f64> {
        self.validate()?;
        let total: f64 = self.points.iter().map(|p| p.rate).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(Error::InvalidLattice(
                "total lattice birth rate is zero (check point rates)",
            ));
        }
        Ok(total)
    }

    /// Total relative source strength: the bare rate sum, scaled by the
    /// field-period count when symmetry is set (module rustdoc).
    pub fn total_strength(&self) -> Result<f64> {
        let base = self.base_total_rate()?;
        match &self.symmetry {
            None => Ok(base),
            Some(symmetry) => Ok(symmetry.copies() * base),
        }
    }

    /// Neutron-branch spectra at species temperatures `(t_d, t_t)` \[keV\]:
    /// `(reaction, weight, mean, sigma)` tuples with weights summing to 1 —
    /// the same rate rule as the parametric kernel (single fuel → one
    /// unit-weight branch; mixture → D-T at the pair effective temperature
    /// and D-D at `t_d`, degenerating to the D-T line where both
    /// reactivities vanish).
    pub fn spectrum_branches(
        &self,
        t_d: f64,
        t_t: f64,
    ) -> Result<Vec<(FusionReaction, f64, f64, f64)>> {
        match &self.fuel_mixture {
            None => {
                let (mu, sigma) = self.reaction.moments_mev(t_d)?;
                Ok(vec![(self.reaction, 1.0, mu, sigma)])
            }
            Some(mixture) => {
                let t_eff = SpeciesIonTemperatures::dt_effective_kev(t_d, t_t);
                let (k_dt, k_dd) = mixture.branch_weights();
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
                Ok(vec![
                    (FusionReaction::Dt, p_dt, mu_dt, sigma_dt),
                    (FusionReaction::Dd, 1.0 - p_dt, mu_dd, sigma_dd),
                ])
            }
        }
    }

    /// D-D neutron-branch share of one node's rate (the proton-accounting
    /// weight): `1` for single-fuel D-D, `0` for single-fuel D-T, the
    /// rate-rule D-D fraction for a mixture (cold nodes contribute `0`).
    fn dd_branch_share(&self, node_ti_kev: f64) -> Result<f64> {
        match &self.fuel_mixture {
            None => Ok(match self.reaction {
                FusionReaction::Dd => 1.0,
                FusionReaction::Dt => 0.0,
            }),
            Some(_) => {
                let (t_d, t_t) = self.species_temperatures_at(node_ti_kev);
                let branches = self.spectrum_branches(t_d, t_t)?;
                if branches.len() == 1 {
                    return Ok(0.0);
                }
                Ok(1.0 - branches[0].1)
            }
        }
    }

    /// Total relative D-D neutron-branch strength (the proton bookkeeping
    /// total — the pinned 50/50 convention: the D(d,p)T proton rate IS this
    /// density). Symmetry-aware like [`Self::total_strength`].
    pub fn total_dd_strength(&self) -> Result<f64> {
        self.validate()?;
        let mut total = 0.0;
        for point in &self.points {
            total += point.rate * self.dd_branch_share(point.ion_temperature_kev)?;
        }
        if !total.is_finite() {
            return Err(Error::InvalidLattice("lattice D-D branch total is empty"));
        }
        match &self.symmetry {
            None => Ok(total),
            Some(symmetry) => Ok(symmetry.copies() * total),
        }
    }

    /// Rate-weighted mean birth energy \[MeV\] over the cloud (branch
    /// mixtures included): the analytic moment the sampler must reproduce.
    pub fn mean_birth_energy(&self) -> Result<f64> {
        self.validate()?;
        let mut num = 0.0;
        let mut den = 0.0;
        for point in &self.points {
            let (t_d, t_t) = self.species_temperatures_at(point.ion_temperature_kev);
            let branches = self.spectrum_branches(t_d, t_t)?;
            let mean: f64 = branches.iter().map(|b| b.1 * b.2).sum();
            num += point.rate * mean;
            den += point.rate;
        }
        if den <= 0.0 || !den.is_finite() {
            return Err(Error::InvalidLattice(
                "total lattice birth rate is zero (check point rates)",
            ));
        }
        Ok(num / den)
    }

    /// Rate-weighted mean birth position \[cm\] over the config points (base
    /// sector when symmetry is set — the full-cloud mean follows by
    /// replication symmetry).
    pub fn mean_birth_position(&self) -> Result<[f64; 3]> {
        self.validate()?;
        let mut num = [0.0f64; 3];
        let mut den = 0.0;
        for point in &self.points {
            for (k, n) in num.iter_mut().enumerate() {
                *n += point.rate * point.position_cm[k];
            }
            den += point.rate;
        }
        if den <= 0.0 || !den.is_finite() {
            return Err(Error::InvalidLattice(
                "total lattice birth rate is zero (check point rates)",
            ));
        }
        Ok([num[0] / den, num[1] / den, num[2] / den])
    }
}

/// Expand a base-sector cloud to the full orbit-major lattice: `n` copies per
/// point (copy `0` first, bit-for-bit the base point; copy `k` rotates the
/// base position by `k·period` through [`rotate_xy`]).
pub fn expand_lattice(base: &[LatticePoint], symmetry: &LatticeSymmetry) -> Vec<LatticePoint> {
    let n = symmetry.field_periods as usize;
    let period = symmetry.period();
    let mut full = Vec::with_capacity(base.len() * n);
    for point in base {
        for k in 0..n {
            if k == 0 {
                full.push(*point);
            } else {
                let (x, y) = rotate_xy(
                    point.position_cm[0],
                    point.position_cm[1],
                    k as f64 * period,
                );
                full.push(LatticePoint {
                    position_cm: [x, y, point.position_cm[2]],
                    rate: point.rate,
                    ion_temperature_kev: point.ion_temperature_kev,
                });
            }
        }
    }
    full
}

/// Fold an orbit-major full lattice back to the base-sector cloud: chunked
/// consecutive groups of `n` (the [`expand_lattice`] order), rates accumulated
/// by orbit mean (exact `f64` sum divided by the period count — the per-sector
/// rate, so the symmetric total `n·Σbase` reproduces the full-cloud total),
/// first position kept. `Ti` must match exactly across each orbit — a
/// non-symmetric input is a loud [`Error::InvalidSymmetry`], as is a point
/// count that is not a multiple of `n`.
pub fn fold_lattice(
    full: &[LatticePoint],
    symmetry: &LatticeSymmetry,
) -> Result<Vec<LatticePoint>> {
    symmetry.validate()?;
    let n = symmetry.field_periods as usize;
    if full.is_empty() {
        return Err(Error::InvalidLattice("lattice needs at least one point"));
    }
    if full.len() % n != 0 {
        return Err(Error::InvalidSymmetry(
            "field-period fold needs a point count divisible by the period count",
        ));
    }
    let mut base = Vec::with_capacity(full.len() / n);
    for orbit in full.chunks(n) {
        let first = orbit[0];
        first.validate()?;
        for other in &orbit[1..] {
            other.validate()?;
            if other.ion_temperature_kev != first.ion_temperature_kev {
                return Err(Error::InvalidSymmetry(
                    "field-period fold needs matching ion temperatures across each orbit",
                ));
            }
        }
        let rate: f64 = orbit.iter().map(|p| p.rate).sum();
        base.push(LatticePoint {
            position_cm: first.position_cm,
            rate: rate / n as f64,
            ion_temperature_kev: first.ion_temperature_kev,
        });
    }
    Ok(base)
}

/// Deterministic sampler for a lattice source.
pub struct LatticeSampler {
    config: LatticeSourceConfig,
    /// Cumulative rate over the config points (base sector when symmetric).
    cumulative: Vec<f64>,
    /// Bare rate sum (one sector).
    total: f64,
    rng: Rng,
}

impl LatticeSampler {
    /// New sampler over `config` pinned to `seed`.
    pub fn new(config: LatticeSourceConfig, seed: u64) -> Result<Self> {
        config.validate()?;
        let mut cumulative = Vec::with_capacity(config.points.len() + 1);
        cumulative.push(0.0);
        let mut acc = 0.0;
        for point in &config.points {
            acc += point.rate;
            cumulative.push(acc);
        }
        if acc.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || !acc.is_finite() {
            return Err(Error::InvalidLattice(
                "total lattice birth rate is zero (check point rates)",
            ));
        }
        Ok(Self {
            config,
            cumulative,
            total: acc,
            rng: Rng::new(seed),
        })
    }

    /// Borrow the underlying configuration.
    pub fn config(&self) -> &LatticeSourceConfig {
        &self.config
    }

    /// Pick a point index and the within-bin fraction from one uniform draw.
    /// A single node without symmetry (or with the trivial period) skips the
    /// draw (the point-source recovery anchor: no index draw on the one-node
    /// path, so the stream matches [`crate::SourceSampler`] bit-for-bit for
    /// single-fuel configs). A single node with a non-trivial symmetry draws
    /// the field-period copy fraction directly — still exactly one selection
    /// draw per particle, so the folded stream stays in lockstep with the
    /// expanded full-cloud stream.
    fn pick_index(&mut self) -> (usize, f64) {
        let last = self.config.points.len() - 1;
        if last == 0 {
            let frac = match &self.config.symmetry {
                Some(symmetry) if symmetry.field_periods > 1 => self.rng.uniform(),
                _ => 0.0,
            };
            return (0, frac);
        }
        let u = self.rng.uniform();
        let target = u * self.total;
        let idx = match self.cumulative.binary_search_by(|v| v.total_cmp(&target)) {
            Ok(i) => i.min(last),
            Err(i) => i.saturating_sub(1).min(last),
        };
        let width = (self.cumulative[idx + 1] - self.cumulative[idx]).max(1e-300);
        let frac = ((target - self.cumulative[idx]) / width).min(1.0);
        (idx, frac)
    }

    /// Draw the birth spectrum moments `(mean, sigma)` \[MeV\] at species
    /// temperatures `(t_d, t_t)`: the config reaction's Ballabio moments, or
    /// — for a fuel mixture — the branch roulette weighted by the rate rule
    /// (module rustdoc). The single-fuel path consumes no RNG, preserving the
    /// landed sampling stream bit-for-bit on the recovery path.
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

    /// Sample one particle.
    ///
    /// The `Err(_) => 0.0` arm is unreachable-in-practice (the config is
    /// validated at construction, so `moments_mev` cannot fail here); the
    /// zero fallback keeps the infallible sampling stream total. A fallible
    /// caller should use [`LatticeSampler::try_sample`].
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

    /// Fallible single-particle sample: like [`LatticeSampler::sample`] but
    /// surfaces a spectrum-moment failure as a loud [`Error`].
    pub fn try_sample(&mut self) -> Result<Particle> {
        let (idx, frac) = self.pick_index();
        let point = self.config.points[idx];
        // With symmetry the config points are the base sector: the
        // field-period copy comes from the within-bin fraction (no extra RNG
        // draw), so the folded stream reproduces the expanded stream
        // bit-for-bit (module rustdoc gate).
        let position_cm = match &self.config.symmetry {
            None => point.position_cm,
            Some(symmetry) => {
                let n = symmetry.field_periods as usize;
                if n == 1 {
                    point.position_cm
                } else {
                    let copy = (frac * n as f64).floor() as usize;
                    let copy = copy.min(n - 1);
                    if copy == 0 {
                        point.position_cm
                    } else {
                        let (x, y) = rotate_xy(
                            point.position_cm[0],
                            point.position_cm[1],
                            copy as f64 * symmetry.period(),
                        );
                        [x, y, point.position_cm[2]]
                    }
                }
            }
        };
        let (t_d, t_t) = self
            .config
            .species_temperatures_at(point.ion_temperature_kev);
        let (mean, sigma) = self.sample_spectrum_moments(t_d, t_t)?;
        let energy_mev = if sigma > 0.0 {
            mean + sigma * self.rng.standard_normal()
        } else {
            mean
        };
        Ok(Particle {
            position_cm,
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

/// Marginal histograms of the lattice source for card emission:
/// cylindrical-`R` birth profile, vertical birth profile, and the global
/// rate-weighted marginal energy spectrum (`bins` controls the output binning
/// of each marginal).
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeEmissionHistograms {
    /// Cylindrical-`R` marginal over `[0, R_max]` \[cm\].
    pub radial: BinnedDistribution,
    /// Vertical marginal over the cloud `z` span \[cm\].
    pub vertical: BinnedDistribution,
    /// Global marginal energy spectrum \[MeV\].
    pub energy: BinnedDistribution,
    /// Captured energy probability mass (1 − tail truncation).
    pub energy_coverage: f64,
    /// Half the L1 distance between the true `(R, z)` birth joint and the
    /// product of the two marginals (0 = independent).
    pub joint_correlation: f64,
    /// Node count of the emitted cloud.
    pub node_count: usize,
    /// Total relative source strength (symmetry-aware).
    pub total_rate: f64,
}

/// Energy tabulation half-width in sigma (the ring/point convention).
const ENERGY_WIDTH_SIGMA: f64 = 4.0;

/// Build the marginal histograms for one lattice configuration.
pub fn lattice_emission_histograms(
    config: &LatticeSourceConfig,
    bins: usize,
) -> Result<LatticeEmissionHistograms> {
    use crate::spectrum::normal_cdf;

    config.validate()?;
    let bins = bins.max(4);
    let total = config.base_total_rate()?;
    let copies = match &config.symmetry {
        None => 1.0,
        Some(symmetry) => symmetry.copies(),
    };

    let r_max = config
        .points
        .iter()
        .map(|p| (p.position_cm[0] * p.position_cm[0] + p.position_cm[1] * p.position_cm[1]).sqrt())
        .fold(0.0f64, f64::max);
    let r_span = if r_max > 0.0 { r_max } else { 1.0 };
    let z_min = config
        .points
        .iter()
        .map(|p| p.position_cm[2])
        .fold(f64::INFINITY, f64::min);
    let z_max = config
        .points
        .iter()
        .map(|p| p.position_cm[2])
        .fold(f64::NEG_INFINITY, f64::max);
    let z_span = if z_max > z_min { z_max - z_min } else { 1.0 };
    let z_lo = if z_max > z_min { z_min } else { z_min - 0.5 };

    let mut energy_lo = f64::INFINITY;
    let mut energy_hi = f64::NEG_INFINITY;
    for point in &config.points {
        let (t_d, t_t) = config.species_temperatures_at(point.ion_temperature_kev);
        for &(_, _, mu, sigma) in &config.spectrum_branches(t_d, t_t)? {
            energy_lo = energy_lo.min(mu - ENERGY_WIDTH_SIGMA * sigma);
            energy_hi = energy_hi.max(mu + ENERGY_WIDTH_SIGMA * sigma);
        }
    }
    // A fully monoenergetic cloud (all Ti = 0) has a degenerate window —
    // widen symmetrically so the line lands in one bin at coverage 1.
    if !energy_lo.is_finite() || !energy_hi.is_finite() {
        return Err(Error::InvalidLattice("lattice energy window is empty"));
    }
    if energy_hi <= energy_lo {
        let mid = 0.5 * (energy_lo + energy_hi);
        let half = mid.abs() * 1e-6 + 1e-9;
        energy_lo = mid - half;
        energy_hi = mid + half;
    }
    // Spatial marginals plus the (R, z) joint for the correlation distance.
    let mut radial_m = vec![0.0f64; bins];
    let mut vertical_m = vec![0.0f64; bins];
    let mut joint = vec![vec![0.0f64; bins]; bins];
    for point in &config.points {
        let r = (point.position_cm[0] * point.position_cm[0]
            + point.position_cm[1] * point.position_cm[1])
            .sqrt();
        let ri = ((r / r_span) * bins as f64).floor() as usize;
        let ri = ri.min(bins - 1);
        let zi = (((point.position_cm[2] - z_lo) / z_span) * bins as f64).floor() as usize;
        let zi = zi.min(bins - 1);
        radial_m[ri] += point.rate;
        vertical_m[zi] += point.rate;
        joint[ri][zi] += point.rate / total;
    }
    let radial: Vec<f64> = radial_m.iter().map(|m| m / total).collect();
    let vertical: Vec<f64> = vertical_m.iter().map(|m| m / total).collect();
    let mut corr = 0.0;
    for i in 0..bins {
        for j in 0..bins {
            corr += (joint[i][j] - radial[i] * vertical[j]).abs();
        }
    }
    corr *= 0.5;

    // Rate-weighted marginal energy spectrum over the global window.
    let de = (energy_hi - energy_lo) / bins as f64;
    let mut energy_m = vec![0.0f64; bins];
    let mut energy_coverage = 0.0f64;
    for point in &config.points {
        let (t_d, t_t) = config.species_temperatures_at(point.ion_temperature_kev);
        let branches = config.spectrum_branches(t_d, t_t)?;
        for (b, mass) in energy_m.iter_mut().enumerate() {
            let lo = energy_lo + b as f64 * de;
            let hi = lo + de;
            let mut p = 0.0;
            for &(_, weight, mu, sigma) in &branches {
                let pb = if sigma > 0.0 {
                    normal_cdf((hi - mu) / sigma) - normal_cdf((lo - mu) / sigma)
                } else if lo < mu && mu <= hi {
                    1.0
                } else {
                    0.0
                };
                p += weight * pb;
            }
            *mass += point.rate * p;
        }
        let mut covered = 0.0;
        for &(_, weight, mu, sigma) in &branches {
            let cb = if sigma > 0.0 {
                normal_cdf((energy_hi - mu) / sigma) - normal_cdf((energy_lo - mu) / sigma)
            } else {
                1.0
            };
            covered += weight * cb;
        }
        energy_coverage += point.rate * covered;
    }
    energy_coverage /= total;

    Ok(LatticeEmissionHistograms {
        radial: BinnedDistribution {
            centers: (0..bins)
                .map(|b| (b as f64 + 0.5) * (r_span / bins as f64))
                .collect(),
            masses: radial,
        },
        vertical: BinnedDistribution {
            centers: (0..bins)
                .map(|b| z_lo + (b as f64 + 0.5) * (z_span / bins as f64))
                .collect(),
            masses: vertical,
        },
        energy: BinnedDistribution {
            centers: (0..bins)
                .map(|b| energy_lo + (b as f64 + 0.5) * de)
                .collect(),
            masses: energy_m.iter().map(|m| m / total).collect(),
        },
        energy_coverage,
        joint_correlation: corr,
        node_count: config.points.len(),
        total_rate: copies * total,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn two_point_lattice() -> LatticeSourceConfig {
        LatticeSourceConfig {
            points: vec![
                LatticePoint {
                    position_cm: [300.0, 0.0, 25.0],
                    rate: 2.0,
                    ion_temperature_kev: 0.0,
                },
                LatticePoint {
                    position_cm: [-300.0, 0.0, 25.0],
                    rate: 1.0,
                    ion_temperature_kev: 0.0,
                },
            ],
            reaction: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        }
    }

    fn ring_lattice(k: usize, radius: f64, height: f64) -> LatticeSourceConfig {
        let points = (0..k)
            .map(|i| {
                let phi = 2.0 * PI * i as f64 / k as f64;
                LatticePoint {
                    position_cm: [radius * phi.cos(), radius * phi.sin(), height],
                    rate: 1.0,
                    ion_temperature_kev: 0.0,
                }
            })
            .collect();
        LatticeSourceConfig {
            points,
            reaction: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        }
    }

    #[test]
    fn hand_vectors_total_mean_position_and_energy() {
        let config = two_point_lattice();
        assert_eq!(config.total_strength().unwrap(), 3.0);
        assert_eq!(config.mean_birth_position().unwrap(), [100.0, 0.0, 25.0]);
        assert_eq!(config.mean_birth_energy().unwrap(), 14.021);
    }

    #[test]
    fn single_point_reproduces_point_source_stream_bit_for_bit() {
        use crate::{PlasmaSourceConfig, SourceSampler};
        let lattice = LatticeSourceConfig {
            points: vec![LatticePoint {
                position_cm: [1.0, -2.0, 3.5],
                rate: 1.0,
                ion_temperature_kev: 20.0,
            }],
            reaction: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            symmetry: None,
            weight: 2.5,
        };
        let point =
            PlasmaSourceConfig::point(1.0, -2.0, 3.5, FusionReaction::Dt, 20.0).with_weight(2.5);
        for seed in [0, 1, 42, 2026] {
            let a = LatticeSampler::new(lattice.clone(), seed)
                .unwrap()
                .sample_n(64);
            let b = SourceSampler::new(point, seed).unwrap().sample_n(64);
            assert_eq!(a, b, "seed {seed}");
        }
    }

    #[test]
    fn expand_is_orbit_major_with_exact_base_copy() {
        let base = two_point_lattice();
        let symmetry = LatticeSymmetry::new(4, 0.0).unwrap();
        let full = expand_lattice(&base.points, &symmetry);
        assert_eq!(full.len(), 8);
        // Copy 0 of each orbit is the untouched base point, bit-for-bit.
        assert_eq!(full[0], base.points[0]);
        assert_eq!(full[4], base.points[1]);
        // Copy k rotates by k·period; quarter-turn of (300, 0) is (0, 300).
        assert!((full[1].position_cm[0]).abs() < 1e-9);
        assert!((full[1].position_cm[1] - 300.0).abs() < 1e-9);
        assert_eq!(full[1].rate, 2.0);
        assert_eq!(full[1].ion_temperature_kev, 0.0);
    }

    #[test]
    fn fold_inverts_expand_with_exact_summed_rates() {
        let base = two_point_lattice();
        let symmetry = LatticeSymmetry::new(4, 0.0).unwrap();
        let full = expand_lattice(&base.points, &symmetry);
        let folded = fold_lattice(&full, &symmetry).unwrap();
        assert_eq!(folded.len(), 2);
        assert_eq!(folded[0].position_cm, base.points[0].position_cm);
        assert_eq!(folded[0].rate, 2.0);
        assert_eq!(folded[1].rate, 1.0);
        // Fold inverts expand on exact-friendly rates, bit-for-bit.
        assert_eq!(folded, base.points);
        // Trivial period is the identity on data.
        let one = LatticeSymmetry::new(1, 0.5).unwrap();
        assert_eq!(expand_lattice(&base.points, &one), base.points);
        assert_eq!(fold_lattice(&base.points, &one).unwrap(), base.points);
    }

    #[test]
    fn folded_sampler_reproduces_expanded_stream_bit_for_bit() {
        let base = two_point_lattice();
        let symmetry = LatticeSymmetry::new(4, 0.0).unwrap();
        let full_points = expand_lattice(&base.points, &symmetry);
        let folded_rates = fold_lattice(&full_points, &symmetry).unwrap();
        let folded = LatticeSourceConfig {
            points: folded_rates,
            symmetry: Some(symmetry),
            ..base.clone()
        };
        // Symmetry-aware total equals the expanded bare sum, exactly.
        let expanded = LatticeSourceConfig {
            points: full_points,
            symmetry: None,
            ..base.clone()
        };
        assert_eq!(
            folded.total_strength().unwrap(),
            expanded.total_strength().unwrap()
        );
        for seed in [0, 1, 7, 42, 2026] {
            let a = LatticeSampler::new(folded.clone(), seed)
                .unwrap()
                .sample_n(256);
            let b = LatticeSampler::new(expanded.clone(), seed)
                .unwrap()
                .sample_n(256);
            assert_eq!(a, b, "seed {seed}");
        }
    }

    #[test]
    fn single_node_symmetry_replicates_and_matches_expanded_stream() {
        let base = LatticeSourceConfig {
            points: vec![LatticePoint {
                position_cm: [300.0, 0.0, 0.0],
                rate: 1.0,
                ion_temperature_kev: 0.0,
            }],
            reaction: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        };
        let symmetry = LatticeSymmetry::new(4, 0.0).unwrap();
        let folded = LatticeSourceConfig {
            symmetry: Some(symmetry),
            ..base.clone()
        };
        let expanded = LatticeSourceConfig {
            points: expand_lattice(&base.points, &symmetry),
            ..base.clone()
        };
        // The lone base node fans out over all four copies, bit-identical to
        // sampling the expanded full cloud (same seed).
        for seed in [0, 5, 42] {
            let a = LatticeSampler::new(folded.clone(), seed)
                .unwrap()
                .sample_n(256);
            let b = LatticeSampler::new(expanded.clone(), seed)
                .unwrap()
                .sample_n(256);
            assert_eq!(a, b, "seed {seed}");
        }
        let cloud = LatticeSampler::new(folded, 5).unwrap().sample_n(4000);
        let mut quadrants = [0usize; 4];
        for p in &cloud {
            let qx = usize::from(p.position_cm[0] > 0.0);
            let qy = usize::from(p.position_cm[1] > 0.0);
            quadrants[qx + 2 * qy] += 1;
        }
        assert!(quadrants.iter().all(|&q| q > 0), "quadrants {quadrants:?}");
    }

    #[test]
    fn ring_lattice_moments_match_closed_form() {
        let config = ring_lattice(64, 300.0, 25.0);
        assert_eq!(config.total_strength().unwrap(), 64.0);
        let mut sampler = LatticeSampler::new(config, 42).unwrap();
        let particles = sampler.sample_n(200_000);
        let n = particles.len() as f64;
        let mut mean = [0.0f64; 3];
        for p in &particles {
            for (m, v) in mean.iter_mut().zip(p.position_cm.iter()) {
                *m += v;
            }
            let r =
                (p.position_cm[0] * p.position_cm[0] + p.position_cm[1] * p.position_cm[1]).sqrt();
            assert!((r - 300.0).abs() < 1e-9, "radius {r}");
            assert!((p.position_cm[2] - 25.0).abs() < 1e-12);
            assert_eq!(p.energy_mev, 14.021);
        }
        for m in mean.iter_mut() {
            *m /= n;
        }
        let se = 300.0 / (2.0 * n).sqrt();
        assert!(mean[0].abs() < 8.0 * se, "mean x {}", mean[0]);
        assert!(mean[1].abs() < 8.0 * se, "mean y {}", mean[1]);
        assert!((mean[2] - 25.0).abs() < 1e-12);
    }

    #[test]
    fn equal_species_pair_recovers_shared_kernel_bit_for_bit() {
        let flat = LatticeSourceConfig {
            points: vec![
                LatticePoint {
                    position_cm: [10.0, 0.0, 0.0],
                    rate: 1.0,
                    ion_temperature_kev: 20.0,
                },
                LatticePoint {
                    position_cm: [0.0, 10.0, 5.0],
                    rate: 3.0,
                    ion_temperature_kev: 20.0,
                },
            ],
            reaction: FusionReaction::Dt,
            fuel_mixture: Some(FuelMixture::new(0.5, 0.5).unwrap()),
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        };
        let pair = LatticeSourceConfig {
            species_temperatures: Some(SpeciesIonTemperatures::new(20.0, 20.0).unwrap()),
            ..flat.clone()
        };
        for seed in [0, 17] {
            let a = LatticeSampler::new(flat.clone(), seed)
                .unwrap()
                .sample_n(128);
            let b = LatticeSampler::new(pair.clone(), seed)
                .unwrap()
                .sample_n(128);
            assert_eq!(a, b, "seed {seed}");
        }
        assert_eq!(
            flat.mean_birth_energy().unwrap(),
            pair.mean_birth_energy().unwrap()
        );
    }

    #[test]
    fn mixture_lattice_fires_both_branches() {
        let config = LatticeSourceConfig {
            points: vec![LatticePoint {
                position_cm: [0.0, 0.0, 0.0],
                rate: 1.0,
                ion_temperature_kev: 20.0,
            }],
            reaction: FusionReaction::Dt,
            fuel_mixture: Some(FuelMixture::new(0.7, 0.3).unwrap()),
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        };
        let particles = LatticeSampler::new(config, 11).unwrap().sample_n(200_000);
        let frac_dd = particles.iter().filter(|p| p.energy_mev < 10.0).count() as f64 / 200_000.0;
        assert!((0.002..0.015).contains(&frac_dd), "D-D share {frac_dd}");
    }

    #[test]
    fn invalid_lattices_are_loud() {
        let empty = LatticeSourceConfig {
            points: vec![],
            reaction: FusionReaction::Dt,
            fuel_mixture: None,
            species_temperatures: None,
            symmetry: None,
            weight: 1.0,
        };
        assert_eq!(
            empty.validate(),
            Err(Error::InvalidLattice("lattice needs at least one point"))
        );
        let mut negative = two_point_lattice();
        negative.points[0].rate = -1.0;
        assert_eq!(
            negative.validate(),
            Err(Error::InvalidLattice("lattice point rate must be >= 0"))
        );
        let mut nan = two_point_lattice();
        nan.points[0].position_cm[0] = f64::NAN;
        assert_eq!(nan.validate(), Err(Error::NonFinite("lattice point x")));
        let mut zero = two_point_lattice();
        zero.points[0].rate = 0.0;
        zero.points[1].rate = 0.0;
        assert!(zero.total_strength().is_err());
        assert_eq!(
            LatticeSymmetry::new(0, 0.0),
            Err(Error::InvalidSymmetry("field periods must be >= 1"))
        );
        assert_eq!(
            LatticeSymmetry::new(2, f64::INFINITY),
            Err(Error::NonFinite("lattice symmetry base angle"))
        );
        // Fold needs a multiple of n and matching Ti.
        let symmetry = LatticeSymmetry::new(4, 0.0).unwrap();
        assert_eq!(
            fold_lattice(&two_point_lattice().points, &symmetry),
            Err(Error::InvalidSymmetry(
                "field-period fold needs a point count divisible by the period count"
            ))
        );
        let mut mismatched = expand_lattice(&two_point_lattice().points, &symmetry);
        mismatched[1].ion_temperature_kev = 5.0;
        assert_eq!(
            fold_lattice(&mismatched, &symmetry),
            Err(Error::InvalidSymmetry(
                "field-period fold needs matching ion temperatures across each orbit"
            ))
        );
        // Species pair without a mixture stays loud.
        let mut no_mixture = two_point_lattice();
        no_mixture.species_temperatures = Some(SpeciesIonTemperatures::new(20.0, 20.0).unwrap());
        assert!(no_mixture.validate().is_err());
    }

    #[test]
    fn emission_histograms_conserve_and_bound_correlation() {
        let hist = lattice_emission_histograms(&two_point_lattice(), 8).unwrap();
        assert_eq!(hist.node_count, 2);
        assert_eq!(hist.total_rate, 3.0);
        for dist in [&hist.radial, &hist.vertical, &hist.energy] {
            let sum: f64 = dist.masses.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "masses sum {sum}");
        }
        assert!((hist.energy_coverage - 1.0).abs() < 1e-12);
        assert!((0.0..=1.0).contains(&hist.joint_correlation));
    }
}
