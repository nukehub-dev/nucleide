//! Seeded sampling of source configurations to particle vectors.
//!
//! [`SourceSampler`] turns a [`PlasmaSourceConfig`] into a stream of
//! [`Particle`] records (species, position, direction, energy, weight) over a small
//! self-contained RNG (SplitMix64 seeding + xoshiro256**), so a pinned seed
//! reproduces the same stream on a given platform. There is deliberately no
//! MCPL (or any file-format) dependency here: MCPL projection stays
//! caller-side, the same layering rule as `vr-tools` `KdeSampler` — the
//! caller writes particle vectors with `nucleide-mcpl-io` / `nucleide-mcnp-io`
//! if it wants a file.
//!
//! Determinism note: the stream is fixed by the seed and the crate version.
//! `sin`/`cos`/`ln`/`sqrt` come from the platform libm, so the last-ulp
//! stream can differ across platforms; moments and card text are unaffected
//! at any meaningful tolerance.
//!
//! The D(d,p)T proton sampler ([`SourceSampler::proton`]) is an independent
//! instance with its own `Rng(seed)`: never a joint neutron+proton stream,
//! never a per-event product roulette on the neutron path. Its draw
//! accounting mirrors the neutron path (ring azimuth 1 uniform, energy
//! Box–Muller ≤ 2 uniforms, direction 2 uniforms), so sampling neutrons
//! then protons from two samplers with the same seed leaves the neutron
//! stream bit-for-bit.

use std::f64::consts::PI;

use nucleide_nuclei::particles::ParticleId;

use crate::{Error, FusionReaction, PlasmaSourceConfig, Result, SourceModel, SpectrumSpec};

/// One sampled source particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
    /// Source particle species: [`ParticleId::Neutron`] on neutron streams,
    /// [`ParticleId::Proton`] on proton streams. Proton streams are
    /// independent sampler instances (own `Rng(seed)`), never a joint
    /// neutron+proton stream and never a per-event roulette on the neutron
    /// path.
    pub species: ParticleId,
    /// Birth position `(x, y, z)` \[cm\].
    pub position_cm: [f64; 3],
    /// Unit direction `(u, v, w)`.
    pub direction: [f64; 3],
    /// Kinetic energy \[MeV\].
    pub energy_mev: f64,
    /// Particle weight.
    pub weight: f64,
}

/// SplitMix64 seed expander and uniform stream (public-domain algorithm,
/// implemented clean-room from the published reference).
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// xoshiro256** — the sampler's core uniform generator (public-domain
/// algorithm, implemented clean-room from the published reference).
struct Xoshiro256 {
    s: [u64; 4],
}

impl Xoshiro256 {
    fn from_seed(seed: u64) -> Self {
        let mut sm = SplitMix64::new(seed);
        Self {
            s: [sm.next_u64(), sm.next_u64(), sm.next_u64(), sm.next_u64()],
        }
    }

    fn next_u64(&mut self) -> u64 {
        let [s0, s1, s2, s3] = self.s;
        let result = s1.wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s1 << 17;
        self.s[2] = s2 ^ s0;
        self.s[3] = s3 ^ s1;
        self.s[1] = s1 ^ self.s[2];
        self.s[0] = s0 ^ self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform in `[0, 1)` (53-bit mantissa granularity).
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Uniform/deviate facade shared by the ring/point and parametric samplers.
pub(crate) struct Rng {
    inner: Xoshiro256,
    /// Cached second Box–Muller deviate.
    normal_spare: Option<f64>,
}

impl Rng {
    pub(crate) fn new(seed: u64) -> Self {
        Self {
            inner: Xoshiro256::from_seed(seed),
            normal_spare: None,
        }
    }

    /// Uniform in `[0, 1)` (53-bit mantissa granularity).
    pub(crate) fn uniform(&mut self) -> f64 {
        self.inner.uniform()
    }

    /// Standard normal deviate (Box–Muller, cached pair).
    pub(crate) fn standard_normal(&mut self) -> f64 {
        if let Some(spare) = self.normal_spare.take() {
            return spare;
        }
        let u1 = 1.0 - self.uniform();
        let u2 = self.uniform();
        let radius = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * PI * u2;
        self.normal_spare = Some(radius * theta.sin());
        radius * theta.cos()
    }

    /// Isotropic unit direction: cos(theta) uniform in [-1, 1], azimuth
    /// uniform in [0, 2*pi).
    pub(crate) fn isotropic_direction(&mut self) -> [f64; 3] {
        let mu = 2.0 * self.uniform() - 1.0;
        let phi = 2.0 * PI * self.uniform();
        let sin_theta = (1.0 - mu * mu).sqrt();
        [sin_theta * phi.cos(), sin_theta * phi.sin(), mu]
    }
}

/// Deterministic sampler for one source configuration.
pub struct SourceSampler {
    config: PlasmaSourceConfig,
    spectrum: SpectrumSpec,
    species: ParticleId,
    rng: Rng,
}

impl SourceSampler {
    /// New neutron sampler over `config` pinned to `seed`.
    ///
    /// Errors: the configuration is validated up front (see
    /// [`PlasmaSourceConfig::validate`]).
    pub fn new(config: PlasmaSourceConfig, seed: u64) -> Result<Self> {
        Self::build(config, seed, ParticleId::Neutron)
    }

    /// New D(d,p)T proton sampler over `config` pinned to `seed`.
    ///
    /// The proton spectrum is the Brysk-Gaussian 3.02 MeV line
    /// ([`SpectrumSpec::from_proton`]); the warm draw accounting mirrors the
    /// neutron path (ring azimuth 1 uniform, energy Box–Muller ≤ 2
    /// uniforms, direction 2 uniforms) on a fresh `Rng(seed)` instance, so
    /// a neutron sampler with the same seed reproduces its stream
    /// bit-for-bit (module rustdoc: separate streams, no neutron-path
    /// draws).
    ///
    /// Errors: the configuration is validated up front, and a D-T-only
    /// config is loud: [`Error::NotYetSupported`]
    /// (`"proton sampling needs a D-D source"`).
    pub fn proton(config: PlasmaSourceConfig, seed: u64) -> Result<Self> {
        Self::build(config, seed, ParticleId::Proton)
    }

    /// Shared construction: neutron (the landed path) and proton (the
    /// D(d,p)T branch) select only the species tag and the spectrum; the
    /// per-particle draw accounting is identical.
    fn build(config: PlasmaSourceConfig, seed: u64, species: ParticleId) -> Result<Self> {
        config.validate()?;
        let spectrum = match species {
            ParticleId::Proton => {
                if config.reaction != FusionReaction::Dd {
                    return Err(Error::NotYetSupported("proton sampling needs a D-D source"));
                }
                SpectrumSpec::from_proton(config.ion_temperature_kev)?
            }
            _ => config.spectrum()?,
        };
        Ok(Self {
            config,
            spectrum,
            species,
            rng: Rng::new(seed),
        })
    }

    /// Borrow the underlying configuration.
    pub fn config(&self) -> &PlasmaSourceConfig {
        &self.config
    }

    /// Standard normal deviate (Box–Muller, cached pair).
    fn standard_normal(&mut self) -> f64 {
        self.rng.standard_normal()
    }

    /// Isotropic unit direction: cos(theta) uniform in [-1, 1], azimuth
    /// uniform in [0, 2*pi).
    fn isotropic_direction(&mut self) -> [f64; 3] {
        self.rng.isotropic_direction()
    }

    /// Sample one particle.
    pub fn sample(&mut self) -> Particle {
        let position_cm = match self.config.model {
            SourceModel::Point(p) => [p.x_cm, p.y_cm, p.z_cm],
            SourceModel::Ring(r) => {
                let phi = 2.0 * PI * self.rng.uniform();
                [
                    r.radius_cm * phi.cos(),
                    r.radius_cm * phi.sin(),
                    r.height_cm,
                ]
            }
        };
        let energy_mev = match self.spectrum {
            SpectrumSpec::Mono { energy_mev } => energy_mev,
            SpectrumSpec::Gaussian {
                mean_mev,
                sigma_mev,
            } => mean_mev + sigma_mev * self.standard_normal(),
        };
        Particle {
            species: self.species,
            position_cm,
            direction: self.isotropic_direction(),
            energy_mev,
            weight: self.config.weight,
        }
    }

    /// Sample `n` particles (simple forward stream; `n` particles consume
    /// `n` samples in order).
    pub fn sample_n(&mut self, n: usize) -> Vec<Particle> {
        (0..n).map(|_| self.sample()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FusionReaction;

    const N: usize = 200_000;

    fn ring_config() -> PlasmaSourceConfig {
        PlasmaSourceConfig::ring(300.0, 25.0, FusionReaction::Dt, 20.0)
    }

    /// Population moments of a ring with uniform azimuth: every spatial
    /// moment is closed form (radius fixed, height fixed, phi uniform).
    #[test]
    fn ring_spatial_moments_match_closed_form() {
        let mut sampler = SourceSampler::new(ring_config(), 42).unwrap();
        let particles = sampler.sample_n(N);
        let mut mean = [0.0f64; 3];
        let mut radius_sum = 0.0;
        let mut cos_phi_sum = 0.0;
        for p in &particles {
            for (m, v) in mean.iter_mut().zip(p.position_cm.iter()) {
                *m += v;
            }
            let r =
                (p.position_cm[0] * p.position_cm[0] + p.position_cm[1] * p.position_cm[1]).sqrt();
            radius_sum += r;
            cos_phi_sum += p.position_cm[0] / r;
            assert!((p.position_cm[2] - 25.0).abs() < 1e-12);
        }
        for m in mean.iter_mut() {
            *m /= N as f64;
        }
        // <x> = <y> = 0 (uniform azimuth); standard error of the mean is
        // R/sqrt(2N) ~ 0.47 cm, so gate at 8 standard errors.
        let se = 300.0 / (2.0_f64 * N as f64).sqrt();
        assert!(mean[0].abs() < 8.0 * se, "mean x {}", mean[0]);
        assert!(mean[1].abs() < 8.0 * se, "mean y {}", mean[1]);
        assert!((mean[2] - 25.0).abs() < 1e-12);
        // <r> = R exactly for every particle.
        assert!((radius_sum / N as f64 - 300.0).abs() < 1e-9);
        // <cos(phi)> = 0 with the same sampling error.
        assert!(cos_phi_sum.abs() / (N as f64) < 8.0 * se / 300.0);
    }

    #[test]
    fn spectrum_moments_match_closed_form() {
        // DT at 20 keV: mean 14.0728742533532 MeV, sigma 0.337721762921614 MeV.
        let mut sampler = SourceSampler::new(ring_config(), 7).unwrap();
        let particles = sampler.sample_n(N);
        let (mean, sigma) = FusionReaction::Dt.moments_mev(20.0).unwrap();
        let mut e_mean = 0.0;
        let mut e_sq = 0.0;
        for p in &particles {
            e_mean += p.energy_mev;
            e_sq += p.energy_mev * p.energy_mev;
        }
        e_mean /= N as f64;
        let e_var = e_sq / N as f64 - e_mean * e_mean;
        // Standard error of the mean ~ sigma/sqrt(N) ~ 7.6e-4 MeV.
        assert!(
            (e_mean - mean).abs() < 5e-3,
            "sampled mean {e_mean} vs {mean}"
        );
        assert!(
            (e_var.sqrt() - sigma).abs() < 5e-3,
            "sampled sigma {} vs {sigma}",
            e_var.sqrt()
        );
    }

    #[test]
    fn directions_are_isotropic() {
        let mut sampler = SourceSampler::new(ring_config(), 3).unwrap();
        let particles = sampler.sample_n(N);
        let mut mean_w = 0.0;
        let mut mean_w2 = 0.0;
        let mut mean_w3 = 0.0;
        for p in &particles {
            let w = p.direction[2];
            mean_w += w;
            mean_w2 += w * w;
            mean_w3 += w * w * w;
            let norm = p.direction[0] * p.direction[0]
                + p.direction[1] * p.direction[1]
                + p.direction[2] * p.direction[2];
            assert!((norm - 1.0).abs() < 1e-12, "unit direction, got {norm}");
        }
        assert!(
            mean_w.abs() / (N as f64) < 0.01,
            "<w> = {}",
            mean_w / N as f64
        );
        // Isotropic: <w^2> = 1/3.
        assert!((mean_w2 / N as f64 - 1.0 / 3.0).abs() < 0.01);
        // Symmetry: <w^3> ~ 0.
        assert!(mean_w3.abs() / (N as f64) < 0.01);
    }

    #[test]
    fn pinned_seed_reproduces_the_stream() {
        let a = SourceSampler::new(ring_config(), 1234)
            .unwrap()
            .sample_n(64);
        let b = SourceSampler::new(ring_config(), 1234)
            .unwrap()
            .sample_n(64);
        assert_eq!(a, b);
        let c = SourceSampler::new(ring_config(), 1235)
            .unwrap()
            .sample_n(64);
        assert_ne!(a, c);
    }

    #[test]
    fn golden_stream_pins_the_rng() {
        // First four particles of a pinned configuration, recorded to pin the
        // RNG stream across refactors (regenerate by printing on failure).
        let particles = SourceSampler::new(ring_config(), 2026).unwrap().sample_n(4);
        let golden: [[f64; 3]; 4] = [
            [-268.378_974_680_236_4, -134.062_395_732_677_5, 25.0],
            [72.383_198_813_085_1, -291.136_862_196_434, 25.0],
            [192.694_047_107_026_68, -229.932_607_973_542_96, 25.0],
            [-33.306_004_347_049_054, -298.145_451_205_337_6, 25.0],
        ];
        let golden_energy = [
            14.178_460_351_981_025,
            13.818_009_560_047_981,
            14.316_129_698_208_458,
            14.020_042_918_011_535,
        ];
        for (p, (want, want_e)) in particles
            .iter()
            .zip(golden.iter().zip(golden_energy.iter()))
        {
            for k in 0..3 {
                assert!(
                    (p.position_cm[k] - want[k]).abs() < 1e-9,
                    "position {want:?} vs {:?}",
                    p.position_cm
                );
            }
            assert!(
                (p.energy_mev - want_e).abs() < 1e-12,
                "energy {want_e} vs {}",
                p.energy_mev
            );
        }
    }

    #[test]
    fn point_and_weight_are_carried() {
        let config =
            PlasmaSourceConfig::point(1.0, -2.0, 3.5, FusionReaction::Dd, 10.0).with_weight(2.5);
        let mut sampler = SourceSampler::new(config, 9).unwrap();
        let p = sampler.sample();
        assert_eq!(p.position_cm, [1.0, -2.0, 3.5]);
        assert_eq!(p.weight, 2.5);
        // DD at 10 keV: mean 2.48245474652507 MeV.
        assert!((p.energy_mev - 2.482_454_746_525_07).abs() < 0.2);
    }

    #[test]
    fn invalid_config_fails_at_construction() {
        let config = PlasmaSourceConfig::ring(0.0, 0.0, FusionReaction::Dt, 10.0);
        assert!(SourceSampler::new(config, 0).is_err());
    }

    #[test]
    fn neutron_stream_tags_neutron_species() {
        let particles = SourceSampler::new(ring_config(), 2026)
            .unwrap()
            .sample_n(16);
        assert!(particles
            .iter()
            .all(|p| p.species == nucleide_nuclei::particles::ParticleId::Neutron));
        assert!(particles.iter().all(|p| (p.weight - 1.0).abs() < 1e-12));
    }

    #[test]
    fn proton_moments_match_brysk_line_at_8_sigma() {
        // D(d,p)T proton line at Ti = 20 keV: mean 3.02 MeV, sigma 0.173782
        // MeV (reaction.rs hand vectors). 8 standard errors at N = 200k.
        let config = PlasmaSourceConfig::ring(300.0, 25.0, FusionReaction::Dd, 20.0);
        let mut sampler = SourceSampler::proton(config, 7).unwrap();
        let particles = sampler.sample_n(N);
        assert!(particles
            .iter()
            .all(|p| p.species == nucleide_nuclei::particles::ParticleId::Proton));
        let (mean, sigma) = (3.02_f64, 0.173782_f64);
        let mut e_mean = 0.0;
        let mut e_sq = 0.0;
        for p in &particles {
            e_mean += p.energy_mev;
            e_sq += p.energy_mev * p.energy_mev;
        }
        e_mean /= N as f64;
        let e_std = (e_sq / N as f64 - e_mean * e_mean).sqrt();
        let se_mean = sigma / (N as f64).sqrt();
        let se_std = sigma / (2.0 * N as f64).sqrt();
        assert!(
            (e_mean - mean).abs() < 8.0 * se_mean,
            "sampled mean {e_mean} vs {mean}"
        );
        assert!(
            (e_std - sigma).abs() < 8.0 * se_std,
            "sampled sigma {e_std} vs {sigma}"
        );
        // Weight is config.weight, exactly like neutrons.
        assert!(particles.iter().all(|p| p.weight == 1.0));
    }

    #[test]
    fn proton_zero_temperature_is_monoenergetic() {
        let config = PlasmaSourceConfig::point(1.0, 2.0, 3.0, FusionReaction::Dd, 0.0);
        let particles = SourceSampler::proton(config, 9).unwrap().sample_n(32);
        assert!(particles.iter().all(|p| p.energy_mev == 3.02));
        assert!(particles
            .iter()
            .all(|p| p.species == nucleide_nuclei::particles::ParticleId::Proton));
    }

    #[test]
    fn proton_dt_only_is_loud() {
        let config = PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, 20.0);
        assert_eq!(
            SourceSampler::proton(config, 1).map(|_| ()),
            Err(crate::Error::NotYetSupported(
                "proton sampling needs a D-D source"
            ))
        );
    }

    #[test]
    fn neutron_stream_is_independent_of_proton_sampler_same_seed() {
        // Separate instances, own Rng(seed): building/sampling a proton
        // sampler with the same seed cannot move the neutron stream.
        let reference = SourceSampler::new(ring_config(), 1234)
            .unwrap()
            .sample_n(64);
        let proton_config = PlasmaSourceConfig::ring(300.0, 25.0, FusionReaction::Dd, 20.0);
        let _protons = SourceSampler::proton(proton_config, 1234)
            .unwrap()
            .sample_n(64);
        let again = SourceSampler::new(ring_config(), 1234)
            .unwrap()
            .sample_n(64);
        assert_eq!(reference, again);
    }
}
