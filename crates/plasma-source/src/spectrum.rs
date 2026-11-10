//! Source energy spectrum: a monoenergetic line or a Brysk Gaussian.
//!
//! Cards cannot carry a continuous Gaussian through the discrete `SI`/`SP`
//! subset this workspace reads, so [`SpectrumSpec::tabulate`] renders the
//! Gaussian as discrete lines at bin centers with probabilities from the
//! exact normal CDF — the same discrete-distribution card shape the
//! spectroscopy decay-source emitter uses (E9). The tabulation truncates the
//! tails beyond a symmetric `±width_sigma` window; the dropped probability is
//! reported as [`SpectrumTable::coverage`] drift in the emission report.

use crate::reaction::proton_moments_mev;
use crate::{Error, FusionReaction, Result};

/// Standard normal CDF `Phi(z) = 0.5 * (1 + erf(z / sqrt(2)))`.
pub(crate) fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

/// Error function via the Abramowitz & Stegun 7.1.26 rational approximation
/// (|epsilon| <= 1.5e-7), implemented clean-room from the published formula.
/// Card probabilities are rendered to 6 significant digits, far above the
/// approximation error; unit tests cross-check against high-accuracy
/// numerical quadrature.
pub(crate) fn erf(x: f64) -> f64 {
    const P: f64 = 0.327_591_1;
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + P * ax);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-ax * ax).exp();
    sign * y
}

/// Energy spectrum of a fusion source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectrumSpec {
    /// Monoenergetic line \[MeV\].
    Mono {
        /// Line energy \[MeV\].
        energy_mev: f64,
    },
    /// Gaussian spectrum (Brysk 1973; Ballabio et al. 1998 coefficients).
    Gaussian {
        /// Mean \[MeV\].
        mean_mev: f64,
        /// Standard deviation \[MeV\]; strictly positive.
        sigma_mev: f64,
    },
}

impl SpectrumSpec {
    /// Spectrum of `reaction` at `ti_kev` \[keV\]. Zero ion temperature
    /// degenerates to the monoenergetic nominal line.
    pub fn from_reaction(reaction: FusionReaction, ti_kev: f64) -> Result<Self> {
        let (mean, sigma) = reaction.moments_mev(ti_kev)?;
        if sigma == 0.0 {
            Ok(SpectrumSpec::Mono { energy_mev: mean })
        } else {
            Ok(SpectrumSpec::Gaussian {
                mean_mev: mean,
                sigma_mev: sigma,
            })
        }
    }

    /// D(d,p)T proton spectrum at `ti_kev` \[keV\] (the Brysk-Gaussian
    /// 3.02 MeV line, [`proton_moments_mev`](crate::reaction::proton_moments_mev)).
    /// Zero ion temperature degenerates to the monoenergetic 3.02 MeV line.
    pub fn from_proton(ti_kev: f64) -> Result<Self> {
        let (mean, sigma) = proton_moments_mev(ti_kev)?;
        if sigma == 0.0 {
            Ok(SpectrumSpec::Mono { energy_mev: mean })
        } else {
            Ok(SpectrumSpec::Gaussian {
                mean_mev: mean,
                sigma_mev: sigma,
            })
        }
    }

    /// Central energy: the line energy, or the Gaussian mean.
    pub fn mean_mev(&self) -> f64 {
        match *self {
            SpectrumSpec::Mono { energy_mev } => energy_mev,
            SpectrumSpec::Gaussian { mean_mev, .. } => mean_mev,
        }
    }

    /// True for the monoenergetic form.
    pub fn is_mono(&self) -> bool {
        matches!(self, SpectrumSpec::Mono { .. })
    }

    /// Render the spectrum as discrete `(energy, probability)` card lines.
    ///
    /// A monoenergetic spectrum returns a single line with probability 1.
    /// A Gaussian returns `n_bins` lines at bin centers over the symmetric
    /// window `mean ± width_sigma * sigma`, with probabilities from the exact
    /// normal CDF. `n_bins >= 2` and `width_sigma > 0` are required for the
    /// Gaussian form.
    pub fn tabulate(&self, n_bins: usize, width_sigma: f64) -> Result<SpectrumTable> {
        match *self {
            SpectrumSpec::Mono { energy_mev } => Ok(SpectrumTable {
                energies_mev: vec![energy_mev],
                probabilities: vec![1.0],
                coverage: 1.0,
                width_sigma: 0.0,
            }),
            SpectrumSpec::Gaussian {
                mean_mev,
                sigma_mev,
            } => {
                if n_bins < 2 {
                    return Err(Error::BadTabulation(
                        "gaussian tabulation needs at least 2 bins",
                    ));
                }
                if !width_sigma.is_finite() || width_sigma <= 0.0 {
                    return Err(Error::BadTabulation("width_sigma must be > 0"));
                }
                let step = 2.0 * width_sigma * sigma_mev / n_bins as f64;
                let lo = mean_mev - width_sigma * sigma_mev;
                let mut energies = Vec::with_capacity(n_bins);
                let mut probabilities = Vec::with_capacity(n_bins);
                for i in 0..n_bins {
                    let a = lo + i as f64 * step;
                    let b = a + step;
                    energies.push((a + b) / 2.0);
                    probabilities.push(
                        normal_cdf((b - mean_mev) / sigma_mev)
                            - normal_cdf((a - mean_mev) / sigma_mev),
                    );
                }
                let coverage: f64 = probabilities.iter().sum();
                Ok(SpectrumTable {
                    energies_mev: energies,
                    probabilities,
                    coverage,
                    width_sigma,
                })
            }
        }
    }
}

/// Discrete rendering of a spectrum for source cards.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumTable {
    /// Line energies \[MeV\], ascending.
    pub energies_mev: Vec<f64>,
    /// Probability per line.
    pub probabilities: Vec<f64>,
    /// Total probability captured by the table (1.0 when nothing is
    /// truncated; `Phi(width_sigma) - Phi(-width_sigma)` for a Gaussian).
    pub coverage: f64,
    /// Half-width of the tabulated window in sigma (0 for monoenergetic).
    pub width_sigma: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// High-accuracy reference for the normal CDF: adaptive Simpson
    /// quadrature of `exp(-t^2/2)/sqrt(2*pi)` — independent of the
    /// A&S 7.1.26 approximation used by `erf`.
    fn reference_normal_cdf(z: f64) -> f64 {
        fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
            let m = (a + b) / 2.0;
            (b - a) / 6.0 * (f(a) + 4.0 * f(m) + f(b))
        }
        fn integrand(t: f64) -> f64 {
            const INV_SQRT_2PI: f64 = 0.398_942_280_401_432_7;
            INV_SQRT_2PI * (-0.5 * t * t).exp()
        }
        // Fixed fine grid; the integrand is smooth and the window modest.
        let n = 4000;
        let step = z / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            sum += simpson(integrand, i as f64 * step, (i + 1) as f64 * step);
        }
        0.5 + sum
    }

    #[test]
    fn erf_matches_quadrature_reference() {
        for &z in &[-4.0, -3.2, -1.0, -0.25, 0.0, 0.5, 1.7, 3.0, 4.0] {
            let got = normal_cdf(z);
            let want = reference_normal_cdf(z);
            assert!((got - want).abs() < 2e-7, "Phi({z}): {got} vs {want}");
        }
    }

    #[test]
    fn gaussian_tabulation_captures_expected_window_mass() {
        // 21 bins over +/-4 sigma: coverage = Phi(4) - Phi(-4).
        let spec = SpectrumSpec::Gaussian {
            mean_mev: 14.0,
            sigma_mev: 0.3,
        };
        let table = spec.tabulate(21, 4.0).unwrap();
        assert_eq!(table.energies_mev.len(), 21);
        assert!((table.coverage - (normal_cdf(4.0) - normal_cdf(-4.0))).abs() < 1e-7);
        assert!(table.coverage < 1.0 && table.coverage > 0.9999);
        let sum: f64 = table.probabilities.iter().sum();
        assert!((sum - table.coverage).abs() < 1e-12);
        // Symmetric window around the mean: probabilities mirror.
        let n = table.probabilities.len();
        for i in 0..n / 2 {
            assert!((table.probabilities[i] - table.probabilities[n - 1 - i]).abs() < 1e-12);
        }
        // Bin centers are ascending and span the window (first center sits
        // half a bin above the lower edge).
        assert!(table.energies_mev.windows(2).all(|w| w[0] < w[1]));
        let half = 4.0 * 0.3;
        let step = 2.0 * half / 21.0;
        assert!((table.energies_mev[0] - (14.0 - half + step / 2.0)).abs() < 1e-12);
        assert!((table.energies_mev[n - 1] - (14.0 + half - step / 2.0)).abs() < 1e-12);
    }

    #[test]
    fn mono_tabulation_is_a_single_unit_line() {
        let table = SpectrumSpec::Mono { energy_mev: 2.4495 }
            .tabulate(21, 4.0)
            .unwrap();
        assert_eq!(table.energies_mev, vec![2.4495]);
        assert_eq!(table.probabilities, vec![1.0]);
        assert_eq!(table.coverage, 1.0);
        assert_eq!(table.width_sigma, 0.0);
    }

    #[test]
    fn bad_tabulation_requests_are_loud() {
        let spec = SpectrumSpec::Gaussian {
            mean_mev: 14.0,
            sigma_mev: 0.3,
        };
        assert_eq!(
            spec.tabulate(1, 4.0),
            Err(Error::BadTabulation(
                "gaussian tabulation needs at least 2 bins"
            ))
        );
        assert_eq!(
            spec.tabulate(21, 0.0),
            Err(Error::BadTabulation("width_sigma must be > 0"))
        );
        assert_eq!(
            spec.tabulate(21, f64::NAN),
            Err(Error::BadTabulation("width_sigma must be > 0"))
        );
    }

    #[test]
    fn zero_temperature_reaction_is_mono() {
        let spec = SpectrumSpec::from_reaction(FusionReaction::Dt, 0.0).unwrap();
        assert_eq!(spec, SpectrumSpec::Mono { energy_mev: 14.021 });
        assert!(spec.is_mono());
    }

    #[test]
    fn proton_spectrum_is_brysk_gaussian_or_mono() {
        // Hand vector: Ti = 20 keV -> (3.02, 0.173782) MeV (reaction.rs).
        let spec = SpectrumSpec::from_proton(20.0).unwrap();
        assert_eq!(
            spec,
            SpectrumSpec::Gaussian {
                mean_mev: 3.02,
                sigma_mev: (1.51 * 0.02_f64).sqrt(),
            }
        );
        let mono = SpectrumSpec::from_proton(0.0).unwrap();
        assert_eq!(mono, SpectrumSpec::Mono { energy_mev: 3.02 });
        assert!(mono.is_mono());
        assert_eq!(
            SpectrumSpec::from_proton(-1.0),
            Err(Error::NegativeIonTemperature(-1.0))
        );
    }
}
