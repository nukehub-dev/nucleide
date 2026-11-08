//! Gaussian kernel-density source sampling (KDSource-class, clean-room).
//!
//! [`KdeSampler`] fits an axis-aligned Gaussian KDE over caller-supplied
//! particle vectors (energy/position/direction rows, e.g. projected from
//! MCPL records by the caller) and resamples synthetic particles with the
//! same density to boost downstream statistics. The KDE equations
//! (Gaussian product kernel, Silverman bandwidth rule) are textbook
//! statistics implemented here from scratch; no upstream code or data is
//! read, ported, or vendored.
//!
//! The adaptive per-sample rule ([`Bandwidth::AdaptiveKnn`], KDSource
//! §2/§3.2 item 3-class) replaces the single global width per axis with
//! one width per fitted sample: `h_i0 = c · dist_k(x_i)` is the
//! `k`-th-nearest-neighbour (leave-one-out) distance in `σ_j`-scaled
//! coordinates, scaled by the caller factor `c`, and the per-sample
//! per-dimension width is `h_ij = h_i0 · σ_j` with `σ_j` the unbiased
//! per-dimension standard deviation. The density is then
//! `pdf(x) = (1/N) Σ_i ∏_j N(x_j; p_ij, h_ij²)` and [`KdeSampler::draw`]
//! perturbs the selected centre row by `h_ij · z_j` per dimension — the
//! same deterministic shape as the Silverman path, with per-sample widths.
//! Hand vector: samples `0, 0.2, 0.4, 10` (`d = 1`, `n = 4`), `k = 1`,
//! `scale = 1` gives `h = 0.2, 0.2, 0.2, 9.6` (the cluster is tight, the
//! tail is smooth) with `σ = 4.9027203333115645`, hence scaled base
//! widths `h_i0 = 0.04079367910119179` (cluster) and `1.9580965968572057`
//! (tail); `pdf(0) = 0.8746688805236177` (term by term
//! `1.9947114020071635 + 1.2098536225957168 + 0.2699548325659403 +
//! 0.024155664925650286`, all over `n = 4`) and
//! `draw(0.5, [1.0]) = 0.4 + 0.2 ≈ 0.6` (centre index 2).
//!
//! Draws are deterministic in the house style ([`crate::sampling`]:
//! caller-supplied randoms, no RNG inside): [`KdeSampler::draw`] takes one
//! uniform `u` selecting the kernel centre plus one standard normal per
//! dimension. MCPL projection stays caller-side (a ~10-line map from
//! particle fields to rows); this module never depends on `mcpl-io`.

use crate::{Error, Result};

/// sqrt(2π): Gaussian normalization factor (no std const exists).
const SQRT_2PI: f64 = 2.506_628_274_631_000_2;

/// Bandwidth rule for [`KdeSampler::fit`].
#[derive(Debug, Clone, PartialEq)]
pub enum Bandwidth {
    /// Silverman's rule per dimension,
    /// `h = σ (4 / (d + 2) / n)^(1 / (d + 4))` with the unbiased
    /// per-dimension standard deviation `σ`. Zero-variance dimensions are
    /// an [`Error`] (a zero bandwidth is a delta spike, never a density);
    /// use [`Bandwidth::Fixed`] with an explicit width instead.
    Silverman,
    /// Caller-supplied per-dimension widths (all finite and positive).
    Fixed(Vec<f64>),
    /// Adaptive per-sample k-nearest-neighbour widths (KDSource §2/§3.2
    /// item 3-class): `h_i0 = scale · dist_k(x_i)` in `σ_j`-scaled
    /// coordinates, `h_ij = h_i0 · σ_j` per dimension. `dist_k(x_i)` is
    /// the Euclidean distance from sample `i` to its `k`-th nearest
    /// neighbour, excluding itself — `k` is 1-based, so `dist_1` is the
    /// nearest neighbour and `dist_k` is undefined unless
    /// `1 <= k <= n - 1` (loud [`Error::BadAdaptiveKnn`] otherwise).
    /// `scale` (the factor `c`, conventionally `1`) must be finite and
    /// positive. Zero-variance dimensions (loud, as in
    /// [`Bandwidth::Silverman`]) and a zero `dist_k` from duplicate
    /// points (loud [`Error::ZeroAdaptiveBandwidth`]) are both rejected:
    /// never a delta spike.
    AdaptiveKnn {
        /// 1-based neighbour rank (`dist_1` = nearest neighbour, self
        /// excluded).
        k: usize,
        /// Caller scale factor `c` (must be finite and positive).
        scale: f64,
    },
}

/// Axis-aligned Gaussian KDE over `n` samples in `d` dimensions.
///
/// `widths` holds the fit-rule parameters: the global per-dimension
/// vector (length `d`) for [`Bandwidth::Silverman`] /
/// [`Bandwidth::Fixed`], or the per-sample scaled base widths `h_i0`
/// (length `n`) for [`Bandwidth::AdaptiveKnn`]. `sample_widths` is the
/// expanded `n × d` matrix of physical per-sample bandwidths `h_ij`
/// for every rule; Silverman/Fixed rows are all identical to
/// `widths`, while AdaptiveKnn rows carry the per-sample widths.
#[derive(Debug, Clone, PartialEq)]
pub struct KdeSampler {
    samples: Vec<Vec<f64>>,
    widths: Vec<f64>,
    sample_widths: Vec<Vec<f64>>,
}

/// Unbiased per-dimension standard deviation `σ_j` over `samples`.
/// A zero-variance dimension is an [`Error`] (a zero bandwidth is a
/// delta spike, never a density).
fn sigma_per_dim(samples: &[Vec<f64>], dims: usize, n: f64) -> Result<Vec<f64>> {
    let mut sigmas = Vec::with_capacity(dims);
    for j in 0..dims {
        let mean = samples.iter().map(|row| row[j]).sum::<f64>() / n;
        let var = samples
            .iter()
            .map(|row| (row[j] - mean).powi(2))
            .sum::<f64>()
            / (n - 1.0).max(1.0);
        let sigma = var.sqrt();
        if sigma <= 0.0 {
            return Err(Error::ZeroVarianceDim { dim: j });
        }
        sigmas.push(sigma);
    }
    Ok(sigmas)
}

impl KdeSampler {
    /// Fit a KDE over `samples` (non-empty, rectangular, all finite).
    pub fn fit(samples: &[Vec<f64>], bandwidth: Bandwidth) -> Result<Self> {
        if samples.is_empty() {
            return Err(Error::EmptyPdf);
        }
        let dims = samples[0].len();
        if dims == 0 {
            return Err(Error::LengthMismatch {
                expected: 1,
                got: 0,
            });
        }
        for (i, row) in samples.iter().enumerate() {
            if row.len() != dims {
                return Err(Error::LengthMismatch {
                    expected: dims,
                    got: row.len(),
                });
            }
            if row.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFiniteTally {
                    field: "kde sample",
                    index: i,
                });
            }
        }
        let n_usize = samples.len();
        let n = n_usize as f64;
        let d = dims as f64;
        let (widths, sample_widths) = match bandwidth {
            Bandwidth::Fixed(w) => {
                if w.len() != dims {
                    return Err(Error::LengthMismatch {
                        expected: dims,
                        got: w.len(),
                    });
                }
                if w.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                    return Err(Error::NonFiniteTally {
                        field: "kde bandwidth",
                        index: w
                            .iter()
                            .position(|v| !v.is_finite() || *v <= 0.0)
                            .unwrap_or(0),
                    });
                }
                let sample_widths = vec![w.clone(); n_usize];
                (w, sample_widths)
            }
            Bandwidth::Silverman => {
                let factor = (4.0 / (d + 2.0) / n).powf(1.0 / (d + 4.0));
                let sigmas = sigma_per_dim(samples, dims, n)?;
                let widths = sigmas.iter().map(|s| s * factor).collect::<Vec<_>>();
                let sample_widths = vec![widths.clone(); n_usize];
                (widths, sample_widths)
            }
            Bandwidth::AdaptiveKnn { k, scale } => {
                if k == 0 || k > n_usize.saturating_sub(1) || !scale.is_finite() || scale <= 0.0 {
                    return Err(Error::BadAdaptiveKnn {
                        n: n_usize,
                        k,
                        scale,
                    });
                }
                let sigmas = sigma_per_dim(samples, dims, n)?;
                let mut base_widths = Vec::with_capacity(n_usize);
                let mut sample_widths = Vec::with_capacity(n_usize);
                for (i, row) in samples.iter().enumerate() {
                    let mut dists = Vec::with_capacity(n_usize - 1);
                    for (m, other) in samples.iter().enumerate() {
                        if m == i {
                            continue;
                        }
                        let mut d2 = 0.0;
                        for (j, (&a, &b)) in row.iter().zip(other).enumerate() {
                            let z = (a - b) / sigmas[j];
                            d2 += z * z;
                        }
                        dists.push(d2.sqrt());
                    }
                    dists.sort_by(|a, b| a.total_cmp(b));
                    let dist_k = dists[k - 1];
                    if dist_k <= 0.0 {
                        return Err(Error::ZeroAdaptiveBandwidth { sample: i });
                    }
                    let base = scale * dist_k;
                    base_widths.push(base);
                    sample_widths.push(sigmas.iter().map(|s| base * s).collect());
                }
                (base_widths, sample_widths)
            }
        };
        Ok(Self {
            samples: samples.to_vec(),
            widths,
            sample_widths,
        })
    }

    /// Number of fitted samples.
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    /// Sample dimensionality.
    pub fn dims(&self) -> usize {
        self.samples[0].len()
    }

    /// Fitted bandwidths in the fit-rule parameterization: the global
    /// per-dimension vector (length [`KdeSampler::dims`]) for
    /// [`Bandwidth::Silverman`] / [`Bandwidth::Fixed`], or the per-sample
    /// scaled base widths `h_i0` (length [`KdeSampler::n_samples`]) for
    /// [`Bandwidth::AdaptiveKnn`]. Use
    /// [`KdeSampler::sample_bandwidths`] for the physical per-sample `h_ij`.
    pub fn bandwidths(&self) -> &[f64] {
        &self.widths
    }

    /// Per-sample physical bandwidths `h_ij`: one row per fitted sample,
    /// one entry per dimension. Silverman/Fixed rows are all identical to
    /// [`KdeSampler::bandwidths`]; AdaptiveKnn rows carry the per-sample
    /// widths `h_i0 · σ_j`.
    pub fn sample_bandwidths(&self) -> &[Vec<f64>] {
        &self.sample_widths
    }

    /// KDE density at `point` (Gaussian product kernel with per-sample
    /// widths, `pdf(x) = (1/N) Σ_i ∏_j N(x_j; p_ij, h_ij²)`,
    /// length-checked).
    pub fn pdf(&self, point: &[f64]) -> Result<f64> {
        if point.len() != self.dims() {
            return Err(Error::LengthMismatch {
                expected: self.dims(),
                got: point.len(),
            });
        }
        if point.iter().any(|v| !v.is_finite()) {
            return Err(Error::NonFiniteTally {
                field: "kde point",
                index: point.iter().position(|v| !v.is_finite()).unwrap_or(0),
            });
        }
        let mut density = 0.0;
        for (row, widths) in self.samples.iter().zip(&self.sample_widths) {
            let norm = widths.iter().map(|h| h * SQRT_2PI).product::<f64>();
            let mut z2 = 0.0;
            for ((x, c), h) in point.iter().zip(row).zip(widths) {
                let z = (x - c) / h;
                z2 += z * z;
            }
            density += (-0.5 * z2).exp() / norm;
        }
        Ok(density / self.samples.len() as f64)
    }

    /// Resample one synthetic particle: `u` in `[0, 1)` selects kernel
    /// centre `floor(u * n)`; `normals` (one finite standard normal per
    /// dimension) perturbs it by `width * z` per dimension, using the
    /// centre row's per-sample widths.
    pub fn draw(&self, u: f64, normals: &[f64]) -> Result<Vec<f64>> {
        if !(0.0..1.0).contains(&u) {
            return Err(Error::BadDraw { value: u });
        }
        if normals.len() != self.dims() {
            return Err(Error::LengthMismatch {
                expected: self.dims(),
                got: normals.len(),
            });
        }
        if normals.iter().any(|v| !v.is_finite()) {
            return Err(Error::NonFiniteTally {
                field: "kde normal",
                index: normals.iter().position(|v| !v.is_finite()).unwrap_or(0),
            });
        }
        let idx = (u * self.samples.len() as f64) as usize;
        let centre = &self.samples[idx];
        let widths = &self.sample_widths[idx];
        Ok(centre
            .iter()
            .zip(widths)
            .zip(normals)
            .map(|((c, h), z)| c + h * z)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic test stream (LCG + Box-Muller); test-only, never shipped.
    struct TestRng(u64);

    impl TestRng {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }

        fn normal(&mut self) -> f64 {
            let (u1, u2) = (self.uniform().max(1e-300), self.uniform());
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }
    }

    fn gaussian_samples(mean: f64, sigma: f64, n: usize) -> Vec<Vec<f64>> {
        let mut rng = TestRng(0x1234_5678_9abc_def0);
        (0..n).map(|_| vec![mean + sigma * rng.normal()]).collect()
    }

    #[test]
    fn input_errors_are_loud() {
        assert!(KdeSampler::fit(&[], Bandwidth::Silverman).is_err());
        assert!(KdeSampler::fit(&[vec![]], Bandwidth::Silverman).is_err());
        assert!(KdeSampler::fit(&[vec![1.0], vec![1.0, 2.0]], Bandwidth::Silverman).is_err());
        assert!(KdeSampler::fit(&[vec![f64::NAN]], Bandwidth::Silverman).is_err());
        assert!(KdeSampler::fit(&[vec![1.0]], Bandwidth::Fixed(vec![])).is_err());
        assert!(KdeSampler::fit(&[vec![1.0]], Bandwidth::Fixed(vec![0.0])).is_err());
        assert!(KdeSampler::fit(&[vec![1.0]], Bandwidth::Fixed(vec![-2.0])).is_err());
        // Zero variance under Silverman is a delta spike, not a density.
        assert!(KdeSampler::fit(&vec![vec![3.0]; 8], Bandwidth::Silverman).is_err());
        // ...but an explicit fixed width is the caller's choice.
        assert!(KdeSampler::fit(&vec![vec![3.0]; 8], Bandwidth::Fixed(vec![0.5])).is_ok());
    }

    #[test]
    fn draw_is_exact_and_deterministic() {
        let kde = KdeSampler::fit(
            &[vec![1.0, 2.0], vec![3.0, 4.0]],
            Bandwidth::Fixed(vec![0.5, 2.0]),
        )
        .unwrap();
        // u = 0.75 selects centre 1; perturbation is width * z per dim.
        let got = kde.draw(0.75, &[1.0, -0.5]).unwrap();
        assert_eq!(got, vec![3.5, 3.0]);
        assert_eq!(kde.draw(0.75, &[1.0, -0.5]).unwrap(), got);
        assert!(kde.draw(1.0, &[0.0, 0.0]).is_err());
        assert!(kde.draw(-0.1, &[0.0, 0.0]).is_err());
        assert!(kde.draw(0.5, &[0.0]).is_err());
        assert!(kde.draw(0.5, &[0.0, f64::INFINITY]).is_err());
    }

    #[test]
    fn gaussian_recovery_and_normalization() {
        // 20000 draws from N(5, 2²): KDE mean within 5 SE, pdf at the mode
        // within 2% of the closed form, integral within 1e-3 of 1.
        let samples = gaussian_samples(5.0, 2.0, 20_000);
        let kde = KdeSampler::fit(&samples, Bandwidth::Silverman).unwrap();
        let h = kde.bandwidths()[0];
        assert!(h > 0.0 && h < 1.0, "silverman width {h}");
        let mut draws = Vec::with_capacity(4096);
        let mut rng = TestRng(0xabcd);
        for _ in 0..4096 {
            draws.push(kde.draw(rng.uniform(), &[rng.normal()]).unwrap()[0]);
        }
        let mean = draws.iter().sum::<f64>() / draws.len() as f64;
        let se = 2.0 / (draws.len() as f64).sqrt();
        assert!((mean - 5.0).abs() < 5.0 * se, "mean {mean}");
        let closed = (-0.5 * ((5.0f64 - 5.0) / 2.0).powi(2)).exp() / (2.0 * SQRT_2PI);
        let got = kde.pdf(&[5.0]).unwrap();
        assert!(
            (got - closed).abs() / closed < 0.02,
            "pdf {got} vs {closed}"
        );
        // Trapezoid integral over ±8σ.
        let (lo, hi, m) = (-11.0, 21.0, 2048);
        let mut area = 0.0;
        let mut prev = kde.pdf(&[lo]).unwrap();
        for i in 1..=m {
            let x = lo + (hi - lo) * i as f64 / m as f64;
            let cur = kde.pdf(&[x]).unwrap();
            area += 0.5 * (prev + cur) * (hi - lo) / m as f64;
            prev = cur;
        }
        assert!((area - 1.0).abs() < 1e-3, "integral {area}");
    }

    #[test]
    fn bandwidth_is_deterministic() {
        let samples = gaussian_samples(0.0, 1.0, 512);
        let a = KdeSampler::fit(&samples, Bandwidth::Silverman).unwrap();
        let b = KdeSampler::fit(&samples, Bandwidth::Silverman).unwrap();
        assert_eq!(a.bandwidths(), b.bandwidths());
    }

    /// Sparse-cluster hand vector: `{0, 0.2, 0.4, 10}`, `k = 1`,
    /// `scale = 1`. The unbiased `σ = sqrt(72.11/3) = 4.9027203333115645`
    /// cancels in one dimension, so the physical widths are the
    /// leave-one-out nearest-neighbour distances `0.2, 0.2, 0.2, 9.6`
    /// (tight cluster smooths locally, the tail smooths broadly); the
    /// scaled base widths are `0.2/σ` and `9.6/σ`. At `x = 0` the density
    /// is the four-term mean
    /// `(1.9947114020071635 + 1.2098536225957168 + 0.2699548325659403 +
    /// 0.024155664925650286) / 4 = 0.8746688805236177`, and `u = 0.5`
    /// selects centre index 2 (`0.4`), so `draw(0.5, [1.0])` is
    /// `0.4 + h_2 · 1.0 ≈ 0.6` bit-exactly against the stored row.
    #[test]
    fn adaptive_knn_hand_vector() {
        let samples = [vec![0.0], vec![0.2], vec![0.4], vec![10.0]];
        let kde = KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 }).unwrap();
        assert_eq!(kde.n_samples(), 4);
        assert_eq!(kde.dims(), 1);
        let sw = kde.sample_bandwidths();
        assert_eq!(sw.len(), 4);
        assert!((sw[0][0] - 0.2).abs() < 1e-12, "h0 = {}", sw[0][0]);
        assert!((sw[1][0] - 0.2).abs() < 1e-12, "h1 = {}", sw[1][0]);
        assert!((sw[2][0] - 0.2).abs() < 1e-12, "h2 = {}", sw[2][0]);
        assert!((sw[3][0] - 9.6).abs() < 1e-12, "h3 = {}", sw[3][0]);
        assert_ne!(sw[0][0], sw[3][0]);
        // Scaled base widths `h_i0 = 0.2/σ` (cluster) and `9.6/σ` (tail).
        let base = kde.bandwidths();
        assert!(
            (base[0] - 0.04079367910119179).abs() < 1e-12,
            "base0 = {base:?}"
        );
        assert!(
            (base[3] - 1.9580965968572057).abs() < 1e-12,
            "base3 = {base:?}"
        );
        let pdf0 = kde.pdf(&[0.0]).unwrap();
        assert!((pdf0 - 0.8746688805236177).abs() < 1e-9, "pdf(0) = {pdf0}");
        // Centre selection + per-sample row: bit-exact against the stored
        // `h_2`, and ≈ 0.6 nominally.
        let got = kde.draw(0.5, &[1.0]).unwrap();
        assert_eq!(got[0], 0.4 + sw[2][0] * 1.0);
        assert!((got[0] - 0.6).abs() < 1e-12, "draw = {}", got[0]);
    }

    /// `k = 2` picks the second-nearest neighbour (self excluded), and
    /// `scale` multiplies every `h_i0` uniformly (hence every `h_ij`).
    #[test]
    fn adaptive_knn_rank_and_scale() {
        let samples = [vec![0.0], vec![0.2], vec![0.4], vec![10.0]];
        let k2 = KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 2, scale: 1.0 }).unwrap();
        let sw = k2.sample_bandwidths();
        assert!((sw[0][0] - 0.4).abs() < 1e-12, "h0 = {}", sw[0][0]);
        assert!((sw[1][0] - 0.2).abs() < 1e-12, "h1 = {}", sw[1][0]);
        assert!((sw[2][0] - 0.4).abs() < 1e-12, "h2 = {}", sw[2][0]);
        assert!((sw[3][0] - 9.8).abs() < 1e-12, "h3 = {}", sw[3][0]);
        let scaled =
            KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 1, scale: 2.0 }).unwrap();
        let sw2 = scaled.sample_bandwidths();
        assert!((sw2[0][0] - 0.4).abs() < 1e-12, "h0 = {}", sw2[0][0]);
        assert!((sw2[3][0] - 19.2).abs() < 1e-12, "h3 = {}", sw2[3][0]);
    }

    /// 2-D correlated pairs on the line `y = x`: the joint centre row
    /// carries the correlation, so every `h_ij = √2` (scaled distance
    /// `√2/σ` with `σ = sqrt(17/3)` per dimension), the density is higher
    /// on the line than off it, and a joint draw stays on the line.
    #[test]
    fn adaptive_knn_2d_correlated_pairs() {
        let samples = [
            vec![0.0, 0.0],
            vec![1.0, 1.0],
            vec![4.0, 4.0],
            vec![5.0, 5.0],
        ];
        let kde = KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 }).unwrap();
        for row in kde.sample_bandwidths() {
            assert!(
                (row[0] - std::f64::consts::SQRT_2).abs() < 1e-12,
                "h_x = {row:?}"
            );
            assert!(
                (row[1] - std::f64::consts::SQRT_2).abs() < 1e-12,
                "h_y = {row:?}"
            );
        }
        assert!(
            kde.pdf(&[2.5, 2.5]).unwrap() > kde.pdf(&[2.5, 5.0]).unwrap(),
            "density must be elevated on the correlation line"
        );
        let got = kde.draw(0.0, &[1.0, 1.0]).unwrap();
        assert_eq!(got[0], got[1]);
        assert!((got[0] - std::f64::consts::SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn adaptive_knn_loud_errors() {
        let samples = [vec![0.0], vec![0.2], vec![0.4], vec![10.0]];
        // k = 0, k > n - 1, and non-positive / non-finite scales.
        assert!(matches!(
            KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 0, scale: 1.0 }),
            Err(Error::BadAdaptiveKnn { .. })
        ));
        assert!(matches!(
            KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 4, scale: 1.0 }),
            Err(Error::BadAdaptiveKnn { .. })
        ));
        assert!(matches!(
            KdeSampler::fit(&[vec![1.0]], Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 }),
            Err(Error::BadAdaptiveKnn { .. })
        ));
        for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    KdeSampler::fit(&samples, Bandwidth::AdaptiveKnn { k: 1, scale }),
                    Err(Error::BadAdaptiveKnn { .. })
                ),
                "scale {scale} must be loud"
            );
        }
        // Duplicate points: `dist_1 = 0` is a delta spike, never loud-free.
        assert!(matches!(
            KdeSampler::fit(
                &[vec![0.0], vec![0.0], vec![1.0]],
                Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 }
            ),
            Err(Error::ZeroAdaptiveBandwidth { .. })
        ));
        // Zero-variance dimensions stay loud, as in the Silverman path.
        assert!(matches!(
            KdeSampler::fit(
                &[vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]],
                Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 }
            ),
            Err(Error::ZeroVarianceDim { dim: 0 })
        ));
    }

    #[test]
    fn adaptive_knn_is_deterministic() {
        let samples = [vec![0.0], vec![0.2], vec![0.4], vec![10.0]];
        let rule = || Bandwidth::AdaptiveKnn { k: 1, scale: 1.0 };
        let a = KdeSampler::fit(&samples, rule()).unwrap();
        let b = KdeSampler::fit(&samples, rule()).unwrap();
        assert_eq!(a.bandwidths(), b.bandwidths());
        assert_eq!(a.sample_bandwidths(), b.sample_bandwidths());
    }
}
