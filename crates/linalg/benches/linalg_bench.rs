//! Criterion benchmarks for the dense-real kernels added in 0.9.0:
//! `weighted_lstsq` (faer column-pivoted QR) and the UQ-lite samplers
//! (`sample_mvn`, `sample_lognormal`). All inputs are deterministic
//! synthetic values with recorded provenance (no evaluated data).

use criterion::{criterion_group, criterion_main, Criterion};
use nucleide_linalg::lstsq::weighted_lstsq;
use nucleide_linalg::sample::{sample_lognormal, sample_mvn};

/// AR(1) correlation block: `C[i][j] = rho^|i-j|`, positive-definite.
fn ar1_cov(dim: usize, rho: f64, scale: f64) -> Vec<Vec<f64>> {
    (0..dim)
        .map(|i| {
            (0..dim)
                .map(|j| scale * rho.powi((i as i32 - j as i32).abs()))
                .collect()
        })
        .collect()
}

/// Vandermonde-ish design over `xs`: `X[i][j] = xs[i]^j`.
fn vandermonde(xs: &[f64], ncols: usize) -> Vec<Vec<f64>> {
    xs.iter()
        .map(|&x| (0..ncols).map(|j| x.powi(j as i32)).collect())
        .collect()
}

fn bench_linalg(c: &mut Criterion) {
    // Small fit (spectroscopy E7 shape): 12 points, 3 coefficients.
    let xs_small: Vec<f64> = (1..=12).map(|i| (i as f64) * 0.25).collect();
    let x_small = vandermonde(&xs_small, 3);
    let y_small: Vec<f64> = xs_small.iter().map(|&x| 1.0 + 2.0 * x - x * x).collect();
    let w_small = vec![1.0; xs_small.len()];

    // Larger fit: 2000 points, 6 coefficients.
    let xs_big: Vec<f64> = (1..=2000).map(|i| (i as f64) * 0.01).collect();
    let x_big = vandermonde(&xs_big, 6);
    let y_big: Vec<f64> = xs_big.iter().map(|&x| x.sin() + x).collect();
    let w_big = vec![1.0; xs_big.len()];

    let mut group = c.benchmark_group("linalg_weighted_lstsq");
    group.sample_size(20);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.bench_function("small_12x3", |b| {
        b.iter(|| weighted_lstsq(&x_small, &y_small, &w_small).expect("small lstsq"))
    });
    group.bench_function("big_2000x6", |b| {
        b.iter(|| weighted_lstsq(&x_big, &y_big, &w_big).expect("big lstsq"))
    });
    group.finish();

    // UQ-lite samplers: 8-dim AR(1) block, 2000 draws, pinned seed.
    let dim = 8;
    let mean = vec![0.0; dim];
    let cov = ar1_cov(dim, 0.5, 0.04);

    let mut group = c.benchmark_group("linalg_sample");
    group.sample_size(20);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.bench_function("mvn_8d_2000", |b| {
        b.iter(|| sample_mvn(&mean, &cov, 2000, 0xC0FFEE).expect("mvn"))
    });
    group.bench_function("lognormal_8d_2000", |b| {
        b.iter(|| sample_lognormal(&mean, &cov, 2000, 0xC0FFEE).expect("lognormal"))
    });
    group.finish();
}

criterion_group!(benches, bench_linalg);
criterion_main!(benches);
