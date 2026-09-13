//! Criterion benchmark for the E7-fit efficiency coefficient fit
//! (`fit_efficiency`) through the shared `linalg` kernel.
//! Inputs are deterministic synthetic calibration points (no lab data).

use criterion::{criterion_group, criterion_main, Criterion};
use nucleide_spectroscopy::fit_efficiency;

/// Synthetic log-polynomial efficiency law: `eff = exp(-2.8 - 0.73*lnE)`.
fn synthetic_points(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let energies: Vec<f64> = (1..=n).map(|i| 0.1 + (i as f64) * 0.25).collect();
    let effs: Vec<f64> = energies
        .iter()
        .map(|&e| (-2.8 - 0.73 * e.ln()).exp())
        .collect();
    let weights = vec![1.0; n];
    (energies, effs, weights)
}

fn bench_spectroscopy(c: &mut Criterion) {
    let (e12, f12, w12) = synthetic_points(12);
    let (e200, f200, w200) = synthetic_points(200);

    let mut group = c.benchmark_group("spectroscopy_fit_efficiency");
    group.sample_size(20);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.bench_function("fit1_order2_12pt", |b| {
        b.iter(|| fit_efficiency(&e12, &f12, &w12, 2, 1).expect("12pt fit"))
    });
    group.bench_function("fit1_order3_200pt", |b| {
        b.iter(|| fit_efficiency(&e200, &f200, &w200, 3, 1).expect("200pt fit"))
    });
    group.finish();
}

criterion_group!(benches, bench_spectroscopy);
criterion_main!(benches);
