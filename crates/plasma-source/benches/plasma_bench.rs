//! Criterion benchmarks for the parametric plasma source.
//!
//! ITER-ish synthetic H-mode configuration (hand round numbers, same as the
//! crate's gate fixtures): one-time weight-table build plus steady-state
//! particle sampling.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use nucleide_plasma_source::miller::MillerGeometry;
use nucleide_plasma_source::parametric::{ParametricPlasmaConfig, ParametricSampler};
use nucleide_plasma_source::profile::{DensityProfile, ProfileMode, TemperatureProfile};
use nucleide_plasma_source::reaction::FusionReaction;

/// ITER-ish synthetic H-mode case (mirrors the crate's gate fixture).
fn iter_h_mode() -> ParametricPlasmaConfig {
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

fn bench_plasma(c: &mut Criterion) {
    let config = iter_h_mode();

    let mut group = c.benchmark_group("plasma_source_parametric");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(8));

    // One-time O(R x THETA) table build (65 536 volume elements).
    group.bench_function("weight_table_build", |b| {
        b.iter(|| ParametricSampler::new(config, 20260915).expect("bench sampler"))
    });

    // Steady-state sampling throughput.
    let mut sampler = ParametricSampler::new(config, 20260915).expect("bench sampler");
    group.throughput(Throughput::Elements(10_000));
    group.bench_function("sample_10k", |b| b.iter(|| sampler.sample_n(10_000)));

    group.finish();
}

criterion_group!(benches, bench_plasma);
criterion_main!(benches);
