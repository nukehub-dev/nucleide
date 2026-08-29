//! Criterion benchmarks for MAGIC weight-window generation and source sampling.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use mcnp_io::meshtal::{MeshTallyData, ParticleKind};
use vr_tools::{
    magic, magic_with,
    sampling::{AliasTable, MeshSourceSampler, Mode},
    MagicSelection,
};

/// Build a large synthetic structured mesh tally: 50 x 50 x 20 voxels,
/// 10 energy groups, for a total of 500 000 per-group entries.
fn synthetic_tally() -> MeshTallyData {
    let nx = 50;
    let ny = 50;
    let nz = 20;
    let n_eg = 10;
    let n_ves = nx * ny * nz;

    let x_bounds: Vec<f64> = (0..=nx).map(|i| i as f64 * 2.0).collect();
    let y_bounds: Vec<f64> = (0..=ny).map(|i| i as f64 * 2.0).collect();
    let z_bounds: Vec<f64> = (0..=nz).map(|i| i as f64 * 2.0).collect();
    let e_bounds: Vec<f64> = (0..=n_eg).map(|i| 1e-1 * (i as f64).powi(2)).collect();

    let mut result = Vec::with_capacity(n_ves);
    let mut rel_error = Vec::with_capacity(n_ves);
    let mut total_result = Vec::with_capacity(n_ves);
    let mut total_rel_error = Vec::with_capacity(n_ves);

    for ve in 0..n_ves {
        let mut per_group = Vec::with_capacity(n_eg);
        let mut per_err = Vec::with_capacity(n_eg);
        let mut total = 0.0;
        for g in 0..n_eg {
            // Smooth spatial + energy dependence with a non-zero floor.
            let v = 1.0 + ((ve * (g + 1)) % 97) as f64;
            per_group.push(v);
            per_err.push(0.05 + ((ve * (g + 3)) % 13) as f64 * 0.01);
            total += v;
        }
        result.push(per_group);
        rel_error.push(per_err);
        total_result.push(total);
        total_rel_error.push(0.05 + ((ve * 5) % 7) as f64 * 0.01);
    }

    MeshTallyData {
        tally_number: 1,
        particle: ParticleKind::Neutron,
        dose_response: false,
        x_bounds,
        y_bounds,
        z_bounds,
        e_bounds,
        column_idx: Default::default(),
        result,
        rel_error,
        total_result,
        total_rel_error,
    }
}

/// Deterministic LCG returning values in (0, 1).
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = (16807 * self.0) % 2147483647;
        self.0 as f64 / 2147483647.0
    }
}

fn bench_vr_tools(c: &mut Criterion) {
    let tally = synthetic_tally();

    let mut group = c.benchmark_group("vr_tools_magic");
    group.sample_size(20);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.bench_function("magic_total_mode", |b| {
        b.iter(|| magic(&tally).expect("magic total"))
    });
    group.bench_function("magic_per_group_mode", |b| {
        b.iter(|| {
            magic_with(&tally, MagicSelection::PerGroup, Default::default()).expect("magic pg")
        })
    });
    group.finish();

    // Alias-table construction and sampling benchmarks.
    let n = tally.num_ves();
    let analog_pdf: Vec<f64> = tally
        .total_result
        .iter()
        .enumerate()
        .map(|(i, q)| {
            q * (tally.x_bounds[(i / (tally.dims()[1] * tally.dims()[2])) + 1]
                - tally.x_bounds[i / (tally.dims()[1] * tally.dims()[2])])
                * (tally.y_bounds[(i / tally.dims()[2]) % tally.dims()[1] + 1]
                    - tally.y_bounds[(i / tally.dims()[2]) % tally.dims()[1]])
                * (tally.z_bounds[(i % tally.dims()[2]) + 1] - tally.z_bounds[i % tally.dims()[2]])
        })
        .collect();

    let mut group_table = c.benchmark_group("vr_tools_alias_table");
    group_table.sample_size(30);
    group_table.measurement_time(std::time::Duration::from_secs(5));
    group_table.bench_function("build_alias_table", |b| {
        b.iter(|| AliasTable::new(&analog_pdf).expect("alias table"))
    });
    group_table.finish();

    let sampler_analog =
        MeshSourceSampler::new(&tally, Mode::Analog, None).expect("analog sampler");
    let sampler_uniform =
        MeshSourceSampler::new(&tally, Mode::Uniform, None).expect("uniform sampler");
    let user_pdf: Vec<f64> = (0..n).map(|i| 1.0 + (i % 5) as f64).collect();
    let sampler_user =
        MeshSourceSampler::new(&tally, Mode::User, Some(&user_pdf)).expect("user sampler");

    let mut group_sample = c.benchmark_group("vr_tools_sampling_1m");
    group_sample.throughput(Throughput::Elements(1_000_000));
    group_sample.sample_size(10);
    group_sample.measurement_time(std::time::Duration::from_secs(8));

    group_sample.bench_function("analog_1m_samples", |b| {
        b.iter(|| {
            let mut rng = Lcg(123_456_789);
            for _ in 0..1_000_000 {
                let _ = sampler_analog.sample(rng.next_f64(), rng.next_f64());
            }
        })
    });

    group_sample.bench_function("uniform_1m_samples", |b| {
        b.iter(|| {
            let mut rng = Lcg(123_456_789);
            for _ in 0..1_000_000 {
                let _ = sampler_uniform.sample(rng.next_f64(), rng.next_f64());
            }
        })
    });

    group_sample.bench_function("user_1m_samples", |b| {
        b.iter(|| {
            let mut rng = Lcg(123_456_789);
            for _ in 0..1_000_000 {
                let _ = sampler_user.sample(rng.next_f64(), rng.next_f64());
            }
        })
    });

    group_sample.finish();
}

criterion_group!(benches, bench_vr_tools);
criterion_main!(benches);
