//! Criterion benchmarks for the activation-to-transport source contract.
//!
//! Synthetic regular meshes with 24 groups and irregular per-voxel spectra.
//! Measures source construction/validation, structured reporting, density
//! conversion and joint probabilities separately at 1k, 10k and 100k
//! voxels (2.4M values at the top size). No I/O and no OpenMC objects.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use nucleide_r2s::{PhotonSource, PhotonSourceDraft, SourceReport};

const GROUPS: usize = 24;

fn source(voxels: usize) -> PhotonSource {
    let mut strengths = Vec::with_capacity(voxels);
    for voxel in 0..voxels {
        let mut row = Vec::with_capacity(GROUPS);
        for group in 0..GROUPS {
            let value = if (voxel + group) % 7 == 0 {
                0.0
            } else {
                ((voxel * 31 + group * 17) % 1000) as f64 + 0.25
            };
            row.push(value);
        }
        strengths.push(row);
    }
    let mut edges = Vec::with_capacity(GROUPS + 1);
    for group in 0..=GROUPS {
        edges.push(group as f64 * (group as f64 + 1.0));
    }
    PhotonSource::new(
        [voxels, 1, 1],
        [0.0, 0.0, 0.0],
        [voxels as f64, 1.0, 1.0],
        vec![1.0; voxels],
        edges,
        0.0,
        strengths,
    )
    .expect("synthetic source")
}

fn bench_source(c: &mut Criterion) {
    let mut group = c.benchmark_group("r2s_source_contract");
    group.sample_size(10);

    for voxels in [1_000usize, 10_000, 100_000] {
        let built = source(voxels);
        group.throughput(Throughput::Elements((voxels * GROUPS) as u64));

        group.bench_function(format!("new_{voxels}"), |b| b.iter(|| source(voxels)));
        group.bench_function(format!("report_{voxels}"), |b| {
            let draft = PhotonSourceDraft::new(
                built.mesh_shape,
                built.lower_left_cm,
                built.upper_right_cm,
                built.voxel_volumes_cm3.clone(),
                built.edges_ev.clone(),
                built.cooling_s,
                built.strengths_photons_per_s.clone(),
            );
            b.iter(|| SourceReport::from_draft(&draft))
        });
        group.bench_function(format!("density_{voxels}"), |b| {
            b.iter(|| built.density_photons_per_s_cm3().expect("density"))
        });
        group.bench_function(format!("joint_probabilities_{voxels}"), |b| {
            b.iter(|| built.joint_probabilities().expect("probabilities"))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_source);
criterion_main!(benches);
