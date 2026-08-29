//! Criterion benchmarks for enrichment cascade solves.

use std::collections::BTreeMap;

use criterion::{criterion_group, criterion_main, Criterion};
use enrichment::{
    default_uranium_cascade, multicomponent, solve_numeric, Cascade, Stream, DEFAULT_MAX_ITER,
    DEFAULT_TOLERANCE,
};
use nuclei::NuclideId;

fn vision_feed_cascade() -> Cascade {
    let mut orig = default_uranium_cascade();
    orig.x_prod_j = 0.055;
    orig.mat_feed = Stream::from_comp(BTreeMap::from([
        (NuclideId::from_name("U234").unwrap(), 0.000183963025893197),
        (NuclideId::from_name("U235").unwrap(), 0.00818576605617839),
        (NuclideId::from_name("U236").unwrap(), 0.00610641667100979),
        (NuclideId::from_name("U238").unwrap(), 0.985523854246919),
    ]));
    orig
}

fn bench_enrichment(c: &mut Criterion) {
    let mut group = c.benchmark_group("enrichment_default_uranium");
    group.sample_size(20);
    group.measurement_time(std::time::Duration::from_secs(5));

    let orig = default_uranium_cascade();
    group.bench_function("solve_numeric_fixed_mstar", |b| {
        b.iter(|| solve_numeric(&orig, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER).expect("solve"))
    });

    group.bench_function("multicomponent_optimize_mstar", |b| {
        b.iter(|| multicomponent(&orig, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER).expect("multi"))
    });

    group.finish();

    let mut group_multi = c.benchmark_group("enrichment_multicomponent_feeds");
    group_multi.sample_size(10);
    group_multi.measurement_time(std::time::Duration::from_secs(10));

    let vision = vision_feed_cascade();
    group_multi.bench_function("vision_4component_feed", |b| {
        b.iter(|| multicomponent(&vision, DEFAULT_TOLERANCE, DEFAULT_MAX_ITER).expect("vision"))
    });

    group_multi.finish();
}

criterion_group!(benches, bench_enrichment);
criterion_main!(benches);
