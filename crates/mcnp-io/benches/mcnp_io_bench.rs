//! Criterion benchmark for MCNP meshtal end-to-end parsing.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use mcnp_io::meshtal::Meshtal;

fn fixture_path(name: &str) -> String {
    format!(
        "{}/../../fixtures/mcnp/meshtal/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn bench_mcnp_io(c: &mut Criterion) {
    let single_path = fixture_path("mcnp_meshtal_single_meshtal.txt");
    let multi_path = fixture_path("mcnp_meshtal_multiple_meshtal.txt");

    let single_text = std::fs::read_to_string(&single_path).expect("read single meshtal");
    let multi_text = std::fs::read_to_string(&multi_path).expect("read multiple meshtal");

    let mut group = c.benchmark_group("mcnp_io_parse_meshtal");
    group.sample_size(30);
    group.measurement_time(std::time::Duration::from_secs(5));

    group.throughput(Throughput::Bytes(single_text.len() as u64));
    group.bench_function("single_tally", |b| {
        b.iter(|| Meshtal::parse(&single_text).expect("parse single"))
    });

    group.throughput(Throughput::Bytes(multi_text.len() as u64));
    group.bench_function("multiple_tallies", |b| {
        b.iter(|| Meshtal::parse(&multi_text).expect("parse multiple"))
    });

    group.bench_function("single_from_file", |b| {
        b.iter(|| Meshtal::from_file(&single_path).expect("from file"))
    });

    group.finish();
}

criterion_group!(benches, bench_mcnp_io);
criterion_main!(benches);
