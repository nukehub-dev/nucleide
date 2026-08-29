//! Criterion benchmarks for nuclide identification and data lookups.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use nuclei::{
    data::{atomic_mass, atomic_mass_by_name, half_life, half_life_by_name},
    dialects::{alara_to_id, fluka_to_id, from_cinder, from_serpent, from_zaid, nist_to_id},
    NuclideId,
};

const NAMES: &[&str] = &[
    "H1", "He4", "C12", "O16", "Fe56", "U235", "U238", "Pu239", "Am242_m1", "Xe135", "Cs137",
    "Ba137m", "Co60", "Ni58", "I135", "Gd157", "Mo95", "Zr90", "Nb95", "Ru106",
];

const ZAID_IDS: &[u32] = &[92235, 92238, 94239, 95242, 1001, 2004];

const CINDER_IDS: &[u32] = &[2350920, 2380920, 2390940, 2420951];

const SERPENT_NAMES: &[&str] = &["U-235", "Pu-239", "Am-242m", "Xe-135", "Cs-137"];

const NIST_NAMES: &[&str] = &["235U", "239Pu", "242Am", "137Cs"];

const ALARA_NAMES: &[&str] = &["u:235", "pu:239", "am:242", "cs:137"];

const FLUKA_NAMES: &[&str] = &["235-U", "238-U", "HELIUM-4", "LITHIU-7", "BORON-10"];

fn bench_nuclei(c: &mut Criterion) {
    let mut group = c.benchmark_group("nuclei_from_name");
    group.throughput(Throughput::Elements(NAMES.len() as u64));
    group.bench_function("parse_20_names", |b| {
        b.iter(|| {
            for name in NAMES {
                let _ = NuclideId::from_name(name).expect("parse");
            }
        })
    });
    group.finish();

    let mut group_dialect = c.benchmark_group("nuclei_dialects");
    group_dialect.throughput(Throughput::Elements(
        (ZAID_IDS.len()
            + CINDER_IDS.len()
            + SERPENT_NAMES.len()
            + NIST_NAMES.len()
            + ALARA_NAMES.len()
            + FLUKA_NAMES.len()) as u64,
    ));
    group_dialect.bench_function("parse_all_dialects", |b| {
        b.iter(|| {
            for &z in ZAID_IDS {
                let _ = from_zaid(z).expect("zaid");
            }
            for &c in CINDER_IDS {
                let _ = from_cinder(c).expect("cinder");
            }
            for name in SERPENT_NAMES {
                let _ = from_serpent(name).expect("serpent");
            }
            for name in NIST_NAMES {
                let _ = nist_to_id(name).expect("nist");
            }
            for name in ALARA_NAMES {
                let _ = alara_to_id(name).expect("alara");
            }
            for name in FLUKA_NAMES {
                let _ = fluka_to_id(name).expect("fluka");
            }
        })
    });
    group_dialect.finish();

    let mut group_lookups = c.benchmark_group("nuclei_data_lookups");
    group_lookups.throughput(Throughput::Elements(NAMES.len() as u64));

    group_lookups.bench_function("atomic_mass_by_name", |b| {
        b.iter(|| {
            for name in NAMES {
                let _ = atomic_mass_by_name(name);
            }
        })
    });

    group_lookups.bench_function("atomic_mass_by_nucid", |b| {
        let ids: Vec<u32> = NAMES
            .iter()
            .map(|n| NuclideId::from_name(n).unwrap().nucid())
            .collect();
        b.iter(|| {
            for &id in &ids {
                let _ = atomic_mass(id);
            }
        })
    });

    group_lookups.bench_function("half_life_by_name", |b| {
        b.iter(|| {
            for name in NAMES {
                let _ = half_life_by_name(name);
            }
        })
    });

    group_lookups.bench_function("half_life_by_nucid", |b| {
        let ids: Vec<u32> = NAMES
            .iter()
            .map(|n| NuclideId::from_name(n).unwrap().nucid())
            .collect();
        b.iter(|| {
            for &id in &ids {
                let _ = half_life(id);
            }
        })
    });

    group_lookups.finish();
}

criterion_group!(benches, bench_nuclei);
criterion_main!(benches);
