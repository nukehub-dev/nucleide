//! Criterion benchmarks for depletion matrix construction and CRAM solves.

use std::collections::BTreeMap;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use depletion::{
    chain::{Chain, ChainNuclide, DecayMode, Reaction},
    matrix::{DepletionSystem, ReactionRates},
    Order,
};

fn fixture_path(name: &str) -> String {
    format!(
        "{}/../../fixtures/depletion/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// `chain_ni.xml` is a fragment with a few `<decay>` entries lacking a
/// `target` attribute.  The benchmark patches in realistic daughters and
/// appends them as stable nuclides so the chain can be built into a
/// [`DepletionSystem`] without modifying the library parser or the fixture.
fn patched_chain_ni() -> Chain {
    let text = std::fs::read_to_string(fixture_path("chain_ni.xml")).expect("read chain_ni.xml");

    let daughters: std::collections::HashMap<&str, &str> = [
        ("Fe55", "Mn55"),
        ("Fe59", "Co59"),
        ("Ni57", "Co57"),
        ("Ni59", "Co59"),
        ("Ni63", "Cu63"),
        ("Ni65", "Cu65"),
    ]
    .into_iter()
    .collect();

    let mut current = String::new();
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        if let Some(start) = line.find("<nuclide name=\"") {
            if let Some(end) = line[start + 15..].find('"') {
                current = line[start + 15..start + 15 + end].to_string();
            }
        }
        if line.contains("<decay") && !line.contains("target=") {
            let daughter = daughters
                .get(current.as_str())
                .copied()
                .unwrap_or("Nothing");
            out.push_str(&line.replace("/>", &format!(" target=\"{daughter}\" />")));
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }

    // Append stable daughter nuclides that are missing from the fragment.
    let mut extra = String::new();
    let mut seen = std::collections::HashSet::new();
    for &daughter in daughters.values() {
        if daughter != "Nothing"
            && !out.contains(&format!("name=\"{daughter}\""))
            && seen.insert(daughter)
        {
            extra.push_str(&format!(
                "  <nuclide name=\"{daughter}\" reactions=\"0\"/>\n"
            ));
        }
    }
    out = out.replace("</depletion_chain>", &format!("{extra}</depletion_chain>"));

    Chain::from_xml(&out).expect("parse patched chain_ni.xml")
}

/// Generate a synthetic linear depletion chain of `n` nuclides where each
/// nuclide decays to the next with a small branching loss, plus a few
/// capture side-reactions.  This is large enough to show scaling but simple
/// enough to stay physically plausible.
fn synthetic_chain(n: usize) -> Chain {
    let mut nuclides: Vec<ChainNuclide> = Vec::with_capacity(n + 2);

    // Light-particle sinks that may be produced by alpha / proton decay.
    nuclides.push(ChainNuclide {
        name: "H1".into(),
        ..Default::default()
    });
    nuclides.push(ChainNuclide {
        name: "He4".into(),
        ..Default::default()
    });

    for i in 0..n {
        let name = format!("A{i}");
        let target = if i + 1 < n {
            format!("A{}", i + 1)
        } else {
            "A0".into()
        };
        let half_life = 1e3_f64 + i as f64 * 10.0;
        let mut reactions = Vec::new();
        if i % 7 == 0 {
            reactions.push(Reaction {
                kind: "(n,gamma)".into(),
                target: Some(format!("A{}", (i + 2) % n.max(1))),
                q: 0.0,
                branching_ratio: 1.0,
            });
        }
        if i % 11 == 0 {
            reactions.push(Reaction {
                kind: "(n,2n)".into(),
                target: Some(format!("A{}", i.saturating_sub(1))),
                q: 0.0,
                branching_ratio: 1.0,
            });
        }
        nuclides.push(ChainNuclide {
            name,
            half_life: Some(half_life),
            decay_energy: 1e3,
            decay_modes: vec![DecayMode {
                kind: "beta-".into(),
                target,
                branching_ratio: 0.99,
            }],
            reactions,
            neutron_fission_yields: Vec::new(),
        });
    }

    Chain::from_nuclides(nuclides).expect("synthetic chain is valid")
}

fn make_rates(chain: &Chain) -> ReactionRates {
    let mut rates = ReactionRates::new();
    for (i, nuc) in chain.nuclides.iter().enumerate() {
        let mut map = std::collections::HashMap::new();
        for r in &nuc.reactions {
            if r.kind == "(n,gamma)" {
                map.insert("(n,gamma)".into(), 1e-7);
            } else if r.kind == "(n,2n)" {
                map.insert("(n,2n)".into(), 1e-8);
            }
        }
        if !map.is_empty() {
            rates.insert(i, map);
        }
    }
    rates
}

fn bench_depletion(c: &mut Criterion) {
    let chain_ni = patched_chain_ni();
    let rates_ni = ReactionRates::new();

    // Pre-build the system so CRAM benchmarks measure only the solver.
    let sys_ni = DepletionSystem::build(chain_ni.clone(), &rates_ni).expect("build Ni system");
    let n0_ni: Vec<f64> = (0..sys_ni.chain.len())
        .map(|i| if i == 0 { 1e24 } else { 0.0 })
        .collect();
    let dt_ni = 2.592e6; // 30 days in seconds

    let mut group = c.benchmark_group("depletion_ni_chain");
    group.sample_size(30);
    group.measurement_time(std::time::Duration::from_secs(5));

    group.bench_function("build_system", |b| {
        b.iter(|| DepletionSystem::build(chain_ni.clone(), &rates_ni).expect("build"))
    });

    group.bench_function("cram16_solve", |b| {
        b.iter(|| depletion::cram(&sys_ni, Order::Order16, &n0_ni, dt_ni).expect("cram16"))
    });

    group.bench_function("cram48_solve", |b| {
        b.iter(|| depletion::cram(&sys_ni, Order::Order48, &n0_ni, dt_ni).expect("cram48"))
    });

    group.bench_function("deplete_end_to_end", |b| {
        let mut n0 = BTreeMap::new();
        n0.insert("Ni58".into(), 1e24);
        b.iter(|| depletion::deplete(&sys_ni, Order::Order48, &n0, dt_ni).expect("deplete"))
    });

    group.finish();

    // Synthetic large-chain scaling benchmarks.
    let chain_syn = synthetic_chain(300);
    let rates_syn = make_rates(&chain_syn);
    let sys_syn = DepletionSystem::build(chain_syn.clone(), &rates_syn).expect("build synthetic");
    let n0_syn: Vec<f64> = (0..sys_syn.chain.len())
        .map(|i| if i == 2 { 1e24 } else { 0.0 })
        .collect();

    let mut group_syn = c.benchmark_group("depletion_synthetic_300");
    group_syn.sample_size(20);
    group_syn.measurement_time(std::time::Duration::from_secs(8));
    group_syn.bench_function("build_system", |b| {
        b.iter(|| DepletionSystem::build(chain_syn.clone(), &rates_syn).expect("build synthetic"))
    });
    group_syn.bench_function("cram48_solve", |b| {
        b.iter(|| depletion::cram(&sys_syn, Order::Order48, &n0_syn, dt_ni).expect("cram48 syn"))
    });
    group_syn.finish();

    // Very fast name-parsing / indexing micro-benchmark with throughput.
    let mut group_micro = c.benchmark_group("depletion_indexing");
    group_micro.throughput(Throughput::Elements(sys_ni.chain.len() as u64));
    group_micro.bench_function("index_lookup_all_nuclides", |b| {
        b.iter(|| {
            for nuc in &sys_ni.chain.nuclides {
                let _ = sys_ni.chain.index_of(&nuc.name);
            }
        })
    });
    group_micro.finish();
}

criterion_group!(benches, bench_depletion);
criterion_main!(benches);
