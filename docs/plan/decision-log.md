---
title: Decision Log
sidebar:
  order: 1
---

Architecture and process decisions for Nucleide, recorded with enough rationale
that future contributors can understand the constraints.

## AD-1: Rust core with Python facade

**Decision:** Implement the core in Rust and expose it through a typed Python
API.

**Rationale:** Nuclear-engineering workflow glue spends most of its time moving
and transforming data between legacy codes. A memory-safe Rust core removes
whole classes of bugs from hand-written parsers, while Python remains the
familiar user-facing layer.

## AD-2: One crate per capability area

**Decision:** Split the workspace into focused crates (`nuclei`, `material`,
`mcnp-io`, etc.) rather than one large crate.

**Rationale:** Faster builds, clearer boundaries, and the ability to publish and
version crates independently. It also keeps `enrichment` independent of
`material` by design.

## AD-3: PyO3 mixed layout with pure-Python facade

**Decision:** Use maturin's mixed layout: `bindings/python/` builds the compiled
extension, and `python/nucleide/` provides typed stubs and re-exports.

**Rationale:** Users import a normal Python package, but the heavy lifting
happens in Rust. The facade lets us evolve the Python surface without changing
Rust symbols.

## AD-4: Workspace version as single source of truth

**Decision:** Keep the release version in `[workspace.package] version` in the
root `Cargo.toml`; all crates inherit it with `version.workspace = true`.

**Rationale:** One place to bump, no drift between crate versions, and maturin
reads the version from the workspace manifest.

## AD-5: Golden-byte fixtures as parser oracles

**Decision:** Validate parsers against byte-exact fixture files and treat parity
with legacy output as intentional.

**Rationale:** Downstream workflows depend on stable file formats. Changing
parser output requires a written reason and updated assertions.

## AD-6: `linalg` as isolation facade

**Decision:** Route linear-algebra needs through `crates/linalg` so the rest of
the workspace never depends on the numeric backend directly.

**Rationale:** Makes it possible to swap or upgrade the backend without touching
parser or depletion code.

## AD-7: License-free decay-data pack instead of vendored ICRP-107

**Decision:** Decline vendoring ICRP-107 decay data. Build the branching-ratio /
progeny / mode store, isomer masses, and free-form name normalizer from
permissive evaluations only: ENDF/B-VIII.0 decay tapes (MF8/MT457 NDK records
for branches, File-1 MT451 ELIS for isomer excitation energies) plus AME2020
ground-state masses.

**Rationale:** ICRP-107 carries redistribution terms incompatible with
vendoring into this repository, while ENDF/B-VIII.0 plus AME2020 supply the
needed branching fractions, daughter states, and excitation energies under
terms that permit redistribution. Keeping one permissive basis also keeps the
half-life, branch, and isomer-mass tables mutually consistent (the VII.1-based
decay-energy table stays as-is and records its own basis in its header).

## Open questions

- Whether to enable `abi3-py311` or stay on `abi3-py310` as the minimum Python
  version.
