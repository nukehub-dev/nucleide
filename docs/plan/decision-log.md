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

## AD-8: ENDF/B-VIII.0 half-life writer as source of truth

**Decision:** Emit `half_life.tsv` from the VIII.0 decay checkout
(`ndk_modes` half-lives, `repr` values, same stable-absent rule as the
branch table) in `run_decay8`; the legacy VII.1 `gen_decay` half-life dict
stays intentionally unused. The writer reproduces the committed data rows
byte-identically (3,821 tapes → 3,561 half-lives; verified 2026-09-10
against `zips/ENDF-B-VIII.0_decay.zip`, MD5
`aa80cd0a880d9d7e0905940b868370c3` as served). Committed values are never
silently rewritten: header-only provenance changes only.

**Rationale:** One permissive basis keeps the half-life, branch, and
isomer-mass tables mutually consistent (AD-7); byte-identical reproduction
makes the TSV a checkable artifact of its pinned upstream zip.

## AD-9: Natural-abundance source marked TO BE SOURCED

**Decision:** Record `natural_abundance.tsv` as vendored with **no**
generator path in this repo and mark the exact transcription path
**TO BE SOURCED** in the TSV header and `docs/theory/nuclear-data.mdx`
rather than claiming an evaluation basis. The only checkable statement
kept is the cross-code agreement in `validation/results.md` (289 isotopes,
max abs diff 0.0 vs OpenMC `NATURAL_ABUNDANCE`).

**Rationale:** Real provenance only: the previous "ENDF/B-VIII.0" label had
no build attribution anywhere in the repo, docs, or validation harness, so
it was removed instead of repeated.

## AD-10: CRAM coefficient provenance documented, values frozen

**Decision:** Document both CRAM sets as IPF (incomplete-partial-fraction)
coefficients — order-16 per Pusa–Leppänen (2010), order-48 per the
higher-order treatment in Pusa (2016) — in `docs/theory/depletion.mdx`,
with pinned anchors transcribed from the code (order-16/48 `alpha0`, the
order-16 k=1 pole). No published table number is claimed for order-48
until confirmed against the coefficients. Coefficient values must not
change; the `cram.rs` "Pusa & Li" comment correction and the order-48
pinned-values transcription-guard test belong to the solver owner.

**Rationale:** Checkability without touching solver code: the IPF variant
is what the code itself implements, and the two references already exist in
`paper.bib` (`pusa2010cram`, `pusa2016cram`).

## AD-11: Julian-year production convention with a documented split

**Decision:** Keep the Julian year (31,557,600 s) in production
(`depletion` inventory, `alara-io` schedule/deck) and document the two
deliberate neighbors: `fispact-io` uses a 365-day year for ALARA
file-format fidelity, and some test anchors divide by the Gregorian mean
year (gap ≈ 6e-4 yr against a 0.01 yr tolerance, e.g. the Cs-137 anchor).
Test anchors migrate to the production convention only where the current
tolerance already hides the drift; tolerances and validation gates stay
untouched.

**Rationale:** Unifying the year by code edit would move validation gates;
documenting the split with rationale keeps the physics convention explicit
at zero gate risk.

## AD-12: Reserved equation labels and literal-value cites

**Decision:** Leave (E6)–(E7) intentionally unassigned in both
`crates/depletion/src/bateman.rs` and `docs/theory/depletion.mdx` (a
code-side close-the-gap renumber belongs to the depletion owner; no
formula changes either way). Literal cites live at their use sites:
`THERMAL_EV`/`FAST_EV` in `scripts/gen-nuclear-data.py` cite the NIST
2200 m/s page / ENDF 14-MeV MF3/MT1 reference, and the enrichment cascade
tuning constants already carry empirical-rationale comments in
`crates/enrichment/src/cascade.rs` (any relabeling there belongs to the
enrichment owner). Every TSV `Regenerate` line names its exact flags.

**Rationale:** Numbering and constants stay consistent across code and docs
without cross-owner edits; explicit regen flags make every table
reproducible from its pinned inputs.

## Open questions

- Whether to enable `abi3-py311` or stay on `abi3-py310` as the minimum Python
  version.
