---
title: Roadmap
sidebar:
  order: 0
---

Nucleide is pre-alpha. APIs may change without notice.

## Current status

The workspace is bootstrapped with nine crates, PyO3 and WASM bindings, a
typed Python facade, and golden-byte fixtures. The canonical CI checks (format,
clippy, workspace tests, maturin build, pytest, ruff, mypy) run on every PR.

## Recently landed

- Canonical `nucid` representation and cross-code naming dialects (`nuclei`).
- Material compositions, compendium loading, and XML export (`material`).
- MCNP-family readers for xsdir, meshtal, SSW, PTRAC, WWINP, and MCTAL
  (`mcnp-io`).
- Serpent `_res.m`, `_dep.m`, `_det.m` readers (`serpent-io`).
- FLUKA USRBIN reader and material/compound card generation (`fluka-io`).
- CRAM depletion solver and chain XML parsing (`depletion`).
- Multicomponent enrichment cascade solver (`enrichment`).
- MAGIC weight windows and mesh source sampling (`vr-tools`).
- Typed Python facade and `.pyi` stubs (`python/nucleide/`).
- Criterion benchmarks for CRAM, cascade solving, and parser throughput
  (`crates/*/benches/`).
- Cross-code validation harness against PyNE and OpenMC (`validation/`).
- Documentation website with interactive WASM tutorials (`website/`).

## Upcoming priorities

- Stabilize the Rust public API across all crates.
- Expand parser coverage and add more golden-byte fixtures.
- Add ndarray/NumPy zero-copy bridges where it improves Python ergonomics.
- Cut the first tagged release: `vX.Y.Z` tags publish Python wheels to PyPI;
  crates.io publishing is available as an opt-in release-workflow input.

## JOSS publication milestone (~6 months of public history)

A JOSS paper draft (`paper.md`, `paper.bib`) and the cross-code validation
harness (`validation/`) are in place. JOSS pre-review gates require the public
repository to show more than six months of active, iterative development and
demonstrated research use before submission, so submission waits while the
project matures in the open. Until then:

- Keep commits steady and incremental; tag real releases (0.1.0 → 0.2.0 → …)
  with matching `CHANGELOG.md` sections and PyPI wheels.
- Make research use visible (e.g. Nucleide as a dependency inside the NukeHub
  ecosystem / NukeIDE, public examples or notebooks) — this is the
  research-impact evidence the submission form asks for.
- Preserve per-version validation numbers (committed `validation/results/*.json`
  archive) so the paper can show correctness across versions.
- Prepare the mandatory JOSS AI-usage disclosure (tools/models, scope, human
  review statement) at submission time.
- Final steps at submission: draft-PDF proofread, Zenodo DOI into
  `CITATION.cff` and the paper, then the JOSS submission form.

## Out of scope

Transport solvers, Fortran discrete-ordinates ports, ENSDF evaluators,
MOAB-dependent meshing, and GUIs. Nucleide complements transport codes; it does
not replace them.
