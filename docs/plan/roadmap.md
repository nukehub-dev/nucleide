---
title: Roadmap
sidebar:
  order: 0
---

Nucleide is pre-alpha. APIs may change without notice.

## Current status

The workspace is bootstrapped with fourteen crates, PyO3 and WASM bindings, a
typed Python facade, and golden-byte fixtures. The canonical CI checks (format,
clippy, workspace tests, maturin build, pytest, ruff, mypy) run on every PR.

## Recently landed

- Canonical `nucid` representation and cross-code naming dialects (`nuclei`).
- Material compositions, compendium loading, and XML export (`material`).
- MCNP-family readers for xsdir, meshtal, SSW, PTRAC, WWINP, and MCTAL
  (`mcnp-io`).
- Serpent `_res.m`, `_dep.m`, `_det.m` readers (`serpent-io`).
- FLUKA USRBIN reader and material/compound card generation (`fluka-io`).
- ALARA Phase 1 interop (`alara-io`): deck/flux/libs/output/photon/schedule
  glue (solver out of scope).
- CRAM depletion solver and chain XML parsing (`depletion`).
- Multicomponent enrichment cascade solver (`enrichment`).
- MAGIC weight windows and mesh source sampling (`vr-tools`).
- CCCC text-subset parsers + PARTISN writer (`cccc-io`): ISOTXS/RTFLUX
  readers and deck validation (no solver).
- FISPACT-II inventory output parser (`fispact-io`): output-only, reusing the
  ALARA response frame.
- Scoped ORIGEN 2.2 TAPE5/6/9 readers (`origen-io`).
- R2S orchestrator landed as a scoped workflow builder (`r2s`): zone-to-flux
  linking, schedule expansion, and uniform-split photon assembly (transport
  and activation solving stay out of scope).
- Typed Python facade and `.pyi` stubs (`python/nucleide/`).
- Criterion benchmarks for CRAM, cascade solving, and parser throughput
  (`crates/*/benches/`).
- Cross-code validation harness against PyNE and OpenMC (`validation/`).
- Activation-code I/O validation (`validation/activation_vs_refs.py`):
  ALARA/CCCC/FISPACT/ORIGEN/R2S checks on committed fixtures vs PyNE oracle
  probes + synthetic self-consistency.
- Documentation website with interactive WASM tutorials (`website/`).

## Upcoming priorities

- Stabilize the Rust public API across all crates.
- Add ndarray/NumPy zero-copy bridges where it improves Python ergonomics.
- Cut releases: `vX.Y.Z` tags publish Python wheels to PyPI and workspace
  crates to crates.io.

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
