# Roadmap

Nucleide is pre-alpha. APIs may change without notice.

## Current status

The workspace is bootstrapped with ten crates, a PyO3 mixed-layout binding, a
Python facade, and golden-byte fixtures. The canonical CI checks (format, clippy,
workspace tests, maturin build, pytest, ruff, mypy) run on every PR.

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

## Upcoming priorities

- Stabilize the Rust public API across all crates.
- Expand parser coverage and add more golden-byte fixtures.
- Add ndarray/NumPy zero-copy bridges where it improves Python ergonomics.
- Publish Rust crates to crates.io and Python wheels to PyPI.
- Add benchmarking for CRAM, cascade solving, and parser throughput.

## Out of scope

Transport solvers, Fortran discrete-ordinates ports, ENSDF evaluators,
MOAB-dependent meshing, and GUIs. Nucleide complements transport codes; it does
not replace them.
