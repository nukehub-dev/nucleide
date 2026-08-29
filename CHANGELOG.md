# Changelog

All notable changes to Nucleide are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Releases are cut with `scripts/bump-version.sh X.Y.Z`, which stamps the
`[Unreleased]` section below and updates `[workspace.package] version` in the
root `Cargo.toml`. Git tags (`vX.Y.Z`) are the release source of truth; CI
publishes Python wheels (and optionally Rust crates) from tags.

## [Unreleased]

### Added

- `depletion`: reaction `branching_ratio` support in chain XML (default 1.0),
  with branched gains in the burnup matrix.
- `depletion`: light-particle production in the burnup matrix — He-4/H-1 from
  alpha/proton decays and from reaction secondaries such as (n,α) and (n,p),
  when the product nuclide is in the chain.
- `depletion`: decay branching ratios are renormalized to sum to 1 (largest
  branch adjusted), matching OpenMC.
- `nuclei`: `NuclideId::from_name` now accepts PyNE-normalized forms such as
  `"U-235"`, `"u235"`, and uppercase-`M` metastable markers.
- Python API: `Cascade.solve_multicomponent()` (M\*-optimizing solve),
  `DepletionSystem` builder and `solve_vec()` for repeated solves without
  per-call overhead, and `MeshTally.total_rel_error`.
- Criterion benchmarks in `crates/*/benches/` (`cargo bench`).
- `validation/`: runnable cross-code validation harness against PyNE 0.7.5 and
  OpenMC 0.16.0 (containerized via `validation/Containerfile` +
  `validation/run_container.sh`) with committed results
  (`validation/results.md`).
- JOSS submission materials: `paper.md`, `paper.bib`, `CITATION.cff`,
  `CONTRIBUTING.md`, and a draft-PDF workflow.

### Changed

- Python API reorganized into domain submodules (`nucleide.nuclei`,
  `nucleide.material`, `nucleide.mcnp`, `nucleide.serpent`, `nucleide.fluka`,
  `nucleide.vr`, `nucleide.enrichment`, `nucleide.depletion`) mirroring the
  workspace crates; the flat `nucleide.*` namespace is reduced to `version()`
  and `__version__` before the first release.
- `vr-tools`: `MeshSourceSampler` now rejects negative or non-finite tally /
  user-density values with an error instead of silently absolutizing them.
- `depletion`: reaction loss is subtracted once per reaction type per nuclide;
  duplicate entries only add branched gains.
- `depletion`: CRAM pole solves reuse scratch buffers across poles.
- `nuclei`: FLUKA name lookups use lazily-built maps instead of linear scans.

### Fixed

- `nuclei`: `NuclideId::new` rejects mass numbers above 999, fixing a debug
  overflow panic from inputs like `"U999999"`.
- `depletion`: zero/negative/non-finite half-lives and non-positive `dt` now
  return errors instead of producing inf/NaN.
- `linalg`: LU solve validates matrix/RHS dimensions and symbolic-pattern
  identity instead of panicking inside faer.
- `enrichment`: validates `alpha > 1`, assay ordering, and `M*` bounds; guards
  secant divide-by-zero; golden-section M* polish no longer trusts unconverged
  probe solves; `recompute_nm` reset path no longer uses stale right-hand sides.
- `vr-tools`: `magic` no longer panics on tallies with empty energy bounds and
  validates array lengths and finiteness.
- `bindings/wasm`: fraction and result maps are now serialized as plain JS
  objects instead of JS `Map`s, fixing empty tables in all interactive
  tutorials.

### Removed

- `vr-tools`: unused `nuclei` dependency.
