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
  branch adjusted) for programmatically built chains, matching OpenMC's
  chain-generation behavior.
- `depletion`: CASL-style fission-yield borrowing in chain XML
  (`<neutron_fission_yields parent="X"/>` resolves X's yields, transitively),
  matching OpenMC.
- `nuclei`: `NuclideId::from_name` now accepts PyNE-normalized forms such as
  `"U-235"`, `"u235"`, and uppercase-`M` metastable markers.
- Python API: `Cascade.solve_multicomponent()` (M\*-optimizing solve),
  `DepletionSystem` builder and `solve_vec()` for repeated solves without
  per-call overhead, and `MeshTally.total_rel_error`.
- Python API: `nucleide.data` module (`fetch`, `fetch_compendium`,
  `default_ref`) for downloading repo data files — the Materials Compendium,
  sample depletion chains — pinned to the installed version's tag, since the
  wheel bundles no data files.
- Criterion benchmarks in `crates/*/benches/` (`cargo bench`).
- `validation/`: runnable cross-code validation harness against PyNE 0.7.5 and
  OpenMC 0.16.0 (containerized via `validation/Containerfile` +
  `validation/run_container.sh`) with committed results
  (`validation/results.md`).
- `validation/`: full-chain depletion validation on the CASL/VERA simplified
  chain (downloaded to the git-ignored `validation/.cache/`), parser
  cross-validation against PyNE and serpentTools oracles
  (`parsers_vs_refs.py`), and generated paper figures
  (`validation/figures/`, via `make_figures.py`).
- JOSS submission materials: `paper.md`, `paper.bib`, `CITATION.cff`,
  `CONTRIBUTING.md`, and a draft-PDF workflow.
- `bindings/wasm`: `wasm-bindgen` crate exposing a subset of the workspace to
  the browser, including a `WasmMaterialsCompendium` API that parses the
  DOE/PNNL Materials Compendium from its JSON text.
- Documentation website (`website/`, Astro + `@nukehub/docs-kit`) with content
  synced from `docs/`, theory pages, and interactive WASM tutorials (including
  a compendium browser that fetches the staged `MaterialsCompendium.json` and
  charts compositions with Plotly), deployed to GitHub Pages.
- Release workflow: `vX.Y.Z` tags build and publish Python wheels for Linux,
  macOS, and Windows to PyPI and draft a GitHub release from the matching
  changelog section; crates.io publishing is opt-in via `workflow_dispatch`.

### Changed

- Documentation website moved from `https://nukehub-dev.github.io/nucleide` to
  a custom domain, `https://nucleide.nukehub.org`, hosted on GitHub Pages with
  Cloudflare DNS; the Astro `base` path changed from `/nucleide` to `/`.

- Documentation theory pages (`depletion`, `enrichment`, `nuclear-data`,
  `variance-reduction`) now use the citation support from `@nukehub/docs-kit`:
  references are declared in frontmatter, cited inline with `<Citation />`, and
  rendered as a linked bibliography with copy-to-clipboard export.

- Python API reorganized into domain submodules (`nucleide.nuclei`,
  `nucleide.material`, `nucleide.mcnp`, `nucleide.serpent`, `nucleide.fluka`,
  `nucleide.vr`, `nucleide.enrichment`, `nucleide.depletion`) mirroring the
  workspace crates; the top level re-exports the domain submodules alongside
  `version()` and `__version__`.
- `vr-tools`: `MeshSourceSampler` now rejects negative or non-finite tally /
  user-density values with an error instead of silently absolutizing them.
- `depletion`: reaction loss is subtracted once per reaction type per nuclide;
  duplicate entries only add branched gains.
- `depletion`: CRAM pole solves reuse scratch buffers across poles.
- `nuclei`: FLUKA name lookups use lazily-built maps instead of linear scans.

### Fixed

- Documentation website deployment now builds the `bindings/wasm` package in
  CI, so the interactive tutorials can load `/wasm/nucleide_wasm.js` instead of
  receiving the SPA fallback HTML response.

- Interactive tutorials: file/text inputs now use the docs-kit `Textarea`
  auto-resize (content-fitted height, no manual resize grip); kit bumped for a
  `Select` dropdown fix so short lists flipped above the trigger no longer
  float with a gap.
- `depletion`: `Chain::from_xml` no longer renormalizes decay branching
  ratios — file values are used verbatim, matching OpenMC's `Chain.from_xml`
  (renormalization only happens at chain *generation*). Found by the CASL
  full-chain validation (e.g. I-128 β⁻ branching ratio 0.931).
- `depletion`: fission production now uses the yield set at the lowest
  incident neutron energy, matching OpenMC's `get_default_fission_yields`.

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
