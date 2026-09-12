# Nucleide

Nucleide is a modern Rust toolkit for nuclear-engineering data, measurement, and
 workflow glue: legacy transport-code I/O, nuclide identification, materials,
 CRAM depletion, enrichment analytics, point kinetics, gamma-ray spectroscopy,
 and code-card emission — exposed through a typed Python API.

The project is a fresh Rust implementation of capabilities pioneered by
[PyNE](https://github.com/pyne/pyne), focused on memory safety, fast builds,
and `pip install`-able wheels. Scope is intentionally narrow today and will
expand as more parsers and workflow pieces land.

## Installation

```bash
pip install nucleide
```

Prebuilt wheels cover Linux, macOS, and Windows for Python >= 3.10 (abi3). To
build from source instead, see the Development section below.

## Why

Nuclear-engineering workflows spend most of their time moving data between
codes rather than solving transport itself. The established tooling for that
glue layer carries a heavy build chain (CMake + Fortran + Cython) and
hand-written parsers that are hard to extend and harder to embed. `Nucleide`
rebuilds the high-value subset in memory-safe Rust with one-command
`pip install` wheels, keeping Python as the user-facing API.

## Features

| Area | Capabilities |
| --- | --- |
| Nuclide core (`nucleide-nuclei`) | Canonical nucid representation, particle registry, reaction-name registry (labels, MT mapping, hashes), name-dialect conversions (ZZAAAMM, ZAID/MCNP, Serpent, FLUKA, NIST, CINDER, ALARA, SZA, ARMI/MCC3), AME2020 masses (incl. isomer masses), natural abundances, half-lives, screening cross sections / scattering lengths / prompt decay energies (generated from ENDF/B + NIST), ENDF/B-VIII.0 decay branches, free-form name normalization, dose factors |
| Materials (`nucleide-material`) | Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium loading, materials XML export, activity/decay-heat/dose-per-gram analytics, label-collision checks and conservation audits, mass-efficiency separator / fixed-ratio blender, Page CUSUM change detector |
| MCNP I/O (`nucleide-mcnp-io`) | xsdir, meshtal, SSW/SURFSRC, PTRAC, WWINP, MCTAL (headers, kcode, standard tally bodies), ENDL readers; NumPy `result_array()` / `totals_array()` meshtal and `tally_vals_array()` MCTAL bridges; material extraction from input decks; full-deck parse/edit/write round-trip (cells, surfaces, materials); L3 semantic views (MODE/TRn/universes/lattices/FILL/tallies) with validation; mesh-to-geometry deck generation |
| MCPL I/O (`nucleide-mcpl-io`) | Monte Carlo Particle List interchange reader/writer (format versions 2/3, single/double precision, gzip-transparent) plus neutron/gamma-only SSW↔MCPL conversion |
| Serpent I/O (`nucleide-serpent-io`) | `_res.m`, `_dep.m`, `_det.m` readers producing structured records |
| FLUKA I/O (`nucleide-fluka-io`) | USRBIN tally reader, material/compound card generation |
| ALARA I/O (`nucleide-alara-io`) | Deck/flux/matlib-elelib-WDR/output/photon/schedule-expansion glue; solver out of scope |
| Depletion (`nucleide-depletion`) | CRAM (orders 16/48) matrix exponential, analytic Bateman fast path (`method=` selector with CRAM-48 fallback), depletion-chain XML parsing, Predictor/CECM/CF4 time-series integrators with activity/decay-heat observables, unit-aware decay inventories, cumulative decays and chain-lineage queries |
| Enrichment (`nucleide-enrichment`) | Multicomponent cascade solver (numeric), SWU closed-form helpers |
| Point kinetics (`nucleide-kinetics`) | Prescribed-reactivity PKE solver, inhour roots, prompt-jump factor |
| Spectroscopy (`nucleide-spectroscopy`) | Spectrum smoothing, gross/net counting, energy/efficiency calibration, X-ray lines, SPE parsing, decay-line SDEF source cards (E9) fed from caller lists or the runtime decay-lines TSV interchange |
| Variance reduction (`nucleide-vr-tools`) | MAGIC weight-window generation, mesh source sampling with alias tables |
| CCCC I/O (`nucleide-cccc-io`) | ISOTXS/RTFLUX text-subset parsers + PARTISN deck writer (no solver) |
| FISPACT I/O (`nucleide-fispact-io`) | FISPACT-II inventory output parser reusing the ALARA response frame (output-only) |
| ORIGEN I/O (`nucleide-origen-io`) | Scoped ORIGEN 2.2 TAPE5 input-echo, TAPE6 inventory, and TAPE9 decay readers |
| R2S (`nucleide-r2s`) | Scoped R2S workflow builder: zone-to-flux linking, schedule expansion, photon assembly (uniform-split placeholder plus per-voxel source tags and `.photonSrc` group spectra), ARMI database-snapshot adapter |
| Emit (`nucleide-emit`) | Single-material emission to MCNP/Serpent/FLUKA/ALARA/PARTISN cards with mass-drift report, plus an ARMI blueprint-key bridge |
| Python bindings | PyO3 extension module behind a typed pure-Python facade (`nucleide._internal`, `.pyi` stubs, `py.typed`) |

## Out of scope

Transport solvers, Fortran discrete-ordinates ports, ENSDF evaluators, MOAB-dependent
meshing, and GUIs. `Nucleide` complements transport codes; it does not replace them.

## Layout

```text
nucleide/
├── crates/
│   ├── nuclei/        # nuclide ids, naming conventions, physical data
│   ├── material/      # compositions, mixing, libraries, XML export
│   ├── mcnp-io/       # xsdir/meshtal/SSW/MCTAL/PTRAC/WWINP
│   ├── mcpl-io/       # MCPL interchange read/write + SSW conversion
│   ├── serpent-io/    # res/dep/det readers
│   ├── fluka-io/      # usrbin reader, material cards
│   ├── alara-io/      # ALARA deck/flux/libs/output/photon/schedule glue (no solver)
│   ├── cccc-io/       # ISOTXS/RTFLUX text-subset parsers + PARTISN writer (no solver)
│   ├── fispact-io/    # FISPACT-II inventory output parser (output-only)
│   ├── origen-io/     # scoped ORIGEN TAPE5/6/9 readers
│   ├── r2s/           # scoped R2S workflow builder (photon tags + spectra)
│   ├── vr-tools/      # MAGIC weight windows, source sampling
│   ├── enrichment/    # cascades, SWU
│   ├── depletion/     # CRAM + chain files
│   ├── kinetics/      # prescribed-reactivity point kinetics + inhour
│   ├── spectroscopy/  # smoothing, counting, calibration, X-ray, SPE
│   ├── emit/          # five-dialect card emission + mass-drift reports
│   └── linalg/        # isolation facade over the linear-algebra backend
├── bindings/python/   # PyO3 crate -> nucleide._internal
├── python/nucleide/   # typed pure-Python facade (maturin mixed layout)
├── fixtures/          # golden-byte test data
├── validation/        # cross-code validation harness vs PyNE/OpenMC
└── tests/             # Python-side tests
```

## Development

```bash
git clone https://github.com/nukehub-dev/nucleide.git
cd nucleide

# Rust side
cargo test                       # workspace unit tests
cargo clippy --all-targets -- -D warnings

# Python side (needs: rustup, pip install maturin)
pip install maturin pytest pytest-cov ruff mypy
maturin develop                  # build + install into current venv
pytest tests/

# Lint / type-check / format the Python surface
ruff format python tests
ruff check python tests
mypy                             # strict type-check against .pyi stubs

# Rust coverage (needs llvm-tools-preview component)
cargo llvm-cov --workspace       # or --lcov for CI upload
```

## Tooling

| Layer | Format | Lint | Types | Coverage |
| --- | --- | --- | --- | --- |
| Rust | rustfmt (`cargo fmt`) | clippy `-D warnings` | — | cargo-llvm-cov (CI) |
| Python | ruff format | ruff check | mypy `--strict` via `.pyi` stubs | pytest-cov |

Wheels are built with **maturin** (PyO3 mixed layout). One wheel serves all
Python >= 3.10 via abi3 — the same stack used by pydantic-core, polars, and ruff.

## Validation strategy

1. Parsers are validated against **golden-byte fixtures** in `fixtures/`;
   parser output must match recorded snapshots before any release.
2. Numeric kernels (CRAM, cascade solving) are checked against published
   analytic vectors and cross-code results on shared inputs; the runnable
   cross-code harness in [`validation/`](https://github.com/nukehub-dev/nucleide/tree/main/validation) compares
   Nucleide against PyNE and OpenMC and commits its measured results.
3. Behavioral compatibility with legacy tool output is asserted wherever a
   fixture exists, so downstream workflows see identical data.

Criterion benchmarks (`cargo bench`) cover the numeric kernels and parsers.

## Citing

If you use Nucleide in research, see
[`CITATION.cff`](https://github.com/nukehub-dev/nucleide/blob/main/CITATION.cff)
and the JOSS paper draft in
[`paper.md`](https://github.com/nukehub-dev/nucleide/blob/main/paper.md).

## Status

Pre-alpha. APIs may change without notice.

## Documentation

Additional tutorials, reference pages, and developer guides live in the
[`docs/`](https://github.com/nukehub-dev/nucleide/tree/main/docs) tree.

## Acknowledgments

Nucleide is a fresh Rust implementation of workflow-glue capabilities pioneered
by [PyNE](https://github.com/pyne/pyne) ("Python for Nuclear Engineering",
BSD-3-Clause). Some reference data and golden test fixtures — notably the
DOE/PNNL Materials Compendium — are vendored directly from PyNE; see
`fixtures/data/MaterialsCompendium.LICENSE` for its terms.

## License

[BSD-2-Clause](https://github.com/nukehub-dev/nucleide/blob/main/LICENSE).
