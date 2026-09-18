# Nucleide

Nucleide is a modern Rust toolkit for nuclear-engineering data, measurement, and
 workflow glue: legacy transport-code I/O, nuclide identification, materials,
 CRAM depletion, enrichment analytics, point kinetics, gamma-ray spectroscopy,
 tritium transport, TBR/blanket bookkeeping, equilibrium data readers,
 an OpenMC statepoint tally bridge, fusion neutron sources, spectrum
 unfolding, variance
 reduction, UQ sampling, CSG translation, and code-card emission — exposed
 through a typed Python API.

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
| Nuclide core (`nucleide-nuclei`) | Canonical nucid representation, particle registry, reaction-name registry (labels, MT mapping, hashes), name-dialect conversions (ZZAAAMM, ZAID/MCNP, Serpent, FLUKA, NIST, CINDER, ALARA, SZA, ARMI/MCC3), AME2020 masses (incl. isomer masses), natural abundances, half-lives, screening cross sections / scattering lengths / prompt decay energies (generated from ENDF/B + NIST), ENDF/B-VIII.0 decay branches, free-form name normalization, dose factors, EPA FGR 15 external-dosimetry coefficients (runtime hash-pinned download), IRDFF-II dosimetry-response pack (34 named reactions; runtime hash-pinned download, parsed into unfolding response rows) |
| Materials (`nucleide-material`) | Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium loading, materials XML export, activity/decay-heat/dose-per-gram analytics, label-collision checks and conservation audits, mass-efficiency separator / fixed-ratio blender, Page CUSUM change detector |
| MCNP I/O (`nucleide-mcnp-io`) | xsdir, meshtal, SSW/SURFSRC, PTRAC, WWINP, MCTAL (headers, kcode, standard tally bodies), ENDL readers; typed legacy SDEF fixed-source reader (round-trips the spectroscopy emitter byte-identically); NumPy `result_array()` / `totals_array()` meshtal and `tally_vals_array()` MCTAL bridges; material extraction from input decks; full-deck parse/edit/write round-trip (cells, surfaces, materials); L3 semantic views (MODE/TRn/universes/lattices/FILL/tallies) with validation; mesh-to-geometry deck generation |
| CSG translation (`nucleide-csg-xlate`) | Scoped MCNP CSG translation to OpenMC `geometry.xml`, Serpent, PHITS, and GDML (Geant4) output (surfaces, cells, universes/fills, rectangular `LAT=1` lattices, material stub) with a drift report; macrobodies expand, unsupported constructs fail with a clear error |
| MCPL I/O (`nucleide-mcpl-io`) | Monte Carlo Particle List interchange reader/writer (format versions 2/3, single/double precision, gzip-transparent) plus SSW↔MCPL conversion and the merge/extract/stats/repair particle-list utilities |
| Serpent I/O (`nucleide-serpent-io`) | `_res.m`, `_dep.m`, `_det.m` readers producing structured records |
| FLUKA I/O (`nucleide-fluka-io`) | USRBIN tally reader, material/compound card generation |
| ALARA I/O (`nucleide-alara-io`) | Deck/flux/matlib-elelib-WDR/output/photon/schedule-expansion glue; solver out of scope. Clearance / waste-classification analytics: clearance index CI = Σ Aᵢ/CLᵢ and the sum-of-fractions screening rule over parsed inventories, caller-supplied limit tables plus the EU 2013/59/Euratom Annex VII Table A vendored default and the Spanish CSN conditional NORM landfill tables as opt-in tables (screening arithmetic, never a compliance decision). Sublet S1–S7 totals: total activity with the IRT α/β/γ split, decay heat over caller-supplied decay energies (each with the excluding-tritium companion), committed ingestion/inhalation hazards over caller 50-year dose coefficients (each with the companion), the transport Bq/A₂ ratio with its effective A₂, the IAEA clearance-index variant screened at ≤ 1, and the slab/point gamma dose over caller groups with the loud 0.3 m clamp |
| Depletion (`nucleide-depletion`) | CRAM (orders 16/48) matrix exponential, analytic Bateman fast path (`method=` selector with CRAM-48 fallback), depletion-chain XML parsing, Predictor/CECM/CF4 time-series integrators with activity/decay-heat observables, unit-aware decay inventories, cumulative decays and chain-lineage queries |
| Enrichment (`nucleide-enrichment`) | Multicomponent cascade solver (numeric), SWU closed-form helpers |
| Point kinetics (`nucleide-kinetics`) | Prescribed-reactivity PKE solver, inhour roots, prompt-jump factor |
| Tritium transport (`nucleide-tritium`) | 1D Fick + McNabb–Foster diffusion-trapping kernel, Dirichlet/Sieverts/Henry/zero-flux/recombination surfaces (steady state and transient), permeation breakthrough and time lag, multi-layer series stacks with Sieverts, Henry, or vented-sink recombination internal interfaces |
| TBR/blanket books (`nucleide-blanket`) | Raw TBR from caller tallies, multiplicative per-port coverage haircuts, blanket energy multiplication, tritium burn rate, and breeding-margin / net-surplus fuel-cycle metrics (pure arithmetic — no transport, no geometry optimization, no coupling into `tritium`) |
| Fusion sources (`nucleide-plasma-source`) | Tokamak neutron sources — ring/point (D-D 2.45 MeV, D-T 14.1 MeV) and a parametric Miller-geometry plasma with caller-supplied L/H/A-mode profiles (Fausser 2012), reactivity-weighted emission (Bosch–Hale 1992) over equimolar D-T, pure D-D, or D/T fuel mixtures at a shared ion temperature (Eriksson 2016), per-species ion temperatures, one pinned deuterium hot-tail shape, toroidal sectors, an arbitrary-3D birth-rate lattice with field-period symmetry, and closed-form ECRH accessibility (cold resonance, relativistic shift, O1/X1 cut-offs, beamline crossings), ion-temperature-broadened Gaussian spectra (Brysk/Ballabio), seeded sampling to particle vectors, MCNP SDEF + Serpent source-card emission with drift report |
| Damage metrics (`nucleide-damage`) | NRT-dpa and arc-dpa displacement functions (NRT 1975; Nordlund 2018), He/H appm production and He/dpa ratios by spectral folding of caller flux with caller response functions, UQ on the folds over caller MVN blocks (SPECTER is the validation oracle behind an opt-in vendored Table VII fallback), plus coil fast-fluence / lifetime bookkeeping (fast-flux sums, history accumulation, weakest-link life over caller limit tables) |
| Spectroscopy (`nucleide-spectroscopy`) | Spectrum smoothing, gross/net counting, energy/efficiency calibration, X-ray lines, SPE parsing, decay-line SDEF source cards (E9) fed from caller lists or the runtime decay-lines TSV interchange |
| Spectrum unfolding (`nucleide-unfold`) | SAND-II, STAYSL-class least-squares, GRAVEL chi-square-weighted, and MAXED maximum-entropy adjustment of a guess neutron spectrum against measured activation rates (caller-supplied response matrix, or the IRDFF-II runtime pack), with convergence diagnostics; non-convergence is a hard error |
| Equilibrium data (`nucleide-equilib-io`) | Classic-netCDF VMEC `wout` reader (CDF-1/CDF-2 only; netCDF-4/CDF-5 convert facade-side) plus the `&INDATA` input-text grammar, with flux-surface Jacobian helpers for volume weighting and wall-load mapping in flux coordinates (reads data, never solves equilibria) |
| OpenMC tally bridge (`nucleide.openmc`, pure Python) | Statepoint mesh/cell tally extraction through the caller-side OpenMC API into plain arrays the damage folds accept, plus a CSV interchange for the same dict |
| Variance reduction (`nucleide-vr-tools`) | MAGIC weight-window generation, OpenMC/Serpent weight-window emission, mesh source sampling with alias tables |
| UQ sampling (`nucleide-linalg`) | Seeded MVN + log-normal + LHS draws over caller-supplied covariance blocks, SANDY-compatible estimators, decay-data and fission-yield perturbation consumers |
| CCCC I/O (`nucleide-cccc-io`) | ISOTXS/RTFLUX text-subset parsers + PARTISN deck writer (no solver) |
| FISPACT I/O (`nucleide-fispact-io`) | FISPACT-II inventory output parser reusing the ALARA response frame, plus the clearance-bearing wide inventory table printed with the `HAZARDS` + `CLEAR` keywords (output-only) |
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
│   ├── csg-xlate/     # scoped MCNP CSG -> OpenMC/Serpent/PHITS/GDML translation
│   ├── mcpl-io/       # MCPL interchange read/write + SSW conversion + merge/extract/stats/repair
│   ├── serpent-io/    # res/dep/det readers
│   ├── fluka-io/      # usrbin reader, material cards
│   ├── alara-io/      # ALARA deck/flux/libs/output/photon/schedule glue + clearance + Sublet S1–S7 totals (no solver)
│   ├── cccc-io/       # ISOTXS/RTFLUX text-subset parsers + PARTISN writer (no solver)
│   ├── fispact-io/    # FISPACT-II inventory output + CLEAR-keyword clearance table (output-only)
│   ├── origen-io/     # scoped ORIGEN TAPE5/6/9 readers
│   ├── r2s/           # scoped R2S workflow builder (photon tags + spectra)
│   ├── vr-tools/      # MAGIC weight windows + OpenMC/Serpent emission, source sampling
│   ├── enrichment/    # cascades, SWU
│   ├── depletion/     # CRAM + chain files
│   ├── kinetics/      # prescribed-reactivity point kinetics + inhour
│   ├── tritium/       # 1D diffusion-trapping kernel + single/multi-layer permeation checks
│   ├── blanket/       # TBR + port penalties + power multiplication + burn/fuel-cycle margin
│   ├── plasma-source/ # tokamak ring/point/parametric-plasma fusion sources + SDEF/Serpent cards
│   ├── damage/        # NRT/arc-dpa + He/H appm + He/dpa spectral folds + UQ
│   ├── unfold/        # neutron spectrum unfolding: SAND-II, STAYSL-class, GRAVEL, MAXED
│   ├── equilib-io/    # classic-netCDF wout reader + INDATA grammar + Jacobian helpers (no solver)
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
