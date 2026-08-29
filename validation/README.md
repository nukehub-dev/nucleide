# Nucleide cross-code validation harness

This directory contains runnable scripts that compare Nucleide against the
reference codes PyNE 0.7.5 and OpenMC 0.16.0 for a JOSS publication. The
committed `results.md` holds the measured numbers that the paper quotes.

## Files

| File | Purpose |
| --- | --- |
| `depletion_vs_openmc.py` | CRAM-48 depletion on `fixtures/depletion/chain_ni.xml` + 3-nuclide analytic Bateman check |
| `depletion_casl_vs_openmc.py` | CRAM-48 depletion on the full CASL/VERA chain (fission yields, branching) vs OpenMC |
| `enrichment_vs_pyne.py` | Uranium and tungsten enrichment cascades vs PyNE `multicomponent` |
| `magic_vs_pyne.py` | MAGIC weight-window generation vs PyNE (or equivalent formula when PyMOAB is missing) |
| `nuclear_data_vs_refs.py` | Atomic masses, abundances, half-lives, name-dialect conversions vs PyNE/OpenMC |
| `parsers_vs_refs.py` | Serpent/MCNP/FLUKA parser cross-validation vs serpentTools and PyNE oracles |
| `timings.py` | Coarse wall-time comparisons (Python overhead included) |
| `make_figures.py` | Generates the paper figures in `figures/` from the JSON reports |
| `common.py` | Shared helpers and `Report` class used by the comparison scripts |
| `render_results.py` | Renders `results.md` from the JSON reports in `results/` |
| `environment.json` | Hand-maintained source annotations for the environment table |
| `run_all.sh` | Run the whole harness with a configurable Python binary |
| `run_container.sh` | One-command container run: rebuilds the image (layer-cached), then the wheel, then the harness |
| `Containerfile` | Validation environment definition (PyNE 0.7.5 + OpenMC 0.16.0 + serpentTools 0.11.0, Python 3.12) |
| `magic_tally.txt` | Small synthetic MCNP meshtal used by the MAGIC comparison |
| `results/` | Per-run JSON reports produced by each comparison script |
| `results.md` | Committed measured results (generated; do not hand-edit) |
| `figures/` | Generated paper figures (`timings.png`, `depletion_agreement.png`); committed, never hand-edited |
| `.cache/` | Download cache for the CASL chain (git-ignored; see below) |

## Environment

The canonical environment is the container built from `Containerfile`:

- Python 3.12 on `condaforge/miniforge3`
- PyNE 0.7.5 (`nomoab_noopenmc` build from conda-forge, so `pyne.mcnp.Meshtal`
  is unavailable). Note `pyne.__version__` reports `0.7.1` — upstream's Python
  metadata lags the actual release; the conda package version (0.7.5) is the
  truth. 0.7.5 is numerically identical to PyNE 0.7.8 for the modules used
  here — the 0.7.5..0.7.8 diff in `src/enrichment.cpp`, `src/data.cpp` and
  `pyne/dbgen/atomic_mass.py` is quote style and `std::isnan` qualification
  only. PyNE 0.7.8 itself cannot build on Python 3.12 (its `setup.py` still
  does `import imp`, removed in 3.12), and conda-forge tops out at 0.7.5.
- OpenMC 0.16.0 built from the upstream release tag. conda-forge still ships
  0.15.3, and the `docker.io/openmc/openmc:v0.16.0` image actually contains
  0.15.3 (its Dockerfile defaults `openmc_branch=master` and the published tag
  was built from the 0.15.3 release commit), so neither binary channel carries
  0.16.0 yet.
- Nucleide abi3 wheel built in release mode by the `ghcr.io/pyo3/maturin`
  container (manylinux2014; release mode is required for meaningful timing
  numbers).
- serpentTools 0.11.0 from PyPI, as the independent oracle for the Serpent
  parser comparison (conda-forge's `serpent-tools` build lags).

A conda env with `pyne` and `openmc` from conda-forge also works (see
`run_all.sh`), but tracks the older conda-forge versions.

## Running

One command (requires `podman`; the first image build takes ~20 min, later
builds reuse the layer cache and are near-instant unless `Containerfile`
changed):

```bash
./validation/run_container.sh
```

Or against any existing Python that has PyNE, OpenMC and the Nucleide release
wheel installed:

```bash
./validation/run_all.sh /path/to/python
```

Each script exits nonzero only on an unexpected runtime failure; reported
quantitative differences are printed, not treated as failures.

## Adding a new comparison

1. Copy an existing script as a template.
2. Import `common.Report` and `common.fmt` and use `report.heading()`,
   `report.prose()`, and `report.table()` to capture every section you want in
   `results.md`.
3. Call `report.emit()` at the end of `main()` to write
   `validation/results/<name>.json` and print a plain-text summary.
4. Name the file `*_vs_*.py` — `run_all.sh` auto-discovers every such script
   and runs it (sorted) before `timings.py`, `render_results.py`, and
   `make_figures.py`.
5. Add the new report name to `render_results.py`'s `SECTION_ORDER`.
6. Regenerate `results.md` and the figures by running
   `./validation/run_container.sh` or `./validation/run_all.sh [python]`.

## Download cache and provenance

`depletion_casl_vs_openmc.py` downloads the CASL simplified depletion chain
(thermal-spectrum branching, 228 nuclides; VERA Depletion Benchmark chain from
CASL-U-2015-1014-000) from the official OpenMC pregenerated-chain archive at
`https://anl.box.com/shared/static/3nvnasacm2b56716oh5hyndxdyauh5gs.xml`. It is
cached at `validation/.cache/chain_casl_thermal.xml` on first use and
re-downloaded only if absent. `.cache/` is git-ignored — never commit the
chain file. Everything else runs on committed inputs from `fixtures/` or this
directory.

## Parser oracle notes

`parsers_vs_refs.py` cross-checks Nucleide's readers against independent
oracles on the committed fixtures:

- **Serpent** — serpentTools 0.11.0 (pip-pinned in `Containerfile`). Compared
  on `sample2_dep.m` (dep) and `serp2_det.m` (det): `ZAI`, `DAYS`, `BU`,
  `NAMES`, per-material `ADENS`/`MDENS`/`VOLUME`, detector tally and
  relative-error columns, and shared bin grids. serpentTools cannot read the
  Serpent-1-format fixtures (`sample_res.m`, `sample1_dep.m`, `sample_det.m`)
  and its `ResultsReader` postprocessing crashes on both `*_res.m` files, so
  those comparisons skip loudly.
- **MCNP** — PyNE `pyne.mcnp`. Working oracles: `Xsdir` (table entries),
  `SurfSrc` (headers and track payloads on all four `.w` fixtures), and
  `PtracReader` (problem title and variable counts on three of four `.ptrac`
  files; PyNE has no `PtracFile` class). `Wwinp` and `Meshtal` require PyMOAB
  (absent in the nomoab build) and skip loudly.
- **FLUKA** — no working oracle: `pyne.fluka.Usrbin` needs PyMOAB and reads
  only binary USRBIN, while our fixtures are ASCII `.lis`. Skipped with an
  explicit note.

## Known limitations

- `chain_ni.xml` contains `<decay>` entries without `target` attributes. The
  harness patches the missing targets identically for both Nucleide and OpenMC
  before loading (mirroring `crates/depletion/benches/depletion_bench.rs`).
- The CASL chain exercises two parser behaviors that differ across codes:
  fission-yield borrowing (`<neutron_fission_yields parent="..."/>`) and decay
  branching ratios that do not sum to 1 (OpenMC reads them verbatim from XML;
  Nucleide's `Chain::from_xml` now does the same). Both are covered by unit
  tests in `crates/depletion`.
- Nucleide's Python `Cascade` exposes both fixed-`M*` `solve()` and the
  `M*`-optimizing `solve_multicomponent()` path. `enrichment_vs_pyne.py` reports
  both the fixed-vs-optimizing comparison and the like-for-like optimizing
  comparison.
- The CRAM timing loop uses `nucleide.depletion.build_depletion_system()` and
  `DepletionSystem.solve_vec()` so that only the linear-algebra solve is timed.
- PyNE in this environment is built without PyMOAB, so `pyne.variancereduction.magic`
  cannot be invoked directly. `magic_vs_pyne.py` falls back to a pure-Python
  reimplementation of PyNE's documented MAGIC formula for the element-wise
  comparison.
