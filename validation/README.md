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
| `magic_vs_pyne.py` | MAGIC weight-window generation vs PyNE (or equivalent formula when PyMOAB is missing), plus always-run OpenMC/Serpent weight-window emission probes (container OpenMC load cross-check; Serpent load is a loud SKIP) |
| `nuclear_data_vs_refs.py` | Atomic masses, abundances, half-lives, name-dialect conversions vs PyNE/OpenMC |
| `dose_vs_pyne.py` | Dose coefficients (`dose_factor`/`dose_per_g`) vs PyNE accessors |
| `emit_vs_self.py` | Five-dialect emission + drift self-consistency (incl. ARMI-key equivalence; no external oracle exists) |
| `parsers_vs_refs.py` | Serpent/MCNP/FLUKA parser cross-validation vs serpentTools and PyNE oracles |
| `activation_vs_refs.py` | ALARA/CCCC/FISPACT/ORIGEN/R2S I/O checks, analytic regrouping/volume vectors, optional PyNE overlap comparison + synthetic self-consistency |
| `decay_vs_radioactivedecay.py` | Single-nuclide decay (H-3, Co-60, Cs-137) vs the radioactivedecay oracle (ICRP-107, container-only) |
| `unfold_vs_analytic.py` | SAND-II neutron spectrum unfolding: synthetic forward-fold round-trips + IRDFF-II analytical benchmark shapes (published facts, cited; no external unfolding oracle exists) |
| `equilib_vs_desc.py` | Equilibrium data readers: synthetic classic-`wout` magic/format/Jacobian gates + `&INDATA` grammar gates (always run) with DESC/simsopt container probes (loud SKIP outside) |
| `actinv_handoff_probe.py` | Optional pinned ACTINV CLI synthetic export/readback, with binary/source hashes and conservation receipts; no activation solve |
| `tritium_vs_festim.py` | 1D tritium kernel parity probe vs FESTIM 2.1 (black-box; probe-only image; see probe section below) |
| `timings.py` | Coarse wall-time comparisons (Python overhead included) |
| `make_figures.py` | Generates the paper figures in `figures/` from the JSON reports |
| `common.py` | Shared helpers and `Report` class used by the comparison scripts |
| `render_results.py` | Renders `results.md` from the JSON reports in `results/` |
| `environment.json` | Hand-maintained source annotations for the environment table |
| `run_all.sh` | Run the whole harness with a configurable Python binary |
| `run_container.sh` | One-command container run: rebuilds the image (layer-cached), then the wheel, then the harness |
| `Containerfile` | Validation environment definition (PyNE 0.7.5 + OpenMC 0.16.0 + serpentTools 0.11.0 + radioactivedecay 0.6.1, Python 3.12) |
| `magic_tally.txt` | Small synthetic MCNP meshtal used by the MAGIC comparison |
| `results/` | Per-run JSON reports produced by each comparison script |
| `results.md` | Committed measured results (generated; do not hand-edit) |
| `figures/` | Generated paper figures (`timings.png`, `depletion_agreement.png`, `kinetics_transient.png`, `spectroscopy_overlay.png`); committed, never hand-edited |
| `.cache/` | Runtime data, source-pinned probe caches, and local run/benchmark receipts (git-ignored) |

## Environment

The canonical environment is the container built from `Containerfile`:

- Python 3.12 in a dedicated `/opt/validation` prefix on `condaforge/miniforge3`
- PyNE 0.7.5 (`nomoab_noopenmc` build from conda-forge, so `pyne.mcnp.Meshtal`
  is unavailable). Note `pyne.__version__` reports `0.7.1` — upstream's Python
  metadata lags the actual release; the conda package version (0.7.5) is the
  truth. 0.7.5 is numerically identical to PyNE 0.7.8 for the modules used
  here — the 0.7.5..0.7.8 diff in `src/enrichment.cpp`, `src/data.cpp` and
  `pyne/dbgen/atomic_mass.py` is quote style and `std::isnan` qualification
  only. PyNE 0.7.8 itself cannot build on Python 3.12 (its `setup.py` still
  does `import imp`, removed in 3.12), and conda-forge tops out at 0.7.5.
- OpenMC 0.16.0 built from the upstream release tag. conda-forge now ships
  0.16.0 as well (previously 0.15.3) — use 0.16.0 everywhere: 0.15.3's
  `Chain.form_matrix` drops decay constants on some diagonal entries (52 on
  the CASL chain, up to 0.26% relative), which fails the 1e-8 CRAM gate;
  0.16.0 includes them and matches to 8.9e-15. (The
  `docker.io/openmc/openmc:v0.16.0` image actually contains 0.15.3 — its
  Dockerfile defaults `openmc_branch=master` and the published tag was built
  from the 0.15.3 release commit — so that image is not a substitute.)
- Nucleide abi3 wheel built in release mode by the `ghcr.io/pyo3/maturin`
  container (manylinux2014; release mode is required for meaningful timing
  numbers).
- serpentTools 0.11.0 from PyPI, as the independent oracle for the Serpent
  parser comparison (conda-forge's `serpent-tools` build lags).
- radioactivedecay 0.6.1 from PyPI, as the independent oracle for the
  single-nuclide decay comparison. Its ICRP-107 decay data are read only at
  harness runtime inside the container — never vendored into fixtures, code,
  or docs (licensing boundary).

A conda env with `pyne` and `openmc` from conda-forge also works (see
`run_all.sh`), but must carry OpenMC 0.16.0 — 0.15.3 fails the depletion
gate (see Environment). Caller-side oracle extras (`radioactivedecay`,
`pypact`, `openmc-plasma-source` + NeSST from PyPI) additionally cover the
decay, clearance-overlap, and plasma-moment oracle legs locally; DESC and
simsopt stay container-absent optionals (import-surface probes only).

## Running

One command with Docker or Podman (the first build installs the oracle
packages and compiles OpenMC; later builds reuse the layer cache):

```bash
./validation/run_container.sh
# To choose Docker explicitly:
CONTAINER_ENGINE=docker ./validation/run_container.sh
```

Automatic selection prefers installed Podman, then Docker. An explicit missing
or failing engine raises an error. Both use the same Containerfile and
comparison scripts. Runtime jobs have 6 GiB memory, no swap, 128 tasks and
two CPUs; Rust and OpenMC compilation use two jobs. Numerical library threads
are capped at two, so keep these settings with timing comparisons. The image
build sends only `validation/` as context; the wheel build bind-mounts the checkout.

Or against any existing Python that has PyNE, OpenMC and the Nucleide release
wheel installed:

```bash
./validation/run_all.sh /path/to/python
```

Scripts with fixed acceptance gates exit nonzero when those gates fail;
other quantitative differences are reported as comparisons rather than
universal equality claims. Read individual table statuses and explicit skips.
Native timing rows identify Criterion artifacts and their hashes. If raw
artifacts are absent, historical fallback estimates are labeled as such;
they are not presented as measurements made by the container run.

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

`nuclear_data_vs_refs.py` additionally downloads the EPA Federal Guidance
Report No. 15 external-dosimetry coefficient zip (EPA 402-R-25-001, July
2025; data file dated 2025-05-28) from
`https://www.epa.gov/system/files/other-files/2025-07/fgr15_data_2025_05_28.zip`,
cached at `validation/.cache/fgr15_data_2025_05_28.zip` on first use and
re-downloaded only if absent. The download is pinned to SHA-256
`71314b3f1d73c197da8b589e74c450f61befce48c66ace6ac1f6e0597560bd91`
(retrieved 2026-09-14); a hash mismatch fails loudly because EPA may have
revised the file. As a U.S. federal government work the FGR 15 tables carry
no copyright, but they are still not vendored — users and the harness fetch
from EPA directly.

`damage_vs_specter.py` additionally downloads the SPECTER manual (Greenwood
& Smither, ANL/FPP/TM-197, January 1985 — a U.S. government work, public
domain) from OSTI at `https://www.osti.gov/servlets/purl/6022143`, cached at
`validation/.cache/specter_anl_fpp_tm_197.pdf` on first use and re-downloaded
only if absent. The download is pinned to SHA-256
`489bca1482ffb33c160329283c0bbb2245841d2808458749ff431caf8d89ce44`
(retrieved 2026-09-15); a hash mismatch fails loudly. The oracle uses the
container's `pdftotext` (poppler-utils) to anchor a handful of transcribed
output spots against the report text; the 41-isotope displacement tables are
never vendored — SPECTER stays an oracle, not a data source.

`parsers_vs_refs.py` additionally downloads the GDML schema (Geant4
`v11.4.2` tag, schema version 3.1.7 — the Geant4 Software License is
BSD-style) from the Geant4 source mirror at
`https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/persistency/gdml/schema/`,
cached at `validation/.cache/gdml-schema-v11.4.2/` on first use and
re-downloaded only if absent. Each of the eight schema files is pinned to a
SHA-256 hash recorded in the script (retrieved 2026-09-15); a hash mismatch
fails loudly. The schema is a format definition used only to validate
Nucleide's emitted GDML documents — never vendored into the repository. The
GDML format is described by the CHEP 2005 paper (CERN-CDS-1023367).

The 0.3.0 oracle additions need no new downloads: `nuclear_data_vs_refs.py`
reads PyNE's bundled `nuc_data.h5` (`/neutron/simple_xs` via
`pyne.xs.data_source.SimpleDataSource`) in place, checks scattering lengths
against hardcoded NIST NCNR anchors, and checks decay energies against
hardcoded chain/ENDF-derived spot values; `depletion_vs_openmc.py` builds its
synthetic A→B→C chain inline in a temp file (no new fixtures).

## Optional ACTINV synthetic exchange probe

`actinv_handoff_probe.py` checks the released ACTINV 1.7.2 Linux exporter
against a freshly rebuilt Nucleide extension. It requires an explicit
`--actinv` path and verifies the pinned binary SHA-256 before execution;
it downloads neither software nor evaluated data. This is source-format
interoperability evidence, not independent activation-physics validation,
so it is separate from `run_all.sh`'s cross-code comparisons.

After `maturin develop`, run on a systemd host with a verified local binary:

```bash
systemd-run --user --scope --quiet \
    -p MemoryMax=6G -p MemorySwapMax=0 -p TasksMax=128 -p CPUQuota=200% \
    .venv/bin/python validation/actinv_handoff_probe.py --actinv /path/to/actinv
```

The pinned exporter is a leaf process; its 30 s timeout kills and reaps that
child. Explicit missing paths, digest mismatches, malformed exports, timeouts,
and conservation errors fail. Each invocation creates a fresh directory under
`validation/.cache/actinv-handoff` and retains its JSON success/failure receipt,
including software/source hashes, actual cgroup v2 limits when available,
and normalization/time selections.

The synthetic cell has 2 cm³ volume, groups `[4, 8]` photons/s, source time
3900 s and shutdown 300 s. Exported densities `[2, 4]` at 3600 s cooling
reconstruct 12 photons/s; volume splitting over `[0.5, 1.5]` cm³ returns
`[3, 9]` photons/s. A decoy cooling row checks time selection. Centroid labels
remain labels, not transport-bin edges. This synthetic probe remains separate
from the evaluated-data calculation below.

## Fresh ACTINV activation comparison

The external example in [ACTINV PR 6](https://github.com/AvilaLabs/ACTINV/pull/6),
revision `6f8f893fdfc1dbe76257705d54f77b1ee68c633b`, was executed with the
hash-pinned 1.7.2 binary above, the rebuilt Nucleide release wheel,
NumPy 2.5.3, Python 3.14.8, and all six verified catalog 1.2.0 artifacts.
The adapter remains upstream-owned; no activation driver or evaluated data
is vendored into Nucleide. It is an optional interoperability check, separate
from the independent physics comparisons in `run_all.sh`.

To reproduce, obtain `contrib/nucleide_r2s/` from that pinned revision into
`validation/.cache/`, install NumPy 2.5.3 and the rebuilt Nucleide wheel in
an isolated environment, and install the external bundle with
`actinv data fetch --output /absolute/data-root` followed by
`actinv data verify --output /absolute/data-root`.
Keep the original example and patch a separate copy of `e2e.py`: change
only its catalog guard and corresponding error message from `1.1.0` to
`1.2.0`. Its four required activation/index/decay hashes match the new catalog.
The original script SHA-256 is
`8911f4f7e7ab8742258e2e9c4a88a262627559eeed4dda5a8579362b5d16d689`;
the patched script SHA-256 is
`4dbddd670cf6b55751737f18122adcc17ae612077a27e0aaeac5b2f101e4e528`.
Retain the diff and all data hashes alongside the report.

```bash
ACTINV_CACHE_DIR=/absolute/prepared-cache \
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 RAYON_NUM_THREADS=2 \
systemd-run --user --scope --quiet \
    -p MemoryMax=6G -p MemorySwapMax=0 -p TasksMax=128 -p CPUQuota=200% \
    timeout --kill-after=5s 360s /absolute/env/bin/python /absolute/example/e2e.py \
    --actinv-bin /absolute/actinv --data-root /absolute/data-root \
    --inputs /absolute/original-example/inputs --out /absolute/new-run
```

The upstream script bounds each solver/export subprocess to 60 s and checks
the full native 709-group input grid. The Fe-56 case (15.748 g, 2 cm³,
10¹² n/(cm²·s), 300 s irradiation and 3600 s cooling) produced
463393807.7641661 photons/s and 20622094.12329952 Bq/g of Mn-56.
The translated and direct solves agreed exactly across all 24 photon groups
and Mn-56 activity (relative tolerance `1e-9`, absolute `1e-12`).
ALARA text readback contains 18 groups because six trailing zero groups are
omitted; density × volume conserved the total to relative error `1.2863e-16`
(tolerance `1e-10`). This is same-solver interoperability evidence, not
experimental validation or a test of arbitrary regrouped spectra.

Local receipts live under `validation/.cache/actinv-solve-018/`:
`fresh-run-01/report.json` records inputs, data hashes, results and versions;
`example-provenance.json` and `catalog-compat.patch` record the compatibility
change. The wheel extension SHA-256 was
`6f108cc47492f3eaf89ffa3553f87578647548c4000535f50754cc8370bd7f12`.

## Tritium/FESTIM parity probe (probe-only image)

`tritium_vs_festim.py` compares the `nucleide.tritium` 1D kernel against
FESTIM 2.1 (Apache-2.0, `festim==2.1`, v2.x API) black-box only: the script
drives FESTIM's documented surface (`HydrogenTransportProblem`,
`Mesh1D`, `Material`, `Species`/`Reaction`, BC classes, `Settings`,
exports) and reads only the exported flux and profile numbers. No FESTIM
code is read, copied, or vendored, and no synthetic-input rule is bent —
every probe parameter is a hand-picked value in the script.

FESTIM 2.1 needs DOLFINx, so the probe never runs in the main
`Containerfile` image (which stays FESTIM-free by contract) and never in
the Rust workspace. It runs in the separate probe-only image
`dolfinx/dolfinx@sha256:2ae4bfbc0d9077268880faf04c72750528bee986c94ab223a2c159969bd56fa8`.
Docker and Podman also support the separate probe below. Build the release
wheel with `run_container.sh` first. The probe image has no Rust toolchain;
install that wheel so every comparison runs the actual rebuilt artifact:

```bash
FESTIM_IMAGE='dolfinx/dolfinx@sha256:2ae4bfbc0d9077268880faf04c72750528bee986c94ab223a2c159969bd56fa8'
docker pull "$FESTIM_IMAGE"
FESTIM_IMAGE_ID="$(docker image inspect --format '{{.Id}}' "$FESTIM_IMAGE")"
docker run --rm --memory=6g --memory-swap=6g --pids-limit=128 --cpus=2 \
    -e NUCLEIDE_FESTIM_IMAGE_ID="$FESTIM_IMAGE_ID" \
    -v "$PWD:/work" -w /work \
    "$FESTIM_IMAGE" \
    bash -c 'pip install --no-deps target/container-wheels/nucleide-*-manylinux*x86_64.whl && pip install festim==2.1 scifem==0.27.0 io4dolfinx==1.5.0 fenicsx-compat==0.4.0 numpy==2.4.6 && python validation/tritium_vs_festim.py'
```

The probe prefers an installed package. Only when Nucleide is absent does it
load the mounted checkout's extension, which must first be rebuilt using
`maturin develop`. A broken installed package is never silently replaced.
The runner supplies the inspected immutable image ID for cache provenance;
the expected image name alone cannot establish the running environment.
For Podman, replace `docker` in the pull, inspect and run commands.
After the probe, regenerate the combined tables and figures in the main image:

```bash
docker run --rm --memory=6g --memory-swap=6g --pids-limit=128 --cpus=2 \
    -e OMP_NUM_THREADS=2 -e OPENBLAS_NUM_THREADS=2 -e MKL_NUM_THREADS=2 \
    -v "$PWD:/io" -w /io localhost/nucleide-validation:latest \
    bash -c 'pip install --no-deps target/container-wheels/nucleide-*-manylinux*x86_64.whl && python validation/render_results.py && python validation/make_figures.py'
```

Progression: F1 trap-free Dirichlet-Dirichlet → F2 one McNabb-Foster trap
→ F3 Sieverts ends → F4 two-layer stack → F5 trap-free recombination outlet.
F1–F4 gates: J_ss 2%, J(t) 3% above a `1e-3·J_ss` floor excluding
t < 0.2·t_lag, t_lag 5%, steady c(x) 2% max-norm; each side self-converges
to < 1/3 of its gate first, otherwise the leg is inconclusive, never a
pass. Always-run legs VV Val-1da (analytic permeation) and M1
(manufactured-solution order) gate the Nucleide side in any environment.
Outside the probe image the FESTIM legs SKIP loudly with their reason and
are recorded as executed probes, never as silent passes. `results.md`
numbers for this section are regenerated only by a probe-image run. The
image digest is pinned; record the tested git commit and actual package
versions with the probe results. The tested image provides DOLFINx
0.11.0.post0 and Python 3.12.3.

F5 uses a 1 mm slab, diffusivity `1e-9 m²/s`, inlet `1 mol/m³`, zero
initial concentration and outlet law `J = Kr c²` with
`Kr = 1e-6 m⁴/(mol·s)`. The scalar positive-root steady solution is
compared against each code. FESTIM receives atomic concentrations and
fluxes: multiply molar values by exact SI `NA = 6.02214076e23 mol⁻¹`.
Its two-identical-reactant `SurfaceReactionBC` consumes two atoms per
reaction, so the reaction coefficient is `Kr/(2 NA)`, with dissociation
zero. The outward-positive derivative-based `SurfaceFlux` is compared
directly; it is never replaced by a concentration-derived flux.

F5 uses quadratic elements and atomic residual tolerance `NA*1e-14`
(relative tolerance `1e-10`). Mesh refinement is 50/100/200 cells at
a nominal `dt=0.25 s` cap; time refinement uses caps
`dt=0.5/0.25/0.125 s` at 200 cells, separately in each code. Both begin
at `1e-6 s`, growing steps by `1.02` up to that cap to resolve the
initial Dirichlet boundary layer. The report records the actual grid. Transient comparison times are
`t/(L²/D) = {0.2,0.5,1,2,5}`. Flux tolerance is 3% with a
`1e-3 J_ss` denominator floor; steady and common-grid profile tolerances
are 2%. Refinement changes must be less than one-third of each gate
before parity can pass. The last 10% of the trajectory must also show
less than `2%/3` relative flux change; no Dirichlet time-lag formula is
used for F5.

Mobile inventory change is checked against integrated inlet minus outlet
flux on the actual exported time grid (relative residual at most `1e-3`,
normalized by the larger of inventory change and absolute throughput).
The report records the grid, initial exported inventory, integrals and
independent stride-2/stride-4 flux-quadrature changes; each quadrature
change must be below `1e-3/3`. Nucleide inventory integrates cell averages;
FESTIM inventory uses the public `TotalVolume` integral export for its
quadratic finite-element field. Its final-time proximity rule is handled by
solving two nominal timesteps beyond the requested end, so actual exports
bracket every requested sample. Profile/scalar checks use the requested
end; inventory balance uses the complete actual trajectory, including this
guard interval. No missing final value is extrapolated. Profile-export sample positions are not
assumed to be integration nodes.
Unconverged refinement is INCONCLUSIVE and exits nonzero.
Calculation/export failures and broken installed environments are FAIL for
all legs; missing FESTIM is a named SKIP.
Successful reduced F5 calculations can be cached under `.cache/f5-runs`.
Reuse requires `NUCLEIDE_FESTIM_IMAGE_ID` from the runner's inspect command,
installed package versions (including PETSc's Python bindings, UFL, Basix
and FFCx), relevant solver/thread settings and hashes of loaded native
libraries, including PETSc and DOLFINx. These augment the probe/cache-helper/facade/extension
hashes and model parameters. Keep the environment fixed during each run.
Missing or incomplete runtime identity disables both persistent cache reads
and writes; calculations still run normally. Old cache entries cannot match
the new identity. Acceptance verdicts are recomputed on reuse.
Adding the probe does not establish a measured F5 pass: regenerate
the report in the pinned image before quoting its results.

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
- **MCNP CSG translation** — OpenMC 0.16.0 fully loads every translated GO
  deck `geometry.xml` (`Geometry.from_xml` with stub materials for the emitted
  material numbers); lattice decks are additionally verified geometrically
  against the source deck (element↔`FILL`-matrix mapping in MCNP `k, j, i`
  order, element-centroid containment, `pitch * dimension` span vs the lattice
  cell RPP bounds). When `openmc` is not importable these probes skip loudly;
  structural Serpent/PHITS probes always run.
- **FLUKA** — no working oracle: `pyne.fluka.Usrbin` needs PyMOAB and reads
  only binary USRBIN, while our fixtures are ASCII `.lis`. Skipped with an
  explicit note.

## Activation oracle notes

`activation_vs_refs.py` checks the `alara-io`, `cccc-io`, `fispact-io`,
`origen-io` and `r2s` crates on committed `fixtures/` inputs:

- **ALARA** — upstream sample decks/listings (`sample2`, `sample3`,
  `sample2.out`, `fluxin2`) with a container-only `pyne.alara` deck block-set
  probe (API-guarded: it compares block sets only where PyNE exposes a usable
  deck entry point, and skips loudly otherwise).
- **FISPACT / CCCC / ORIGEN / R2S** — no independent oracle exists (PyNE ships
  no matching readers), so the synthetic fixtures are self-consistency checks:
  row/variable counts, spot values, totals, plus PARTISN render/validate and
  R2S workflow smoke tests.
- Approximations recorded in the report: the uniform photon-group split in
  `r2s_assemble` (total strength only), the FISPACT `half_life_s == -1.0`
  missing-data sentinel, the CCCC text-analog subset (production files are
  binary), and the ALARA year conventions (365.25 d in decks, 365 d in the
  output listing) kept verbatim.

## Decay oracle notes

`decay_vs_radioactivedecay.py` compares Nucleide depletion decay analytics
(CRAM-48 parent atoms + `Inventory` parent activity) against the
`radioactivedecay` oracle (default ICRP-107 dataset, pip-pinned in
`Containerfile`) on synthetic single-nuclide cases (H-3, Co-60, Cs-137)
built inline in a temp file — no fixtures, no downloads:

- **Licensing boundary**: ICRP-107 decay data are read only at harness
  runtime inside the container. No half-life or branching value from the
  oracle is vendored into fixtures, code, or docs; every oracle number in
  the report is a live query.
- **Honest tolerances**: Nucleide uses ENDF/B-VIII.0 half-lives, the oracle
  ICRP-107, so exact agreement is not expected. The gate is on the
  table-corrected residual (measured ratio vs the ratio predicted from the
  two runtime-read half-lives); raw gaps are findings, and the half-life
  deltas themselves are reported.

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
