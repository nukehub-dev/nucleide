# Changelog

All notable changes to Nucleide are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Releases are cut with `scripts/bump-version.sh X.Y.Z`, which stamps the
`[Unreleased]` section below and updates `[workspace.package] version` and
workspace crate dependency versions in the root `Cargo.toml`. Git tags
(`vX.Y.Z`) are the release source of truth; CI publishes Python wheels and
workspace crates from tags.

## [Unreleased]

### Added

- Serpent and FLUKA output-parsing Python tutorials
  (`docs/tutorials/python/parse-serpent-output.md`,
  `docs/tutorials/python/parse-fluka-output.md`): read `_res.m`/`_dep.m`/
  `_det.m` files via `nucleide.serpent.read_serpent` and USRBIN `.lis`
  tallies via `nucleide.fluka.read_usrbin`, mirroring the MCNP parsing
  tutorial. Both are listed in the Python tutorial index and the docs index,
  and cross-linked from the tutorials "Finding more examples" section.
- SVG figures on all six theory pages (`docs/theory/figures/`): burnup-matrix
  column anatomy (depletion), MARC cascade schematic (enrichment), Walker
  alias-table construction (variance reduction), nucid digit layout (nuclear
  data), step-reactivity prompt-jump transient (kinetics), and peak/background
  anatomy (spectroscopy). Figures are staged into the website by
  `sync-data.mjs` and embedded with the `<ImageFigure>` shortcode.
- Theme-adaptive inline-SVG capability map on the documentation home page
  (`docs/README.md`): `currentColor` follows the site theme and
  `var(--primary)` follows the accent picker. The home intro was also rewritten
  as a reader-facing welcome with a "Where to start" section.
- Decay-line SDEF source emission (`nucleide-spectroscopy`, equation E9):
  `nucleide.spectroscopy.sdef_decay_source` takes caller-supplied decay lines
  (energy, intensity), merges duplicates, sorts, and normalizes them to
  probability bins rendered as an MCNP point-source card — PyNE field order
  (`POS=` always; `VEC=`/`DIR=1` only when directed; `WGT=`; `PAR=` via the
  `nuclei` particle dialects, MCNP versions 5/6), inline `ERG=<E>` for a
  single line, and discrete-distribution `ERG=D1` + `SI1 L`/`SP1 D` cards
  wrapped at 80 columns for several. Distribution-card syntax is
  parser-verified surface only (MCNP sampling semantics stay caller-side);
  no decay data is vendored. Validation gains E9 normalization/card gates
  plus container-gated byte-exact diffs against PyNE `PointSource.mcnp`.
- Interactive legacy-reader tutorials: Serpent output (`_res`/`_dep`/`_det`)
  and FLUKA USRBIN browser demos, and ORIGEN TAPE5/6/9 tabs in the activation
  demo — `serpent-io`, `fluka-io`, and `origen-io` are now WASM-exposed with
  live demos and E2E coverage.
- "Common workflows" section in the Python API reference (five task-oriented
  snippets in the submodule idiom), an "Optional runtime dependencies" note in
  getting started (when `numpy`/`h5py` are needed), and "See also" links from
  the depletion, enrichment, kinetics, and spectroscopy tutorials to the
  committed cross-code validation results.

### Changed

- Binding error messages for unsupported CRAM orders and unknown dose
  pathways/sources now name the accepted values (Python and WASM wordings
  aligned).

### Fixed

- Interactive demos: the Serpent parser demo no longer crashes on a
  one-column `IMP_KEFF` row, and the FLUKA/MCNP mesh views no longer
  overflow the browser argument limit on meshes above ~65k bins (min/max now
  computed with loops instead of spread).

## [0.5.0] - 2026-09-11

### Added

- Legacy I/O completion (`nucleide-mcnp-io` + `nucleide.mcnp`):
  ENDL EEDL/EPDL table reader (`endl` module: `END_OF_TABLE` framing, header
  slices, `NFIELDS_RPROP` body widths, `get_rx` selector lookup, and an
  `endftod` matching the PyNE C++ semantics) exposed as
  `nucleide.mcnp.read_endl` / `EndlLibrary.get_rx` / `endl_endftod`, with a
  synthetic `fixtures/endl/` oracle; shared Fortran unformatted-record
  framework (`fortran` module: typed `BadRecordMarker` instead of the
  upstream `AttributeError` path, little-endian only) now backing the SSW
  reader/writer; SSW combining (`combine_files` port of
  `scripts/ssw_combine.py`: signed-`orignp1` sum, plain-`nrss` sum,
  sign-preserving `nps` shift, typed incompatibility errors) exposed as
  `nucleide.mcnp.combine_ssw_files` (+ `ssw_combine_main` CLI); Python-side
  PTRAC event-row stream (`ptrac_event_rows` over `read_ptrac`, pinned
  against the 19-column `PtracEvent` schema) with an `h5py` table writer
  (`write_ptrac_hdf5`, `ptrac_to_hdf5_main` CLI; HDF5 bytes never asserted)
  and a MOAB-free meshtal mesh-data extractor (`meshtal_mesh_data`, MOAB
  tagging stays caller-side). Validation gains an ENDL-vs-`pyne.endl` oracle
  section (loud skips preserved). Known parity quirks stay documented, not
  fixed: SSW `SF_00001` split and negative-`np1`⇒table-2, `ncrd` sign⇒abs
  width, PTRAC echo-repack to 8-byte mode, meshtal single-group total mirror
  and `Rel_` column map, WWINP `{0:13.5E}` fields. Big-endian real-bytes
  input stays unsupported (no fixture oracle).
- Prescribed-reactivity point kinetics (new `nucleide-kinetics` crate):
  PKE solver for 1+ delayed-neutron groups under constant/step/impulse/
  ramp/polyline reactivity schedules in Δk, equilibrium initials, the
  inhour relation with stable-period solve, and the prompt-jump
  approximation. The integrator is an adaptive implicit θ-method
  (trapezoidal default, backward Euler on request) with exact stepping to
  schedule knots — no explicit-solver parity. All delayed data are
   caller-supplied (`from_ifp` documents the OpenMC provenance note);
   thermal feedback and flux coupling are out of scope (the WASM
   `kineticsTransient` tutorial facade is covered in its own entry below).
- Python API: `nucleide.kinetics.solve` / `equilibrium` / `initial_rate` /
   `inhour_rho` / `stable_period` / `prompt_jump` thin wrappers over the new
   core, plus `tests/test_kinetics.py` and the `validation/kinetics_vs_pyrk.py`
   two-tier oracle (synthetic-fixture analytic gates O1–O4 plus a container
   PyRK ramp cross-check at 1e-3).
- Gamma-ray spectroscopy and measurement (new `nucleide-spectroscopy`
  crate): E1 rectangular / E2 five-point smoothing (edges copied), E3
  `m == 1` background / E4 half-open gross / E5 net counts with positional
  indexing, E6 quadratic energy bins, E7 log-polynomial efficiency
  evaluation in MeV (caller coefficients only, no fitting), E8 X-ray line
  algebra over caller-supplied atomic constants (no vendored tables), and
  the dollar/plain `.spe` readers with the upstream quirks pinned
  (duplicate-tag first-wins, missing-tag errors, live-real `$MEAS_TIM:`,
  last-plus-one `$DATA:`, `keV`-suffixed triplets, ignored `$ROI:`/
  `$PRESETS:`/`$ENER_FIT:`, positional 0-based dollar labels). FWHM fits
  are parsed, never evaluated; peak search/fit, activities, decay spectra,
  and plotting stay out of scope.
- Python API: `nucleide.spectroscopy.rect_smooth` / `five_point_smooth` /
  `calc_bg` / `gross_count` / `net_counts` / `energy_bins` /
  `detector_efficiency` / `xray_lines` / `parse_dollar_spe` / `parse_spe` /
  `read_dollar_spe` / `read_spe` thin wrappers over the new core, plus
  `tests/test_spectroscopy.py` and the
  `validation/spectroscopy_vs_pyne.py` two-tier oracle (synthetic E1–E8
   gates plus a container PyNE cross-check at 1e-12).
- Browser-interactive tutorials for point kinetics and spectroscopy:
  `kineticsTransient` (step-reactivity PKE solve returning the `n(t)`
  series plus the E4 prompt-jump value) and `spectroscopySmooth` (E1
  rectangular / E2 five-point smoothing plus E3–E5 gross/background/net
  counting) thin WASM facades over the verified crate APIs, with
  `KineticsTransient` / `SpectroscopyDemo` demos and
  `tutorials/interactive/kinetics.mdx` / `spectroscopy.mdx` pages. Demo
  presets stay synthetic (`fixtures/kinetics/` + `fixtures/spectroscopy/`
  values, hand-picked small numbers otherwise).
- NumPy bridges for MCTAL/WWINP/PTRAC (bindings-only, same
  `result_array` pattern as the meshtal bridge: flatten → `Vec::into_pyarray`
  → reshape view, float64 C-order, `ValueError` on ragged/out-of-range):
  `Wwinp.ww_row_array` / `ww_column_array` / `ww_particle_array` (`(nft,)`,
  `(n_groups,)`, `(n_groups, nft)` with `nft = nf[0] * nf[1] * nf[2]`),
  `Mctal.k_arrays` (five `(n_cycles,)` series) / `averages_array`
  (`(n_cycles, 14)`, `(0, 14)` when empty), `PtracFile.events_array`
  (`(n_events, 19)` in `ptrac_event_columns` order, absent variables as 0.0)
  / `event_field_array` (`(n_events,)`). All arrays are owned, writable, and
  decoupled from the file data; the wheel stays NumPy-free (NumPy required
  only at call time, `pytest.importorskip` in `tests/test_mcnp_io.py`).
  Plain-copy helpers (`ww_row`, `ww_column`, the `k_col`/`averages` getters,
  `events`) are unchanged.

### Fixed

- Website build: the `phillips-1978` entry in
  `docs/theory/spectroscopy.mdx` carried no URL, which the docs content
  schema requires (`references.0.url`), so `npm run build` failed. The
  entry is now plain prose in the References section (citation text only,
  no journal-verification claim, no invented link).
- Stable-as-zero analytics (`nucleide-material`): a known nuclide with mass
  data but no decay constant is now stable (λ = 0) and contributes exactly
  0.0 to `activity` / `decay_heat` / `dose_per_g` instead of raising
  `MissingDecay`. Stable nuclides skip the dose-factor lookup entirely, so
  fully stable compositions such as H2O return 0 on every pathway/source.
  Nuclides with no mass data still raise `MissingMass`, and radioactive
  nuclides without a dose row (including `-1` air sentinels) still raise
  `MissingDose`. Python/WASM signatures are unchanged.

## [0.4.0] - 2026-09-10

### Added

- Python API: `MeshTally.result_array()` / `totals_array()` zero-copy NumPy
  bridge over the meshtal `result`/`rel_error` tables (supersedes the 0.3.0
  deferral): `result_array()` returns owned writable C-order float64
  `(ve, group)` arrays with `ve = (i * ny + j) * nz + k`, moved out of one
  row-major flatten via `Vec::into_pyarray` + a reshape view (no second
  copy); `totals_array()` returns `(num_ves,)` totals pairs directly.
  Requires NumPy installed at runtime (`numpy>=1.26` in the `test` extra;
  base wheel stays dependency-free). `to_list()` / `totals_list()` stay as
  the NumPy-free plain-copy path.
- License-free decay-data pack (`nucleide-nuclei`, AD-7): per-branch
  daughters from ENDF/B-VIII.0 decay tapes
  (`crates/nuclei/src/data/decay_branches.tsv`, 5 068 rows over 3 541
  parents; MF8/MT457 NDK RTYP/RFS/BR with SF/fission branches dropped and
  zero-half-life tapes absent), exposed as `decay_branches` /
  `branching_fraction` lookups plus the `DecayData` facade and the
  `DecayBranch` / `DecayBranchMode` types. Isomer masses extend
  `ame2020.tsv` (738 rows: `m_ground + ELIS/931.49410242 u` from each
  isomer tape's File-1 MT451 record; tapeless isomers fall back to the
  ground-state mass, Q-values stay ground-state-only). New free-form
  `normalize_nuclide_name` in `crates/nuclei/src/dialects.rs`
  (symbol-first, then mass-first so `N15` stays nitrogen and `92235` stays
  a ZAID; bare symbols rejected; `Ir-192n` resolves to state 2).
  `scripts/gen-nuclear-data.py --endf-decay8` regenerates both tables with
  hard ENDF spot-checks (K-40 beta-/EC pair summing to 1, Es-254 members);
  `decay_energy.tsv` stays on ENDF/B-VII.1 with its basis in its header.
- Python API: `nucleide.nuclei.decay_branches` /
  `nucleide.nuclei.decay_branch_fraction` thin wrappers over the new
  tables. WASM `decay_branches` / `decay_branch_fraction` table-level
  mirrors.
- `validation/nuclear_data_vs_refs.py`: `DECAY_BRANCH_SPOTS` oracle
  (K-40 pins 3.93839e16 s plus its beta-/EC pair; Es-254 members per the
  repo tables) following the decay-energy spot pattern; committed results
  untouched, container rerun deferred.
- `nucleide-depletion` analytic Bateman decay fast path
  (`crates/depletion/src/bateman.rs`): cached `C`/`C⁻¹` eigendecomposition
  closed form (Bateman 1910; Amaku–Pascholati–Vanin CPC 181 (2010)) for
  decay-only lower-triangular chains, generated at runtime from the crate's
  own chain data, plus an f64-careful `BatemanHp` variant (magnitude-sorted
  terms, Neumaier compensation, `exp_m1` for small `λt`; no new
  dependencies). New `Method::{Cram(Order), Bateman, BatemanHp}` selector
  (`"cram16"`/`"cram48"`/`"bateman"`/`"bateman_hp"`) threading
  `deplete_with_method` / `integrate_with_method` /
  `DecayInventory::decay_with_method` (existing `Order`-based entry points
  stay as CRAM shims; `dt == 0.0` returns the input exactly). Non-decay
  systems fall back to CRAM-48 inside the dispatcher: near-degenerate
  half-lives (`1e-12` relative gap with a live coupling path), cyclic or
  out-of-order topology, and reactions/fission on; stable nuclides use limit
  forms inline.
- Python API: `method=` on `nucleide.depletion.deplete` / `deplete_series` /
  `DepletionSystem.solve` / `solve_vec` / `Inventory.decay` (default
  `"cram48"`; an explicitly non-default `method` overrides `order`).
  WASM `deplete` / `depleteSeries` gain an optional trailing `method`.
- `validation/depletion_vs_openmc.py`: in-memory Bateman-vs-CRAM-48
  cross-check column (single-step vs CRAM-48 within `1e-6`, vs analytic
  within `1e-8`, series band within `1e-6`); committed results untouched,
  container rerun deferred.

- `nucleide-mcnp-io` L3 semantic objects (`crates/mcnp-io/src/semantic.rs`):
  typed `MODE` (37 particle shorthands, default `{N}`), `TRn` cards
  (degrees flag, 3-entry displacement, 5–9 entry rotation, main-to-aux flag,
  hidden inline `FILL` transforms), auto-created universes (`U`/`-U`,
  data-block lists with `J`/`nJ`/`nR`/`nM` expansion), `LAT` (1|2 only),
  cell `FILL` (single universe or 3-D matrix with transform reference or
  hidden transform; data-block lists are simple per-cell only), per-cell
  `IMP`/`VOL`, periodic surface pointers, and a minimal typed tally model
  (`Fn` number + particle classifier + entries with grouped `FMn`/`En`).
  `DeckProblem::validate()` centralizes duplicate-number conflicts,
  dangling material/surface/complement/transform/periodic/fill links,
  redundant cell+data definitions, write-time state checks, and
  nucleide-defined lattice/fill cross-checks (`LAT`-without-`FILL` and
  `FILL`-matrix-without-`LAT` are errors); particle/mode mismatches are
  notes via `validation_notes()`. Synthetic fixture
  `fixtures/mcnp/inp/deck_l3.txt` (no license needed).
- Python API: `DeckProblem` `mode`/`transforms`/`universes`/`lattices`/
  `fills`/`importances`/`volumes`/`tallies` getters, `validate()` /
  `validation_notes()`, and `set_mode`/`set_cell_universe`/
  `set_cell_lattice`/`set_cell_fill` setters.
- `nucleide-emit` single-material emission (`crates/emit/`): one `Material`
  renders through five code dialects — MCNP `m` cards (mass fractions,
  configurable xs suffix, 128-column packing), Serpent `mat` cards
  (`{zaid}.{lib}` ids), FLUKA `COMPOUND` cards, ALARA `mixture` blocks,
  single-zone PARTISN decks — plus a mass-drift report (`emit_drift`) with
  per-nuclide drop reasons and re-parse verification where readers exist
  (MCNP, ALARA).
- Python API: `nucleide.emit.emit_cards` / `emit_drift_table` thin wrappers
  over the new crate.
- `nucleide-emit` ARMI blueprint bridge (`crates/emit/src/armi.rs`):
  `from_armi_mass_fracs` builds a `Material` from ARMI-side post-expansion
  mass-fraction keys (database names, bare names, ZAIDs-as-strings, AAAZZZS
  ids, unambiguous MC2-3 labels via the nuclei bridge) for unchanged
  MCNP/Serpent/FLUKA/ALARA/PARTISN emission. One-way, dict-in only:
  elemental keys (caller passes expanded `massFrac`), number fractions,
  enrichment shorthands, `balance`, and temperatures stay caller-resolved;
  bare `AM242` is rejected (pass `AM242M`/`AM242G` explicitly).
- Python API: `nucleide.emit.emit_armi_cards` / `emit_armi_drift_table`
  thin wrappers over the bridge.
- `nucleide-nuclei` dose factors (`crates/nuclei/src/data/dose_factors.tsv`,
  1,116 rows: 93 folded nuclides × 4 pathways × 3 sources from the BSD-3
  PyNE `dbgen/dosefactors*.csv` tables, HNF-SD-WM-TI-707 Rev.1 / HNF-5636
  App. O; `+D` folds into the parent, GENII/DOE air are `-1` sentinels):
  `DosePathway`/`DoseSource`/`DoseEntry` lookups plus a `DoseData` facade,
  generated by the committed stdlib-only `scripts/gen-nuclear-data.py`
  (`--dose-*.csv` flags over local copies).
- `nucleide-material` dose per gram: `DoseProvider` trait (`DoseFactors`,
  `NoDoses`, plus `impl` for `nucleide_nuclei::data::DoseData`),
  `Material::dose_per_g` / `total_dose_per_g` following PyNE's
  `Material::dose_per_g` equations (per-gram map semantics; screening-level
  only, not for safety decisions).
- Python API: `nucleide.nuclei.dose_factor` and
  `nucleide.material.dose_per_g` thin wrappers over the vendored tables +
  material analytics.
- WASM `WasmDeckProblem` (`bindings/wasm/src/lib.rs`) + `deck-editor` demo
  (`website/src/components/interactive/DeckEditor.tsx`,
  `docs/tutorials/interactive/deck-editor.mdx`): `fromText` deck parsing
  (text only, no filesystem in the browser), typed `cells`/`surfs`/`mode`/
  `transforms`/`universes`/`lattices`/`fills`/`importances`/`volumes`/
  `tallies` views, `materialNumbers`/`dataNames`, `cellInventory`
  (cell→material map), `validate`/`validationNotes`, and `setCellDensity`/
  `setCellMaterial`/`setMode`/`setCellUniverse`/`setCellLattice`/
  `setCellFill` setters over the existing Tier 1 `DeckProblem` API.
  `WasmInventory` over the existing `DecayInventory` API (`units` default
  `"atoms"`, `decay` honoring rates + solver `method` via the predictor
   single-step, `activities`/`masses`/`moles`, fractions, `halfLivesReadable`,
   `add`/`sub`/`mul`/`div`, `toCsv`/`fromCsv`) plus `cumulativeDecays`/
   `progeny`/`branchingFraction`/`decayMode`/`chainEdges` free functions.
   The demo parses the inline synthetic `deck_minimal.txt`, fetches the
   synthetic `deck_l3.txt` sample staged by `sync-data.mjs`, edits cells
   through the setters, and shows validation plus byte-identical `dumps`.
- `nucleide-r2s` ARMI DB-snapshot adapter (`crates/r2s/src/snapshot.rs`):
  `deck_from_snapshot` / `R2sWorkflow::from_snapshot` /
  `snapshot_workflow` build a validated volumes-method ALARA template deck
  from versionless caller-dumped dicts (opaque zone ids 1:1 from block
  names, 1:1 zone→mixture with atoms/barn-cm number densities as element
  `vol_fraction`, caller-supplied fluxes/cooling/schedule; empty
  compositions map to `void` and are skipped). Dict-in only: no HDF5
  dependency, no ARMI layout versioning mirrored, synthetic test data only.
  Composition keys follow the emit ARMI-input rule (post-expansion nuclide
  keys; elemental keys, bare `AM242`, and unknown names are errors).
- Python API: `nucleide.r2s.r2s_from_snapshot` thin wrapper returning
  `{workflow, deck, decks}` (workflow summary, canonical template deck,
  one deck per step).
- Interactive `mcnp-io` demo (`website/src/components/interactive/McnpParser.tsx`,
  `docs/tutorials/interactive/mcnp-io.mdx`): the xsdir histogram now pairs the
  existing per-element table-count bar with a table-count bar by library suffix
  and an AWR-distribution histogram over the staged `xsdir_sample.txt`.
- Interactive `variance-reduction` demo
  (`website/src/components/interactive/VrDemo.tsx`,
  `docs/tutorials/interactive/variance-reduction.mdx`): the MAGIC energy
  upper-bound step line (`hv` shape, log energy axis) gains a group-count title.
- Interactive tutorial E2E coverage (`website/e2e/smoke.spec.ts`): follow-up
  clicks with per-click output assertions on every demo page — deck-editor
  `Validate` + a setter flow + `Load L3 sample`, variance-reduction `Sample
  index` / `Sample voxel`, materials `Mix` + `To XML`, enrichment `Optimize
  M*`, activation ALARA-output / FISPACT / R2S mode tabs, mcnp-io `wwinp`
  parse, depletion `Load sample chain`.
- WASM mirrors for the new Rust APIs (`bindings/wasm/src/lib.rs`, typed in
  `website/src/types/nucleide-wasm.d.ts`): `emitCards` / `emitDriftTable`
  over `nucleide-emit` (all five dialects, `EmitOpts` overrides, drift rows
  with per-nuclide drop reasons), `emitArmiCards` / `emitArmiDriftTable`
  over the ARMI bridge (elemental keys, bare `AM242`, and unknown names
  throw), `doseFactor` / `dosePerGram` over the vendored HNF-5636 tables
  (pathway `air`/`soil`/`ingest`/`inhale`, source `EPA`/`DOE`/`GENII`,
  screening-level only), and `r2sFromSnapshot` over the snapshot adapter
  (typed `Snapshot*Json` inputs returning a `{workflow, deck, decks}`
  bundle).
- Interactive `emitter` demo (`website/src/components/interactive/
  EmitterDemo.tsx`, `docs/tutorials/interactive/emitter.mdx`): GNDS/ARMI key
  toggle with inline presets, material name/density/per-dialect options, all
  five dialect cards plus the drift table, with inline errors for
  elemental/AM242/bad keys.
- Interactive `materials` demo dose-per-gram section (`website/src/
  components/interactive/MaterialBuilder.tsx`): pathway/source selects
  scoring the current formula via `dosePerGram`, with a screening-only note.
- Interactive `activation` demo `r2s-snapshot` mode (`website/src/
  components/interactive/ActivationDemo.tsx`): an inline two-zone JSON
  snapshot parsed via `r2sFromSnapshot` into workflow steps plus the template
  and per-step deck count.
- Interactive tutorial E2E coverage (`website/e2e/smoke.spec.ts`): emitter
  GNDS emit plus ARMI re-emit with drift assertion, materials `Compute dose`
  value, activation `R2S snapshot` tab with steps and deck-count outputs.

### Fixed

- Interactive tutorial input guards (inline `WASM error:` before any WASM
  call instead of NaN passthrough): deck-editor density / material / universe /
  lattice / fill / mode setters reject non-finite input
  (`website/src/components/interactive/DeckEditor.tsx`,
  `docs/tutorials/interactive/deck-editor.mdx`); variance-reduction tally /
  tolerance / null value plus alias-table PDF tokens reject non-finite input
  instead of silently dropping bad tokens
  (`website/src/components/interactive/VrDemo.tsx`,
  `docs/tutorials/interactive/variance-reduction.mdx`); depletion atom counts
  / time step / burnup steps require finite (and positive where applicable)
  values (`website/src/components/interactive/DepletionStep.tsx`,
  `docs/tutorials/interactive/depletion.mdx`).
- Invalid nuclide-id fallback labels instead of panics (`crates/nuclei/src/
  lib.rs`, `crates/nuclei/src/armi.rs`): `NuclideId::from_nucid` on a raw
  integer outside the validated `(Z, A, state)` domain previously panicked
  downstream — out-of-bounds `ELEMENTS[z]` indexing in `to_name`, or an
  `.expect()` on the element-symbol lookup in `nucid_to_armi_label`. Both
  are now total: such ids render the diagnostic fallback `Z{z}A{a}[m{s}]`
  (with the same `n`-prefix/capitalization treatment on the ARMI side),
  which `from_name` does not parse. New `try_from_nucid` / `is_valid`
  validate untrusted integers (`1 <= Z <= 118`, `Z <= A <= 999`, canonical
  tail) with the matching `Error` instead of producing a fallback id.
- Poisoned-lock getters now return errors instead of panicking
  (`bindings/python/src/lib.rs`): every `PyCascade` getter and
  `PyDeckProblem` getter/setter previously called `.lock().unwrap()`, so a
  poisoned `Mutex` aborted the interpreter. Each now maps the poison to a
  `ValueError` (`"cascade lock poisoned"` / `"deck lock poisoned"`).
- Cell-`FILL` matrix cap at 1 000 000 cells
  (`crates/mcnp-io/src/semantic.rs`): a hostile `FILL` range such as
  `-2147483648:2147483647` previously overflowed the `i32` width
  computation (debug panic / release wrap) and could drive a gigabyte
  allocation. Widths are now computed in `i64` with checked multiplication,
  and any matrix above the cap fails with a cap error before allocation;
  larger lattices must be built programmatically, not spelled per cell.
- Non-finite `n0` validation at the Bateman solve entry
  (`crates/depletion/src/bateman.rs`): `BatemanCache::solve` previously
  propagated NaN/infinite/negative initial counts into the triangular
  solves. It now rejects them with a `Linalg` error (`n0` must hold finite
  atom counts `>= 0`), matching the existing `dt`/length-mismatch errors.
- Snapshot-mixture emission returns an error instead of panicking
  (`crates/r2s/src/snapshot.rs`, `crates/r2s/src/error.rs`): the
  `deck_from_snapshot` mixture writer previously ended its entry match with
  `unreachable!`, so a future non-element variant would panic. It now
  returns `Error::Invalid` (the `error.rs:16` arm) with the offending
  entry attached.

## [0.3.0] - 2026-09-08

### Added

- `nucleide-depletion` time-series integrators (`crates/depletion/src/integrate.rs`):
  `Integrator::{Predictor, Cecm, Cf4}` with `integrate(sys, n0, steps, integrator,
  order)` returning per-node atoms, activity (`λ·N`), and decay heat
  (`λ·N·E`, chain `decay_energy` → ENDF/B-VII.1 table → zero).
  Synthetic A→B→C fixture `fixtures/depletion/chain_abc.xml` plus
  `tests/test_depletion_timeseries.py` covering multi-step Bateman accuracy.
- `nucleide-nuclei` data breadth (no HDF5, no hand values): `simple_xs.tsv`
  (241 rows; thermal totals from NIST NCNR bound XS as free-atom
  `xs·(A/(A+1))² + xs_a`, fast totals from ENDF/B-VII.1 MF3/MT1 at 14 MeV),
  `scattering_lengths.tsv` (267 rows from NIST NCNR Sears-1992, complex
  lengths by real part, monoisotopic element attributions),
  `decay_energy.tsv` (3,557 rows of mean *prompt* recoverable MeV/decay from
  ENDF/B-VII.1 MF8/MT457; daughter gammas live on daughter rows, so Cs137
  reports prompt 0.179 MeV and Ba137_m1 0.661 MeV — chain codes sum members).
  All three generated by the committed stdlib-only
  `scripts/gen-nuclear-data.py` (run `--help` for upstream URLs); resonance
  absorbers without NIST rows and isomers stay absent by construction.
- `nucleide-material` decay heat: `DecayEnergyProvider` trait
  (`DecayEnergies`, `NoDecayEnergies`, plus `impl` for
  `nucleide_nuclei::data::DecayData` so depletion gets λ and MeV
  from one provider with no circular deps), `Material::decay_heat` /
  `total_decay_heat` in watts, and `MEV_TO_JOULES`. Dose coefficients stay a
  documented deferral (PyNE's `dbgen/dosefactors*.csv` are the future source).
- `nucleide-depletion`: single decay-energy source of truth — the hand-kept
  fallback table is removed and `decay_energy_mev` resolves chain →
  ENDF/B-VII.1 table → `0.0`; new shared `decay_energies_by_name` helper.
- Python API: `nucleide.depletion.deplete_series` — multi-step depletion
  driver over the core `integrate` series (`predictor`/`cecm`/`cf4`;
  per-step `rates`/`rates_list`), returning `times`/`atoms`/`activity`/
  per-nuclide `decay_heat` maps through the shared resolver (identical to
  core `integrate` heats; the core `t = 0` row is omitted).
- Python API: `nucleide.nuclei` (`simple_xs`, `scattering_length`,
  `decay_energy`) and `nucleide.material` (`decay_heat`) thin wrappers over
  the vendored tables + material analytics.
- Python API: `MeshTally.to_list()` / `totals_list()` plain-copy helpers;
  the zero-copy NumPy bridge (`result_array()`) stays deferred — no `numpy`
  dependency introduced.
- `nucleide-mcnp-io` full-deck round-trip (`crates/mcnp-io/src/cell.rs`,
  `surf.rs`, `problem.rs`): typed cell/surface cards (CSG expressions with
  `:` union, juxtaposition intersection, parentheses, `#` complement, `*`
  reflecting markers; 40 surface kinds incl. macrobodies), `DeckProblem`
  parse/edit/write with format-preserving write-back (byte-identical when
  unedited; only touched cards re-render canonically) plus
  `set_cell_density` / `set_cell_material`. Synthetic fixtures
  `fixtures/mcnp/inp/deck_{minimal,macro,params}.txt` (no license needed).
- `nucleide-depletion` unit-aware inventories (`crates/depletion/src/
  inventory.rs`): `DecayInventory` over atom counts with activity/mass/mole
  units (`Bq`…`Ci`, `g`/`kg`, `mol`, time in `s`/`m`/`h`/`d`/`y`),
  fractions, readable half-lives, arithmetic, CSV round-trip, plus
  `cumulative_decays` and chain observables (`progeny`,
  `branching_fraction`, `decay_mode`, `chain_edges`).
- `nucleide-nuclei` ARMI dialect bridge (`crates/nuclei/src/armi.rs`):
  `armi_name_to_nucid` (`nU235`, ZAIDs, MCC3, AAAZZZS),
  `nucid_to_armi_label`, `mcc3_to_nucid`.
- `nucleide-material` composition checks (`crates/material/src/check.rs`):
  `check_labels` (truncated-label collisions at DIF3D/MC2 6/8-char widths
  across GNDS/ZAID/Serpent/ALARA/ARMI labels) and conservation `audit`
  (non-positive totals, negative masses, duplicates, unknown masses).
- Python API: `nucleide.mcnp.DeckProblem` (`parse_deck`/`read_deck`,
  `dumps`, setters, cell/surf/data accessors), `nucleide.depletion.
  Inventory` + `cumulative_decays`/`progeny`/`branching_fraction`/
  `decay_mode`/`chain_edges`, `nucleide.nuclei.armi_to_nucid`/
  `nucid_to_armi`/`mcc3_to_nucid`, `nucleide.material.check_labels`/
  `audit_material`.

## [0.2.0] - 2026-09-06

### Added

- `cccc-io` crate: CCCC text-subset parsers (ISOTXS multigroup libraries,
  RTFLUX/ATFLUX/RZFLUX flux files) + PARTISN deck writer with ISOTXS nuclide
  mapping and validation (no solver). Synthetic fixtures under
  `fixtures/cccc/` (no license needed).
- `fispact-io` crate: FISPACT-II inventory-output parser producing
  ALARA-compatible response frames (output-only, reusing `ResponseFrame`;
  no solver). Synthetic fixture `fixtures/fispact/inventory.fis`
  (no license needed).
- `origen-io` crate: scoped ORIGEN 2.2 TAPE readers — `TAPE5` input echo,
  `TAPE6` output inventory, `TAPE9`-style decay constants. Synthetic fixtures
  under `fixtures/origen/` (no license needed).
- `r2s` crate: scoped R2S workflow builder — zone-to-flux linking from ALARA
  decks, schedule expansion, per-zone deck emission, and photon-source
  assembly with a documented uniform-split approximation (total strength
  preserved; group-wise emission lines stay in ALARA `.photonSrc` spectra).
  Transport and activation solving stay out of scope.
- Python API: `nucleide.cccc` (`isotxs_parse`, `rtflux_parse`,
  `partisn_render`, `partisn_validate`), `nucleide.fispact`
  (`fispact_parse_output`), `nucleide.origen` (`origen_parse_tape5/6/9`),
  and `nucleide.r2s` (`r2s_from_deck`, `r2s_validate`, `r2s_expand`,
  `r2s_assemble`) exposing the four new crates as plain dicts/lists.
- `alara-io` crate: ALARA input-deck, group-flux, material/element/WDR library,
  activation-output, photon-source, and schedule-expansion parsers (glue only;
  solver out of scope). Fixtures vendored verbatim from UW ALARA BSD samples
  under `fixtures/alara/` (terms in `fixtures/alara/LICENSE.ALARA`).
- Python API: `nucleide.alara` (`alara_parse_deck`, `alara_parse_flux`,
  `alara_parse_output`, `alara_expand_schedule`) exposing the `alara-io`
  deck/flux/output/schedule glue as plain dicts/lists.

### Changed

- Release workflow: removed the non-functional `cargo publish --dry-run` step
  for workspace crates and added `set -euo pipefail` to the publish loop so a
  failure in one crate stops the job.

### Fixed

- Workspace crate dependencies now specify a `version` requirement so
  `cargo publish` accepts them for crates.io.
- `scripts/bump-version.sh` now updates workspace crate dependency versions
  alongside `[workspace.package]` version.

## [0.1.0] - 2026-08-30

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
- Interactive tutorials: MCNP meshtal z-slice heatmap with log colorscale and
  relative-error toggle, xsdir table-count histogram by element, depletion
  burnup-curve plot, material atom-vs-weight-fraction bar chart, enrichment
  cascade stage-profile line chart, and MAGIC energy upper-bound step line.
- `enrichment`: `Cascade::stage_profile()` returns per-stage assays of the
  enriching key via the ideal-cascade recurrence, exposed through WASM as
  `WasmCascade.stageProfile()`.
- Release workflow: `vX.Y.Z` tags build and publish Python wheels for Linux,
  macOS, and Windows to PyPI, publish workspace crates to crates.io, and draft
  a GitHub release from the matching changelog section.

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

- Documentation website "Edit this page on GitHub" links now point at real
  repo files (docs-kit sync injects the source path as `editPath`
  frontmatter); previously every link 404'd — all were missing the file
  extension, and the changelog link pointed at `docs/changelog` instead of
  the repo-root `CHANGELOG.md`.
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
