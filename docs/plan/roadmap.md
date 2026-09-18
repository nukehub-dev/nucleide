---
title: Roadmap
sidebar:
  order: 0
---

Nucleide is pre-alpha. APIs may change without notice.

## Current status

The workspace is bootstrapped with twenty-five crates, PyO3 and WASM bindings, a
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
  glue (solver out of scope), plus clearance / waste-classification
  analytics: clearance index CI = Σ Aᵢ/CLᵢ and sum-of-fractions screening
  over parsed inventories, caller-supplied limit tables with the EU
  2013/59/Euratom Annex VII Table A vendored default, plus the Spanish CSN
  conditional NORM landfill tables as explicit opt-in tables (screening
  arithmetic, never a compliance decision).
- CRAM depletion solver and chain XML parsing (`depletion`).
- Multicomponent enrichment cascade solver (`enrichment`).
- MAGIC weight windows with OpenMC (`settings.xml`) and Serpent (WWINP
  via `wwin ... wf 2`) emission, plus mesh source sampling (`vr-tools`).
- CCCC text-subset parsers + PARTISN writer (`cccc-io`): ISOTXS/RTFLUX
  readers and deck validation (no solver).
- FISPACT-II inventory output parser (`fispact-io`): output-only, reusing the
  ALARA response frame, plus the clearance-bearing wide inventory table
  printed with the `HAZARDS` + `CLEAR` keywords.
- Scoped ORIGEN 2.2 TAPE5/6/9 readers (`origen-io`).
- R2S orchestrator landed as a scoped workflow builder (`r2s`): zone-to-flux
  linking, schedule expansion, and uniform-split photon assembly (transport
  and activation solving stay out of scope).
- Typed Python facade and `.pyi` stubs (`python/nucleide/`).
- `MeshTally.result_array()` / `totals_array()` zero-copy NumPy bridge
  (`mcnp-io`): owned C-order float64 `(ve, group)` views over one row-major
  flatten (plus `(num_ves,)` totals); `to_list()` stays the NumPy-free path.
- Criterion benchmarks for CRAM, cascade solving, and parser throughput
  (`crates/*/benches/`).
- Cross-code validation harness against PyNE and OpenMC (`validation/`).
- Activation-code I/O validation (`validation/activation_vs_refs.py`):
  ALARA/CCCC/FISPACT/ORIGEN/R2S checks on committed fixtures vs PyNE oracle
  probes + synthetic self-consistency.
- Documentation website with interactive WASM tutorials (`website/`).
- Five-dialect card emission with mass-drift report (`emit`).
- ARMI blueprint-key bridge to emission (`emit`).
- Dose factors (`nuclei`) and dose per gram (`material`).
- L3 semantic deck objects with validation (`mcnp-io`).
- Analytic Bateman fast path with `method=` selector (`depletion`).
- Deck-editor demo and WASM decay inventory (`bindings/wasm`, `website`).
- ARMI database-snapshot adapter (`r2s`).
- License-free ENDF/B-VIII.0 decay-data pack: branches, isomer masses,
  free-form name normalizer (`nuclei`).
- `MeshTally.result_array()` / `totals_array()` NumPy bridge (`mcnp-io`).
- Library/inventory charts in the interactive demos (`website`).
- Hardened guards: invalid-id fallback labels, poisoned-lock errors,
  cell-`FILL` cap, Bateman `n0` validation, snapshot-mixture error.
- Rust API stability pass (all crates): every public error enum is
  `#[non_exhaustive]` (new variants are no longer breaking), crate-root
  error re-exports and `Result` aliases complete the surface, and all 25
  crates carry `#![warn(missing_docs)]`.
- Point-kinetics solver (`kinetics`): prescribed-reactivity PKE with inhour
  and prompt-jump analyses (no transport, no feedback).
- Spectroscopy toolkit (`spectroscopy`): smoothing, counting, calibration,
  X-ray, and SPE line-list readers.
- UQ-lite sampling kernel (`linalg`): seeded MVN/lognormal/LHS draws over
  caller-supplied blocks, plus `decay` branch/energy/fission-yield perturbers
  (`perturb_fission_yields` closes the last named-open hook).
- License-free ENDF/B-VIII.0 fission-yield pack (`nuclei`): independent and
  cumulative yields verbatim, 36 parents, chain-XML depletion fallback.
- Tritium 1D diffusion-trapping kernel (`tritium`): T1–T2 mobile/trap
  transport with the full surface taxonomy; recombination ends closed in
  steady state and transient (G5/G6); multi-layer series stacks with
  Sieverts or Henry internal interfaces (G7/G8/G9) plus vented-sink
  recombination gaps (G10/G11).
- Scoped MCNP→OpenMC/Serpent/PHITS/GDML CSG translation (`csg-xlate`): surfaces,
  cells, nested universes, and rectangular `LAT=1` lattices with per-item
  drift reports (GDML lattices expand to per-element placements).

- Gaussian KDE source resampling (`vr-tools`): `KdeSampler` fits over caller
  particle vectors with deterministic draw/pdf.
- Typed legacy SDEF fixed-source reader (`mcnp-io`): keyword and discrete
  `SI`/`SP`/`SB` subset with canonical re-emission.
- EPA FGR 15 runtime-download dosimetry pack (`nuclei` + `nucleide.data`):
  hash-pinned fetch of the official EPA zip; the `fgr15` module parses the
  seven scenario tables (nothing vendored).
- Tokamak fusion neutron sources (`plasma-source`): ring/point D-D and D-T
  sources with ion-temperature-broadened spectra (Brysk 1973; Ballabio 1998),
  seeded sampling to particle vectors, and MCNP SDEF + Serpent source-card
  emission with drift reports (MCPL projection stays caller-side).
- Parametric tokamak plasma source (`plasma-source`): Miller-geometry flux
  surfaces (Fausser 2012) with caller-supplied L/H/A-mode profiles,
  reactivity-weighted emission (Bosch–Hale 1992), Miller-Jacobian volume
  gates, and marginal-histogram source cards with joint-correlation drift.
  Fuel is equimolar D-T, pure D-D, or a D/T mixture at the shared ion
  temperature (Eriksson/DRESS two-branch rate); toroidal sectors, per-species
  ion temperatures, the D(d,p)T proton branch, and the pinned T-T
  normalization also landed (sectors scale totals by rotation/2π with a PHI
  marginal; sampler and cards stay neutron-only).
- Damage and gas-production metrics (`damage`): NRT-dpa and arc-dpa
  displacement functions (NRT 1975; Nordlund 2018), He/H appm and He/dpa
  ratios by spectral folding of caller flux with caller response functions,
  and UQ on the folds over caller MVN blocks; the SPECTER report
  (ANL/FPP/TM-197, US-gov PD) is the validation oracle behind an opt-in
  vendored Table VII fallback for callers with no damage-data pipeline
  (ASTM E693/E521 stay designation-only; PKA-spectrum solving stays out).
- Neutron spectrum unfolding (`unfold`): SAND-II iterative spectral
  adjustment (McElroy et al., AFWL-TR-67-41, 1967), STAYSL-class damped
  least-squares (Perey, ORNL/TM-6062, 1977), GRAVEL chi-square-weighted
  adjustment (Matzke, PTB-N-19, 1994), and MAXED maximum-entropy adjustment
  (Reginatto & Goldhagen, Health Phys. 77 (1999) 579 — journal equations
  only, the closed UMG package never touched) of a caller-supplied guess
  spectrum against measured activation rates over a caller-supplied
  response matrix, with per-group relative-change convergence diagnostics
  and a hard `NotConverged` past the explicit iteration cap (MAXED adds a
  caller chi-square target to the acceptance leg); the IRDFF-II response
  pack (v1 foil subset + v2 extension) ships as a runtime hash-pinned download, never vendored.
- TBR and blanket power bookkeeping (`blanket`): raw TBR from caller
  tallies, multiplicative per-port coverage haircuts, blanket energy
  multiplication, the tritium burn rate from pinned constants, and the
  breeding-margin / net-surplus fuel-cycle metrics — pure arithmetic over
  caller transport tallies (no transport solving, no geometry
  optimization, no coupling into `tritium`).
- Equilibrium data readers (`equilib`): classic-netCDF VMEC `wout`
  files (CDF-1/CDF-2 only; netCDF-4/CDF-5 convert facade-side) plus the
  `&INDATA` input-text grammar, with flux-surface Jacobian helpers for
  volume weighting — reads data, never solves equilibria.
- Wall-load mapping in flux coordinates (`equilib`): closed-form per-cell
  accumulation of a caller birth-rate density field onto the caller wall
  surface sharing the same flux coordinates, plus the angle-weighted
  one-field-period total that conserves the discrete births (pure glue —
  no transport, no shadowing).
- Arbitrary-3D birth-rate lattice source (`plasma-source`): the caller
  supplies the full 3D birth distribution as a point cloud with optional
  field-period symmetry (base-sector cloud replicated uniformly, totals
  scaled by the period count), sampled and emitted through the landed
  mixture/spectrum/card machinery.
- Deuterium hot-tail fraction (`plasma-source`): the one pinned
  single-tail-temperature non-Maxwellian shape on a D/T fuel mixture (bulk
  plus hot tail reacting at their own effective temperatures, five
  sub-branches through strengths, sampler, cards, and drift row).
- Closed-form ECRH accessibility kernel (`plasma-source`): cold
  electron-cyclotron resonance, the weakly-relativistic thermal shift, and
  the O1/X1 cut-off densities located by interpolation along a caller
  beamline — the per-port penalty inputs for blanket bookkeeping (no ray
  tracing, no launcher design).
- Coil fast-fluence / lifetime bookkeeping (`damage`): fast-flux sums above
  a caller threshold, fluence as flux × time, piecewise-constant history
  accumulation, and weakest-link life over caller-supplied limit tables
  (dpa weighting reuses the folds).
- OpenMC statepoint tally bridge (`nucleide.openmc`, pure Python — no Rust
  changes): mesh/cell tally extraction through the caller-side OpenMC API
  into the plain arrays the damage folds accept, plus a CSV interchange
  for the same dict (HDF5 never read in Rust).
- Sublet S1–S7 radiological totals (`alara`): total activity with the IRT
  α/β/γ split and decay heat over caller-supplied decay energies (each
  with the excluding-tritium companion), committed ingestion/inhalation
  hazards over caller 50-year dose coefficients (each with the companion),
  the transport Bq/A₂ ratio with its effective A₂, the IAEA
  clearance-index variant screened at ≤ 1, and the slab/point gamma dose
  over caller groups with the loud 0.3 m clamp — all pinned to the open
  CCFE-PR(16)53 preprint and FISPACT-II open docs.
- He/dpa ratio uncertainty propagation (`damage`): per-draw ratios over the
  seeded joint MVN block with k-SE gates against the bias-corrected
  expectation and the delta-propagated spread, plus the distribution-free
  68% interval gated against the Fieller-construction quantiles (holds at
  honestly-sized blocks); zero-dpa draws fail loudly.
- Per-species ion temperatures on single-fuel parametric configs
  (`plasma-source`): single-fuel D-T reacts at the mass-weighted T_DT,
  single-fuel D-D at T_D, with per-branch Ballabio lines and bit-for-bit
  equal-T recovery.

## Upcoming priorities

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
