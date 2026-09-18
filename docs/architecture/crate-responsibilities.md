---
title: Crate Responsibilities
sidebar:
  order: 2
---

Each Nucleide crate has a single responsibility and a well-defined place in the
workspace dependency graph.

## Foundation crates

### `nucleide-linalg`

Isolation facade over the linear-algebra backend. Today it pulls in `faer`,
`num-complex`, `rand` (OS-entropy features off: seeded ChaCha only, so the
crate still builds for wasm), and `roxmltree` so numeric dependencies do not
leak into other crates. Besides the complex sparse LU core it owns the
UQ-lite sampling kernel: the `sample` module (seeded multivariate-normal
sampling over caller-supplied covariance blocks — Cholesky primary,
eigen-clipping fallback, relative/absolute conventions, SANDY-style
convergence diagnostics) and the `decay` consumer (branch/energy perturbers
preserving the `1 − BR(SF)` deficit, plus `perturb_fission_yields` over
caller-supplied yield blocks — `raw = base * (1 + rel)`, negatives clamped to
zero, rescaled to the incoming block sum; U7 gate in
`validation/uq_lite_vs_sandy.py`). The 0.9.0 `sample` extension adds
`sample_lognormal` plus the closed-form `lognormal_mean`/`lognormal_cov`
helpers. The 0.9.0 `lstsq` module owns the shared dense-real weighted
least-squares kernel. The 0.10.0 `sample` extension adds `sample_lhs`
stratified sampling (theory U7, separate G1/G2 gate).
Other workspace crates depend on `nucleide-linalg`, not on the backend directly.

### `nucleide-nuclei`

Canonical nuclide identification. Owns:

- `NuclideId`: compact `u32` representation compatible with legacy integer ids.
- Element symbol/number tables.
- Name parsing (`U235`, `Am242_m1`, `Ba137m`).
- Naming dialects: MCNP ZAID, Serpent, FLUKA, NIST, CINDER, ALARA, SZA.
- Particle registry and reaction-name registry (labels, MT mapping, hashes).
- Physical data access: AME2020 masses, natural abundances, half-lives, plus
  screening-level `simple_xs`, `scattering_length`, and `decay_energy_mev`
  TSV tables.
- EPA FGR 15 external-dosimetry coefficients: the runtime-download `fgr15`
  module parses the seven `Table_4_*.DAT` scenario tables (fetched and
  hash-pinned by `nucleide.data.fetch_fgr15`; nothing vendored).
- IRDFF-II dosimetry-response pack (v1 foil subset + v2 extension): the
  runtime-download `irdff` module parses
  the `MF=3` sections of the 34 named reactions from `IRDFF-II.g725` text
  (fetched and hash-pinned by `nucleide.data.fetch_irdff`; nothing
  vendored) into caller-ready response rows over the SAND-II 725-group
  structure for spectrum unfolding.

## Capability crates

### `nucleide-material`

Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium
loading, and materials XML export. Radioanalytics (activity, decay heat) take
injected `DecayProvider`/`DecayEnergyProvider` traits; `DecayEnergies` reads
the screening-level `decay_energy.tsv` table. Depends on `nucleide-nuclei`.

### `nucleide-mcnp-io`

MCNP-family file I/O: `xsdir`, `meshtal`, SSW/SURFSRC, PTRAC, WWINP, MCTAL
readers; material extraction from input decks; mesh-to-geometry deck
generation; the typed legacy `SDEF` fixed-source reader (discrete
`SI`/`SP`/`SB` subset, canonical re-emission round-tripping the
spectroscopy decay-source emitter). Depends on `nucleide-nuclei`.

### `nucleide-mcpl-io`

MCPL particle-list interchange read/write (format versions 2 and 3,
single/double precision, gzip-transparent paths), the neutron/gamma-only
SSW↔MCPL conversion v1 (`ssw` module: explicit per-track surface/kind
parameters, reference-header cloning for the return leg, named errors past
n/γ), and the particle-list utilities (`merge_mcpl` first-file header wins
plus provenance comment, `stat:sum` never synthesized, precision promotes
to double on mixed merge; `extract_mcpl` range/predicate subset with the
header preserved verbatim; `mcpl_stats` counts/energy moments/PDG
histogram; `repair_mcpl` paper-pinned count repair for interrupted writes).
Depends on `nucleide-mcnp-io` for the SSW header/track types.

### `nucleide-csg-xlate`

Scoped MCNP CSG translation (v3 scope: surfaces, cells, simple nested
universes, and rectangular `LAT=1` lattices with a full matrix `FILL`) to
OpenMC `geometry.xml`, Serpent input, PHITS input, and GDML (Geant4
geometry, schema version 3.1.7 pinned from the Geant4 `v11.4.2` tag), with a
per-cell/per-surface drift report in the `emit` drift pattern. Maps axis
planes, spheres, on-axis cylinders, `SPH`, `RPP`, and axis-aligned `RCC`
(half-space expansion in the OpenMC direction; native spellings in the
Serpent/PHITS directions; named boolean solids in GDML, where infinite
half-spaces are bounded by a per-deck cutoff recorded as drift); cones,
quadrics, tori, general planes, other macrobodies, complements beyond flat
intersections, reflecting/periodic boundaries (GDML has no boundary
spelling), transforms, hexagonal lattices, matrix fills without `LAT=1`,
tallies, and sources are loud errors. Depends on `nucleide-mcnp-io` only
among workspace crates; never the reverse, never on bindings.

### `nucleide-serpent-io`

Parsers for Serpent MATLAB-style output files (`_res.m`, `_dep.m`, `_det.m`).
No internal crate dependencies.

### `nucleide-fluka-io`

FLUKA interface: USRBIN tally reader and MATERIAL/COMPOUND card generation.
Depends on `nucleide-nuclei`.

### `nucleide-alara-io`

ALARA activation-code interop: input-deck, group-flux, material/element/WDR
library, activation-output, photon-source, and schedule-expansion glue, plus
clearance / waste-classification analytics (clearance index and the
sum-of-fractions rule over parsed inventories, with the EU 2013/59/Euratom
Annex VII Table A vendored as the default limit table and the Spanish CSN
conditional NORM landfill tables as opt-in tables) and Sublet S1–S7
radiological totals (total activity with the IRT α/β/γ split, decay heat
over caller-supplied decay energies, committed ingestion/inhalation hazards,
the transport Bq/A₂ ratio, the IAEA clearance-index variant, and the
slab/point gamma dose — all over caller inventories with caller-supplied
coefficients, never vendored). Depends
on `nucleide-nuclei` only among workspace crates; the solver stays inside ALARA.

### `nucleide-enrichment`

Multicomponent enrichment cascades and SWU analytics. Depends on `nucleide-nuclei` and
stays independent of `nucleide-material`.

### `nucleide-depletion`

CRAM matrix exponential (orders 16 and 48), depletion-chain XML parsing, and
multi-step `Integrator::{Predictor, Cecm, Cf4}` time series with
activity/decay-heat output. Depends on `nucleide-linalg` and `nucleide-nuclei`.

### `nucleide-kinetics`

Prescribed-reactivity point-kinetics transients: one PKE solve over
caller-supplied delayed-neutron data (1+ groups, constant/step/impulse/ramp/
polyline insertions in Δk, equilibrium initials) through an implicit
θ-method integrator with exact schedule-knot stepping, plus the inhour
relation with stable-period solve and the prompt-jump approximation. No
transport, no thermal feedback, no flux coupling, no tabulated data. No
internal crate dependencies.

### `nucleide-spectroscopy`

Gamma-ray spectroscopy and measurement: rectangular/five-point smoothing,
gross/background/net peak counting, quadratic energy bins, log-polynomial
energy/efficiency calibration (including a weighted least-squares
coefficient fit through the shared `linalg::lstsq` kernel), X-ray line
algebra over caller-supplied constants, caller-line SDEF decay-source cards,
and the dollar/plain `.spe` text readers with upstream quirks pinned. Peak
search/fit, activities, and plotting stay out. Depends on
`nucleide-linalg` and `nucleide-nuclei`.

### `nucleide-unfold`

Neutron spectrum unfolding from activation-type measurements: the SAND-II
iterative spectral adjustment (McElroy et al., AFWL-TR-67-41, 1967 — US
government work, clean-room from the report), the STAYSL-class damped
least-squares adjustment (Perey, ORNL/TM-6062, 1977) on the shared
`linalg::lstsq` kernel, the GRAVEL chi-square-weighted adjustment
(Matzke, PTB-N-19, 1994), and the MAXED maximum-entropy adjustment
(Reginatto & Goldhagen, Health Phys. 77 (1999) 579 — journal equations
only, the closed UMG package never touched) — one method per cycle, each
an iterator over the caller-supplied response matrix, measured rates, and
guess spectrum, with per-group relative-change convergence diagnostics and
a hard `NotConverged` past the explicit iteration cap (never a silent
partial spectrum; MAXED additionally accepts at or below a caller
chi-square target, defaulting to the detector count). Every response value
and group bound is caller-supplied: IRDFF and other IAEA-copyright
libraries are never vendored (the IRDFF-II pack ships as a runtime
hash-pinned download, parsed into caller-ready rows). Depends on
`nucleide-linalg` only
among workspace crates; bindings depend on it, never the reverse.

### `nucleide-vr-tools`

MAGIC weight-window generation, OpenMC/Serpent weight-window emission over
`MagicOutput` (pinned `<weight_windows>`/WWINP spellings; Serpent text
re-parses through the `nucleide-mcnp-io` WWINP reader), mesh source sampling
with alias tables, and Gaussian KDE source sampling (`KdeSampler`). Depends
on `nucleide-mcnp-io` (`nucleide-nuclei` comes in transitively).

### `nucleide-cccc-io`

CCCC text-subset readers (ISOTXS multigroup libraries, RTFLUX/ATFLUX/RZFLUX
flux files) plus a minimal PARTISN deck writer with ISOTXS nuclide mapping.
Depends on `nucleide-nuclei` only among workspace crates; the solver stays out of
scope.

### `nucleide-fispact-io`

FISPACT-II inventory-output parser producing ALARA-compatible response
frames (`nucleide-alara-io` `ResponseFrame` rows), plus the clearance-bearing
wide inventory table printed with the `HAZARDS` + `CLEAR` keywords. Depends on `nucleide-alara-io` and `nucleide-nuclei`;
activation solving stays inside FISPACT-II.

### `nucleide-origen-io`

Scoped ORIGEN 2.2 TAPE readers: `TAPE5` input echo, `TAPE6` output
inventories, `TAPE9`-style decay constants. Depends on `nucleide-nuclei` only among
workspace crates; burnup driving stays inside ORIGEN.

### `nucleide-r2s`

Scoped rigorous two-step (R2S) workflow builder: zone-to-flux linking from
ALARA decks, schedule expansion, uniform-split photon-source assembly,
WATTS-class sweep expansion (`expand_sweep`), and facility-flow accounting
(`snapshot_inventory`) for differencing facility snapshots. Also owns the
versionless ARMI DB-snapshot → deck adapter (`snapshot.rs`):
dict-in only, no HDF5, no ARMI layout versioning, opaque zone ids,
volumes-method decks with per-zone mixtures carrying atoms/barn-cm number
densities. Depends on `nucleide-alara-io`, `nucleide-mcnp-io`,
`nucleide-vr-tools`, `nucleide-material`, `nucleide-depletion`, and
`nucleide-nuclei`; transport and activation solving stay inside their respective
codes. The uniform split preserves only the total shutdown strength until
group-wise emission data is wired through per nuclide.

### `nucleide-tritium`

1D tritium diffusion-trapping kernel: Fickian mobile transport (T1) coupled
to N extrinsic McNabb–Foster trap species (T2) on a slab with
caller-supplied temperature, solved by cell-centred finite volumes with
implicit theta-stepping through the shared `linalg::tridiag` Thomas solver.
Owns the Dirichlet/Sieverts/Henry/zero-flux surface taxonomy plus
recombination ends (`J = K_r c²`, closed in steady state and transient by
the face-response construction, G5/G6), and multi-layer series stacks with
Sieverts or Henry internal interface conditions (`c/K` and flux continuous;
the linear interface flux folds into the tridiagonal step matrix —
G7/G8/G9) plus vented-sink recombination gaps (single face concentration
with the `K_r x²` desorption jump, closed in closed form on the CUT matrix
— G10/G11). The flux-continuous product law has no spelling by
construction; multi-D/FEM,
heat coupling, and vendored property tables stay out. Depends on
`nucleide-linalg` only among workspace crates; bindings depend on it, never
the reverse.

### `nucleide-plasma-source`

Tokamak fusion-neutron source creation. Ring and point sources over the D-D
(2.45 MeV) and D-T (14.1 MeV) reactions with ion-temperature-broadened
Gaussian spectra (Brysk 1973; Ballabio et al. 1998 coefficients), plus a
parametric Miller-geometry plasma: caller-supplied L/H/A-mode density and
temperature profiles (Fausser et al. 2012) over closed-form flux surfaces,
reactivity-weighted emission (Bosch & Hale 1992), and a seeded sampler to
particle vectors (position, direction, energy, weight). Fuel is equimolar
D-T, pure D-D, a D/T mixture at the shared profile ion temperature
(`S = n²·[f_D·f_T·⟨σv⟩_DT + (f_D²/2)·⟨σv⟩_DD]`, the Eriksson/DRESS rate
rule with exact equimolar/pure recovery anchors), at per-species ion
temperatures (D-T at the mass-weighted `T_DT`, D-D at `T_D`), or with the
one pinned deuterium hot-tail fraction (bulk plus tail sub-pairs at their
own effective temperatures, five sub-branches through strengths, sampler,
cards, and drift row; `eta = 0` recovers exactly). Toroidal sectors
restrict births with totals scaled by `rotation/2π` (full rotation
recovers the full torus bit-for-bit; partial sectors add a PHI marginal).
An arbitrary-3D birth-rate lattice drops axisymmetry for
stellarator-class callers (caller point cloud with field-period symmetry
reduction, same machinery per node). The D(d,p)T proton *rate* is
accounted alongside the neutron source (sampler and cards stay
neutron-only); the T-T neutron term is pinned as the third mixture term
but not transported. Closed-form ECRH accessibility (cold resonance,
relativistic shift, O1/X1 cut-offs, beamline locators) reports
heating-access positions for blanket bookkeeping — no ray tracing, no
launcher design. MCNP `SDEF` +
Serpent `src` card emission with a drift report: ring/point cards
round-trip byte-identically through the typed `nucleide-mcnp-io` reader
(whose accepted subset carries the ring's `AXS`/`RAD`/`EXT` keywords);
parametric cards carry the radial/vertical/energy marginals as histograms
and the drift report quantifies the tabulation truncation and the
joint-correlation distance a product-form card cannot carry; Serpent rows
are analytic by design. Profiles are caller inputs — nothing computes them.
Out of scope (loud `NotYetSupported`): T-T neutron transport, proton
transport, reactant distributions beyond the one pinned tail shape (the
rest of the full Eriksson generalization; no tail spelling on lattice
configs), and per-species temperatures without a fuel mixture.
MCPL projection stays caller-side (`vr-tools` KDE layering rule). Depends on
`nucleide-mcnp-io` and `nucleide-nuclei`; never on `mcpl-io`, never on
bindings.

### `nucleide-damage`

Damage and gas-production metrics by spectral folding: the NRT-dpa
displacement function (Norgett–Robinson–Torrens 1975, modified
Kinchin–Pease with the Lindhard damage-energy partition via the Robinson
fit) and the arc-dpa efficiency correction (Nordlund et al. 2018,
CC BY 4.0) as closed forms over caller-supplied material constants; He/H
production in appm and He/dpa ratios from piecewise-constant-per-group
folds of caller `(flux, response, bounds)` slices; UQ on the folds over
caller MVN blocks through `linalg::sample` (pinned-seed, k-SE gates),
including per-draw He/dpa ratio UQ over the joint flux/He/dpa block with
the bias-corrected expectation gate (zero-dpa draws fail loudly).
Zero-flux groups contribute exactly 0; the He/dpa ratio at zero dpa is a
named error, never `inf`. Caller-supplied response functions stay the core;
the one opt-in exception is the vendored SPECTER Table VII fallback
(spectrum-averaged damage-energy cross sections for 24 elements plus the
Table II `E_d` column, displacement XS only — never consulted implicitly).
ASTM E693/E521 are designation-only, SPECTER (ANL/FPP/TM-197, US-gov PD) is
the validation oracle behind that fallback, and PKA-spectra solving stays
out (the fispact-org PKA evaluator is GPL-3.0, never read). Coil
fast-fluence / lifetime bookkeeping lives here too (fast-flux sums above a
caller threshold, history accumulation, weakest-link life over
caller-supplied limit tables; dpa weighting reuses the folds) —
magnetics/quench/structural analysis and vendored limit tables stay out.
Depends on `nucleide-linalg` and `nucleide-nuclei`; bindings depend on it,
never the reverse.

### `nucleide-emit`

Single-material emission to MCNP/Serpent/FLUKA/ALARA/PARTISN cards plus a
mass-drift report. Pure glue: depends on `nucleide-material` and the five
`*-io` crates, never the reverse, and never on bindings.

### `nucleide-blanket`

TBR and blanket power bookkeeping over caller transport tallies: the raw
TBR ratio from caller `(bred, source)` pairs, multiplicative per-port
coverage haircuts, blanket energy multiplication, the tritium burn rate
from the pinned `17.6` MeV / `3.0160492` u constants, and the
breeding-margin / net-surplus fuel-cycle metrics. Pure arithmetic —
Stellaris Point-A values (raw TBR `1.1070`, `1.074` after the 3% ECRH-port
haircut, multiplication `1.20`, burn `416.6` g/day) gate the kernel as
regression vectors. Transport solving, TBR target solving, and any
coupling into `tritium` stay out (that crate keeps its
no-breeding-coupling scope line; this companion owns the breeding-side
arithmetic). No workspace dependencies; bindings depend on it, never the
reverse.

### `nucleide-equilib-io`

Equilibrium data readers (never a solver): a pure classic-netCDF `wout`
reader (CDF-1/CDF-2 accepted per magic probe; HDF5-backed netCDF-4 and
CDF-5 rejected loudly as `WrongVariant` for facade-side conversion) over
the documented variable lists (dims `radius`/`mn_mode`/`mn_mode_nyq`;
`nfp`, `ns`, `xm`, `xn`, `rmnc`, `zmns`, `lmns`, `gmnc` minimum plus
later-use `bmnc`/`bsubumnc`/`bsubvmnc`/`bsubsmns`/`currumnc`/`currvmnc`
and optional `mpol`/`ntor`/`phiedge`/`volume_p`), the VMEC `&INDATA`
namelist text grammar (power-series profiles, `NCURR = 1` as the `I'(s)`
profile, fixed-boundary stance, scalar indices only), and flux-surface
Jacobian helpers (J1 point sum, J2 tensor grid over one field period,
J3 surface average) for volume weighting, plus closed-form wall-load
mapping in flux coordinates (per-cell `q[j][k]` accumulation plus the
discrete-conserving angle-weighted total — no transport, no shadowing).
Synthetic fixtures are built
from the published variable lists only — never solver output; the
on-disk variant of DESC-saved files and VMEC++-era outputs stays out
until re-probed. No workspace dependencies; bindings depend on it,
never the reverse.

## Binding crates

### `nucleide-bindings`

PyO3 extension module that exposes workspace crates to Python as
`nucleide._internal`. It is the only crate allowed to depend on `pyo3` and the
only crate allowed to know about the Python API surface. It depends on most
capability crates.

### `nucleide-wasm`

`wasm-bindgen` crate that exposes a subset of the workspace to the browser for
the interactive tutorials on the docs site. It is not published to crates.io
(the release workflow's publish list excludes it); the site builds it with
`wasm-pack`.

## Dependency rules

- Workspace crates may depend on other workspace crates.
- Workspace crates must not depend on `nucleide-bindings`, `nucleide-wasm`, or `pyo3`.
- `nucleide-bindings` may depend on workspace crates.
- `nucleide-enrichment` must not depend on `nucleide-material`.

## Release order

When publishing to crates.io, publish in dependency order:

1. `nucleide-linalg`
2. `nucleide-nuclei`
3. `nucleide-material`
4. `nucleide-mcnp-io`
5. `nucleide-csg-xlate`
6. `nucleide-mcpl-io`
7. `nucleide-serpent-io`
8. `nucleide-fluka-io`
9. `nucleide-enrichment`
10. `nucleide-depletion`
11. `nucleide-kinetics`
12. `nucleide-spectroscopy`
13. `nucleide-unfold`
14. `nucleide-tritium`
15. `nucleide-plasma-source`
16. `nucleide-damage`
17. `nucleide-blanket`
18. `nucleide-equilib-io`
19. `nucleide-vr-tools`
20. `nucleide-alara-io`
21. `nucleide-cccc-io`
22. `nucleide-fispact-io`
23. `nucleide-origen-io`
24. `nucleide-r2s`
25. `nucleide-emit`
26. `nucleide-bindings`

(`nucleide-wasm` is cdylib-only and never published; keep this list in
sync with the publish list in `.github/workflows/release.yml`.)
