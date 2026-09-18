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

## [0.16.0] - 2026-09-18

### Added

- Per-species ion temperatures on single-fuel parametric configs in
  `nucleide-plasma-source` (exposed through the `species_temperatures`
  dict without a `fuel` dict): single-fuel D-T reacts at the landed
  mass-weighted `T_DT = T_D + (2/5)·(T_T − T_D)`, single-fuel D-D at `T_D`,
  with the Ballabio lines following per branch (D-T at `T_DT`, D-D at
  `T_D`). Equal pair temperatures `(T, T)` reproduce the
  shared-temperature single-fuel kernel bit-for-bit (strengths, sampled
  stream, cards). Anchors: `(T_D, T_T) = (20, 30)` keV gives `T_DT = 24`
  keV with `⟨σv⟩_DT = 5.414667327922193e-22` m³/s and the D-T line at mean
  14.077934372230146 MeV, sigma 0.3700390508551786 MeV. Mixture paths are
  untouched; the `None`-pair path keeps the landed expression trees
  verbatim (no new RNG draws); the deuterium tail still requires a fuel
  mixture.

- Sublet S4 ingestion and S5 inhalation hazard kernels in
  `nucleide-alara-io` (one shared `HAZARDS` shape): each folds per-nuclide
  activities with the caller-supplied 50-year committed dose coefficients
  (`e_ing` / `e_inh`, Sv/Bq, never vendored) into the `TOTAL ... HAZARD
  FOR ALL MATERIALS` dose (Sv) plus the excluding-tritium companion, with
  no split key. Pure arithmetic over caller inventories; negative /
  non-finite inputs and overflow are loud errors. Exposed as
  `nucleide.alara.alara_ingestion_hazard` / `alara_inhalation_hazard`
  (entry keys `nuclide`, `activity_bq`, `e_ing_sv_per_bq` /
  `e_inh_sv_per_bq`), with hand-vector gates in `tests/test_alara.py` and
  C9 rows in `validation/clearance_vs_pypact.py`.

- Sublet S3 gamma dose-rate kernel in `nucleide-alara-io` (new `dose`
  module): the slab dose `C·B/2·Σ μa/μm·Sγ` (`B = 2`) and the point dose
  `C·Σ μa/(4πr²)·e^(−μr)·mₛ·Sγ` (Sv/h, `C = 3.6e9·|e|`, `Sγ = I·A`) over
  caller gamma groups — specific activity, group yields, and air/mixture
  attenuation all caller-supplied, never vendored (`mixture_mu` folds
  elemental values with mass fractions). Point distances below 0.3 m clamp
  to 0.3 m and report it loudly via the `clamped` flag, never silently.
  Exposed as `nucleide.alara.alara_dose_slab` / `alara_dose_point` /
  `alara_dose_mixture_mu`, with hand-vector gates in `tests/test_alara.py`
  and C8 rows in `validation/clearance_vs_pypact.py`.

- Sublet S6 transport-ratio and S7 IAEA clearance-index kernels in
  `nucleide-alara-io`: S6 folds per-nuclide activities with the
  caller-supplied `A2` limits (TBq, never vendored) into the dimensionless
  `Total Bq/A2` ratio plus the effective A2 it defines
  (`nucleide.alara.alara_transport_ratio`, entry keys `nuclide`,
  `activity_bq`, `a2_tbq`); S7 folds activities with the caller-supplied
  IAEA levels (Bq/kg, never vendored) and the total mass into the
  dimensionless clearance index with its `<= 1` screening class (boundary
  included) and dominant contributor
  (`nucleide.alara.alara_iaea_clearance_index`, entry keys `nuclide`,
  `activity_bq`, `limit_bq_per_kg`). Pure arithmetic; bad inputs and
  overflow are loud errors. Hand-vector gates in `tests/test_alara.py` and
  C10 rows in `validation/clearance_vs_pypact.py`.

- He/dpa ratio uncertainty propagation in `nucleide-damage`
  (`he_dpa_ratio_uq`, exposed as `nucleide.damage.he_dpa_ratio_uq`):
  per-draw ratios over the seeded joint `[flux, he, dpa]` MVN block (the
  landed `linalg` engine, no new sampling machinery), reporting the draw
  mean ± draw standard deviation with `k`-SE gates against the
  second-order bias-corrected expectation and the first-order
  delta-propagated standard deviation, plus the distribution-free 68%
  interval (draw 16th/50th/84th percentiles) gated against the
  Fieller-construction quantiles — the all-regime verdict that holds at
  honestly-sized blocks where the symmetric-spread gate cannot. Any draw
  at non-positive dpa fails loudly with the landed `ZeroDpa` vocabulary,
  never `inf`/`NaN` with a spread. The single-response
  `fold_uq("he_dpa_ratio")` spelling stays a loud error naming the
  two-response function. Gates in `tests/test_damage.py` and G6–G7 rows in
  `validation/damage_vs_specter.py`.

- Runnable workflow notebooks in `notebooks/` (first-wall damage + UQ,
  activation screening + radiological totals, tokamak source sampling +
  card emission): synthetic inputs and committed fixtures only, fixed
  seeds, no outputs vendored in the files, executed end to end in CI
  (`tests/test_notebooks.py`, needs `nbclient`/`nbformat`/`ipykernel`),
  with Colab badges for one-click runs once the PyPI release is live.
- Tutorial code-execution gate (`tests/test_tutorial_code.py`): every
  runnable ```python block in `docs/tutorials/python/*.md` executes in CI
  in document order, so the tutorials double as verified examples.
  Placeholder-path snippets now point at real fixtures (or carry an
  explicit skip marker for user-side/h5py/scratch-file demos), and the
  `fetch_compendium` tutorial snippet reads the committed compendium file
  instead of downloading.

- Interactive demos for the new surfaces (new tabbed radiological-totals
  demo covering S3–S7, He/dpa ratio-UQ block in the damage demo,
  single-fuel species-pair inputs in the fusion-source demo, all backed by
  new thin WASM endpoints), plus theory-page coverage of the new equation
  sets, tutorial prose, and citations — including verified references for
  the previously uncited equilibrium (VMEC, DESC, simsopt) and blanket
  (Stellaris concept paper) theory pages.
- Tutorial sidebar reorganization: one workflow-grouped reading order
  shared by the sidebar, the index pages, and the docs index (duplicate
  and missing position values fixed).

### Fixed

- Target-less `<decay>` chain entries now parse as decay out of the modeled
  chain (pure diagonal loss, no gain term) — the decay analogue of the
  long-standing target-less `<reaction>` convention, including
  `target="nothing"`. Previously the committed
  `fixtures/depletion/chain_ni.xml` was rejected outright even though its
  target-less decays (Fe55→Mn55 and friends) are physically just truncated
  chains. An explicit `target` naming a nuclide absent from the chain stays
  a loud `UnknownNuclide` error; CRAM and Bateman agree exactly on the
  pure-loss exponential.

## [0.15.0] - 2026-09-18

### Added

- OpenMC statepoint tally bridge (`nucleide.openmc`, pure Python — no Rust
  changes). Statepoint files are HDF5 and are never read in Rust; the
  facade drives the caller-side OpenMC Python API (call sequence pinned
  against OpenMC 0.16.0, the validation-container build:
  `StatePoint(path, autolink=False)` with an explicit `close()`, then
  `tally.get_reshaped_data(value=..., expand_dims=True)` for `"mean"`
  and `"std_dev"`) and returns plain nested lists the landed folds
  already accept: `tally_arrays` / `read_statepoint` yield `[voxel][group]`
  `flux` plus matching `rel_err` (`std_dev / mean`, 0.0 where the mean is
  0) and the `G+1` MeV energy edges (`[]` unbinned — a single total
  group), with `cell_flux` / `voxel_flux` selecting the group vector the
  `damage` folds take alongside the edges as `bounds`. `tally_to_csv` /
  `tally_from_csv` move the same dict through a commented CSV with stock
  stdlib only, so an OpenMC machine exports and any machine reads back
  bit-identically. Mesh tallies are `RegularMesh`-only (C-order voxels,
  centimetre bounds); a missing file, missing OpenMC, missing
  tally/score, multi-nuclide tallies, multi-bin non-energy filters, a
  second mesh/energy filter, and non-regular meshes are all loud errors
  naming the remedy — never silent empties. OpenMC stays a caller-side
  lazy import (the base package never depends on it). Gates are
  synthetic-statepoint unit tests in `tests/test_openmc_bridge.py`
  (pinned fetch sequence, mesh C-order and unit conversion, singleton
  cell-filter transparency, CSV round-trips, loud-error vectors, and the
  `cell_flux` → `nrt_dpa` feed check).

- Wall-load mapping in flux coordinates in `nucleide-equilib-io` (exposed
  through `nucleide.equilib` as `wall_load`): closed-form per-cell
  accumulation `q[j][k] = Σ_i S·J·Δs` of a caller birth-rate density
  field onto the caller wall surface sharing the same flux coordinates
  (radial-major voxels, one-field-period toroidal convention, Jacobian
  voxels from the landed `fourier_jacobian` kernel) plus the
  angle-weighted one-field-period total, which conserves the discrete
  births to roundoff. Pure glue — no transport, no shadowing, no
  FEM/thermal, no CAD, no HDF5; shape mismatches and non-finite or
  negative densities/Jacobians are loud named errors reusing the crate
  vocabulary. Gates are hand-computed toy geometries at 1e-12 relative
  (axisymmetric uniform field vs the analytic `S·J0` wall flux with the
  `12π²` total, a cosine-Jacobian vector, radial weighting, and discrete
  conservation) in `wall.rs` unit tests and `tests/test_equilib.py`, with
  always-run E9–E10 rows in the same `equilib` validation report.

- Closed-form ECRH accessibility kernel in `nucleide-plasma-source`
  (exposed through `nucleide.plasma_source` as `ecrh_scalars` /
  `ecrh_accessibility`): cold electron-cyclotron resonance
  (`B_res = 2π·m_e·f/(n·e)`), the weakly-relativistic Maxwell–Jüttner
  shift (`B_res(T_e) = γ·B_res`, `γ = 1 + T_e/511 keV`), and the O1/X1
  cut-off densities (`n_c,O1 = ε₀·m_e·(2π·f)²/e²`,
  `n_c,X1 = n_c,O1·(f − f_ce)/f`, evanescent `None` where `f ≤ f_ce`)
  from the published ECRH/dispersion facts (Bornatici et al., Nucl. Fusion
  23 (1983) 1153; Stix, Waves in Plasmas (1992) Ch. 1), located by linear
  interpolation along a caller beamline (`s` [m] strictly increasing,
  `B` [T], `n_e` [m⁻³]) into plain-float resonance/cut-off positions —
  the per-port penalty inputs for blanket bookkeeping. Pure functions
  with named loud errors (empty/ragged/non-monotonic beamlines, negative
  densities, bad frequency/harmonic/temperature); no ray tracing, no
  launcher design, no HDF5. Gates are hand-computed vectors at 1e-9
  relative (170 GHz ITER-class resonance/shift/cut-offs plus beamline
  crossings) in `ecrh.rs` unit tests and `tests/test_ecrh_access.py`.

- Sublet S1+S2 radiological totals in `nucleide-alara-io` (exposed through
  `nucleide.alara` as `alara_total_activity` / `alara_decay_heat`). S1 sums
  per-nuclide activities `Ai = Ni λi` (Bq) with the FISPACT-II IRT
  α/β/γ split (IRT 4 → alpha; IRT 1, 2, 11, 14, 16, 17, 19, 20 → beta;
  IRT 3 → gamma; IRT 12, 13 split α/β and IRT 15 splits α/γ by
  caller-supplied branch fractions); S2 folds caller-supplied average decay
  energies into per-class decay heat `Ai·E·C1` (kW, C1 = eV→kJ). Both carry
  the excluding-tritium companion and conserve parts-to-total exactly;
  unmapped IRTs and bad values are loud named errors. Pinned to the open
  CCFE-PR(16)53 preprint (Table X, Tables VI–VII) and the open FISPACT-II
  `output_interpretation` reference — never the paywalled journal pages;
  decay energies stay caller-supplied, never vendored. S3–S7, ICRP/IAEA
  coefficients, bremsstrahlung, DPA/KERMA/gas, and fission counts stay out.
  `validation/clearance_vs_pypact.py` gains always-run C7 gates
  (hand-computed vectors at exact equality) in the same `clearance` report.

- TBR and blanket power bookkeeping (`nucleide-blanket`,
  `nucleide.blanket`): raw TBR from caller `(bred, source)` tallies,
  multiplicative per-port coverage haircuts, blanket energy
  multiplication, the tritium burn rate from pinned `17.6` MeV /
  `3.0160492` u constants, and the breeding-margin / net-surplus
  fuel-cycle metrics. Pure arithmetic over caller transport tallies —
  no transport solving, no geometry optimization, and no coupling into
  the `tritium` permeation kernel (which keeps its no-breeding-coupling
  scope line). Stellaris Point-A values gate the kernel as regression
  vectors: raw TBR `1.1070`, `1.074` after the 3% ECRH-port haircut,
  multiplication `1.20`, burn `416.6` g/day.

- Coil fast-fluence / lifetime analytics in `nucleide-damage` (exposed
  through `nucleide.damage` as `coil_fast_flux` / `coil_fast_fluence` /
  `coil_accumulate` / `coil_lifetime` / `coil_remaining` plus
  `FPY_SECONDS`): magnet lifetime bookkeeping over caller spectra — the
  fast-flux sum above a caller threshold (groups count iff their upper edge
  clears the threshold; a threshold cutting through a group includes the
  whole group), fluence as flux × time, piecewise-constant history
  accumulation, and weakest-link life (`min` of limit/rate, or
  limit-minus-accumulated over rate clamped at zero). dpa spectral
  weighting stays in the landed folds (the NRT/arc fold at one second is
  the rate the life kernel consumes). Zero rates never breach and are
  skipped (all zero is infinite life with a `None` limiting channel); all
  shapes, signs, and orderings are loud named errors reusing the fold
  vocabulary. Limit tables are always caller-supplied — published design
  numbers (3e22 / 1.5e23 m⁻² fast-fluence limits at a 99th-percentile coil
  fast flux of 9.5e13 m⁻²s⁻¹ → ~10 FPY, limited by the first channel) are
  gates at pinned tolerance, never defaults. Magnetics, quench,
  structural analysis, and transport solving stay out.

- Equilibrium data readers in `nucleide-equilib-io` (exposed through
  `nucleide.equilib` as `probe_variant` / `read_wout` / `jacobian` /
  `jacobian_grid` / `parse_indata` plus typed `indata_*` helpers). A
  pure classic-netCDF `wout` reader accepts CDF-1 (`43 44 46 01`) and
  CDF-2 (`43 44 46 02`) per magic probe (`ncdump -k` must report
  `classic` / `64-bit offset`) and rejects HDF5-backed netCDF-4 and
  CDF-5 loudly for facade-side conversion (no HDF5 in Rust); it reads
  the documented variable lists (`nfp`, `ns`, `xm`, `xn`, `rmnc`,
  `zmns`, `lmns`, `gmnc` minimum with dims
  `radius`/`mn_mode`/`mn_mode_nyq`, later-use `bmnc`/`bsubumnc`/
  `bsubvmnc`/`bsubsmns`/`currumnc`/`currvmnc`, optional
  `mpol`/`ntor`/`phiedge`/`volume_p`). The `&INDATA` namelist grammar
  covers scalar switches, power-series profiles, and boundary Fourier
  tables (sliced indices like `RBC(0:4,2)` are a loud error;
  `NCURR = 1` reads as the `I'(s)` profile; absent `LFREEB` reads
  fixed boundary). Flux-surface Jacobian helpers evaluate the `gmnc`
  Fourier sum at points and over one-field-period grids for volume
  weighting. Reads data, never solves equilibria; synthetic fixtures
  are built from the published variable lists only (never solver
  output). Gates are magic-byte accept/reject vectors, Jacobian hand
  vectors, and INDATA grammar vectors in-crate plus
  `tests/test_equilib.py` and the `validation/equilib_vs_desc.py`
  oracle (DESC/simsopt legs SKIP loudly outside the container).

- Arbitrary-3D / stellarator birth-rate lattice source in
  `nucleide-plasma-source` (exposed through `nucleide.plasma_source`
  with `kind="lattice"`): the caller supplies the full 3D birth
  distribution as a point list — `points` of `{"position": [x, y, z]
  [cm], "rate": w >= 0, "ion_temperature_kev": Ti}` (relative
  birth-rate weights, arbitrary global scale; `Ti = 0` is the
  monoenergetic nominal line) — with the landed machinery applied per
  node (D/T mixture rates, Ballabio spectra, seeded sampling, card
  emission with drift reports; MCPL projection stays caller-side).
  The optional `field_periods`/`base_angle` pair declares field-period
  symmetry (both keys or neither): the points are the base-sector cloud,
  replicated uniformly around the machine axis, with the symmetric total
  scaled by the period count. A two-point hand vector pins the contract
  (total 3, mean birth position `(100, 0, 25)`, mono 14.021 MeV line); a
  single node reproduces the landed point source bit-for-bit; a uniform
  ring cloud converges to the analytic ring moments; the folded
  base-sector stream reproduces the expanded full-cloud stream bit-for-bit
  (including the single-node case: one base node fans out over all
  copies); cards carry the cloud's cylindrical-`R`/vertical/energy
  marginals with the same drift rows plus a `lattice discretization` row
  and still round-trip byte-identically through the typed SDEF reader.
  Non-finite positions, negative rates, empty or zero-total clouds, and
  non-symmetric fold inputs are loud named errors (`InvalidLattice` /
  `InvalidSymmetry`). The Stellaris-class precedent (30 radial × 50
  poloidal × 100 toroidal point sources, 14.06 MeV monoenergetic) is the
  named design point; equilibrium solving, transport solving, CAD/DAGMC,
  tabulated plasma data, and HDF5 stay out. `validation/
  plasma_source_vs_openmc.py` adds always-run P12 gates (hand vectors,
  single-point recovery, ring-limit closure, fold replication and totals,
  card round trip, mixture branch fire) and a container-only O12 probe on
  a dense monoenergetic ring cloud vs the upstream ring geometry and
  Ballabio helpers.

- Deuterium hot-tail fraction for the D/T fuel mixture in
  `nucleide-plasma-source` (exposed through `nucleide.plasma_source`
  as the `deuterium_tail={"fraction": eta, "temperature_kev": T_tail}`
  dict [keV] on a parametric `fuel` blend): the one pinned
  single-tail-temperature non-Maxwellian shape within the Eriksson et
  al., Comput. Phys. Commun. 199 (2016) 40 arbitrary-distribution
  framework — the bulk `(1 − eta)` deuterium at `T_D` plus a hot tail
  `eta` at `T_tail`, each sub-pair reacting at its own mass-weighted
  relative temperature (D-T at `T_DT`/`T_DTt`, D-D at `T_D`/`T_mix`/
  `T_tail`) with the Ballabio lines following per sub-branch through the
  strength rule, the sampler's branch roulette, the card energy
  marginals, and a `deuterium tail` drift row carrying the tail neutron
  share. The pinned 70/30 hand vector (`T_D = 20`, `T_T = 30` keV,
  `eta = 0.05` at `T_tail = 60` keV) gates the rate at 1e-12 relative;
  `eta = 0` reproduces the no-tail mixture kernel bit-for-bit (strengths,
  sampled stream, cards — regression gates); the tail requires a `fuel`
  mixture and composes with species temperatures and toroidal sectors.
  Non-finite, out-of-range, or negative tail parameters are loud named
  errors (`NonFinite` / `InvalidTail` / `NegativeIonTemperature`). The
  effective-temperature line convention carries the same error stance as
  the species pair (the residual against the full Eriksson numerical
  integration is out of scope); the rest of the full Eriksson
  generalization stays a loud boundary. `validation/
  plasma_source_vs_openmc.py` adds always-run P13 gates (sub-rate-rule
  moments and branch share, the drift tail-neutron share, zero-fraction
  recovery, loud errors) and a container-only O13 probe against the
  upstream quadrature with the NeSST reactivities evaluated at the five
  sub-pair temperatures.
- Interactive demos for the new surfaces: the fusion-source demo gains a
  parametric + hot-tail tab, the activation demo a Sublet S1+S2 tab, plus
  new ECRH, blanket, equilibrium, coil (damage), and wall-load demos.
  Fixed the WASM `sampleFusionSource` parametric path, which rejected the
  `kind`/`n`/`seed` sampling envelope its own dispatcher requires
  (`FusionParametricSpecJson` denies unknown fields to keep stray
  `fuel_mixture`/`species_temperatures`/`sector` keys loud — only the
  envelope keys are stripped before parsing).

## [0.14.0] - 2026-09-17

### Added

- Vented-sink recombination internal interfaces for multi-layer tritium
  stacks in `nucleide-tritium` (exposed through `nucleide.tritium`):
  `LayerStack::new` now accepts `Interface::recombination(rate)` alongside
  the linear Sieverts/Henry laws. A vented gap carries a single face
  concentration `x ≥ 0` with half-cell fluxes `J_L = g_L·(c_L − x)`,
  `J_R = g_R·(x − c_R)` (`g = 2D/dx` per side, no solubility involved) and
  the desorption sink `K_r·x²` venting from the face, so
  `J_L − J_R = K_r·x²`. The bulk matrix is CUT at the gap
  (block-diagonal); a lone gap closes in closed form
  `x = 2P/(√(Q²+4K_rP)+Q)`, coupled faces (several gaps, or recombination
  outer ends alongside) by Newton iteration with the analytic Jacobian —
  same FV θ-stepper and `linalg::tridiag` solve, with the per-step close
  fused into the trap Picard loop exactly like the G6 outer ends. The
  `K_r → 0` limit is a continuous-concentration joint (not a Sieverts
  law); a non-positive permeation drive is a loud error in the steady
  state (the transient follows the G6 clamping stance per step). The
  flux-continuous product law has no spelling by construction (spurious
  insulated root, symmetry breaking — permanently rejected). Gates
  G10a–G10g/G11a–G11d (in-crate unit tests, synthetic stacks with
  hand-derived closed forms, no `fixtures/tritium/` changes): 2-layer
  closed form at `1e-12` with gap residual and desorption-carrying outer
  balance, `K_r` limit pins, a 4-layer mixed-law stack (one interface of
  each law) with flux continuity across both linear gaps, coupled-Newton
  residuals at the Newton contract, traps on the owning layer's isotherm,
  transient θ-balance with the desorption sink to roundoff, `dt`-halving
  order bands, late-time asymptote, and the trap-coupled Picard-fusion
  balance. The Python `steady_layers`/`transient_layers` facades accept a
  `{"kind": "recombination", "rate": Kr}` dict per gap (bare
  `"recombination"` strings fail loudly — the rate has nowhere to go) and
  report `interface_faces` (face value per gap, `None` at linear gaps;
  `[time][gap]` in transients) for closing the discrete balance with the
  desorption term; single-slab and default-stack signatures are unchanged.

- IRDFF-II dosimetry-response registry expansion (`nucleide-nuclei`,
  `nucleide._internal.parse_irdff_g725`): the `V2_REACTIONS` extension adds
  the 26 further named `MF=3` dosimetry sections of the same pinned
  `IRDFF-II_g725.zip` (same URL + SHA-256 pins, same cache, same loud
  offline/mismatch errors) — the remaining threshold (n,p)/(n,a)/(n,2n)
  monitors, capture foils, fission chambers, and high-threshold bismuth
  monitors, for 34 reactions total. `MF=10` isomer sections stay out, the
  default pack stays v1, and the `unfold` iterators are unchanged. Parser
  review: the pinned file carries `L2=99` (not 0) in three `MF=3` head cards
  (notably v1's `56Fe(n,p)`), now accepted as 0-or-99 and loudly pinned.
  Gates: synthetic mocked-fetch unit tests plus a container live-fetch gate
  (loud SKIP when offline) with hand-checked first-response groups and
  `unfold` round-trips at pinned tolerance. Nothing IAEA-copyright is
  vendored.

- MAXED maximum-entropy unfolding (`nucleide-unfold`, `nucleide.unfold.maxed`):
  the fourth and last adjustment family in the crate — Shannon relative
  entropy against the guess-as-prior over the exponential Lagrange family,
  driven by an equilibrated Levenberg-Marquardt dual step with a chi-square
  target (`target_chi2`, default: the detector count) and a per-iteration
  trust-region cap. Caller sigmas are the `1/σ²` weights, zero measurements
  carry zero weight (the GRAVEL house rule), and a relatively converged run
  whose chi-square still exceeds the target is a hard `NotConverged` (never
  a partial spectrum). Validation adds always-run U19-U23 gates (determined
  recovery, counting-statistics underdetermined round-trip, bracketed
  IRDFF-II analytical shape probes, SAND-II cross fixed points, contract
  errors) to the same `unfold` report.

- D(d,p)T proton bookkeeping and the pinned T-T branch normalization for
  the parametric tokamak plasma source (`nucleide-plasma-source`,
  `nucleide.plasma_source`): `ParametricPlasmaConfig::proton_strength_density`
  / `total_proton_strength` (plus `FuelMixture::tt_neutron_coefficient`)
  report the D(d,p)T proton production rate alongside the neutron source
  under the pinned 50/50 convention — the proton density IS the D-D
  neutron-branch density bit-for-bit (pure D-D gives exactly one proton per
  neutron; D-T-only gives exactly none) — while the sampler draws no proton
  energies, the cards carry no proton distributions, and proton transport
  stays out of scope. The T-T neutron term is pinned as the third mixture
  term `S = n²·[f_D·f_T·⟨σv⟩_DT + (f_D²/2)·⟨σv⟩_DD + f_T²·⟨σv⟩_TT]` (neutron
  multiplicity 2 folded in, reacting at `T_T`, matching the upstream
  `openmc-plasma-source` rule) but is NOT transported: no Bosch–Hale T-T
  fit exists (Bosch & Hale 1992 covers only D(d,n), D(d,p), D-T, D-³He),
  no Ballabio-class T-T line exists (three-body continuum), and the
  oracle's vendored tables stay out of scope — a pure-tritium mixture stays
  a loud error. A mixture with `f_T = 0` reproduces the landed two-branch
  kernel bit-for-bit (regression gates). New `proton_accounting(spec)`
  facade (per-neutron ratio plus sector-aware totals for parametric specs;
  `None` totals for ring/point specs, which carry no density model).
  `validation/plasma_source_vs_openmc.py` adds always-run P9 gates (exact
  pure-fuel ratios, the 70/30 share vs quadrature, the no-T recovery
  anchor, the loud pure-T error) and a container-only O9 probe on a
  tritium-rich 10/90 blend showing the sampled stream sits on the
  two-branch upstream quadrature while the three-branch reference pulls
  away by exactly the T-T-predicted scale (protons have no oracle —
  upstream models neutrons only).

- Toroidal sectors for the parametric tokamak plasma source
  (`nucleide-plasma-source`, `nucleide.plasma_source`): the additive
  `start_angle`/`rotation_angle` pair [rad] (`ToroidalSector`) restricts
  births to `[start_angle, start_angle + rotation_angle)` with uniform birth
  angles and emission totals scaled by `rotation_angle/2π` (half torus ×0.5
  exactly, quarter torus ×0.25 exactly); the exact full-rotation spelling
  reproduces the landed full-torus kernel bit-for-bit (same seeded stream,
  same cards — regression gates). A partial sector adds a uniform
  angle-bin marginal to the cards (`PHI=D4` on SDEF, `phi d4` on Serpent)
  plus a `toroidal sector` drift row; the SDEF card still round-trips
  byte-identically through the typed reader (the `mcnp-io` SDEF subset
  gains the `PHI` keyword). Non-finite or out-of-range angles are loud
  errors (`NonFinite`/`InvalidSector`; both keys or neither in the Python
  spec). `validation/plasma_source_vs_openmc.py` adds always-run P10 gates
  (in-sector uniform births, full-rotation stream/card recovery, the PHI
  marginal round-trip with its drift row) and a container-only O10 probe
  showing a partial sector leaves the upstream-quadrature (r, z, E) moments
  untouched.

- Per-species ion temperatures for the parametric tokamak plasma source
  (`nucleide-plasma-source`, `nucleide.plasma_source`): the additive
  `species_temperatures={"D": T_D, "T": T_T}` dict [keV]
  (`SpeciesIonTemperatures`) reacts a D/T fuel mixture at distinct
  Maxwellian species temperatures — the D-T branch at the mass-weighted
  relative temperature `T_DT = T_D + (2/5)·(T_T − T_D)`, the D-D branch at
  `T_D` (Eriksson et al. 2016, the full distinct-Maxwellian
  generalization) — through the emission strength, the sampler's branch
  roulette, the card energy marginals, and the axis spectrum summary. The
  D-T Ballabio line at `T_DT` is an effective-temperature convention
  (Table III was fitted to single-temperature Maxwellians; the residual
  against the full Eriksson numerical integration is out of scope, and
  non-Maxwellian reactants stay a loud boundary). A uniform `T_D = T_T`
  pair reproduces the shared-temperature mixture kernel bit-for-bit
  (strengths, sampled stream, cards — regression gates). The pair
  requires a `fuel` mixture (a single-fuel config with the pair set is a
  loud error); the profile ion temperature is then unused for rate and
  spectrum (still validated), while the density profile keeps shaping the
  emission. Non-finite or negative temperatures are loud errors
  (`NonFinite`/`NegativeIonTemperature`). `validation/
  plasma_source_vs_openmc.py` adds always-run P11 gates (pair-temperature
  sampled moments and branch share vs quadrature, equal-T stream/card
  recovery) and a container-only O11 probe against the upstream quadrature
  with the NeSST reactivities evaluated at the pair temperatures.

## [0.13.0] - 2026-09-16

### Added

- D/T fuel mixtures for the parametric tokamak plasma source
  (`nucleide-plasma-source`, `nucleide.plasma_source`): a `fuel` dict
  `{"D": f_D, "T": f_T}` on the parametric source spec (the
  openmc-plasma-source spelling; then `reaction` is optional and unused)
  switches the emission model to the Eriksson et al., Comput. Phys. Commun.
  199 (2016) 40 mixture rate at the shared profile ion temperature,
  `S = n²·[f_D·f_T·⟨σv⟩_DT + (f_D²/2)·⟨σv⟩_DD]` — the standard
  `n_a·n_b·⟨σv⟩` rate with the `1/(1+δ_ab)` same-species guard, the DRESS
  rate construction upstream `openmc-plasma-source` implements. Birth
  energies run a per-particle branch roulette between the two Ballabio
  lines; emitted cards carry the same product-form marginals, and the
  magnetic-axis spectrum summary reports the two-branch Gaussian-mixture
  moments. The equimolar (`f_D = f_T = 1/2`) and pure-D-D (`f_D = 1`)
  limits reproduce the landed single-fuel kernels bit-for-bit (regression
  gates), and non-finite, negative, or non-summing fractions are loud named
  errors. Toroidal sectors stay out of scope. `validation/plasma_source_vs_openmc.py`
  adds always-run analytic gates P8 (pinned 70/30 blend: sampled moments vs
  in-script mixture quadrature, D-D branch share) and a container-only O8
  cross-check against the upstream map/profiles with the NeSST D-T and D-D
  reactivities.

- Spanish CSN conditional-clearance tables for NORM landfill disposal in
  `nucleide-alara-io` (exposed through `nucleide.alara` as
  `alara_clearance_es_table`). The Consejo de Seguridad Nuclear draft
  technical opinion CSN/PDT/AICD/TGE/2503/02 (TGE/VAR/2025/1, hosted on
  csn.es, accessed 2026-09-16 — RD 1029/2022 and RD 1217/2024 (RINR), both
  transposing Directive 2013/59/Euratom) publishes Tables 1-3, per-landfill
  type (inert / non-hazardous / hazardous) clearance levels in Bq/g for each
  NORM material nature (rocks / ashes / sands / slags / oil & gas), keyed by
  the Table 4 natural decay-chain keys. All three tables are vendored as
  attributed embedded TSVs in the EU-table pattern (`try_` + cache split,
  row-count-pinned tests); the Table 4 chain keys expand to per-member
  entries at the parent value through the `nuclei` dialect machinery (EU
  Part-2 precedent), with the last claiming key supplying a shared member's
  level, as the source lists each secular-equilibrium chain before its
  subchains. The caller selects the table and material column explicitly —
  no cross-table logic, no "most permissive wins"; the EU Annex VII table
  stays the vendored default. Screening arithmetic only, never a compliance
  decision. `validation/clearance_vs_pypact.py` gains always-run C6 gates
  (transcription spots, chain-expansion members, hand-computed screening
  vector at exact equality).

- Henry internal interfaces for multi-layer tritium stacks in
  `nucleide-tritium` (exposed through `nucleide.tritium`): `LayerStack::new`
  now accepts `Interface::Henry` alongside `Interface::Sieverts`. The
  linear Henry law (`c/K_H` continuous, flux continuous) shares the
  Sieverts resistance form — face flux `J = (c_i/K_i − c_{i+1}/K_{i+1})/R_f`
  with `R_f` the two half-cell resistances `dx/(2·D·K)` and the layer
  `solubility` playing `K_H` — so it folds into the same tridiagonal
  θ-step matrix with no new solver machinery; recombination internal
  interfaces stay loud `UnsupportedInterface` errors (recorded limitation,
  deferred to its own spike cycle). Gates G9a/G9b (in-crate unit tests,
  synthetic stacks with hand-derived series-resistance closed forms, no
  `fixtures/tritium/` changes): a 2-layer all-Henry stack against the
  series-resistance closed form at `1e-12` relative on the fluxes (`1e-18`
  absolute floor) with interface flux continuity and half-cell-corrected
  interface potentials to roundoff, and a 3-layer mixed Sieverts/Henry
  stack holding flux continuity across both interface laws against the
  same closed form; the G7c one-layer recovery anchor is unchanged. The
  Python `steady_layers`/`transient_layers` facades gain an additive
  `interfaces` keyword (`"sieverts"` per gap by default, `"henry"`
  supported, `"recombination"` fails with a clear error); single-slab and
  default-stack signatures are unchanged.

- Vendored SPECTER Table VII fallback in `nucleide-damage` (exposed
  through `nucleide.damage` as `specter_table` / `specter_damage_energy` /
  `specter_ed` / `specter_spectra`): spectrum-averaged damage-energy cross
  sections in keV-b for 24 elements (Be through Pb) across the report's 7
  spectra (thermal / fission / 14 MeV / HFIR / EBR-II / FFTF / fusion),
  transcribed from Greenwood & Smither, ANL/FPP/TM-197 (January 1985,
  US-government public domain — OSTI `purl/6022143`, the same hash-pinned
  PDF the validation oracle downloads) by reading the scan pages directly,
  plus the Table II `E_d` column (eV) per element. Displacement cross
  sections in barns follow the report's own `0.8/2E_d` rule with the
  vendored `E_d` (or a caller-chosen value). Scope is displacement XS only
  (no gas columns); the group-wise Appendix A printouts are image-only
  pages in the scan and stay untranscribed, so the table covers the 24
  Table VII elements rather than the 41 Table I entries. The table is an
  opt-in caller input only — the fold kernels are unchanged and never
  consult it implicitly; the `try_` + cache constructor split, row-count
  pin (24), and spot tests follow the EU-table pattern. The HFIR column
  reproduces the oracle's Table VI Fe/Ti/Cu spots at Table VII print
  precision (3 significant figures, 5e-3 relative); the oracle script
  itself is untouched and keeps passing.

- STAYSL-class damped least-squares spectral adjustment in
  `nucleide-unfold` (exposed as `nucleide.unfold.staysl`, second method in
  the one-method-per-cycle order after SAND-II): the Perey ORNL/TM-6062
  (1977) anchored-damped objective
  `φ† = argmin_φ Σ_i (c_i(φ) − N_i)²/σ_i² + λ·Σ_j (φ_j − φ_guess,j)²`,
  re-solved each cycle as one weighted least-squares system through the
  shared `linalg::lstsq` kernel (no second least-squares implementation)
  on the current active set, with non-positive groups clipped to zero and
  pinned NNLS-style. Sigmas are caller-supplied, one finite positive value
  per detector (`w_i = 1/σ_i²`); zero, non-finite, or missing sigmas are
  loud named errors. Damping `λ` is absolute Tikhonov damping on the
  distance to the guess; convergence is the shared contract (per-group
  relative-change tolerance with an explicit iteration cap; exhausting the
  cap is a hard `NotConverged`, reusing the landed variant). Gates are
  synthetic forward-fold-then-recover round-trips at pinned tolerances
  (determined recovery, sigma-weighted precise-reading dominance,
  underdetermined rate reproduction, IRDFF-II analytical benchmark-field
  shapes) plus the cross-method fixed-point check (a SAND-II fixed point
  is a least-squares fixed point) and the named-contract errors; the
  always-run `validation/unfold_vs_analytic.py` oracle gains gates U9-U13
  in the same `unfold` report (no `SECTION_ORDER` change).

- GRAVEL chi-square-weighted spectral adjustment in `nucleide-unfold`
  (exposed as `nucleide.unfold.gravel`, third method in the
  one-method-per-cycle order after SAND-II and STAYSL-class): the Matzke
  PTB-N-19 (1994) multiplicative iteration
  `φ_j ← φ_j·exp(Σ_i W_ij·ln(N_i/c_i) / Σ_i W_ij)` with
  `W_ij = (R_ij·φ_j / c_i)·(N_i²/σ_i²)` — the SAND-II base rate share
  times the caller-sigma measurement-weight factor — as an iterator with
  the same diagnostics shape (`spectrum`/`rates`/`rate_factors`/
  `iterations`/`tolerance`/`max_rel_change`). Sigmas are caller-supplied,
  one finite positive value per detector; zero, non-finite, or missing
  sigmas are loud named errors, never silent uniform weighting.
  Convergence is the shared contract (per-group relative-change tolerance
  with an explicit iteration cap; exhausting the cap is a hard
  `NotConverged`, reusing the landed variant). Pinned divergences from
  SAND-II: zero measurements carry zero weight and are simply not fitted
  (they pin nothing to zero); at a fixed point every rate ratio is 1, so
  the weights are irrelevant and SAND-II, STAYSL-class, and GRAVEL fixed
  points coincide where the formulations agree (cross-checked both
  directions). Gates are synthetic forward-fold-then-recover round-trips
  at pinned tolerances (determined recovery, counting-statistics-sigma
  underdetermined rate reproduction plus error reduction) plus the
  IRDFF-II analytical benchmark-field shapes (Trkov et al., Nucl. Data
  Sheets 163 (2020) 1 — published facts used with citation, not vendored
  data); the always-run `validation/unfold_vs_analytic.py` oracle gains
  gates U14-U18 in the same `unfold` report (no `SECTION_ORDER` change).
  MAXED stays recorded for a later cycle.

- IRDFF-II runtime-download response pack in `nucleide-nuclei`
  (`nucleide.data.fetch_irdff` + `nucleide._internal.parse_irdff_g725`):
  the official IAEA `IRDFF-II_g725.zip` distribution is fetched at runtime
  and hash-pinned (SHA-256 recorded in `nucleide.data`, FGR-15 precedent —
  nothing IAEA-copyright is vendored), and the `IRDFF-II.g725` member text
  is parsed into caller-ready response rows over the SAND-II 725-group
  structure (`groups` 726 eV bounds plus one cross-section vector per
  named reaction) for the v1 foil-activation subset (`au197_ng`,
  `in115_ng`, `u235_nf`, `u238_nf`, `fe56_np`, `ni58_np`, `al27_na`,
  `na23_n2n`; MF=3 sections only — the pointwise MF=10 isomer reactions
  stay out). The rows feed `unfold.sandii`/`staysl`/`gravel` unchanged;
  structural problems, missing/duplicate sections, off-structure
  energies, and row-shape violations are loud errors. Offline or
  hash-mismatch fetches are loud `RuntimeError`s naming the URL or both
  digests.

## [0.12.0] - 2026-09-16

### Added

- Clearance / waste-classification analytics in
  `nucleide-alara-io` and `nucleide-fispact-io` (exposed through
  `nucleide.alara` and `nucleide.fispact`). The clearance index
  CI = Σ Aᵢ/CLᵢ and the sum-of-fractions screening rule (RS-G-1.7 §5,
  referenced by designation) run as pure arithmetic over parsed
  `(nuclide, activity)` inventories — Sublet et al., NDS 139 (2017) 77 is
  the methodology reference; further indices from that paper stay
  recorded, not implemented. Limit tables are caller-supplied, with the EU
  2013/59/Euratom Annex VII Table A vendored as the attributed default
  (official legal text transcribed from EUR-Lex CELEX:32013L0059, OJ L 13,
  17.1.2014, accessed 2026-09-15; reusable per Decision (EU) 2011/833;
  Bq/g == kBq/kg; IAEA tables never vendored). Nuclides key through the
  `nuclei` dialect machinery; missing limits and bad values are loud named
  errors; the boundary CI == 1 classifies as satisfied ("does not
  exceed"). Screening arithmetic only — never a compliance decision.
  `fispact-io` gains the clearance-block reader for the real FISPACT-II
  wide inventory table printed with the `HAZARDS` + `CLEAR` keywords
  (grammar cross-checked against the Apache-2.0 `fispact/workshops`
  reference outputs): per-step `(interval, time_s, cooling)` rows with
  activity, per-nuclide clearance index, flags, and half-life (`Stable` →
  `-1`). Validation adds `validation/clearance_vs_pypact.py`: always-run
  analytic gates C1-C5 (exact CI vectors, boundary probes both sides, EU
  table spots, fixture parse, end-to-end screening) plus a pypact
  (Apache-2.0) cross-check of the overlapping inventory columns with loud
  SKIP outside the oracle environment.

- GDML (Geant4) emission direction for `nucleide-csg-xlate`: scoped MCNP
  CSG decks now translate to GDML documents (`deck_csg_to_gdml` in Rust,
  `nucleide.mcnp.parse_csg_to_gdml` / `read_csg_to_gdml` on the Python side)
  alongside the OpenMC/Serpent/PHITS directions, over the same v3 scope
  (surfaces, cells, nested universes, rectangular `LAT=1` lattices,
  material stub). Cells become named boolean solids inside a `world`
  volume, universes become assemblies, lattices expand to one placement per
  element, and infinite half-spaces are bounded by a per-deck cutoff
  recorded in the drift report (new actions `halfspace-bounded`,
  `lattice-expanded`, `material-stub`). Emitted documents pin the GDML
  schema version 3.1.7 (Geant4 `v11.4.2` tag); the validation harness
  validates them against the runtime-fetched, SHA-256-pinned XSD wherever
  `lxml` is available. Reflecting and periodic boundaries stay loud errors
  in this direction (Geant4 expresses boundaries through wrapper code, not
  geometry markup). New `fixtures/mcnp/inp/deck_csg_cylinders.txt` extends
  the `deck_csg_*` family.
- Parametric tokamak plasma source in `nucleide-plasma-source`
  (exposed through `nucleide.plasma_source` with
  `kind="parametric"`): Miller-geometry flux surfaces (Fausser et al., Fus.
  Eng. Des. 87 (2012) 787 — `major_radius`/`minor_radius`/`elongation`/
  `triangularity`/`shafranov_factor`) carrying caller-supplied L/H/A-mode
  density and temperature profiles (`mode`, `pedestal_radius`, and the
  `ion_density_*`/`ion_temperature_*` parameter sets; profiles are inputs,
  never computed), with reactivity-weighted emission (Bosch & Hale, Nucl.
  Fusion 32 (1992) 611, Atzeni–Meyer-ter-Vehn parametrization; new
  `nucleide.plasma_source.reactivity` exposes ⟨σv⟩ in m³/s). Birth
  positions are drawn ∝ `S·R·|J|` — closed-form Miller Jacobian with the
  toroidal `R` factor — and energies from the local-ion-temperature
  Ballabio Gaussian, all through the seeded sampler. Source cards for
  parametric sources carry the radial/vertical/energy marginals as
  discrete histograms (`RAD`/`EXT`/`ERG` distributions on `SDEF`, `rad`/
  `ext`/`erg` on Serpent `src`) and the drift report quantifies the
  tabulation truncation plus the joint-correlation distance a product-form
  card cannot carry; SDEF cards still round-trip byte-identically through
  the typed reader. Analytic gates: `J == r` torus limit, `J == κ·r·(1 −
  2·esh·r·cosθ/a²)`, area `π·κ·a²`, volume `2π²·κ·a²·R₀` (any Shafranov
  shift), finite-difference Jacobian agreement at general triangularity,
  profile pedestal continuity, reactivity transcription goldens, and
  sampled-moment gates vs fine quadrature. The validation oracle adds
  always-run P5–P7 and container-only O4–O7 legs against
  openmc-plasma-source's Miller map, Fausser profile functions, and NeSST
  reactivities (loud SKIP outside the container). Loud boundary: fuel
  mixtures (Eriksson-weighted reactant distributions), toroidal sectors,
  and the D(d,p)T proton branch stay named `NotYetSupported` errors.
- Tokamak fusion neutron sources in the new `nucleide-plasma-source`
  crate (exposed as `nucleide.plasma_source`): ring and point sources over
  the D-D (2.45 MeV) and D-T (14.1 MeV) reactions with ion-temperature-
  broadened Gaussian spectra (Brysk, Plasma Phys. 15 (1973) 611; Ballabio
  et al., Nucl. Fusion 38 (1998) 1723 Table III coefficients), a seeded
  sampler to particle vectors (position/direction/energy/weight arrays,
  deterministic under a pinned seed), and MCNP `SDEF` + Serpent `src`
  card emission with a drift report (`rel_drift` of the tabulated
  spectrum, `reparsed` flags). Emitted SDEF cards round-trip
  byte-identically through the typed `nucleide-mcnp-io` reader, whose
  accepted subset now also carries the ring keywords `AXS`/`RAD`/`EXT`;
  Serpent rows are analytic by design (no Serpent source reader). MCPL
  projection stays caller-side (`vr-tools` KDE layering rule): the crate
  has no `mcpl-io` dependency. Parametric Miller-geometry plasma profiles,
  mixed-fuel spectra, toroidal sectors, and the D(d,p)T proton branch are
  loud named `NotYetSupported` errors — the second landing's scope.
  New-crate checklist applied (workspace/release wiring, `gen-reference`
  entries, crate-responsibilities section, README feature table); the
  two-part `validation/plasma_source_vs_openmc.py` oracle adds always-run
  analytic gates (ring/spectrum closed-form moments, determinism, card
  round trips) plus container-only moment cross-checks against the
  upstream MIT `openmc-plasma-source` package (loud SKIP outside the
  container; the `Containerfile` gains the pinned oracle layer).

- OpenMC/Serpent weight-window emission in `nucleide-vr-tools` (exposed as
  `nucleide.vr.emit_openmc_weight_windows` / `emit_serpent_wwin`, taking the
  same `MeshTally` + `MagicOutput` pair as `nucleide.vr.magic`): MAGIC lower
  bounds are formatted into an OpenMC `settings.xml` fragment (`<mesh>` +
  `<weight_windows>` — spelling pinned against the public OpenMC docs
  §3.29/§3.66 with the sub-element layout the C++ reader parses; energies
  converted MeV → eV; flat bounds emitted x-fastest per energy group) or into
  the MCNP WWINP text spelling that Serpent reads via
  `wwin <name> wf "<file>" 2` (user guide §2.2.8.3 — the only `wwin` FMT
  whose layout is publicly pinned; the Serpent-native FMT=1 layout is never
  emitted). OpenMC upper bounds are synthesized as `5 × lower` and drift
  notes report synthesized or null-cell structure; ill-formed inputs
  (negative/non-finite windows, unsorted energy or mesh bounds, OpenMC
  tuning parameters outside the reader's ranges) are loud named errors.
  Golden-text fixtures and re-parse assertions cover both emitters — the
  Serpent file round-trips through the existing `nucleide-mcnp-io` WWINP
  reader — and always-run structural probes extend the existing `magic`
  validation report (container OpenMC load cross-check; Serpent load is a
  loud SKIP, proprietary and not installed).

- Damage and gas-production metrics in the new `nucleide-damage` crate
  (exposed as `nucleide.damage`): NRT-dpa (Norgett, Robinson & Torrens,
  Nucl. Eng. Des. 33 (1975) 50–54 — modified Kinchin–Pease with the
  Lindhard damage-energy partition via the Robinson fit) and arc-dpa
  (Nordlund et al., Nat. Commun. 9 (2018) 1084, CC BY 4.0 — piecewise
  efficiency ξ(T_d)) as closed forms over caller-supplied material
  constants (`E_d`, `b_arc`, `c_arc`; no fitted table is vendored); He/H
  production in appm and He/dpa ratios by piecewise-constant-per-group
  spectral folding of caller `(flux, response, bounds)` slices
  (`nrt_dpa` / `arc_dpa` / `gas_appm` / `he_dpa_ratio`). Zero-flux groups
  contribute exactly 0; negative inputs, non-increasing bounds, and the
  He/dpa ratio at zero dpa are loud named errors — never `inf` (the
  `ZeroMaxFlux` precedent). `fold_uq` propagates caller-block uncertainty
  through the landed `linalg` MVN machinery (pinned seed, `k`-standard-
  error gates against the exact bilinear expectation and the first-order
  propagated standard deviation; UQ on the ratio is a named-open).
  SPECTER (Greenwood & Smither, ANL/FPP/TM-197, US-gov PD) is the
  validation oracle only — `validation/damage_vs_specter.py` folds
  report-transcribed HFIR-CTR32 spots (Fe/Ti/Cu dpa; C12/Li7/B10/N14 gas)
  against the runtime hash-pinned PDF at report print precision, loud
  SKIP outside the container; no SPECTER table is vendored and ASTM
  E693/E521 stay designation-only. PKA-spectra solving stays out (the
  fispact-org PKA evaluator is GPL-3.0, never read). New-crate checklist
  applied (workspace/release wiring, `gen-reference` entries,
  crate-responsibilities section, README feature table); the
  `Containerfile` gains `poppler-utils` for the oracle's PDF text
  extraction.
- MCPL particle-list utilities in `nucleide-mcpl-io` (exposed as
  `nucleide.mcpl.merge_mcpl` / `extract_mcpl` / `mcpl_stats` /
  `repair_mcpl`, all file-based and `.gz`-transparent): `merge_mcpl`
  concatenates compatible files with the first-file header plus a
  provenance comment (`stat:sum` comments ride along verbatim; sums are
  never synthesized or updated; mixed single/double input promotes to
  double, the lossless direction; any other header-option disagreement is
  a loud error); `extract_mcpl` writes an index-range (`start`/`stop`) or
  predicate subset with the source header preserved verbatim;
  `mcpl_stats` returns record counts, energy moments, total weight, and a
  per-PDG-code histogram; `repair_mcpl` implements the paper-pinned
  `mcpl_repair` semantics, recomputing the particle count of a file whose
  writer never patched it and ignoring a partially written trailing
  record. Synthetic golden fixtures in `fixtures/mcpl/utils/` pin merge
  concat order, precision promotion, and extract subsets; new always-run
  U1-U9 gates in `validation/mcpl_vs_refs.py` cover the same contracts.

- Neutron spectrum unfolding in the new `nucleide-unfold` crate (exposed as
  `nucleide.unfold`): the SAND-II iterative spectral adjustment (McElroy et
  al., AFWL-TR-67-41, 1967 — US government work, clean-room from the report)
  as an iterator over a caller-supplied response matrix, measured
  activation rates, and guess spectrum, with per-group relative-change
  convergence diagnostics (iteration count, final max change, per-detector
  measured/folded rate factors). The guess is adjusted by a positivity-
  preserving weighted geometric mean of per-detector correction factors
  until the fold reproduces the rates; non-convergence past the explicit
  iteration cap is a named hard `NotConverged` error — never a silent
  partial spectrum (the `tritium` face-Newton precedent) — and unreachable
  measurements (zero response rows or rates wiped by pinned-to-zero groups)
  are a named `RatesUnreachable`. Every response value and energy-group
  bound is caller-supplied: IRDFF and other IAEA-copyright libraries are
  never vendored (a runtime-download pack stays a later decision under the
  FGR-15 precedent). STAYSL-class least-squares adjustment (on the shared
  `linalg::lstsq` kernel), GRAVEL, and MAXED are recorded for later
  one-method-per-cycle landings. Gates are synthetic
  forward-fold-then-recover round-trips at pinned tolerances plus the
  IRDFF-II analytical benchmark-field shapes (Trkov et al., Nucl. Data
  Sheets 163 (2020) 1 — published facts used with citation, not vendored
  data); the new always-run `validation/unfold_vs_analytic.py` oracle
  replays them (no external unfolding oracle exists, so no container
  dependency is added). New-crate checklist applied (workspace/release
  wiring, `gen-reference` entries, crate-responsibilities section, README
  feature table).
- Multi-layer tritium permeation in `nucleide-tritium` (exposed
  through `nucleide.tritium` as `steady_layers`/`transient_layers`):
  series stacks of caller-specified layers (thickness, cells, Arrhenius
  `D`, solubility `K_S`, per-layer traps/temperature/source — e.g.
  W/Cu/CuCrZr first walls) with Sieverts internal interface conditions
  (`c/K_S` continuous, flux continuous). The interface fluxes are linear in
  the cell values and fold directly into the tridiagonal θ-step matrix
  (same finite-volume stepper and `linalg::tridiag` Thomas solve), so no
  interface iteration is needed; outer ends reuse the landed surface
  taxonomy including recombination ends through the same face-response
  construction. Henry/recombination internal interface laws are loud
  `UnsupportedInterface` errors in v1 (recorded limitation). A one-layer
  stack dispatches to the landed single-slab kernel and reproduces it
  exactly. Gates are analytic and synthetic, in unit tests with recorded
  provenance (no `fixtures/tritium/` changes): 2/3-layer series-resistance
  closed forms at `1e-12` with interface flux continuity to roundoff,
  property-continuous 2-layer ≡ landed kernel, layered recombination-outlet
  closed form, transient `dt`-halving order bands and late-time asymptote
  (Dirichlet and recombination outlets), and the discrete mass balance to
  roundoff.

### Fixed

- Robustness hardening across the 0.12.0 surface (all loud errors, never
  silent partials or panics): GDML `#n` complement of an empty region is
  `ComplementTooComplex` instead of a panic; the FISPACT-II clearance
  reader rejects non-finite/negative inventory values and step times and
  `total_clearance_index` reports overflow as an error (now returning
  `Result<f64>`); `alara-io` fraction sums report overflow instead of
  `Ok(inf)` and the EU Annex VII loader gains a fallible
  `try_eu_annex_vii` constructor with a row-count-pinned test;
  `mcpl-io` decoding rejects non-finite particle fields and universal
  weights on read (mirroring the write path) so corrupt files fail loudly
  instead of poisoning `mcpl_stats`; parametric sampling uses a panic-free
  total-order CDF search plus a fallible `try_sample`; SAND-II replaces its
  last `expect` with a loud error; layered-tritium `unreachable!` sites
  become named errors.
- Efficiency: `merge_mcpl` pre-sizes from declared counts; layered-tritium
  per-cell layer lookup is a cached-offset binary search; Criterion
  coverage added for every 0.12.0 kernel and reworked emission path (new
  `unfold`, `tritium`, `mcpl`, `plasma` benches; `windows_emit` and
  `wwinp_to_text` cases in the existing `vr-tools`/`mcnp-io` benches).

## [0.11.0] - 2026-09-14

### Added

- EPA FGR 15 external-dosimetry coefficients (EPA 402-R-25-001, July 2025)
  as a runtime-download feature — nothing from the EPA file is vendored.
  `nucleide.data.fetch_fgr15` downloads the official coefficient zip into the
  per-user platform cache (`%LOCALAPPDATA%\nucleide\Cache` on Windows,
  `~/Library/Caches/nucleide` on macOS, `$XDG_CACHE_HOME/nucleide` or
  `~/.cache/nucleide/` elsewhere), hash-pinned to SHA-256
  `71314b3f1d73c197da8b589e74c450f61befce48c66ace6ac1f6e0597560bd91`
  (mismatch fails loudly — EPA may have revised the file). The new
  `nucleide-nuclei` `fgr15` module parses the seven `Table_4_*.DAT` scenario
  tables (1,252 nuclides × 6 age columns; units string recorded from the
  header), and `nucleide.nuclei` gains `fgr15_dose_rate` /
  `load_fgr15_table` / `parse_fgr15_table` / `fgr15_age_index` (FGR 15 name
  spellings `H-3`, `Ba-137m`, `Sb-124n`). FGR 15 is external exposure only
  and complements — never silently overrides — `dose_per_g`
  (HNF-5636/PyNE). `validation/nuclear_data_vs_refs.py` gains an FGR 15
  section (row-count/column gates plus ten hand-transcribed spots at exact
  equality). Screening-level only — not for safety decisions.
- Fission-yield perturbation consumer: `nucleide-linalg`
  `perturb_fission_yields` now perturbs caller-supplied yield blocks
  (`raw = base * (1 + rel)`, negatives clamped to zero, rescaled to the
  incoming block sum — independent blocks sum to 2, cumulative blocks
  higher; zero-base rows stay zero), closing the last named-open hook in
  the UQ kernel (the `FissionYieldsOpen` error variant is removed). The
  thin Python `nucleide.uq.perturb_fission_yields` wrapper starts working
  (signature unchanged); the caller supplies the block and builds `rel`
  (tape `dY/Y` sigmas or correlated draws forwarded from the sampling
  kernel). New U7 validation gate in `validation/uq_lite_vs_sandy.py`
  (theory U8) over a caller-selected U235-thermal subset: perturbed sums
  exact at the pinned seed, draw moments within the IID k-SE bound.
- Scoped MCNP→OpenMC CSG translation v1 (new `nucleide-csg-xlate` crate):
  `deck_csg_to_openmc_xml` maps deck surfaces, cells, and a material stub
  (`material="void"` or the MCNP number) to OpenMC `geometry.xml` with a
  per-cell/per-surface drift report (macrobody expansions, complement
  inlines, reflecting/periodic links, dropped params/cards, each with
  reasons). Axis planes, spheres, on-axis cylinders, `SPH`, `RPP`, and
  axis-aligned `RCC` map; cones, quadrics, tori, general planes, other
  macrobodies, non-flat complements, transforms, lattices, matrix or
  transformed fills, tallies, and sources are loud named errors, while
  simple nested universes (`U=k`, single-universe `FILL n`) map to
  `universe=`/`fill=`. Thin Python
  `nucleide.mcnp.parse_csg_to_openmc` / `read_csg_to_openmc` facades plus
  synthetic `fixtures/mcnp/inp/deck_csg_*.txt` decks and a CSG section in
  `validation/parsers_vs_refs.py` (structural probes always run; the OpenMC
   `Region.from_expression` cross-check runs in the container).
- Scoped MCNP→Serpent CSG translation (same `nucleide-csg-xlate` crate):
  `deck_csg_to_serpent_input` emits `surf`/`cell` cards over the same v2
  scope with Serpent-native simplifications (`RPP` to `cuboid`,
  axis-aligned `RCC` to truncated cylinders, `#n` passed through as the
  native cell complement, empty regions via a synthesized `inf` surface;
  materials as `m<n>`/`void` names for caller-supplied `mat` cards).
  Reflecting/periodic boundaries are a loud `SerpentBoundaryOutOfScope`
  (Serpent `set bc` is global; per-surface mapping unverified). Thin
  Python `nucleide.mcnp.parse_csg_to_serpent` / `read_csg_to_serpent`
  facades plus structural `serpent structure` probe rows in
  `validation/parsers_vs_refs.py`.
- Scoped MCNP→PHITS CSG translation (same `nucleide-csg-xlate` crate):
  `deck_csg_to_phits_input` emits `[Surface]`/`[Cell]` sections over the
  same v2 scope with manual-verified identical symbols (planes, spheres,
  on-axis cylinders, `SPH`/`RPP`/`RCC`-any-orientation/axis-aligned
  `BOX`, coefficients verbatim), native `#` complements, `U=`/`FILL=`
  params, `*` reflective surfaces, verbatim densities, and an outer-void
  `-1` heuristic for void union/complement cells (drift-noted).
  Periodic pointers are a loud `PhitsBoundaryOutOfScope`. Thin Python
  `nucleide.mcnp.parse_csg_to_phits` / `read_csg_to_phits` facades plus
  structural `phits structure` probe rows in
  `validation/parsers_vs_refs.py`.
- Rectangular-lattice (`LAT=1`) CSG translation in all three directions
  (same `nucleide-csg-xlate` crate): a `LAT=1` cell with a full matrix
  `FILL` (every element filled, no transform) now emits a rectangular
  lattice — OpenMC `<lattice type="rectangular">` (`dimension`,
  `lower_left`, `pitch`, `universes`; MCNP `k, j, i` order maps to
  `z`-ascending / `y`-descending / `x`-ascending row-major), a Serpent
  cuboidal `lat` card (type 11) filled via `fill <lattice id>`, and PHITS
  `LAT=1` with a matrix `FILL` (ranges plus the universe list in MCNP
  order verbatim) — with pitch and lower-left derived only from the
  lattice cell's `RPP` or axis-plane-box bounds and lattice ids allocated
  outside both cell and universe id spaces (each noted as
  `lattice-emitted` drift). Hexagonal `LAT=2` lattices, `0`-holes,
  non-`RPP`/axis-plane-bounded lattice cells, single-universe fills of
  lattice type, and fill transforms stay loud named errors. New synthetic
  `fixtures/mcnp/inp/deck_csg_lattice_rect.txt` deck, per-direction
  universe-order tests, and lattice cross-reference probe rows in
  `validation/parsers_vs_refs.py`.
- Tritium-transport analytic-gate spec (no kernel yet):
  `docs/theory/tritium.mdx` pins the T1–T2 equation set and the G1–G5
  gate contract (steady linear, permeation time-lag, single-trap limits,
  Sieverts steady, recombination steady; recombination transient stays
  named-open) plus
   `fixtures/tritium/` oracles for the analytic-gate replay at kernel landing.
- 1D tritium-transport kernel v1 (new `nucleide-tritium` crate): T1 mobile
  diffusion with N extrinsic McNabb–Foster trap species (T2), caller
  supplied Arrhenius data and steady temperature profile (no tables, no
  heat solve), Dirichlet/Sieverts/Henry/zero-flux/recombination surface
  taxonomy (recombination steady state closed by the exact face-response
  construction with G5a/G5b gates; recombination transient stays a
  named-open G5 error), and a cell-centred
  finite-volume theta-stepper (Crank–Nicolson default, backward Euler on
  request) solving through the new shared `nucleide-linalg` `tridiag`
   Thomas-solver module. Analytic-gate replay only (Rust fixture tests
   plus `tests/test_tritium.py` on G1–G5 steady at the pinned tolerances; no
   `validation/*_vs_*.py`, no FESTIM oracle). Thin Python
  `nucleide.tritium` facade, WASM `tritiumBreakthrough` facade with a
  breakthrough-curve interactive tutorial.
- Gaussian KDE source sampling (KDSource-class, clean-room) in
  `nucleide-vr-tools`: `KdeSampler` fits an axis-aligned Gaussian KDE over
  caller particle vectors (Silverman or fixed bandwidths; zero-variance
  dims are loud) with deterministic `draw`/`pdf` (caller randoms, no RNG
  inside), gated by Gaussian recovery, normalization, and determinism
  tests. Thin Python `nucleide.vr.KdeSampler` facade.
- WATTS-class sweep expansion in `nucleide-r2s`: `SweepAxis` validation,
  deterministic cartesian `expand_sweep` to named case bundles, and
  `assemble_results` requiring exact per-case coverage, plus thin Python
  `nucleide.r2s.r2s_expand_sweep` facade.
- Facility-flow accounting in `nucleide-r2s`: `snapshot_inventory` totals
  atoms per ARMI bare name (`N × V × 1e-24` over snapshot zones) for
  differencing facility snapshots, plus thin Python
  `nucleide.r2s.r2s_snapshot_inventory` facade.
- JADE-class library assessment oracle
  (`validation/library_vs_refs.py`, auto-discovered): data hygiene,
  independent-block sums to 2.0, and cumulative coverage over the
  committed fission-yield pack (synthetic thresholds, no tapes).
- Legacy `SDEF` fixed-source reader in `nucleide-mcnp-io` (new `sdef`
  module): typed `SdefCard`/`SdefDist`/`SdefProblem` model over the `POS`,
  `CELL`, `SURF`, `VEC`, `DIR`, `ERG`, `NRM`, `PAR`, `WGT`, `TME` keywords
  (inline literals or `Dn` references) plus the discrete `SIn L`/`SPn D`/
  `SBn D` distribution forms; every other keyword or option letter is a loud
  drift note or a clean `SdefError`, never silently misread. Validation
  covers duplicate cards/keywords, dangling `Dn` references, orphan
  `SPn`/`SBn` cards, entry-count mismatches, and empty tables; canonical
  re-emission round-trips the spectroscopy `sdef_decay_source` card dialect
  byte-identically (80-column wrapping, C++ defaultfloat precision-6).
  `DeckProblem.sdef` typed view (enforced in `validate`) and thin Python
  `nucleide.mcnp.parse_sdef` facade; E9 reader round-trip rows in
  `validation/spectroscopy_vs_pyne.py` plus synthetic `tests/test_sdef.py`
  round-trip and error-case coverage.
- Tritium recombination transient closure in `nucleide-tritium` (the last
  named-open gate in the kernel): `solve` now accepts recombination ends
  (`J = K_r c_m²`) in the transient, closing each implicit θ-step by the
  G5 affine face-response construction reused per step — one base Thomas
  solve plus one sensitivity column per recombination end per `dt`,
  closed-form face for one end, analytic-Jacobian Newton for two, fused
  into the trap-Picard loop and sharing `rtol`/`atol`. The
  `RecombinationOpen` error variant is removed. Since no closed form
  exists for the recombination transient, asymptotic + self-convergence
  gates replace the algebraic style: G6a late-time asymptote to the G5a
  steady flux (1e-6 relative at 60 permeation lags), G6b K_r→∞/K_r→0
  recovery of the Dirichlet/zero-flux transients, G6c discrete mass
  balance with the recombination leak to roundoff, G6d dt-halving
  θ-method order (Crank–Nicolson ≈ 2, backward Euler ≈ 1), G6e
  independent method-of-lines cross-check (node-centred central FD +
  explicit RK4, ghost-node recombination face) of the mid-transient
  trajectory at 0.5/1/2 permeation lags against the production transient,
  trap-free and trap-coupled (1.5e-4 relative on the mobile-profile
  max-norm, ≈5x measured headroom), plus
  positivity under the face clamp. Python `nucleide.tritium.transient`
  facade unchanged (dict shape and signatures stay); the theory page
  gains the G6 section.

## [0.10.0] - 2026-09-14

### Added

- Fission product yields from the ENDF/B-VIII.0 neutron-induced and
  spontaneous fission-yield sublibraries: committed
  `crates/nuclei/src/data/fission_yields.tsv` (151,490 rows over 36 parents
  and 122 incident-energy sets; MF8/MT454 independent + MF8/MT459 cumulative,
  values verbatim — both sublibraries are unchanged carries of ENDF/B-VII.1
  evaluations per the tapes' README.txt, with the Mattera–Sonzogni
  cumulative-yield correction not applied), generated by the new
  `--endf-nfpy`/`--endf-sfpy` flags of `scripts/gen-nuclear-data.py`.
  `nucleide-nuclei` gains `FissionYieldOrigin`/`FissionYieldKind` enums plus
  `fission_yield_table`/`fission_yields`/`default_fission_yields`/
  `fission_yield` lookups (lowest-energy independent set = the OpenMC
  `get_default_fission_yields` depletion convention); thin Python
  `nucleide.nuclei.fission_yields`/`fission_yield` facades. The `depletion`
  chain XML parser falls back to this library (lowest-energy independent
  set, products filtered to chain members) when a fissionable parent has no
  `<neutron_fission_yields>`, raising the existing `BadStructure` error when
  the library lacks the parent. Spot gates in
  `validation/nuclear_data_vs_refs.py` pin the table to tape values, with a
  loud-skip OpenMC cross-check when the nfy tapes are absent.
- Seeded Latin-hypercube sampling (`nucleide-linalg` `sample_lhs` over
  caller-supplied blocks: per-dimension one jittered draw per stratum via
  Fisher–Yates permutation, hand-rolled inverse-normal CDF, shared
  Cholesky/eigen-clip factor path with `x = μ + Bz`; thin Python
  `nucleide.uq.sample_lhs` + WASM `sampleLhs` facades). Ships with its
   separate validation gate U6 (theory U7) (`validation/uq_lite_vs_sandy.py`: G1
  stratification-exact plus G2 moments within the IID bound as an upper
  bound only) over synthetic `fixtures/uq/lhs_2x2.json`. A draw mode, not
  a perturbation convention (`perturb_energies(..., "lhs")` stays an error).
- Interactive UQ demo: MVN/LHS draw-mode toggle calling the WASM
  `sampleLhs` facade with the same (mean, cov, n, seed) inputs and result
  shape, with E2E coverage.
- Interactive deterministic demo: the RTFLUX section now renders a chart-only
  per-point flux profile (one grouped-bar trace per energy group over the
  capped `parseRtflux` values, following the ISOTXS grouped-bar precedent)
  with E2E coverage.
- Interactive deterministic demo: PARTISN writer section rendering and
  validating structured deck dicts in the browser via new WASM
  `partisnRender`/`partisnValidate` facades (exact `title`/`dim`/`zones`
  (`id`/`material`/`isotxs_labels`/`density`)/`source` keys mirroring the
  Python `partisn_render`/`partisn_validate` shape, no camelCase aliases;
  zone labels validated against a pasted ISOTXS library via `parseIsotxs`)
  with render preview, inline validation-error display, and E2E coverage.
- Interactive activation demo: the ORIGEN TAPE6 per-step comparison now takes
  a dynamic list of N pasted snapshots (add/remove, each re-parsed with the
  same `parseOrigenTape6` reader) into a per-nuclide series chart over all
  steps in the existing grouped log-bar bundle, with E2E coverage.
  Multi-snapshot file grammar stays out — pasted snapshots only.
- Interactive activation demo: the R2S snapshot tab hosts an R2S voxels
  section over a drift-accepted `VoxelTags` copy-port in `nucleide-wasm`
  (copy/split zone-total tagging plus `.photonSrc` group select-and-sum over
  `nucleide-alara-io` photon types only — no `nucleide-r2s` dependency, so
  the Rayon chain stays untouched). Zone totals and voxel maps tag in the
  browser (200-voxel demo cap) with a scatter/bar/heatmap/histogram chart
  bundle, with E2E coverage.
- Crate-root error re-exports and `Result` aliases (API-stability pass,
  additive half): `nucleide-linalg` re-exports `DecayError` /
  `LstsqError` / `SampleError`, `nucleide-nuclei` re-exports
  `DialectError` / `ParticlesError` / `RxnameError`, `nucleide-depletion`
  re-exports `CramError`, `nucleide-mcpl-io` re-exports `SswError`, and
  `nucleide-fluka-io` re-exports `MaterialError` / `UsrbinError`
  (`nucleide-mcnp-io` stays module-path-only by design); the nine crates
  lacking one gain `pub type Result<T>` (crate-root in `nucleide-depletion`,
  `nucleide-kinetics`, `nucleide-linalg`, `nucleide-mcpl-io`,
  `nucleide-nuclei`, `nucleide-spectroscopy`, `nucleide-vr-tools`,
  per-module in `nucleide-fluka-io` and all nine `nucleide-mcnp-io` error
  modules). `#![warn(missing_docs)]` now covers 17 of 18 crates —
  `nucleide-mcnp-io` (~132 missing docs) is the recorded deferral, owned
  by a future doc-volume cycle.

### Changed

- `#[non_exhaustive]` on all 38 public error enums across the 18 workspace
  crates (`nucleide-mcnp-io`'s nine per-module errors included): the
  API-stability sweep breaks exhaustive downstream `match`es on these
  enums — a wildcard arm is now required — and new variants become
  non-breaking additions. No data struct or non-error enum was touched,
  and 154 internal `Result<T, Error>` respellings to the new aliases
  change no public signature.

### Fixed

- Sampling kernel draw-count cap (`nucleide-linalg` `sample::MAX_SAMPLES` = 10M
  enforced in `validate_block` as `SampleError::TooManySamples` → `ValueError`,
  covering `sample_mvn`/`sample_lhs`/`sample_lognormal`).
- Corrected `inv_normal_cdf` docs: inputs at/beyond `[0, 1]` yield `NaN`, not
  `∓inf` (unreachable — `sample_lhs` clamps into `(0, 1)`).
- WASM UQ facades (`uqSample`/`sampleLhs`) reject `n < 2` with a
  moments-need-`n >= 2` message instead of the confusing `sample_cov`
  `DimensionMismatch`, and reject `seed >= 2^53` (f64 integer precision —
  the saturating `as u64` cast would collide streams).
- WASM `voxelTagsFromTotals` rejects empty `zone_of_voxel` (previously a
  vacuous success); `voxelPhotonSums` rejects non-finite `timeS` and empty
  `nuclides`; `partisnValidate` rejects empty `zones` (previously a vacuous
  pass) while `partisnRender` keeps rendering zoneless decks by design
  (callers validate separately).
- UQ demo clears stale results when switching MVN/LHS and rejects
  `seed >= 2^53` client-side; activation demo treats a blank voxel
  cooling-time field as an error instead of silently querying `0.0`.
- U6 validation gate reports FAIL instead of crashing with
  `ZeroDivisionError` on a hypothetical `n < 2` fixture.

## [0.9.0] - 2026-09-13

### Added

- MCTAL mesh tallies (`nucleide-mcnp-io` `mctal` headline): parses mesh-tally bodies (`detector_type <= -1`) — the 4-int
  mesh `f` line (`unknown ni nj nk`, bare-`f` or `f<tally>` spellings),
  `(ni+1)+(nj+1)+(nk+1)` `cora`/`corb`/`corc` bounds, the shared
  `d`/`u`/`s`/`m`/`c`/`e`/`t` cards, and `vals` pairs over the
  `ni*nj*nk` mesh cells (writer-loop order, `i` fastest). Synthetic
  `fixtures/mcnp/mctal/synthetic_mesh.mctal` oracle with closed-form
  val/err pairing, Rust tests, and Python `Mctal.mesh_tallies` /
  `mesh_tally_vals_array` views. Radiograph (`detector_type >= 3`),
  point-detector (tally names ending in 5), and perturbation bodies stay
  named-open errors.
- MCTAL S-tails: optional `tfc` blocks after standard
  tallies (`jtf` line plus 3–4-float data rows, exposed as
  `tallies[].tfc`), total/cumulative card variants (`ut`/`uc`/…, stored
  verbatim in `tallies[]` card `variant`), and third-token flags on
  `c`/`e`/`t` cards (stored verbatim as card `flag`), with synthetic
  `fixtures/mcnp/mctal/synthetic_tfc_variants.mctal` coverage. WWINP
  cylindrical (`nr=16`), Serpent `UnsupportedExpr`, SDEF versions beyond
  5/6, and SPE gating stay as-is (loud errors or unchanged surfaces —
  no second field-order source on disk).
- `validation/parsers_vs_refs.py`: MCTAL section comparing header
  scalars plus per-cycle keff/lifetime series against PyNE `Mctal` on
  the kcode-only fixtures; body-bearing files (standard, mesh,
  tfc/variant) are loud SKIPs — PyNE skips tally bodies without
   advancing past them, so no oracle exists there.
- Interactive completion bundle (WASM thin facades only —
  no new math): MCPL `readMcpl` / `writeMcpl` / `ssw2mcpl` / `mcpl2ssw`
  (bytes-based, gzip magic-sniffed, 200-row cap) with a file-upload
  `mcpl-io` tutorial staging the golden `ssw2mcpl_expected.mcpl` +
  `reference.w` pair; material `materialSeparate` / `materialBlend` /
  `cusumDetect`; spectroscopy `parseLinesTsv` / `energyBins` /
  `detectorEfficiency` / `parseDollarSpe` / `parsePlainSpe`; kinetics
  `inhourRho` / `stablePeriod` / `promptJump`; deterministic `parseRtflux`
  (kind-switch mirrors Python); ORIGEN per-step comparison over
  `parseOrigenTape6` plus chart-only ISOTXS totals, TAPE5 flux-vs-step,
  and TAPE6 per-nuclide views. Named-open and recorded in AD-13: R2S
  `VoxelTags` port (rayon feature unification), PARTISN writer, RTFLUX
  profile / ORIGEN per-step full views, multi-snapshot TAPE6 grammar.

- Shared dense-real weighted least-squares kernel (`nucleide-linalg`
  `lstsq` module, zero dependency change — faer 0.20 column-pivoted QR
  `solve_lstsq` route): `weighted_lstsq(X, y, w)` over caller
  `Vec<Vec<f64>>` inputs (non-negative caller weights, `sqrt(w)` row
  scaling, named errors for shape/finiteness/weight/underdetermined
  inputs) with closed-form synthetic gates. One module, two consumers —
  the UQ-lite sampler keeps its Cholesky/eigen factor path.
- Efficiency-coefficient fit (`nucleide-spectroscopy`, E7-fit extension,
  not a new equation number): `fit_efficiency` solves the log-space linear
  least squares (`y = ln eff`, `X[i][j] = (ln E_i)^j` for `fit == 1` or
  `(1/E_i)^j` for `fit == 2`, caller-supplied weights) through the shared
  kernel, returning `order + 1` coefficients for `detector_efficiency`.
  Peak search/fit and FWHM stay out; every other coefficient vector stays
  a caller input. Synthetic `fixtures/spectroscopy/efficiency_fit.json`
  oracle (recovery + round-trip at 1e-9).
- Python API: `nucleide.spectroscopy.fit_efficiency` thin wrapper over the
  new fit.
- `validation/spectroscopy_vs_pyne.py`: E7-fit synthetic recovery +
   round-trip gates (always run) plus a loud oracle-check SKIP — the upstream
  module ships no efficiency-fitting routine, so there is no container
   cross-check.
- Python facade bundle over existing Rust (thin wrappers
  only — no new math or data): `nucleide.nuclei` gains `rxname_label` /
  `rxname_doc` / `rxname_reaction` / `rxname_id_from_nucdelta` /
  `rxname_child` / `rxname_parent`, `particle_is_valid` /
  `particle_is_valid_pdc` / `particle_is_hydrogen` /
  `particle_is_heavy_ion`, and `dose_f1` / `dose_lung_model`;
  `nucleide.material` gains `mix_by_mass` / `mix_by_volume` /
  `specific_activity` / `materials_doc_to_xml` / `expand_elements` /
  `collapse_elements` (bare element symbols map to natural-element
  placeholders); `nucleide.fluka` gains `fluka_material_str` /
  `fluka_compound_str` / `fluka_builtin_set`; `nucleide.alara` gains
  `alara_validate_deck` / `alara_check_block` / `alara_flux_total` /
  `alara_flux_len` / `alara_output_totals` /
  `alara_output_total_activity` / `alara_photon_total_strength` /
  `alara_schedule_total_time`; `nucleide.origen` gains
  `origen_tape6_find` / `origen_tape6_total_activity` /
  `origen_tape9_find`; `nucleide.cccc` gains `cccc_rtflux_npoints` /
  `cccc_rtflux_point` / `cccc_rtflux_total` / `cccc_isotxs_find` /
  `cccc_isotxs_len`; `nucleide.fispact` gains `fispact_is_output`;
  `nucleide.enrichment` gains the six assay mass ratios
  (`prod_per_feed`, `tail_per_feed`, `tail_per_prod`, `feed_per_prod`,
  `feed_per_tail`, `prod_per_tail`) plus `alphastar_i`;
  `nucleide.kinetics` gains `from_ifp`; `nucleide.vr` gains `magic_with`
  plus `MeshSourceSampler` `mode` / `num_voxels` / `table_len`
  accessors; `nucleide.mcnp` is unchanged (dict surface complete);
  `nucleide.mcpl` gains `mcpl_statsum_validate` /
  `mcpl_statsum_comment`. `depletion`/`serpent` stay internals by design
  (no wrappers). Covered by `tests/test_facade_bundle.py` on
  synthetic inputs and committed fixtures.
- Python tutorials: new `docs/tutorials/python/
  vr-magic.md` (MAGIC + alias-table sampling) and
  `docs/tutorials/python/rxname.md` (registry, graph ops, particle
  helpers), a standalone-SWU section in `enrichment-cascade.md`, and
  worked mixing/specific-activity/`MaterialsDoc`/element-expansion
  sections in `build-materials.md` (replacing the "not yet exposed"
  note); all registered in `docs/tutorials/python/index.md` and
   `docs/README.md`.
- UQ-lite log-normal sampling (log-normal scope only; LHS stays out):
  `nucleide-linalg` `sample` gains a `LogNormal` perturbation convention
  (`parse("lognormal"/"log_normal"/"log-normal")`,
  `apply` as `nominal * exp(delta)`) plus `sample_lognormal` (shared
  Cholesky/eigen-clip factor path and ChaCha RNG as `sample_mvn`, then a
  `y = exp(x)` post-transform with log-space mean/covariance semantics)
  and the closed-form (U6) `lognormal_mean` / `lognormal_cov` helpers.
  Synthetic `fixtures/uq/lognormal_2x2.json` oracle with pinned seed and
  separate U5 validation gate (`validation/uq_lite_vs_sandy.py`: `ln(y)`
  MVN recovery plus closed-form mean; existing k-SE gates untouched).
  Python API: `nucleide.uq.sample_lognormal` / `lognormal_mean` /
  `lognormal_cov` thin wrappers (`perturb_energies(..., "lognormal")`
  now parses).
- SSW↔MCPL fidelity tail (`nucleide-mcpl-io` `ssw`; bitarray words stay
  out — no `isurf`/`rawtype` decode, no type-word reverse-engineering
  without a licensed type-table source; synthetic fixtures only, no
  transport semantics):
  `mcpl2ssw` gains `force_cs_to_one` (default `false` keeps the true cosine;
  opt-in forces the stored `cs` slot to `1.0` with `u`/`v` verbatim,
  removing the stored-cosine delta vs the 2.2.8 binary), `niss_override` (`None` =
  reference passthrough, the upstream-2.2.8 spelling verified over the
  synthetic pair; `Some(v >= 0)` stamps the tally-convention value, negative
  is `NissOutOfRange`), and `allow_polarisation` (default rejects non-zero
  input vectors with `PolarisationPresent`; opt-in drops them since SSW has
  no slot); `ssw2mcpl` gains opt-in `polarisation` (uniform vector +
  `has_polarisation`, non-finite is `InvalidPolarisation`), `universal_pdg`
  (single-kind inputs collapse file-wide, mixed is `MixedPdgForUniversal`),
  and `universal_weight` (equal weights collapse file-wide, mixed is
  `MixedWeightForUniversal`); the SSW-PDG table widens beyond neutron/gamma
  to electron 11 / positron −11 / proton 2212 with identical verbatim
  geometry (anything else stays loud `UnsupportedPdg`, no particle framework
  beyond the table). Each leg has named errors, synthetic both-direction
  Rust tests, Python facade options (`ssw2mcpl` `polarisation` /
  `universal_pdg` / `universal_weight`; `mcpl2ssw` `force_cs_to_one` /
  `niss` / `allow_polarisation`; five accepted kind spellings), committed
  `fixtures/mcpl/ssw_conversion/` goldens (extended, polarised, universal),
  and `validation/mcpl_vs_refs.py` gates S5–S8 in the same `mcpl` report
  (no `SECTION_ORDER` change). WASM/demo untouched (no cheap surface).
- Python API: `nucleide.mcpl.mcpl2ssw` keyword args `force_cs_to_one` /
  `niss` / `allow_polarisation` and `ssw2mcpl` options `polarisation` /
  `universal_pdg` / `universal_weight` over the widened kind table.

## [0.8.0] - 2026-09-13

### Added

- SSW↔MCPL conversion, neutron/gamma-only v1 (`nucleide-mcpl-io` `ssw`
  module, new `nucleide-mcnp-io` dependency for the SSW header/track
  types): `ssw2mcpl` maps kinetic energy verbatim (MeV), time shakes→ms
  (×1e-5), surface id→`userflags` (opt-out), and kind→PDG 2112/22, with
  single-precision + gzip-compressed-bytes output, an optional ≤100 MiB
  deck-embed blob, and named errors for other particle kinds;
  `mcpl2ssw` clones a reference SSW header (code/version/deck
  passthrough, `nrss`/`np1`/`orignp1` patched to the converted tally with
  the reference's table-2 sign) with surface id from `userflags` or an
  explicit `[1, 999999]` override, time ms→shakes (×1e5), and direction
  cosines propagated verbatim. Surface id and particle kind ride in as
  explicit caller parameters — the local SSW reader has no `isurf`/`rawtype`
  decode, so v1 never guesses them from `bitarray`. Deferred, all
  named-open: non-neutron/gamma particles (`UnsupportedPdg`), polarisation
  and universal codes, upstream `bitarray` type-word round-trips, and any
  transport semantics.
- Python API: `nucleide.mcpl.ssw2mcpl` / `mcpl2ssw` thin file-based wrappers
  over the new module (per-track surf+kind pairing, options dict, `.gz`
  transparent), plus a synthetic `fixtures/mcpl/ssw_conversion/` oracle pair
  (hand-framed two-track SSW reference with provenance notes + golden
  `ssw2mcpl` output).
- `validation/mcpl_vs_refs.py`: SSW round-trip gates S1–S4 on the synthetic
  pair plus a live upstream `ssw2mcpl`/`mcpl2ssw` (`mcpl-extra` 2.2.8)
  converter cross-check T3–T5 — count and closed-form energies agree both
  directions; the scripts SKIP loudly when absent, never fail without the
  oracle.
- UQ-lite sampling kernel, decay-only sub-scope (`nucleide-linalg`
  `sample` + `decay` modules, new `rand` 0.9 dependency with OS-entropy
  features off so the crate still builds for wasm): seeded
  multivariate-normal sampling over caller-supplied covariance blocks
  (`sample_mvn` — Cholesky primary with named `cholesky`/`eigen_clip`
  reporting, eigen-clipping fallback at 1e-12 relative floor, relative /
  absolute perturbation conventions, sample mean/unbiased-covariance
  convergence diagnostics with caller-supplied tolerances), plus a decay consumer
  (`perturb_branches` preserving the `1-BR(SF)` deficit by renormalisation,
  `perturb_energies`, finiteness-checked `passthrough`). Fission-yield
  perturbation stays a named-open hook (`FissionYieldsOpen`, waits on
  ENDF fission-yield tapes); ERRORR/NJOY, MF32/40, transport coupling, and
  vendored covariance stores stay out. Closed-form covariance-recovery
  gates on synthetic blocks only (pinned seeds, k-SE statistical
  tolerances, `fixtures/uq/`).
- Python API: `nucleide.uq` thin wrappers (`sample_mvn` / `sample_mean` /
  `sample_cov` / `check_convergence` / `perturb_branches` /
  `perturb_energies` / `passthrough` / `perturb_fission_yields`) over the
  new kernel.
- `validation/uq_lite_vs_sandy.py`: synthetic covariance-recovery gates
  U1–U4 (always run) plus a live SANDY `Samples.get_mean`/`get_cov` moment
  cross-check over the nucleide draws (tape-free, no NJOY) with a loud
  SKIP when SANDY/pandas are absent; tape-driven `sandy.sampling`/ERRORR
  comparisons stay NJOY-gated skips.
- UQ-lite tutorial + interactive demo: `docs/tutorials/python/uq-sampling.md`
  (seeded `sample_mvn` over the synthetic `fixtures/uq/` blocks, convergence
  diagnostics, decay perturbation) and a thin `uqSample` WASM facade with an
  interactive `tutorials/interactive/uq` sampling demo.

## [0.7.0] - 2026-09-12

### Added

- Fuel-cycle micro-adds (`nucleide-material`, `nucleide-enrichment`):
  `Material::separate` splits a composition into product/tails streams by
  per-nuclide efficiency (finite values in `[0, 1]`, unlisted nuclides go
  wholly to tails, per-nuclide mass conserved), `Material::blend` mixes
  streams at fixed ratios with explicit normalization (errors on empty,
  all-zero, or negative recipes — never a silent uniform split), and a
  dependency-free one-sided Page CUSUM change detector (`Cusum` with
  Welford running mean/variance, `new(ref_shift_k=0.5, alarm_h=4.0,
  startup=10)` tuning plus `update`/`status`/`statistic`/`reset`).
- Python API: `nucleide.material.separate_material` / `blend_material` /
  `Cusum` and `nucleide.enrichment.value_func` / `swu_per_feed` /
  `swu_per_prod` / `swu_per_tail` thin wrappers over the new core.
- `validation/enrichment_swu_vs_cyclus.py`: two-tier SWU oracle —
  hand-recomputed closed-form gates on the cyclus assay ladder (feed
  0.0072, product 0.05, tails 0.002, 10 units of product) plus a live
  `cyclus.toolkit.enrichment` cross-check with a loud SKIP when cyclus is
  absent.
- Full MCTAL tally bodies (`nucleide-mcnp-io`): the standard (non-mesh,
  non-radiograph) per-tally layout — `f`/`d`/`u`/`s`/`m`/`c`/`e`/`t` bin
  cards plus `vals` `(value, rel_error)` pairs — with a synthetic
  `fixtures/mcnp/mctal/synthetic_tally_bodies.mctal` oracle (closed-form
  val/err pairing). Mesh tallies, radiograph/point-detector specials,
  `tfc` blocks, total/cumulative variants, and perturbation bodies stay
  named-open errors, never silent skips.
- Python API: `Mctal.tallies` / `tally_nums` / `npert` views plus
  `Mctal.tally_vals_array` NumPy bridge (`(n_pairs, 2)`, float64 C-order)
  over the new bodies.
- MCPL particle-interchange reader/writer (new `nucleide-mcpl-io` crate):
  format versions 2 (read-only octahedral directions) and 3 (read/write,
  adaptive-projection directions), single/double precision, polarisation,
  user flags, universal PDG/weight, header comments/blobs with upstream
  `stat:sum` syntax validation, gzip-transparent paths, and byte-exact
  synthetic round-trips. Record-level interop is cross-checked both
  directions against the upstream 2.2.8 implementation. SSW conversion is
  deferred.
- Python API: `nucleide.mcpl.read_mcpl` / `write_mcpl` / `McplFile` thin
  wrappers over the new crate, plus a `docs/tutorials/python/
  mcpl-interchange.md` tutorial.
- `validation/mcpl_vs_refs.py`: synthetic round-trip gates plus an upstream
  `mcpl`-package record-level cross-check (count, energies, PDG codes) with
  a loud SKIP when absent.
- R2S per-voxel photon-source tags (`nucleide-r2s`, no MOAB/HDF5):
  `VoxelTags` mapping zone totals onto native-mesh voxels
  (`tag_zone_totals` copies, `split_zone_totals` conserves), plus
  `.photonSrc` group selection (`photon_groups_at`) and element-wise sums
  (`sum_group_strengths`) feeding real group spectra where the uniform-split
  placeholder stood.
- Python API: `nucleide.r2s.tag_zone_strength` / `photon_group_sums` thin
  wrappers over the new tags.
- SDEF runtime line tables (`nucleide-spectroscopy`, E9 input path): the
  documented decay-lines TSV interchange (`parse_lines_tsv` /
  `read_decay_lines` over `energy_MeV intensity` rows) feeding the existing
  normalizer — no vendored lines, no legacy reader, synthetic fixtures only.
- Python API: `nucleide.spectroscopy.parse_lines_tsv` / `read_decay_lines`
  thin wrappers over the new reader.
- `validation/spectroscopy_vs_pyne.py`: runtime TSV interchange oracle —
  Cs-137 ENSDF lines through `pyne.data` round-trip the reader (energy
  column within 1e-9, E9 normalization to 1.0), loud SKIP when PyNE is
  absent or the energy/intensity pairing stays open.
- Interactive MCTAL demo: `parseMctal` WASM function over the new tally
  bodies (header plus per-tally pair counts and totals) with an `mctal`
  mode in the MCNP parser demo (inline sample plus the staged
  `synthetic_tally_bodies.mctal` fixture) and E2E coverage.

### Changed

- Release workflow (`release.yml`): the crates.io publish list now covers
  every publishable workspace member (`nucleide-kinetics` and
  `nucleide-spectroscopy` were never listed, which failed the 0.6.0
  release at `nucleide-bindings` with `no matching package named
  nucleide-kinetics`), plus the new `nucleide-mcpl-io`, and a fail-fast
  guard rejects future omissions instead of retrying lag-waits against a
  deterministic failure.

### Fixed

- Hostile-input hardening in the MCTAL reader (`nucleide-mcnp-io`):
  allocations now follow the bytes actually present instead of reserving
  from declared bin/cycle counts (same class as the 0.5.0 `FILL`-matrix
  cap), the eight-card bin product saturates instead of wrapping, and
  fractional/negative bin counts are malformed-input errors rather than
  silent truncations. Valid files parse identically.

## [0.6.0] - 2026-09-12

### Added

- Plotly charts in the interactive Serpent and activation tutorials: the
  Serpent `_res.m` demo plots the `IMP_KEFF` mean ± uncertainty across burnup
  blocks (new `keff_history` field on `parseSerpentRes`), the `_det.m` demo
  plots per-detector tally spectra against energy-bin midpoints (new
  `spectra` field on `parseSerpentDet`, with Serpent 1 vs 2 column layouts),
  and the ALARA/FISPACT output tabs plot the reported `total` specific
  activity against cooling time. The activation demo stages the ALARA
  `output/sample2.out` and FISPACT `inventory.fis` fixtures behind new
  **Load sample** buttons. Website E2E smoke tests assert each chart renders.
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
  `sync-data.mjs` and embedded with the kit's theme-aware `<SvgFigure>`
  shortcode.
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

- Theory-page SVG figures no longer render as white boxes in the site's dark
  theme: the six figures in `docs/theory/figures/` are theme-aware (no opaque
  background, near-black strokes/text follow `currentColor`, muted grays follow
  `var(--muted-foreground)`, panel fills follow `var(--muted)`), and the pages
  inline them with the kit's `<SvgFigure>` shortcode instead of `<ImageFigure>`.
  Labels that sit on hardcoded pastel fills keep hardcoded dark text (the nucid
  Z/A digits no longer turn near-white in dark theme), and the peak-anatomy
  sideband/background labels gained pastel pill backdrops for contrast in both
  themes.
- `nucleide.serpent.read_serpent` no longer wraps matrix variables in a
  spurious one-element outer list: matrices now map directly to 2-D lists of
  row lists (`r["ABS_KEFF"][cycle]` instead of `r["ABS_KEFF"][0][cycle]`,
  same for `_dep.m`/`_det.m` matrices), matching the scalar/vector mapping.
  A matrix holding non-numeric values now raises `ValueError` instead of
  silently reading as an empty list. The Serpent tutorial and the parser
  cross-validation script were updated to the faithful shapes.
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
