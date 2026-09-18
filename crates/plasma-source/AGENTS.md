# `crates/plasma-source` AGENTS.md

## Purpose

Tokamak fusion-neutron source creation: ring and point sources plus a
parametric Miller-geometry plasma, over the D-D (2.45 MeV) and D-T
(14.1 MeV) reactions — singly or as arbitrary D/T fuel mixtures at a shared
ion temperature, at distinct Maxwellian species temperatures
(`SpeciesIonTemperatures`: D-T at the mass-weighted `T_DT`, D-D at `T_D`),
or with the one pinned deuterium hot-tail fraction (`DeuteriumTail`: bulk
plus hot-tail sub-pairs at their own effective temperatures) — with
ion-temperature-broadened Gaussian spectra, a seeded
sampler to particle vectors, and MCNP `SDEF` + Serpent `src` card emission
with a drift report. An arbitrary-3D birth-rate lattice (`lattice.rs`)
drops axisymmetry for stellarator-class callers: a caller-supplied point
cloud with field-period symmetry reduction, sampled and emitted through the
same machinery.

## Ownership

Owns `crates/plasma-source/src/` (`reaction.rs`, `reactivity.rs`,
`spectrum.rs`, `source.rs`, `sample.rs`, `miller.rs`, `profile.rs`,
`parametric.rs`, `lattice.rs`, `ecrh.rs`, `emit_sdef.rs`, `emit_serpent.rs`, `report.rs`,
`error.rs`), the Python surface (`nucleide.plasma_source`,
`plasma_source_*` in `_internal`), `tests/test_plasma_source.py` plus
`tests/test_ecrh_access.py` (ECRH hand vectors + loud-error gates), and the
`plasma_source_vs_openmc.py` validation oracle.

## Local Contracts

- Provenance: physics is clean-room from published closed forms — Brysk,
  Plasma Phys. 15 (1973) 611 + Ballabio et al., Nucl. Fusion 38 (1998) 1723
  Table III (spectra); Fausser et al., Fus. Eng. Des. 87 (2012) 787
  (Miller map, L/H/A profiles); Bosch & Hale, Nucl. Fusion 32 (1992) 611
  in the Atzeni–Meyer-ter-Vehn parametrization (reactivity). Coefficients
  are transcribed as plain published facts and pinned by golden tests. The
  MIT `openmc-plasma-source`/`NeSST` packages are runtime oracles only
  (never ported); GPL upstreams (KDSource) stay never-read.
- Layering: depends on `nucleide-mcnp-io` (the `SDEF` dialect) and
  `nucleide-nuclei` (`PAR=` designators) only — never on `mcpl-io` or
  bindings. MCPL projection stays caller-side (the `vr-tools` KDE rule):
  the crate outputs particle vectors and card strings; the caller writes
  files with `nucleide-mcpl-io` / `nucleide-mcnp-io`.
- The MCNP `SDEF` accepted subset carries the ring keywords `AXS`/`RAD`/
  `EXT` (owned by the `mcnp-io` reader); emitted cards must round-trip
  byte-identically through `parse_sdef_text` (the E9 precedent).
  Serpent `src` rows are analytic by design (no Serpent source reader in
  the workspace).
- Parametric model: profiles are caller inputs (never computed, no
  equilibrium); emission strength is `f_fuel·n²·⟨σv⟩` with equimolar D-T
  (`f=1/4`) or pure D-D (`f=1/2`, the Fausser/openmc-plasma-source
  convention), or — via `FuelMixture` on the config / a `fuel={"D": f_D,
  "T": f_T}` dict in the Python spec — the Eriksson/DRESS mixture rule
  `S = n²·[f_D·f_T·⟨σv⟩_DT + (f_D²/2)·⟨σv⟩_DD]` at the shared profile `T_i`.
  The T-T neutron term is pinned as the third mixture term
  `+ f_T²·⟨σv⟩_TT` (neutron multiplicity 2 folded in, reacting at `T_T`,
  `FuelMixture::tt_neutron_coefficient`) but not transported — no
  publishable T-T reactivity fit (Bosch & Hale 1992 covers only D(d,n),
  D(d,p), D-T, D-³He) or Ballabio-class line (three-body continuum) exists,
  and the oracle's vendored tables stay out of scope. The D(d,p)T proton
  *rate* is accounted alongside the neutron source
  (`proton_strength_density` / `total_proton_strength`, 50/50 with the D-D
  neutron branch, bit-for-bit; `proton_accounting` facade) while the
  sampler and the cards stay neutron-only; proton transport stays out.
  The normalization is pinned with hand vectors in the `parametric` module
  rustdoc, and the recovery anchors are regression gates: `f_D=f_T=1/2`
  reproduces the equimolar kernel's D-T branch bit-for-bit, `f_D=1`
  reproduces pure D-D bit-for-bit; the landed single-fuel expression trees
  are untouched. Birth positions weight `S·R·|J|` (the `R` toroidal factor
  is mandatory — without it the core is over-represented); birth energies
  run a per-particle branch roulette between the two Ballabio lines. Card
  emission for parametric sources is product-form: radial/vertical/energy
  *marginals* as discrete histograms, and the drift report must carry the
  joint-correlation distance (half the L1 distance between the true
  `(r, z)` joint and the product of marginals). An optional
  `ToroidalSector` (`start_angle`/`rotation_angle` pair, or the same keys in
  the Python spec) restricts births to
  `[start_angle, start_angle + rotation_angle)` with uniform birth angles
  and totals scaled by `rotation_angle/2π`; the exact full-rotation spelling
  recovers the full torus bit-for-bit (regression gate). A partial sector
  adds a uniform angle-bin marginal to the cards (`PHI=D4` on SDEF, `phi d4`
  on Serpent) plus a `toroidal sector` drift row; sector angles carry their
  own loud errors (`Error::NonFinite`, `Error::InvalidSector`) at
  construction.
- Loud boundary (`Error::NotYetSupported` / `ValueError`, never a guess):
  T-T neutron transport
  (third-term normalization pinned, `tt_neutron_coefficient`; no
  publishable fit or line exists) and proton transport (the D(d,p)T proton
  *rate* is accounted: `proton_strength_density` /
  `total_proton_strength`, `proton_accounting` facade; sampler and cards
  stay neutron-only), reactant distributions beyond the one pinned
  deuterium hot-tail fraction on a D/T mixture (the rest of the full
  Eriksson generalization; the per-species-temperature form is supported via
  `SpeciesIonTemperatures` — including on single-fuel configs, which react
  at the pair temperatures — the single-tail-temperature form via
  `DeuteriumTail` — both parametric-only, there is no tail spelling on
  lattice configs), per-species temperatures without a fuel mixture on
  lattice configs, and profile
  self-consistency (zero total strength). Mixture fractions, species
  temperatures, tail parameters, lattice clouds/symmetry, and sector angles carry their own loud errors
  (`Error::NonFinite`, `Error::InvalidFuelMixture`,
  `Error::NegativeIonTemperature`, `Error::InvalidTail`,
  `Error::InvalidLattice`, `Error::InvalidSymmetry`, `Error::InvalidSector`) at
  construction.
- Units: cm, MeV, keV (card convention); profile density is m⁻³ (Fausser
  convention, converted internally); reactivity is exposed in m³/s. The
  validation oracle converts at the openmc-plasma-source boundary (m, eV).
- Sampler determinism: pinned seed reproduces the stream on a given
  platform (libm functions); golden vectors in `sample.rs` pin the RNG
  across refactors, and the parametric inverse-CDF tables are built at
  construction (no per-particle quadrature).

## Work Guidance

- Miller/Jacobian changes must keep the pinned gates green: `J == r` (torus
  limit), `J == κ·r·(1 − 2·esh·r·cosθ/a²)` at δ=0, area `π·κ·a²`, volume
  `2π²·κ·a²·R₀` (any esh), finite-difference agreement at general δ, and
  the strength ratio gates (`<σv>` scaling).
- New emission dialects follow `emit_serpent.rs`: golden card tests plus
  drift rows (`reparsed` false when no reader exists).
- Fuel-model changes keep the convention-first rule: the normalization
  (module rustdoc hand vectors) and the exact recovery anchors move before
  any sampler change, and the single-fuel sampling stream must stay
  bit-for-bit (no new RNG draws on the `fuel_mixture: None` path).
- Generalizing the mixture model further (T-T branch transport, tail
  shapes beyond the one pinned `DeuteriumTail`) moves the
  `NotYetSupported` boundary only after its own pinned normalization brief.
- Lattice changes keep the sampler draw accounting: no new RNG draws on any
  path (the symmetry copy comes from the within-bin fraction), so the
  folded/expanded bit-identity and the single-point recovery anchors hold.

## Verification

- `cargo test -p nucleide-plasma-source` (map/Jacobian closed forms,
  profile/reactivity goldens, sampler moment gates, golden cards, reader
  round trips, mixture hand vectors + exact recovery anchors + loud
  fraction errors, branch gates: T-T coefficient vectors, proton
  bit-identity with the D-D branch, the no-T recovery anchor, sector gates:
  rotation/2π hand vectors, full-rotation bit-for-bit recovery, loud angle
  errors, species gates: T_DT hand vectors, pair-temperature strength and
  per-branch spectra, equal-T bit-for-bit recovery, loud temperature
  errors, tail gates: sub-rate hand vectors at 1e-12, five-sub-branch
  spectra, eta=0 bit-for-bit recovery, loud tail errors, lattice gates:
  two-point hand vectors, single-point bit-for-bit point-source recovery,
  expand/fold exactness, folded-vs-expanded stream bit-identity (incl. the
  single-node regression), ring-limit moment closure, equal-pair recovery,
  card round trips with the lattice row, loud lattice/symmetry errors).
- `pytest tests/test_plasma_source.py` after `maturin develop`.
- `validation/plasma_source_vs_openmc.py` runs inside `run_all.sh`
  (two-part: P1–P13 always, O1–O13 vs openmc-plasma-source/NeSST
  container-only with loud SKIP outside; the Containerfile layer installs
  `openmc-plasma-source` + `NeSST` and shims `scipy.integrate.cumtrapz`
  in the oracle process for scipy ≥ 1.14).

## Child NAD Index

None.
