# `crates/tritium` AGENTS.md

## Purpose

1D tritium diffusion-trapping kernel: Fickian mobile transport (T1) coupled
to N extrinsic McNabb–Foster trap species (T2) on a slab with
caller-supplied temperature, plus the Dirichlet/Sieverts/Henry/zero-flux
surface taxonomy. Multi-layer series stacks (e.g. W/Cu/CuCrZr first-wall
stacks) extend the same kernel per layer with linear internal interface
conditions — Sieverts (`u = c_m/K_S`) or Henry (`u = c_m/K_H`), both
folding into the tridiagonal θ-step matrix — or the vented-sink
recombination law (R-S: a single face concentration with the `K_r·x²`
desorption jump, closed on the CUT matrix in closed form for a lone gap
and by Newton for coupled faces). Pure 1D
finite-volume PDE; no multi-D, no FEM, no heat solve, no property tables.

## Ownership

Owns `crates/tritium/src/` (`params.rs`, `bc.rs`, `solve.rs`, `layers.rs`,
`error.rs`), the Python surface (`nucleide.tritium`, `tritium_*` in
`_internal`), `tests/test_tritium.py`, the `tritiumBreakthrough` WASM
facade, the `TritiumBreakthrough` demo, and
`docs/tutorials/interactive/tritium.mdx`. Replays (never rewrites)
`fixtures/tritium/`; shares `linalg::tridiag` with no other consumer yet.

## Local Contracts

- Equation set (T1–T2) is pinned: mobile balance, McNabb–Foster trap
  kinetics with the Langmuir equilibrium, Arrhenius caller-data stance,
  surface taxonomy. Derivations live in `docs/theory/tritium.mdx`; code
  comments cite equation labels, never external paths.
- Boundary taxonomy: `dirichlet` (default surface), `sieverts`
  (`c = K_S sqrt(p)`), `henry` (linear variant), `zero_flux` (symmetry /
  impermeable wall); `recombination` (`J = K_r c²`) is closed in the
  steady state by the exact face-response construction (G5) and in the
  transient by the per-step face Newton sharing the affine face-response
  construction (G6a–G6e) — never a silent pass. The same taxonomy governs
  the outer ends of multi-layer stacks (`layers.rs`).
- Multi-layer contract (G7/G8/G9/G10/G11, `layers.rs`): a `LayerStack` chains
  `N ≥ 1` caller-specified layers (thickness, cells, Arrhenius `D`,
  solubility `K`, per-layer traps/temperature/source). Internal interfaces
  are linear local-equilibrium laws — Sieverts (`u = c_m / K_S` continuous)
  or Henry (`u = c_m / K_H` continuous), flux continuous in both; the layer
  `solubility` carries `K_S` or `K_H` per the adjacent interface's law —
  or the vented-sink recombination law (R-S): `Interface::Recombination`
  carries the desorption rate `K_r` (finite, `> 0`), the gap holds a single
  face concentration `x ≥ 0` with half-cell fluxes `J_L = g_L·(c_L − x)`,
  `J_R = g_R·(x − c_R)` (`g = 2D/dx` per side, no solubility involved) and
  the sink `K_r·x²` (`J_L − J_R = K_r·x²`). The `K_r → 0` limit is a
  continuous-concentration joint, not a Sieverts law; a non-positive
  permeation drive `P ≤ 0` fails loudly in the steady state (the transient
  follows the G6 clamping stance per step). The flux-continuous product
  law has no spelling by construction (permanently rejected: spurious
  insulated root, symmetry breaking, blocked transient).
  Both linear laws share the face flux
  `J = (c_i/K_i − c_{i+1}/K_{i+1})/R_f` with `R_f` the sum of the two
  half-cell resistances `dx/(2·D·K)` (`K_S` or `K_H` per the adjacent
  law); the flux is linear in the cell values
  and folds directly into the tridiagonal θ-step matrix — no interface
  Newton (the spike outcome; the nonlinear G5/G6 face machinery is needed
  only for recombination *outer ends*, reused unchanged on the layered
  matrix). A vented gap instead CUTs the matrix (block-diagonal): adjacent
  values stay affine in the face value, so a lone gap closes in closed
  form `x = 2P/(√(Q²+4K_rP)+Q)` while coupled faces (several gaps, or
  recombination outer ends alongside) close by Newton with the analytic
  Jacobian; the per-step close fuses into the trap Picard loop in the G6
  slot, and the discrete balance closes against the reported boundary
  fluxes plus the desorption sinks from the recorded gap faces
  (`interface_faces`). A one-layer stack dispatches to the landed single-slab kernel
  and reproduces it exactly (G7c regression anchor). Goldens for the
  layered gates live in unit tests with recorded provenance (synthetic
  stacks, hand-derived series-resistance or vented-sink oracles), never in
  `fixtures/tritium/`.
- Discretization rule: cell-centred finite volume with implicit
  theta-stepping (Crank–Nicolson default, backward Euler on request);
  diffusion implicit, traps via the exact per-cell backward-Euler map with
  Picard coupling to rtol/atol. Method changes re-run the G1–G4 gates.
- Synthetic fixtures only: hand-built params plus closed-form values with
  recorded provenance; never evaluated-library data. `fixtures/tritium/`
  is read-only for this crate (no new fixture files; analytic oracles
  already exist).
- Oracle tolerances are gates, not claims: 1e-12 algebraic (G1/G3a/G3b
  isotherm/G3c/G4/G5a steady, G5b recovery at 1e-6), 1e-6 transient
  breakthrough curve (G2), 2% time-lag
  intercept, roundoff mass conservation. G6 has no algebraic oracle
  (asymptotic + self-convergence, plus the G6e independent cross-check):
  G6a steady-asymptote at 1e-6 relative on fluxes and flux balance at
  t = 60·t_lag with the G5a
  linear profile at 1e-9, G6b K_r→∞ flux recovery at 1e-6 and K_r→0
  zero-flux profile at 1e-5 absolute, G6c discrete mass balance at 1e-12
  relative against the reported K_r cf² faces, G6d dt-halving order
  bands 1.5–2.5 (Crank–Nicolson) and 0.7–1.3 (backward Euler) on two
  successive halvings in the resolved band (t = 1200 s, dt ≤ 3 s), G6e
  mid-transient cross-check at t = 0.5/1/2·t_lag against an independent
  method-of-lines solver (node-centred central FD + explicit RK4,
  ghost-node recombination end) at 1.5e-4 relative on the mobile-profile
  max-norm (6e-4 trapped profile, outlet flux 5e-3 relative with a
  1e-3·J_ss floor), trap-free and trap-coupled. G7/G8 multi-layer gates
  (synthetic stacks, in-test closed forms): G7a/G7b series-resistance flux
  at 1e-12 relative with a 1e-18 absolute floor (fluxes ~1e-7), interface
  flux continuity and half-cell-corrected interface potentials to roundoff
  (1e-12), G7c single-layer recovery exact (bit-for-bit against the landed
  kernel via dispatch), G7d property-continuous 2-layer ≡ landed at 1e-12,
  G7e layered recombination outer end at 1e-12 relative, G8a dt-halving
  bands as G6d (measured ≈2.0 CN / ≈1.0 BE), G8b layered steady asymptote
  at 1e-6 relative on fluxes and 1e-9 on the profile (Dirichlet and
  recombination outlets), G8c layered discrete mass balance at 1e-12
  relative. G9 Henry/mixed-interface gates (synthetic stacks, in-test
  closed forms, same resistance form with `K_H` in place of `K_S`):
  G9a 2-layer all-Henry stack series-resistance flux at 1e-12 relative
  with a 1e-18 absolute floor, interface flux continuity and half-cell-
  corrected interface potentials to roundoff (1e-12), profile the exact
  piecewise-linear Henry potential at the nodes; G9b 3-layer mixed
  Sieverts/Henry stack flux continuity across both interface laws to
  roundoff against the same closed form. G10/G11 vented-sink gates
  (synthetic stacks, in-test scalar-quadratic oracles): G10a 2-layer
  closed form at 1e-12 relative with a 1e-18 absolute floor on face value,
  fluxes, profile, gap residual, and the desorption-carrying outer
  balance; G10b `K_r→∞` upstream-block pinning and `K_r→0`
  continuous-joint recovery (uncut slab at 1e-12, no `K`-ratio jump);
  G10c 4-layer one-of-each-law stack with linear-gap continuity to
  roundoff and the super-block quadratic at 1e-12; G10d/G10e coupled
  Newton residuals at the Newton contract with the balance closed modulo
  the residuals; G10f traps on the owning layer's isotherm; G10g zero
  drive loud in the steady state (zero trajectory in the transient, G6
  clamp stance); G11a trap-free θ-balance with desorption at 1e-12
  relative; G11b dt-halving bands as G6d on the pinned pairs (CN on the
  spike's pair, BE on both); G11c vented steady asymptote at 1e-6 on
  fluxes and 1e-9 on the profile; G11d trap-coupled (backward Euler)
  mobile-plus-trapped balance with desorption at 1e-9.
- Out of scope (do not expand here): multi-D/FEM, heat coupling,
  plasma-facing implantation models, TBR coupling, FESTIM-file I/O, any
  dolfinx linkage, vendored D/K tables (caller-supplied only);
  the flux-continuous product-law recombination interface (permanently
  rejected — no spelling exists, never approximated).
- WASM/tutorial surface (owned): `tritiumBreakthrough` in `bindings/wasm`
  (trap-free permeation solve returning the `J/J_ss` series plus `t_lag`
  and `J_ss`; thin facade, fixed 200-cell grid), the
  `TritiumBreakthrough` demo in
  `website/src/components/interactive/TritiumBreakthrough.tsx`, and the
  `docs/tutorials/interactive/tritium.mdx` page. Demo presets stay
  synthetic (`fixtures/tritium/` values only).

## Work Guidance

- New analyses add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate depends on `nucleide-linalg` (tridiag)
  only; bindings depend on it, never the reverse.
- No `validation/*_vs_*.py` script: no external oracle exists (analytic-gate
  replay only, per the theory validation section). A FESTIM
  cross-check, if ever attempted, stays a loud-SKIP probe and never adds
  FESTIM/FEniCS to the validation container.

## Verification

- `cargo test -p nucleide-tritium` (includes fixture replay).
- `pytest tests/test_tritium.py` after `maturin develop`.

## Child NAD Index

None.
