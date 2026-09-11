# `crates/kinetics` AGENTS.md

## Purpose

Prescribed-reactivity point-kinetics transients: one PKE solve over
caller-supplied delayed-neutron data, plus the inhour and prompt-jump
analyses. Pure ODEs; no transport, no feedback, no tabulated data.

## Ownership

Owns `crates/kinetics/src/` (`params.rs`, `reactivity.rs`, `solve.rs`,
`inhour.rs`, `error.rs`), the Python surface (`nucleide.kinetics`,
`kinetics_*` in `_internal`), `tests/test_kinetics.py`, the
`fixtures/kinetics/` synthetic oracles, and the `kinetics_vs_pyrk.py`
validation oracle.

## Local Contracts

- Equation set (E1–E4) is pinned: PKE right-hand sides, equilibrium
  initials, inhour relation, prompt-jump formula. Derivations live in
  `docs/theory/kinetics.mdx`; code comments cite equation labels, never
  external paths.
- Insertion taxonomy: `constant` (default), `step`, `impulse`, `ramp`
  mirror the upstream reactivity-insertion classes with right-continuous
  knot semantics; `polyline` is the Nucleide extension for table input.
- Stiffness rule: no explicit-solver parity. The integrator stays implicit
  (A-stable θ-method) with exact knot stepping; method changes re-run the
  O1–O5 gates.
- Synthetic fixtures only: hand-built params plus closed-form values with
  recorded provenance; never evaluated-library data. New oracles extend
  `fixtures/kinetics/` and the fixture tests in the same change.
- Oracle tolerances are gates, not claims: 1e-12 algebraic, 1e-10 inhour
  residual, 1e-6 1-group, 1–2% 6-group tail, 1e-3 PyRK cross-check.
- Out of scope (do not expand here): thermal-hydraulic feedback,
  flux-coupled transients, time-dependent tallies.
- WASM/tutorial surface (owned): `kineticsTransient` in `bindings/wasm`
  (step-reactivity PKE solve returning the `n(t)` series plus the E4
  prompt-jump value; thin facade, default solver options), the
  `KineticsTransient` demo in
  `website/src/components/interactive/KineticsTransient.tsx`, and the
  `docs/tutorials/interactive/kinetics.mdx` page. Demo presets stay
  synthetic (`fixtures/kinetics/` values only).

## Work Guidance

- New analyses add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate is independent; bindings depend on it,
  never the reverse.
- Validation oracle needs the upstream package at run time
  (`Containerfile` PyRK layer); a missing oracle is a recorded SKIP, never
  a silent pass.

## Verification

- `cargo test -p nucleide-kinetics` (includes fixture replay).
- `pytest tests/test_kinetics.py` after `maturin develop`.
- `validation/kinetics_vs_pyrk.py` runs inside `run_all.sh` (two-tier:
  analytic gates always, PyRK cross-check when importable).

## Child NAD Index

None.
