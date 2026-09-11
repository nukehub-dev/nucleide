# fixtures/kinetics/ — synthetic point-kinetics inputs and oracles

Hand-built synthetic data only — no laboratory or evaluated-library values.
Params use round numbers chosen for this repo; expected values are closed-form
consequences of the PKE system (E1–E4), computed once in float64 with the
provenance recorded in each file.

- `step_oracle.json` — 1-group params + a `0 → 0.002 Δk` step at `t = 1 s`:
  algebraic `initial_rate`, the E4 `prompt_jump_factor`, and the
  two-exponential closed form `n(t)` at seven times since the step.
- `ramp_table.csv` — shared transient *input*: `t_s,rho_dk` sampling a
  `0 → 0.002 Δk` ramp over `[1, 3] s` with a hold to `10 s` (exact linear
  values, 0.25 s spacing). Both the Rust tests (as a polyline) and the
  `validation/kinetics_vs_pyrk.py` PyRK cross-check drive from this table,
  so the two codes see identical reactivity histories.
- `inhour_check.json` — 6-group synthetic params + the stable period at
  `ρ = 0.002 Δk` from Newton iteration on the inhour equation (spelled
  `inhour`, not the plan draft's `ihour` typo).
