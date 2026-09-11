---
title: Run Kinetics
sidebar:
  order: 8
---

Nucleide point kinetics solves prescribed-reactivity transients with the
`nucleide-kinetics` crate. This tutorial covers the Python API; the
implementation lives in `crates/kinetics`. For the physics and math, see the
[Point kinetics theory](../../theory/kinetics.mdx) page.

## Solve a step transient

```python
from nucleide.kinetics import solve

betas = [0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032]
lambdas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
out = solve(
    betas,
    lambdas,
    Lambda=1e-5,
    rho={"kind": "step", "t_step": 1.0, "rho_init": 0.0, "rho_final": 0.002},
    t=[0.5, 1.0, 1.05, 2.0, 10.0],
    n0=1.0,
)
print(out["n"])  # neutron level at each output time
print(out["C"][0])  # precursor groups at the first output time
```

Reactivities are in Δk and times in seconds. Precursor initials default to
equilibrium; pass `C0=[...]` to override. `method` selects `"trapezoidal"`
(default, second order) or `"backward_euler"` (extra damping); `rtol` /
`atol` control the adaptive stepper.

## Reactivity schedules

`rho` is a spec dict with `kind` selecting the schedule:

```python
step = {"kind": "step", "t_step": 1.0, "rho_init": 0.0, "rho_final": 0.002}
ramp = {
    "kind": "ramp",
    "t_start": 1.0,
    "t_end": 3.0,
    "rho_init": 0.0,
    "rho_rise": 0.002,
    "rho_final": 0.002,
}
pulse = {"kind": "impulse", "t_start": 1.0, "t_end": 2.0, "rho_init": 0.0, "rho_max": 0.002}
table = {"kind": "polyline", "times": [0.0, 1.0, 3.0], "values": [0.0, 0.0, 0.002]}
flat = {"kind": "constant", "rho": 0.0}
```

Steps are right-continuous (`t >= t_step` reads the post-step level) and the
solver lands exactly on every schedule knot.

## Inhour and prompt-jump helpers

```python
from nucleide.kinetics import equilibrium, inhour_rho, initial_rate, prompt_jump, stable_period

c0 = equilibrium(betas, lambdas, 1e-5, n0=1.0)
print(initial_rate(betas, lambdas, 1e-5, step, n0=1.0))  # == n0 * rho(0) / Lambda
print(stable_period(betas, lambdas, 1e-5, rho=0.002))  # asymptotic period in seconds
print(inhour_rho(betas, lambdas, 1e-5, omega=0.056))  # reactivity for a growth rate
print(prompt_jump(1.0, rho_before=0.0, rho_after=0.002, beta_total=sum(betas)))
```

`prompt_jump` needs `rho_after < beta_total` (no finite prompt equilibrium
past prompt critical); `stable_period` needs `0 < rho < beta_total`.

## OpenMC IFP data note

OpenMC's iterated-fission-probability estimator reports effective delayed
fractions and the generation time but no precursor decay constants: pass its
`betas`/`Lambda` through unchanged and supply `lambdas` from the same
delayed-neutron library the IFP run used (the `from_ifp` provenance note in
`crates/kinetics/src/params.rs`).

## See also

- [`crates/kinetics/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/kinetics/src/lib.rs)
  for the Rust API.
- `tests/test_kinetics.py` for worked examples.
