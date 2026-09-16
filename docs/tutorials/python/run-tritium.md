---
title: Run Tritium
sidebar:
  order: 17
---

Nucleide's 1D tritium kernel solves diffusion-trapping permeation problems in a
slab with the `nucleide-tritium` crate. This tutorial covers the Python API;
the implementation lives in `crates/tritium`. For the physics — the
McNabb–Foster trap kinetics, the surface-law taxonomy, and the closed-form
permeation checks (G1–G6 on the theory page) — see the
[Tritium transport theory](../../theory/tritium.mdx) page.

## Steady-state permeation

```python
from nucleide.tritium import steady

out = steady(
    1e-3,  # slab length [m]
    32,  # cells
    1e-9,  # mobile diffusivity D [m²/s]
    {"kind": "dirichlet", "value": 1.0},  # left face concentration [mol/m³]
    {"kind": "dirichlet", "value": 0.0},  # right face
)
print(out["mobile"])  # mobile concentration at each cell centre
print(out["flux_left"], out["flux_right"])  # outward-positive face fluxes
print(out["inventory_mobile"], out["inventory_trapped"])
```

Every coefficient is caller data — the kernel ships no material tables. `D`
turns Arrhenius (`D0 * exp(-E_D / R / T)`) when `E_D` is passed; `temperature`
[K] and volumetric `source` [mol/m³/s] take one value or one per cell.
McNabb–Foster traps are one spec dict per species (`k0`, `p0`, `site_density`
required; `e_k`/`e_p` default to 0), returned per cell as `trapped[cell][trap]`:

```python
out = steady(
    1e-3,
    64,
    1e-9,
    {"kind": "dirichlet", "value": 1.0},
    {"kind": "dirichlet", "value": 0.0},
    traps=[{"k0": 1.0, "p0": 1.0, "site_density": 100.0}],
)
print(out["inventory_trapped"])  # total trapped inventory [mol/m²]
```

## Surface taxonomy

Each face is a spec dict whose `kind` picks the surface law; the two faces may
mix kinds freely:

- `"dirichlet"` — prescribed surface concentration, `"value"` [mol/m³].
- `"sieverts"` — diatomic-gas equilibrium `c = K_S √p`, `solubility` and
  `pressure` [Pa].
- `"henry"` — the linear variant `c = K_H p`, same keys as Sieverts.
- `"zero_flux"` — sealed wall or symmetry plane, no extra keys.
- `"recombination"` — molecular recombination outlet `J = K_r c²`,
  `"rate"` [m⁴/mol/s] (covered below).

```python
out = steady(
    1e-3,
    32,
    1e-9,
    {"kind": "sieverts", "solubility": 0.1, "pressure": 1e5},
    {"kind": "sieverts", "solubility": 0.1, "pressure": 1e3},
)
print(out["flux_right"])  # D·K_S·(√p₁ − √p₂)/L
```

## Transient breakthrough

`transient` integrates the same slab over a caller-supplied output grid `t`
[s], returning the mobile/trapped profiles and the outward flux series at each
output time. `time_lag` gives the permeation time lag `t_lag = L²/6D` and
`breakthrough` the analytic normalized outlet-flux series — the two standard
ways to read a breakthrough curve:

```python
from nucleide.tritium import breakthrough, time_lag, transient

D, L = 1e-9, 1e-3
times = [0.5 * time_lag(L, D), time_lag(L, D), 2.0 * time_lag(L, D), 5.0 * time_lag(L, D)]

out = transient(
    L,
    256,
    D,
    {"kind": "dirichlet", "value": 1.0},
    {"kind": "dirichlet", "value": 0.0},
    times,
)
print(out["flux_right"])  # outlet flux at each entry of times
print(breakthrough(D, L, times))  # analytic J(L,t)/J_ss at the same times
```

Initial profiles default to zero; pass `mobile0` / `trapped0` to override.
`method` selects `"crank_nicolson"` (default) or `"backward_euler"`, and
`rtol`/`atol`/`dt_max` control the adaptive stepper.

## Recombination boundaries

A recombination face makes the problem nonlinear: the outlet flux is
`K_r c²` at the surface concentration, not a fixed value. In steady state the
kernel closes each recombination face in closed form from the slab's affine
face response (the G5 construction); the transient reuses that face response
inside every implicit step — one Newton solve per step, warm-started (G6). See
the [theory page](../../theory/tritium.mdx) for the derivations; here is the
call shape:

```python
from nucleide.tritium import recombination_rate, steady, time_lag, transient

kr = recombination_rate(1e-6, 0.0, 500.0)  # K_r = kr0 * exp(-e_r / R / T)
ss = steady(
    1e-3,
    128,
    1e-9,
    {"kind": "dirichlet", "value": 1.0},
    {"kind": "recombination", "rate": kr},
)
print(ss["flux_right"])  # K_r * c_s² at the recombination face

t_end = 60.0 * time_lag(1e-3, 1e-9)
out = transient(
    1e-3,
    128,
    1e-9,
    {"kind": "dirichlet", "value": 1.0},
    {"kind": "recombination", "rate": 1e-6},
    [t_end],
    dt_max=1.0,
)
print(out["flux_right"])  # lands on the steady value above
```

Both ends may be recombination faces; a large `rate` recovers the Dirichlet
limit and a tiny one the sealed slab.

## Closed-form helpers

The closed-form helpers behind those checks are exported directly:

```python
from nucleide.tritium import irreversible_fill, langmuir, oriani, sieverts

print(sieverts(0.1, 1e5))  # c = K_S * sqrt(p)
print(oriani(1e-9, 3.0, 1.0))  # D_eff = D / (1 + K * N)
print(langmuir(100.0, 1.0, 1.0))  # c_t = N * K * c / (1 + K * c)
print(irreversible_fill(0.1, 1.0, 100.0, [1.0, 10.0]))  # c_t = N * (1 - e^{-kct})
```

## See also

- [Tritium transport theory](../../theory/tritium.mdx) for the T1–T2
  equation set and the G1–G6 closed-form checks.
- `tests/test_tritium.py` for replays of those checks and the input-error
  cases.
- [Interactive tritium permeation](../interactive/tritium.mdx) to run
  breakthrough curves in the interactive demo.
