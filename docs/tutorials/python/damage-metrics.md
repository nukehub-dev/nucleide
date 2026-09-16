---
title: Damage metrics
sidebar:
  order: 18
---

Nucleide's damage metrics fold a caller-supplied multigroup neutron flux
against caller-supplied response cross sections to give the standard
first-wall engineering numbers: NRT-dpa and arc-dpa displacement damage,
helium or hydrogen gas production in atomic parts per million (appm), and
the He/dpa ratio. `fold_uq` propagates uncertainty through any of the
linear folds. This tutorial covers the Python API; the implementation
lives in the `nucleide-damage` crate. Every flux and cross section is
caller data — the kernel ships no nuclear-data tables.

## A first fold

`nrt_dpa` takes three lists and an exposure time:

- `bounds` — the `G + 1` group boundaries in MeV, strictly increasing.
  They document the grid; the fold never uses the widths.
- `flux[g]` — the per-group integrated flux (n/cm²/s through group g),
  not a per-unit-energy density.
- `response[g]` — the group damage cross section in barns, condensed by
  your own data pipeline (an NJOY- or SPECTER-class code, for example).
- `seconds` — the exposure time at that flux.

```python
from nucleide.damage import arc_dpa, nrt_dpa

bounds = [0.0, 0.1, 1.0, 20.0]  # three groups: G+1 MeV boundaries
flux = [1.0e12, 2.0e12, 4.0e12]  # per-group integrated flux [n/cm²/s]
dpa_xs = [100.0, 200.0, 50.0]  # group NRT dpa cross section [barns]

print(nrt_dpa(flux, dpa_xs, bounds, 2.0))  # two seconds of exposure
print(arc_dpa(flux, dpa_xs, bounds, 2.0))
```

```text
1.4e-09
1.4e-09
```

The fold is piecewise-constant per group: `dpa = seconds·1e-24·Σ_g
flux[g]·response[g]`, where `1e-24` converts barns to cm². Hand-check:
`1e12·100 + 2e12·200 + 4e12·50 = 7e14`, so `2·1e-24·7e14 = 1.4e-9` dpa.
`arc_dpa` is the same fold applied to arc-corrected cross sections — the
correction lives in how the response was condensed, not in the fold, so
both calls agree on the synthetic response above.

Because the fold never re-weights inside a group, moving the boundaries
changes nothing:

```python
wide = [0.0, 5.0, 10.0, 20.0]  # same groups, stretched
print(nrt_dpa(flux, dpa_xs, wide, 2.0))  # identical: widths never enter
```

```text
1.4e-09
```

Zero-flux groups contribute exactly 0, and nothing ever divides by a group
flux, so a dead group cannot produce `NaN` or `inf`:

```python
sparse = [0.0, 2.0e12, 0.0]
resp = [1.0e3, 1.5e2, 7.0e3]
print(nrt_dpa(sparse, resp, bounds, 1.0))  # dead groups contribute exactly 0
```

```text
3e-10
```

## Seconds or fluence

`seconds` multiplies a per-second flux. To fold a time-integrated fluence
instead, park the fluence in `flux` and pass `seconds=1.0`:

```python
fluence = [4.78373e22]  # time-integrated fluence [n/cm²], one group
fe56_xs = [191.18]  # Fe-56 dpa cross section [barns]
print(nrt_dpa(fluence, fe56_xs, [0.0, 20.0], 1.0))
```

```text
9.145535014
```

This is the single-group Fe-56 spot the repo's validation harness
cross-checks against the SPECTER report. Either way, `seconds` must be
strictly positive.

## Gas production and the He/dpa ratio

`gas_appm` is the same fold with the per-atom normalization — the scale is
`1e-18` instead of `1e-24`, adding the parts-per-million factor. The
response counts whichever gas your cross sections describe (He or H, one
nuclide's production per target atom). `he_dpa_ratio` folds the He
production and damage responses over the same flux and divides; the
seconds cancel, leaving appm per dpa.

```python
from nucleide.damage import gas_appm, he_dpa_ratio

flux = [1.0e13, 1.0e12, 1.0e11]
he_xs = [0.2, 0.4, 0.6]  # group He production cross section [barns]
damage_xs = [50.0, 100.0, 200.0]

print(f"{gas_appm(flux, he_xs, bounds, 3.0):.6g}")  # He [appm]
print(f"{he_dpa_ratio(flux, he_xs, damage_xs, bounds, 3.0):.6g}")  # appm/dpa
```

```text
7.38e-06
3967.74
```

A damage fold of exactly zero raises a clear error — never `inf`.

## Closed-form damage functions

The functions a condensation pipeline would tabulate into group cross
sections are exported directly, so you can spot-check a response or build
one by hand. `lindhard_partition` gives the fraction of recoil energy left
for displacements after electronic losses; `damage_energy` is the recoil
energy times that fraction; `nrt_displacements` is the
Norgett–Robinson–Torrens count (0 below `E_d`, 1 up to `2·E_d/0.8`, then
`0.8·T_dam/(2·E_d)`); and `arc_displacements` scales the high branch by
the Nordlund efficiency `arc_efficiency`.

```python
from nucleide.damage import (
    arc_displacements,
    arc_efficiency,
    damage_energy,
    lindhard_partition,
    nrt_displacements,
)

print(f"{lindhard_partition(1.0e5, 'Fe56', 'Fe56'):.6g}")  # fraction
print(f"{damage_energy(1.0e5, 'Fe56', 'Fe56'):.6g}")  # damage energy [eV]
print(f"{nrt_displacements(1.0e5, 40.0, 'Fe56'):.6g}")  # NRT count
print(f"{arc_efficiency(1.0e5, 40.0, -0.55, 0.3):.6g}")  # arc efficiency
print(f"{arc_displacements(1.0e5, 40.0, 'Fe56', -0.55, 0.3):.6g}")  # arc count
```

```text
0.618651
61865.1
618.651
0.315671
198.221
```

`E_d` (here 40 eV) and the arc constants (`b_arc`, `c_arc`) are
caller-supplied material data; `b_arc` must be negative and `c_arc` must
lie in `(0, 1)`.

## Error cases

Malformed input raises a clear error naming the cause — never a silent
clamp or a non-finite result:

```python
try:
    nrt_dpa([-1.0, 1.0, 1.0], dpa_xs, bounds, 1.0)  # negative group flux
except ValueError as exc:
    print("ValueError:", exc)

try:
    nrt_dpa(flux, dpa_xs, bounds, 0.0)  # non-positive exposure time
except ValueError as exc:
    print("ValueError:", exc)

try:
    he_dpa_ratio([0.0, 0.0, 0.0], he_xs, damage_xs, bounds, 1.0)  # zero dpa
except ValueError as exc:
    print("ValueError:", exc)
```

```text
ValueError: damage: flux must be non-negative
ValueError: damage: seconds must be strictly positive
ValueError: damage: He/dpa ratio at zero dpa is undefined (zero flux or zero damage response)
```

The same validation covers mismatched list lengths (`expected 3, got 2`),
non-increasing bounds, `NaN`/`inf` entries ("non-finite"), and negative
response values.

## Propagating uncertainty

`fold_uq` perturbs the stacked `[flux, response]` vector with seeded
multivariate-normal (MVN) draws and refolds the metric per draw. `mean`
and `cov` describe relative perturbations of that vector (dimension `2G`,
flux first); the seed pins the stream, so identical inputs reproduce
bit-for-bit. The returned dict carries the sample moments, the exact
expectation of the fold, the first-order propagated standard deviation,
and a correctness check: `passed` is `True` when the sample moments land
within `k` standard errors of the analytic values.

```python
from nucleide.damage import fold_uq

flux = [1.0e13, 3.0e12]
resp = [80.0, 160.0]
bounds = [0.0, 0.5, 20.0]
mean = [0.0, 0.0, 0.0, 0.0]
cov = [  # relative variances: 2% and 3% on flux, 5% and 1% on response
    [0.0004, 0.0, 0.0, 0.0],
    [0.0, 0.0009, 0.0, 0.0],
    [0.0, 0.0, 0.0025, 0.0],
    [0.0, 0.0, 0.0, 0.0001],
]
out = fold_uq("nrt_dpa", flux, resp, bounds, 2.0, mean, cov, 4000, 20260915, 5.0)
print(out["mean"], out["std"])
print(out["analytic_std"], out["expected"])
print(out["passed"])
```

```text
2.559163388924067e-09 9.264502319685633e-11
9.135425551116927e-11 2.56e-09
True
```

Uncertainty on the He/dpa ratio is not yet supported — the ratio is a
nonlinear statistic. Keep covariances in the small-perturbation regime: a
draw that pushes a flux or response negative surfaces the same clear error
as above, never a silent clip.

## See also

- `tests/test_damage.py` for replays of these checks and the input-error
  cases.
- [UQ sampling](uq-sampling.md) for the seeded MVN engine `fold_uq`
  reuses.
- The repo validation harness folds published SPECTER report values for
  comparison (`validation/damage_vs_specter.py`).
- No damage-data pipeline of your own? `nucleide.damage.specter_table`
  serves the vendored SPECTER Table VII fallback (one spectrum-averaged
  cross section per element for seven named spectra) as an explicit
  one-group response.
