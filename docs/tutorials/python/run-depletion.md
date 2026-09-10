---
title: Run Depletion
sidebar:
  order: 3
---

Nucleide depletion uses CRAM (Chebyshev Rational Approximation Method) to solve
the burnup matrix exponential. This tutorial covers the Python API; the
implementation lives in `crates/depletion`. For the physics and math, see the
[Depletion theory](../../theory/depletion.mdx) page.

## Load a depletion chain

```python
from nucleide.depletion import read_chain

chain = read_chain("path/to/chain.xml")
print(chain.nuclides[:10])
```

## Run a CRAM solve

```python
from nucleide.depletion import deplete

n0 = {"U235": 1.0, "U238": 20.0}
rates = {"U235:(n,gamma)": 1e-20}  # optional reaction rates, keyed "Name:reaction"
result = deplete(chain, n0, dt=1e6, rates=rates, order=48)
```

`order` selects CRAM-16 or CRAM-48.

## Pick the solver kernel with `method=`

`deplete`, `deplete_series`, `DepletionSystem.solve` / `solve_vec`, and
`Inventory.decay` all accept `method=` (`"cram16"`, `"cram48"`,
`"bateman"`, `"bateman_hp"`; default `"cram48"`). An explicitly
non-default `method` overrides `order`; omit it (or pass `"cram48"`) to
keep the legacy `order` behavior. The Bateman arms evaluate the analytic
closed form per step and fall back to CRAM-48 — never an error — outside
decay-only triangular topology: (D1) near-degenerate half-lives
(`1e-12` relative gap with a live coupling path), (D2) cyclic or
out-of-order topology, (D3) stable nuclides via inline limit forms, and
(D4) live reactions or fission. The closed form sums alternating-sign
terms, so long chains (tens of members) or extreme half-life spreads can
suffer catastrophic cancellation in float64 — `bateman_hp` (magnitude-
sorted terms, Neumaier compensation, `exp_m1` for small `λt`) extends the
usable range, but ill-conditioned chains should use CRAM-48, which stays
the default. `n0` entries must be finite atom counts `>= 0` on the
Bateman path (other spellings are errors, as are non-positive or
non-finite `dt` values other than exactly `0.0`, which echoes the input):

```python
from nucleide.depletion import deplete, read_chain

chain = read_chain("fixtures/depletion/chain_abc.xml")
n0 = {"A": 1.0e15}
for method in ("cram48", "cram16", "bateman", "bateman_hp"):
    out = deplete(chain, n0, dt=1.0e5, method=method)
    print(method, f"{out['B']:.4e}")
```

For the kernel math and the D1–D4 limits see the
[Depletion theory](../../theory/depletion.mdx) page.

## Nuclide names, decay branches, and isomer masses

Chain lookups key on GNDS names, and the `nuclei` tables behind them
accept free-form spellings. `normalize_nuclide` resolves symbol-first
(`N15` stays nitrogen) then mass-first (`92235` stays a ZAID); bare
symbols are rejected and `Ir-192n` resolves to state 2:

```python
from nucleide.nuclei import (
    atomic_mass,
    decay_branches,
    decay_branch_fraction,
    normalize_nuclide,
)

print(normalize_nuclide("241Pu"))  # Pu241
print(normalize_nuclide("Ba-137m"))  # Ba137_m1
print(normalize_nuclide("Ir-192n"))  # Ir192_m2
print(decay_branches("K40"))
print(decay_branch_fraction("K40", "Ca40"))
print(atomic_mass("Ba137_m1"), atomic_mass("Ba137"))
```

`decay_branches` lists per-branch `(daughter, fraction, mode)` rows from
the ENDF/B-VIII.0 decay tapes (SF/fission branches dropped); isomer
masses extend the AME2020 ground-state table with each isomer tape's
File-1 MT451 excitation energy (`m_ground + ELIS/931.49410242 u`),
falling back to the ground-state mass where no tape exists. Q-values
stay ground-state-only.

## Run a multi-step series

`deplete_series` threads single-step CRAM solves forward over a list of step
lengths (predictor in the multi-step sense: rates held constant within each
step). With a synthetic decay chain on disk:

```python
from nucleide.depletion import deplete_series, read_chain

chain = read_chain("path/to/chain_synth.xml")
out = deplete_series(chain, {"A": 1.0e15}, [1.0e5, 1.0e5, 1.0e5])
print(out["times"])  # cumulative seconds: [1e5, 2e5, 3e5]
print(out["atoms"][0])  # atom counts after the first step
print(out["activity"][0])  # Bq per nuclide (A = λN)
```

Per-step rates use `rates` (every step) or `rates_list` (one entry per step,
`None` for decay-only), keyed `"Name:reaction"` as in `deplete`. The
`integrator` argument selects `"predictor"`, `"cecm"`, or `"cf4"` (same
core series as Rust `integrate`, minus the `t = 0` row); per-nuclide
`decay_heat` maps use the same chain → ENDF/B-VII.1 → `0.0` energy
resolution as the core (synthetic nuclides report `0.0`). For the exact
call shape see
`deplete_series` in the
[Python API](../../reference/python-api.mdx#nucleidedepletion) and
`tests/test_03_depletion_series.py`.

## Unit-aware inventories

`Inventory` wraps atom counts with activity/mass/mole units over a chain,
plus fractions, readable half-lives, arithmetic, and CSV round-trip:

```python
from nucleide.depletion import Inventory

inv = Inventory(chain, {"Co60": 1.0}, units="Ci")
aged = inv.decay(1.0, time_unit="y")
print(aged.activities("Bq"), aged.masses("g"))
print(inv.half_lives_readable())
```

`cumulative_decays` integrates decays over one step, and `progeny` /
`branching_fraction` / `decay_mode` / `chain_edges` expose chain lineage
as plain data (see `tests/test_inventory.py`).

## See also

- [`crates/depletion/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/depletion/src/lib.rs)
  for the Rust API.
- `tests/test_depletion.py` for worked examples.
