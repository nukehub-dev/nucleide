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
`None` for decay-only), keyed `"Name:reaction"` as in `deplete`. Only
`integrator="predictor"` is exposed; per-nuclide `decay_heat` maps use the
same chain → ENDF/B-VII.1 → `0.0` energy resolution as the core (synthetic
nuclides report `0.0`). For the exact call shape see
`deplete_series` in the
[Python API](../../reference/python-api.mdx#nucleidedepletion) and
`tests/test_03_depletion_series.py`.

## See also

- [`crates/depletion/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/depletion/src/lib.rs)
  for the Rust API.
- `tests/test_depletion.py` for worked examples.
