---
title: Run Depletion
sidebar:
  order: 5
---

Nucleide depletion uses CRAM (Chebyshev Rational Approximation Method) to solve
the burnup matrix exponential. This tutorial covers the Python API; the
implementation lives in `crates/depletion`. For the physics and math, see the
[Depletion theory](../theory/depletion.mdx) page.

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

## See also

- [`crates/depletion/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/depletion/src/lib.rs)
  for the Rust API.
- `tests/test_depletion.py` for worked examples.
