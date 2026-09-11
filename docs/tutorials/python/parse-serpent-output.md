---
title: Parse Serpent Output
sidebar:
  order: 2
---

Nucleide reads Serpent 1 and Serpent 2 MATLAB-style output files without
running Serpent itself. This tutorial shows the Python API; equivalent Rust
types live in the `nucleide-serpent-io` crate.

## `_res.m`: run results

A `*_res.m` file collects run metadata, k-eigenvalues, and cycle statistics
as `NAME = [...]` assignments.

```python
from nucleide.serpent import read_serpent

r = read_serpent("fixtures/serpent/serp2_res.m", "res")
print(r["VERSION"][0], r["POP"][0])

keff = r["ABS_KEFF"]  # one [mean, stdev] row per cycle block
print(keff[0])
```

`read_serpent(path, kind)` returns a plain dict keyed by variable name.
Vectors (`VERSION`, `POP`, ...) are flat lists. Matrix-valued variables
come back as 2-D lists of row lists (one row per Serpent block), so index
`keff[cycle]` directly to reach a row.

## `_dep.m`: depletion inventories

A `*_dep.m` file stores per-step inventories: `ZAI`/`NAMES` identify the
nuclides, `DAYS`/`BU` the depletion steps, and `TOT_*`/`MAT_<name>_*`
hold the large density matrices.

```python
d = read_serpent("fixtures/serpent/sample1_dep.m", "dep")
print(len(d["ZAI"]), d["DAYS"], d["BU"])

adens = d["TOT_ADENS"]  # rows align with ZAI; one column per step
print(adens[0])
```

Per-material blocks follow the same layout (`MAT_fuel_ADENS`,
`MAT_fuel_H`, `MAT_fuel_BURNUP`, ...), with totals in `TOT_ADENS`,
`TOT_MASS`, and friends.

## `_det.m`: detector tallies

A `*_det.m` file holds one flat array per detector, reshaped by the reader
into rows of bin indices plus a tally value and relative error.

```python
det = read_serpent("fixtures/serpent/serp2_det.m", "det")
bins = det["DET1"]  # one 13-column row per detector bin
print(bins[2][11], bins[2][12])  # tally value, relative error
print(det["DET1E"][0])  # energy bin: [lower, upper, midpoint]
```

The bin-grid arrays keep Serpent's naming: `DET1E` energies, `DET1T`
times, `DET2X`/`DET2Y` spatial edges. Serpent 1 detector files (the
`sample_det.m` fixture) expose their shapes through `DET<name>_VALS` and
`DET<name>_EBINS` scalars instead.

## Fixtures

Golden-byte reference files live under `fixtures/serpent/` (`sample_*` for
Serpent 1, `serp2_*` for Serpent 2). Tests assert that Nucleide reproduces
them byte-for-byte where parity is intended.

See the
[Python API reference](../../reference/python-api.mdx#nucleideserpent) for
the full signature.
