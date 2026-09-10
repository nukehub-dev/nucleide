---
title: Parse MCNP Output
sidebar:
  order: 1
---

Nucleide reads several common MCNP-family files without running MCNP itself.
This tutorial shows the Python API; equivalent Rust types live in the `mcnp-io`
crate.

## xsdir

An `xsdir` file indexes the cross-section tables available to MCNP.

```python
from nucleide.mcnp import read_xsdir

xs = read_xsdir("path/to/xsdir")
print(xs.datapath)
for t in xs.tables[:5]:
    print(t.name, t.zaid())
```

## meshtal

A `meshtal` file contains one or more FMESH tallies.

```python
from nucleide.mcnp import read_meshtal

mt = read_meshtal("path/to/meshtal")
print(mt.version, mt.histories)

t4 = mt.tallies[4]
print(t4.dims(), t4.num_ves())
```

`result_array()` / `totals_array()` expose the same tables as NumPy
arrays (requires NumPy at runtime; `to_list()` / `totals_list()` stay
the NumPy-free plain-copy path). Each returns a `(result, rel_error)`
pair: `result_array()` yields owned writable C-order float64
`(ve, group)` arrays with `ve = (i * ny + j) * nz + k`, while
`totals_array()` yields `(num_ves,)` totals pairs:

```python
import numpy as np

from nucleide.mcnp import read_meshtal

mt = read_meshtal("fixtures/mcnp/meshtal/mcnp_meshtal_single_meshtal.txt")
t = mt.tallies[4]
result, rel_err = t.result_array()
print(result.shape, result.dtype, result.flags["C_CONTIGUOUS"])
totals, totals_err = t.totals_array()
print(totals.shape, totals.dtype)
```

Two hardening notes for deck-adjacent work: cell-`FILL` universe
matrices are capped at 1 000 000 cells (hostile `i:j` ranges fail with
a cap error before any allocation — larger lattices must be built
programmatically), and raw nucid integers outside the validated
`(Z, A, state)` domain render the diagnostic `Z{z}A{a}[m{s}]` fallback
in `to_name` / ARMI labels instead of panicking (validate untrusted
integers with `try_from_nucid` / `is_valid` on the Rust side).

## MCTAL, WWINP, PTRAC, and SSW

```python
from nucleide.mcnp import read_mctal, read_ptrac, read_ssw, read_wwinp

k = read_mctal("path/to/mctal")
ww = read_wwinp("path/to/wwinp")
pt = read_ptrac("path/to/ptrac")
ss = read_ssw("path/to/ssw")
```

## Writing SSW files

Use `write_ssw` to write a modified surface-source file back to disk.

```python
from nucleide.mcnp import write_ssw

write_ssw(ss, "path/to/output.ssw")
```

## Input decks: parse, edit, write back

`DeckProblem` parses a full input deck (message/title/cell/surface/data
blocks) into typed cell/surface cards plus material and passthrough data
cards, and writes it back byte-identical when unedited — only cards touched
through the setters re-render canonically:

```python
from nucleide.mcnp import read_deck

deck = read_deck("path/to/model.i")
print(deck.title, [c["num"] for c in deck.cells])
deck.set_cell_density(1, -10.0)
open("path/to/model_edited.i", "w").write(deck.dumps())
```

## L3 semantics: universes, lattices, tallies, and validation

`DeckProblem` also exposes typed semantic views over the raw cards —
`MODE` particles, `TRn` transforms, auto-created universes (`U`/`-U`),
`LAT`/`FILL` assignments, importances, volumes, and `F`/`FM`/`E` tallies —
plus `validate()`, which centralizes duplicate-number, dangling-link,
redundant-definition, and write-time-state checks (lattice/fill
cross-checks are nucleide-defined: `LAT` without `FILL`, and a `FILL`
matrix without `LAT`, are errors). Particle/mode mismatches are notes,
not errors:

```python
from nucleide.mcnp import read_deck

deck = read_deck("fixtures/mcnp/inp/deck_l3.txt")
print(deck.mode)  # {'particles': 'N P'}
print([(u["number"], u["cells"]) for u in deck.universes])
print(deck.fills[0]["universes"], deck.tallies[0]["e_bins"])
deck.validate()
print(deck.validation_notes())
deck.set_cell_universe(4, 5, True)  # writes U=-5
deck.validate()
```

## Fixtures

Golden-byte reference files live under `fixtures/mcnp/`. Tests assert that
Nucleide reproduces them byte-for-byte where parity is intended.
