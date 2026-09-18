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

xs = read_xsdir("fixtures/mcnp/xsdir/dummy_xsdir")
print(xs.datapath)
for t in xs.tables[:5]:
    print(t.name, t.zaid())
```

## meshtal

A `meshtal` file contains one or more FMESH tallies.

```python
from nucleide.mcnp import read_meshtal

mt = read_meshtal("fixtures/mcnp/meshtal/mcnp_meshtal_single_meshtal.txt")
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

k = read_mctal("fixtures/mcnp/mctal/synthetic_kcode5.mctal")
ww = read_wwinp("fixtures/mcnp/wwinp/mcnp_wwinp_wwinp_np.txt")
pt = read_ptrac("fixtures/mcnp/ptrac/mcnp_ptrac_i4_little.ptrac")
ss = read_ssw("fixtures/mcnp/ssw/mcnp_surfsrc_onetrack.w")
```

`read_mctal` returns the header, the standard tally bodies (`tallies`, each
with per-card `{count, values, variant, flag}` bins plus `(value, rel_error)`
`vals` pairs, an optional `tfc` fluctuation-chart block, and a `total`), the
rectangular/cylindrical/spherical mesh tallies (`mesh_tallies`, with the
`ni`/`nj`/`nk` mesh counts and `cora`/`corb`/`corc` bounds), and the `kcode`
cycles. Total/cumulative card variants (`ut`/`uc`/…) and the `c`/`e`/`t`
third-token flags are stored verbatim, never interpreted. Radiograph and
point-detector specials plus perturbation bodies are not yet supported —
reading them raises a clear error instead of parsing silently. `tally_vals_array(number)` exposes one tally's
pairs as a NumPy `(n_pairs, 2)` array (`mesh_tally_vals_array` covers mesh
tallies).

## Writing SSW files

Use `write_ssw` to write a modified surface-source file back to disk
(writes scratch files, so this snippet is excluded from code execution;
`tests/test_mcnp.py` covers the writer):

<!-- code-test: skip -->
```python
from nucleide.mcnp import write_ssw

write_ssw(ss, "path/to/output.ssw")
```

## Input decks: parse, edit, write back

`DeckProblem` parses a full input deck (message/title/cell/surface/data
blocks) into typed cell/surface cards plus material and passthrough data
cards, and writes it back byte-identical when unedited — only cards touched
through the setters re-render canonically (paths are placeholders, so this
snippet is excluded from code execution):

<!-- code-test: skip -->
```python
from nucleide.mcnp import read_deck

deck = read_deck("path/to/model.i")
print(deck.title, [c["num"] for c in deck.cells])
deck.set_cell_density(1, -10.0)
open("path/to/model_edited.i", "w").write(deck.dumps())
```

## SDEF source cards

`parse_sdef` reads one legacy `SDEF` source definition plus its `SIn`/`SPn`/
`SBn` distribution cards into typed dicts. The accepted keyword set is
`POS`/`CELL`/`SURF`/`VEC`/`DIR`/`ERG`/`NRM`/`PAR`/`WGT`/`TME`, and only the
discrete distribution forms (`SIn L`, `SPn D`, `SBn D`) parse — non-discrete
options, dangling distribution references, and length mismatches all raise a
`ValueError`:

```python
from nucleide.mcnp import parse_sdef

parsed = parse_sdef("SDEF POS=0 0 0\n     ERG=D1\nSI1 L 0.662 1.17\nSP1 D 0.5 0.5")
print(parsed["pos"], parsed["erg"])  # 0 0 0, D1
print(parsed["distributions"][0]["si"])  # 0.662 1.17
print(parsed["card"])  # canonical re-render of the full card
```

Unknown keywords are not dropped silently: they land in the `ignored` list as
drift notes — entries in the drift report, the running list of every judgment
call the parser made. A deck read through `read_deck`/`parse_deck` exposes the
same view as the `sdef` property (`None` when the deck has no source card),
and cards emitted by the spectroscopy SDEF source round-trip byte-identical:

```python
from nucleide import mcnp, spectroscopy

bins, card = spectroscopy.sdef_decay_source([(0.662, 2.0), (1.33, 1.0)])
assert mcnp.parse_sdef(card)["card"] == card

deck = mcnp.parse_deck("msg\ntitle\n1 0 -1\n\n1 so 1.0\n\n" + card + "\n")
print(deck.sdef["erg"])  # D1
```

See [Emit an SDEF decay source](run-spectroscopy.md#emit-an-sdef-decay-source)
for the emitter side.

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

## Sample files

Sample reference files live under `fixtures/mcnp/`. The test suite checks
that Nucleide reproduces them byte-for-byte where exact reproduction is
intended; see the
[sample data files](../../reference/fixtures.mdx) index for the full list.

## ENDL, SSW combining, PTRAC export, and mesh data

The same module covers the remaining legacy touchpoints (all Python-side
except where noted):

```python
import tempfile
from nucleide.mcnp import (
    combine_ssw_files,
    meshtal_mesh_data,
    ptrac_event_rows,
    read_endl,
    write_ptrac_hdf5,
)

tmp = tempfile.mkdtemp(prefix="nucleide-mcnp-")
ref_w = "fixtures/mcnp/ssw/mcnp_surfsrc_onetrack.w"

lib = read_endl("fixtures/endl/synthetic_eedl.txt")
print(lib.nuclides())  # [820000000]
print(lib.get_rx(820000000, 9, 10, 0)[:1])  # integrated table rows

combine_ssw_files(f"{tmp}/merged.w", [ref_w, ref_w])

rows = ptrac_event_rows("fixtures/mcnp/ptrac/mcnp_ptrac_i4_little.ptrac")  # 19 PtracEvent columns, no HDF5 needed
```

The HDF5 writer needs your `h5py` (bytes never asserted here; covered by
`tests/test_mcnp_io.py`), so it is excluded from code execution:

<!-- code-test: skip -->
```python
from nucleide.mcnp import write_ptrac_hdf5

write_ptrac_hdf5("fixtures/mcnp/ptrac/mcnp_ptrac_i4_little.ptrac", "run.h5")
```python
from nucleide.mcnp import meshtal_mesh_data

mesh = meshtal_mesh_data("fixtures/mcnp/meshtal/mcnp_meshtal_single_meshtal.txt", as_numpy=True)
print(mesh["tallies"][4]["result"].shape)  # (ve, groups)
```

`meshtal_mesh_data` returns plain dicts/arrays only — MOAB tagging stays
caller-side. See the
[Python API reference](../../reference/python-api.mdx#nucleidemcnp) for the
full signatures.
