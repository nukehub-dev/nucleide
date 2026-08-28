---
title: Parse MCNP Output
sidebar:
  order: 3
---

Nucleide reads several common MCNP-family files without running MCNP itself.
This tutorial shows the Python API; equivalent Rust types live in the `mcnp-io`
crate.

## xsdir

An `xsdir` file indexes the cross-section tables available to MCNP.

```python
from nucleide import read_xsdir

xs = read_xsdir("path/to/xsdir")
print(xs.datapath)
for t in xs.tables[:5]:
    print(t.name, t.zaid())
```

## meshtal

A `meshtal` file contains one or more FMESH tallies.

```python
from nucleide import read_meshtal

mt = read_meshtal("path/to/meshtal")
print(mt.version, mt.histories)

t4 = mt.tallies[4]
print(t4.dims(), t4.num_ves())
```

## MCTAL, WWINP, PTRAC, and SSW

```python
from nucleide import read_mctal, read_wwinp, read_ptrac, read_ssw

k = read_mctal("path/to/mctal")
ww = read_wwinp("path/to/wwinp")
pt = read_ptrac("path/to/ptrac")
ss = read_ssw("path/to/ssw")
```

## Writing SSW files

Use `write_ssw` to write a modified surface-source file back to disk.

```python
from nucleide import write_ssw

write_ssw(ss, "path/to/output.ssw")
```

## Fixtures

Golden-byte reference files live under `fixtures/mcnp/`. Tests assert that
Nucleide reproduces them byte-for-byte where parity is intended.
