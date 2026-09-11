---
title: Parse FLUKA Output
sidebar:
  order: 3
---

Nucleide reads FLUKA USRBIN mesh tallies from `.lis` output files without
running FLUKA itself. This tutorial shows the Python API; equivalent Rust
types live in the `nucleide-fluka-io` crate.

## USRBIN tallies

`read_usrbin` parses every USRBIN block in a `.lis` file into one
`UsrbinTally` each: detector metadata, mesh bounds, track-length data, and
percentage errors as plain lists — no meshing layer.

```python
from nucleide.fluka import read_usrbin

tallies = read_usrbin("fixtures/fluka/fluka_usrbin_single.lis")
t = tallies[0]
print(t.name, t.particle)

nx, ny, nz = t.dims()
print(nx * ny * nz == len(t.data) == len(t.error))
print(t.x_bounds)  # nx + 1 bin edges
```

`data` and `error` are flat in the file's `A(ix,iy,iz)` order (first index
fastest), so `len(data) == nx * ny * nz` and each of `x_bounds`,
`y_bounds`, `z_bounds` carries one more edge than bins.

## Multiple detectors per file

A `.lis` file can hold several USRBIN blocks; each becomes one tally. The
degenerate fixture shows that per-axis bin counts need not be ordered
largest-first:

```python
tallies = read_usrbin("fixtures/fluka/fluka_usrbin_multiple.lis")
print([t.name for t in tallies])
```

## Fixtures

Golden-byte reference files live under `fixtures/fluka/` (single, multiple,
and degenerate binning cases plus the generating input deck). Tests assert
that Nucleide reproduces them byte-for-byte where parity is intended.

See the
[Python API reference](../../reference/python-api.mdx#nucleidefluka) for
the full signatures.
