---
title: Getting Started
sidebar:
  order: 1
---

Install Nucleide and run your first Python snippet.

## Install from PyPI

```bash
pip install nucleide
```

Prebuilt wheels cover Linux, macOS, and Windows for Python >= 3.10 (abi3: one
wheel per platform serves every supported Python version).

## Optional runtime dependencies

The base wheel declares no runtime dependencies: `pip install nucleide` alone
covers everything on this page. Only two call sites reach for third-party
packages, each lazily and with a clear error if missing:
`meshtal_mesh_data(..., as_numpy=True)` imports `numpy` to return meshtal
tallies as arrays (the default `as_numpy=False` path returns plain lists), and
`write_ptrac_hdf5` imports `h5py` (and `numpy`) to export PTRAC events to
HDF5 — its error message points at `pip install h5py`, which brings numpy
along. The `test` extra (`pip install nucleide[test]`) adds pytest and numpy
for running the test suite.

## Verify the Python surface

```python
import nucleide as nuc

u = nuc.nuclei.Nuclide("U235")
print(u.nucid, u.zaid, u.serpent)

print("version:", nuc.__version__)
```

## Build from source

Use this path to try unreleased changes or to contribute. You need Rust stable
via `rustup` (workspace pins `rust-version = "1.83"`) and Python >= 3.10.

```bash
git clone https://github.com/nukehub-dev/nucleide.git
cd nucleide

cargo test --workspace   # Rust workspace
pip install maturin
maturin develop          # build + install the Python package
```

The full contributor toolchain (pytest, ruff, mypy, WASM, website) is covered
in [Local development](../development/local-dev.md).

## Next steps

- Read the [crate overview](../reference/crate-overview.mdx).
- Try [parsing an MCNP output file](python/parse-mcnp-output.md).
