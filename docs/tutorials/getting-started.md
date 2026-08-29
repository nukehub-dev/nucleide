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

- Read the [crate overview](../reference/crate-overview.md).
- Try [parsing an MCNP output file](parse-mcnp-output.md).
