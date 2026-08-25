# Getting Started

Install Nucleide, build the Rust workspace, and exercise the Python bindings.

## What you need

- Rust stable via `rustup` (workspace pins `rust-version = "1.83"`).
- Python >= 3.10.
- `pip install maturin pytest pytest-cov ruff mypy`.

## Clone and build

```bash
git clone https://github.com/nukehub-dev/nucleide.git
cd nucleide

cargo build --workspace
cargo test --workspace
```

## Install the Python package locally

```bash
maturin develop
```

This compiles the PyO3 extension and installs the `nucleide` package into your
current virtual environment.

## Verify the Python surface

```python
import nucleide as nuc

u = nuc.Nuclide("U235")
print(u.nucid, u.zaid, u.serpent)

print("version:", nuc.__version__)
```

## Next steps

- Read the [crate overview](../reference/crate-overview.md).
- Try [parsing an MCNP output file](parse-mcnp-output.md).
