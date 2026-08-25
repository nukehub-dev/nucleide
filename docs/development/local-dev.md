# Local Development

Set up a Nucleide development environment and run the canonical checks.

## Required tooling

Install once before making changes:

- **Rust** stable via `rustup` (workspace pins `rust-version = "1.83"`).
- **Python** >= 3.10 with:
  `pip install maturin pytest pytest-cov ruff mypy`

The compiled extension must be rebuilt (`maturin develop`) after any Rust change
before running Python tests.

## Rust side

```bash
cargo test                       # workspace unit tests
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

## Python side

```bash
pip install maturin pytest pytest-cov ruff mypy
maturin develop                  # build + install into current venv
pytest tests/
```

## Lint, format, and type-check

```bash
ruff format python tests
ruff check python tests
mypy                             # strict; stubs in *.pyi
```

## Coverage

- Python: pytest runs with coverage enabled by default (`pyproject.toml`
  addopts); keep new code covered.
- Rust: `cargo llvm-cov --workspace` when touching numeric kernels or parsers;
  CI uploads an lcov report.

## Before committing

Run these from the repo root. They are the canonical "did I break anything"
checks:

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --workspace
maturin develop
pytest tests/
ruff format --check python tests && ruff check python tests
mypy
```

CI (`.github/workflows/ci.yml`) runs exactly these checks plus coverage; keep
it green.
