---
title: Local Development
sidebar:
  order: 0
---

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

## Documentation website

The site lives in `website/` and consumes `docs/` via `nukehub-sync-docs`.
Theory pages use LaTeX math rendered by KaTeX.

```bash
cd website
npm install                 # one-time; pulls remark-math and rehype-katex
npm run sync-docs           # copy ../docs into src/content/docs
npm run build               # generate static site in dist/
npm run preview             # optional: serve locally
```

### Interactive WASM tutorials

Pages under `docs/tutorials/interactive/` run Nucleide in the browser through
`bindings/wasm/`. Build the WASM module before serving or building the site:

```bash
cd website
npm run build:wasm          # runs wasm-pack into public/wasm/
```

`website/public/wasm/` is git-ignored; regenerate it after any Rust change that
affects the WASM API.

The docs CI (`.github/workflows/docs.yml`) runs the website build plus
markdown lint and link checks.
