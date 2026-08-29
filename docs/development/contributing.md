---
title: Contributing to Nucleide
sidebar:
  order: 1
---

Thank you for contributing to Nucleide. This document describes the workflow,
conventions, and checks expected for code changes.

## Before you start

1. Read the root `AGENTS.md` and any `AGENTS.md` files in directories you plan
   to touch.
2. Open an issue first if your change is large, architectural, or introduces
   new dependencies.
3. Make sure you can run the local development stack:
   [Local development](local-dev.md).

## Development workflow

1. Create a feature branch from `main`:

   ```bash
   git checkout main
   git pull
   git checkout -b feature/your-feature-name
   ```

2. Make your changes following the conventions below.

3. Add or update tests for new behavior.

4. Run the canonical checks:

   ```bash
   cargo fmt --check
   cargo clippy --all-targets -- -D warnings
   cargo test --workspace
   maturin develop
   pytest tests/
   ruff format --check python tests && ruff check python tests
   mypy
   ```

5. Commit with a clear message explaining what changed and why.

6. Push and open a pull request against `main`.

## Code conventions

### Rust

- Workspace pins `rust-version = "1.83"`; avoid unstable features.
- Format with `cargo fmt` (rustfmt defaults).
- Lint with `clippy --all-targets -- -D warnings`; zero warnings tolerated.
- Keep crate layering clean: workspace crates must not depend on
  `bindings/python` or on Python.
- Add tests next to changed code (`#[cfg(test)]` modules or `tests/` files).

### Python

- Target Python >= 3.10.
- Format with `ruff format`.
- Lint with `ruff check`.
- Type-check with `mypy --strict` against the `.pyi` stubs.
- The compiled `_internal` module is checked via its stub; do not add business
  logic to the pure-Python facade.

### Shell scripts

- Run `shellcheck` and `shfmt` where available.
- Prefer `#!/usr/bin/env bash`.
- Use `set -euo pipefail`.

## Documentation

Documentation is a first-class deliverable. Update docs when your change affects:

- Public Rust or Python API → `docs/reference/`
- Tutorials or worked examples → `docs/tutorials/`
- Architecture, component boundaries, or layering → `docs/architecture/`
- Developer workflow → `docs/development/`

Do not duplicate information that already lives in `README.md`, `AGENTS.md`,
generated API docs, or fixture READMEs. Link instead.

## Testing

### Rust tests

```bash
cargo test --workspace
cargo test -p nuclei
```

### Python tests

```bash
maturin develop
pytest tests/ -v
```

### Parser parity regression tests

Parsers validated against golden-byte fixtures must continue to reproduce those
fixtures byte-for-byte. If you intentionally change output, update the fixture
and all consumers in the same PR and explain why in the commit message.

## Commit messages

Use clear, imperative commit messages:

```text
Add support for custom mesh tally formats

- Adds UsrbinTally parser for Cartesian and cylindrical meshes
- Adds golden-byte fixture for cylindrical USRBIN
- Updates tutorial with cylindrical example
```

## Pull request checklist

- [ ] Branch is based on the latest `main`
- [ ] `cargo fmt --check` passes
- [ ] `cargo clippy --all-targets -- -D warnings` passes
- [ ] `cargo test --workspace` passes
- [ ] `maturin develop` succeeds
- [ ] `pytest tests/` passes
- [ ] `ruff format --check python tests && ruff check python tests` passes
- [ ] `mypy` passes
- [ ] Documentation updated for user-facing or architectural changes
- [ ] No secrets, credentials, or personal data committed
- [ ] Commit messages explain the change

## Getting help

- Open a discussion for questions.
- Open an issue for bugs or feature requests.
- Tag maintainers on security-related changes.

## License

By contributing, you agree that your contributions will be licensed under the
BSD-2-Clause license.
