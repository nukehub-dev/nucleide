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

## API stability (pre-1.0 deprecation policy)

All workspace crates are `0.x`, so minor releases may break the Rust and
Python APIs. Breaking changes are signaled, never silent:

1. Record every breaking change under `CHANGELOG.md → [Unreleased] →
   Changed` in the same PR that makes it.
2. Prefer deprecating over removing: `#[deprecated(since = "…", note =
   "…")]` on Rust items (removed no earlier than the next minor) and
   `warnings.warn(..., DeprecationWarning)` on the Python facade for the
   same window.
3. Pure additions (new functions, modules, error variants, struct fields)
   need no deprecation period, but error enums stay exhaustive within a
   minor line so downstream `match`es keep compiling — a new variant is a
   breaking change and gets a changelog entry like any other.
4. Stabilized crates carry `#![warn(missing_docs)]`: every public item
   ships documented, types with invariants construct through validating
   constructors (record-style structs may keep public fields, validated at
   use), and fallible paths return `Result` — no user-reachable panics
   (private `expect`/`unwrap` on validated invariants must cite why the
   input cannot occur).

The 0.5.0 API freeze covered (`nucleide-kinetics`,
`nucleide-spectroscopy`, the `nucleide-mcnp-io` `endl`/`fortran` modules,
and the numpy bridges): doc-comment and attribute edits only, no signature
or behavior changes.

## Documentation

Documentation is a first-class deliverable. Update docs when your change affects:

- Public Rust or Python API → `docs/reference/`
- Tutorials or worked examples → `docs/tutorials/`
- Architecture, component boundaries, or layering → `docs/architecture/`
- Developer workflow → `docs/development/`

Do not duplicate information that already lives in `README.md`, `AGENTS.md`,
generated API docs, or fixture READMEs. Link instead.

### Maintaining the docs

1. **Keep docs in sync with code.** A PR that changes a public API, parser
   output, or crate boundary must update the matching tutorial or reference page.
2. **Regenerate the reference tables.** After changing `python/nucleide/_internal.pyi`,
   the `python/nucleide/*.py` facades, crate manifests/descriptions, or adding
   fixtures (document the area README), run
   `python3 scripts/gen-reference.py --write` and commit the result (CI enforces
   freshness).
3. **Prefer deletion over stale historical notes.** If a section no longer
   reflects current behavior, delete it or move it to an explicit "Historical"
   appendix with a removal date.
4. **Do not duplicate details that live elsewhere.** Link to the
   [project README](../../README.md), API stubs, and fixture READMEs instead of
   copying them.
5. **Use relative links.** Internal links must be relative so documentation
   stays usable offline and in branches.
6. **Use `.mdx` for component-heavy pages.** Pages that use `@nukehub/docs-kit`
   shortcodes such as `<Mermaid>`, `<Callout>`, or `<DataTable>` must have an
   `.mdx` extension. Plain `.md` is fine for prose-only pages.
7. **Every page sets `title` and `sidebar.order` in frontmatter.** Without an
   explicit order the site sidebar falls back to alphabetical sorting. Keep the
   order values aligned with the reading order in the section index tables, and
   do not repeat the title as an in-body `#` heading (the site renders the
   frontmatter title as the page heading).

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
