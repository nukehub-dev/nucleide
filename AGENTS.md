# Nuke Agent Doc (NAD) Framework

## Purpose

Binding work contract for AI agents and human contributors working on Nucleide.

## Ownership

This root `AGENTS.md` owns the NAD hierarchy, project-wide workflow rules, and
cross-domain standards. Domain-specific guidance lives in child `AGENTS.md`
files listed in the Child NAD Index.

## NAD Core Contract

- `AGENTS.md` files are binding work contracts for their subtrees.
- Work products, source materials, instructions, records, assets, and durable
  docs must stay understandable from the nearest applicable `AGENTS.md` plus
  every parent `AGENTS.md` above it.

### Read Before Editing

1. Read this root `AGENTS.md`.
2. Identify every file or folder you expect to touch.
3. Walk from the repository root to each target path.
4. Read every `AGENTS.md` found along each route.
5. If a parent `AGENTS.md` lists a child `AGENTS.md` whose scope contains the
   path, read that child and continue from there.
6. Use the nearest `AGENTS.md` as the local contract and parent docs for
   repo-wide rules.
7. If docs conflict, the closer doc controls local work details, but no child
   doc may weaken NAD.

### Update After Editing

Every meaningful change requires a NAD pass before the task is done.

Update the closest owning `AGENTS.md` when a change affects:

- purpose, scope, ownership, or responsibilities
- durable structure, contracts, workflows, or operating rules
- required inputs, outputs, permissions, constraints, side effects, or artifacts
- user preferences about behavior, communication, process, organization, or quality
- `AGENTS.md` creation, deletion, move, rename, or index contents

Update parent docs when parent-level structure, ownership, workflow, or child
index changes. Update child docs when parent changes alter local rules. Remove
stale or contradictory text immediately. Small edits that do not change
behavior or contracts may leave docs unchanged, but the NAD pass still must
happen.

### Docs Pass

`AGENTS.md` updates do not cover user/dev documentation. In the same change,
also update `README.md` when a change alters user-visible behavior — features,
public API (Rust or Python), file formats, install/dev workflows.

Internal refactors, bug fixes with no behavior change, and test-only work need
neither.

## Hierarchy

- Root `AGENTS.md` is the NAD rail: project-wide instructions, global
  preferences, durable workflow rules, and the top-level Child NAD Index.
- Child `AGENTS.md` files own domain-specific instructions and their own Child
  NAD Index.
- Each parent explains what its direct children cover and what stays owned by
  the parent.
- The closer a doc is to the work, the more specific and practical it must be.

## Child Doc Shape

Create a child `AGENTS.md` when a folder becomes a durable boundary with its
own purpose, rules, responsibilities, workflow, materials, or quality standards.

Default section order:

- Purpose
- Ownership
- Local Contracts
- Work Guidance
- Verification
- Child NAD Index

## Style

- Keep docs concise, current, and operational.
- Document stable contracts, not diary entries.
- Put broad rules in parent docs and concrete details in child docs.
- Prefer direct bullets with explicit names.
- Do not duplicate rules across many files unless each scope needs a local version.
- Delete stale notes instead of explaining history.
- Trim obvious statements, repeated rules, misplaced detail, and warnings for
  risks that no longer exist.

## Closeout

1. Re-check changed paths against the NAD chain.
2. Update nearest owning docs and any affected parents or children.
3. Refresh every affected Child NAD Index.
4. Remove stale or contradictory text.
5. Run existing verification when relevant.
6. Report any docs intentionally left unchanged and why.

## User Preferences

When the user requests a durable behavior change, record it here or in the
relevant child `AGENTS.md`.

---

## Nucleide Project Guidance

## Required tooling

Install once before making changes:

- **Rust** stable via `rustup` (workspace pins `rust-version = "1.83"`).
- **Python** >= 3.10 with:
  `pip install maturin pytest pytest-cov ruff mypy`
- The compiled extension must be rebuilt (`maturin develop`) after any Rust
  change before running Python tests.
- **WASM** (only for interactive tutorials):
  `rustup target add wasm32-unknown-unknown` and `cargo install wasm-pack`.
- **Website E2E tests** (only when touching interactive demos):
  `cd website && npx playwright install chromium` after `npm install`.

## Before committing

Run these from the repo root. They are the canonical "did I break anything"
checks:

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings   # zero warnings tolerated
cargo test --workspace
maturin develop                              # rebuild Python extension
pytest tests/
ruff format --check python tests && ruff check python tests
mypy                                         # strict; stubs in *.pyi
cd website && npm run build:wasm && npm run check && npm run build && npm run test:e2e:ci
```

Notes:

- Formatter/linter configs live in `pyproject.toml` (ruff line-length 100,
  mypy strict). Rust formatting is rustfmt defaults.
- CI (`.github/workflows/ci.yml`) runs exactly these checks plus coverage;
  keep it green.

## Coverage

- Python: pytest runs with coverage enabled by default (`pyproject.toml`
  addopts); keep new code covered.
- Rust: `cargo llvm-cov --workspace` when touching numeric kernels or parsers;
  CI uploads an lcov report.

## Architecture pointer

High-level layout; see the Child NAD Index below for domain-specific details.

- `crates/` — Rust workspace members (one crate per capability area):
  `nuclei`, `material`, `mcnp-io`, `serpent-io`, `fluka-io`, `vr-tools`,
  `enrichment`, `depletion`, `linalg`.
- `bindings/python/` — PyO3 crate exposing `nucleide._internal`; thin facade,
  no business logic.
- `bindings/wasm/` — `wasm-bindgen` crate that lets tutorials run Nucleide in
  the browser without Python.
- `python/nucleide/` — typed pure-Python package surface (`.pyi` stubs +
  `py.typed`); re-exports only.
- `fixtures/` — golden-byte test data (see `fixtures/README.md`).
- `tests/` — Python-side tests (run after `maturin develop`).
- `.research/` — local-only working notes; git-ignored, never referenced by
  committed docs or code.

## Documentation

The `docs/` tree owns durable user and contributor documentation.

- Audience-based layout: `tutorials/`, `reference/`, `development/`, `architecture/`, `plan/`.
- Index and maintenance rules live in `docs/README.md`.
- Internal links must be relative and must not duplicate details already in
  `README.md`, `AGENTS.md`, generated API stubs, or fixture READMEs.
- Docs changes trigger `.github/workflows/docs.yml` for markdown lint and link
  checking.

The `website/` directory is the Astro-based documentation website. It consumes
`@nukehub/docs-kit`, pulls content from `../docs/`, and deploys to GitHub Pages via
`.github/workflows/docs-deploy.yml`. Keep site branding and routing in
`website/src/data/` and `website/src/pages/`; shared UI lives in the kit.

Interactive tutorials use the WASM build. Run `npm run build:wasm` from
`website/` to regenerate `website/public/wasm/` (git-ignored) after any Rust
change that affects `bindings/wasm`. Changes that affect interactive demos
must pass the E2E smoke tests (`npm run test:e2e:ci` from `website/`).

## Release workflow

Git tags of the form `vX.Y.Z` trigger `.github/workflows/release.yml`, which:

1. Runs the canonical verify checks.
2. Builds and publishes Python wheels for Linux, macOS, and Windows via
   `maturin publish` (requires `PYPI_API_TOKEN`).
3. Optionally publishes workspace crates to crates.io in dependency order
   (requires `CARGO_REGISTRY_TOKEN`; enabled via a `workflow_dispatch` input).
4. Drafts a GitHub release with the wheels and the matching `CHANGELOG.md`
   section.

Use `scripts/bump-version.sh X.Y.Z` to bump the workspace version and stamp
`CHANGELOG.md` before tagging.

## Common pitfalls

- **Never edit generated/build artifacts**: `target/`, `*.so`, `__pycache__/`,
  `coverage.xml`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`.
- **Golden fixtures are byte-exact oracles**: changing one requires a written
  reason in the commit message and updated assertions in the same change.
- **Keep the layering**: `bindings/python` may depend on workspace crates;
  workspace crates must never depend on bindings or on Python.
  `enrichment` stays independent of `material` by design.
- **Parser parity quirks are intentional**: where a reader reproduces legacy
  output byte-for-byte (padding, exponent formats, sentinel values), do not
  "fix" it without updating fixtures and their consumers together.
- **Commit `Cargo.lock`** (binary workspace); do not hand-edit it.
- **Version bumps** happen in `[workspace.package] version` in the root
  `Cargo.toml`; crates inherit via `version.workspace = true`. Run
  `scripts/bump-version.sh X.Y.Z` to update the workspace version and stamp
  `CHANGELOG.md`; then commit, tag `vX.Y.Z`, and push.

## Child NAD Index

- `website/AGENTS.md` — website build, preview, sync, and E2E test workflow.

Create additional child `AGENTS.md` files under `crates/<name>/` or other
folders once they grow their own durable contracts (e.g. fixture policy
details, parser-parity rules), and list them here.
