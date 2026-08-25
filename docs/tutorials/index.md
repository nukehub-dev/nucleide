# Tutorials

Hands-on guides for Nucleide. Each tutorial is short, self-contained, and
assumes you have already installed the project (see
[getting-started.md](getting-started.md)).

## Suggested order

1. [getting-started.md](getting-started.md) — install Nucleide and verify the
   Rust and Python surfaces.
2. [parse-mcnp-output.md](parse-mcnp-output.md) — read common MCNP output files.
3. [build-materials.md](build-materials.md) — build, mix, and serialize materials.
4. [run-depletion.md](run-depletion.md) — run a CRAM depletion solve.
5. [enrichment-cascade.md](enrichment-cascade.md) — solve a multicomponent
   enrichment cascade.

## Finding more examples

- Rust unit tests live under `crates/<name>/src/` and `crates/<name>/tests/`.
- Python tests live under `tests/`.
- Golden-byte fixtures and their descriptions live under `fixtures/`.
