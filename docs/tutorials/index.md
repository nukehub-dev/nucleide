---
title: Tutorials
sidebar:
  order: 0
---

Hands-on guides for Nucleide. Each tutorial is short, self-contained, and
assumes you have already installed the project (see
[Getting started](getting-started.md)).

## Suggested order

1. [Getting started](getting-started.md) — install Nucleide and verify the
   Rust and Python surfaces.
2. [Interactive tutorials](interactive/index.mdx) — run Nucleide in your browser,
   no installation required.
3. [Parse MCNP output](parse-mcnp-output.md) — read common MCNP output files.
4. [Build materials](build-materials.md) — build, mix, and serialize materials.
5. [Run depletion](run-depletion.md) — run a CRAM depletion solve.
6. [Enrichment cascade](enrichment-cascade.md) — solve a multicomponent
   enrichment cascade.

## Finding more examples

- Rust unit tests live in inline `#[cfg(test)]` modules under `crates/<name>/src/`.
- Python tests live under `tests/`.
- Golden-byte fixtures and their descriptions live under `fixtures/`.
