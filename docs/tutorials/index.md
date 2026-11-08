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
3. [Python tutorials](python/index.md) — then follow its suggested order:
   twenty-two short guides covering parsing, materials, solvers, emission,
   CSG translation, and UQ.

## Finding more examples

- More output-parsing walkthroughs: [Parse Serpent output](python/parse-serpent-output.md)
  and [Parse FLUKA output](python/parse-fluka-output.md) mirror the MCNP one.
- Rust unit tests live in inline `#[cfg(test)]` modules under `crates/<name>/src/`.
- Python tests live under `tests/`.
- Sample files used by the test suite (with per-file descriptions) live under
  `fixtures/`; see the
  [sample data files](../reference/fixtures.mdx) index for the full list.
