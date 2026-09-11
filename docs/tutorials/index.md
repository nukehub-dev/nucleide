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
3. [Python tutorials](python/index.md) — hands-on Python guides for parsing,
   materials, depletion, enrichment, activation, and deterministic I/O.
4. [Parse MCNP output](python/parse-mcnp-output.md) — read common MCNP output files.
5. [Build materials](python/build-materials.md) — build, mix, and serialize materials.
6. [Run depletion](python/run-depletion.md) — run a CRAM depletion solve.
7. [Enrichment cascade](python/enrichment-cascade.md) — solve a multicomponent
   enrichment cascade.
8. [Activation analysis](python/activation-analysis.md) — read ALARA, FISPACT-II,
   and ORIGEN files and assemble an R2S workflow.
9. [Deterministic I/O](python/deterministic-io.md) — read ISOTXS and flux files
   and write PARTISN decks.
10. [Emit code cards](python/emit-cards.md) — render one material to MCNP,
    Serpent, FLUKA, ALARA, and PARTISN cards with a mass-drift report.

## Finding more examples

- More output-parsing walkthroughs: [Parse Serpent output](python/parse-serpent-output.md)
  and [Parse FLUKA output](python/parse-fluka-output.md) mirror the MCNP one.
- Rust unit tests live in inline `#[cfg(test)]` modules under `crates/<name>/src/`.
- Python tests live under `tests/`.
- Golden-byte fixtures and their descriptions live under `fixtures/`.
