# Nucleide Documentation

This directory contains durable documentation for Nucleide: a Rust toolkit for
nuclear-engineering workflow glue, exposed through a typed Python API. The docs
are organized by audience so you can find what you need without reading
everything else.

Start with the [project README](../README.md) for the feature list, status, and
license.

## How to use this index

- **New users** → start with [Getting started](tutorials/getting-started.md)
- **Analysts parsing legacy code output** → see the [tutorials](tutorials/index.md)
- **Readers who want the physics/math** → see the [Theory](theory/index.mdx) pages
- **Developers contributing code** → read [Local development](development/local-dev.md),
  then [Contributing](development/contributing.md)
- **Architects reviewing boundaries** → start with [Architecture overview](architecture/overview.mdx)
- **Release managers** → see [Roadmap](plan/roadmap.md) and the [changelog](../CHANGELOG.md)

## Documentation structure

### Tutorials

| Document | Purpose |
| --- | --- |
| [Tutorials](tutorials/index.md) | Tutorial index and suggested reading order |
| [Getting started](tutorials/getting-started.md) | Install Nucleide from PyPI and run your first Python snippet |
| [Python tutorials](tutorials/python/index.md) | Python tutorial index and suggested reading order |
| [Parse MCNP output](tutorials/python/parse-mcnp-output.md) | Read xsdir, meshtal, MCTAL, WWINP, PTRAC, and SSW files |
| [Build materials](tutorials/python/build-materials.md) | Build materials from formulae, mix compositions, and export XML |
| [Run depletion](tutorials/python/run-depletion.md) | Load a depletion chain and run single-step and multi-step CRAM solves |
| [Enrichment cascade](tutorials/python/enrichment-cascade.md) | Set up and solve a multicomponent enrichment cascade |
| [Activation analysis](tutorials/python/activation-analysis.md) | Read ALARA, FISPACT-II, and ORIGEN files and assemble an R2S workflow |
| [Deterministic I/O](tutorials/python/deterministic-io.md) | Read ISOTXS and flux files and write PARTISN decks |
| [Interactive tutorials](tutorials/interactive/index.mdx) | Run Nucleide in the browser through the WASM build |
| [Interactive — nuclides](tutorials/interactive/nuclides.mdx) | Nuclide identifiers and nuclear data |
| [Interactive — materials](tutorials/interactive/materials.mdx) | Formulas, fractions, mixing, and XML export |
| [Interactive — enrichment](tutorials/interactive/enrichment.mdx) | Solve MARC cascades live |
| [Interactive — depletion](tutorials/interactive/depletion.mdx) | One-step CRAM depletion |
| [Interactive — MCNP I/O](tutorials/interactive/mcnp-io.mdx) | Parse MCNP file snippets |
| [Interactive — variance reduction](tutorials/interactive/variance-reduction.mdx) | MAGIC bounds and alias-table sampling |
| [Interactive — activation](tutorials/interactive/activation.mdx) | ALARA decks/outputs, FISPACT, R2S workflows |
| [Interactive — deterministic I/O](tutorials/interactive/deterministic.mdx) | ISOTXS fluxes and PARTISN decks |

### Theory

| Document | Purpose |
| --- | --- |
| [Theory](theory/index.mdx) | Theory index and suggested reading order |
| [Depletion](theory/depletion.mdx) | Burnup matrices, the Bateman equation, CRAM, and time-series integrators |
| [Enrichment cascades](theory/enrichment.mdx) | MARC cascades, separation factors, and SWU |
| [Variance reduction](theory/variance-reduction.mdx) | MAGIC weight windows and alias-table source sampling |
| [Nuclear data](theory/nuclear-data.mdx) | Nuclide IDs, name dialects, masses, half-lives, and screening data |

### Reference

| Document | Purpose |
| --- | --- |
| [Reference](reference/index.md) | Reference index and quick links |
| [Crate overview](reference/crate-overview.mdx) | One-line responsibilities for every workspace crate |
| [Python API](reference/python-api.mdx) | Python facade overview and module map |
| [Fixtures](reference/fixtures.mdx) | Golden-byte test data index with per-file source links |

### Development

| Document | Purpose |
| --- | --- |
| [Local development](development/local-dev.md) | Toolchain, local build, test, and lint commands |
| [Contributing](development/contributing.md) | Branch workflow, commit style, and PR checklist |
| [Cross-code validation](development/validation.md) | Reproduce the PyNE/OpenMC validation harness and measured tables |

### Architecture

| Document | Purpose |
| --- | --- |
| [Architecture overview](architecture/overview.mdx) | High-level system overview and layer boundaries |
| [Crate responsibilities](architecture/crate-responsibilities.md) | Crate-level responsibilities and dependency rules |

### Plan

| Document | Purpose |
| --- | --- |
| [Roadmap](plan/roadmap.md) | Current status and upcoming priorities |
| [Decision log](plan/decision-log.md) | Architecture and process decisions with rationale |

See the project [AGENTS.md](../AGENTS.md) for ownership and contract details.
Rules for maintaining these docs live in
[Contributing](development/contributing.md#maintaining-the-docs).
