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
| [Parse MCNP output](tutorials/parse-mcnp-output.md) | Read xsdir, meshtal, MCTAL, WWINP, PTRAC, and SSW files |
| [Build materials](tutorials/build-materials.md) | Build materials from formulae, mix compositions, and export XML |
| [Run depletion](tutorials/run-depletion.md) | Load a depletion chain and run a short CRAM solve |
| [Enrichment cascade](tutorials/enrichment-cascade.md) | Set up and solve a multicomponent enrichment cascade |
| [Interactive tutorials](tutorials/interactive/index.mdx) | Run Nucleide in the browser through the WASM build |
| [Interactive — nuclides](tutorials/interactive/nuclides.mdx) | Nuclide identifiers and nuclear data |
| [Interactive — materials](tutorials/interactive/materials.mdx) | Formulas, fractions, mixing, and XML export |
| [Interactive — enrichment](tutorials/interactive/enrichment.mdx) | Solve MARC cascades live |
| [Interactive — depletion](tutorials/interactive/depletion.mdx) | One-step CRAM depletion |
| [Interactive — MCNP I/O](tutorials/interactive/mcnp-io.mdx) | Parse MCNP file snippets |
| [Interactive — variance reduction](tutorials/interactive/variance-reduction.mdx) | MAGIC bounds and alias-table sampling |

### Theory

| Document | Purpose |
| --- | --- |
| [Theory](theory/index.mdx) | Theory index and suggested reading order |
| [Depletion](theory/depletion.mdx) | Burnup matrices, the Bateman equation, and CRAM |
| [Enrichment cascades](theory/enrichment.mdx) | MARC cascades, separation factors, and SWU |
| [Variance reduction](theory/variance-reduction.mdx) | MAGIC weight windows and alias-table source sampling |
| [Nuclear data](theory/nuclear-data.mdx) | Nuclide IDs, name dialects, masses, and half-lives |

### Reference

| Document | Purpose |
| --- | --- |
| [Reference](reference/index.md) | Reference index and quick links |
| [Crate overview](reference/crate-overview.md) | One-line responsibilities for every workspace crate |
| [Python API](reference/python-api.md) | Python facade overview and module map |

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

## Maintenance rules

1. **Keep docs in sync with code.** A PR that changes a public API, parser
   output, or crate boundary must update the matching tutorial or reference page.
2. **Prefer deletion over stale historical notes.** If a section no longer
   reflects current behavior, delete it or move it to an explicit "Historical"
   appendix with a removal date.
3. **Do not duplicate details that live elsewhere.** Link to the
   [project README](../README.md), API stubs, and fixture READMEs instead of
   copying them.
4. **Use relative links.** Internal links must be relative so documentation
   stays usable offline and in branches.
5. **Use `.mdx` for component-heavy pages.** Pages that use `@nukehub/docs-kit`
   shortcodes such as `<Mermaid>`, `<Callout>`, or `<DataTable>` must have an
   `.mdx` extension. Plain `.md` is fine for prose-only pages.
6. **Every page sets `title` and `sidebar.order` in frontmatter.** Without an
   explicit order the site sidebar falls back to alphabetical sorting. Keep the
   order values aligned with the reading order in the section index tables, and
   do not repeat the title as an in-body `#` heading (the site renders the
   frontmatter title as the page heading).

See the project [AGENTS.md](../AGENTS.md) for ownership and contract details.
