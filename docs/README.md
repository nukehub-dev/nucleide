# Nucleide Documentation

This directory contains durable documentation for Nucleide: a Rust toolkit for
nuclear-engineering workflow glue, exposed through a typed Python API. The docs
are organized by audience so you can find what you need without reading
everything else.

Start with the project overview in [`../README.md`](../README.md) for the
feature list, status, and license.

## How to use this index

- **New users** → start with [tutorials/getting-started.md](tutorials/getting-started.md)
- **Analysts parsing legacy code output** → see the tutorials under
  [tutorials/](tutorials/index.md)
- **Developers contributing code** → read
  [development/local-dev.md](development/local-dev.md), then
  [development/contributing.md](development/contributing.md)
- **Architects reviewing boundaries** → start with
  [architecture/overview.md](architecture/overview.md)
- **Release managers** → see [plan/roadmap.md](plan/roadmap.md) and
  [`../CHANGELOG.md`](../CHANGELOG.md)

## Documentation structure

### Tutorials

| Document | Purpose |
|---|---|
| [tutorials/index.md](tutorials/index.md) | Tutorial index and suggested reading order |
| [tutorials/getting-started.md](tutorials/getting-started.md) | Install Nucleide and run your first Rust or Python snippet |
| [tutorials/parse-mcnp-output.md](tutorials/parse-mcnp-output.md) | Read xsdir, meshtal, MCTAL, WWINP, PTRAC, and SSW files |
| [tutorials/build-materials.md](tutorials/build-materials.md) | Build materials from formulae, mix compositions, and export XML |
| [tutorials/run-depletion.md](tutorials/run-depletion.md) | Load a depletion chain and run a short CRAM solve |
| [tutorials/enrichment-cascade.md](tutorials/enrichment-cascade.md) | Set up and solve a multicomponent enrichment cascade |

### Reference

| Document | Purpose |
|---|---|
| [reference/index.md](reference/index.md) | Reference index and quick links |
| [reference/crate-overview.md](reference/crate-overview.md) | One-line responsibilities for every workspace crate |
| [reference/python-api.md](reference/python-api.md) | Python facade overview and module map |

### Development

| Document | Purpose |
|---|---|
| [development/local-dev.md](development/local-dev.md) | Toolchain, local build, test, and lint commands |
| [development/contributing.md](development/contributing.md) | Branch workflow, commit style, and PR checklist |

### Architecture

| Document | Purpose |
|---|---|
| [architecture/overview.md](architecture/overview.md) | High-level system overview and layer boundaries |
| [architecture/crate-responsibilities.md](architecture/crate-responsibilities.md) | Crate-level responsibilities and dependency rules |

### Plan

| Document | Purpose |
|---|---|
| [plan/roadmap.md](plan/roadmap.md) | Current status and upcoming priorities |
| [plan/decision-log.md](plan/decision-log.md) | Architecture and process decisions with rationale |

## Maintenance rules

1. **Keep docs in sync with code.** A PR that changes a public API, parser
   output, or crate boundary must update the matching tutorial or reference page.
2. **Prefer deletion over stale historical notes.** If a section no longer
   reflects current behavior, delete it or move it to an explicit "Historical"
   appendix with a removal date.
3. **Do not duplicate details that live elsewhere.** Link to
   [`../README.md`](../README.md), API stubs, and fixture READMEs instead of
   copying them.
4. **Use relative links.** Internal links must be relative so documentation
   stays usable offline and in branches.

See [`../AGENTS.md`](../AGENTS.md) for ownership and contract details.
