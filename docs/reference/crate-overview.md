# Crate Overview

Nucleide is a Cargo workspace. Each crate owns one capability area and exposes a
thin, focused API.

| Crate | Path | Responsibility |
| --- | --- | --- |
| `linalg` | `crates/linalg` | Isolation facade over the linear-algebra backend so numeric dependencies stay in one place |
| `nuclei` | `crates/nuclei` | Canonical nucid representation, element tables, naming dialects (MCNP/Serpent/FLUKA/NIST/CINDER/ALARA/SZA), particles, reaction names |
| `material` | `crates/material` | Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium loading, materials XML export |
| `mcnp-io` | `crates/mcnp-io` | xsdir, meshtal, SSW/SURFSRC, PTRAC, WWINP, MCTAL readers; material extraction from input decks; mesh-to-geometry deck generation |
| `serpent-io` | `crates/serpent-io` | `_res.m`, `_dep.m`, `_det.m` readers producing structured records |
| `fluka-io` | `crates/fluka-io` | USRBIN tally reader, MATERIAL/COMPOUND card generation |
| `enrichment` | `crates/enrichment` | Multicomponent cascade solver (numeric + assignment), SWU closed-form helpers |
| `depletion` | `crates/depletion` | CRAM matrix exponential (orders 16/48), depletion-chain XML parsing |
| `vr-tools` | `crates/vr-tools` | MAGIC weight-window generation, mesh source sampling with alias tables |
| `nucleide-bindings` | `bindings/python` | PyO3 extension module exposing the workspace to Python as `nucleide._internal` |

## Dependency rules

- `bindings/python` may depend on workspace crates.
- Workspace crates must never depend on `bindings/python` or on Python.
- `enrichment` stays independent of `material` by design.

See [Crate responsibilities](../../architecture/crate-responsibilities.md) for the dependency rationale and layering rules.
