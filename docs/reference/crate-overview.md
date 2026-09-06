---
title: Crate Overview
sidebar:
  order: 1
---

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
| `alara-io` | `crates/alara-io` | ALARA deck/flux/libs/output/photon/schedule glue (NOT a solver) |
| `cccc-io` | `crates/cccc-io` | ISOTXS/RTFLUX text-subset parsers + PARTISN deck writer (NOT a solver) |
| `fispact-io` | `crates/fispact-io` | FISPACT-II inventory output parser reusing the ALARA response frame (output-only) |
| `origen-io` | `crates/origen-io` | Scoped ORIGEN 2.2 TAPE5 input-echo, TAPE6 inventory, TAPE9 decay readers |
| `r2s` | `crates/r2s` | Scoped R2S workflow builder: zone-to-flux linking, schedule expansion, uniform-split photon assembly |
| `enrichment` | `crates/enrichment` | Multicomponent cascade solver (numeric), SWU closed-form helpers |
| `depletion` | `crates/depletion` | CRAM matrix exponential (orders 16/48), depletion-chain XML parsing |
| `vr-tools` | `crates/vr-tools` | MAGIC weight-window generation, mesh source sampling with alias tables |
| `nucleide-bindings` | `bindings/python` | PyO3 extension module exposing the workspace to Python as `nucleide._internal` |
| `nucleide-wasm` | `bindings/wasm` | `wasm-bindgen` crate powering the browser-based interactive tutorials |

## Dependency rules

- `bindings/python` may depend on workspace crates.
- Workspace crates must never depend on `bindings/python` or on Python.
- `enrichment` stays independent of `material` by design.

See [Crate responsibilities](../architecture/crate-responsibilities.md) for the dependency rationale and layering rules.
