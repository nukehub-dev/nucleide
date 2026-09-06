---
title: Crate Overview
sidebar:
  order: 1
---

Nucleide is a Cargo workspace. Each crate owns one capability area and exposes a
thin, focused API.

<!-- GEN:crate-table:START -->

| Crate | Path | Responsibility | Depends on |
| --- | --- | --- | --- |
| `linalg` | `crates/linalg` | Isolation facade over the numeric backend (planned: faer) so the rest of the workspace never depends on it directly | - |
| `nuclei` | `crates/nuclei` | Nuclide identification, naming conventions, and reaction names | - |
| `material` | `crates/material` | Nuclear material compositions: mixing, conversions, materials XML export | `nuclei` |
| `mcnp-io` | `crates/mcnp-io` | MCNP-family file I/O: xsdir, meshtal, SSW, PTRAC readers/writers | `nuclei` |
| `serpent-io` | `crates/serpent-io` | Parsers for Serpent Monte Carlo MATLAB-style output files (\_res.m, \_dep.m, \_det.m) | - |
| `fluka-io` | `crates/fluka-io` | FLUKA Monte Carlo interface: USRBIN tally reading and MATERIAL/COMPOUND card generation | `nuclei` |
| `alara-io` | `crates/alara-io` | ALARA activation-code interop: input-deck, flux, schedule, and output glue (no solver) | `nuclei` |
| `cccc-io` | `crates/cccc-io` | CCCC binary-standard readers (ISOTXS, RTFLUX/ATFLUX, RZFLUX) and PARTISN deck writer | `nuclei` |
| `fispact-io` | `crates/fispact-io` | FISPACT-II output parser producing ALARA-compatible response frames (no solver) | `alara-io`, `nuclei` |
| `origen-io` | `crates/origen-io` | ORIGEN 2.2 TAPE readers (scoped: TAPE5 input echo, TAPE6 output, TAPE9 decay constants) | `nuclei` |
| `enrichment` | `crates/enrichment` | Multicomponent enrichment cascades and SWU analytics | `nuclei` |
| `depletion` | `crates/depletion` | CRAM matrix exponential (orders 16/48), depletion-chain XML parsing | `linalg`, `nuclei` |
| `vr-tools` | `crates/vr-tools` | MAGIC weight-window generation, mesh source sampling with alias tables | `mcnp-io` |
| `r2s` | `crates/r2s` | Rigorous two-step (R2S) shutdown-dose-rate orchestration over mcnp-io, alara-io, and vr-tools | `alara-io`, `depletion`, `material`, `mcnp-io`, `nuclei`, `vr-tools` |
| `nucleide-bindings` | `bindings/python` | PyO3 bindings exposing Nucleide to Python as nucleide.\_internal | `nuclei`, `material`, `mcnp-io`, `depletion`, `serpent-io`, `fluka-io`, `vr-tools`, `enrichment`, `alara-io`, `cccc-io`, `fispact-io`, `origen-io`, `r2s` |
| `nucleide-wasm` | `bindings/wasm` | wasm-bindgen facade exposing Nucleide core capabilities to the browser | `nuclei`, `material`, `enrichment`, `depletion`, `mcnp-io`, `alara-io`, `cccc-io`, `fispact-io`, `vr-tools`, `linalg` |

<!-- GEN:crate-table:END -->

## Dependency rules

- `bindings/python` may depend on workspace crates.
- Workspace crates must never depend on `bindings/python` or on Python.
- `enrichment` stays independent of `material` by design.

See [Crate responsibilities](../architecture/crate-responsibilities.md) for the dependency rationale and layering rules.
