---
title: Crate Responsibilities
sidebar:
  order: 2
---

Each Nucleide crate has a single responsibility and a well-defined place in the
workspace dependency graph.

## Foundation crates

### `nucleide-linalg`

Isolation facade over the linear-algebra backend. Today it pulls in `faer`,
`num-complex`, and `roxmltree` so numeric dependencies do not leak into other
crates. Other workspace crates depend on `nucleide-linalg`, not on the backend directly.

### `nucleide-nuclei`

Canonical nuclide identification. Owns:

- `NuclideId`: compact `u32` representation compatible with legacy integer ids.
- Element symbol/number tables.
- Name parsing (`U235`, `Am242_m1`, `Ba137m`).
- Naming dialects: MCNP ZAID, Serpent, FLUKA, NIST, CINDER, ALARA, SZA.
- Particle registry and reaction-name registry (labels, MT mapping, hashes).
- Physical data access: AME2020 masses, natural abundances, half-lives.

## Capability crates

### `nucleide-material`

Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium
loading, and materials XML export. Depends on `nucleide-nuclei`.

### `nucleide-mcnp-io`

MCNP-family file I/O: `xsdir`, `meshtal`, SSW/SURFSRC, PTRAC, WWINP, MCTAL
readers; material extraction from input decks; mesh-to-geometry deck
generation. Depends on `nucleide-nuclei`.

### `nucleide-serpent-io`

Parsers for Serpent MATLAB-style output files (`_res.m`, `_dep.m`, `_det.m`).
No internal crate dependencies.

### `nucleide-fluka-io`

FLUKA interface: USRBIN tally reader and MATERIAL/COMPOUND card generation.
Depends on `nucleide-nuclei`.

### `nucleide-alara-io`

ALARA activation-code interop: input-deck, group-flux, material/element/WDR
library, activation-output, photon-source, and schedule-expansion glue. Depends
on `nucleide-nuclei` only among workspace crates; the solver stays inside ALARA.

### `nucleide-enrichment`

Multicomponent enrichment cascades and SWU analytics. Depends on `nucleide-nuclei` and
stays independent of `nucleide-material`.

### `nucleide-depletion`

CRAM matrix exponential (orders 16 and 48) and depletion-chain XML parsing.
Depends on `nucleide-linalg` and `nucleide-nuclei`.

### `nucleide-vr-tools`

MAGIC weight-window generation and mesh source sampling with alias tables.
Depends on `nucleide-mcnp-io` (`nucleide-nuclei` comes in transitively).

### `nucleide-cccc-io`

CCCC text-subset readers (ISOTXS multigroup libraries, RTFLUX/ATFLUX/RZFLUX
flux files) plus a minimal PARTISN deck writer with ISOTXS nuclide mapping.
Depends on `nucleide-nuclei` only among workspace crates; the solver stays out of
scope.

### `nucleide-fispact-io`

FISPACT-II inventory-output parser producing ALARA-compatible response
frames (`nucleide-alara-io` `ResponseFrame` rows). Depends on `nucleide-alara-io` and `nucleide-nuclei`;
activation solving stays inside FISPACT-II.

### `nucleide-origen-io`

Scoped ORIGEN 2.2 TAPE readers: `TAPE5` input echo, `TAPE6` output
inventories, `TAPE9`-style decay constants. Depends on `nucleide-nuclei` only among
workspace crates; burnup driving stays inside ORIGEN.

### `nucleide-r2s`

Scoped rigorous two-step (R2S) workflow builder: zone-to-flux linking from
ALARA decks, schedule expansion, and uniform-split photon-source assembly.
Depends on `nucleide-alara-io`, `nucleide-mcnp-io`, `nucleide-vr-tools`, `nucleide-material`, `nucleide-depletion`, and
`nucleide-nuclei`; transport and activation solving stay inside their respective
codes. The uniform split preserves only the total shutdown strength until
group-wise emission data is wired through per nuclide.

## Binding crates

### `nucleide-bindings`

PyO3 extension module that exposes workspace crates to Python as
`nucleide._internal`. It is the only crate allowed to depend on `pyo3` and the
only crate allowed to know about the Python API surface. It depends on most
capability crates.

### `nucleide-wasm`

`wasm-bindgen` crate that exposes a subset of the workspace to the browser for
the interactive tutorials on the docs site. It is not published to crates.io
(the release workflow's publish list excludes it); the site builds it with
`wasm-pack`.

## Dependency rules

- Workspace crates may depend on other workspace crates.
- Workspace crates must not depend on `nucleide-bindings`, `nucleide-wasm`, or `pyo3`.
- `nucleide-bindings` may depend on workspace crates.
- `nucleide-enrichment` must not depend on `nucleide-material`.

## Release order

When publishing to crates.io, publish in dependency order:

1. `nucleide-linalg`
2. `nucleide-nuclei`
3. `nucleide-material`
4. `nucleide-mcnp-io`
5. `nucleide-serpent-io`
6. `nucleide-fluka-io`
7. `nucleide-enrichment`
8. `nucleide-depletion`
9. `nucleide-vr-tools`
10. `nucleide-alara-io`
11. `nucleide-cccc-io`
12. `nucleide-fispact-io`
13. `nucleide-origen-io`
14. `nucleide-r2s`
15. `nucleide-bindings`
