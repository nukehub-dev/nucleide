# Crate Responsibilities

Each Nucleide crate has a single responsibility and a well-defined place in the
workspace dependency graph.

## Foundation crates

### `linalg`

Isolation facade over the linear-algebra backend. Today it pulls in `faer`,
`num-complex`, and `roxmltree` so numeric dependencies do not leak into other
crates. Other workspace crates depend on `linalg`, not on the backend directly.

### `nuclei`

Canonical nuclide identification. Owns:

- `NuclideId`: compact `u32` representation compatible with legacy integer ids.
- Element symbol/number tables.
- Name parsing (`U235`, `Am242_m1`, `Ba137m`).
- Naming dialects: MCNP ZAID, Serpent, FLUKA, NIST, CINDER, ALARA, SZA.
- Particle registry and reaction-name registry (labels, MT mapping, hashes).
- Physical data access: AME2020 masses, natural abundances, half-lives.

## Capability crates

### `material`

Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium
loading, and materials XML export. Depends on `nuclei`.

### `mcnp-io`

MCNP-family file I/O: `xsdir`, `meshtal`, SSW/SURFSRC, PTRAC, WWINP, MCTAL
readers; material extraction from input decks; mesh-to-geometry deck
generation. Depends on `nuclei`.

### `serpent-io`

Parsers for Serpent MATLAB-style output files (`_res.m`, `_dep.m`, `_det.m`).
No internal crate dependencies.

### `fluka-io`

FLUKA interface: USRBIN tally reader and MATERIAL/COMPOUND card generation.
Depends on `nuclei`.

### `enrichment`

Multicomponent enrichment cascades and SWU analytics. Depends on `nuclei` and
stays independent of `material`.

### `depletion`

CRAM matrix exponential (orders 16 and 48) and depletion-chain XML parsing.
Depends on `linalg` and `nuclei`.

### `vr-tools`

MAGIC weight-window generation and mesh source sampling with alias tables.
Depends on `mcnp-io` and `nuclei`.

## Binding crate

### `nucleide-bindings`

PyO3 extension module that exposes workspace crates to Python as
`nucleide._internal`. It is the only crate allowed to depend on `pyo3` and the
only crate allowed to know about the Python API surface. It depends on most
capability crates.

## Dependency rules

- Workspace crates may depend on other workspace crates.
- Workspace crates must not depend on `nucleide-bindings` or on `pyo3`.
- `nucleide-bindings` may depend on workspace crates.
- `enrichment` must not depend on `material`.

## Release order

When publishing to crates.io, publish in dependency order:

1. `linalg`
2. `nuclei`
3. `material`
4. `mcnp-io`
5. `serpent-io`
6. `fluka-io`
7. `enrichment`
8. `depletion`
9. `vr-tools`
10. `nucleide-bindings`
