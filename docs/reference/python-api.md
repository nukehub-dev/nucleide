---
title: Python API
sidebar:
  order: 2
---

Nucleide exposes a typed pure-Python package at `python/nucleide/` that re-exports
a compiled PyO3 extension built from `bindings/python/`. The extension module is
named `nucleide._internal`; users import from the `nucleide` domain submodules,
which mirror the Rust workspace crates.

## Module layout

```python
import nucleide as nuc
```

The top level carries only `version()` and `nucleide.__version__`. Everything
else lives in a domain submodule, each re-exporting its symbols from
`nucleide._internal`:

| Submodule | Backing crate | Contents |
| --- | --- | --- |
| `nucleide.nuclei` | `nuclei` | Nuclide/particle identifiers, nuclear data, reaction names |
| `nucleide.material` | `material` | Compositions, compendium, activity, XML export |
| `nucleide.mcnp` | `mcnp-io` | MCNP file readers and writers |
| `nucleide.serpent` | `serpent-io` | Serpent output readers |
| `nucleide.fluka` | `fluka-io` | FLUKA USRBIN readers |
| `nucleide.vr` | `vr-tools` | MAGIC weight windows, source sampling |
| `nucleide.enrichment` | `enrichment` | Enrichment cascades |
| `nucleide.depletion` | `depletion` | Depletion chains and CRAM solves |

## `nucleide.nuclei`

- `Nuclide(name)` — canonical nuclide identifier with cross-code name conversions.
- `Particle(spec)` — particle species with MCNP/FLUKA/Geant4 translations.
- `atomic_mass(key)` — by nucid integer or name.
- `natural_abundance(key)`
- `half_life(key)`
- `decay_constant(key)`
- `q_value_capture(key)`
- `q_value_alpha(key)`
- `rxname_id(name)` → reaction id
- `rxname_name(id)` → name
- `rxname_mt(id)` → ENDF MT number

## `nucleide.material`

- `from_formula(formula)` → composition `dict`
- `MaterialsCompendium.load(path)` → `MaterialsCompendium`
- `activity(comp)` → activity `dict`
- `to_xml(comp, name, density, units="g/cm3")` → materials XML string

## `nucleide.mcnp`

- `read_xsdir(path)` → `Xsdir`
- `read_meshtal(path)` → `Meshtal`
- `read_wwinp(path)` → `Wwinp`
- `read_mctal(path)` → `Mctal`
- `read_ssw(path)` → `SurfSrc`
- `read_ptrac(path)` → `PtracFile`
- `write_ssw(ssw, path, tracks=None)`
- `read_inp(path)` → list of material `dict`s
- `mesh_to_geom(...)` → MCNP geometry deck string

## `nucleide.serpent`

- `read_serpent(path, kind)` → structured `dict`

## `nucleide.fluka`

- `read_usrbin(path)` → list of `UsrbinTally`

## `nucleide.vr`

- `magic(tally, per_group=False, tolerance=0.5)` → `MagicOutput`
- `AliasTable(pdf)` — Walker alias table for discrete sampling
- `MeshSourceSampler(tally, mode, user_pdf=None)` — mesh source sampler

## `nucleide.enrichment`

- `Cascade.default_uranium()` → `Cascade`
- `c.solve(tolerance=None, max_iterations=None)`

## `nucleide.depletion`

- `read_chain(path)` → `Chain`
- `build_depletion_system(chain, rates)` → `DepletionSystem`
- `deplete(chain, n0, dt, rates=None, order=48)` → result `dict`

## Version

- `version()` → workspace version string
- `nucleide.__version__` — same value

For full type signatures, see `python/nucleide/_internal.pyi`.
