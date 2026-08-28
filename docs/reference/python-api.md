---
title: Python API
sidebar:
  order: 2
---

Nucleide exposes a typed pure-Python package at `python/nucleide/` that re-exports
a compiled PyO3 extension built from `bindings/python/`. The extension module is
named `nucleide._internal`; users import from `nucleide` directly.

## Module layout

```python
import nucleide as nuc
```

All public symbols are re-exported from `nucleide._internal` in
`python/nucleide/__init__.py`.

## Core identifiers

- `Nuclide(name)` — canonical nuclide identifier with cross-code name conversions.
- `Particle(spec)` — particle species with MCNP/FLUKA/Geant4 translations.

## MCNP I/O

- `read_xsdir(path)` → `Xsdir`
- `read_meshtal(path)` → `Meshtal`
- `read_wwinp(path)` → `Wwinp`
- `read_mctal(path)` → `Mctal`
- `read_ssw(path)` → `SurfSrc`
- `read_ptrac(path)` → `PtracFile`
- `write_ssw(ssw, path, tracks=None)`
- `mesh_to_geom(...)` → MCNP geometry deck string

## Serpent and FLUKA

- `read_serpent(path, kind)` → structured `dict`
- `read_usrbin(path)` → list of `UsrbinTally`

## Materials

- `from_formula(formula)` → composition `dict`
- `MaterialsCompendium.load(path)` → `MaterialsCompendium`
- `activity(comp)` → activity `dict`

## Depletion

- `read_chain(path)` → `Chain`
- `deplete(chain, n0, dt, rates=None, order=48)` → result `dict`

## Enrichment

- `Cascade.default_uranium()` → `Cascade`
- `c.solve(tolerance=None, max_iterations=None)`

## Variance reduction

- `magic(tally, per_group=False, tolerance=0.5)` → `MagicOutput`

## Data lookups

- `atomic_mass(key)` — by nucid integer or name.
- `natural_abundance(key)`
- `half_life(key)`
- `decay_constant(key)`
- `q_value_capture(key)`
- `q_value_alpha(key)`

## Reaction names

- `rxname_id(name)` → reaction id
- `rxname_name(id)` → name
- `rxname_mt(id)` → ENDF MT number

## Version

- `version()` → workspace version string
- `nucleide.__version__` — same value

For full type signatures, see `python/nucleide/_internal.pyi`.
