---
title: Build Materials
sidebar:
  order: 4
---

Nucleide models materials as maps from canonical nuclide identifiers to masses,
with optional density and metadata. Materials live in the `material` crate and
are exposed through `nucleide.material.MaterialsCompendium`, `from_formula`, and
related helpers.

## Build from a chemical formula

```python
from nucleide.material import from_formula

comp = from_formula("UO2")
print(comp)
```

## Load the PNNL Materials Compendium

The compendium ships as a standalone JSON file; it is not bundled in the
Python wheel. Download it with `nucleide.data`, which pins the download to the
installed version's tag (in a source checkout, the same file is already under
`fixtures/data/`):

```python
from nucleide.data import fetch_compendium
from nucleide.material import MaterialsCompendium

path = fetch_compendium()  # add ref="main" or a commit SHA to override
lib = MaterialsCompendium.load(path)
print(len(lib), "materials")
print(lib.names()[:5])
```

The equivalent one-liner with curl:

```bash
curl -LO https://raw.githubusercontent.com/nukehub-dev/nucleide/main/fixtures/data/MaterialsCompendium.json
```

## Mix materials

On the Rust side, `Material::mix_by_mass` and `Material::mix_by_volume` combine
compositions by mass or volume fractions. Mixing is not yet exposed through the
Python facade.

## Export materials XML

A composition dictionary can be serialized to an OpenMC-style `<material>` XML
fragment. `to_xml` takes a name, density, and optional density units:

```python
from nucleide.material import from_formula, to_xml

comp = from_formula("UO2")
xml = to_xml(comp, "fuel", 10.0, "g/cm3")
print(xml)
```

## See also

- Crate docs in [`crates/material/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/material/src/lib.rs).
- Fixture license in `fixtures/data/MaterialsCompendium.LICENSE`.
