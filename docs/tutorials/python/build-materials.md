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

`mix_by_mass` combines `(composition, weight)` streams by relative mass;
`mix_by_volume` combines `(composition, volume, density)` triples through
each stream's density:

```python
from nucleide.material import mix_by_mass, mix_by_volume

mass = mix_by_mass([({"U235": 19.0}, 1.0), ({"U238": 1.0}, 1.0)])
vol = mix_by_volume([({"U235": 19.0}, 1.0, 19.1), ({"U238": 1.0}, 1.0, 19.1)])
print(mass, vol)
```

## Specific activity and multi-material documents

`specific_activity` returns the whole-material Bq/g; `materials_doc_to_xml`
bundles named materials (each with its density) into one `<materials>`
document, and `expand_elements` / `collapse_elements` convert between
natural-element placeholders (bare symbols such as `"U"`) and isotopic
breakdowns:

```python
from nucleide.material import (
    collapse_elements,
    expand_elements,
    materials_doc_to_xml,
    specific_activity,
)

print(specific_activity({"U235": 1.0}))
xml = materials_doc_to_xml([("fuel", {"U235": 19.0}, 10.0)], cross_sections="xs.xml")
print(expand_elements({"U": 20.0}))
print(collapse_elements({"U235": 19.0, "U238": 1.0}))  # {"U": 20.0}
```

## Export materials XML

A composition dictionary can be serialized to an OpenMC-style `<material>` XML
fragment. `to_xml` takes a name, density, and optional density units:

```python
from nucleide.material import from_formula, to_xml

comp = from_formula("UO2")
xml = to_xml(comp, "fuel", 10.0, "g/cm3")
print(xml)
```

## Dose per gram

`dose_per_g` maps a composition in grams through the vendored dose-factor
tables (`nucleide.nuclei.dose_factor` over the HNF-5636 / PyNE
`dbgen/dosefactors` basis: 93 folded nuclides × 4 pathways × 3 sources;
`+D` daughters fold into the parent, GENII/DOE air cells are `-1`
sentinels and error). The pathway is one of `air`, `soil`, `ingest`, or
`inhale`; the source defaults to `EPA` (`DOE` and `GENII` also resolve).
Units follow the table — air in mrem/h per g per m³, soil in mrem/h per
g per m², ingest/inhale in mrem per g — following PyNE's
`Material::dose_per_g` equations with per-gram map semantics. Results
are screening-level only, never for safety decisions:

```python
from nucleide.material import dose_per_g

print(dose_per_g({"Co60": 1.0}, "ingest"))  # mrem per g
print(dose_per_g({"Co60": 1.0}, "air"))  # mrem/h per g per m^3
```

Nuclides lacking mass data are errors, as are radioactive nuclides
lacking dose data (including the `-1` sentinels). Stable nuclides — known
mass but no decay constant — contribute exactly 0.

## See also

- Crate docs in [`crates/material/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/material/src/lib.rs).
- Fixture license in `fixtures/data/MaterialsCompendium.LICENSE`.
