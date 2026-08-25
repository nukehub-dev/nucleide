# Build Materials

Nucleide models materials as maps from canonical nuclide identifiers to masses,
with optional density and metadata. Materials live in the `material` crate and
are exposed through `nucleide.MaterialsCompendium`, `from_formula`, and related
helpers.

## Build from a chemical formula

```python
from nucleide import from_formula

comp = from_formula("UO2")
print(comp)
```

## Load the PNNL Materials Compendium

```python
from nucleide import MaterialsCompendium

lib = MaterialsCompendium.load("fixtures/data/MaterialsCompendium.txt")
print(len(lib), "materials")
print(lib.names()[:5])
```

## Mix materials

On the Rust side, `Material::mix` combines compositions by mass fractions. The
Python facade exposes equivalent helpers through the `_internal` module.

## Export materials XML

The `material` crate can serialize a library to the materials XML format used by
several transport codes. From Python:

```python
# see nucleide._internal for the current XML export helpers
```

## See also

- Crate docs in `crates/material/src/lib.rs`.
- Fixture license in `fixtures/data/MaterialsCompendium.LICENSE`.
