---
title: Reaction Names
sidebar:
  order: 15
---

Nucleide maps reaction names, ENDF MT numbers, and numeric ids through one
registry (`crates/nuclei/src/rxname.rs`, a port of the PyNE `rxname` map).
The Python entry points live in `nucleide.nuclei`.

## Names, ids, and MT numbers

```python
from nucleide.nuclei import rxname_id, rxname_mt, rxname_name

fid = rxname_id("fission")
print(fid, rxname_name(fid), rxname_mt(fid))  # MT 18
print(rxname_id("18") == fid)  # numeric MT strings resolve too
print(rxname_id("absorption") == rxname_id("27"))
```

## Labels, docs, and registry rows

```python
from nucleide.nuclei import rxname_doc, rxname_label, rxname_reaction

print(rxname_label(fid))  # "(z,fission)"
print(rxname_doc(fid))  # long description
print(rxname_reaction(fid))  # {"id", "name", "mt", "label", "doc"}
print(rxname_reaction(0xFFFFFFFF))  # None for unknown ids
```

Unknown ids return `""` for the label/doc and `0` for the MT number, matching
the upstream map default.

## Reaction graph: channels, daughters, parents

Channels connect nuclides under a projectile flag (`"n"`, `"p"`, `"d"`,
`"t"`, `"He3"`, `"a"`, `"gamma"`, `"decay"`; neutron-induced is the default):

```python
from nucleide.nuclei import rxname_child, rxname_id_from_nucdelta, rxname_parent

u235 = 922350000
u236 = 922360000
absorption = rxname_id("absorption")
print(rxname_id_from_nucdelta(u235, u236, "n") == absorption)
print(rxname_child("U235", "absorption"))  # "U236"
print(rxname_parent("U236", "absorption"))  # "U235"
```

`rx` accepts a canonical name or a numeric id.

## Particle helpers

Validity and classification helpers cover both particle aliases and nuclides:

```python
from nucleide.nuclei import (
    particle_is_heavy_ion,
    particle_is_hydrogen,
    particle_is_valid,
    particle_is_valid_pdc,
)

print(particle_is_valid("n"))  # True (alias, PDC number, or nuclide)
print(particle_is_valid_pdc(2112))  # True (neutron)
print(particle_is_hydrogen("Proton"))  # True
print(particle_is_heavy_ion("U235"))  # True
```

## See also

- `tests/test_facade_bundle.py` for runnable examples.
