---
title: Emit Code Cards
sidebar:
  order: 7
---

Nucleide renders one material through five legacy transport-code dialects —
MCNP, Serpent, FLUKA, ALARA, PARTISN — and reports how much mass survives
each translation. This tutorial shows the Python API; equivalent Rust types
live in the `emit` crate.

## Cards for all five codes

```python
from nucleide import emit

comp = {"U235": 5.0, "U238": 95.0}
cards = emit.emit_cards(comp, "umetal", density=19.1)
print(sorted(cards))
print(cards["MCNP"])
```

The MCNP card uses mass fractions with a configurable library suffix, the
Serpent card uses `{zaid}.{lib}` ids, and long MCNP cards wrap within the
128-column limit:

```python
print(cards["Serpent"])
print(cards["ALARA"])
```

## Mass-drift report

`emit_drift_table` returns one row per dialect with the represented mass,
the relative drift, dropped nuclides with reasons, and whether the text was
re-parse verified (MCNP and ALARA only):

```python
rows = emit.emit_drift_table(comp, "umetal", density=19.1)
for row in rows:
    print(row["code"], row["rel_drift"], row["dropped"], row["reparsed"])
```

Nuclides a dialect cannot represent — O16 has no FLUKA isotope-table entry —
show up as reported drift instead of silent loss:

```python
rows = emit.emit_drift_table({"U235": 1.0, "O16": 2.0}, "uox", density=10.0)
fluka = next(r for r in rows if r["code"] == "FLUKA")
print(fluka["mass_out"], fluka["dropped"])
```
