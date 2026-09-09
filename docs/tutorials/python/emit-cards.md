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

## ARMI blueprint keys

`emit_armi_cards` / `emit_armi_drift_table` accept ARMI-side mass-fraction
dicts — database names (`nU235`), ZAIDs (`92235`), AAAZZZS ids, and
unambiguous MC2-3 labels (`U-2355`, `U235_7`) — and emit the same five cards
as GNDS names:

```python
armi_cards = emit.emit_armi_cards({"nU235": 5.0, "nU238": 95.0}, "umetal", density=19.1)
assert armi_cards == emit.emit_cards({"U235": 5.0, "U238": 95.0}, "umetal", density=19.1)
```

The bridge is one-way (dict in, never ARMI YAML or `.h5`) and keeps
ARMI-side resolution with the caller:

- Pass post-expansion nuclide fractions (`material.massFrac`). Elemental
  keys (`ZR`, `FE`, `O`, …) are rejected with an "expand first" error.
- Values are masses in grams; number fractions/densities must be converted
  first. Enrichment shorthands, `balance`, and temperatures are resolved by
  the caller, which passes the hot density in g/cm³ (`None` keeps the usual
  `MissingDensity` behavior for Serpent/FLUKA/PARTISN).
- Bare `AM242` is rejected as ambiguous — pass `AM242M` for Am-242m or
  `AM242G` for ground explicitly.
