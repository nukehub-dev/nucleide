---
title: Activation Analysis
sidebar:
  order: 5
---

Nucleide reads activation-code inputs and outputs without running the codes
themselves. This tutorial shows the Python API; equivalent Rust types live in
the `alara-io`, `fispact-io`, `origen-io`, and `r2s` crates.

## ALARA decks, flux, and output

Parse an input deck, a group-flux file, and the activation-output listing into
a shared response frame, then expand the irradiation schedule.

```python
from pathlib import Path

from nucleide.alara import (
    alara_expand_schedule,
    alara_parse_deck,
    alara_parse_flux,
    alara_parse_output,
)

deck_text = Path("fixtures/alara/decks/sample2").read_text()
deck = alara_parse_deck(deck_text)
print(deck["geometry"], [m["name"] for m in deck["mixtures"]])

flux = alara_parse_flux(Path("fixtures/alara/flux/fluxin2").read_text(), "fluxin2")
print(flux["num_intervals"], flux["groups_per_interval"], flux["total"])

rows = alara_parse_output(
    Path("fixtures/alara/output/sample2.out").read_text(), "sample2"
)
print(len(rows), sorted({r["variable"] for r in rows}))

steps = alara_expand_schedule(deck_text)
print([(s["flux"], s["duration_s"], s["is_cooling"]) for s in steps])
```

## FISPACT output in the same frame shape

FISPACT-II inventory tables parse into the same row shape as ALARA output, so
both codes share one analysis frame. The inventory table carries no half-life
column, so every row reports the `half_life_s` sentinel value of `-1.0`.

```python
from pathlib import Path

from nucleide.fispact import fispact_parse_output

rows = fispact_parse_output(
    Path("fixtures/fispact/inventory.fis").read_text(), "synth"
)
print(len(rows), sorted({r["time_s"] for r in rows}))
print(sorted({r["variable"] for r in rows}))
print({r["half_life_s"] for r in rows})
```

## ORIGEN TAPE5, TAPE6, and TAPE9

Read the TAPE5 input echo (titles, irradiation steps, materials), the TAPE6
output inventories (per-nuclide grams and activity plus the total activity),
and the TAPE9-style decay constants.

```python
from pathlib import Path

from nucleide.origen import (
    origen_parse_tape5,
    origen_parse_tape6,
    origen_parse_tape9,
)

tape5 = origen_parse_tape5(Path("fixtures/origen/tape5_sample").read_text())
print(tape5["titles"], tape5["irradiation_steps"])

tape6 = origen_parse_tape6(Path("fixtures/origen/tape6_sample").read_text())
print(len(tape6["records"]), tape6["total_activity"])

tape9 = origen_parse_tape9(Path("fixtures/origen/tape9_sample").read_text())
print([(e["nuclide"], e["decay_const"]) for e in tape9[:2]])
```

## R2S shutdown-dose-rate workflow

Build a rigorous two-step workflow from an ALARA deck, validate it against the
same deck, expand its schedule, and assemble a decay photon source for one
zone. Assembly sums shutdown `Specific Activity` over the zone and splits the
total uniformly over the requested groups, so it preserves only the total
shutdown strength — not the nuclide- and energy-dependent lines of a
`.photonSrc` spectrum.

```python
from pathlib import Path

from nucleide.r2s import r2s_assemble, r2s_expand, r2s_from_deck, r2s_validate

deck_text = Path("fixtures/alara/decks/sample2").read_text()
workflow = r2s_from_deck(deck_text)
r2s_validate(workflow, deck_text)
print(workflow["steps"], workflow["top_schedule"])

print(len(r2s_expand(deck_text)))

source = r2s_assemble(
    Path("fixtures/alara/output/sample2.out").read_text(),
    "sample2",
    "inner_zone",
    4,
)
print(source["zone"], source["total"], source["groups"])
```

## Fixtures

Reference files live under `fixtures/alara` (decks, flux files, output
listings), `fixtures/fispact` (inventory tables), and `fixtures/origen`
(TAPE5/6/9 samples).
