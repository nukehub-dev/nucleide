---
title: Deterministic I/O
sidebar:
  order: 8
---

Nucleide reads deterministic-code data files and writes minimal PARTISN
inputs without running a transport solve. This tutorial shows the Python API;
equivalent Rust types live in the `cccc-io` crate.

## ISOTXS libraries

Parse a multigroup library, then find one nuclide entry by label.

```python
from pathlib import Path

from nucleide.cccc import isotxs_parse

lib = isotxs_parse(Path("fixtures/cccc/isotxs_sample").read_text())
u235 = next(n for n in lib["nuclides"] if n["label"] == "U235")
print(u235["zaid"], u235["groups"], u235["total_xs"])
```

## Flux files

Parse an RTFLUX file (ATFLUX and RZFLUX use the same reader with
`kind="atflux"` or `kind="rzflux"`), then read out one spatial point's group
values.

```python
from pathlib import Path

from nucleide.cccc import rtflux_parse

flux = rtflux_parse(Path("fixtures/cccc/rtflux_sample").read_text())
print(flux["kind"], flux["npoints"], flux["groups"], flux["total"])

point = flux["values"][: flux["per_point"]]
print(point)
```

## PARTISN decks

Render a minimal PARTISN deck and validate its ISOTXS labels against library
text. A label with no library entry fails validation.

```python
from pathlib import Path

from nucleide.cccc import partisn_render, partisn_validate

deck = {
    "title": "synthetic slab",
    "dim": 1,
    "zones": [
        {"id": 1, "material": "fuel", "isotxs_labels": ["U235"], "density": 10.0},
        {
            "id": 2,
            "material": "blanket",
            "isotxs_labels": ["PU239"],
            "density": 5.0,
        },
    ],
    "source": "isotropic",
}
print(partisn_render(deck))

isotxs_text = Path("fixtures/cccc/isotxs_sample").read_text()
partisn_validate(deck, isotxs_text)

deck["zones"][0]["isotxs_labels"].append("U238")
try:
    partisn_validate(deck, isotxs_text)
except ValueError as exc:
    print("dangling label:", exc)
```

## Fixtures

Reference files live under `fixtures/cccc`. Both samples are synthetic
fixtures authored for parser tests, not measured data.
