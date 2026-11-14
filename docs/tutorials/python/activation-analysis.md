---
title: Activation Analysis
sidebar:
  order: 7
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

## Zone mappings, explicit volumes, and regrouping

The deck summary includes ordered `mat_loading` entries (`zone`, `mixture`)
and explicit `volumes` entries (`zone`, `volume_cm3`). Missing blocks give
empty lists. Join by zone name, since the blocks may use different orders:

```python
mapping_text = """geometry rectangular
mixture Steel
material external_steel 1 1
end
mat_loading
ZoneA Steel
Empty void
end
volumes
1 Empty
2 ZoneA
end
"""
mapping = alara_parse_deck(mapping_text)
volumes = {entry["zone"]: entry["volume_cm3"] for entry in mapping["volumes"]}
joined = [(entry["zone"], entry["mixture"], volumes[entry["zone"]])
          for entry in mapping["mat_loading"]]
assert joined == [("ZoneA", "Steel", 2.0), ("Empty", "void", 1.0)]
```

A mixture reference is not an evaluated material composition. Load and resolve
external compositions yourself. Parsing does not check cross-references;
`alara_validate_deck` performs that separate check. Volumes are cm³ and are
never inferred from geometry.

`regroup_flux` moves group-integrated flux between caller-supplied ascending
eV boundaries with exactly matching endpoints. Its required model assumes
constant differential flux per eV inside each source group:

```python
from nucleide.alara import regroup_flux

regrouped = regroup_flux([0, 2, 5], [4, 9], [0, 1, 3, 5], model="flat_energy")
assert regrouped["flux"] == [2.0, 5.0, 6.0]
assert regrouped["source_total"] == regrouped["target_total"] == 13.0
assert regrouped["split_source_groups"] == [0, 1]  # zero-based source indices
```

The report names the model and identifies every split group, including groups
with zero flux. Identity grids preserve the input values. Finite nonnegative
flux and strictly ascending finite nonnegative edges are required. There is
no clipping or extrapolation. Totals and regrouping contributions use compensated summation to retain
small contributions beside large ones. Values below the representable
floating-point range can still round to zero. Total flux is conserved up
to roundoff;
arbitrary reaction rates are not. Splitting 175 groups into 709 cannot recover
unmeasured fine structure, so activation comparisons still require matching
native grids. Photon centroids cannot supply bin boundaries.

## Split integrated sources over unequal voxel volumes

Multiply a parsed photon density by its activation-zone volume once to get
photons/s. The volume option then splits that integrated total at uniform
source density within each zone, preserving caller voxel order:

```python
from nucleide.r2s import tag_zone_strength

zone_totals = [30.0 * 4.0, 20.0 * 2.0]  # density × zone cm³ -> photons/s
voxel_volumes = [1.0, 3.0, 2.0]
tags = tag_zone_strength(zone_totals, [0, 0, 1], split=True,
                         voxel_volumes_cm3=voxel_volumes)
assert tags["source_strength"] == [30.0, 90.0, 40.0]
assert tags["total"] == 160.0
```

Outputs are integrated strengths, not densities. Zone sums and the returned
volume-path total use compensated summation; outputs too small to represent
can round to zero. Volumes must be finite and
positive; every positive-source zone must have voxels. Supplying volumes
requires `split=True`. Omitting volumes preserves the existing equal-count
split; omitting `split` copies the zone total to each voxel. Tag decay-time
metadata remains shutdown (`0.0`), even for strengths selected at another time.

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

## ARMI database snapshots without ARMI in the loop

`r2s_from_snapshot` builds the same workflow from a versionless snapshot dict
— no HDF5 dependency, no ARMI checkout, and no layout versioning on this
side. You run your own version-adapted dump script (your `h5py` plus your
ARMI checkout) and pass plain dicts; if a new ARMI release changes location
packing, only your script changes, never this schema. Zone ids are opaque
strings (1:1 from block names such as `B0120-005`); merging blocks into
homogenized zones is your script's job and is never inferred here.

Units the dict carries:

| Field | Unit |
| --- | --- |
| `volume_cm3` | cm³ |
| `composition` values | atoms/barn-cm |
| `zbottom_cm` / `ztop_cm` | cm (informational, kept for round-trip) |
| `temperature_C` | °C (informational, kept for round-trip) |
| `cooling_s` | seconds |
| flux files named by `flux_defs` | n/cm²/s per group (see note) |

ARMI `mgFlux` is volume-integrated (n·cm/s): divide by the homogenized zone
volume when writing each per-zone flux file, or fold the normalization into
`flux_defs[].scale`. Composition keys must be post-expansion nuclide names
(`U235`, `nPu239`, `AM242M`); elemental keys are rejected with "expand
first", bare `AM242` is rejected (pass `AM242M`/`AM242G` explicitly), and an
empty composition means a `void` zone that the workflow skips. Cooling times
and the irradiation schedule are caller-supplied: empty `cooling_s` fails
validation, and `schedule_text` may carry `schedule`/`pulsehistory` blocks
only (omit it for a 1-day-per-flux default).

Adapt the dump script per ARMI release: pin the ARMI commit, discover
`*Block*` groups dynamically instead of hardcoding concrete class names, and
join the `layout`/`type` + `name` + `indexInData` rows with your version's
location packers.

<!-- code-test: skip -->
```python
import h5py
import json

# USER-SIDE script: needs your h5py + your ARMI checkout, adapted to your
# database layout version. Never import this from Nucleide.
with h5py.File("snap016.h5", "r") as db:
    core = db["c01n01"]  # cycle/node group; names vary by model
    zones = []
    for key, grp in core.items():
        if "Block" not in key:
            continue  # discover concrete block groups dynamically
        loc = grp["layout"][...]  # row-address with type+name+indexInData
        vol = float(grp["volume"][...].sum())
        comp = {}
        for nuclide, ndens in zip(
            grp["componentNuclides"][...], grp["numberDensities"][...]
        ):
            comp.setdefault(nuclide.decode(), 0.0)
            comp[nuclide.decode()] += float(ndens) * vol
        zones.append(
            {
                "id": key,  # e.g. B0120-005; opaque downstream
                "volume_cm3": vol,
                "zbottom_cm": float(grp["zbottom"][...]),
                "ztop_cm": float(grp["ztop"][...]),
                "xs_type": str(grp["xsType"][...]),
                "composition": {k: v / vol for k, v in comp.items()},
            }
        )

snapshot = {
    "zones": zones,
    "flux_defs": [{"name": "snap_flux", "file": "data/flux", "scale": 1.0}],
    "cooling_s": [86400.0],
}
Path("snapshot.json").write_text(json.dumps(snapshot))
```

<!-- code-test: skip -->
```python
import json
from pathlib import Path

from nucleide.r2s import r2s_expand, r2s_from_snapshot, r2s_validate

bundle = r2s_from_snapshot(json.loads(Path("snapshot.json").read_text()))
r2s_validate(bundle["workflow"], bundle["deck"])
print(bundle["workflow"]["steps"], bundle["workflow"]["top_schedule"])
print(len(r2s_expand(bundle["deck"])), len(bundle["decks"]))
```

## Sample files

Sample reference files live under `fixtures/alara` (decks, flux files, output
listings), `fixtures/fispact` (inventory tables), and `fixtures/origen`
(TAPE5/6/9 samples); see the
[sample data files](../../reference/fixtures.mdx) index for the full list.
