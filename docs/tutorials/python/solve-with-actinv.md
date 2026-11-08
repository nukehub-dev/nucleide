---
title: Solve activation with ACTINV
sidebar:
  order: 23
---

Nucleide builds R2S workflows, reads activation-code inputs and outputs,
and assembles shutdown photon sources — but it does not solve activation
itself. [ACTINV](https://github.com/AvilaLabs/ACTINV) is an optional
external solver that fills that seat: you drive it from its command
line, with problem files written in the `actinv-spec-1` JSON format.
Nucleide never requires it — everything in Nucleide works without
it. The ACTINV example files live in the ACTINV repository (they
maintain them against their releases); this page shows where Nucleide
hands off.

## Check whether ACTINV is installed

```python
import shutil

print(shutil.which("actinv"))
```

```text
None
```

`None` simply means the external solver is not on your `PATH` yet —
skip to the install step below. Your Nucleide workflow is unaffected.

## Install the solver and its data

ACTINV software and its nuclear data are separate downloads. The
software installs with pip; the data (the activation library, decay
data, and their index files) is fetched by the CLI into a versioned
folder and checked by SHA-256 before use:

```bash
python -m pip install actinv==1.4.0
actinv data fetch
actinv data verify
```

The default data bundle is `tendl-2025-neutron` from data catalog
`1.1.0` (about 140 MiB to download). The default neutron library
artifact `tendl-2025-neutron-709g` is pinned to SHA-256
`ec4c72bf598dc8ad3d533d9cfafdcf493e2d1f949a3e4db6251495659b68cc44`;
`actinv data manifest` prints the catalog with every artifact hash, and
`actinv data verify` re-checks installed files against it. Record the
hashes you actually used with your calculation — a different bundle
changes the answer.

## Map your Nucleide inputs

The smallest mapping is one material, one spectrum, and one schedule,
written as an `actinv-spec-1` document. The three Nucleide inputs come
from the same surfaces used in
[Activation analysis](activation-analysis.md): deck mixtures,
group fluxes, and the expanded irradiation schedule.

```python
from pathlib import Path

from nucleide.alara import alara_parse_deck, alara_parse_flux

deck_text = Path("fixtures/alara/decks/sample2").read_text()
deck = alara_parse_deck(deck_text)
print(deck["geometry"], [m["name"] for m in deck["mixtures"]])

flux = alara_parse_flux(Path("fixtures/alara/flux/fluxin2").read_text(), "fluxin2")
print(flux["groups_per_interval"], flux["num_intervals"], flux["total"])
```

```text
rectangular ['inner_mix', 'outer_mix']
175 3 4177383924000000.0
```

Expanding the deck's hierarchical schedule gives flat steps — each a
duration in seconds plus the flux name it runs at (an empty flux name
means cooling):

```python
from nucleide.alara import alara_expand_schedule

steps = alara_expand_schedule(deck_text)
print(steps)
```

```text
[{'duration_s': 31557600.0, 'flux': 'flux_1', 'is_cooling': False}]
```

Each input maps onto the spec as follows (all six points are stated,
not assumed):

| # | Nucleide side | `actinv-spec-1` side |
| --- | --- | --- |
| 1 | Group fluxes in ALARA group order (175 groups per interval here) | `spectrum.structure: "fispact-709"` with exactly 709 group-integrated fluxes, or a `custom` structure whose boundaries match the library. An ALARA-order file cannot be pasted as-is — supply a 709-group spectrum (for example a transport spectrum) or rebin onto the library boundaries while conserving the total flux. |
| 2 | Deck `flux <name> <file> <scale>` scaling on top of the file values | Physical group-integrated fluxes in n/cm²/s. Fold the deck scale into your numbers first, then write `spectrum.total` for the sum you intend; ACTINV rescales the vector to that total while preserving its shape. |
| 3 | ALARA output spellings such as `mn-56` | Canonical `Mn56` (isomers such as `Ba137m1`; a bare `m` means `m1`). Composition keys accept an element (`Fe`) or an explicit nuclide (`Fe56`), but never both together. Nucleide's own reader accepts both spellings, so the conversion happens at this boundary. A wrong key raises a clear error. |
| 4 | `FlatStep {duration_s, flux}` with delay steps flux-off | `schedule` entries `{"dt": ..., "flux": ...}` — irradiation at `flux: 1.0`, cooling (any flux-off step) at `flux: 0.0`. Seconds are accepted, as are minutes, hours, days, and years. |
| 5 | — | Data provisioning: run `actinv data fetch`, keep the library hash with the calculation record, and check hashes with `actinv validate ... --hashes` before running. |
| 6 | The photon-source handoff below | Return path carries `TOTAL` rows only (no per-nuclide spectra), on a centroid grid, one volume element per zone in the demo. Treat it as a source input, not as nuclide-resolved data. |

Building the spec for the iron case the ACTINV example ships —
one gram of pure iron, the measured 709-group spectrum stored
highest-energy-first (`descending: true`) at its recorded total, and
five minutes of irradiation followed by cooling — looks like this:

```python
schedule = [
    {"dt": f"{step['duration_s']} s", "flux": 0.0 if step["is_cooling"] else 1.0}
    for step in steps
]
spec = {
    "spec": "actinv-spec-1",
    "title": "FNS Fe, 5-minute irradiation, 1996 campaign",
    "projectile": "neutron",
    "library": {"path": "catalog:tendl-2025-neutron-709g"},
    "decay": {
        "primary": "catalog:endfb-viii-0-decay",
        "fallback": "catalog:jeff-3-3-decay",
    },
    "material": {"mass_g": 1.0, "basis": "wt_percent", "composition": {"Fe": 100.0}},
    "spectrum": {
        "structure": "fispact-709",
        "descending": True,
        "total": 1.116e10,
        "flux_per_group": "<709 group-integrated fluxes>",
    },
    "schedule": [{"dt": "300.0 s", "flux": 1.0}, {"dt": "3600.0 s", "flux": 0.0}],
}
print(spec["material"]["composition"], spec["spectrum"]["structure"])
print(schedule)
```

```text
{'Fe': 100.0} fispact-709
[{'dt': '31557600.0 s', 'flux': 1.0}]
```

In practice you start from `actinv new problem.json` (a complete,
editable iron example with the full spectrum — you do not type a
709-group vector) and edit the title, composition, and schedule toward
your case. Reuse the shipped spectrum only when it represents your
irradiation; changing the material does not make that spectrum
appropriate for a different facility.

## Validate and solve

```bash
actinv new problem.json
actinv validate problem.json --hashes
actinv run problem.json result.json
```

Or shell out from Python:

<!-- code-test: skip -->
```python
import subprocess

subprocess.run(["actinv", "validate", "problem.json", "--hashes"], check=True)
subprocess.run(["actinv", "run", "problem.json", "result.json"], check=True)
```

Read the final step of the result:

<!-- code-test: skip -->
```python
import json
from pathlib import Path

result = json.loads(Path("result.json").read_text(encoding="utf-8"))
last = result["steps"][-1]
print("Time (s):", last["t_s"])
print("Activity (Bq/g):", sum(last["activity_Bq_per_g"].values()))
print("Decay heat (W/g):", last["heat_W_per_g"]["total"])
```

The result's `ledger` lists data and calculation notes — which data
files were used, and which reactions or products were not covered.
Keep it with the result: missing coverage is reported there, never
silently filled in.

## Read the photon source back with Nucleide

For an R2S calculation, export the decay photon source and read it
with Nucleide's own ALARA photon-source readers. A multi-cell mesh
calculation exports through `actinv export-r2s` and then one ALARA
`.photonSrc` file per cell plus a provenance index:

```bash
actinv export-r2s mesh-result.ndjson STEP r2s-source.ndjson
actinv export-source alara r2s-source.ndjson photon-src --shutdown-t-s T
```

<!-- code-test: skip -->
```python
from pathlib import Path

from nucleide.alara import alara_photon_total_strength
from nucleide.r2s import photon_group_sums, tag_zone_strength

index_text = Path("photon-src/actinv-alara-index.json").read_text()
# For each cell entry in the index, load its .photonSrc file with
# Nucleide's readers and tag the voxel through r2s.tag_zone_strength
# (one voxel per zone, since mesh cells are already independent
# volume elements).
```

A single-cell calculation instead exports a point source with
`actinv export-openmc result.json STEP source.py` or
`actinv export-mcnp result.json STEP source.sdef`. The ACTINV
repository's output-side demo (`contrib/nucleide_r2s`) runs exactly
this readback and prints a conservation table: each cell's reconstructed
total against the declared photons per second.

## Caveats

- **Research-grade.** ACTINV is calculation software for research, not
  a licensing or safety basis. Record the exact software version,
  problem, data versions, and hashes; apply your own applicability,
  uncertainty, and independent-review process.
- **Known data defects.** The default TENDL-2025 neutron corpus carries
  known upstream defects, including affected product-state data for 17
  IRDFF-II benchmark targets. Check target coverage for your
  reactions before trusting an inventory, and read the result ledger.
- **TOTAL-only photons.** The return path carries total photon rows,
  not per-nuclide spectra, on a centroid grid rather than energy-bin
  bounds. Supply the matching bounds yourself if you run the transport
  leg downstream.
- **Data is separate from software.** Installing `actinv` does not
  install data, and updating the solver does not update your data.
  Verify hashes; for offline machines, transfer the files listed by
  `actinv data manifest` and check them with `actinv data verify`.

## See also

- [Activation analysis](activation-analysis.md) for the Nucleide-side
  deck, flux, schedule, and R2S workflow surfaces.
- ACTINV's [handbook and specification](https://github.com/AvilaLabs/ACTINV/blob/master/docs/guide/specification.md),
  its [iron example](https://github.com/AvilaLabs/ACTINV/tree/master/examples/fns_iron),
  and its [Nucleide R2S demo](https://github.com/AvilaLabs/ACTINV/tree/master/contrib/nucleide_r2s).
- [Python API reference](../../reference/python-api.mdx#nucleidealara) for
  the full `nucleide.alara` and `nucleide.r2s` signatures.
