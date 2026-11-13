---
title: Solve activation with ACTINV
sidebar:
  order: 23
---

Nucleide reads activation inputs and outputs and builds rigorous two-step
(R2S) workflows. [ACTINV](https://github.com/AvilaLabs/ACTINV) is an
optional external activation solver. This page explains the exchange
boundaries, records a fresh activation comparison, and demonstrates source readback; ACTINV
and evaluated nuclear data are not Nucleide dependencies.

## Install software separately from nuclear data

Use the [ACTINV 1.7.2 release](https://github.com/AvilaLabs/ACTINV/releases/tag/v1.7.2)
for the commands below. Its released Linux binary was tested with
`mesh`, `run`, `export-r2s`, and `export-source alara` using verified nuclear data. Download the
`actinv-1.7.2-linux-x86_64.tar.gz` asset and `SHA256SUMS`, verify the archive,
and extract it into a separate directory:

```bash
sha256sum --ignore-missing --check SHA256SUMS
tar -xzf actinv-1.7.2-linux-x86_64.tar.gz
export PATH="$PWD/target/x86_64-unknown-linux-gnu/release:$PATH"
actinv --version
```

The expected archive SHA-256 is
`ad39446999ef79575b77ac6daafa2d720304073b0e0700d1a962b7b1de393a1d`;
the extracted binary's SHA-256 is
`7d52a7b234f3125cc48ce222b32a1c22235acb756d89adbff386260387cf28d1`.
The release tag resolves to source commit
`d3dc4e463346c1a6d3dc426432048022932e14c5`.

The older published **1.4.0** release predates `export-source alara`.
Reproducing the older upstream example requires its specified source
commit `0353d5e3288a90cfefaabe7d2263b22fbaeeecad` and a source build
(`cargo build --release -p actinv-cli` in that checkout), rather than
installing the published 1.4.0 package and assuming the exporter exists.
The end-to-end adapter is recorded in
[upstream PR 6](https://github.com/AvilaLabs/ACTINV/pull/6), revision
[`6f8f893fdfc1dbe76257705d54f77b1ee68c633b`](https://github.com/AvilaLabs/ACTINV/tree/6f8f893fdfc1dbe76257705d54f77b1ee68c633b/contrib/nucleide_r2s).
That revision was unmerged when inspected on 2026-11-12; its adapter and
reported calculation results remain maintained by ACTINV.

Nuclear data is a second download:

```bash
actinv data manifest
actinv data fetch
actinv data verify
```

Release 1.7.2 defaults to catalog **1.2.0**, distinct from the older
example's 1.1.0. Its default `tendl-2025-neutron` bundle reports about
185.5 MiB download and 276.2 MiB installed. The
[pinned catalog](https://github.com/AvilaLabs/ACTINV/blob/d3dc4e463346c1a6d3dc426432048022932e14c5/crates/actinv-data/data/actinv-data-catalog-v1.2.0.json)
expects activation artifact `tendl-2025-neutron-709g` to have SHA-256
`ec4c72bf598dc8ad3d533d9cfafdcf493e2d1f949a3e4db6251495659b68cc44`.
Record all installed artifact hashes, not just that library's hash.
All six catalog artifacts were installed and hash-verified for the fresh
comparison below. The data and prepared solver cache stay outside the package.

## Fresh activation comparison

The pinned upstream example was executed with ACTINV 1.7.2, the rebuilt
Nucleide release wheel, NumPy 2.5.3, and catalog 1.2.0. Its compatibility
patch changes only the example's catalog-version check from `1.1.0` to
`1.2.0`; the activation, index, and two decay artifact hashes are identical
between these catalogs. No physical inputs or acceptance thresholds changed.

The case uses 15.748 g of Fe-56 in 2 cm³, a native 709-group neutron grid
with total flux 10¹² n/(cm²·s), 300 s irradiation, and 3600 s cooling.
The Nucleide-translated mesh calculation and a separately specified direct
ACTINV calculation produced identical 24-group photon strengths and Mn-56
activity: **463393807.7641661 photons/s** and **20622094.12329952 Bq/g**.
The comparison requires relative agreement within `1e-9`, with `1e-12`
absolute tolerance. ALARA export/readback reconstructed the integrated
photon total with relative error `1.2863e-16` (required: `1e-10`).
The text export contains 18 groups: its six trailing zero groups are omitted.

This establishes a fresh activation solve and agreement between two input
paths through the same solver. It does not independently validate the
evaluated nuclear data or experimental activation accuracy. See the
[validation instructions](../../../validation/README.md#fresh-actinv-activation-comparison)
for reproduction and provenance.

## Keep six exchange boundaries explicit

| Boundary | Required mapping |
| --- | --- |
| Zone and mixture | `alara_parse_deck(...)["mat_loading"]` preserves zone/mixture associations. Join explicit `volumes` by zone name; their list order need not match. Parsing exposes declarations; call `alara_validate_deck` separately to check references. |
| Composition and mass | A mixture reference is not an evaluated isotopic composition. Resolve external material libraries yourself, choose nuclide fractions and mass, and convert names such as `mn-56` to ACTINV's `Mn56`. Do not supply both an element and its explicit isotopes in the same composition. |
| Neutron flux | Supply physical group-integrated flux in n/(cm²·s), with energy edges and ordering matching the activation library. Select the interval after the deck's `skip` count and apply its `scale`. ACTINV's spectrum total must equal the intended physical sum. |
| Time | Expand irradiation/pulse schedules into durations and flux-on/off steps. ALARA cooling requests are times after shutdown; convert successive requests into incremental cooling durations. |
| Evaluated data | Provision and verify the actual activation/decay files separately. Retain software pins, data hashes, problem and result ledger with each calculation. |
| Photon readback | Select the requested cooling time and `TOTAL` rows, then multiply photon density by activation-zone volume exactly once to obtain integrated photons/s before voxel splitting. |

A 175-group ALARA spectrum is not a 709-group spectrum. The
`nucleide.alara.regroup_flux(..., model="flat_energy")` helper can
conserve integrated flux using explicit ascending eV edges and constant
differential flux within each source interval. Splitting a coarse group
can change reaction rates despite conserving the total. Use matching
native neutron edges for the activation comparison; a conservative
regrouping alone does not establish equivalent activation input.

## Parse a small synthetic input

This example uses a two-group spectrum solely to show interval selection,
scaling and time conversion. It is not a spectrum for the 709-group library.

```python
from nucleide.alara import (
    alara_expand_schedule,
    alara_parse_deck,
    alara_parse_flux,
    alara_validate_deck,
)

# FE56 is an external material-library reference, not a resolved composition.
deck_text = """geometry rectangular
volume
2 zone1
end
mixture iron
material FE56 1 1
end
mat_loading
zone1 iron
end
flux neutron synthetic.flux 2.5 1 default
schedule exposure
300 s neutron steady 0 s
end
pulsehistory steady
1 0 s
end
cooling
600 s
3600 s
end
"""
deck = alara_parse_deck(deck_text)
alara_validate_deck(deck_text)
mixture_by_zone = {row["zone"]: row["mixture"] for row in deck["mat_loading"]}
volume_by_zone = {row["zone"]: row["volume_cm3"] for row in deck["volumes"]}
assert mixture_by_zone["zone1"] == "iron"
assert volume_by_zone["zone1"] == 2.0

spectra = alara_parse_flux("100 200\n\n4 8\n", "synthetic")
flux_definition = deck["fluxes"][0]
selected = spectra["intervals"][flux_definition["skip"]]
physical_flux = [value * flux_definition["scale"] for value in selected]
assert physical_flux == [10.0, 20.0]

steps = alara_expand_schedule(deck_text)
schedule = [
    {"dt": step["duration_s"], "flux": 0.0 if step["is_cooling"] else 1.0}
    for step in steps
]
previous = 0.0
for cooling_time in deck["cooling_times_s"]:
    schedule.append({"dt": cooling_time - previous, "flux": 0.0})
    previous = cooling_time
assert [step["dt"] for step in schedule] == [300.0, 600.0, 3000.0]
assert sum(step["dt"] for step in schedule) == 3900.0
```

The two requested cooling times are 600 s and 3600 s after the same
shutdown at absolute time 300 s. Adding cooling durations of 600 s and
3600 s would select the wrong final time. For multiple activation zones,
select each zone's corresponding interval after the skip offset; do not
reuse the first selected interval for every zone.

For a real `actinv-spec-1` problem, start from `actinv new problem.json`
and replace its composition, mass, matching-grid spectrum and schedule
with reviewed inputs. The pinned upstream adapter describes the full
translation. Nucleide does not own an activation driver or resolve
ACTINV's data paths.

## Validate, solve and export externally

These commands require the optional solver and verified nuclear data;
they are not run by the synthetic example above:

```bash
actinv new problem.json
actinv validate problem.json --hashes
actinv run problem.json result.json
```

Keep the result's `ledger`, including missing reaction/product coverage.
For a mesh calculation, export the selected step and pass the physical
shutdown time (replace `STEP` with its step index):

```bash
actinv export-r2s mesh-result.ndjson STEP r2s-source.ndjson
actinv export-source alara r2s-source.ndjson photon-src --shutdown-t-s 300
```

An absolute source time of 3900 s then becomes 3600 s cooling in
`actinv-alara-index.json` and the emitted `.photonSrc` rows. Read each
index entry's file and `volume_cm3` together. The released CLI's
synthetic 2 cm³ export reconstructed 12 photons/s from densities `[2, 4]`;
this checks the exchange format, not an activation solve.

## Read and split a synthetic photon source

This runs without ACTINV or nuclear data. Distinct cooling rows make
accidental time selection visible. The selected density is
6 photons/(s·cm³); multiplying by 2 cm³ gives 12 photons/s:

```python
from nucleide.r2s import photon_group_sums, tag_zone_strength

photon_text = "TOTAL 600 s 8 4\nTOTAL 3600 s 2 4\n"
source = photon_group_sums(photon_text, ["TOTAL"], 3600.0)
assert source["sums"] == [2.0, 4.0]
assert source["total"] == 6.0
zone_total_photons_s = source["total"] * 2.0

tags = tag_zone_strength(
    [zone_total_photons_s],
    [0, 0],
    split=True,
    voxel_volumes_cm3=[0.5, 1.5],
)
assert tags["source_strength"] == [3.0, 9.0]
assert tags["total"] == 12.0
```

Volume-weighted splitting assumes uniform source density within each
activation zone and returns integrated voxel strengths. Volumes must be
finite and positive; every positive-strength zone must have a voxel.
Without volumes, `split=True` retains the existing equal-count split.
Tag decay-time metadata remains shutdown `0.0`; retain the selected
cooling time separately when using these tags downstream.

The exporter supplies **TOTAL-only photon spectra** on energy-centroid
labels. Centroids do not determine bin edges. A transport source needs
reviewed photon edges and an explicit mapping that preserves the groups'
meaning. An export containing 18 nonzero centroids is not automatically
a complete 24-group photon representation; matching its total is
insufficient to establish a dose result.

## See also

- [Activation analysis](activation-analysis.md) for the deck, flux,
  schedule and R2S workflow APIs.
- ACTINV's [pinned specification](https://github.com/AvilaLabs/ACTINV/blob/d3dc4e463346c1a6d3dc426432048022932e14c5/docs/guide/specification.md)
  and [iron example](https://github.com/AvilaLabs/ACTINV/tree/d3dc4e463346c1a6d3dc426432048022932e14c5/examples/fns_iron).
- [Python API reference](../../reference/python-api.mdx#nucleidealara)
  for the ALARA and R2S signatures.
