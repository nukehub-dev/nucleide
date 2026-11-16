---
title: Activation to transport
sidebar:
  order: 24
---

This tutorial connects three things that are easy to mix up: a raw OpenMC
tally, a physical group flux, and an activation-derived photon source that a
transport code can sample. Nucleide keeps the raw data, the physical
conversion and the source contract separate so a downstream fold cannot
silently use the wrong units.

The example is synthetic. No activation solve, transport run or evaluated
data is required.

## 1. Declare the raw tally normalization

OpenMC's default track-length `flux` tally is volume-integrated. In a
fixed-source run it is normalized by the simulated particle count **and**
scaled by the sum of the configured source strengths. With the default source
strength (`1.0`, or normalized selection weights), the result is per source
particle; with source strengths written as absolute particles/s or photons/s,
the source rate is already inside the tally. A random-ray tally run with
`settings.random_ray["volume_normalized_flux_tallies"] = True` is
volume-averaged instead. The statepoint bridge cannot infer which combination
a caller produced, so it records `normalization = "unknown"` and the
conversion asks for a declaration.

```python
import json
import tempfile
from pathlib import Path

import nucleide.openmc as omc
import nucleide.r2s as r2s

raw = {
    "tally_id": 7,
    "score": "flux",
    "kind": "cell",
    "n_voxels": 1,
    "n_groups": 2,
    "energy_bounds_mev": [0.0, 1.0, 2.0],
    "flux": [[4.0, 9.0]],
    "rel_err": [[0.05, 0.10]],
    "mesh_shape": None,
    "mesh_bounds": None,
    "openmc_version": "0.16.0",
    "normalization": "unknown",
    "units": "unknown",
}

physical = omc.normalize_flux(
    raw,
    source_rate_s=10.0,
    volumes_cm3=[2.0],
    input_normalization="per-source-volume-integrated",
)
assert physical["flux"] == [[20.0, 45.0]]
assert physical["units"] == "1/(cm2 s)"
assert raw["flux"] == [[4.0, 9.0]]
```

The arithmetic is `raw × source_rate / volume` for per-source
volume-integrated input. Volume-averaged input applies only `source_rate`.
When OpenMC already applied the total source strength, declare
`rate-scaled-volume-integrated` (divide by volume only) or
`rate-scaled-volume-averaged` (no arithmetic); multiplying those by the rate
again would double-count it. Physical input is rejected rather than scaled
twice. Relative errors are unchanged, and source-rate or volume uncertainty is
not folded in.

```python
rate_scaled = omc.normalize_flux(
    raw,
    source_rate_s=10.0,
    volumes_cm3=[2.0],
    input_normalization="rate-scaled-volume-integrated",
)
assert rate_scaled["flux"] == [[2.0, 4.5]]
```

A cell tally uses one volume. For a mesh tally, pass one positive volume per
voxel in the bridge's C-order (x slowest, then y, then z).

## 2. Keep the metadata through a CSV roundtrip

Legacy `nucleide-openmc-tally-v1` files carry no normalization state and read
back as `unknown`; they are never relabelled as physical flux by default. New
files use v2 and retain the normalization, units and optional source rate.

```python
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "flux.csv"
    omc.tally_to_csv(physical, path)
    readback = omc.tally_from_csv(path)
assert readback["normalization"] == "physical"
assert readback["source_rate_s"] == 10.0
assert readback["flux"] == [[20.0, 45.0]]
```

## 3. Build an explicit photon source

The R2S source contract takes a regular mesh, ascending eV group edges, one
cooling time and **integrated** photons/s per voxel and group. The mesh shape
`(2, 1, 1)` below means two voxels in x; each voxel has its own spectrum.

```python
source = r2s.photon_source(
    mesh_shape=[2, 1, 1],
    lower_left_cm=[0.0, 0.0, 0.0],
    upper_right_cm=[2.0, 1.0, 1.0],
    voxel_volumes_cm3=[1.0, 1.0],
    edges_ev=[0.0, 1.0e6, 2.0e6, 5.0e6],
    cooling_s=0.0,
    strengths_photons_per_s=[[100.0, 0.0, 20.0], [0.0, 300.0, 0.0]],
)

assert source["total_photons_per_s"] == 420.0
assert source["voxel_totals_photons_per_s"] == [120.0, 300.0]
assert r2s.source_density(source) == [[100.0, 0.0, 20.0], [0.0, 300.0, 0.0]]
```

A density applies each voxel volume exactly once. Do not multiply an
integrated source by volume again when handing it to a transport code. The
exportable example uses uniform 1 cm³ regular-mesh cells; the contract also
accepts caller volumes for density conversion, but OpenMC export requires them
to match the regular-mesh cell volume.

The joint voxel/group probabilities are derived from the same strengths, not
from independent spatial and energy marginals:

```python
probabilities = r2s.source_probabilities(source)
assert probabilities == [[100.0 / 420.0, 0.0, 20.0 / 420.0], [0.0, 300.0 / 420.0, 0.0]]
```

## 4. Inspect validation and assumptions

The report uses stable issue codes and deterministic ordering. Assumptions
are labelled separately from errors, and invalid data never gets repaired.

```python
report = r2s.source_report(source)
assert report["valid"] is True
codes = [issue["code"] for issue in report["issues"]]
assert "source-assumption-uniform-within-bin" in codes
assert "source-assumption-uniform-within-voxel" in codes
```

The two modelling assumptions are deliberate: photons are uniform in energy
density inside each group (per eV) and uniform in position inside each voxel.
They conserve group and voxel strengths; they do not reconstruct unresolved
emission lines or sub-voxel structure. A zero-source voxel is valid inside a
nonzero source; an all-zero distribution is valid data but cannot be sampled
as a probability distribution.

## 5. Restore omitted groups explicitly

Activation exports sometimes drop zero groups. Pass the original full grid and
the dense indices of the exported columns; omitted groups become exactly zero.
Centroids are never interpreted as bin edges.

```python
restored = r2s.photon_source_exported(
    mesh_shape=[1, 1, 1],
    lower_left_cm=[0.0, 0.0, 0.0],
    upper_right_cm=[1.0, 1.0, 1.0],
    voxel_volumes_cm3=[1.0],
    full_edges_ev=[0.0, 1.0e6, 2.0e6, 5.0e6],
    cooling_s=0.0,
    exported_group_indices=[1, 2],
    exported_strengths_photons_per_s=[[5.0, 7.0]],
)
assert restored["strengths_photons_per_s"] == [[0.0, 5.0, 7.0]]
```

Duplicate, out-of-range or unsorted indices raise clear errors.

## 6. Export a cooling-time batch and verify it

Batch selection is exact. There is no interpolation, no nearest-time
substitution and no silent duplicate collapse. Distinct float spellings keep
close-but-distinct offsets in separate files. The batch manifest is written
only after every requested output succeeds.

```python
def write_json_source(item: dict, path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "cooling_s": item["cooling_s"],
                "total_photons_per_s": item["total_photons_per_s"],
            }
        ),
        encoding="utf-8",
    )

hour = r2s.photon_source(
    mesh_shape=[2, 1, 1],
    lower_left_cm=[0.0, 0.0, 0.0],
    upper_right_cm=[2.0, 1.0, 1.0],
    voxel_volumes_cm3=[1.0, 1.0],
    edges_ev=[0.0, 1.0e6, 2.0e6, 5.0e6],
    cooling_s=3600.0,
    strengths_photons_per_s=[[7.0, 0.0, 1.0], [0.0, 13.0, 0.0]],
)

with tempfile.TemporaryDirectory() as tmp:
    batch = omc.export_cooling_batch(
        [source, hour],
        requested_times_s=[0.0, 3600.0],
        out_dir=tmp,
        writer=write_json_source,
    )
    manifest = omc.read_manifest(Path(tmp) / "batch-manifest.json")
    check = omc.verify_manifest(manifest, base_dir=tmp)

assert batch["manifest"]["status"] == "complete"
assert [entry["time_s"] for entry in batch["outputs"]] == [0.0, 3600.0]
assert check["ok"] is True
```

A failed writer raises `nucleide.openmc.BatchExportError` with the partial
outputs attached; no complete manifest is claimed. Existing outputs are
refused unless `overwrite=True`. An all-zero time is an error unless
`skip_zero=True`, which records the skipped time in the manifest. The batch
manifest records its schema, requested offsets, units, voxel order,
assumptions, and each output's SHA-256, total photons/s, energy-grid hash and
validation summary; skipped entries are listed separately.

`source_manifest` records the measured input/output hashes, the energy-grid
hash, voxel order, cooling offset and validation summary in a compact JSON
document that can be checked offline:

```python
manifest = omc.source_manifest(source=source, input_files=[])
assert manifest["source"]["voxel_order"] == "c-order-(i,j,k)"
assert len(manifest["source"]["energy_grid_sha256"]) == 64
```

Caller-supplied evaluated-library identity is recorded as such, without an
invented hash. Timestamps and environment variables are deliberately absent
from the deterministic content.

## 7. Optional OpenMC mesh source

If caller-side OpenMC 0.16.0 is installed, the same source becomes a
`MeshSource` with one energy distribution per mesh element. The adapter maps
the C-order voxel list into OpenMC's structured-mesh order explicitly and
stores absolute photons/s as the element strengths. OpenMC then scales its own
tallies by that total strength, so a subsequent flux conversion declares a
`rate-scaled-*` input; do not multiply the same rate again.

<!-- code-test: skip -->

```python
if omc.openmc_version() is not None:
    mesh_source = omc.mesh_source(source)
    assert mesh_source.strength == 420.0
    with tempfile.TemporaryDirectory() as tmp:
        omc.write_mesh_source_xml(source, Path(tmp) / "settings.xml")
```

OpenMC is optional and imported lazily. The generic tutorial above runs
without it. Exporting a source is not transport validation: source
conservation, XML/API compatibility and spatial-energy sampling are separate
claims from shutdown dose rate, which needs geometry, evaluated photon data,
a dose response and an independent reference.

The validation harness carries the external checks separately:
`validation/openmc_handoff_probe.py` pins the OpenMC 0.16 source/tally/XML
semantics, and `validation/openmc_real_data_probe.py` runs a bounded
fixed-source water/iron sphere with a caller-mounted NNDC HDF5 library and
checks the same per-source and rate-scaled conversions through
`normalize_flux`. The latter is a normalization and data-loading check, not a
dose or transport benchmark.
