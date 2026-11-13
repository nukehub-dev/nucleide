# `crates/r2s` AGENTS.md

## Purpose

Rigorous two-step (R2S) shutdown-dose-rate orchestration: ALARA deck
workflows, versionless ARMI snapshot bridge, photon assembly, voxel tags,
parameter-sweep expansion (WATTS-class), and facility-flow accounting.
Builds decks, links zones, assembles sources; transport and activation
solving stay inside their respective codes.

## Ownership

Owns `crates/r2s/src/` (`workflow.rs`, `tags.rs`, `photon.rs`,
`snapshot.rs`, `sweep.rs`), the Python surface (`nucleide.r2s`), and
`tests/test_r2s.py` + `tests/test_r2s_tags.py`.

## Local Contracts

- `workflow`: zone-to-flux linking, schedule expansion, per-step deck
  emission; uniform-split photon assembly.
- `tags`: legacy count-based splitting keeps its signed-tag and unmapped-zone
  semantics. `split_zone_totals_by_volume` splits integrated photons/s by
  caller-supplied positive voxel volumes in cm³, assuming uniform density
  within each zone; returns integrated photons/s in caller voxel order at
  shutdown (`0.0`). Require finite nonnegative zone totals, finite positive
  volumes and per-zone volume sums, matching lengths and valid indices;
  reject unmapped positive sources. Empty meshes allow only zero sources.
  Python `tag_zone_strength(..., voxel_volumes_cm3=...)` requires `split=True`.
  Compensate source/volume sums and the returned Python total; preserve legacy
  total getters; `checked_total_strength` provides an additive finite aggregate.
  Underflowing strengths can round to zero, so relative
  conservation cannot be promised below the representable range.
- `snapshot`: ARMI dict-in bridge (never touches HDF5, never mirrors
  ARMI internals); `snapshot_inventory` totals atoms per bare name
  (`N × V × 1e-24`) for facility snapshot differencing.
- `sweep`: named-axis validation, deterministic cartesian `expand_sweep`,
  `assemble_results` requiring exact per-case coverage. Pure data
  plumbing: no I/O, no template engine, no code execution.
- Snapshot compositions follow the emit ARMI-input rule (post-expansion
  nuclide keys); invalid densities/volumes are loud errors.
- Synthetic fixtures only; bindings stay thin.

## Work Guidance

- New analyses add a module plus Python/tests in the same change.
- Keep the layering: depends on `alara-io`/`mcnp-io`/`nuclei`/`vr-tools`;
  never the reverse, never on bindings.

## Verification

- `cargo test -p nucleide-r2s`.
- `pytest tests/test_r2s.py tests/test_r2s_tags.py` after `maturin develop`.
- `validation/activation_vs_refs.py` runs inside `run_all.sh`.

## Child NAD Index

None.
