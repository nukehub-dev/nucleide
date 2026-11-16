# `crates/r2s` AGENTS.md

## Purpose

Rigorous two-step (R2S) shutdown-dose-rate orchestration: ALARA deck
workflows, versionless ARMI snapshot bridge, photon assembly, flux
normalization, activation-to-transport source handoff, voxel tags,
parameter-sweep expansion (WATTS-class), and facility-flow accounting.
Builds decks, links zones, assembles sources; transport and activation
solving stay inside their respective codes.

## Ownership

Owns `crates/r2s/src/` (`workflow.rs`, `tags.rs`, `photon.rs`, `source.rs`,
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
- `source`: activation-to-transport handoff. `normalize_flux` takes an
  explicit `TallyNormalization`; per-source volume-integrated tallies
  multiply by the source rate and divide by caller voxel volumes,
  per-source volume-averaged tallies multiply by the source rate only, and
  the `*RateScaled` variants skip the rate multiplication because OpenMC
  fixed-source tallies already include the configured total source strength.
  `Physical` input is rejected. All tally/volume/rate inputs and
  intermediate/result values must be finite; inputs are never mutated.
  `PhotonSource` is the dense structured-mesh
  contract: C-order voxels (`(i*ny+j)*nz+k`), finite ascending mesh bounds,
  positive finite voxel volumes, `G+1` strictly ascending finite bin edges
  with first `>= 0` (centroids are never accepted), finite `cooling_s >= 0`,
  finite nonnegative integrated photons/s with one row per voxel and `G`
  entries per row. Totals and densities use Neumaier compensation;
  `joint_probabilities` is the joint voxel/group distribution and rejects a
  non-finite or `<= 0` total; `energy_probabilities` is the per-voxel
  spectrum (zero rows for zero-source voxels). `from_exported_groups` requires
  strictly increasing unique in-range indices and zero-fills omitted bins.
  `select_cooling_indices` matches exact `==` only, preserves requested
  order, and rejects duplicates/missing offsets loudly.
- `source` diagnostics: `PhotonSourceDraft` is unvalidated; `SourceReport`
  and `PhotonSource::from_draft` share one check order (mesh, bounds,
  volumes, edges, cooling, strengths, totals, warnings, assumptions), never
  repair, and use stable kebab-case codes documented in `source.rs`.
  `SOURCE_ISSUE_LIMIT` caps stored examples while `issue_count` and
  `truncated` keep the true count. Assumptions (`source-assumption-*`) and
  `source-all-zero` do not invalidate a source; any `Error` severity does.
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
- `cargo bench -p nucleide-r2s --bench source_bench` for the 1k/10k/100k
  voxel source-contract baseline (24 groups; no I/O).
- `pytest tests/test_r2s.py tests/test_r2s_tags.py tests/test_openmc_handoff.py`
  after `maturin develop`.
- `validation/activation_vs_refs.py` runs inside `run_all.sh`.

## Child NAD Index

None.
