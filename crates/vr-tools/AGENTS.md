# `crates/vr-tools` AGENTS.md

## Purpose

Monte Carlo variance-reduction tools: MAGIC weight windows, mesh source
sampling (alias tables), and Gaussian KDE source resampling
(KDSource-class, clean-room). No transport, no MOAB.

## Ownership

Owns `crates/vr-tools/src/` (`magic.rs`, `windows.rs`, `sampling.rs`,
`kde.rs`), the Python surface (`nucleide.vr`, `magic`/`magic_with`/
`emit_openmc_weight_windows`/`emit_serpent_wwin`/`AliasTable`/
`MeshSourceSampler`/`KdeSampler` in `_internal`), and the MAGIC, sampling,
and KDE sections of `tests/test_serpent_fluka_vr.py`.

## Local Contracts

- `magic`: flux/`(2*max)` weight-window bounds over meshtal data.
- `windows`: pure formatting of `MagicOutput` for other transport codes —
  OpenMC `settings.xml` (`<mesh>` + `<weight_windows>` sub-element spelling
  pinned against the public OpenMC docs §3.29/§3.66 and the C++ reader) and
  Serpent via the MCNP WWINP text spelling (`wwin <name> wf "<file>" 2`,
  user guide §2.2.8.3; the Serpent-native FMT=1 layout is not publicly
  documented and is never emitted). Bounds are reordered z-fastest
  (meshtal) → x-fastest (target codes), energies MeV → eV for OpenMC only.
  Drift notes cover synthesized upper bounds and null cells; negative,
  non-finite, unsorted, or out-of-range inputs are loud named errors. The
  Serpent text re-parses through `nucleide-mcnp-io`'s WWINP reader, which
  pins its layout — do not hand-edit that stream format here.
- `sampling`: Walker/Vose alias tables + `MeshSourceSampler`
  (ANALOG/UNIFORM/USER bias modes, birth weights); draws take
  caller-supplied uniforms, no RNG inside.
- `kde`: axis-aligned Gaussian KDE over caller particle vectors
  (Silverman, fixed, or adaptive per-sample kNN bandwidths; zero-variance
  dims and duplicate-point zero distances are loud); deterministic
  `draw(u, normals)`/`pdf`, no RNG inside, no `mcpl-io` dependency
  (MCPL projection stays caller-side).
- KDE equations are textbook statistics implemented from scratch; no
  upstream code or data read, ported, or vendored (GPL upstream).
- Synthetic fixtures only; golden values live in unit tests with
  recorded provenance (test-only LCG seed).
- Bindings stay thin; error enum stays `#[non_exhaustive]`.

## Work Guidance

- New samplers add a module plus Python/tests in the same change.
- Keep the layering: depends on `nucleide-mcnp-io` (meshtal) only.

## Verification

- `cargo test -p nucleide-vr-tools` (includes KDE recovery gates).
- `pytest tests/test_serpent_fluka_vr.py` after `maturin develop`.

## Child NAD Index

None.
