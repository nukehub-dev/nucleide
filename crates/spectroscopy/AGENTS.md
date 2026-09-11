# `crates/spectroscopy` AGENTS.md

## Purpose

Gamma-ray spectroscopy and measurement: smoothing, peak counting,
energy/efficiency calibration, X-ray line algebra, and the two text `.spe`
readers. Pure algorithms over caller vectors; no tabulated data, no I/O
beyond `&str` parsing.

## Ownership

Owns `crates/spectroscopy/src/` (`spectrum.rs`, `smooth.rs`, `counts.rs`,
`calib.rs`, `spe.rs`, `xray.rs`, `error.rs`), the Python surface
(`nucleide.spectroscopy`, `spectroscopy_*` in `_internal`),
`tests/test_spectroscopy.py`, the `fixtures/spectroscopy/` synthetic
oracles, and the `spectroscopy_vs_pyne.py` validation oracle.

## Local Contracts

- Equation set (E1–E8) is pinned: rectangular/five-point smoothing,
  `m == 1` background, half-open gross counts, net counts, quadratic
  energy bins, two log-polynomial efficiency fits, caller-constant X-ray
  algebra, dollar/plain `.spe` grammars. Derivations live in
  `docs/theory/spectroscopy.mdx`; code comments cite equation labels,
  never external paths.
- Pin the upstream quirks, never "fix" them: dollar duplicate-tag
  first-wins, missing-tag errors, `$MEAS_TIM:` live-real reversal,
  `$DATA:` last-plus-one channel count, triplet-plus-`keV` lines,
  `$ROI:`/`$PRESETS:`/`$ENER_FIT:` ignored, plain unknown-key ignore,
  double-space `Energy Fit:` indices 0/2/4, positional 0-based dollar
  labels under nonzero `start_chan_num`, positional E3–E5 indexing.
- Caller constants only: atomic yields/ratios/energies, calibration fits,
  and efficiency coefficients are inputs. No vendored atomic tables, no
  coefficient fitting, no FWHM evaluation.
- Synthetic fixtures only: hand-built counts plus closed-form values with
  recorded provenance; never laboratory or evaluated-library data. New
  oracles extend `fixtures/spectroscopy/` and the fixture tests in the
  same change.
- Out of scope (do not expand here): peak search/fit, activities, decay
  spectra/SDEF vectors, plotting, WASM exposure.

## Work Guidance

- New analyses add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate is independent; bindings depend on it,
  never the reverse.
- Validation oracle needs the upstream package at run time
  (container PyNE layer); a missing oracle is a recorded SKIP, never a
  silent pass.

## Verification

- `cargo test -p nucleide-spectroscopy` (includes fixture replay).
- `pytest tests/test_spectroscopy.py` after `maturin develop`.
- `validation/spectroscopy_vs_pyne.py` runs inside `run_all.sh` (two-tier:
  synthetic gates always, PyNE cross-check when importable).

## Child NAD Index

None.
