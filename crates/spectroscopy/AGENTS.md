# `crates/spectroscopy` AGENTS.md

## Purpose

Gamma-ray spectroscopy and measurement: smoothing, peak counting,
energy/efficiency calibration, X-ray line algebra, caller-line SDEF
decay-source cards, and the two text `.spe` readers. Pure algorithms over
caller vectors; no tabulated data, no I/O beyond `&str` parsing.

## Ownership

Owns `crates/spectroscopy/src/` (`spectrum.rs`, `smooth.rs`, `counts.rs`,
`calib.rs`, `spe.rs`, `xray.rs`, `sdef.rs`, `lines.rs`, `error.rs`), the Python surface
(`nucleide.spectroscopy`, `spectroscopy_*` in `_internal`),
`tests/test_spectroscopy.py`, the `fixtures/spectroscopy/` synthetic
oracles, and the `spectroscopy_vs_pyne.py` validation oracle.

## Local Contracts

- Equation set (E1–E9) is pinned: rectangular/five-point smoothing,
  `m == 1` background, half-open gross counts, net counts, quadratic
  energy bins, two log-polynomial efficiency fits, caller-constant X-ray
  algebra, caller-line SDEF decay-source cards, dollar/plain `.spe`
  grammars. Derivations live in `docs/theory/spectroscopy.mdx`; code
  comments cite equation labels, never external paths.
- Pin the upstream quirks, never "fix" them: dollar duplicate-tag
  first-wins, missing-tag errors, `$MEAS_TIM:` live-real reversal,
  `$DATA:` last-plus-one channel count, triplet-plus-`keV` lines,
  `$ROI:`/`$PRESETS:`/`$ENER_FIT:` ignored, plain unknown-key ignore,
  double-space `Energy Fit:` indices 0/2/4, positional 0-based dollar
  labels under nonzero `start_chan_num`, positional E3–E5 indexing.
- Caller constants only: atomic yields/ratios/energies, calibration fits,
  efficiency coefficients, and decay lines (energy + intensity pairs, from
  lists or the runtime TSV interchange in `lines.rs`) are inputs. No
  vendored atomic tables or decay-line libraries, no coefficient fitting,
  no FWHM evaluation, no legacy line-list reader.
- SDEF decay sources (E9): caller lines merge at duplicate energies, sort
  ascending, and normalize to probabilities summing to 1.0; cards keep the
  upstream monoenergetic field order (single line → inline `ERG=<E>`,
  multiple → `ERG=D1` + `SI1 L`/`SP1 D`, wrapped at 80 columns). The
  distribution syntax is parser-verified surface only — MCNP sampling
  semantics stay the caller's responsibility. `PAR=` designators come from
  `nucleide-nuclei` (this crate's one workspace dependency); versions 5/6
  only.
- Synthetic fixtures only: hand-built counts plus closed-form values with
  recorded provenance; never laboratory or evaluated-library data. New
  oracles extend `fixtures/spectroscopy/` and the fixture tests in the
  same change.
- Out of scope (do not expand here): peak search/fit, activities, plotting.
- WASM/tutorial surface (owned): `spectroscopySmooth` in `bindings/wasm`
  (E1 rectangular / E2 five-point smoothing plus E3–E5 gross/background/net
  counting; thin facade, positional channels), the `SpectroscopyDemo` demo
  in `website/src/components/interactive/SpectroscopyDemo.tsx`, and the
  `docs/tutorials/interactive/spectroscopy.mdx` page. Demo defaults stay
  synthetic (`fixtures/spectroscopy/` values or hand-picked small numbers).

## Work Guidance

- New analyses add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate depends only on `nucleide-nuclei` inside
  the workspace (particle dialect for E9); bindings depend on it, never the
  reverse.
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
