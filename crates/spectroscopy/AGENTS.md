# `crates/spectroscopy` AGENTS.md

## Purpose

Gamma-ray spectroscopy and measurement: smoothing, peak counting,
energy/efficiency calibration, X-ray line algebra, caller-line SDEF
decay-source cards, and the two text `.spe` readers. Pure algorithms over
caller vectors; no tabulated data, no I/O beyond `&str` parsing.

## Ownership

Owns `crates/spectroscopy/src/` (`spectrum.rs`, `smooth.rs`, `counts.rs`,
`calib.rs` including the E7-fit coefficient fit, `spe.rs`, `xray.rs`,
`sdef.rs`, `lines.rs`, `error.rs`), the Python surface
(`nucleide.spectroscopy` including `fit_efficiency`, `spectroscopy_*` in
`_internal`), `tests/test_spectroscopy.py`, the `fixtures/spectroscopy/`
synthetic oracles (including `efficiency_fit.json`), and the
`spectroscopy_vs_pyne.py` validation oracle.

## Local Contracts

- Equation set (E1–E9) is pinned: rectangular/five-point smoothing,
  `m == 1` background, half-open gross counts, net counts, quadratic
  energy bins, two log-polynomial efficiency evaluations plus the E7-fit
  log-space coefficient fit (an E7 extension, not E10 unless the equation
  set is amended), caller-constant X-ray algebra, caller-line SDEF
  decay-source cards, dollar/plain `.spe` grammars. Derivations live in
  `docs/theory/spectroscopy.mdx`; code comments cite equation labels, never
  external paths.
- Pin the upstream quirks, never "fix" them: dollar duplicate-tag
  first-wins, missing-tag errors, `$MEAS_TIM:` live-real reversal,
  `$DATA:` last-plus-one channel count, triplet-plus-`keV` lines,
  `$ROI:`/`$PRESETS:`/`$ENER_FIT:` ignored, plain unknown-key ignore,
  double-space `Energy Fit:` indices 0/2/4, positional 0-based dollar
  labels under nonzero `start_chan_num`, positional E3–E5 indexing.
- Caller constants only: atomic yields/ratios/energies, calibration fits,
  decay lines (energy + intensity pairs, from lists or the runtime TSV
  interchange in `lines.rs`), and efficiency-fit points/weights are inputs.
  No vendored atomic tables or decay-line libraries, no FWHM evaluation, no
  legacy line-list reader. Coefficient fitting exists on exactly one path:
  the E7-fit weighted least-squares fit (`fit_efficiency`, targets
  `y = ln eff` over the `(ln E)^j` / `(1/E)^j` basis with caller-supplied
  weights, solved through the `nucleide-linalg` `lstsq` kernel); every other
  coefficient vector stays a caller input with no fitting routine.
- SDEF decay sources (E9): caller lines merge at duplicate energies, sort
  ascending, and normalize to probabilities summing to 1.0; cards keep the
  upstream monoenergetic field order (single line → inline `ERG=<E>`,
  multiple → `ERG=D1` + `SI1 L`/`SP1 D`, wrapped at 80 columns). The
  distribution syntax is parser-verified surface only — MCNP sampling
  semantics stay the caller's responsibility. `PAR=` designators come from
  `nucleide-nuclei` (particle dialect for E9); versions 5/6 only.
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
- Keep the layering: this crate depends only on `nucleide-nuclei` (particle
  dialect for E9) and `nucleide-linalg` (least-squares kernel for the E7
  fit) inside the workspace; bindings depend on it, never the reverse.
  Never implement a second least-squares route here — the shared
  `nucleide-linalg` `lstsq` module is the one dense-real path.
- Validation oracle needs the upstream package at run time
  (container PyNE layer); a missing oracle is a recorded SKIP, never a
  silent pass. The E7 fit has no upstream counterpart (no fitting routine
  upstream), so its tier-2 row is always a loud SKIP pinned by synthetic
  recovery + round-trip gates.

## Verification

- `cargo test -p nucleide-spectroscopy` (includes fixture replay, among
  them `efficiency_fit.json` recovery + round-trip at 1e-9).
- `pytest tests/test_spectroscopy.py` after `maturin develop`.
- `validation/spectroscopy_vs_pyne.py` runs inside `run_all.sh` (two-tier:
  synthetic gates always, PyNE cross-check when importable; E7-fit tier-2
  always SKIP).

## Child NAD Index

None.
