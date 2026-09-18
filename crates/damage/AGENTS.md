# `crates/damage` AGENTS.md

## Purpose

Damage and gas-production metrics by spectral folding: the closed-form
NRT-dpa and arc-dpa displacement functions, He/H production in appm, He/dpa
ratios, and MVN uncertainty propagation through the folds, plus coil
fast-fluence / lifetime bookkeeping (fast-flux sums, history accumulation,
weakest-link life over caller-supplied limit tables) reusing those folds
as the dpa spectral weight. The first-wall and magnet engineering metric
set for the fusion workstream.

## Ownership

Owns `crates/damage/src/` (`physics.rs`, `fold.rs`, `coil.rs`, `uq.rs`,
`error.rs`, `specter.rs`, `data/specter_table_vii.tsv`), the Python surface
(`nucleide.damage`, `damage_*` in `_internal`), and the
`validation/damage_vs_specter.py` oracle.

## Local Contracts

- Equation set is pinned and cited in rustdoc: NRT-dpa (Norgett, Robinson &
  Torrens, NED 33 (1975) 50–54) with the Lindhard partition via the
  Robinson fit (`k_L = 0.0793…`, `E_L = 30.724… eV`, `g(ε) = ε + 3.4008
  ε^{1/6} + 0.40244 ε^{3/4}` — constants as documented in Griffin
  SAND2016-2269 §3.2, which corrects the NJOY manual's `k_L` exponent);
  arc-dpa (Nordlund et al., Nat. Commun. 9 (2018) 1084, CC BY 4.0,
  Eqs. (5)–(7)). Piecewise boundaries sit on the PKA energy at `E_d` and
  `2E_d/κ`; the high branch evaluates `κ·T_dam/(2E_d)` on the Lindhard
  damage energy. Equations are facts: cite papers, never external paths.
- Clean-room provenance: no ported code or tables. The fispact-org PKA
  evaluator is GPL-3.0 — never read. ASTM E693/E521 are paywalled:
  designation-string references only, never transcribed.
- SPECTER (ANL/FPP/TM-197, US-gov PD) is the validation oracle AND the
  source of one opt-in vendored fallback: `validation/damage_vs_specter.py`
  transcribes a handful of output spots (Tables V/VI, HFIR-CTR32) and folds
  them at report print precision (5 digits, tol 1e-4 relative), while
  `specter.rs` + `data/specter_table_vii.tsv` vendor Table VII
  (spectrum-averaged damage-energy XS in keV-b for 24 elements × 7 spectra)
  plus the Table II `E_d` column — displacement XS only, EU-table pattern
  (`try_` + cache split, row-count pin, spot tests). The HFIR column
  reproduces the oracle's Fe/Ti/Cu spots at Table VII print precision
  (3 s.f., tol 5e-3). Scope narrowed deliberately: the group-wise
  Appendix A printouts are image-only pages in the report scan and cannot
  meet the transcription bar, so the table covers the 24 Table VII
  elements, not the 41 Table I entries; gas-production columns stay out.
  The fallback is never consulted implicitly — fold kernels only see
  caller slices.
- Fold conventions: piecewise-constant per group — group flux is the
  per-group integrated flux, response the caller-condensed group cross
  section (barns); bounds are validated (`G+1`, strictly increasing MeV)
  and document the grid but never enter the sum. Zero-flux groups
  contribute exactly 0.0; negative flux/response is a loud error; the
  He/dpa ratio at zero dpa is `Error::ZeroDpa`, never `inf` (the
  `ZeroMaxFlux` precedent).
- Material constants (`E_d`, `b_arc`, `c_arc`) are caller-supplied;
  `ArcParams::new` enforces the MD-fitted domain (`b_arc < 0`,
  `0 < c_arc < 1`). This crate ships no fitted parameter table.
- UQ reuses `linalg::sample` exclusively (seeded MVN, relative
  perturbations, `k`-SE gates against the exact bilinear expectation and
  the first-order propagated standard deviation). No new sampling
  machinery; UQ on the He/dpa ratio is a named-open (`NotYetSupported`).
- Out of scope (do not expand here): PKA-spectra solving, transport
  solving, group-wise displacement-table vendoring (Appendix A stays
  image-only), gas-production columns, IAEA/IRDFF/TENDL consumption
  (caller-supplied or runtime-download only), magnetics/quench/structural
  analysis, and vendored coil limit tables (limits stay caller-supplied;
  published design numbers are gates, never defaults).
- Coil bookkeeping (`coil.rs`) is arithmetic over caller spectra reusing
  the fold conventions: fast fluence sums groups whose upper edge clears
  the caller threshold (whole-group inclusion when the threshold cuts a
  group), dpa rates come from the landed folds at one second, and life is
  the weakest-link `min` of limit/rate (remaining life clamps at zero;
  all-zero rates are infinite life, never `inf` from a division).

## Work Guidance

- Keep the layering: depends on `nucleide-linalg` + `nucleide-nuclei`
  only; bindings depend on it, never the reverse.
- New metrics add a fold fn plus Python facade/tests/validation rows in
  the same change (bindings thin, no logic).
- New transcribed oracle spots must cite report table/page, verify
  against the runtime-fetched PDF text (OCR-tolerant matching), and keep
  the internal arithmetic check (`appm = Φ·σ·1e-18`) green.

## Verification

- `cargo test -p nucleide-damage` (closed-form limits, piecewise gates,
  fold conventions, UQ k-SE gates, Table VII row-count/spot gates plus the
  HFIR-column oracle reproduction).
- `pytest tests/test_damage.py` after `maturin develop`.
- `validation/damage_vs_specter.py` runs inside `run_all.sh` (two-part:
  analytic gates always, SPECTER cross-check when the hash-pinned PDF and
  `pdftotext` are available — loud SKIP otherwise).

## Child NAD Index

None.
