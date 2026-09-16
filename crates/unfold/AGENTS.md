# `crates/unfold` AGENTS.md

## Purpose

Neutron spectrum unfolding from activation-type measurements: adjust a
caller-supplied guess spectrum against measured detector/reaction rates
until the forward fold through the caller-supplied response matrix
reproduces them. Inverse-problem iteration, not a transport solve.

## Ownership

Owns `crates/unfold/src/` (`sandii.rs` SAND-II iterator, `staysl.rs`
STAYSL-class iterator, `gravel.rs` GRAVEL iterator, `lib.rs` forward
operator + `error.rs`), the Python surface
(`nucleide.unfold`, `unfold_*` in `_internal`), `tests/test_unfold.py`,
and the `validation/unfold_vs_analytic.py` oracle.

## Local Contracts

- Equation sets are pinned: SAND-II (S1–S3: forward fold, detector
  rate-share weights, multiplicative weighted-geometric-mean adjustment),
  STAYSL-class (T1–T4: fold, anchored-damped objective, increment form,
  clip + active-set pin), GRAVEL (G1–G3: fold, chi-square weights,
  adjustment). Derivations live in the crate-level rustdoc (one module doc
  per method); code comments cite equation labels, never external paths.
- One method per cycle: this crate ships SAND-II (McElroy et al.,
  AFWL-TR-67-41, 1967 — US government work, clean-room from the report),
  STAYSL-class (Perey, ORNL/TM-6062, 1977), and GRAVEL (Matzke, PTB-N-19,
  1994), each in its own module with the same iterator + diagnostics
  shape; do not add a second method to `sandii.rs`. MAXED lands in a later
  cycle the same way.
- Convergence contract: per-group relative change `|Δφ|/φ` strictly below
  an explicit tolerance, with an explicit iteration cap. Exhausting the
  cap is `Error::NotConverged` — a hard fail, never a silent partial
  spectrum (the `tritium` face-Newton precedent).
- Degenerate-input policy (all named): zero response row against a nonzero
  rate, or a rate whose support was pinned to zero mid-iteration, is
  `Error::RatesUnreachable`; detectors folding to zero contribute no
  weight; unconstrained groups keep the guess exactly (factor 1). Pinned
  divergences: SAND-II zero measurements pin their groups to zero and only
  the guess must be strictly positive (the update is multiplicative);
  GRAVEL zero measurements carry zero weight and are not fitted; STAYSL-class
  clips non-positive groups to zero and pins them (active set).
- Caller constants only: response matrix, measured rates, guess spectrum,
  and energy-group bounds are inputs. IRDFF and other IAEA-copyright
  libraries are never vendored; the IRDFF-II v1 pack ships as a runtime
  hash-pinned download (`nucleide.data.fetch_irdff` plus
  `parse_irdff_g725`, owned by `nucleide-nuclei`), parsed into caller-ready
  response rows. The IRDFF-II analytical
  benchmark-field shapes (Trkov et al., Nucl. Data Sheets 163 (2020) 1)
  are published facts and may be re-derived as test inputs with citation;
  the tabulated IAEA group spectra are never used.
- Any least-squares inside this crate routes through the shared
  `nucleide-linalg` `lstsq` kernel (the STAYSL-class consumer); the crate
  declares its `linalg` workspace edge up front.
- Synthetic fixtures only: hand-built response matrices plus closed-form
  spectra with recorded provenance; never evaluated-library data.
- Out of scope (do not expand here): uncertainty propagation beyond the
  landed UQ-lite conventions, outlier-foil rejection (the original code's
  σ-discard feature), cover-material corrections, MAXED maximum entropy.

## Work Guidance

- New methods add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate depends on `nucleide-linalg` only among
  workspace crates; bindings depend on it, never the reverse.
- No external unfolding oracle exists (PyNE/OpenMC ship no SAND-II/STAYSL/
  GRAVEL counterpart), so the validation oracle is analytic-only: synthetic
  round-trips plus IRDFF-II analytical shapes, always run, nothing can
  SKIP.

## Verification

- `cargo test -p nucleide-unfold` (forward-fold identity, fixed point,
  determined/underdetermined round-trips, IRDFF-II shape recovery,
  contract gates).
- `pytest tests/test_unfold.py` after `maturin develop`.
- `validation/unfold_vs_analytic.py` runs inside `run_all.sh`
  (auto-discovered; report name `unfold` is registered in
  `render_results.py` `SECTION_ORDER`).

## Child NAD Index

None.
