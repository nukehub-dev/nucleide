# `crates/unfold` AGENTS.md

## Purpose

Neutron spectrum unfolding from activation-type measurements: adjust a
caller-supplied guess spectrum against measured detector/reaction rates
until the forward fold through the caller-supplied response matrix
reproduces them. Inverse-problem iteration, not a transport solve.

## Ownership

Owns `crates/unfold/src/` (`sandii.rs` SAND-II iterator, `staysl.rs`
STAYSL-class iterator, `gravel.rs` GRAVEL iterator, `maxed.rs` MAXED
maximum-entropy iterator, `lib.rs` forward operator + `error.rs`), the
Python surface (`nucleide.unfold`, `unfold_*` in `_internal`),
`tests/test_unfold.py`, and the `validation/unfold_vs_analytic.py` oracle.

## Local Contracts

- Equation sets are pinned: SAND-II (S1–S3: forward fold, detector
  rate-share weights, multiplicative weighted-geometric-mean adjustment),
  STAYSL-class (T1–T4: fold, anchored-damped objective, increment form,
  clip + active-set pin), GRAVEL (G1–G3: fold, chi-square weights,
  adjustment), MAXED (M1–M6: fold, relative entropy, chi-square, exponential
  family, dual curvature, equilibrated damped dual step with trust-region
  cap). Derivations live in the crate-level rustdoc (one module doc
  per method); code comments cite equation labels, never external paths.
- One method per cycle: this crate ships SAND-II (McElroy et al.,
  AFWL-TR-67-41, 1967 — US government work, clean-room from the report),
  STAYSL-class (Perey, ORNL/TM-6062, 1977), GRAVEL (Matzke, PTB-N-19,
  1994), and MAXED (Reginatto & Goldhagen, Health Phys. 77 (1999) 579 +
  NIM A 476 (2002) 242 — journal equations only, the closed UMG package is
  never touched), each in its own module with the same iterator +
  diagnostics shape; do not add a second method to any method module, and
  no fifth method is scoped after MAXED (the recorded order closes here).
- Convergence contract: per-group relative change `|Δφ|/φ` strictly below
  an explicit tolerance, with an explicit iteration cap. Exhausting the
  cap is `Error::NotConverged` — a hard fail, never a silent partial
  spectrum (the `tritium` face-Newton precedent). MAXED additionally takes
  a caller chi-square target (default: the detector count): converging in
  relative change to a fit whose chi-square still exceeds the target is
  `Error::NotConverged` too.
- Degenerate-input policy (all named): zero response row against a nonzero
  rate, or a rate whose support was pinned to zero mid-iteration, is
  `Error::RatesUnreachable`; detectors folding to zero contribute no
  weight; unconstrained groups keep the guess exactly (factor 1). Pinned
  divergences: SAND-II zero measurements pin their groups to zero and only
  the guess must be strictly positive (the update is multiplicative);
  GRAVEL zero measurements carry zero weight and are not fitted; STAYSL-class
  clips non-positive groups to zero and pins them (active set); MAXED zero
  measurements carry zero weight like GRAVEL (an all-zero measurement set
  returns the guess unchanged after one confirming no-op), and only the
  guess must be strictly positive (the update is exponential).
- IRDFF-II analytical benchmark-shape probes use group bounds bracketing
  each shape's support for solve-based methods (STAYSL-class, MAXED): with
  bounds extending decades past a shape's decay the trailing groups sit at
  the 1e-30 positivity floor, and a pinned-to-zero recovered group reads a
  spurious relative error of exactly 1.0 against that floor. Solve-free
  multiplicative methods (SAND-II, GRAVEL) probe the full range.
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
  σ-discard feature), cover-material corrections, a fifth unfolding method.

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
  contract gates — U1-U23 in `validation/unfold_vs_analytic.py` mirror the
  Rust gates).
- `pytest tests/test_unfold.py` after `maturin develop`.
- `validation/unfold_vs_analytic.py` runs inside `run_all.sh`
  (auto-discovered; report name `unfold` is registered in
  `render_results.py` `SECTION_ORDER`).

## Child NAD Index

None.
