# `crates/blanket` AGENTS.md

## Purpose

TBR and blanket power bookkeeping over caller transport tallies: raw TBR
from caller `(bred, source)` pairs, multiplicative per-port coverage
haircuts, blanket energy multiplication, tritium burn rate, and the
breeding-margin / net-surplus fuel-cycle metrics. Pure arithmetic —
no transport solving, no geometry optimization, no evaluated data.

## Ownership

Owns `crates/blanket/src/` (`lib.rs` kernel + gates, `error.rs`), the
Python surface (`nucleide.blanket`, `blanket_*` in `_internal`), and
`tests/test_blanket.py`. No fixtures, no validation oracle script.

## Local Contracts

- Equation set (B1–B4, P1–P2, F1–F2) is pinned in the crate-level rustdoc:
  `(B1)` raw TBR ratio, `(B2)` multiplicative port haircut
  `raw * Prod_i (1 - f_i)`, `(B3)` margin `eff - 1`, `(B4)` requirement
  check, `(P1)` multiplication ratio, `(P2)` blanket power product,
  `(F1)` burn `P * K` with the pinned `17.6` MeV / `3.0160492` u
  constants, `(F2)` net surplus `(eff - 1) * burn`. Comments cite
  equation labels, never external paths.
- Stellaris Point-A values are regression gates, not inputs: raw TBR
  `1.1070` at exact equality, `1.07379` penalized (report print `1.074`),
  multiplication `1.20`, burn `416.6` g/day at the implied `2714.82` MW
  through `(F1)`. Sub-breakeven results stay negative (deficit reported,
  never clamped); malformed inputs fail loudly through the named error
  set.
- Synthetic hand vectors only: the only "data" in this crate is the
  pinned `(F1)` constant pair plus the SI exact conversions, all spelled
  as `pub const` with derivations in rustdoc. Never vendor nuclear data.
- Out of scope (do not expand here): transport solving, TBR target
  solving / geometry optimization, any coupling into `tritium` (that
  crate keeps its no-breeding-coupling scope line — this companion owns
  the breeding-side arithmetic), HDF5 paths, evaluated-data consumption.

## Work Guidance

- New metrics add a kernel fn plus Python facade/tests in the same change
  (bindings thin, no logic).
- Keep the layering: this crate has no workspace dependencies; bindings
  depend on it, never the reverse.
- No `validation/*_vs_*.py` script: no external oracle exists for
  bookkeeping arithmetic (analytic-gate replay only, per the crate tests
  plus `tests/test_blanket.py` with its independent burn-constant
  transcription).

## Verification

- `cargo test -p nucleide-blanket` (Stellaris hand-vector gates,
  constant identity, loud-input gates).
- `pytest tests/test_blanket.py` after `maturin develop`.

## Child NAD Index

None.
