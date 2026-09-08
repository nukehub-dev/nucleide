# `crates/emit` AGENTS.md

## Purpose

Single-material emission to legacy transport-code cards: one `Material`
renders to MCNP/Serpent/FLUKA/ALARA/PARTISN dialects plus a mass-drift
report. Pure glue over the `*-io` crates; no new physics.

## Ownership

Owns `crates/emit/src/` (one module per dialect plus shared
`EmitOptions`/`DriftTable` in `lib.rs`), the Python surface
(`nucleide.emit`, `emit_cards`/`emit_drift_table` in `_internal`),
`tests/test_emit.py`, and the `emit_vs_self.py` validation oracle.

## Local Contracts

- Dialect conventions (deviations are bugs, not style):
  - MCNP: `m<number>` card, `{zaid}.{xs_suffix}` ids (default `80c`),
    negative mass fractions, greedy packing within the 128-column cap
    (LA-UR-18-20808 §2.6.2), 5-blank continuation indent.
  - Serpent: `mat <name> <-density>` with `{zaid}.{lib}` ids (default
    `03c`) and negative per-nuclide mass fractions, one nuclide per line
    (no continuation markers needed).
  - FLUKA: whole material as one `COMPOUND` card with mass fractions via
    existing `compound_str`; nuclides without a table name/mass are
    `dropped`, never fatal.
  - ALARA: `mixture <name>` block with `element <name> 1.0 <massfrac>`
    entries; text matches `AlaraDeck` display formatting.
  - PARTISN: single-zone deck, uppercase canonical ISOTXS labels; `validate()`
    is never called (no `IsotxsLib` available at emission time).
- Drift semantics: `rel_drift = (mass_in - mass_out) / mass_in` from
  `accounted` masses; `reparsed` is true only after a successful reader
  round-trip (MCNP, ALARA). Serpent/FLUKA/PARTISN have no material readers
  in this workspace — their drift is analytic by design.
- Density rule: Serpent/FLUKA/PARTISN require mass density
  (`EmitOptions::density` overrides `Material::density`, else
  `MissingDensity`); MCNP/ALARA cards carry none.
- Synthetic fixtures only; golden strings live in unit tests, not files.

## Work Guidance

- New dialects add a module plus a `Code` variant and keep
  `emit_all`/`Code::all` order in sync; Python/tests/docs follow in the
  same change (bindings thin, no logic).
- Keep the layering: this crate depends on `material` + `*-io` crates;
  never the reverse, never on bindings.

## Verification

- `cargo test -p nucleide-emit` (includes MCNP/ALARA re-parse round-trips).
- `pytest tests/test_emit.py` after `maturin develop`.
- `validation/emit_vs_self.py` runs inside `run_all.sh` (self-consistency
  oracle; no external code offers this comparison).

## Child NAD Index

None.
