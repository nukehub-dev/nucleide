# Nuclei NAD

## Purpose

`nucleide-nuclei`: canonical nucid representation, element/naming tables,
name dialects, particle and reaction-name registries, and static nuclear
reference data (`src/data/`) — plus two runtime-parsed packs: the EPA FGR 15
external-dosimetry tables (`src/fgr15.rs`) and the IAEA IRDFF-II v1
foil-response pack (`src/irdff.rs`).

## Ownership

This file owns the FGR 15 and IRDFF-II modules' distribution and parsing
contracts.
Nucid conventions, dialects, and the vendored TSV tables stay documented in
`src/lib.rs`, `src/dialects.rs`, and `src/data.rs` (change those docs with
the code).

## Local Contracts

- **Nothing from the EPA FGR 15 file is committed to the repository**: no
  TSV extracts, no fixture copies, no generated tables. The zip is fetched
  at runtime by `python/nucleide/data.py::fetch_fgr15`, pinned to the
  SHA-256 in `python/nucleide/data.py` (`FGR15_SHA256`), and member text is
  passed to `fgr15::parse_table` (no network or HTTP in this crate).
- **Tests use synthetic hand-built tables only** (see `fgr15.rs` tests and
  `tests/test_fgr15.py`), in the exact `Table_4_*.DAT` layout. Never copy
  real EPA content into a test.
- **Parsing is strict and loud**: data rows are `Symbol-A[letter]` +
  exactly six finite floats; element separators, dashed separators, and
  blank lines are structural; BOM/zero-width characters are stripped; rows
  are never silently skipped, duplicates error, and the row total must equal
  the caller's `expected_rows` (1,252 for the published tables). The units
  string is recorded verbatim from the header — never hardcode a units
  wording.
- **Isomer names** follow the repo-wide `mnopqrstuvxyz` letter convention
  (`m` → state 1, `n` → state 2); the FGR 15 spelling (`Ba-137m`,
  `Sb-124n`) is the module's public name dialect.
- **Screening-level only** — not for safety decisions (EPA screening
  context; say so in user-facing docs and docstrings).
- **Nothing from the IAEA IRDFF-II files is committed to the repository**
  (same distribution contract as FGR 15): the `IRDFF-II_g725.zip` is fetched
  at runtime by `python/nucleide/data.py::fetch_irdff`, pinned to the
  SHA-256 in `python/nucleide/data.py` (`IRDFF_SHA256`), and member text is
  passed to `irdff::parse_g725` (no network or HTTP in this crate). The v1
  pack covers the eight named foil reactions (`V1_REACTIONS`, `MF=3`
  sections only) over the SAND-II 725-group structure; rows feed
  `unfold.sandii`/`staysl`/`gravel` unchanged.
- **IRDFF tests use synthetic hand-built sections only** (see `irdff.rs`
  tests), in the exact `MF=3` layout. Never copy real IAEA content into a
  test.

## Work Guidance

- Extend `fgr15.rs` (parser, scenario/age enums, `Fgr15Table`) rather than
  `data.rs` — FGR 15 is runtime-parsed, not a vendored TSV.
- A new EPA revision means: update `FGR15_URL`/`FGR15_SHA256` in
  `python/nucleide/data.py`, bump the row-count expectation if the nuclide
  census changed, and refresh the spots in
  `validation/nuclear_data_vs_refs.py` together, in one change.
- Out of scope for v1: `Nuclide_Coefficients/` 28-tissue tables,
  `Mono_Coefficients/` monoenergetic photons, biokinetic modeling.

## Verification

- `cargo test -p nucleide-nuclei` (synthetic fgr15/irdff tests live
  in-module).
- `ruff format --check`, `ruff check`, and `mypy` clean for the Python
  facade (`python/nucleide/nuclei.py`, `python/nucleide/data.py`) and
  `tests/test_fgr15.py`.
- `python3 scripts/gen-reference.py --check` after touching the facade.

## Child NAD Index

None.
