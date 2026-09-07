# fixtures/ — golden test data

Parsers are validated against byte-exact golden files; any change that alters
parsed output must update a fixture deliberately, never accidentally.

1. Golden files are vendored verbatim under format-specific directories.
2. Parser output is pinned by unit-test assertions and snapshots; CI fails on
   any drift.
3. Numeric kernels (CRAM) are checked against analytic vectors on shared
   chain inputs.

Contents:

Each `fixtures/<area>/README.md` describes its area; the index below is
generated from those files — run `python3 scripts/gen-reference.py --write`
after adding fixtures (CI enforces freshness via `--check`).

<!-- GEN:fixture-index:START -->

| Area | Files | Size | Contents |
| --- | --- | --- | --- |
| `alara/` | 10 | 97,611 bytes | alara — UW ALARA samples vendored verbatim (terms in `alara/LICENSE.ALARA`) |
| `cccc/` | 2 | 285 bytes | cccc — Synthetic ISOTXS/RTFLUX samples authored for Nucleide (no license needed) |
| `data/` | 2 | 7,859,273 bytes | data — DOE/PNNL Materials Compendium JSON (+ its license) |
| `depletion/` | 4 | 37,843 bytes | depletion — Depletion-chain XML files (simple chains, Ni chain) |
| `fispact/` | 1 | 1,058 bytes | fispact — Synthetic FISPACT-II-style inventory authored for Nucleide (no license needed) |
| `fluka/` | 4 | 7,685 bytes | fluka — USRBIN `.lis` files (single/multiple/degenerate) + test input |
| `mcnp/` | 20 | 468,802 bytes | mcnp — MCNP input decks, tallies, and binary-format oracles |
| `origen/` | 3 | 681 bytes | origen — Synthetic ORIGEN TAPE samples authored for Nucleide (no license needed) |
| `serpent/` | 6 | 653,120 bytes | serpent — Serpent 1 & 2 res/dep/det outputs |

<!-- GEN:fixture-index:END -->

Still to add in later phases: truncated PTRAC samples for fuzzing; recorded
CRAM input/output pairs for regression pinning (currently validated
analytically).
