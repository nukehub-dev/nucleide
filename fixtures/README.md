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

Runtime consumption: the interactive website tutorials load a subset of these
fixtures in the browser. `npm run sync-data` (from `website/`) stages them
into the git-ignored `website/public/data/`; see `website/AGENTS.md` for the
contract. Currently staged: the MCNP meshtal/xsdir samples, the
`deck_minimal`/`deck_l3` decks, the simple depletion chain, the Serpent
`sample_res.m`/`sample2_dep.m`/`sample_det.m` samples, the ALARA
`output/sample2.out` and FISPACT `inventory.fis` activation outputs, and the
FLUKA `fluka_usrbin_single.lis` USRBIN sample.

<!-- GEN:fixture-index:START -->

| Area | Files | Size | Contents |
| --- | --- | --- | --- |
| `alara/` | 10 | 97,611 bytes | alara — UW ALARA samples vendored verbatim (terms in `alara/LICENSE.ALARA`) |
| `cccc/` | 2 | 285 bytes | cccc — Synthetic ISOTXS/RTFLUX samples authored for Nucleide (no license needed) |
| `data/` | 2 | 7,859,273 bytes | data — DOE/PNNL Materials Compendium JSON (+ its license) |
| `depletion/` | 4 | 37,843 bytes | depletion — Depletion-chain XML files (simple chains, Ni chain) |
| `endl/` | 1 | 762 bytes | endl — Synthetic EEDL-style tables authored for Nucleide (no license needed) |
| `fispact/` | 1 | 1,058 bytes | fispact — Synthetic FISPACT-II-style inventory authored for Nucleide (no license needed) |
| `fluka/` | 4 | 7,685 bytes | fluka — USRBIN `.lis` files (single/multiple/degenerate) + test input |
| `kinetics/` | 3 | 2,493 bytes | fixtures/kinetics/ — synthetic point-kinetics inputs and oracles |
| `mcnp/` | 27 | 471,097 bytes | mcnp — MCNP input decks, tallies, and binary-format oracles |
| `mcpl/` | 6 | 4,981 bytes | mcpl — Synthetic MCPL interchange fixtures authored for Nucleide (no license needed) |
| `origen/` | 3 | 681 bytes | origen — Synthetic ORIGEN TAPE samples authored for Nucleide (no license needed) |
| `serpent/` | 6 | 653,120 bytes | serpent — Serpent 1 & 2 res/dep/det outputs |
| `spectroscopy/` | 6 | 6,380 bytes | fixtures/spectroscopy/ — synthetic gamma-spectroscopy inputs and oracles |
| `uq/` | 4 | 1,127 bytes | UQ-lite synthetic oracle inputs |

<!-- GEN:fixture-index:END -->

Still to add in later phases: truncated PTRAC samples for fuzzing; recorded
CRAM input/output pairs for regression pinning (currently validated
analytically).
