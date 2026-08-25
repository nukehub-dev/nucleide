# fixtures/ — golden test data

Parsers are validated against byte-exact golden files; any change that alters
parsed output must update a fixture deliberately, never accidentally.

1. Golden files are vendored verbatim under format-specific directories.
2. Parser output is pinned by unit-test assertions and snapshots; CI fails on
   any drift.
3. Numeric kernels (CRAM) are checked against analytic vectors on shared
   chain inputs.

Contents:

- `data/` — DOE/PNNL Materials Compendium JSON (+ its license)
- `mcnp/xsdir` — xsdir parser oracle
- `mcnp/meshtal` — single + multiple meshtal files
- `mcnp/ssw` — mcnp5/mcnp6/mcnpx surfsrc + one-track (SSW round-trip oracles)
- `mcnp/ptrac` — i4/i8 + MCNP6 variants + input card
- `mcnp/wwinp` — n/p/np weight-window files
- `mcnp/inp` — material-bearing input decks (+ commented variant)
- `mcnp/mctal/synthetic_*` — generated kcode decks (no public corpus exists)
- `depletion/` — depletion-chain XML files (simple chains, Ni chain)
- `serpent/` — Serpent 1 & 2 res/dep/det outputs
- `fluka/` — USRBIN `.lis` files (single/multiple/degenerate) + test input

Still to add in later phases: truncated PTRAC samples for fuzzing; recorded
CRAM input/output pairs for regression pinning (currently validated
analytically).
