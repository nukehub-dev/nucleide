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
- `mcnp/xsdir` — xsdir parser oracle; `dummy_xsdir` is staged for the website
- `mcnp/meshtal` — single + multiple meshtal files; `mcnp_meshtal_single_meshtal.txt`
  is staged for the website
- `mcnp/ssw` — mcnp5/mcnp6/mcnpx surfsrc + one-track (SSW round-trip oracles)
- `mcnp/ptrac` — i4/i8 + MCNP6 variants + input card
- `mcnp/wwinp` — n/p/np weight-window files
- `mcnp/inp` — material-bearing input decks (+ commented variant)
- `mcnp/mctal/synthetic_*` — generated kcode decks (no public corpus exists)
- `depletion/` — depletion-chain XML files (simple chains, Ni chain);
  `chain_simple.xml` is staged for the website
- `serpent/` — Serpent 1 & 2 res/dep/det outputs
- `fluka/` — USRBIN `.lis` files (single/multiple/degenerate) + test input
- `alara/` — UW ALARA samples vendored verbatim (terms in `alara/LICENSE.ALARA`):
  input-deck oracles `decks/sample2` + `decks/sample3`; group-flux files
  `flux/fluxin2` + all-zero edge case `flux/fluxin_zeros`; support libraries
  `libs/sampleMatlib` (material lib), `libs/myElelib` (element lib),
  `libs/NRCA` + `libs/NRCC` (WDR libs); activation-output oracle
  `output/sample2.out`
- `cccc/` — synthetic ISOTXS/RTFLUX samples authored for Nucleide (no license
  needed): `isotxs_sample` (2 nuclides x 3 groups), `rtflux_sample`
  (2 points x 3 groups)
- `fispact/` — synthetic FISPACT-II-style inventory authored for Nucleide (no
  license needed): `inventory.fis` (3 cooling steps x 4 rows x 3 variables)
- `origen/` — synthetic ORIGEN TAPE samples authored for Nucleide (no license
  needed): `tape5_sample` (input echo), `tape6_sample` (inventory),
  `tape9_sample` (decay constants)

Still to add in later phases: truncated PTRAC samples for fuzzing; recorded
CRAM input/output pairs for regression pinning (currently validated
analytically).
