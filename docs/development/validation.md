---
title: Cross-code Validation
sidebar:
  order: 2
---

The `validation/` harness compares Nucleide against the reference codes PyNE
0.7.5 and OpenMC 0.16.0 on shared inputs, plus coarse timing comparisons. It
produces the measured tables quoted by the JOSS paper.

Comparisons cover:

- **Depletion** — CRAM-48 solve on `fixtures/depletion/chain_ni.xml` and a
  3-nuclide analytic Bateman check, vs OpenMC 0.16.0.
- **Enrichment cascades** — default uranium and von-Halle tungsten feeds, vs
  PyNE 0.7.5 `multicomponent` (numeric and symbolic solvers).
- **MAGIC weight windows** — total-mode and per-group MAGIC on a synthetic MCNP
  meshtal, vs PyNE's documented formula.
- **Nuclear data** — atomic masses, natural abundances, half-lives, and
  name-dialect conversions, vs PyNE and OpenMC.
- **Timings** — mean wall-time over 20 repeats for the depletion, enrichment,
  and MAGIC solves, plus native Rust Criterion figures.

The canonical environment is the container built from
`validation/Containerfile` (Python 3.12, PyNE 0.7.5 from conda-forge, OpenMC
0.16.0 built from the upstream release tag, and the Nucleide abi3 release wheel
built by the PyO3 maturin container).

## Reproduce the results

One command from the repository root:

```bash
./validation/run_container.sh
```

This builds the image once, rebuilds the Nucleide wheel, runs every comparison,
and regenerates `validation/results.md` from the machine-readable JSON reports
in `validation/results/`.

## Full measured tables

See the committed results for the complete correctness and timing tables:

<https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md>
