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
- **Depletion (CASL chain)** — full 228-nuclide CASL/VERA simplified chain
  (fission yields, decay branching) on a fresh-UO2 inventory, vs OpenMC
  0.16.0. The chain file is downloaded once to the git-ignored
  `validation/.cache/` (see `validation/README.md` for provenance).
- **Enrichment cascades** — default uranium and von-Halle tungsten feeds, vs
  PyNE 0.7.5 `multicomponent` (numeric and symbolic solvers).
- **MAGIC weight windows** — total-mode and per-group MAGIC on a synthetic MCNP
  meshtal, vs PyNE's documented formula, plus always-run OpenMC/Serpent
  weight-window emission re-parse probes (the container also loads the
  emitted fragment through OpenMC's own XML readers; Serpent load is a loud
  SKIP, proprietary and not installed).
- **Point kinetics** — prescribed-reactivity transients vs analytic gates
  (prompt-jump, 1-group closed form, inhour stable period) plus a ramp
  cross-check against the upstream PyRK neutronics block.
- **Spectroscopy** — smoothing, background/gross/net counting, and energy
  bins on a synthetic 7-channel vector vs hand values, vs PyNE
  `spectanalysis`/`gammaspec`.
- **Nuclear data** — atomic masses, natural abundances, half-lives, and
  name-dialect conversions, vs PyNE and OpenMC.
- **Parsers** — Serpent `res`/`dep`/`det` readers vs serpentTools, and MCNP
  `xsdir`/surface-source/PTRAC readers vs PyNE, on the committed fixtures
  (PyMOAB-dependent oracles skip with an explicit note).
- **Activation workflow arithmetic** — flat-energy spectral regrouping against
  independent hand values and PyNE's public overlap API when available;
  volume-weighted source splitting against per-zone conservation vectors.
  These checks do not establish activation-rate or transport-dose equivalence.
- **Timings** — mean wall-time over 20 repeats for the depletion, enrichment,
  and MAGIC solves, plus native Rust Criterion figures.

`validation/make_figures.py` also renders the paper figures
(`validation/figures/timings.png`, `depletion_agreement.png`,
`kinetics_transient.png`, and `spectroscopy_overlay.png`) from the
machine-readable results; the figures are committed but generated — never
hand-edit them.

The canonical environment is the container built from
`validation/Containerfile` (Python 3.12, PyNE 0.7.5 from conda-forge, OpenMC
0.16.0 built from the upstream release tag, and the Nucleide abi3 release wheel
built by the PyO3 maturin container).

## Reproduce the results

One command from the repository root:

```bash
./validation/run_container.sh
# Select Docker explicitly when desired:
CONTAINER_ENGINE=docker ./validation/run_container.sh
```

The script prefers installed Podman and otherwise uses Docker. An explicit
engine selection fails clearly if unavailable. Both run the same environment
and comparisons. Runtime jobs use 6 GiB memory, no swap, 128 tasks and two
CPUs; OpenMC and Rust builds use two jobs. Numerical library threads are
limited to two so timing results reflect that execution budget.

This builds the image once, rebuilds the Nucleide wheel, runs every comparison,
and regenerates `validation/results.md` from the machine-readable JSON reports
in `validation/results/`.

## Full measured tables

See the committed results for the complete correctness and timing tables:

<https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md>
