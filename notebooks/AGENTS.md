# `notebooks/` AGENTS.md

## Purpose

Runnable end-to-end workflow notebooks (Jupyter `.ipynb`) that double as
adoption material: first-wall damage + UQ, activation screening + Sublet
totals, and tokamak source sampling + card emission. Each notebook runs top
to bottom in CI and reads only committed fixtures with fixed seeds.

## Ownership

Owns `notebooks/*.ipynb`, executed by `tests/test_notebooks.py` (fresh
kernel per notebook, repo root as cwd) and by nothing else.

## Local Contracts

- Synthetic inputs only; committed fixtures by repo-relative paths; fixed
  seeds everywhere. No downloads, no network, no vendored data, no machine
  outputs, no writes outside `tempfile` scratch. New notebooks follow the
  same contract or they do not land.
- Committed notebooks carry no outputs and no execution counts
  (`test_notebooks_are_stripped` fails otherwise). Results are recomputed
  by the reader, never vendored in the file — diffs stay reviewable.
- Stdlib + `nucleide` (+ `numpy`, already a test dependency) only. No
  matplotlib, no extra kernel packages: the CI install line in
  `.github/workflows/ci.yml` (`nbclient nbformat ipykernel`) is the full
  notebook dependency set — extend it only with a written reason.
- Colab badges point at `main` and require the PyPI release to be live;
  the first markdown cell of each notebook says so (`%pip install`
  stays a reader-side step, never a committed cell, so CI stays offline).
- Asserts inside notebooks are gates, not decoration: a stale number
  fails `test_notebooks_execute` instead of rotting silently.

## Verification

- `pytest tests/test_notebooks.py -q` after `maturin develop`.
- `pytest tests/test_tutorial_code.py -q` for the markdown tutorials,
  which the same "docs are verified examples" rule covers.

## Child NAD Index

None.
