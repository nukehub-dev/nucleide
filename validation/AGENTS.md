# Validation AGENTS.md

## Purpose

Cross-code validation harness comparing Nucleide against PyNE and OpenMC on
shared inputs, plus coarse timing comparisons. Its committed results are quoted
by the JOSS paper (`paper.md`).

## Ownership

This file owns the validation harness workflow: scripts, the committed
reference tally (`magic_tally.txt`), and the committed results document.

## Local Contracts

- `results.md` and `validation/results/*.json` are **generated artifacts**:
  never hand-edit them. Regenerate by running `run_all.sh` (or
  `run_container.sh`). Cosmetic fixes that do not alter numbers belong in
  `render_results.py` or `common.py`.
- `validation/environment.json` is the **only hand-maintained input** for the
  environment table; keep it up to date when the container or dependency
  channels change.
- Every number quoted in `paper.md` must come from a recorded `run_all.sh` run.
- Scripts must pass `ruff format --check` and `ruff check` (repo
  `pyproject.toml` settings, line-length 100).
- Do not add fixtures here that duplicate `fixtures/` policy; the harness reads
  from `fixtures/` and patches only at runtime, in memory.

## Work Guidance

- Canonical environment is the container built from `Containerfile`
  (`./validation/run_container.sh`, requires podman): Python 3.12,
  PyNE 0.7.5 from conda-forge, OpenMC 0.16.0 built from the upstream tag, and
  the Nucleide abi3 release wheel built by the `ghcr.io/pyo3/maturin`
  container. Version-pin rationale is documented in the `Containerfile`
  header; update it when conda-forge catches up (PyNE > 0.7.5 or
  OpenMC >= 0.16.0).
- A conda env fallback also works:
  (`mamba create -p ~/.conda/envs/nuke-validation -c conda-forge pyne openmc`)
  with the Nucleide release wheel installed into that env
  (`maturin build --release`, then `pip install --force-reinstall target/wheels/nucleide-*.whl`).
- Run: `./validation/run_container.sh`, or `./validation/run_all.sh [python-binary]`
  against any prepared Python.
- PyNE is typically built without PyMOAB; the MAGIC comparison falls back to a
  formula-equivalent pure-Python reference and must say so in `results.md`.

## Verification

- `run_all.sh` exits nonzero on failure; all scripts must complete.
- Repo-wide markdown lint (`npx markdownlint-cli2 "validation/*.md"`) must be
  clean.

## Child NAD Index

None.
