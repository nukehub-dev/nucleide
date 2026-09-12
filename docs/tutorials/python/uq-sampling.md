---
title: UQ Sampling
sidebar:
  order: 13
---

Nucleide UQ-lite draws seeded multivariate-normal samples over
caller-supplied covariance blocks — the SANDY role (seeded MVN draws over
caller blocks) with no vendored covariance stores. Engine note: SANDY
itself factorises with SVD behind NumPy's PCG64 stream, while this kernel
uses Cholesky-with-eigen-clip behind ChaCha8 — same target distribution,
so draws are not interchangeable, but the moment estimators agree exactly
(`sample_mean`/`sample_cov` match `Samples.get_mean`/`get_cov` at 1e-9 in
the Tier-2 oracle gate). This tutorial covers the Python API; the
kernel lives in the `linalg` crate (`sample` + `decay` modules).

## Sample a covariance block

`sample_mvn` draws `n` samples `x ~ N(mean, cov)` reproducibly from `seed`:
identical inputs always yield identical samples. The block below is the
synthetic 2x2 from `fixtures/uq/cov_2x2.json` (variances 0.25/0.16,
covariance 0.10 — round numbers, no evaluated data):

```python
from nucleide.uq import check_convergence, sample_cov, sample_mean, sample_mvn

mean, cov = [1.0, 2.0], [[0.25, 0.10], [0.10, 0.16]]
out = sample_mvn(mean, cov, 2000, 20260913)
print(out["method"])  # "cholesky" (positive-definite path)
print(sample_mean(out["samples"]))
print(sample_cov(out["samples"]))
```

`method` names the factorisation path: `"cholesky"` for positive-definite
blocks, `"eigen_clip"` when a positive-semidefinite but singular block
falls back to eigen-clipping (then `min_eigen`/`max_eigen` carry the
unclipped extremes). It is reported, never silent.

## Check convergence

`check_convergence` recomputes the sample mean and unbiased sample
covariance (`1/(n-1)`, SANDY `get_cov` style) and compares them against the
inputs that generated the draws. Tolerances are caller-supplied:

```python
rep = check_convergence(mean, cov, out["samples"], 0.05, 0.05)
print(rep["passed"], rep["mean_err_max"], rep["cov_err_fro"])
```

## Perturb decay data

`perturb_branches` applies relative deltas to one parent's kept branch
fractions and renormalises to preserve the incoming `1 - BR(SF)` deficit —
the evaluated store drops spontaneous-fission branches, so the synthetic
`fixtures/uq/decay_perturb.json` vector sums to 0.90 and the output sums to
0.90 too. `perturb_energies` perturbs per-nuclide energies with no sum
constraint (`"relative"` or `"absolute"`, negatives clamped to zero):

```python
from nucleide.uq import perturb_branches, perturb_energies

print(perturb_branches([0.5, 0.3, 0.1], [0.1, -0.2, 0.0]))
print(perturb_energies([0.5, 1.5], [0.2, -0.1], "relative"))
```

Fission-yield perturbation stays a named-open hook (`perturb_fission_yields`
always raises — it waits on FY tapes), and there are no vendored
covariance, branch, or yield stores anywhere in this path.

## See also

- `crates/linalg/src/sample.rs` and `crates/linalg/src/decay.rs`
  for the Rust API.
- `tests/test_uq.py` for worked examples.
- [Cross-code validation results](https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md)
  for the SANDY moment cross-check (`validation/uq_lite_vs_sandy.py`).
