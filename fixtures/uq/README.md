# UQ-lite synthetic oracle inputs

Hand-built sampling-kernel inputs for the UQ-lite covariance-recovery gates
(Rust `nucleide-linalg` tests, `tests/test_uq.py`, and
`validation/uq_lite_vs_sandy.py`). All numbers are synthetic round values
chosen for closed-form checks — no evaluated nuclear data, no covariance
library content, nothing read from ENDF tapes.

- `cov_2x2.json`: 2-D block, variances 0.25/0.16 with covariance 0.10
  (correlation 0.5), nonzero mean, pinned seed, `n`, and the gate width `k`
  in standard errors.
- `cov_3x3.json`: 3-D block with mixed-sign correlations (positive-definite
  by leading principal minors), pinned seed, `n`, and `k`.
- `decay_perturb.json`: synthetic kept-branch vector summing to 0.90 (the
  0.10 deficit stands in for a dropped SF branch per the `1 − BR(SF)` store
  convention) plus a relative perturbation, and a synthetic two-entry energy
  vector with its delta and convention.
- `lognormal_2x2.json`: synthetic log-space 2-D block (`mean_log`, `cov`)
  for log-normal sampling (`y = exp(x)`, `x ~ N(mean_log, cov)`), with
  pinned seed, `n`, and `k`. Gates check `ln(y)` recovers the log-space
  block within `k` MVN standard errors plus the sample mean of `y` against
  the closed-form `E[y_i] = exp(mu_i + C_ii/2)`.

Schema per covariance file: `{"mean": [...], "cov": [[...]],
"seed": u64, "n": draws, "k": gate width in standard errors, "note": str}`.
(`lognormal_2x2.json` uses `"mean_log"` in place of `"mean"` for the
log-space MVN parameters; same types.)
Gates assert sample mean/covariance within `k` standard errors of the inputs
(`std(mean_i) = sqrt(C[i][i]/n)`,
`var(S[i][j]) = (C[i][i]*C[j][j] + C[i][j]^2)/(n-1)` for MVN draws).
