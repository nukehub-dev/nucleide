"""Python-side tests for the UQ-lite sampling kernel (covariance recovery + decay)."""

import json
import math
from pathlib import Path
from typing import Any

import pytest

import nucleide.uq as uq

FIX = Path(__file__).parent.parent / "fixtures" / "uq"


def _load(name: str) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((FIX / name).read_text())
    return data


def test_sample_mvn_is_seeded_reproducible() -> None:
    fx = _load("cov_2x2.json")
    a = uq.sample_mvn(fx["mean"], fx["cov"], 64, fx["seed"])
    b = uq.sample_mvn(fx["mean"], fx["cov"], 64, fx["seed"])
    assert a["samples"] == b["samples"]
    assert a["method"] == "cholesky"
    assert a["min_eigen"] is None and a["max_eigen"] is None
    c = uq.sample_mvn(fx["mean"], fx["cov"], 64, fx["seed"] + 1)
    assert c["samples"] != a["samples"]


def _check_recovery(fx: dict[str, Any]) -> None:
    # Honest statistical gate (same derivation as the Rust tests): for n MVN
    # draws, std(mean_i) = sqrt(C[i][i]/n) and
    # var(S[i][j]) = (C[i][i]*C[j][j] + C[i][j]^2)/(n-1); assert within k
    # standard errors at the fixture-pinned seed.
    mean, cov, n, k, seed = fx["mean"], fx["cov"], fx["n"], fx["k"], fx["seed"]
    out = uq.sample_mvn(mean, cov, n, seed)
    assert out["method"] == "cholesky"
    samples = out["samples"]
    assert len(samples) == n and all(len(s) == len(mean) for s in samples)
    assert all(math.isfinite(v) for s in samples for v in s)
    sm = uq.sample_mean(samples)
    sc = uq.sample_cov(samples)
    dim = len(mean)
    for i in range(dim):
        se = math.sqrt(cov[i][i] / n)
        assert abs(sm[i] - mean[i]) <= k * se
        for j in range(dim):
            se_cov = math.sqrt((cov[i][i] * cov[j][j] + cov[i][j] ** 2) / (n - 1))
            assert abs(sc[i][j] - cov[i][j]) <= k * se_cov


def test_cov_2x2_recovers() -> None:
    _check_recovery(_load("cov_2x2.json"))


def test_cov_3x3_recovers() -> None:
    _check_recovery(_load("cov_3x3.json"))


def test_rank_deficient_block_uses_eigen_clip() -> None:
    out = uq.sample_mvn([0.0, 0.0], [[1.0, 1.0], [1.0, 1.0]], 20000, 20260913)
    assert out["method"] == "eigen_clip"
    assert out["min_eigen"] == pytest.approx(0.0, abs=1e-9)
    assert out["max_eigen"] == pytest.approx(2.0, rel=1e-9)


def test_check_convergence_report() -> None:
    fx = _load("cov_2x2.json")
    out = uq.sample_mvn(fx["mean"], fx["cov"], 512, fx["seed"])
    rep = uq.check_convergence(fx["mean"], fx["cov"], out["samples"], 1.0, 10.0)
    assert rep["passed"] is True
    assert rep["mean_tol"] == 1.0 and rep["cov_tol"] == 10.0
    assert rep["mean_err_max"] >= 0.0 and rep["cov_err_fro"] >= 0.0
    strict = uq.check_convergence(fx["mean"], fx["cov"], out["samples"], 1e-12, 1e-12)
    assert strict["passed"] is False


def test_perturb_branches_preserves_deficit() -> None:
    fx = _load("decay_perturb.json")
    out = uq.perturb_branches(fx["base_branches"], fx["rel"])
    assert sum(out) == pytest.approx(sum(fx["base_branches"]), rel=1e-15)
    assert all(v >= 0.0 for v in out)


def test_perturb_energies_and_conventions() -> None:
    fx = _load("decay_perturb.json")
    out = uq.perturb_energies(fx["base_energies"], fx["delta"], fx["convention"])
    assert out == [pytest.approx(0.5 * 1.2), pytest.approx(1.5 * 0.9)]
    assert uq.perturb_energies([1.0], [0.5], "absolute") == [pytest.approx(1.5)]
    assert uq.perturb_energies([1.0], [-2.0], "absolute") == [0.0]
    # Log-normal: nominal * exp(delta); exp(0) is the identity, stays positive.
    assert uq.perturb_energies([2.0], [0.0], "lognormal") == [pytest.approx(2.0)]
    assert uq.perturb_energies([2.0], [math.log(2.0)], "log_normal") == [pytest.approx(4.0)]
    assert uq.perturb_energies([2.0], [math.log(2.0)], "log-normal") == [pytest.approx(4.0)]
    with pytest.raises(ValueError, match="unknown perturbation convention"):
        uq.perturb_energies([1.0], [0.1], "lhs")


def test_sample_lognormal_reuses_mvn_stream() -> None:
    fx = _load("lognormal_2x2.json")
    mean_log, cov, seed = fx["mean_log"], fx["cov"], fx["seed"]
    mvn = uq.sample_mvn(mean_log, cov, 64, seed)
    logn = uq.sample_lognormal(mean_log, cov, 64, seed)
    assert logn["method"] == mvn["method"] == "cholesky"
    for x, y in zip(mvn["samples"], logn["samples"], strict=True):
        for a, b in zip(x, y, strict=True):
            assert b == pytest.approx(math.exp(a))
            assert math.isfinite(b) and b > 0.0
    again = uq.sample_lognormal(mean_log, cov, 64, seed)
    assert again["samples"] == logn["samples"]
    other = uq.sample_lognormal(mean_log, cov, 64, seed + 1)
    assert other["samples"] != logn["samples"]


def test_lognormal_fixture_moments_match_closed_form() -> None:
    fx = _load("lognormal_2x2.json")
    mean_log, cov, n, k, seed = fx["mean_log"], fx["cov"], fx["n"], fx["k"], fx["seed"]
    out = uq.sample_lognormal(mean_log, cov, n, seed)
    assert out["method"] == "cholesky"
    samples = out["samples"]
    assert len(samples) == n
    assert all(v > 0.0 for s in samples for v in s)
    # L1: ln(y) recovers the log-space block within the MVN k-SE gate.
    ln = [[math.log(v) for v in s] for s in samples]
    sm = uq.sample_mean(ln)
    sc = uq.sample_cov(ln)
    dim = len(mean_log)
    for i in range(dim):
        se = math.sqrt(cov[i][i] / n)
        assert abs(sm[i] - mean_log[i]) <= k * se
        for j in range(dim):
            se_cov = math.sqrt((cov[i][i] * cov[j][j] + cov[i][j] ** 2) / (n - 1))
            assert abs(sc[i][j] - cov[i][j]) <= k * se_cov
    # L2: sample mean of y against the closed-form E[y_i].
    expected = uq.lognormal_mean(mean_log, cov)
    var_y = uq.lognormal_cov(mean_log, cov)
    assert expected == [
        pytest.approx(math.exp(mean_log[0] + cov[0][0] / 2)),
        pytest.approx(math.exp(mean_log[1] + cov[1][1] / 2)),
    ]
    sm_y = uq.sample_mean(samples)
    for i in range(dim):
        se = math.sqrt(var_y[i][i] / n)
        assert abs(sm_y[i] - expected[i]) <= k * se


def test_passthrough_and_fy_hook() -> None:
    assert uq.passthrough([0.1, -0.2]) == [pytest.approx(0.1), pytest.approx(-0.2)]
    with pytest.raises(ValueError, match="non-finite"):
        uq.passthrough([float("inf")])
    with pytest.raises(ValueError, match="named-open"):
        uq.perturb_fission_yields([0.5], [0.1])


def test_malformed_inputs_raise() -> None:
    with pytest.raises(ValueError, match="empty"):
        uq.sample_mvn([], [], 4, 1)
    with pytest.raises(ValueError, match="zero draws"):
        uq.sample_mvn([0.0], [[1.0]], 0, 1)
    with pytest.raises(ValueError, match="dimension mismatch"):
        uq.sample_mvn([0.0, 0.0], [[1.0]], 4, 1)
    with pytest.raises(ValueError, match="asymmetric"):
        uq.sample_mvn([0.0, 0.0], [[1.0, 0.0], [0.5, 1.0]], 4, 1)
    with pytest.raises(ValueError, match="no positive eigenvalue"):
        uq.sample_mvn([0.0, 0.0], [[-1.0, 0.0], [0.0, -2.0]], 4, 1)
    with pytest.raises(ValueError, match="length mismatch"):
        uq.perturb_branches([0.5], [0.1, 0.2])
