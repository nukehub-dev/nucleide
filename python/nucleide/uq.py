"""Seeded UQ-lite sampling kernel and decay perturber (backed by the `nucleide-linalg` crate)."""

from typing import Any

from nucleide._internal import (
    uq_check_convergence,
    uq_lognormal_cov,
    uq_lognormal_mean,
    uq_passthrough,
    uq_perturb_branches,
    uq_perturb_energies,
    uq_perturb_fission_yields,
    uq_sample_cov,
    uq_sample_lognormal,
    uq_sample_mean,
    uq_sample_mvn,
)

__all__ = [
    "sample_mvn",
    "sample_lognormal",
    "lognormal_mean",
    "lognormal_cov",
    "sample_mean",
    "sample_cov",
    "check_convergence",
    "perturb_branches",
    "perturb_energies",
    "passthrough",
    "perturb_fission_yields",
]


def sample_mvn(mean: list[float], cov: list[list[float]], n: int, seed: int) -> dict[str, Any]:
    """Draw ``n`` samples ``x ~ N(mean, cov)`` reproducibly from ``seed``.

    ``cov`` is caller-supplied (no vendored stores). Returns ``samples``
    (list of ``n`` row lists), ``method`` (``"cholesky"`` or
    ``"eigen_clip"``), and the unclipped ``min_eigen``/``max_eigen``
    (``None`` on the Cholesky path). Identical inputs always yield
    identical samples.
    """
    return uq_sample_mvn(mean, cov, n, seed)


def sample_lognormal(
    mean_log: list[float], cov: list[list[float]], n: int, seed: int
) -> dict[str, Any]:
    """Draw ``n`` log-normal samples ``y = exp(x)``, ``x ~ N(mean_log, cov)``.

    ``mean_log``/``cov`` are log-space MVN parameters (never the moments of
    ``y``); the closed-form moments are ``E[y_i] = exp(mu_i + C_ii/2)`` (see
    :func:`lognormal_mean` / :func:`lognormal_cov`). Same return shape as
    :func:`sample_mvn`; identical inputs always yield identical samples.
    """
    return uq_sample_lognormal(mean_log, cov, n, seed)


def lognormal_mean(mean_log: list[float], cov: list[list[float]]) -> list[float]:
    """Closed-form log-normal mean ``E[y_i] = exp(mu_i + C_ii/2)``."""
    return uq_lognormal_mean(mean_log, cov)


def lognormal_cov(mean_log: list[float], cov: list[list[float]]) -> list[list[float]]:
    """Closed-form log-normal covariance."""
    return uq_lognormal_cov(mean_log, cov)


def sample_mean(samples: list[list[float]]) -> list[float]:
    """Sample mean over draws (one entry per dimension)."""
    return uq_sample_mean(samples)


def sample_cov(samples: list[list[float]]) -> list[list[float]]:
    """Unbiased sample covariance (``1/(n-1)``, matching SANDY ``get_cov``)."""
    return uq_sample_cov(samples)


def check_convergence(
    mean: list[float],
    cov: list[list[float]],
    samples: list[list[float]],
    mean_tol: float,
    cov_tol: float,
) -> dict[str, Any]:
    """Sample mean/covariance diagnostics: ``mean_err_max``, ``cov_err_fro``,
    the echoed tolerances, and ``passed``."""
    return uq_check_convergence(mean, cov, samples, mean_tol, cov_tol)


def perturb_branches(base: list[float], rel: list[float]) -> list[float]:
    """Perturb kept branch fractions with relative deltas, preserving the
    incoming ``1 - BR(SF)`` deficit by renormalisation."""
    return uq_perturb_branches(base, rel)


def perturb_energies(
    base: list[float], delta: list[float], convention: str = "relative"
) -> list[float]:
    """Perturb decay energies under ``convention`` (``"relative"``,
    ``"absolute"``, or ``"lognormal"`` (``nominal * exp(delta)``); negative
    results clamp to zero (a no-op for lognormal draws)."""
    return uq_perturb_energies(base, delta, convention)


def passthrough(delta: list[float]) -> list[float]:
    """Finiteness-checked copy of a perturbation vector."""
    return uq_passthrough(delta)


def perturb_fission_yields(base: list[float], rel: list[float]) -> list[float]:
    """Fission-yield perturbation (named-open: waits on ENDF fission-yield tapes).

    Always raises; fission-yield blocks stay out of the decay-only sub-scope.
    """
    return uq_perturb_fission_yields(base, rel)
