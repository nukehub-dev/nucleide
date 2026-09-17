"""Neutron spectrum unfolding (backed by the `nucleide-unfold` crate)."""

from typing import Any

from nucleide._internal import (
    unfold_forward_fold,
    unfold_gravel,
    unfold_maxed,
    unfold_sandii,
    unfold_staysl,
)

__all__ = ["sandii", "staysl", "gravel", "maxed", "forward_fold"]


def sandii(
    response: list[list[float]],
    rates: list[float],
    guess: list[float],
    *,
    tolerance: float = 1e-3,
    max_iterations: int = 200,
) -> dict[str, Any]:
    """SAND-II iterative spectral adjustment (McElroy et al., AFWL-TR-67-41, 1967).

    ``response`` holds one row per detector/reaction (all rows one value per
    energy group), ``rates`` the measured rate per detector, and ``guess``
    one strictly positive value per energy group. The guess is adjusted
    until folding it through the response matrix reproduces the rates.
    ``tolerance`` is the largest per-group relative change between
    successive adjustments the run converges under; ``max_iterations`` is
    the explicit adjustment cap — exhausting it raises ``ValueError``
    (non-convergence is a hard fail, never a silent partial spectrum).
    Returns ``spectrum``/``rates``/``rate_factors``/``iterations``/
    ``tolerance``/``max_rel_change``.
    """
    return unfold_sandii(response, rates, guess, tolerance, max_iterations)


def staysl(
    response: list[list[float]],
    rates: list[float],
    sigmas: list[float],
    guess: list[float],
    *,
    tolerance: float = 1e-3,
    max_iterations: int = 200,
    damping: float = 1e-3,
) -> dict[str, Any]:
    """STAYSL-class damped least-squares adjustment (Perey, ORNL/TM-6062, 1977).

    ``response`` holds one row per detector/reaction (all rows one value per
    energy group), ``rates`` the measured rate per detector, ``sigmas`` one
    strictly positive measurement sigma per detector (the weight is
    ``1/σ²``), and ``guess`` one strictly positive value per energy group —
    the prior the damping pulls toward. The least-squares fit of the rates
    is computed with Tikhonov damping ``damping`` (``λ``): directions the
    data constrains much stronger than ``λ`` are fitted to the rates,
    weaker ones keep the guess. Groups the solve drives non-positive are
    pinned to zero and the system is re-solved on the reduced set.
    ``tolerance`` is the largest per-group relative change between
    successive solves the run converges under; ``max_iterations`` is the
    explicit cycle cap — exhausting it raises ``ValueError``
    (non-convergence is a hard fail, never a silent partial spectrum).
    Unlike :func:`sandii` the update is additive, so the rates are
    reproduced up to the damping slack rather than exactly. Returns
    ``spectrum``/``rates``/``rate_factors``/``iterations``/``tolerance``/
    ``max_rel_change``.
    """
    return unfold_staysl(response, rates, sigmas, guess, tolerance, max_iterations, damping)


def gravel(
    response: list[list[float]],
    rates: list[float],
    sigmas: list[float],
    guess: list[float],
    *,
    tolerance: float = 1e-3,
    max_iterations: int = 200,
) -> dict[str, Any]:
    """GRAVEL chi-square-weighted adjustment (Matzke, PTB-N-19, 1994).

    ``response`` holds one row per detector/reaction (all rows one value per
    energy group), ``rates`` the measured rate per detector, ``sigmas`` one
    strictly positive measurement sigma per detector (the per-detector
    weight factor is ``N_i²/σ_i²``, so a precisely measured rate pulls
    harder than a sloppy one), and ``guess`` one strictly positive value
    per energy group. The guess is adjusted by the same
    positivity-preserving weighted geometric mean as :func:`sandii` until
    the fold reproduces the rates. ``tolerance`` is the largest per-group
    relative change between successive adjustments the run converges under;
    ``max_iterations`` is the explicit adjustment cap — exhausting it raises
    ``ValueError`` (non-convergence is a hard fail, never a silent partial
    spectrum). Unlike :func:`sandii`, zero measurements carry zero weight
    and are simply not fitted — they pin nothing to zero. Returns
    ``spectrum``/``rates``/``rate_factors``/``iterations``/``tolerance``/
    ``max_rel_change``.
    """
    return unfold_gravel(response, rates, sigmas, guess, tolerance, max_iterations)


def maxed(
    response: list[list[float]],
    rates: list[float],
    sigmas: list[float],
    guess: list[float],
    *,
    target_chi2: float | None = None,
    tolerance: float = 1e-3,
    max_iterations: int = 200,
) -> dict[str, Any]:
    """MAXED maximum-entropy adjustment (Reginatto & Goldhagen, Health Phys. 77 (1999) 579).

    ``response`` holds one row per detector/reaction (all rows one value per
    energy group), ``rates`` the measured rate per detector, ``sigmas`` one
    strictly positive measurement sigma per detector (the chi-square weight
    is ``1/σ²``), and ``guess`` one strictly positive value per energy group —
    the default model the relative entropy is measured against. The guess is
    adjusted through the exponential Lagrange family until the fold
    reproduces the rates. ``target_chi2`` is the chi-square the run accepts
    at or below (``None`` selects the detector count, the chi-square
    expectation — pass an explicit smaller target down to 0.0 to demand a
    tighter fit); ``tolerance`` is the largest per-group relative change
    between successive adjustments the run converges under (jointly with the
    target); ``max_iterations`` is the explicit adjustment cap — exhausting
    it, or converging in relative change to a fit whose chi-square still
    exceeds the target, raises ``ValueError`` (non-convergence is a hard
    fail, never a silent partial spectrum). Unlike :func:`sandii`, zero
    measurements carry zero weight and are simply not fitted — they pin
    nothing to zero. Returns ``spectrum``/``rates``/``rate_factors``/
    ``iterations``/``tolerance``/``max_rel_change``.
    """
    return unfold_maxed(response, rates, sigmas, guess, target_chi2, tolerance, max_iterations)


def forward_fold(response: list[list[float]], spectrum: list[float]) -> list[float]:
    """Fold a spectrum through a response matrix (one rate per detector row)."""
    return unfold_forward_fold(response, spectrum)
