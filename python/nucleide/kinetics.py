"""Prescribed-reactivity point kinetics (backed by the `nucleide-kinetics` crate)."""

from typing import Any

from nucleide._internal import (
    kinetics_equilibrium,
    kinetics_from_ifp,
    kinetics_inhour_rho,
    kinetics_initial_rate,
    kinetics_prompt_jump,
    kinetics_solve,
    kinetics_stable_period,
)

__all__ = [
    "solve",
    "equilibrium",
    "initial_rate",
    "inhour_rho",
    "stable_period",
    "prompt_jump",
    "from_ifp",
]


def from_ifp(betas: list[float], lambda_gen: float, lambdas: list[float]) -> dict[str, Any]:
    """Validated delayed-neutron data from OpenMC IFP kinetics data.

    OpenMC's IFP estimator reports effective delayed fractions (``betas``)
    and the generation time (``lambda_gen``) but no precursor decay
    constants: the caller supplies ``lambdas`` from the same data library.
    Returns ``betas``/``lambdas``/``lambda_gen`` plus ``beta_total`` and
    ``groups``.
    """
    return kinetics_from_ifp(betas, lambda_gen, lambdas)


def solve(
    betas: list[float],
    lambdas: list[float],
    Lambda: float,
    rho: dict[str, Any],
    t: list[float],
    n0: float,
    C0: list[float] | None = None,
    method: str = "trapezoidal",
    rtol: float = 1e-9,
    atol: float = 1e-12,
    dt_min: float = 1e-14,
    dt_max: float | None = None,
    max_steps: int = 1000000,
) -> dict[str, Any]:
    """Solve a prescribed-reactivity point-kinetics transient.

    ``rho`` is a spec dict with ``kind`` selecting ``"constant"`` (``rho``),
    ``"step"`` (``t_step``, ``rho_init``, ``rho_final``), ``"impulse"``
    (``t_start``, ``t_end``, ``rho_init``, ``rho_max``), ``"ramp"``
    (``t_start``, ``t_end``, ``rho_init``, ``rho_rise``, ``rho_final``), or
    ``"polyline"`` (``times``, ``values``); reactivities in Δk, times in
    seconds. Returns ``times``/``n``/``C`` plus the echoed initials.
    """
    return kinetics_solve(
        betas,
        lambdas,
        Lambda,
        rho,
        t,
        n0,
        C0,
        method,
        rtol,
        atol,
        dt_min,
        dt_max,
        max_steps,
    )


def equilibrium(betas: list[float], lambdas: list[float], Lambda: float, n0: float) -> list[float]:
    """Equilibrium precursors ``C_i = beta_i/(lambda_i*Lambda)*n0``."""
    return kinetics_equilibrium(betas, lambdas, Lambda, n0)


def initial_rate(
    betas: list[float],
    lambdas: list[float],
    Lambda: float,
    rho: dict[str, Any],
    n0: float,
    C0: list[float] | None = None,
) -> float:
    """Initial rate ``dn/dt`` at ``t = 0`` for the schedule and initials."""
    return kinetics_initial_rate(betas, lambdas, Lambda, rho, n0, C0)


def inhour_rho(betas: list[float], lambdas: list[float], Lambda: float, omega: float) -> float:
    """Inhour right-hand side ``rho(omega)`` [Δk]."""
    return kinetics_inhour_rho(betas, lambdas, Lambda, omega)


def stable_period(betas: list[float], lambdas: list[float], Lambda: float, rho: float) -> float:
    """Stable period ``T = 1/omega`` [s] at reactivity ``rho`` [Δk]."""
    return kinetics_stable_period(betas, lambdas, Lambda, rho)


def prompt_jump(n_before: float, rho_before: float, rho_after: float, beta_total: float) -> float:
    """Prompt-jump estimate across a reactivity step (needs ``rho_after < beta``)."""
    return kinetics_prompt_jump(n_before, rho_before, rho_after, beta_total)
