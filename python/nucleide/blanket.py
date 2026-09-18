"""TBR and blanket power bookkeeping (backed by the `nucleide-blanket` crate).

Pure arithmetic over the caller's transport tallies — TBR with per-port
coverage penalties, blanket energy multiplication, tritium burn rate, and
the breeding-margin / net-surplus fuel-cycle metrics. Tallies come from
the caller's OpenMC/DAGMC runs; this module reports margins and never
optimizes geometry, solves transport, or touches the `tritium`
permeation kernel.
"""

from nucleide._internal import (
    blanket_apply_port_penalty,
    blanket_blanket_power,
    blanket_breeding_margin,
    blanket_energy_multiplication,
    blanket_meets_requirement,
    blanket_net_surplus,
    blanket_tbr_from_tallies,
    blanket_tritium_burn,
)

__all__ = [
    "tbr_from_tallies",
    "apply_port_penalty",
    "breeding_margin",
    "meets_requirement",
    "energy_multiplication",
    "blanket_power",
    "tritium_burn_g_per_day",
    "net_surplus_g_per_day",
]


def tbr_from_tallies(tritons_bred: float, source_neutrons: float) -> float:
    """Raw TBR from caller tallies: ``tritons_bred / source_neutrons``.

    ``tritons_bred`` is the caller-tallied tritons bred (non-negative),
    ``source_neutrons`` the source neutrons (strictly positive). Returns
    the dimensionless ratio.
    """
    return blanket_tbr_from_tallies(tritons_bred, source_neutrons)


def apply_port_penalty(raw_tbr: float, port_fractions: list[float]) -> float:
    """Effective TBR after per-port coverage penalties.

    Each entry of ``port_fractions`` is one port's fractional coverage
    loss in ``[0, 1)``; the haircut is multiplicative,
    ``raw * Prod_i (1 - f_i)``. An empty list returns ``raw_tbr``
    unchanged.
    """
    return blanket_apply_port_penalty(raw_tbr, port_fractions)


def breeding_margin(effective_tbr: float) -> float:
    """Breeding margin: ``effective_tbr - 1`` (dimensionless).

    A negative return is a sub-breakeven deficit and is reported as-is,
    never clamped to zero.
    """
    return blanket_breeding_margin(effective_tbr)


def meets_requirement(effective_tbr: float, required_tbr: float) -> bool:
    """Return whether ``effective_tbr >= required_tbr``.

    ``required_tbr`` must be strictly positive (breakeven itself is
    spelled ``1.0``).
    """
    return blanket_meets_requirement(effective_tbr, required_tbr)


def energy_multiplication(blanket_power_mw: float, fusion_power_mw: float) -> float:
    """Blanket energy multiplication: ``blanket_power_mw / fusion_power_mw``.

    Both powers are in MW; blanket power is non-negative, fusion power
    strictly positive.
    """
    return blanket_energy_multiplication(blanket_power_mw, fusion_power_mw)


def blanket_power(fusion_power_mw: float, multiplication: float) -> float:
    """Blanket thermal power in MW: ``multiplication * fusion_power_mw``."""
    return blanket_blanket_power(fusion_power_mw, multiplication)


def tritium_burn_g_per_day(fusion_power_mw: float) -> float:
    """Tritium burn rate in g/day: one triton mass per 17.6 MeV of D-T
    fusion energy, scaled to ``fusion_power_mw`` MW of fusion power."""
    return blanket_tritium_burn(fusion_power_mw)


def net_surplus_g_per_day(effective_tbr: float, fusion_power_mw: float) -> float:
    """Net tritium surplus in g/day: ``(effective_tbr - 1) * burn``.

    Positive is a breeding surplus, negative a net deficit (reported
    as-is).
    """
    return blanket_net_surplus(effective_tbr, fusion_power_mw)
