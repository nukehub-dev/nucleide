"""Damage and gas-production metrics (backed by the `nucleide-damage` crate).

First-wall engineering metric set by spectral folding of a caller-supplied
multigroup flux with caller-supplied response functions: NRT-dpa
(Norgett–Robinson–Torrens 1975, modified Kinchin–Pease with the Lindhard
damage-energy partition via the Robinson fit), arc-dpa (Nordlund et al.
2018 efficiency correction), He/H production in atomic parts per million,
and He/dpa ratios. The closed-form damage functions
(:func:`nrt_displacements`, :func:`arc_displacements`,
:func:`lindhard_partition`, :func:`damage_energy`,
:func:`arc_efficiency`) are exposed so callers can build per-group response
cross sections with their own nuclear-data pipeline; the folds themselves
are piecewise-constant per group and never re-weight within a group.
Uncertainty propagation (:func:`fold_uq`) reuses the seeded MVN machinery
over caller blocks (the ``nucleide.uq`` engine) with pinned-seed,
``k``-standard-error gates. Caller-supplied response functions stay the
core: the one opt-in exception is the vendored SPECTER Table VII fallback
(:func:`specter_table`, :func:`specter_damage_energy`, :func:`specter_ed`),
spectrum-averaged displacement cross sections transcribed from Greenwood &
Smither, ANL/FPP/TM-197 (US-gov public domain) for callers with no
NJOY/SPECTER-class pipeline of their own — selected explicitly per
spectrum, never consulted implicitly. ASTM E693/E521 are referenced by
designation only, and the SPECTER report doubles as the validation oracle
in the repo harness.
"""

from typing import Any

from nucleide._internal import (
    damage_arc_displacements,
    damage_arc_dpa,
    damage_arc_efficiency,
    damage_coil_accumulate,
    damage_coil_fast_fluence,
    damage_coil_fast_flux,
    damage_coil_lifetime,
    damage_coil_remaining,
    damage_damage_energy,
    damage_fold_uq,
    damage_gas_appm,
    damage_he_dpa_ratio,
    damage_he_dpa_ratio_uq,
    damage_lindhard_partition,
    damage_nrt_displacements,
    damage_nrt_dpa,
    damage_specter_damage_energy,
    damage_specter_ed,
    damage_specter_spectra,
    damage_specter_table,
)

__all__ = [
    "nrt_dpa",
    "arc_dpa",
    "gas_appm",
    "he_dpa_ratio",
    "lindhard_partition",
    "damage_energy",
    "nrt_displacements",
    "arc_efficiency",
    "arc_displacements",
    "fold_uq",
    "he_dpa_ratio_uq",
    "specter_table",
    "specter_damage_energy",
    "specter_ed",
    "specter_spectra",
    "coil_fast_flux",
    "coil_fast_fluence",
    "coil_accumulate",
    "coil_lifetime",
    "coil_remaining",
]

#: Seconds in one full-power year (365.25 d), the life unit of the coil gates.
#: (Module attribute, intentionally outside __all__: the reference generator
#: covers callable facade entries; constants stay importable but unlisted.)
FPY_SECONDS: float = 365.25 * 24.0 * 3600.0


def nrt_dpa(flux: list[float], response: list[float], bounds: list[float], seconds: float) -> float:
    """NRT-dpa over the group flux for ``seconds`` of exposure.

    ``flux[g]`` is the per-group integrated flux (n/cm²/s through group g),
    ``response[g]`` the group NRT dpa cross section in barns (caller-
    condensed, e.g. by a SPECTER/NJOY-class pipeline), and ``bounds`` the
    G+1 MeV group boundaries, strictly increasing. The fold is
    piecewise-constant per group: ``dpa = seconds·1e-24·Σ flux[g]·response[g]``.
    Zero-flux groups contribute exactly 0; negative flux/response and
    non-increasing bounds raise a clear error.
    """
    return damage_nrt_dpa(flux, response, bounds, seconds)


def arc_dpa(flux: list[float], response: list[float], bounds: list[float], seconds: float) -> float:
    """arc-dpa fold: same convention as :func:`nrt_dpa` with arc-corrected
    (``arc_displacements``-based) dpa cross sections."""
    return damage_arc_dpa(flux, response, bounds, seconds)


def gas_appm(
    flux: list[float], response: list[float], bounds: list[float], seconds: float
) -> float:
    """Gas production in atomic parts per million (He or H, whichever gas
    the caller's ``response`` counts), by the same fold with the per-atom
    normalization: ``appm = seconds·1e-18·Σ flux[g]·response[g]``."""
    return damage_gas_appm(flux, response, bounds, seconds)


def he_dpa_ratio(
    flux: list[float],
    he_response: list[float],
    damage_response: list[float],
    bounds: list[float],
    seconds: float,
) -> float:
    """He/dpa ratio in appm per dpa from one fold of the He production and
    damage cross sections over the same flux. A damage fold of exactly zero
    raises a clear error — never ``inf``."""
    return damage_he_dpa_ratio(flux, he_response, damage_response, bounds, seconds)


def lindhard_partition(t_ev: float, recoil: str | int, lattice: str | int) -> float:
    """Lindhard partition fraction ``P(ε) = 1/(1 + k_L·g(ε))``: the fraction
    of a recoil's energy available for displacements after electronic losses.

    ``t_ev`` is the PKA energy in eV; ``recoil`` and ``lattice`` accept a
    nuclide name (``"Fe56"``) or an int nucid (self-recoil: pass the same
    value twice). Robinson-fit constants (Lindhard 1963; Robinson 1970).
    """
    return damage_lindhard_partition(t_ev, recoil, lattice)


def damage_energy(t_ev: float, recoil: str | int, lattice: str | int) -> float:
    """Lindhard damage energy ``T_dam = T·P(ε)`` in eV for a recoil of energy
    ``t_ev`` stopped in ``lattice``."""
    return damage_damage_energy(t_ev, recoil, lattice)


def nrt_displacements(t_ev: float, ed_ev: float, target: str | int) -> float:
    """NRT displacement function ``N_d(T)`` for a self-recoil ``target`` with
    average threshold displacement energy ``ed_ev`` (eV): 0 below ``E_d``,
    1 on ``E_d <= T < 2·E_d/0.8``, and ``0.8·T_dam/(2·E_d)`` above — the
    Norgett–Robinson–Torrens (1975) standard."""
    return damage_nrt_displacements(t_ev, ed_ev, target)


def arc_efficiency(t_dam_ev: float, ed_ev: float, b_arc: float, c_arc: float) -> float:
    """arc-dpa efficiency ``ξ(T_d) = (1−c)·(T_d/(2·E_d/0.8))^b + c``
    (Nordlund et al. 2018, Eq. (7)) at damage energy ``t_dam_ev``. ``b_arc``
    must be negative and ``c_arc`` in (0, 1) — the MD-fitted regime; the
    constants are caller-supplied material data, not vendored tables."""
    return damage_arc_efficiency(t_dam_ev, ed_ev, b_arc, c_arc)


def arc_displacements(
    t_ev: float, ed_ev: float, target: str | int, b_arc: float, c_arc: float
) -> float:
    """arc-dpa displacement function for a self-recoil ``target``: the NRT
    piecewise form with the high branch multiplied by ``ξ_arc(T_dam)``."""
    return damage_arc_displacements(t_ev, ed_ev, target, b_arc, c_arc)


def fold_uq(
    metric: str,
    flux: list[float],
    response: list[float],
    bounds: list[float],
    seconds: float,
    mean: list[float],
    cov: list[list[float]],
    n: int,
    seed: int,
    k: float,
) -> dict[str, Any]:
    """Uncertainty propagation through one spectral fold (the ``nucleide.uq``
    MVN engine, no new sampling machinery).

    ``mean``/``cov`` describe relative perturbations of the stacked
    ``[flux, response]`` vector (dimension 2G). ``metric`` is one of
    ``"nrt_dpa"``, ``"arc_dpa"``, ``"gas_appm"`` (``"he_dpa_ratio"`` is a
    loud error naming :func:`he_dpa_ratio_uq`: the ratio needs both
    responses). Draws are seeded and
    reproducible; the returned dict carries the sample moments (``mean``,
    ``std``), the exact expectation (``expected``), the first-order
    propagated standard deviation (``analytic_std``), and the ``k``-SE gate
    verdict (``passed``) — the same honest-gate pattern as the U1–U4/U7
    samplers.
    """
    return damage_fold_uq(metric, flux, response, bounds, seconds, mean, cov, n, seed, k)


def he_dpa_ratio_uq(
    flux: list[float],
    he_response: list[float],
    damage_response: list[float],
    bounds: list[float],
    seconds: float,
    mean: list[float],
    cov: list[list[float]],
    n: int,
    seed: int,
    k: float,
) -> dict[str, Any]:
    """Uncertainty propagation through the He/dpa ratio (per-draw ratios).

    ``mean``/``cov`` describe relative perturbations of the stacked
    ``[flux, he_response, damage_response]`` vector (dimension 3G). Each
    seeded draw refolds He and dpa and forms the ratio per draw; the
    returned dict carries the draw moments (``mean``, ``std``), the
    second-order bias-corrected expectation (``expected`` — the mean of
    ratios carries the ``Var(dpa)`` bias, so the bare ratio of means is
    reported as ``nominal`` instead), the first-order delta-propagated
    standard deviation (``analytic_std``), and the ``k``-SE gate verdict
    (``passed``). Because a skewed ratio needs asymmetric margins, the dict
    also carries the distribution-free 68% interval — draw percentiles
    (``q16``, ``q50``, ``q84``) with the Fieller-construction quantiles
    (``expected_q16``, ``expected_q50``, ``expected_q84``) and their own
    ``k``-SE gate verdict (``quantiles_passed``), which holds at
    honestly-sized blocks where the symmetric-spread gate cannot. A draw
    landing at non-positive dpa fails loudly, never ``inf``/``NaN`` with a
    spread.
    """
    return damage_he_dpa_ratio_uq(
        flux, he_response, damage_response, bounds, seconds, mean, cov, n, seed, k
    )


def specter_spectra() -> list[str]:
    """Canonical spectrum names of the vendored SPECTER Table VII fallback:
    ``thermal``, ``fission``, ``14mev``, ``hfir``, ``ebr2``, ``fftf``,
    ``fusion``."""
    return damage_specter_spectra()


def specter_table(spectrum: str) -> dict[str, float]:
    """Vendored SPECTER Table VII displacement cross sections in barns for
    one spectrum, keyed by element symbol (``"Fe"``; ``"Ag"`` is natural
    silver, ``"W"`` natural tungsten).

    Spectrum-averaged damage-energy cross sections (keV-b, Greenwood &
    Smither, ANL/FPP/TM-197, US-gov public domain) converted with the
    vendored Table II ``E_d`` via the report's ``0.8/2E_d`` rule — an
    opt-in fallback response for callers with no NJOY/SPECTER-class
    pipeline of their own. The values are whole-spectrum averages, so each
    is a one-group response for its spectrum (fold with a one-group
    fluence, ``seconds=1``). Unknown spectra raise a clear error.
    """
    return damage_specter_table(spectrum)


def specter_damage_energy(spectrum: str) -> dict[str, float]:
    """Verbatim vendored SPECTER Table VII damage-energy cross sections in
    keV-b (as printed) for one spectrum, keyed by element symbol — same
    transcription and spectrum spellings as :func:`specter_table` without
    the ``0.8/2E_d`` conversion."""
    return damage_specter_damage_energy(spectrum)


def specter_ed(element: str) -> float:
    """Vendored SPECTER Table II ``E_d`` (eV) for an element symbol
    (e.g. ``"Fe"``). Unknown elements raise a clear error."""
    return damage_specter_ed(element)


def coil_fast_flux(flux: list[float], bounds: list[float], threshold_mev: float) -> float:
    """Fast flux above ``threshold_mev``: the group sum over groups with
    ``bounds[g+1] > threshold_mev`` (caller area units, e.g. m⁻²s⁻¹).

    A threshold cutting through a group includes the whole group — align the
    grid so the threshold sits on a boundary when that conservatism matters.
    """
    return damage_coil_fast_flux(flux, bounds, threshold_mev)


def coil_fast_fluence(
    flux: list[float], bounds: list[float], threshold_mev: float, seconds: float
) -> float:
    """Fast fluence above ``threshold_mev``: :func:`coil_fast_flux` held for
    ``seconds`` (fluence = flux × time, caller area units)."""
    return damage_coil_fast_fluence(flux, bounds, threshold_mev, seconds)


def coil_accumulate(rates: list[float], durations: list[float]) -> float:
    """Piecewise-constant irradiation-history accumulation:
    ``Σ rates[i]·durations[i]``.

    Fold each interval's spectrum through its metric of choice
    (:func:`coil_fast_flux` or a dpa fold at one second for the rate) and
    accumulate the rate history here."""
    return damage_coil_accumulate(rates, durations)


def coil_lifetime(limits: list[float], rates: list[float]) -> dict[str, Any]:
    """Weakest-link coil life over caller-supplied limit tables:
    ``min`` of ``limits[i] / rates[i]``.

    Limits share their rate's units (fluence or dpa — any accumulated
    metric); published design numbers are validation gates the caller
    supplies, never defaults. Returns ``{"seconds", "limiting"}``: seconds
    to the first breach plus the limiting channel index (all-zero rates give
    infinite seconds with ``limiting`` None — a zero rate never breaches).
    """
    return damage_coil_lifetime(limits, rates)


def coil_remaining(
    limits: list[float], accumulated: list[float], rates: list[float]
) -> dict[str, Any]:
    """Weakest-link remaining life: ``min`` of
    ``(limits[i] − accumulated[i]) / rates[i]``, clamped at zero once a
    limit is reached. Same ``{"seconds", "limiting"}`` shape as
    :func:`coil_lifetime`."""
    return damage_coil_remaining(limits, accumulated, rates)
