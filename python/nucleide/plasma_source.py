"""Tokamak fusion neutron sources (backed by the `nucleide-plasma-source` crate).

Ring, point, and parametric tokamak plasma neutron sources over the D-D
(2.45 MeV) and D-T (14.1 MeV) reactions, with ion-temperature broadening
(Brysk 1973; Ballabio et al. 1998 coefficients). ``particles`` samples a
source spec to particle vectors (seeded, deterministic per platform);
``emit_source_cards`` renders MCNP ``SDEF`` and Serpent ``src`` cards with a
drift report.

``kind="parametric"`` models a Miller-geometry plasma (Fausser et al. 2012)
from caller-supplied profiles: ``major_radius``, ``minor_radius``,
``elongation``, ``triangularity``, ``shafranov_factor`` (cm; shape factors
dimensionless), ``mode`` (``"L"``/``"H"``/``"A"``), ``pedestal_radius`` (cm),
``ion_density_{centre,peaking_factor,pedestal,separatrix}`` (m⁻³ and
dimensionless), and ``ion_temperature_{centre,peaking_factor,beta,pedestal,
separatrix}`` (keV and dimensionless). Emission is reactivity-weighted
(Bosch & Hale 1992). Fuel is the single-fuel ``reaction`` key (``"dt"`` —
equimolar D-T, or ``"dd"`` — pure D-D) or the ``fuel`` dict
``{"D": f_D, "T": f_T}`` (openmc-plasma-source spelling) for a D/T mixture
at the shared profile ion temperature: ``S = n²·[f_D·f_T·⟨σv⟩_DT +
(f_D²/2)·⟨σv⟩_DD]`` with the ``1/(1+δ_ab)`` same-species guard (Eriksson
et al. 2016, the DRESS rate construction). With a ``fuel`` dict the
``reaction`` key is optional and unused. Toroidal sectors are not yet
supported. Lengths are centimetres, energies MeV, ion temperature keV. MCPL
projection stays caller-side: write particle vectors with ``nucleide.mcpl``
when a file is wanted.
"""

from typing import Any

from nucleide._internal import (
    plasma_source_emit_cards,
    plasma_source_particles,
    plasma_source_reactivity,
    plasma_source_spectrum_moments,
)

__all__ = [
    "particles",
    "emit_source_cards",
    "spectrum_moments",
    "reactivity",
]


def particles(spec: dict[str, Any], n: int, seed: int) -> dict[str, Any]:
    """Sample ``n`` source particles into per-field float64 NumPy arrays.

    ``spec`` keys: ``kind`` (``"point"``, ``"ring"``, or ``"parametric"``).
    Point sources take ``position`` [cm] (three-list), ring sources take
    ``radius`` [cm] and ``height`` [cm]; both take ``reaction``
    (``"dt"``/``"dd"``), ``ion_temperature_kev`` (0 for the monoenergetic
    nominal line), and optional ``weight`` (default 1.0). Parametric sources
    take the Miller-geometry and profile keys from the module docstring,
    plus either ``reaction`` (single fuel) or the ``fuel`` fraction dict
    (mixture; then ``reaction`` is optional and unused). Returns ``x``/``y``/
    ``z`` [cm], direction cosines ``u``/``v``/``w`` (unit vectors), ``energy``
    [MeV], and ``weight``. The same ``seed`` reproduces the same stream.
    """
    return plasma_source_particles(spec, n, seed)


def emit_source_cards(spec: dict[str, Any], bins: int = 21) -> dict[str, Any]:
    """Emit MCNP ``SDEF`` and Serpent ``src`` source cards plus drift reports.

    ``spec`` is the source spec from ``particles`` (optional ``mcnp_version``,
    5 or 6, default 5); ``bins`` sets the tabulation bin count. Returns
    ``sdef`` and ``serpent``, each ``{"card": str, "drift": [row dicts]}``
    (rows carry ``quantity``, ``accounted``, ``rel_drift``, ``reparsed``,
    ``note``), plus ``spectrum`` moments (``nominal_mev``, ``mean_mev``,
    ``sigma_mev``, ``mono`` — for a parametric source, the magnetic-axis
    moments; for a fuel mixture, the two-branch Gaussian-mixture summary).
    The SDEF card round-trips through ``nucleide.mcnp.parse_sdef``
    byte-identically; Serpent drift rows are analytic by design. Parametric
    cards carry the radial/vertical/energy *marginals* as histograms; the
    drift rows quantify the tabulation truncation and the joint-correlation
    information a product-form card cannot carry.
    """
    return plasma_source_emit_cards(spec, bins)


def spectrum_moments(reaction: str, ion_temperature_kev: float) -> dict[str, Any]:
    """Closed-form spectrum moments of a fusion reaction at an ion temperature.

    ``reaction`` is ``"dt"`` or ``"dd"``; ``ion_temperature_kev`` is in keV.
    Returns ``reaction``, ``label``, ``nominal_mev`` (the ``T_i = 0`` line),
    ``mean_mev``, and ``sigma_mev`` (0 when monoenergetic).
    """
    return plasma_source_spectrum_moments(reaction, ion_temperature_kev)


def reactivity(reaction: str, ion_temperature_kev: float) -> float:
    """Thermonuclear reactivity ⟨σv⟩ [m³/s] at an ion temperature [keV].

    ``reaction`` is ``"dt"`` or ``"dd"``; the fit is Bosch & Hale (1992) in
    the Atzeni–Meyer-ter-Vehn parametrization. Zero at ``T_i = 0``.
    """
    return plasma_source_reactivity(reaction, ion_temperature_kev)
