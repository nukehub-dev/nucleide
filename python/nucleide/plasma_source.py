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
``reaction`` key is optional and unused. An optional ``species_temperatures``
dict ``{"D": T_D, "T": T_T}`` [keV] reacts at distinct Maxwellian species
temperatures — the D-T arm at the mass-weighted
``T_DT = T_D + (2/5)·(T_T − T_D)``, the D-D arm at ``T_D`` (Eriksson et
al. 2016, the full distinct-Maxwellian generalization); on a mixture the
profile ion temperature is then unused for rate and spectrum (still
validated), while on a single-fuel config the same pair rule applies
(single-fuel D-T at ``T_DT``, single-fuel D-D at ``T_D``). A uniform
``T_D = T_T`` pair reproduces the shared-temperature kernel exactly.
Non-finite or negative temperatures are loud errors. An optional
``deuterium_tail`` dict ``{"fraction": eta, "temperature_kev": T_tail}``
splits the deuterium population into a bulk
``(1 − eta)`` at ``T_D`` and a hot tail ``eta`` at ``T_tail`` [keV] — the one
pinned single-tail-temperature non-Maxwellian shape (Eriksson et al. 2016
arbitrary-distribution framework); each bulk/tail sub-pair reacts at its own
mass-weighted relative temperature (D-T at ``T_DT``/``T_DTt``, D-D at
``T_D``/``T_mix``/``T_tail``) with the Ballabio lines following per
sub-branch, and ``eta = 0`` reproduces the no-tail kernel exactly. Requires a
``fuel`` dict; non-finite, out-of-range, or negative parameters are loud
errors. An optional
``start_angle``/``rotation_angle`` pair [rad] restricts births to the toroidal
sector ``[start_angle, start_angle + rotation_angle)`` (uniform birth angle;
emission totals scale by ``rotation_angle/2π``); both keys or neither, and a
full ``2π`` rotation reproduces the full torus exactly. ``proton_accounting``
reports the D(d,p)T proton rate alongside the neutron source (pinned 50/50
with the D-D neutron branch; sampler and cards stay neutron-only). Lengths are
centimetres, energies MeV, ion temperature keV. MCPL
projection stays caller-side: write particle vectors with ``nucleide.mcpl``
when a file is wanted.

``kind="lattice"`` drops the last tokamak assumption (axisymmetry): the caller
supplies the full 3D birth distribution as ``points``, a list of
``{"position": [x, y, z] [cm], "rate": w >= 0, "ion_temperature_kev": Ti}``
dicts (relative birth-rate weights, arbitrary global scale; ``Ti = 0`` is the
monoenergetic nominal line). Fuel is ``reaction`` or the ``fuel`` dict (then
``reaction`` is optional and unused), with the same per-node branch roulette;
the optional ``species_temperatures`` pair reacts every node at ``(T_D, T_T)``
(requires ``fuel``; node temperatures then unused for rate/spectrum, still
validated). The optional ``field_periods``/``base_angle`` pair declares
field-period symmetry (both or neither; ``field_periods >= 1``): the points
are the base-sector cloud, replicated uniformly around the axis. Cards carry
the cloud's cylindrical-``R``/vertical/energy marginals with drift rows plus
a ``lattice discretization`` row. ``proton_accounting`` reports the D-D branch
share of the caller rates (a pure-tritium lattice still carries caller neutron
strength, unlike the reactivity-derived parametric kernel).
"""

from typing import Any

from nucleide._internal import (
    plasma_source_ecrh_beamline,
    plasma_source_ecrh_scalars,
    plasma_source_emit_cards,
    plasma_source_particles,
    plasma_source_proton_accounting,
    plasma_source_reactivity,
    plasma_source_spectrum_moments,
)

__all__ = [
    "particles",
    "emit_source_cards",
    "spectrum_moments",
    "reactivity",
    "proton_accounting",
    "ecrh_scalars",
    "ecrh_accessibility",
]


def particles(spec: dict[str, Any], n: int, seed: int) -> dict[str, Any]:
    """Sample ``n`` source particles into per-field float64 NumPy arrays.

    ``spec`` keys: ``kind`` (``"point"``, ``"ring"``, ``"parametric"``, or
    ``"lattice"``).
    Point sources take ``position`` [cm] (three-list), ring sources take
    ``radius`` [cm] and ``height`` [cm]; both take ``reaction``
    (``"dt"``/``"dd"``), ``ion_temperature_kev`` (0 for the monoenergetic
    nominal line), and optional ``weight`` (default 1.0). Parametric sources
    take the Miller-geometry and profile keys from the module docstring,
    plus either ``reaction`` (single fuel) or the ``fuel`` fraction dict
    (mixture; then ``reaction`` is optional and unused), plus the optional
    ``species_temperatures`` dict ``{"D": T_D, "T": T_T}`` [keV] for distinct
    species temperatures (both keys; single-fuel D-T reacts at ``T_DT``,
    single-fuel D-D at ``T_D``, mixtures at the same pair temperatures), plus
    the optional ``deuterium_tail`` dict ``{"fraction": eta,
    "temperature_kev": T_tail}`` for the deuterium hot-tail fraction on a
    mixture (both keys; requires ``fuel``), plus
    the optional ``start_angle``/``rotation_angle`` toroidal-sector pair
    [rad] (both or neither). Lattice sources take ``points`` (a list of
    ``{"position", "rate", "ion_temperature_kev"}`` dicts), plus either
    ``reaction`` or the ``fuel`` dict, plus the optional
    ``species_temperatures`` dict, plus the optional
    ``field_periods``/``base_angle`` symmetry pair (both or neither).
    Returns ``x``/``y``/
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
    moments; for a fuel mixture, the two-branch Gaussian-mixture summary;
    for a lattice, the rate-weighted cloud moments with the mixture sigma
    around the mean).
    The SDEF card round-trips through ``nucleide.mcnp.parse_sdef``
    byte-identically; Serpent drift rows are analytic by design. Parametric
    cards carry the radial/vertical/energy *marginals* as histograms (plus a
    uniform toroidal-angle marginal for a partial sector); lattice cards
    carry the cloud's cylindrical-``R``/vertical/energy marginals plus a
    ``lattice discretization`` row; the
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


def proton_accounting(spec: dict[str, Any]) -> dict[str, Any]:
    """D(d,p)T proton bookkeeping for a source spec (accounting only).

    The proton branch shares the D-D total rate 50/50 with the neutron
    branch, so the proton density IS the D-D neutron-branch density.
    Returns ``proton_per_neutron`` (exactly 1.0 for pure D-D, 0.0 for
    D-T-only) plus, for ``kind="parametric"`` only, the integrated
    ``neutron_total`` and ``proton_total`` (sector-aware; ring/point specs
    carry no density model, so their totals are ``None``) and a ``note``
    stating the convention. The sampler and the emitted cards stay
    neutron-only; proton transport is out of scope. A pure-tritium mixture
    (no T-T neutron branch) is a loud error, like sampling it.
    """
    return plasma_source_proton_accounting(spec)


def ecrh_scalars(
    frequency_ghz: float,
    harmonic: int = 1,
    electron_temperature_kev: float | None = None,
    field_t: float | None = None,
) -> dict[str, Any]:
    """Closed-form ECRH resonance and cut-off scalars at a wave frequency.

    ``frequency_ghz`` is the wave frequency [GHz], ``harmonic`` the cyclotron
    harmonic (≥ 1), ``electron_temperature_kev`` (or ``None``) the bulk
    Maxwell–Jüttner shift input [keV], and ``field_t`` (or ``None``) the
    local field [T]. Returns ``f_ce_per_t_ghz`` (≈ 27.992 GHz/T),
    ``f_ce_ghz`` (``None`` without ``field_t``), ``b_cold_t`` (cold
    resonance, ECRH-2), ``b_rel_t`` (``None`` without a temperature, ECRH-3),
    ``n_o1_m3`` (O-mode cut-off density, ECRH-4), and ``n_x1_m3`` (``None``
    without ``field_t`` or where the X1 wave is evanescent, ECRH-5). Plain
    floats; malformed inputs raise ``ValueError``.
    """
    return plasma_source_ecrh_scalars(frequency_ghz, harmonic, electron_temperature_kev, field_t)


def ecrh_accessibility(
    s_m: list[float],
    b_t: list[float],
    ne_m3: list[float],
    frequency_ghz: float,
    harmonic: int = 1,
    electron_temperature_kev: float | None = None,
) -> dict[str, Any]:
    """Locate ECRH resonance and cut-offs along a caller beamline.

    ``s_m``/``b_t``/``ne_m3`` are the caller position [m], field [T], and
    density [m⁻³] profiles at strictly increasing positions;
    ``frequency_ghz``/``harmonic`` select the wave;
    ``electron_temperature_kev`` (or ``None``) selects the relativistically
    shifted resonance. Returns ``b_res_t`` (the located resonant field),
    ``b_cold_t``, ``n_o1_m3``, and the ``resonance_m``/``o1_cutoff_m``/
    ``x1_cutoff_m`` crossing positions [m] — the per-port penalty inputs.
    No ray tracing; ragged or non-monotonic beamlines raise ``ValueError``.
    """
    return plasma_source_ecrh_beamline(
        s_m, b_t, ne_m3, frequency_ghz, harmonic, electron_temperature_kev
    )
