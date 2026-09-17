"""1D tritium diffusion-trapping kernel (backed by the `nucleide-tritium` crate).

Single slabs (:func:`steady`, :func:`transient`) and multi-layer series
stacks with Sieverts/Henry/vented-sink-recombination internal interfaces
(:func:`steady_layers`, :func:`transient_layers`)."""

from typing import Any

from nucleide._internal import (
    tritium_breakthrough,
    tritium_irreversible_fill,
    tritium_langmuir,
    tritium_layers_steady,
    tritium_layers_transient,
    tritium_oriani,
    tritium_recombination_rate,
    tritium_sieverts,
    tritium_steady,
    tritium_time_lag,
    tritium_transient,
)

__all__ = [
    "steady",
    "transient",
    "steady_layers",
    "transient_layers",
    "time_lag",
    "breakthrough",
    "oriani",
    "langmuir",
    "irreversible_fill",
    "sieverts",
    "recombination_rate",
]


def steady(
    length: float,
    cells: int,
    D: float,
    left: dict[str, Any],
    right: dict[str, Any],
    traps: list[dict[str, Any]] | None = None,
    temperature: list[float] | None = None,
    source: list[float] | None = None,
    E_D: float = 0.0,
) -> dict[str, Any]:
    """Trap-free-style steady state of the mobile/trapped slab.

    ``left``/``right`` are boundary-spec dicts with ``kind`` selecting
    ``"dirichlet"`` (``value`` [mol/m³]), ``"sieverts"``/``"henry"``
    (``solubility``, ``pressure`` [Pa]), ``"recombination"`` (``rate`` —
    closed per solve: the G5 face construction in steady state, the G6
    per-step face Newton in :func:`transient`), or ``"zero_flux"``.
    ``traps`` holds one spec dict per species (``k0``, ``p0``,
    ``site_density`` required; ``e_k``/``e_p`` default to 0). ``temperature``
    is one value (uniform, default 500 K) or one per cell; ``source`` is
    ``None`` (zero), one value, or one per cell. Returns ``centres``,
    ``mobile``, ``trapped`` (``[cell][trap]``), ``flux_left``/``flux_right``
    (outward-positive), and the two inventories.
    """
    return tritium_steady(
        length,
        cells,
        D,
        E_D,
        traps if traps is not None else [],
        temperature if temperature is not None else [500.0],
        source,
        left,
        right,
    )


def transient(
    length: float,
    cells: int,
    D: float,
    left: dict[str, Any],
    right: dict[str, Any],
    t: list[float],
    traps: list[dict[str, Any]] | None = None,
    temperature: list[float] | None = None,
    source: list[float] | None = None,
    E_D: float = 0.0,
    mobile0: list[float] | None = None,
    trapped0: list[list[float]] | None = None,
    method: str = "crank_nicolson",
    rtol: float = 1e-9,
    atol: float = 1e-12,
    dt_min: float = 1e-14,
    dt_max: float | None = None,
    max_steps: int = 1000000,
) -> dict[str, Any]:
    """Solve the mobile/trapped transient over the output grid ``t`` [s].

    Same slab/trap/BC arguments as :func:`steady` plus the output times and
    the optional initial profiles (both default to zero). ``method`` is
    ``"crank_nicolson"`` (default) or ``"backward_euler"``. Returns
    ``times``, ``mobile`` (``[time][cell]``), ``trapped``
    (``[time][cell][trap]``), and the outward ``flux_left``/``flux_right``
    series.
    """
    return tritium_transient(
        length,
        cells,
        D,
        E_D,
        traps if traps is not None else [],
        temperature if temperature is not None else [500.0],
        source,
        left,
        right,
        t,
        mobile0,
        trapped0,
        method,
        rtol,
        atol,
        dt_min,
        dt_max,
        max_steps,
    )


def time_lag(length: float, D: float) -> float:
    """Permeation time lag ``t_lag = L²/6D`` [s] (G2-lag)."""
    return tritium_time_lag(length, D)


def breakthrough(D: float, length: float, times: list[float]) -> list[float]:
    """Normalized outlet flux ``J(L,t)/J_ss`` at each time (G2 series)."""
    return tritium_breakthrough(D, length, times)


def oriani(D: float, K: float, N: float) -> float:
    """Oriani effective diffusivity ``D_eff = D/(1 + K N)`` [m²/s] (G3a)."""
    return tritium_oriani(D, K, N)


def langmuir(N: float, K: float, c: float) -> float:
    """Langmuir equilibrium load ``c_t = N K c/(1 + K c)`` [mol/m³] (T2-eq)."""
    return tritium_langmuir(N, K, c)


def irreversible_fill(k: float, c: float, N: float, times: list[float]) -> list[float]:
    """Irreversible-trap fill ``c_t(t) = N(1 − e^{−kct})`` [mol/m³] (G3c)."""
    return tritium_irreversible_fill(k, c, N, times)


def sieverts(K_S: float, p: float) -> float:
    """Sieverts surface concentration ``c = K_S sqrt(p)`` [mol/m³] (G4)."""
    return tritium_sieverts(K_S, p)


def recombination_rate(kr0: float, e_r: float, temp: float) -> float:
    """Recombination rate ``K_r = kr0 * exp(-e_r / R / temp)`` [m⁴/mol/s] (G5)."""
    return tritium_recombination_rate(kr0, e_r, temp)


def steady_layers(
    layers: list[dict[str, Any]],
    left: dict[str, Any],
    right: dict[str, Any],
    interfaces: list[Any] | None = None,
) -> dict[str, Any]:
    """Trap-free-style steady state of a multi-layer series stack (G7/G9/G10).

    ``layers`` holds one spec dict per layer, from the ``left`` face to the
    ``right`` face: ``thickness`` [m], ``cells``, ``D`` [m²/s],
    ``solubility`` (required) — the layer's interface constant ``K``
    (``K_S`` [mol/m³/Pa¹ᐟ²] under a Sieverts law, ``K_H`` [mol/m³/Pa] under a
    Henry law; unused at vented-sink gaps, whose face is a concentration);
    ``E_D`` [J/mol] (default 0), ``traps`` (list of trap-spec
    dicts, default none), ``temperature`` (one value or one per cell of the
    layer, default ``[500]``), and ``source`` (``None``, one value, or one
    per cell). ``interfaces`` optionally holds one internal-interface entry
    per gap — a ``"sieverts"`` (default) or ``"henry"`` (G9) string, or a
    ``{"kind": "recombination", "rate": Kr}`` dict for the vented-sink law
    (G10: single face concentration with the ``Kr·x²`` desorption sink;
    ``Kr → 0`` is a continuous-concentration joint, not a Sieverts law; a
    non-positive permeation drive is a clear error). Linear gaps hold
    ``c/K`` continuous with continuous flux, so a stack whose layers share
    ``D`` and ``K`` is a plain slab with a transparent interface.
    ``left``/``right`` are the same boundary-spec dicts as :func:`steady`.
    A one-layer stack reproduces :func:`steady` exactly. Returns
    ``centres``, ``mobile``, ``trapped`` (``[cell][trap]``, per-layer
    species), ``flux_left``/``flux_right`` (outward-positive),
    ``interface_faces`` (one face value per gap, ``None`` at linear gaps),
    and the two inventories.
    """
    return tritium_layers_steady(layers, left, right, interfaces)


def transient_layers(
    layers: list[dict[str, Any]],
    left: dict[str, Any],
    right: dict[str, Any],
    t: list[float],
    mobile0: list[float] | None = None,
    trapped0: list[list[float]] | None = None,
    method: str = "crank_nicolson",
    rtol: float = 1e-9,
    atol: float = 1e-12,
    dt_min: float = 1e-14,
    dt_max: float | None = None,
    max_steps: int = 1000000,
    interfaces: list[Any] | None = None,
) -> dict[str, Any]:
    """Solve the multi-layer mobile/trapped transient over the grid ``t`` [s] (G8/G11).

    Same layer-stack and boundary arguments as :func:`steady_layers`
    (including the optional ``interfaces`` entry per gap) plus the output
    times and the optional initial profiles (``mobile0`` per cell,
    ``trapped0`` as ``[cell][trap]`` matching each layer's trap count; both
    default to zero). ``method`` is ``"crank_nicolson"`` (default) or
    ``"backward_euler"``. Returns ``times``, ``mobile``
    (``[time][cell]``), ``trapped`` (``[time][cell][trap]``), the outward
    ``flux_left``/``flux_right`` series, and ``interface_faces``
    (``[time][gap]``, ``None`` at linear gaps) for closing the discrete
    balance with the desorption term.
    """
    return tritium_layers_transient(
        layers,
        left,
        right,
        t,
        mobile0,
        trapped0,
        method,
        rtol,
        atol,
        dt_min,
        dt_max,
        max_steps,
        interfaces,
    )
