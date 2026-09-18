"""Equilibrium data readers (backed by the `nucleide-equilib-io` crate)."""

from typing import Any

from nucleide._internal import (
    equilib_jacobian,
    equilib_jacobian_grid,
    equilib_parse_indata,
    equilib_probe_variant,
    equilib_read_wout,
    equilib_wall_load,
)

__all__ = [
    "probe_variant",
    "read_wout",
    "jacobian",
    "jacobian_grid",
    "wall_load",
    "parse_indata",
    "indata_int",
    "indata_float",
    "indata_bool",
    "indata_str",
    "indata_series_1d",
    "indata_table_2d",
    "is_free_boundary",
    "ncurr_is_iprime",
]


def probe_variant(path: str) -> str:
    """Probe a ``wout`` file variant from its magic bytes.

    Returns ``"classic"`` for CDF-1 and ``"64-bit offset"`` for CDF-2
    (the ``ncdump -k`` spellings). HDF5-backed netCDF-4 files and CDF-5
    files raise ``ValueError`` pointing at facade-side conversion: the
    Rust reader never takes HDF5.
    """
    return equilib_probe_variant(path)


def read_wout(path: str) -> dict[str, Any]:
    """Read a classic-netCDF ``wout`` file (CDF-1/CDF-2 only).

    Returns the reader-minimum variables (``nfp``, ``ns``, ``xm``,
    ``xn``, ``rmnc``, ``zmns``, ``lmns``, ``gmnc``), the resolved
    counts, optional scalars (``mpol``, ``ntor``, ``phiedge``,
    ``volume_p``, each ``None`` when absent), the optional Nyquist maps
    (``xm_nyq``/``xn_nyq``), the later-use fields present in the file,
    and the ``variant`` spelling. Matrices arrive as ``{"rows", "cols",
    "values"}`` dicts with nested row lists (outer over radius, inner
    over modes). Reads data, never solves equilibria.
    """
    return equilib_read_wout(path)


def jacobian(
    xm: list[float], xn: list[float], coeffs: list[float], theta: float, zeta: float
) -> float:
    """Evaluate the Jacobian Fourier sum ``Σ c·cos(m·θ − n·ζ)`` (J1).

    ``xm``/``xn`` are the mode maps, ``coeffs`` one half-mesh ``gmnc``
    row, angles in radians with ``zeta`` in the file's toroidal
    convention (one field period — the caller scales for full-torus
    work).
    """
    return equilib_jacobian(xm, xn, coeffs, theta, zeta)


def jacobian_grid(
    xm: list[float],
    xn: list[float],
    coeffs: list[float],
    nfp: int,
    ntheta: int,
    nzeta: int,
) -> list[list[float]]:
    """Evaluate the Jacobian on a tensor grid (J2).

    ``ntheta`` poloidal points over ``[0, 2π)`` by ``nzeta`` toroidal
    points over one field period ``[0, 2π/nfp)``; returns nested row
    lists (outer over θ, inner over ζ).
    """
    return equilib_jacobian_grid(xm, xn, coeffs, nfp, ntheta, nzeta)


def wall_load(
    s_edges: list[float],
    ntheta: int,
    nzeta: int,
    nfp: int,
    birth: list[float],
    jacobian: list[float],
) -> dict[str, Any]:
    """Map a birth-rate density field onto the wall in flux coordinates (W1-W2).

    ``s_edges`` holds the ``nr + 1`` strictly increasing radial edges,
    ``ntheta``/``nzeta`` the angular grid counts, ``nfp`` the field-period
    count; ``birth``/``jacobian`` are the per-voxel birth-rate densities
    and ``sqrt(g)`` values in radial-major order
    (``(i * ntheta + j) * nzeta + k`` — evaluate the Jacobian voxels with
    :func:`jacobian`/:func:`jacobian_grid`). Returns ``{"ntheta",
    "nzeta", "loads", "total"}`` with nested ``loads`` rows (outer over
    theta, inner over zeta; one radial line integral per wall node) and
    the one-field-period ``total``. Closed-form accumulation only — no
    transport, no shadowing. Shape mismatches, empty grids, and
    non-finite/negative densities or Jacobians raise ``ValueError``.
    """
    return equilib_wall_load(s_edges, ntheta, nzeta, nfp, birth, jacobian)


def parse_indata(text: str) -> dict[str, Any]:
    """Parse ``&INDATA ... /`` namelist text.

    Returns ``{"scalars": {NAME: value}, "indexed": {NAME: [[[i, ...],
    value], ...]}}`` with upper-cased names. Sliced indices
    (``RBC(0:4,2)=...``) raise ``ValueError``.
    """
    return equilib_parse_indata(text)


def _scalars(parsed: dict[str, Any]) -> dict[str, Any]:
    return dict(parsed["scalars"])


def indata_int(parsed: dict[str, Any], name: str) -> int:
    """Integer switch by name (``NFP``, ``NCURR``, ``MPOL``, ``NTOR``, ...)."""
    value = _scalars(parsed).get(name.upper())
    if value is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    if isinstance(value, bool):
        raise ValueError(f"`{name.upper()}` is not an integer")
    if isinstance(value, float):
        if not value.is_integer():
            raise ValueError(f"`{name.upper()}` is not an integer")
        return int(value)
    if not isinstance(value, int):
        raise ValueError(f"`{name.upper()}` is not an integer")
    return value


def indata_float(parsed: dict[str, Any], name: str) -> float:
    """Real value by name."""
    value = _scalars(parsed).get(name.upper())
    if value is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"`{name.upper()}` is not a float")
    return float(value)


def indata_bool(parsed: dict[str, Any], name: str) -> bool:
    """Logical switch by name (``LFREEB``, ``LASYM``, ...)."""
    value = _scalars(parsed).get(name.upper())
    if value is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    if not isinstance(value, bool):
        raise ValueError(f"`{name.upper()}` is not a logical")
    return value


def indata_str(parsed: dict[str, Any], name: str) -> str:
    """String value by name (``MGRID_FILE``, ``*_TYPE``, ...)."""
    value = _scalars(parsed).get(name.upper())
    if value is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    if not isinstance(value, str):
        raise ValueError(f"`{name.upper()}` is not a string")
    return value


def indata_series_1d(parsed: dict[str, Any], name: str) -> list[tuple[int, float]]:
    """One-dimensional float series (``AM``, ``RAXIS``, ...): ``(index, value)`` pairs."""
    entries = parsed["indexed"].get(name.upper())
    if entries is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    out = []
    for indices, value in entries:
        if len(indices) != 1 or isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"`{name.upper()}` entry is not a 1-D float")
        out.append((int(indices[0]), float(value)))
    return sorted(out)


def indata_table_2d(parsed: dict[str, Any], name: str) -> list[tuple[tuple[int, int], float]]:
    """Two-dimensional float table (``RBC``, ``ZBS``, ...): ``((m, n), value)`` entries."""
    entries = parsed["indexed"].get(name.upper())
    if entries is None:
        raise ValueError(f"missing variable `{name.upper()}`")
    out = []
    for indices, value in entries:
        if len(indices) != 2 or isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"`{name.upper()}` entry is not a 2-D float")
        out.append(((int(indices[0]), int(indices[1])), float(value)))
    return sorted(out)


def is_free_boundary(parsed: dict[str, Any]) -> bool:
    """Whether ``LFREEB`` is set; absent reads fixed boundary (the only conversion stance)."""
    value = _scalars(parsed).get("LFREEB", False)
    if not isinstance(value, bool):
        raise ValueError("`LFREEB` is not a logical")
    return value


def ncurr_is_iprime(parsed: dict[str, Any]) -> bool | None:
    """Whether ``NCURR = 1`` (the ``I'(s)`` derivative profile); ``None`` when absent."""
    value = _scalars(parsed).get("NCURR")
    if value is None:
        return None
    return indata_int(parsed, "NCURR") == 1
