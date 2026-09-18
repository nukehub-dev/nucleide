"""OpenMC statepoint tally bridge (pure Python; caller-side OpenMC).

Statepoint files are HDF5 and are never read in Rust. This module drives the
caller's OpenMC Python API (call sequence pinned against OpenMC 0.16.0, the
validation-container build) and returns plain nested lists that the Rust
analytics consume directly: per-group ``flux`` vectors feed
:func:`nucleide.damage.nrt_dpa` / ``arc_dpa`` / ``gas_appm`` (plus
``he_dpa_ratio`` and the ``unfold`` iterators, which take the same
caller-supplied group structure), with ``energy_bounds_mev`` as the ``bounds``
argument. Clearance screening consumes derived inventories, not flux, so the
bridge stops at the spectral arrays by design.

Export route (primary): ``openmc.StatePoint(path, autolink=False)`` (explicit
``close()`` in a ``finally`` — 0.16.0 documents no context manager), then
``statepoint.tallies[tally_id]``, then ``tally.get_reshaped_data(value=...,
expand_dims=True)`` for ``"mean"`` and ``"std_dev"``. Energy edges come from
the tally's ``EnergyFilter.bins`` (eV, divided by 1e6 to MeV); mesh geometry
from ``MeshFilter.mesh`` (``RegularMesh`` only: ``dimension`` + ``lower_left``
+ ``upper_right``, centimetres). ``autolink=False`` because tally results live
in the statepoint itself — no ``summary.h5`` hunt, no surprise link failure.

Array contract: ``flux`` is always ``[voxel][group]`` (mesh voxels flattened
in C-order over the expanded ``(nx, ny, nz)`` axes; a cell tally is one
voxel), ``rel_err`` the matching ``std_dev / mean`` table (0.0 where the mean
is 0), ``energy_bounds_mev`` the ``G+1`` MeV edges (``[]`` when the tally has
no energy filter — a single total group), ``mesh_shape`` / ``mesh_bounds``
only for ``kind="mesh"``. ``cell_flux`` unwraps the single voxel for the
damage folds; ``voxel_flux`` selects one mesh voxel.

Export route (fallback): :func:`tally_to_csv` / :func:`tally_from_csv` move
the same dict through a commented CSV with stock :mod:`csv` only, so a
machine with OpenMC can export and any machine can read back bit-identically
— no OpenMC, no pandas, no HDF5 on the reading side.

Failure stance: a missing file, a missing OpenMC, a missing tally/score, a
multi-nuclide tally, a multi-bin non-energy filter (caller slices with
``get_slice`` first), a second mesh/energy filter, or a non-regular mesh are
all loud errors naming the remedy — never silent empties. OpenMC stays a
caller-side import (lazy, inside the functions that need it): the base
package never depends on it.

Synthetic-statepoint protocol (what the unit tests build in-test): a
statepoint exposes ``.tallies`` (dict id -> tally); a tally exposes
``.scores`` / ``.nuclides`` / ``.filters`` plus ``.get_reshaped_data(value=,
expand_dims=)`` returning an object with ``.shape`` / ``.tolist()``; an
energy filter is exactly ``type(f).__name__ == "EnergyFilter"`` with float
``.bins`` (eV); a mesh filter exposes ``.mesh`` whose type name is
``"RegularMesh"`` with ``.dimension`` / ``.lower_left`` / ``.upper_right``;
every other filter must report a single bin via ``.num_bins`` (or a
single-entry ``.bins``).
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path
from typing import Any

__all__ = [
    "openmc_version",
    "tally_arrays",
    "read_statepoint",
    "cell_flux",
    "voxel_flux",
    "tally_to_csv",
    "tally_from_csv",
]

#: OpenMC release the export call sequence is pinned against (validation
#: container build; the duck-typed extraction also works on 0.15.x, which
#: already carries ``get_reshaped_data(expand_dims=...)``).
VERSION_CONTRACT = "0.16.0"

_EV_PER_MEV = 1.0e6
_CSV_FORMAT = "nucleide-openmc-tally-v1"
_CSV_HEADER = ["voxel", "group", "lower_mev", "upper_mev", "flux", "rel_err"]


def openmc_version() -> str | None:
    """Installed caller-side OpenMC version, or ``None`` when not importable."""
    try:
        import openmc as om
    except ImportError:
        return None
    return str(om.__version__)


def _require_openmc() -> Any:
    """Import caller-side OpenMC or raise a loud error with the remedy."""
    try:
        import openmc as om
    except ImportError as exc:
        raise ImportError(
            "nucleide.openmc needs the caller-side OpenMC Python API, which is not "
            "installed in this environment (the base nucleide package never depends "
            f"on it; the call sequence is pinned against OpenMC {VERSION_CONTRACT}, "
            "the validation-container build). Install OpenMC where the statepoint "
            "lives and call read_statepoint() there, or export the tally arrays to "
            "CSV with tally_to_csv() and read them here with tally_from_csv() "
            "(dependency-free)."
        ) from exc
    return om


def tally_arrays(statepoint: Any, tally_id: int, score: str = "flux") -> dict[str, Any]:
    """Extract one tally as plain ``[voxel][group]`` arrays.

    ``statepoint`` is an ``openmc.StatePoint`` (or a synthetic object built
    to the module protocol). ``score`` defaults to ``"flux"``; any other
    score the tally carries works the same way. The returned dict carries
    ``tally_id``, ``score``, ``kind`` (``"mesh"``/``"cell"``), ``n_voxels``,
    ``n_groups``, ``energy_bounds_mev`` (``G+1`` edges, or ``[]`` unbinned),
    ``flux`` / ``rel_err`` (``[V][G]``), ``mesh_shape`` / ``mesh_bounds``
    (mesh only), and ``openmc_version`` (``None`` when OpenMC is not
    importable, e.g. synthetic objects in a bare environment).
    """
    try:
        tallies = statepoint.tallies
    except AttributeError as exc:
        raise TypeError(
            "tally_arrays needs an object exposing .tallies (openmc.StatePoint, "
            "or a synthetic statepoint built to the documented protocol)"
        ) from exc
    if not isinstance(tallies, dict) or tally_id not in tallies:
        available = sorted(tallies) if isinstance(tallies, dict) else tallies
        raise KeyError(f"tally {tally_id!r} not in statepoint; available tally ids: {available}")
    tally = tallies[tally_id]
    scores = [str(s) for s in tally.scores]
    if score not in scores:
        raise ValueError(
            f"tally {tally_id} has no score {score!r}; available scores: {scores} "
            "(score the tally for it, or pass one of the listed scores)"
        )
    score_idx = scores.index(score)
    nuclides = list(tally.nuclides)
    if len(nuclides) > 1:
        raise ValueError(
            f"tally {tally_id} spans {len(nuclides)} nuclides {nuclides}; slice to "
            "one nuclide caller-side (openmc get_slice/summation) first"
        )
    if not nuclides:
        raise ValueError(f"tally {tally_id} carries no nuclide bins; nothing to extract")
    mesh_f, energy_f = _classify_filters(tally.filters, tally_id)
    if energy_f is None:
        energy_bounds: list[float] = []
        n_groups = 1
    else:
        energy_bounds = _energy_bounds_mev(energy_f, tally_id)
        n_groups = len(energy_bounds) - 1
    if mesh_f is None:
        kind = "cell"
        space_dims: list[int] = []
        mesh_shape: list[int] | None = None
        mesh_bounds: dict[str, list[float]] | None = None
    else:
        kind = "mesh"
        mesh_shape, mesh_bounds = _mesh_geometry(mesh_f, tally_id)
        space_dims = list(mesh_shape)
    n_voxels = 1
    for dim in space_dims:
        n_voxels *= dim
    flux, rel_err = _extract_values(
        tally, tally_id, score_idx, space_dims, n_groups, _singleton_count(tally.filters)
    )
    return {
        "tally_id": tally_id,
        "score": score,
        "kind": kind,
        "n_voxels": n_voxels,
        "n_groups": n_groups,
        "energy_bounds_mev": energy_bounds,
        "flux": flux,
        "rel_err": rel_err,
        "mesh_shape": mesh_shape,
        "mesh_bounds": mesh_bounds,
        "openmc_version": openmc_version(),
    }


def read_statepoint(path: str | Path, tally_id: int, score: str = "flux") -> dict[str, Any]:
    """Open a statepoint file with caller-side OpenMC and extract one tally.

    Same returned dict as :func:`tally_arrays`. The file must exist (loud
    ``FileNotFoundError`` otherwise) and OpenMC must be importable (loud
    ``ImportError`` with the CSV-fallback remedy otherwise). The statepoint
    is opened with ``autolink=False`` and always closed.
    """
    state_path = Path(path)
    if not state_path.is_file():
        raise FileNotFoundError(
            f"OpenMC statepoint not found: {state_path} (run OpenMC with tallies "
            "enabled, then point read_statepoint at the statepoint.N.h5 file)"
        )
    om = _require_openmc()
    statepoint = om.StatePoint(str(state_path), autolink=False)
    try:
        return tally_arrays(statepoint, tally_id, score=score)
    finally:
        close = getattr(statepoint, "close", None)
        if callable(close):
            close()


def cell_flux(arrays: dict[str, Any]) -> list[float]:
    """Unwrap the single-voxel group flux for the damage folds.

    Returns ``arrays["flux"][0]`` — the ``flux`` vector
    :func:`nucleide.damage.nrt_dpa` takes alongside
    ``arrays["energy_bounds_mev"]`` as ``bounds``. Tallies with more than
    one voxel are a loud error (select one with :func:`voxel_flux`).
    """
    flux = _flux_table(arrays)
    if len(flux) != 1:
        raise ValueError(
            f"cell_flux needs a single-voxel (cell) tally, got {len(flux)} voxels; "
            "use voxel_flux(arrays, voxel) to select one mesh voxel"
        )
    return [float(v) for v in flux[0]]


def voxel_flux(arrays: dict[str, Any], voxel: int) -> list[float]:
    """Return the group-flux row of one mesh voxel (C-order voxel index)."""
    if isinstance(voxel, bool) or not isinstance(voxel, int):
        raise TypeError(f"voxel index must be an int, got {voxel!r}")
    flux = _flux_table(arrays)
    if not 0 <= voxel < len(flux):
        raise ValueError(f"voxel {voxel} out of range for {len(flux)} voxels")
    return [float(v) for v in flux[voxel]]


def tally_to_csv(arrays: dict[str, Any], path: str | Path) -> None:
    """Write tally arrays to the commented CSV interchange (stock stdlib).

    The ``#``-comment header records the tally metadata (id, score, kind,
    shapes, mesh and energy bounds, OpenMC version); the body holds one row
    per ``(voxel, group)`` with the group bounds inline. The file reads back
    bit-identically through :func:`tally_from_csv`.
    """
    flux = _flux_table(arrays)
    rel_err = arrays.get("rel_err")
    n_voxels = len(flux)
    n_groups = len(flux[0]) if n_voxels else 0
    if (
        not isinstance(rel_err, list)
        or len(rel_err) != n_voxels
        or any(len(row) != n_groups for row in rel_err)
    ):
        raise ValueError("tally_to_csv needs the dict from tally_arrays/tally_from_csv")
    energy_bounds = arrays.get("energy_bounds_mev")
    if not isinstance(energy_bounds, list):
        raise ValueError("tally_to_csv needs the dict from tally_arrays/tally_from_csv")
    if energy_bounds and len(energy_bounds) != n_groups + 1:
        raise ValueError(
            f"energy bounds length {len(energy_bounds)} != groups + 1 ({n_groups + 1})"
        )
    kind = arrays.get("kind")
    mesh_shape = arrays.get("mesh_shape")
    mesh_bounds = arrays.get("mesh_bounds")
    if kind == "cell":
        if n_voxels != 1:
            raise ValueError(f"a cell tally holds one voxel, got {n_voxels}")
    elif kind == "mesh":
        if not isinstance(mesh_shape, list) or not isinstance(mesh_bounds, dict):
            raise ValueError("a mesh tally needs mesh_shape and mesh_bounds")
        expect = 1
        for dim in mesh_shape:
            expect *= int(dim)
        if expect != n_voxels:
            raise ValueError(f"mesh_shape {mesh_shape} holds {expect} voxels, flux has {n_voxels}")
    else:
        raise ValueError(f"unknown tally kind {kind!r}; want 'mesh' or 'cell'")
    version = arrays.get("openmc_version")
    lines: dict[str, str] = {
        "format": _CSV_FORMAT,
        "tally_id": str(arrays.get("tally_id")),
        "score": str(arrays.get("score")),
        "kind": str(kind),
        "n_voxels": str(n_voxels),
        "n_groups": str(n_groups),
        "mesh_shape": "-" if mesh_shape is None else ",".join(str(d) for d in mesh_shape),
        "energy_bounds": "-" if not energy_bounds else ",".join(repr(b) for b in energy_bounds),
        "openmc_version": "-" if version is None else str(version),
    }
    if isinstance(mesh_bounds, dict):
        for axis in ("x", "y", "z"):
            bounds = mesh_bounds.get(axis)
            lines[f"mesh_{axis}_bounds"] = "-" if not bounds else ",".join(repr(b) for b in bounds)
    else:
        for axis in ("x", "y", "z"):
            lines[f"mesh_{axis}_bounds"] = "-"
    with open(Path(path), "w", newline="", encoding="utf-8") as handle:
        for key, value in lines.items():
            handle.write(f"# {key} = {value}\n")
        writer = csv.writer(handle)
        writer.writerow(_CSV_HEADER)
        for voxel in range(n_voxels):
            for group in range(n_groups):
                if energy_bounds:
                    lower = repr(energy_bounds[group])
                    upper = repr(energy_bounds[group + 1])
                else:
                    lower = upper = ""
                writer.writerow(
                    [
                        voxel,
                        group,
                        lower,
                        upper,
                        repr(float(flux[voxel][group])),
                        repr(float(rel_err[voxel][group])),
                    ]
                )


def tally_from_csv(path: str | Path) -> dict[str, Any]:
    """Read back the commented CSV interchange into the tally-arrays dict."""
    state_path = Path(path)
    if not state_path.is_file():
        raise FileNotFoundError(f"tally CSV not found: {state_path}")
    with open(state_path, newline="", encoding="utf-8") as handle:
        raw = handle.read().splitlines()
    meta: dict[str, str] = {}
    header: list[str] | None = None
    rows: list[list[str]] = []
    for line in raw:
        if not line.strip():
            continue
        if line.startswith("#"):
            key, equals, value = line[1:].partition("=")
            if equals:
                meta[key.strip()] = value.strip()
        elif header is None:
            header = next(csv.reader([line]))
        else:
            rows.append(next(csv.reader([line])))
    if meta.get("format") != _CSV_FORMAT:
        raise ValueError(
            f"{state_path} is not a {_CSV_FORMAT} file "
            "(export one with tally_to_csv, or read the statepoint with read_statepoint)"
        )
    for key in ("tally_id", "score", "kind", "n_voxels", "n_groups"):
        if key not in meta:
            raise ValueError(f"{state_path}: tally CSV header is missing {key!r}")
    if header != _CSV_HEADER:
        raise ValueError(f"{state_path}: tally CSV header row must be {_CSV_HEADER}")
    try:
        n_voxels = int(meta["n_voxels"])
        n_groups = int(meta["n_groups"])
    except ValueError as exc:
        raise ValueError(f"{state_path}: tally CSV counts are not ints") from exc
    if n_voxels < 1 or n_groups < 1:
        raise ValueError(f"{state_path}: tally CSV needs positive counts")
    if len(rows) != n_voxels * n_groups:
        raise ValueError(
            f"{state_path}: tally CSV holds {len(rows)} rows, "
            f"want {n_voxels * n_groups} ({n_voxels} voxels x {n_groups} groups)"
        )
    energy_bounds = _parse_float_list(meta.get("energy_bounds", "-"), state_path, "energy_bounds")
    if energy_bounds and len(energy_bounds) != n_groups + 1:
        raise ValueError(f"{state_path}: energy bounds length != groups + 1")
    flux: list[list[float]] = [[0.0] * n_groups for _ in range(n_voxels)]
    rel_err: list[list[float]] = [[0.0] * n_groups for _ in range(n_voxels)]
    for position, row in enumerate(rows):
        if len(row) != 6:
            raise ValueError(f"{state_path}: row {position} has {len(row)} fields, want 6")
        try:
            voxel, group = int(row[0]), int(row[1])
            flux_value, rel_value = float(row[4]), float(row[5])
        except ValueError as exc:
            raise ValueError(f"{state_path}: row {position} has non-numeric entries") from exc
        if voxel != position // n_groups or group != position % n_groups:
            raise ValueError(
                f"{state_path}: row {position} breaks (voxel, group) order "
                "(rows must run voxel-major, groups ascending)"
            )
        if energy_bounds:
            if row[2] == "" or row[3] == "":
                raise ValueError(f"{state_path}: row {position} is missing group bounds")
            if (float(row[2]), float(row[3])) != (
                energy_bounds[group],
                energy_bounds[group + 1],
            ):
                raise ValueError(f"{state_path}: row {position} bounds disagree with header")
        elif row[2] != "" or row[3] != "":
            raise ValueError(f"{state_path}: row {position} carries bounds for an unbinned tally")
        flux[voxel][group] = flux_value
        rel_err[voxel][group] = rel_value
    kind = meta["kind"]
    if kind == "mesh":
        mesh_shape = [
            int(d) for d in _parse_float_list(meta.get("mesh_shape", "-"), state_path, "mesh_shape")
        ]
        mesh_bounds = {
            axis: _parse_float_list(
                meta.get(f"mesh_{axis}_bounds", "-"), state_path, f"mesh_{axis}_bounds"
            )
            for axis in ("x", "y", "z")
        }
        if not mesh_shape or any(not b for b in mesh_bounds.values()):
            raise ValueError(f"{state_path}: mesh tally CSV is missing mesh metadata")
    elif kind == "cell":
        mesh_shape = None
        mesh_bounds = None
    else:
        raise ValueError(f"{state_path}: unknown tally kind {kind!r}")
    version = meta.get("openmc_version", "-")
    return {
        "tally_id": int(meta["tally_id"]),
        "score": meta["score"],
        "kind": kind,
        "n_voxels": n_voxels,
        "n_groups": n_groups,
        "energy_bounds_mev": energy_bounds,
        "flux": flux,
        "rel_err": rel_err,
        "mesh_shape": mesh_shape,
        "mesh_bounds": mesh_bounds,
        "openmc_version": None if version == "-" else version,
    }


def _classify_filters(filters: Any, tally_id: int) -> tuple[Any | None, Any | None]:
    """Split tally filters into (mesh filter|None, energy filter|None).

    Every other filter must select a single bin; anything wider (or
    unparsable) is a loud error telling the caller to slice first.
    """
    mesh_f: Any | None = None
    energy_f: Any | None = None
    for entry in filters:
        name = type(entry).__name__
        mesh = getattr(entry, "mesh", None)
        if mesh is not None and hasattr(mesh, "dimension"):
            if mesh_f is not None:
                raise ValueError(
                    f"tally {tally_id} has two mesh filters; slice to one mesh "
                    "caller-side (openmc get_slice) first"
                )
            mesh_f = entry
            continue
        if name == "EnergyFilter":
            if energy_f is not None:
                raise ValueError(
                    f"tally {tally_id} has two energy filters; slice to one "
                    "caller-side (openmc get_slice) first"
                )
            energy_f = entry
            continue
        if _bin_count(entry, name, tally_id) == 1:
            continue
        raise ValueError(
            f"tally {tally_id} has a multi-bin {name} filter; slice to one bin "
            "caller-side (openmc get_slice) first — the bridge maps one spatial "
            "selector and one energy axis onto the fold inputs"
        )
    return mesh_f, energy_f


def _bin_count(entry: Any, name: str, tally_id: int) -> int:
    """Number of bins of a non-mesh, non-energy filter (loud when unknown)."""
    count = getattr(entry, "num_bins", None)
    if isinstance(count, int) and not isinstance(count, bool):
        return count
    bins = getattr(entry, "bins", None)
    if bins is not None:
        return len(list(bins))
    raise ValueError(
        f"tally {tally_id}: cannot determine the bin count of filter {name!r}; "
        "slice caller-side (openmc get_slice) to a single bin first"
    )


def _singleton_count(filters: Any) -> int:
    """Number of single-bin (non-mesh, non-energy) filters on the tally."""
    count = 0
    for entry in filters:
        mesh = getattr(entry, "mesh", None)
        if mesh is not None and hasattr(mesh, "dimension"):
            continue
        if type(entry).__name__ == "EnergyFilter":
            continue
        count += 1
    return count


def _energy_bounds_mev(entry: Any, tally_id: int) -> list[float]:
    """EnergyFilter bins (eV) converted to strictly increasing MeV edges."""
    try:
        raw = [float(b) for b in entry.bins]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"tally {tally_id}: energy filter bins are not numeric") from exc
    bounds = [b / _EV_PER_MEV for b in raw]
    if len(bounds) < 2 or any(high <= low for low, high in zip(bounds, bounds[1:], strict=False)):
        raise ValueError(
            f"tally {tally_id}: energy filter bins must be strictly increasing "
            f"(got {len(bounds)} edges)"
        )
    return bounds


def _mesh_geometry(entry: Any, tally_id: int) -> tuple[list[int], dict[str, list[float]]]:
    """RegularMesh shape and per-axis centimetre bounds (loud otherwise)."""
    mesh = entry.mesh
    if type(mesh).__name__ != "RegularMesh":
        raise ValueError(
            f"tally {tally_id}: only RegularMesh tallies are bridged "
            f"(got {type(mesh).__name__}); use a regular mesh or a cell tally"
        )
    try:
        shape = [int(d) for d in mesh.dimension]
        lower = [float(v) for v in mesh.lower_left]
        upper = [float(v) for v in mesh.upper_right]
    except (TypeError, ValueError, AttributeError) as exc:
        raise ValueError(f"tally {tally_id}: cannot read RegularMesh dimension/bounds") from exc
    if len(shape) != 3 or len(lower) != 3 or len(upper) != 3:
        raise ValueError(f"tally {tally_id}: RegularMesh geometry must be 3D")
    if any(d < 1 for d in shape) or any(
        high <= low for low, high in zip(lower, upper, strict=True)
    ):
        raise ValueError(f"tally {tally_id}: RegularMesh has a degenerate dimension")
    bounds = {
        axis: [low + i * (high - low) / n for i in range(n + 1)]
        for axis, low, high, n in zip(("x", "y", "z"), lower, upper, shape, strict=True)
    }
    return shape, bounds


def _extract_values(
    tally: Any,
    tally_id: int,
    score_idx: int,
    space_dims: list[int],
    n_groups: int,
    n_singletons: int,
) -> tuple[list[list[float]], list[list[float]]]:
    """Index the reshaped mean/std_dev tables into ``[V][G]`` plain lists."""
    try:
        means = tally.get_reshaped_data(value="mean", expand_dims=True)
        stds = tally.get_reshaped_data(value="std_dev", expand_dims=True)
    except AttributeError as exc:
        raise TypeError(
            f"tally {tally_id} does not expose get_reshaped_data(value=, expand_dims=) "
            f"(needs OpenMC >= 0.13.3; contract {VERSION_CONTRACT})"
        ) from exc
    try:
        mean_shape = tuple(int(d) for d in means.shape)
        std_shape = tuple(int(d) for d in stds.shape)
    except (AttributeError, TypeError, ValueError) as exc:
        raise TypeError(f"tally {tally_id}: reshaped data has no readable .shape") from exc
    roles = ["space"] * len(space_dims)
    sizes = list(space_dims)
    if n_groups > 1 or _has_energy_axis(tally):
        roles.append("energy")
        sizes.append(n_groups)
    roles.extend(["single"] * n_singletons)
    sizes.extend([1] * n_singletons)
    roles.extend(["nuclide", "score"])
    sizes.extend([1, _score_total(tally)])
    if mean_shape != tuple(sizes) or std_shape != tuple(sizes):
        raise RuntimeError(
            f"tally {tally_id}: reshaped data shape {mean_shape} disagrees with the "
            f"filter layout {tuple(sizes)} (OpenMC contract {VERSION_CONTRACT})"
        )
    mean_rows = means.tolist()
    std_rows = stds.tolist()
    n_voxels = 1
    for dim in space_dims:
        n_voxels *= dim
    flux: list[list[float]] = []
    rel_err: list[list[float]] = []
    for spatial in itertools.product(*(range(dim) for dim in space_dims)):
        flux_row: list[float] = []
        rel_row: list[float] = []
        for group in range(n_groups):
            index = _locate(roles, spatial, group, score_idx)
            mean = float(_at(mean_rows, index, tally_id))
            std = float(_at(std_rows, index, tally_id))
            flux_row.append(mean)
            rel_row.append(std / mean if mean != 0.0 else 0.0)
        flux.append(flux_row)
        rel_err.append(rel_row)
    return flux, rel_err


def _has_energy_axis(tally: Any) -> bool:
    """Whether the tally carries an EnergyFilter (single-group axis kept)."""
    return any(type(entry).__name__ == "EnergyFilter" for entry in tally.filters)


def _score_total(tally: Any) -> int:
    """Number of score bins on the tally (loud when unreadable)."""
    try:
        return len(list(tally.scores))
    except TypeError as exc:
        raise ValueError("tally scores are not readable") from exc


def _locate(
    roles: list[str], spatial: tuple[int, ...], group: int, score_idx: int
) -> tuple[int, ...]:
    """Build the full index tuple for one (voxel, group) entry."""
    index: list[int] = []
    axis = 0
    for role in roles:
        if role == "space":
            index.append(spatial[axis])
            axis += 1
        elif role == "energy":
            index.append(group)
        elif role == "score":
            index.append(score_idx)
        else:
            index.append(0)
    return tuple(index)


def _at(nested: Any, index: tuple[int, ...], tally_id: int) -> Any:
    """Index nested lists (loud RuntimeError instead of a bare IndexError)."""
    value = nested
    try:
        for axis in index:
            value = value[axis]
    except (IndexError, TypeError) as exc:
        raise RuntimeError(f"tally {tally_id}: reshaped data indexing failed at {index}") from exc
    return value


def _flux_table(arrays: dict[str, Any]) -> list[list[Any]]:
    """Read and rectangularity-check the ``flux`` table of an arrays dict."""
    try:
        flux = arrays["flux"]
    except (KeyError, TypeError) as exc:
        raise ValueError("needs the dict returned by tally_arrays/tally_from_csv") from exc
    if not isinstance(flux, list) or not flux:
        raise ValueError("needs the dict returned by tally_arrays/tally_from_csv")
    width = len(flux[0])
    if width < 1 or any(not isinstance(row, list) or len(row) != width for row in flux):
        raise ValueError("needs the dict returned by tally_arrays/tally_from_csv")
    return flux


def _parse_float_list(text: str, path: Path, key: str) -> list[float]:
    """Parse a comma-separated ``#``-comment list (``-`` means absent)."""
    if text.strip() == "-" or text.strip() == "":
        return []
    try:
        return [float(part) for part in text.split(",")]
    except ValueError as exc:
        raise ValueError(f"{path}: tally CSV {key!r} is not a float list") from exc
