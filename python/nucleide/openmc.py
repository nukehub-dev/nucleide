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

Normalization and source handoff: raw extraction is intentionally not
relabelled as physical flux. OpenMC 0.16 fixed-source tallies are normalized
by the simulated particle count and scaled by the configured total source
strength, so a default unit-strength source is per source particle while an
absolute-rate source is already scaled. :func:`normalize_flux` requires the
caller to declare the volume state and rate state, and rejects double
scaling. :func:`mesh_source` builds an optional ``openmc.MeshSource`` from an
integrated photon-source dict, with explicit C-order/Fortran-order mapping
and per-element histogram densities; :func:`export_cooling_batch` and the
manifest helpers add deterministic files and offline hash verification.

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
import hashlib
import itertools
import json
import math
import os
import struct
import tempfile
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, NamedTuple

from nucleide import _internal as _core

__all__ = [
    "openmc_version",
    "tally_arrays",
    "read_statepoint",
    "cell_flux",
    "voxel_flux",
    "tally_to_csv",
    "tally_from_csv",
    "normalize_flux",
    "source_manifest",
    "write_manifest",
    "read_manifest",
    "verify_manifest",
    "export_cooling_batch",
    "mesh_source",
    "write_mesh_source_xml",
]

#: OpenMC release the export call sequence is pinned against (validation
#: container build; the duck-typed extraction also works on 0.15.x, which
#: already carries ``get_reshaped_data(expand_dims=...)``).
VERSION_CONTRACT = "0.16.0"

#: Input normalization states for :func:`normalize_flux`. OpenMC's default
#: track-length ``flux`` tally is volume-integrated per source particle; a
#: random-ray tally run with ``volume_normalized_flux_tallies`` is
#: volume-averaged per source particle. The bridge cannot infer which one a
#: caller produced, so conversion requires an explicit declaration.
NORMALIZATION_UNKNOWN = "unknown"
NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE = "per-source-volume-integrated"
NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED = "rate-scaled-volume-integrated"
NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE = "per-source-volume-averaged"
NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED = "rate-scaled-volume-averaged"
NORMALIZATION_PHYSICAL = "physical"

UNITS_VOLUME_INTEGRATED_PER_SOURCE = "cm per source particle"
UNITS_VOLUME_INTEGRATED_RATE_SCALED = "cm/s (source-rate scaled)"
UNITS_VOLUME_AVERAGED_PER_SOURCE = "1/(cm2 source-particle)"
UNITS_VOLUME_AVERAGED_RATE_SCALED = "1/(cm2 s) (source-rate scaled)"
UNITS_PHYSICAL = "1/(cm2 s)"

_MANIFEST_SCHEMA = "nucleide-source-manifest-v1"
_BATCH_SCHEMA = "nucleide-source-batch-v1"
_VOXEL_ORDER = "c-order-(i,j,k)"
_ASSUMPTIONS = [
    "uniform-position-within-voxel",
    "uniform-density-per-ev-within-group",
]
_NORMALIZATION_STATES = {
    NORMALIZATION_UNKNOWN,
    NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
    NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED,
    NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE,
    NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED,
    NORMALIZATION_PHYSICAL,
}

_EV_PER_MEV = 1.0e6
_CSV_FORMAT = "nucleide-openmc-tally-v2"
_CSV_FORMATS = {_CSV_FORMAT, "nucleide-openmc-tally-v1"}
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
        # Raw extraction is per source particle, but whether a flux tally is
        # volume-integrated (default track-length) or already volume-averaged
        # (random ray) is caller-configured. Never guess: normalize_flux
        # requires an explicit declaration.
        "normalization": NORMALIZATION_UNKNOWN,
        "units": "unknown",
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
    normalization = _csv_normalization(arrays)
    units = str(arrays.get("units", _default_units(normalization)))
    rate = arrays.get("source_rate_s")
    if rate is not None:
        rate = float(rate)
        if not math.isfinite(rate) or rate < 0.0:
            raise ValueError("tally_to_csv needs a finite nonnegative source_rate_s")
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
        "normalization": normalization,
        "units": units,
    }
    if rate is not None:
        lines["source_rate_s"] = repr(rate)
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
    file_format = meta.get("format")
    if file_format not in _CSV_FORMATS:
        raise ValueError(
            f"{state_path} is not a {_CSV_FORMAT} (or legacy v1) file "
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
    if file_format == _CSV_FORMAT:
        normalization = meta.get("normalization", NORMALIZATION_UNKNOWN)
        if normalization not in _NORMALIZATION_STATES:
            raise ValueError(
                f"{state_path}: unknown tally normalization {normalization!r}; "
                f"want one of {sorted(_NORMALIZATION_STATES)}"
            )
        units = meta.get("units", _default_units(normalization))
        rate_text = meta.get("source_rate_s")
        source_rate_s = None
        if rate_text not in (None, "-"):
            try:
                source_rate_s = float(rate_text)
            except ValueError as exc:
                raise ValueError(f"{state_path}: source_rate_s is not a float") from exc
            if not math.isfinite(source_rate_s) or source_rate_s < 0.0:
                raise ValueError(f"{state_path}: source_rate_s must be finite and nonnegative")
    else:
        # Legacy v1 carries no normalization metadata. Raw data is not
        # physical flux; callers must declare its state before conversion.
        normalization = NORMALIZATION_UNKNOWN
        units = "unknown"
        source_rate_s = None
    out: dict[str, Any] = {
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
        "normalization": normalization,
        "units": units,
    }
    if source_rate_s is not None:
        out["source_rate_s"] = source_rate_s
    return out


def normalize_flux(
    arrays: dict[str, Any],
    source_rate_s: float,
    volumes_cm3: Sequence[float],
    *,
    input_normalization: str | None = None,
) -> dict[str, Any]:
    """Convert raw OpenMC tally arrays into physical group-integrated flux.

    The conversion is additive: ``arrays`` is not mutated and its raw
    extraction semantics are preserved for every caller that does not opt in.
    ``input_normalization`` declares both axes: volume-integrated
    (:data:`NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE`, OpenMC's default
    track-length flux tally) or volume-averaged
    (:data:`NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE`, random-ray output run
    with ``volume_normalized_flux_tallies``), and per-source-particle or
    rate-scaled. OpenMC 0.16 fixed-source tallies are also scaled by the
    configured total source strength, so the ``*_RATE_SCALED`` variants
    declare that the source rate is already inside the tally and must not be
    multiplied again. The declaration may be omitted only when the arrays dict
    already carries one of those states. Legacy CSV v1 data and fresh
    :func:`tally_arrays` output are ``unknown`` and require a caller
    declaration; :data:`NORMALIZATION_PHYSICAL` is rejected as a double
    normalization.

    For per-source volume-integrated input
    ``phi[v,g] = raw[v,g] * rate / volume[v]``; the rate-scaled variant
    computes only ``raw / volume``. Volume-averaged input applies only the
    declared rate scale (or nothing when already rate-scaled). Output units
    are ``1/(cm2 s)``. Relative errors are copied unchanged (deterministic
    scaling preserves them); source-rate and volume uncertainty, cross-group
    covariance and representability limits are outside this conversion.
    """
    flux = _flux_table(arrays)
    rel_err = arrays.get("rel_err")
    n_voxels, n_groups = len(flux), len(flux[0])
    if (
        not isinstance(rel_err, list)
        or len(rel_err) != n_voxels
        or any(not isinstance(row, list) or len(row) != n_groups for row in rel_err)
    ):
        raise ValueError(
            "normalize_flux needs the dict from tally_arrays/tally_from_csv "
            "(flux and matching rel_err tables)"
        )
    rate = float(source_rate_s)
    if not math.isfinite(rate) or rate < 0.0:
        raise ValueError(f"source_rate_s must be finite and nonnegative, got {source_rate_s!r}")
    volumes = [float(v) for v in volumes_cm3]
    if len(volumes) != n_voxels:
        raise ValueError(
            f"volumes_cm3 has {len(volumes)} entries for {n_voxels} voxels; "
            "pass one positive volume per voxel"
        )
    recorded = str(arrays.get("normalization", NORMALIZATION_UNKNOWN))
    if input_normalization is None:
        if recorded == NORMALIZATION_UNKNOWN:
            raise ValueError(
                "raw tally normalization is unknown; declare input_normalization="
                f"{NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE!r} for OpenMC's default "
                "track-length flux tally with normalized source weights, "
                f"{NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED!r} when the configured "
                "source strengths already sum to the physical rate, or "
                f"{NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE!r}/"
                f"{NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED!r} for random-ray output"
            )
        declared = recorded
    else:
        declared = str(input_normalization)
        if recorded not in (NORMALIZATION_UNKNOWN, declared):
            raise ValueError(
                f"declared input_normalization {declared!r} disagrees with the arrays "
                f"metadata {recorded!r}; refusing to guess which normalization applies"
            )
    if declared == NORMALIZATION_PHYSICAL:
        raise ValueError(
            "these arrays are already physical flux; refusing to normalize twice "
            "(use the input normalization state from before the first conversion)"
        )
    if declared == NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE:
        mode = "volume-integrated"
    elif declared == NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED:
        mode = "volume-integrated-rate-scaled"
    elif declared == NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE:
        mode = "volume-averaged"
    elif declared == NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED:
        mode = "volume-averaged-rate-scaled"
    else:
        raise ValueError(
            f"unsupported input_normalization {declared!r}; want "
            f"{NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE!r}, "
            f"{NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED!r}, "
            f"{NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE!r}, or "
            f"{NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED!r}"
        )
    units = UNITS_PHYSICAL
    converted = _core.r2s_normalize_flux(flux, rate, volumes, mode)
    out = dict(arrays)
    out["flux"] = converted
    out["rel_err"] = [list(row) for row in rel_err]
    out["normalization"] = NORMALIZATION_PHYSICAL
    out["units"] = units
    out["input_normalization"] = declared
    out["source_rate_s"] = rate
    return out


def mesh_source(source: dict[str, Any]) -> Any:
    """Build a caller-side OpenMC ``MeshSource`` from a validated source dict.

    ``source`` is the dict returned by ``nucleide.r2s.photon_source`` or
    ``photon_source_exported``: regular-mesh geometry, ascending eV group
    edges, caller voxel volumes and ``[voxel][group]`` integrated photons/s.
    OpenMC is imported lazily; the base package never depends on it.

    Each mesh element receives its own ``IndependentSource`` with a
    histogram energy density (integrated group strength divided by the group
    width, per eV) and the voxel's absolute photons/s as its strength. The
    joint voxel/group distribution is therefore carried by the element
    strengths and per-element spectra, never reconstructed from independent
    spatial and energy marginals. Voxel order is the bridge's C-order
    ``(i, j, k)``; the nested array passed to OpenMC is indexed ``[i][j][k]``
    so OpenMC's internal Fortran flattening preserves the same elements.
    OpenMC scales its fixed-source tallies by the sum of these strengths, so a
    later flux conversion uses a ``rate-scaled-*`` declaration and must not
    multiply by the same total rate again.

    Zero-source voxels are allowed inside a nonzero source: they receive a
    uniform zero-strength spectrum so OpenMC's sampling probability for them
    is exactly zero. An all-zero source is rejected by the shared validator.
    The regular-mesh cell volumes must equal ``voxel_volumes_cm3``; a
    mismatch is a loud error rather than a bounding-box approximation.
    """
    om = _require_openmc()
    args = _source_args(source)
    mesh_shape, lower, upper, volumes, edges, _cooling, strengths = args
    report = _source_report_from_args(args)
    if not report["valid"]:
        codes = ", ".join(
            issue["code"] for issue in report["issues"] if issue["severity"] == "error"
        )
        raise ValueError(f"cannot export an invalid photon source; errors: {codes}")
    total = float(source.get("total_photons_per_s", 0.0))
    if not math.isclose(total, _sum_table(strengths), rel_tol=1e-12, abs_tol=0.0):
        raise ValueError(
            "source dict total_photons_per_s disagrees with its strengths table; "
            "rebuild it with nucleide.r2s.photon_source"
        )
    if total <= 0.0:
        raise ValueError("an all-zero photon source cannot be exported for sampling")
    probabilities = _core.r2s_source_energy_probabilities(*args)
    widths = [edges[g + 1] - edges[g] for g in range(len(edges) - 1)]
    n_x, n_y, n_z = mesh_shape
    dx = (upper[0] - lower[0]) / n_x
    dy = (upper[1] - lower[1]) / n_y
    dz = (upper[2] - lower[2]) / n_z
    expected_volume = dx * dy * dz
    for voxel, volume in enumerate(volumes):
        if not math.isclose(volume, expected_volume, rel_tol=1e-9, abs_tol=0.0):
            raise ValueError(
                f"voxel {voxel} volume {volume!r} cm3 does not match the regular-mesh "
                f"cell volume {expected_volume!r} cm3; OpenMC export supports only a "
                "regular mesh, not a bounding-box approximation of another geometry"
            )
    mesh = om.RegularMesh()
    mesh.dimension = (n_x, n_y, n_z)
    mesh.lower_left = tuple(lower)
    mesh.upper_right = tuple(upper)
    nested: list[list[list[Any]]] = []
    for i in range(n_x):
        plane: list[list[Any]] = []
        for j in range(n_y):
            column: list[Any] = []
            for k in range(n_z):
                voxel = (i * n_y + j) * n_z + k
                strength = _sum_row(strengths[voxel])
                if strength > 0.0:
                    # OpenMC histogram distributions take a differential PDF
                    # per energy unit; group probabilities must be divided by
                    # the bin width, never passed as integrated bins.
                    density = [
                        probabilities[voxel][group] / widths[group] for group in range(len(widths))
                    ]
                    energy: Any = om.stats.Tabular(edges, density, interpolation="histogram")
                else:
                    energy = om.stats.Uniform(edges[0], edges[-1])
                column.append(
                    om.IndependentSource(
                        energy=energy,
                        strength=strength,
                        particle="photon",
                    )
                )
            plane.append(column)
        nested.append(plane)
    mesh_source_obj = om.MeshSource(mesh=mesh, sources=nested)
    if not math.isclose(float(mesh_source_obj.strength), total, rel_tol=1e-12, abs_tol=0.0):
        raise RuntimeError(
            "OpenMC MeshSource strength does not match the validated absolute total; "
            "refusing to export an inconsistent source"
        )
    return mesh_source_obj


def write_mesh_source_xml(source: dict[str, Any], path: str | Path) -> Path:
    """Write a self-contained OpenMC ``settings.xml`` for one photon source.

    The file contains the regular-mesh definition and the mesh source, so
    caller OpenMC can read it back with ``openmc.Settings.from_xml`` without
    a separate ``meshes.xml``. The parent directory must already exist;
    callers that need atomic replacement should use :func:`export_cooling_batch`
    or write through a temporary path.
    """
    om = _require_openmc()
    target = Path(path)
    if not target.parent.is_dir():
        raise FileNotFoundError(f"output directory does not exist: {target.parent}")
    settings = om.Settings()
    settings.source = [mesh_source(source)]
    settings.export_to_xml(target)
    return target


def source_manifest(
    *,
    source: dict[str, Any],
    input_files: Iterable[str | Path] = (),
    output_files: Iterable[str | Path] = (),
    external_solver: str | None = None,
    external_version: str | None = None,
    adapter: str = "openmc-mesh-source",
    adapter_version: str = "1",
    evaluated_library: str | None = None,
    validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a compact, versioned source reproducibility manifest.

    The manifest records identifiers and measured SHA-256 hashes for files
    that exist, the energy-grid identity, voxel order, cooling offset,
    absolute photons/s and the shared validation summary. It never invents
    hashes for evaluated data: ``evaluated_library`` is a caller-supplied
    identity only. Bulk arrays are referenced by grid hash, not duplicated.
    JSON serialization is deterministic through :func:`write_manifest`; no
    timestamps, environment variables or credentials are recorded.
    """
    args = _source_args(source)
    report = _source_report_from_args(args)
    if not report["valid"]:
        codes = ", ".join(
            issue["code"] for issue in report["issues"] if issue["severity"] == "error"
        )
        raise ValueError(f"cannot manifest an invalid photon source; errors: {codes}")
    grid_hash = _energy_grid_sha256(args[4])
    validation_summary = (
        _report_summary(report) if validation is None else _validation_summary(validation)
    )
    manifest: dict[str, Any] = {
        "schema": _MANIFEST_SCHEMA,
        "created_by": {"nucleide": _core.version()},
        "external": {
            "solver": external_solver,
            "solver_version": external_version,
            "solver_version_source": "caller-supplied" if external_version else None,
            "adapter": adapter,
            "adapter_version": adapter_version,
        },
        "source": {
            "kind": "photons",
            "strength_units": "photons/s",
            "normalization": "integrated-per-voxel-and-group",
            "cooling_s": args[5],
            "mesh_shape": args[0],
            "lower_left_cm": args[1],
            "upper_right_cm": args[2],
            "n_voxels": len(args[3]),
            "n_groups": len(args[4]) - 1,
            "voxel_order": _VOXEL_ORDER,
            "energy_grid_sha256": grid_hash,
            "total_photons_per_s": report["totals"]["total_photons_per_s"],
            "assumptions": list(_ASSUMPTIONS),
        },
        "evaluated_library": (
            None
            if evaluated_library is None
            else {"identity": evaluated_library, "hash": None, "status": "caller-supplied"}
        ),
        "inputs": [_file_record(path, require_exists=True) for path in input_files],
        "outputs": [_file_record(path, require_exists=False) for path in output_files],
        "validation": validation_summary,
    }
    return manifest


def write_manifest(manifest: dict[str, Any], path: str | Path) -> Path:
    """Write a manifest or batch manifest as deterministic, atomic JSON."""
    schema = manifest.get("schema")
    if schema not in {_MANIFEST_SCHEMA, _BATCH_SCHEMA}:
        raise ValueError(
            f"unknown manifest schema {schema!r}; want {_MANIFEST_SCHEMA!r} or {_BATCH_SCHEMA!r}"
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        text = json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False) + "\n"
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "manifest is not JSON-serializable (non-finite number or bad type)"
        ) from exc
    _write_text_atomic(target, text)
    return target


def read_manifest(path: str | Path) -> dict[str, Any]:
    """Read a source or batch manifest, rejecting unknown schemas."""
    manifest_path = Path(path)
    if not manifest_path.is_file():
        raise FileNotFoundError(f"manifest not found: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{manifest_path} is not valid JSON") from exc
    if not isinstance(manifest, dict):
        raise ValueError(f"{manifest_path}: manifest root must be a JSON object")
    schema = manifest.get("schema")
    if schema not in {_MANIFEST_SCHEMA, _BATCH_SCHEMA}:
        raise ValueError(
            f"{manifest_path}: unknown manifest schema {schema!r}; "
            f"want {_MANIFEST_SCHEMA!r} or {_BATCH_SCHEMA!r}"
        )
    return manifest


def verify_manifest(
    manifest_or_path: dict[str, Any] | str | Path,
    *,
    base_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify recorded file hashes and report mismatches, never repair them.

    ``manifest_or_path`` may be a previously built manifest dict or a path.
    Relative paths resolve against the manifest's directory (path input) or
    ``base_dir`` (dict input; default ``"."``). Records without a measured
    hash (explicit unknown/caller-supplied provenance) are not treated as
    verified files; missing files and hash mismatches make ``ok`` false.
    """
    if isinstance(manifest_or_path, (str, Path)):
        manifest_path = Path(manifest_or_path)
        manifest = read_manifest(manifest_path)
        base = manifest_path.parent
    else:
        manifest = manifest_or_path
        base = Path(base_dir if base_dir is not None else ".")
    if manifest.get("schema") not in {_MANIFEST_SCHEMA, _BATCH_SCHEMA}:
        raise ValueError(f"unknown manifest schema {manifest.get('schema')!r}")
    checked = 0
    mismatches: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for section in ("inputs", "outputs"):
        records = manifest.get(section, [])
        if not isinstance(records, list):
            raise ValueError(f"manifest {section!r} must be a list")
        for record in records:
            if not isinstance(record, dict) or "path" not in record:
                raise ValueError(f"manifest {section!r} entry is missing 'path'")
            digest = record.get("sha256")
            if digest is None:
                if record.get("status") == "missing":
                    missing.append({"path": str(record["path"]), "section": section})
                continue
            path = Path(str(record["path"]))
            resolved = path if path.is_absolute() else base / path
            if not resolved.is_file():
                missing.append({"path": str(record["path"]), "section": section})
                continue
            checked += 1
            actual = _sha256_file(resolved)
            if actual != digest:
                mismatches.append(
                    {
                        "path": str(record["path"]),
                        "section": section,
                        "expected_sha256": digest,
                        "actual_sha256": actual,
                    }
                )
    return {
        "ok": not mismatches and not missing,
        "checked": checked,
        "mismatches": mismatches,
        "missing": missing,
    }


class BatchExportError(RuntimeError):
    """A cooling-time batch failed after zero or more outputs were written.

    ``partial`` carries the same fields as a batch manifest with
    ``"status": "partial"``; no success manifest is written for the batch.
    """

    def __init__(self, message: str, *, partial: dict[str, Any]) -> None:
        super().__init__(message)
        self.partial = partial


def export_cooling_batch(
    sources: Sequence[dict[str, Any]],
    requested_times_s: Sequence[float],
    out_dir: str | Path,
    *,
    writer: Callable[[dict[str, Any], Path], None] | None = None,
    overwrite: bool = False,
    skip_zero: bool = False,
) -> dict[str, Any]:
    """Export one source per exact requested cooling offset plus a manifest.

    ``sources`` carry ``cooling_s`` plus the photon-source fields; matching is
    exact (no interpolation, nearest-time substitution or duplicate offsets).
    Outputs are named ``source-t<repr>.xml`` from the round-trippable float
    spelling, so distinct times cannot collide. Every selection, report and
    name is prevalidated before the first write. Files are replaced
    atomically; existing files are refused unless ``overwrite=True``. An
    all-zero source is an error unless ``skip_zero=True``, in which case the
    skipped time is recorded instead of exported. On failure the raised
    :class:`BatchExportError` carries partial outputs and no complete manifest
    is written. The default writer emits OpenMC ``settings.xml`` via
    :func:`write_mesh_source_xml`.
    """
    selected = list(sources)
    if not selected:
        raise ValueError("export_cooling_batch needs at least one source")
    times = [float(t) for t in requested_times_s]
    if not times:
        raise ValueError("requested_times_s must not be empty")
    for t in times:
        if not math.isfinite(t) or t < 0.0:
            raise ValueError(f"cooling offsets must be finite and nonnegative, got {t!r}")
    available: list[float] = []
    for index, source in enumerate(selected):
        cooling = source.get("cooling_s")
        if cooling is None:
            raise ValueError(f"source {index} is missing cooling_s")
        value = float(cooling)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"source {index} cooling_s must be finite and nonnegative")
        available.append(value)
    indices = _core.r2s_select_cooling(available, times)
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    plans: list[tuple[float, dict[str, Any], Path, dict[str, Any], str]] = []
    skipped: list[dict[str, Any]] = []
    for t, index in zip(times, indices, strict=True):
        source = selected[index]
        args = _source_args(source)
        report = _source_report_from_args(args)
        if not report["valid"]:
            codes = ", ".join(
                issue["code"] for issue in report["issues"] if issue["severity"] == "error"
            )
            raise ValueError(f"source at cooling_s={t!r} is invalid; errors: {codes}")
        total = float(report["totals"]["total_photons_per_s"])
        if total <= 0.0:
            if skip_zero:
                skipped.append({"time_s": t, "reason": "all-zero-source"})
                continue
            raise ValueError(
                f"source at cooling_s={t!r} has zero total photons/s; an all-zero "
                "distribution cannot be exported as a sampling source (pass skip_zero=True "
                "to record it as skipped instead)"
            )
        plans.append(
            (
                t,
                source,
                directory / f"source-t{_time_token(t)}.xml",
                report,
                _energy_grid_sha256(args[4]),
            )
        )
    manifest_path = directory / "batch-manifest.json"
    if not overwrite:
        existing = [str(path) for _, _, path, _, _ in plans if path.exists()]
        if manifest_path.exists():
            existing.append(str(manifest_path))
        if existing:
            raise FileExistsError(
                "refusing to overwrite existing batch outputs without overwrite=True: "
                + ", ".join(existing)
            )
    if len({path for _, _, path, _, _ in plans}) != len(plans):
        raise ValueError("deterministic output names collided; distinct cooling offsets required")
    writer_fn = writer if writer is not None else write_mesh_source_xml
    outputs: list[dict[str, Any]] = []
    batch_meta: dict[str, Any] = {
        "schema": _BATCH_SCHEMA,
        "created_by": {"nucleide": _core.version()},
        "requested_times_s": times,
        "normalization": "integrated-per-voxel-and-group",
        "units": "photons/s",
        "voxel_order": _VOXEL_ORDER,
        "assumptions": list(_ASSUMPTIONS),
        "outputs": outputs,
        "skipped": skipped,
    }
    partial: dict[str, Any] = dict(batch_meta, status="partial")
    try:
        for t, source, target, report, grid_hash in plans:
            fd, tmp_name = tempfile.mkstemp(dir=directory, prefix=f".{target.name}.", suffix=".tmp")
            os.close(fd)
            tmp_path = Path(tmp_name)
            try:
                writer_fn(source, tmp_path)
                if not tmp_path.is_file():
                    raise RuntimeError(f"writer did not produce {tmp_path}")
                os.replace(tmp_path, target)
            except BaseException:
                tmp_path.unlink(missing_ok=True)
                raise
            total = float(report["totals"]["total_photons_per_s"])
            outputs.append(
                {
                    "time_s": t,
                    "path": target.name,
                    "sha256": _sha256_file(target),
                    "size_bytes": target.stat().st_size,
                    "total_photons_per_s": total,
                    "energy_grid_sha256": grid_hash,
                    "validation": _report_summary(report),
                }
            )
            partial["outputs"] = outputs
        manifest: dict[str, Any] = dict(batch_meta, status="complete")
        write_manifest(manifest, manifest_path)
    except BaseException as exc:
        raise BatchExportError(
            f"cooling batch failed after {len(outputs)} output(s); no complete manifest written",
            partial=partial,
        ) from exc
    return {
        "manifest": manifest,
        "manifest_path": str(manifest_path),
        "outputs": outputs,
        "skipped": skipped,
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


def _csv_normalization(arrays: dict[str, Any]) -> str:
    """Normalization state recorded on a tally-arrays dict (validated)."""
    normalization = str(arrays.get("normalization", NORMALIZATION_UNKNOWN))
    if normalization not in _NORMALIZATION_STATES:
        raise ValueError(
            f"unknown tally normalization {normalization!r}; want one of "
            f"{sorted(_NORMALIZATION_STATES)}"
        )
    return normalization


def _default_units(normalization: str) -> str:
    """Default units for a declared normalization state."""
    return {
        NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE: UNITS_VOLUME_INTEGRATED_PER_SOURCE,
        NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED: UNITS_VOLUME_INTEGRATED_RATE_SCALED,
        NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE: UNITS_VOLUME_AVERAGED_PER_SOURCE,
        NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED: UNITS_VOLUME_AVERAGED_RATE_SCALED,
        NORMALIZATION_PHYSICAL: UNITS_PHYSICAL,
    }.get(normalization, "unknown")


class _SourceArgs(NamedTuple):
    """Validated photon-source fields passed to the shared Rust kernel."""

    mesh_shape: list[int]
    lower_left_cm: list[float]
    upper_right_cm: list[float]
    volumes_cm3: list[float]
    edges_ev: list[float]
    cooling_s: float
    strengths: list[list[float]]


def _source_args(source: dict[str, Any]) -> _SourceArgs:
    """Extract and type-check the photon-source dict fields."""
    if not isinstance(source, dict):
        raise TypeError("photon source must be the dict returned by nucleide.r2s.photon_source")
    required = (
        "mesh_shape",
        "lower_left_cm",
        "upper_right_cm",
        "voxel_volumes_cm3",
        "edges_ev",
        "cooling_s",
        "strengths_photons_per_s",
    )
    missing = [key for key in required if key not in source]
    if missing:
        raise ValueError(
            "photon source is missing fields: "
            + ", ".join(missing)
            + " (build it with nucleide.r2s.photon_source)"
        )
    try:
        cooling = float(source["cooling_s"])
    except (TypeError, ValueError) as exc:
        raise ValueError("cooling_s must be a finite nonnegative float") from exc
    return _SourceArgs(
        mesh_shape=_int_list(source["mesh_shape"], "mesh_shape"),
        lower_left_cm=_float_list(source["lower_left_cm"], "lower_left_cm"),
        upper_right_cm=_float_list(source["upper_right_cm"], "upper_right_cm"),
        volumes_cm3=_float_list(source["voxel_volumes_cm3"], "voxel_volumes_cm3"),
        edges_ev=_float_list(source["edges_ev"], "edges_ev"),
        cooling_s=cooling,
        strengths=_float_matrix(source["strengths_photons_per_s"], "strengths_photons_per_s"),
    )


def _source_report_from_args(args: _SourceArgs) -> dict[str, Any]:
    """Run the shared Rust report on validated source fields."""
    return _core.r2s_source_report(*args)


def _report_summary(report: dict[str, Any]) -> dict[str, Any]:
    """Compact stable validation summary for a manifest."""
    codes: list[str] = []
    for issue in report.get("issues", []):
        code = str(issue.get("code", ""))
        if code and code not in codes:
            codes.append(code)
    return {
        "valid": bool(report.get("valid", False)),
        "issue_count": int(report.get("issue_count", 0)),
        "truncated": bool(report.get("truncated", False)),
        "codes": codes,
    }


def _validation_summary(validation: dict[str, Any]) -> dict[str, Any]:
    """Normalize a caller-supplied validation summary or report to stable keys."""
    if not isinstance(validation, dict):
        raise TypeError("validation must be a dict with at least a 'valid' flag")
    if "valid" not in validation:
        raise ValueError("validation summary is missing 'valid'")
    if "codes" not in validation and isinstance(validation.get("issues"), list):
        summary = _report_summary(validation)
    else:
        summary = {
            "valid": bool(validation["valid"]),
            "issue_count": int(validation.get("issue_count", 0)),
            "truncated": bool(validation.get("truncated", False)),
            "codes": [str(code) for code in validation.get("codes", [])],
        }
    return summary


def _file_record(path: str | Path, *, require_exists: bool) -> dict[str, Any]:
    """Measure one file (SHA-256 + size) or mark it explicitly missing."""
    file_path = Path(path)
    record: dict[str, Any] = {"path": file_path.as_posix()}
    if not file_path.is_file():
        if require_exists:
            raise FileNotFoundError(f"manifest input not found: {file_path}")
        record.update({"sha256": None, "size_bytes": None, "status": "missing"})
        return record
    record.update({"sha256": _sha256_file(file_path), "size_bytes": file_path.stat().st_size})
    return record


def _sha256_file(path: Path) -> str:
    """Streaming SHA-256 of one file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _energy_grid_sha256(edges_ev: Sequence[float]) -> str:
    """Deterministic identity of an energy grid (little-endian float64 edges)."""
    payload = b"nucleide-energy-grid-ev\x00" + b"".join(
        struct.pack("<d", float(edge)) for edge in edges_ev
    )
    return hashlib.sha256(payload).hexdigest()


def _time_token(time_s: float) -> str:
    """Round-trippable, filename-safe spelling of one cooling offset."""
    token = repr(0.0 if time_s == 0.0 else float(time_s))
    return token.replace("+", "p")


def _write_text_atomic(target: Path, text: str) -> None:
    """Write UTF-8 text through a same-directory temporary file."""
    fd, name = tempfile.mkstemp(dir=target.parent, prefix=f".{target.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
        os.replace(name, target)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def _sum_row(row: Sequence[float]) -> float:
    """Accurate sum of one strength row."""
    return math.fsum(float(value) for value in row)


def _sum_table(table: Sequence[Sequence[float]]) -> float:
    """Accurate sum of a rectangular strength table."""
    return math.fsum(_sum_row(row) for row in table)


def _int_list(value: Any, name: str) -> list[int]:
    """Validate a sequence of integers (bools rejected)."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of ints")
    out: list[int] = []
    for entry in value:
        if isinstance(entry, bool) or not isinstance(entry, int):
            raise TypeError(f"{name} must contain only ints, got {entry!r}")
        out.append(entry)
    return out


def _float_list(value: Any, name: str) -> list[float]:
    """Validate a sequence of finite floats (bools rejected)."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of floats")
    out: list[float] = []
    for entry in value:
        if isinstance(entry, bool):
            raise TypeError(f"{name} must contain only floats, got {entry!r}")
        try:
            number = float(entry)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} must contain only floats, got {entry!r}") from exc
        if not math.isfinite(number):
            raise ValueError(f"{name} must contain only finite values, got {entry!r}")
        out.append(number)
    return out


def _float_matrix(value: Any, name: str) -> list[list[float]]:
    """Validate a rectangular sequence of float sequences."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of float sequences")
    rows: list[list[float]] = []
    width: int | None = None
    for row in value:
        parsed = _float_list(row, name)
        if width is None:
            width = len(parsed)
        elif len(parsed) != width:
            raise ValueError(f"{name} must be rectangular; found rows of {width} and {len(parsed)}")
        rows.append(parsed)
    return rows
