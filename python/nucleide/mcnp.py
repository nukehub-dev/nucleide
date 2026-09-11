"""MCNP file readers and writers (backed by the `nucleide-mcnp-io` crate)."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from typing import Any

from nucleide._internal import (
    DeckProblem,
    EndlLibrary,
    Mctal,
    Meshtal,
    MeshTally,
    PtracFile,
    SurfSrc,
    Wwinp,
    Xsdir,
    XsdirTable,
    endl_endftod,
    mesh_to_geom,
    parse_deck,
    read_deck,
    read_endl,
    read_inp,
    read_mctal,
    read_meshtal,
    read_ptrac,
    read_ssw,
    read_wwinp,
    read_xsdir,
    write_ssw,
)
from nucleide._internal import (
    combine_ssw_files as _combine_ssw_files,
)

__all__ = [
    "Xsdir",
    "XsdirTable",
    "Meshtal",
    "MeshTally",
    "Wwinp",
    "Mctal",
    "SurfSrc",
    "PtracFile",
    "EndlLibrary",
    "DeckProblem",
    "read_xsdir",
    "read_meshtal",
    "read_wwinp",
    "read_mctal",
    "read_ssw",
    "read_ptrac",
    "read_endl",
    "endl_endftod",
    "write_ssw",
    "combine_ssw_files",
    "ssw_combine_main",
    "read_inp",
    "parse_deck",
    "read_deck",
    "mesh_to_geom",
    "ptrac_event_columns",
    "ptrac_event_rows",
    "write_ptrac_hdf5",
    "ptrac_to_hdf5_main",
    "meshtal_mesh_data",
]

#: PTRAC event-table columns in `pyne.mcnp.PtracEvent` order: `event_type`
#: plus the 18 mapped data columns (`node` … `tme`).
_PTRAC_EVENT_COLUMNS: tuple[str, ...] = (
    "event_type",
    "node",
    "nsr",
    "nsf",
    "nxs",
    "ntyn",
    "ipt",
    "ncl",
    "mat",
    "ncp",
    "xxx",
    "yyy",
    "zzz",
    "uuu",
    "vvv",
    "www",
    "erg",
    "wgt",
    "tme",
)


def ptrac_event_columns() -> tuple[str, ...]:
    """PTRAC event-table column names in `PtracEvent` order (19 columns)."""
    return _PTRAC_EVENT_COLUMNS


#: Progress reports while writing PTRAC tables fire every this many events
#: (mirrors upstream `print_progress=1000000`; tests monkeypatch it down).
_PTRAC_PROGRESS_EVERY = 1_000_000


def ptrac_event_rows(path: str) -> list[dict[str, float]]:
    """Stream a PTRAC file into row dicts matching the `PtracEvent` schema.

    Every row carries all 19 columns; variables absent from the file's
    variable list read as 0.0. (Upstream `write_to_hdf5_table` leaves unset
    columns carrying the previous row's values; rows here are independent.)
    """
    ptrac = read_ptrac(path)
    rows: list[dict[str, float]] = []
    for event in ptrac.events():
        rows.append({col: float(event.get(col, 0.0)) for col in _PTRAC_EVENT_COLUMNS})
    return rows


def write_ptrac_hdf5(
    ptrac_path: str,
    hdf5_path: str,
    table_name: str = "ptrac",
    table_title: str = "Ptrac data",
    show_progress: bool = False,
) -> int:
    """Write a PTRAC file's event stream to an HDF5 table (`ptrac_to_hdf5`).

    Rows come from :func:`ptrac_event_rows` (19 `float64` columns in
    `PtracEvent` order); HDF5 bytes are never asserted — the format is not
    byte-stable. Requires `h5py` at runtime. Returns the event count.
    """
    try:
        import h5py as h5py_mod
    except ImportError as exc:
        raise ImportError(
            "write_ptrac_hdf5 requires h5py (pip install h5py); "
            "use ptrac_event_rows() for the dependency-free row stream"
        ) from exc

    import numpy as np

    rows = ptrac_event_rows(ptrac_path)
    dtype = [(col, "<f8") for col in _PTRAC_EVENT_COLUMNS]
    data = np.empty(len(rows), dtype=dtype)
    for i, row in enumerate(rows):
        for col in _PTRAC_EVENT_COLUMNS:
            data[col][i] = row[col]
        if show_progress and (i + 1) % _PTRAC_PROGRESS_EVERY == 0:
            print(f"processing event {i + 1}")
    with h5py_mod.File(hdf5_path, "a") as h5file:
        if table_name in h5file:
            del h5file[table_name]
        dataset = h5file.create_dataset(table_name, data=data)
        dataset.attrs["TITLE"] = table_title
    return len(rows)


def ptrac_to_hdf5_main(argv: Sequence[str] | None = None) -> int:
    """CLI for :func:`write_ptrac_hdf5` (mirrors `pyne.ptrac_to_hdf5`)."""
    parser = argparse.ArgumentParser(
        description="write the contents of a MCNP PTRAC file to a HDF5 table"
    )
    parser.add_argument("ptrac_file", help="MCNP PTRAC file to read from")
    parser.add_argument(
        "hdf5_file", help="HDF5 file to write to (will be created if it does not exist)"
    )
    parser.add_argument(
        "-n",
        "--table-name",
        default="ptrac",
        help='name of the HDF5 table (default is "ptrac")',
    )
    parser.add_argument(
        "-t",
        "--table-title",
        default="Ptrac data",
        help='title of the HDF5 table (default is "Ptrac data")',
    )
    parser.add_argument(
        "-s", "--show-progress", action="store_true", help="show progress indicator"
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    return write_ptrac_hdf5(
        args.ptrac_file,
        args.hdf5_file,
        table_name=args.table_name,
        table_title=args.table_title,
        show_progress=args.show_progress,
    )


def combine_ssw_files(output: str, inputs: Sequence[str]) -> None:
    """Combine SSW surface-source files into one (`ssw_combine` port).

    The output header comes from the first file with the signed `orignp1`
    sum and the plain `nrss` sum; later files' track `nps` shift
    sign-preservingly. Raises `ValueError` on incompatible headers
    (upstream prints and returns `False`).
    """
    if not inputs:
        raise ValueError("combine_ssw_files needs at least one input file")
    _combine_ssw_files(output, list(inputs))


def ssw_combine_main(argv: Sequence[str] | None = None) -> None:
    """CLI for :func:`combine_ssw_files` (mirrors `ssw_combine.py`)."""
    names = list(argv) if argv is not None else sys.argv[1:]
    combine_ssw_files("newssr", names)


def meshtal_mesh_data(path: str, as_numpy: bool = False) -> dict[str, Any]:
    """Extract meshtal tallies as plain dicts/arrays (no MOAB tagging).

    Returns `{"version", "ld", "title", "histories", "tallies"}` where each
    tally carries its bounds, per-cell `result`/`rel_error` (`[ve][group]`,
    x slowest → z fastest), energy-integrated `total_result`/
    `total_rel_error`, `dims`, `num_ves`, and `num_e_groups`. With
    `as_numpy=True` the numeric tables are NumPy arrays instead of lists.
    MOAB tag mapping (`tag_flux_error_from_tally_results`) stays caller-side:
    apply these arrays to caller-owned tags.
    """
    meshtal = read_meshtal(path)
    tallies: dict[int, dict[str, Any]] = {}
    for number, tally in meshtal.tallies.items():
        entry: dict[str, Any] = {
            "tally_number": tally.tally_number,
            "particle": tally.particle,
            "dose_response": tally.dose_response,
            "x_bounds": list(tally.x_bounds),
            "y_bounds": list(tally.y_bounds),
            "z_bounds": list(tally.z_bounds),
            "e_bounds": list(tally.e_bounds),
            "dims": list(tally.dims()),
            "num_ves": tally.num_ves(),
            "num_e_groups": tally.num_e_groups(),
            "result": tally.result,
            "rel_error": tally.rel_error,
            "total_result": tally.total_result,
            "total_rel_error": tally.total_rel_error,
        }
        tallies[int(number)] = entry
    if as_numpy:
        import numpy as np

        for entry in tallies.values():
            entry["result"] = np.asarray(entry["result"], dtype=float)
            entry["rel_error"] = np.asarray(entry["rel_error"], dtype=float)
            entry["total_result"] = np.asarray(entry["total_result"], dtype=float)
            entry["total_rel_error"] = np.asarray(entry["total_rel_error"], dtype=float)
    return {
        "version": meshtal.version,
        "ld": meshtal.ld,
        "title": meshtal.title,
        "histories": meshtal.histories,
        "tallies": tallies,
    }
