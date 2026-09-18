"""Tests for the OpenMC statepoint tally bridge (synthetic statepoints).

OpenMC is never imported here: the fakes below implement the documented
synthetic-statepoint protocol under the real 0.16.0 class names
(``EnergyFilter`` / ``MeshFilter`` / ``RegularMesh``), and ``FakeTally``
asserts the pinned ``get_reshaped_data(value=, expand_dims=True)`` call
sequence. The validation harness replays the same assertions against real
OpenMC objects (loud SKIP without OpenMC).
"""

import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import nucleide.damage as dmg
import nucleide.openmc as omc


class FakeData:
    """Nested-list stand-in for the reshaped tally array (``.shape``/``.tolist``)."""

    def __init__(self, nested: list[Any]) -> None:
        self._nested = nested

    @property
    def shape(self) -> tuple[int, ...]:
        shape: list[int] = []
        node: Any = self._nested
        while isinstance(node, list):
            shape.append(len(node))
            node = node[0] if node else None
        return tuple(shape)

    def tolist(self) -> list[Any]:
        return self._nested


class FakeTally:
    """Fake openmc.Tally recording the pinned numeric-fetch sequence."""

    def __init__(
        self,
        scores: list[str],
        nuclides: list[str],
        filters: list[Any],
        mean: list[Any],
        std: list[Any],
    ) -> None:
        self.scores = scores
        self.nuclides = nuclides
        self.filters = filters
        self._mean = mean
        self._std = std
        self.calls: list[str] = []

    def get_reshaped_data(self, value: str = "mean", expand_dims: bool = False) -> FakeData:
        assert expand_dims, "bridge must request expanded mesh dims"
        self.calls.append(value)
        if value == "mean":
            return FakeData(self._mean)
        if value == "std_dev":
            return FakeData(self._std)
        raise AssertionError(f"unexpected value={value!r}")


class EnergyFilter:
    """Fake openmc.EnergyFilter: float eV bin edges."""

    def __init__(self, bins: list[float]) -> None:
        self.bins = list(bins)
        self.num_bins = len(bins) - 1


class RegularMesh:
    """Fake openmc.RegularMesh (centimetres)."""

    def __init__(
        self,
        dimension: tuple[int, int, int],
        lower: tuple[float, float, float],
        upper: tuple[float, float, float],
    ) -> None:
        self.dimension = dimension
        self.lower_left = lower
        self.upper_right = upper


class MeshFilter:
    """Fake openmc.MeshFilter."""

    def __init__(self, mesh: RegularMesh) -> None:
        self.mesh = mesh


class CellFilter:
    """Fake openmc.CellFilter: int cell ids."""

    def __init__(self, bins: list[int]) -> None:
        self.bins = list(bins)
        self.num_bins = len(bins)


class FakeStatePoint:
    """Fake openmc.StatePoint: a dict of tallies."""

    def __init__(self, tallies: dict[int, FakeTally]) -> None:
        self.tallies = tallies


class _BlockOpenMC:
    """Meta-path finder making ``import openmc`` fail (caller-side stance)."""

    def find_spec(
        self,
        fullname: str,
        path: Any = None,
        target: Any = None,
    ) -> Any:
        if fullname == "openmc" or fullname.startswith("openmc."):
            raise ImportError("blocked for test")
        return None


def _cell_tally() -> FakeTally:
    return FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[EnergyFilter([0.0, 1.0e5, 2.0e7])],
        mean=[[[1.5]], [[2.5]]],
        std=[[[0.15]], [[0.25]]],
    )


def _mesh_nested(values: list[list[float]], nx: int, ny: int, nz: int) -> list[Any]:
    """Nest ``[voxel][group]`` values C-order into ``(x, y, z, group, nuclide, score)``."""
    out: list[Any] = []
    for i in range(nx):
        plane: list[Any] = []
        for j in range(ny):
            column: list[Any] = []
            for _k in range(nz):
                column.append([[[v]] for v in values[i * ny + j]])
            plane.append(column)
        out.append(plane)
    return out


def _cell_statepoint() -> FakeStatePoint:
    return FakeStatePoint({7: _cell_tally()})


def test_cell_energy_tally_exact() -> None:
    out = omc.tally_arrays(_cell_statepoint(), 7)
    assert out["tally_id"] == 7
    assert out["score"] == "flux"
    assert out["kind"] == "cell"
    assert out["n_voxels"] == 1
    assert out["n_groups"] == 2
    assert out["energy_bounds_mev"] == [0.0, 0.1, 20.0]
    assert out["flux"] == [[1.5, 2.5]]
    assert out["rel_err"] == [[pytest.approx(0.1), pytest.approx(0.1)]]
    assert out["mesh_shape"] is None
    assert out["mesh_bounds"] is None
    assert out["openmc_version"] == omc.openmc_version()


def test_numeric_fetch_sequence_pinned() -> None:
    tally = _cell_tally()
    omc.tally_arrays(FakeStatePoint({7: tally}), 7)
    assert tally.calls == ["mean", "std_dev"]


def test_mesh_tally_c_order_and_units() -> None:
    values = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]]
    ordered = _mesh_nested(values, 2, 2, 1)
    zeros = _mesh_nested([[0.0, 0.0]] * 4, 2, 2, 1)
    mesh = RegularMesh((2, 2, 1), (-4.0, -4.0, -4.0), (4.0, 4.0, 4.0))
    tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[MeshFilter(mesh), EnergyFilter([0.0, 1.0e6, 2.0e6])],
        mean=ordered,
        std=zeros,
    )
    out = omc.tally_arrays(FakeStatePoint({3: tally}), 3)
    assert out["kind"] == "mesh"
    assert out["n_voxels"] == 4
    assert out["n_groups"] == 2
    assert out["flux"] == values
    assert out["rel_err"] == [[0.0, 0.0]] * 4
    assert out["energy_bounds_mev"] == [0.0, 1.0, 2.0]
    assert out["mesh_shape"] == [2, 2, 1]
    assert out["mesh_bounds"] == {"x": [-4.0, 0.0, 4.0], "y": [-4.0, 0.0, 4.0], "z": [-4.0, 4.0]}


def test_unbinned_tally_is_single_total_group() -> None:
    tally = FakeTally(scores=["flux"], nuclides=["total"], filters=[], mean=[[5.0]], std=[[0.5]])
    out = omc.tally_arrays(FakeStatePoint({1: tally}), 1)
    assert out["n_groups"] == 1
    assert out["energy_bounds_mev"] == []
    assert out["flux"] == [[5.0]]
    assert out["rel_err"] == [[pytest.approx(0.1)]]


def test_singleton_cell_filter_is_transparent() -> None:
    tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[CellFilter([3]), EnergyFilter([0.0, 2.0e7])],
        mean=[[[[7.0]]]],
        std=[[[[0.0]]]],
    )
    out = omc.tally_arrays(FakeStatePoint({1: tally}), 1)
    assert out["flux"] == [[7.0]]
    assert out["rel_err"] == [[0.0]]


def test_zero_mean_reports_zero_rel_err() -> None:
    tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[EnergyFilter([0.0, 2.0e7])],
        mean=[[[0.0]]],
        std=[[[0.3]]],
    )
    out = omc.tally_arrays(FakeStatePoint({1: tally}), 1)
    assert out["flux"] == [[0.0]]
    assert out["rel_err"] == [[0.0]]


def test_cell_flux_feeds_damage_fold() -> None:
    arrays = omc.tally_arrays(_cell_statepoint(), 7)
    flux = omc.cell_flux(arrays)
    response = [100.0, 200.0]
    hand = 2.0e-24 * (flux[0] * response[0] + flux[1] * response[1])
    assert dmg.nrt_dpa(flux, response, arrays["energy_bounds_mev"], 2.0) == pytest.approx(hand)


def test_voxel_flux_selection_and_errors() -> None:
    values = [[1.0, 10.0], [2.0, 20.0]]
    mesh = RegularMesh((2, 1, 1), (0.0, 0.0, 0.0), (2.0, 1.0, 1.0))
    mean = _mesh_nested(values, 2, 1, 1)
    std = _mesh_nested([[0.0, 0.0]] * 2, 2, 1, 1)
    tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[MeshFilter(mesh), EnergyFilter([0.0, 1.0e6, 2.0e7])],
        mean=mean,
        std=std,
    )
    arrays = omc.tally_arrays(FakeStatePoint({3: tally}), 3)
    assert omc.voxel_flux(arrays, 1) == [2.0, 20.0]
    with pytest.raises(ValueError, match="out of range"):
        omc.voxel_flux(arrays, 2)
    with pytest.raises(ValueError, match="out of range"):
        omc.voxel_flux(arrays, -1)
    with pytest.raises(TypeError, match="must be an int"):
        omc.voxel_flux(arrays, True)
    with pytest.raises(ValueError, match="single-voxel"):
        omc.cell_flux(arrays)


def test_csv_roundtrip_mesh_cell_and_unbinned(tmp_path: Path) -> None:
    mesh = RegularMesh((2, 1, 1), (-4.0, -4.0, -4.0), (4.0, 4.0, 4.0))
    mesh_tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[MeshFilter(mesh), EnergyFilter([0.0, 1.0e6, 2.0e6])],
        mean=_mesh_nested([[1.0, 2.0], [3.0, 4.0]], 2, 1, 1),
        std=_mesh_nested([[0.1, 0.2], [0.3, 0.4]], 2, 1, 1),
    )
    cases = [
        omc.tally_arrays(FakeStatePoint({3: mesh_tally}), 3),
        omc.tally_arrays(_cell_statepoint(), 7),
        omc.tally_arrays(
            FakeStatePoint(
                {
                    1: FakeTally(
                        scores=["flux"], nuclides=["total"], filters=[], mean=[[5.0]], std=[[0.5]]
                    )
                }
            ),
            1,
        ),
    ]
    for index, arrays in enumerate(cases):
        path = tmp_path / f"tally-{index}.csv"
        omc.tally_to_csv(arrays, path)
        assert omc.tally_from_csv(path) == arrays


def test_csv_malformed_is_loud(tmp_path: Path) -> None:
    arrays = omc.tally_arrays(_cell_statepoint(), 7)
    good = tmp_path / "good.csv"
    omc.tally_to_csv(arrays, good)
    text = good.read_text(encoding="utf-8")
    bad_format = tmp_path / "bad-format.csv"
    bad_format.write_text(text.replace("nucleide-openmc-tally-v1", "other"), encoding="utf-8")
    with pytest.raises(ValueError, match="not a nucleide-openmc-tally-v1"):
        omc.tally_from_csv(bad_format)
    bad_header = tmp_path / "bad-header.csv"
    bad_header.write_text(text.replace("voxel,group,lower_mev", "voxel,group,X"), encoding="utf-8")
    with pytest.raises(ValueError, match="header row"):
        omc.tally_from_csv(bad_header)
    short = tmp_path / "short.csv"
    short.write_text("\n".join(text.splitlines()[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="holds 1 rows"):
        omc.tally_from_csv(short)
    shuffled = tmp_path / "shuffled.csv"
    lines = text.splitlines()
    lines[-2:] = [lines[-1], lines[-2]]
    shuffled.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="breaks \\(voxel, group\\) order"):
        omc.tally_from_csv(shuffled)
    with pytest.raises(FileNotFoundError, match="not found"):
        omc.tally_from_csv(tmp_path / "missing.csv")


def test_missing_tally_score_and_shape_are_loud() -> None:
    with pytest.raises(KeyError, match="available tally ids"):
        omc.tally_arrays(_cell_statepoint(), 99)
    with pytest.raises(ValueError, match="available scores"):
        omc.tally_arrays(_cell_statepoint(), 7, score="fission")
    multi = FakeTally(
        scores=["flux"],
        nuclides=["U235", "U238"],
        filters=[EnergyFilter([0.0, 2.0e7])],
        mean=[[[1.0], [2.0]]],
        std=[[[0.0], [0.0]]],
    )
    with pytest.raises(ValueError, match="nuclides"):
        omc.tally_arrays(FakeStatePoint({1: multi}), 1)
    wide_cell = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[CellFilter([1, 2]), EnergyFilter([0.0, 2.0e7])],
        mean=[[[1.0]], [[2.0]]],
        std=[[[0.0]], [[0.0]]],
    )
    with pytest.raises(ValueError, match="multi-bin CellFilter"):
        omc.tally_arrays(FakeStatePoint({1: wide_cell}), 1)
    two_energy = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[EnergyFilter([0.0, 1.0]), EnergyFilter([0.0, 1.0])],
        mean=[[[1.0]]],
        std=[[[0.0]]],
    )
    with pytest.raises(ValueError, match="two energy filters"):
        omc.tally_arrays(FakeStatePoint({1: two_energy}), 1)
    with pytest.raises(TypeError, match="exposing .tallies"):
        omc.tally_arrays(object(), 1)


def test_unsupported_mesh_is_loud() -> None:
    class CylindricalMesh:
        def __init__(self) -> None:
            self.dimension = (2, 3, 1)

    class CylindricalMeshFilter:
        def __init__(self) -> None:
            self.mesh = CylindricalMesh()

    tally = FakeTally(
        scores=["flux"],
        nuclides=["total"],
        filters=[CylindricalMeshFilter(), EnergyFilter([0.0, 2.0e7])],
        mean=[[[1.0]]],
        std=[[[0.0]]],
    )
    with pytest.raises(ValueError, match="only RegularMesh"):
        omc.tally_arrays(FakeStatePoint({1: tally}), 1)


def test_empty_tally_results_propagate() -> None:
    class EmptyTally(FakeTally):
        def get_reshaped_data(self, value: str = "mean", expand_dims: bool = False) -> FakeData:
            raise ValueError("called before the Tally is populated with data")

    tally = EmptyTally(
        scores=["flux"], nuclides=["total"], filters=[EnergyFilter([0.0, 2.0e7])], mean=[], std=[]
    )
    with pytest.raises(ValueError, match="populated with data"):
        omc.tally_arrays(FakeStatePoint({1: tally}), 1)


def test_energy_bins_must_be_strictly_increasing() -> None:
    for bins in ([1.0e6, 1.0e6], [2.0e6, 1.0e6], [5.0]):
        tally = FakeTally(
            scores=["flux"],
            nuclides=["total"],
            filters=[EnergyFilter(bins)],
            mean=[[[1.0]]],
            std=[[[0.0]]],
        )
        with pytest.raises(ValueError, match="strictly increasing"):
            omc.tally_arrays(FakeStatePoint({1: tally}), 1)


def test_read_statepoint_missing_file_needs_no_openmc(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="statepoint not found"):
        omc.read_statepoint(tmp_path / "statepoint.10.h5", 1)


def test_read_statepoint_missing_openmc_is_loud(tmp_path: Path) -> None:
    path = tmp_path / "statepoint.10.h5"
    path.write_bytes(b"fake")
    blocker = _BlockOpenMC()
    sys.meta_path.insert(0, blocker)
    try:
        with pytest.raises(ImportError, match="caller-side OpenMC"):
            omc.read_statepoint(path, 1)
    finally:
        sys.meta_path.remove(blocker)
    assert omc.openmc_version() is None or isinstance(omc.openmc_version(), str)


def test_read_statepoint_pinned_sequence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "statepoint.10.h5"
    path.write_bytes(b"fake")
    seen: dict[str, Any] = {}
    closed: list[bool] = []

    class FakeStatePointFile:
        def __init__(self, filepath: str, autolink: bool = True) -> None:
            seen["path"] = filepath
            seen["autolink"] = autolink
            self.tallies = {7: _cell_tally()}

        def close(self) -> None:
            closed.append(True)

    fake = ModuleType("openmc")
    fake.StatePoint = FakeStatePointFile  # type: ignore[attr-defined]
    fake.__version__ = "0.16.0-fake"  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "openmc", fake)
    out = omc.read_statepoint(path, 7)
    assert seen == {"path": str(path), "autolink": False}
    assert closed == [True]
    assert out["flux"] == [[1.5, 2.5]]
    assert out["openmc_version"] == "0.16.0-fake"
