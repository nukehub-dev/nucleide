"""Activation-to-transport handoff tests (synthetic; no OpenMC import).

Covers physical flux normalization, the v2 CSV metadata interchange, the
shared photon-source contract, the reproducibility manifest and exact
cooling-time batch export. Real pinned-OpenMC index/XML behaviour lives in
``validation/openmc_handoff_probe.py``; the adapter test here uses a fake
``openmc`` module built to the documented 0.16.0 object protocol so the
C-order-to-Fortran-order mapping is checked without the optional dependency.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import nucleide.openmc as omc
import nucleide.r2s as r2s

# ---------------------------------------------------------------------------
# fixtures / helpers


def _tally(flux: list[list[float]], *, normalization: str | None = None) -> dict[str, Any]:
    n_voxels = len(flux)
    n_groups = len(flux[0])
    arrays: dict[str, Any] = {
        "tally_id": 1,
        "score": "flux",
        "kind": "mesh" if n_voxels > 1 else "cell",
        "n_voxels": n_voxels,
        "n_groups": n_groups,
        "energy_bounds_mev": [0.0, 1.0, 2.0][: n_groups + 1],
        "flux": [list(row) for row in flux],
        "rel_err": [[0.1] * n_groups for _ in range(n_voxels)],
        "mesh_shape": [n_voxels, 1, 1] if n_voxels > 1 else None,
        "mesh_bounds": (
            {"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0], "z": [0.0, 1.0]} if n_voxels > 1 else None
        ),
        "openmc_version": "0.16.0",
        "normalization": normalization if normalization is not None else omc.NORMALIZATION_UNKNOWN,
        "units": "unknown",
    }
    return arrays


def _source(
    strengths: list[list[float]] | None = None,
    *,
    cooling_s: float = 0.0,
    mesh_shape: list[int] | None = None,
    volumes: list[float] | None = None,
    edges: list[float] | None = None,
) -> dict[str, Any]:
    strengths = strengths if strengths is not None else [[100.0, 0.0], [0.0, 300.0]]
    mesh_shape = mesh_shape if mesh_shape is not None else [2, 1, 1]
    volumes = volumes if volumes is not None else [1.0, 1.0]
    edges = edges if edges is not None else [0.0, 2.0, 4.0]
    nx, ny, nz = mesh_shape
    return r2s.photon_source(
        mesh_shape,
        [0.0, 0.0, 0.0],
        [float(nx), float(ny), float(nz)],
        volumes,
        edges,
        cooling_s,
        strengths,
    )


# ---------------------------------------------------------------------------
# W1: flux normalization


def test_normalize_flux_hand_case_and_metadata() -> None:
    arrays = _tally([[4.0, 9.0]])
    out = omc.normalize_flux(
        arrays,
        10.0,
        [2.0],
        input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
    )
    assert out["flux"] == [[20.0, 45.0]]
    assert out["rel_err"] == arrays["rel_err"]
    assert out["normalization"] == omc.NORMALIZATION_PHYSICAL
    assert out["units"] == omc.UNITS_PHYSICAL
    assert out["source_rate_s"] == 10.0
    assert out["input_normalization"] == omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE
    # Additive: the raw arrays are untouched.
    assert arrays["flux"] == [[4.0, 9.0]]
    assert arrays["normalization"] == omc.NORMALIZATION_UNKNOWN


def test_normalize_flux_two_unequal_volumes_and_zero_rate() -> None:
    arrays = _tally([[4.0, 9.0], [8.0, 3.0]])
    out = omc.normalize_flux(
        arrays,
        10.0,
        [2.0, 4.0],
        input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
    )
    assert out["flux"] == [[20.0, 45.0], [20.0, 7.5]]
    zero = omc.normalize_flux(
        arrays,
        0.0,
        [2.0, 4.0],
        input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
    )
    assert zero["flux"] == [[0.0, 0.0], [0.0, 0.0]]


def test_normalize_flux_volume_averaged_branch_does_not_divide_again() -> None:
    arrays = _tally([[4.0, 9.0]])
    out = omc.normalize_flux(
        arrays,
        10.0,
        [2.0],
        input_normalization=omc.NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE,
    )
    assert out["flux"] == [[40.0, 90.0]]


def test_normalize_flux_rate_scaled_does_not_multiply_the_rate_again() -> None:
    arrays = _tally([[4.0, 9.0]])
    integrated = omc.normalize_flux(
        arrays,
        10.0,
        [2.0],
        input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_RATE_SCALED,
    )
    assert integrated["flux"] == [[2.0, 4.5]]
    averaged = omc.normalize_flux(
        arrays,
        10.0,
        [2.0],
        input_normalization=omc.NORMALIZATION_VOLUME_AVERAGED_RATE_SCALED,
    )
    assert averaged["flux"] == [[4.0, 9.0]]
    assert averaged["normalization"] == omc.NORMALIZATION_PHYSICAL


def test_normalize_flux_requires_declaration_and_rejects_double_scaling() -> None:
    arrays = _tally([[4.0, 9.0]])
    with pytest.raises(ValueError, match="normalization is unknown"):
        omc.normalize_flux(arrays, 10.0, [2.0])
    with pytest.raises(ValueError, match="already physical"):
        omc.normalize_flux(
            arrays,
            10.0,
            [2.0],
            input_normalization=omc.NORMALIZATION_PHYSICAL,
        )
    recorded = _tally([[4.0, 9.0]], normalization=omc.NORMALIZATION_PHYSICAL)
    with pytest.raises(ValueError, match="already physical"):
        omc.normalize_flux(recorded, 10.0, [2.0])
    with pytest.raises(ValueError, match="disagrees"):
        omc.normalize_flux(
            _tally([[4.0, 9.0]], normalization=omc.NORMALIZATION_VOLUME_AVERAGED_PER_SOURCE),
            10.0,
            [2.0],
            input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
        )
    with pytest.raises(ValueError, match="unsupported input_normalization"):
        omc.normalize_flux(arrays, 10.0, [2.0], input_normalization="mystery")


def test_normalize_flux_rejects_bad_rate_volume_and_values() -> None:
    arrays = _tally([[4.0, 9.0]])
    declared = omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE
    with pytest.raises(ValueError, match="source_rate_s"):
        omc.normalize_flux(arrays, -1.0, [2.0], input_normalization=declared)
    with pytest.raises(ValueError, match="source_rate_s"):
        omc.normalize_flux(arrays, float("nan"), [2.0], input_normalization=declared)
    with pytest.raises(ValueError, match="volumes_cm3"):
        omc.normalize_flux(arrays, 1.0, [2.0, 3.0], input_normalization=declared)
    with pytest.raises(ValueError):
        omc.normalize_flux(_tally([[-1.0, 9.0]]), 1.0, [2.0], input_normalization=declared)
    with pytest.raises(ValueError):
        omc.normalize_flux(_tally([[4.0, 9.0]]), 1.0, [0.0], input_normalization=declared)


def test_csv_v1_reads_as_unknown_and_v2_roundtrips_physical() -> None:
    arrays = _tally([[4.0, 9.0]])
    v1 = "\n".join(
        [
            "# format = nucleide-openmc-tally-v1",
            "# tally_id = 1",
            "# score = flux",
            "# kind = cell",
            "# n_voxels = 1",
            "# n_groups = 2",
            "# mesh_shape = -",
            "# energy_bounds = 0.0,1e-06,2e-06",
            "# openmc_version = 0.16.0",
            "# mesh_x_bounds = -",
            "# mesh_y_bounds = -",
            "# mesh_z_bounds = -",
            "voxel,group,lower_mev,upper_mev,flux,rel_err",
            "0,0,0.0,1e-06,4.0,0.1",
            "0,1,1e-06,2e-06,9.0,0.1",
            "",
        ]
    )
    # The reader accepts v1 metadata but never relabels it as physical.
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as tmp:
        source = Path(tmp) / "legacy.csv"
        source.write_text(v1, encoding="utf-8")
        legacy = omc.tally_from_csv(source)
        assert legacy["normalization"] == omc.NORMALIZATION_UNKNOWN
        assert legacy["units"] == "unknown"
        with pytest.raises(ValueError, match="normalization is unknown"):
            omc.normalize_flux(legacy, 10.0, [2.0])
        converted = omc.normalize_flux(
            legacy,
            10.0,
            [2.0],
            input_normalization=omc.NORMALIZATION_VOLUME_INTEGRATED_PER_SOURCE,
        )
        target = Path(tmp) / "converted.csv"
        omc.tally_to_csv(converted, target)
        back = omc.tally_from_csv(target)
        assert back["flux"] == [[20.0, 45.0]]
        assert back["normalization"] == omc.NORMALIZATION_PHYSICAL
        assert back["units"] == omc.UNITS_PHYSICAL
        assert back["source_rate_s"] == 10.0
        assert back["rel_err"] == arrays["rel_err"]


def test_csv_v2_preserves_raw_metadata(tmp_path: Path) -> None:
    arrays = _tally([[1.0, 2.0]])
    path = tmp_path / "raw.csv"
    omc.tally_to_csv(arrays, path)
    assert "nucleide-openmc-tally-v2" in path.read_text(encoding="utf-8")
    back = omc.tally_from_csv(path)
    assert back == arrays


# ---------------------------------------------------------------------------
# W2: shared photon source contract


def test_photon_source_totals_density_and_joint_probabilities() -> None:
    source = _source()
    assert source["n_voxels"] == 2
    assert source["n_groups"] == 2
    assert source["total_photons_per_s"] == pytest.approx(400.0)
    assert source["voxel_totals_photons_per_s"] == pytest.approx([100.0, 300.0])
    assert source["group_totals_photons_per_s"] == pytest.approx([100.0, 300.0])
    assert r2s.source_density(source) == [[100.0, 0.0], [0.0, 300.0]]
    assert r2s.source_probabilities(source) == [[0.25, 0.0], [0.0, 0.75]]
    report = r2s.source_report(source)
    assert report["valid"] is True
    codes = [issue["code"] for issue in report["issues"]]
    assert "source-assumption-uniform-within-bin" in codes
    assert "source-assumption-uniform-within-voxel" in codes


def test_photon_source_unequal_volumes_and_zero_voxel() -> None:
    source = _source(
        [[100.0, 0.0], [0.0, 300.0]],
        volumes=[1.0, 3.0],
    )
    assert r2s.source_density(source) == [[100.0, 0.0], [0.0, 100.0]]
    # A zero-source voxel is valid inside a nonzero source.
    zero_voxel = _source([[100.0, 0.0], [0.0, 0.0]])
    assert zero_voxel["voxel_totals_photons_per_s"] == pytest.approx([100.0, 0.0])
    assert r2s.source_probabilities(zero_voxel) == [[1.0, 0.0], [0.0, 0.0]]


def test_all_zero_source_is_data_but_not_a_distribution() -> None:
    source = _source([[0.0, 0.0], [0.0, 0.0]])
    assert source["total_photons_per_s"] == 0.0
    assert r2s.source_report(source)["valid"] is True
    with pytest.raises(ValueError):
        r2s.source_probabilities(source)


def test_sparse_group_export_fills_omitted_bins_with_zero() -> None:
    full_edges = [0.0, 1.0, 2.0, 3.0, 4.0]
    source = r2s.photon_source_exported(
        [1, 1, 1],
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [1.0],
        full_edges,
        0.0,
        [1, 3],
        [[5.0, 7.0]],
    )
    assert source["strengths_photons_per_s"] == [[0.0, 5.0, 0.0, 7.0]]
    assert source["total_photons_per_s"] == pytest.approx(12.0)
    with pytest.raises(ValueError):
        r2s.photon_source_exported(
            [1, 1, 1],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0],
            full_edges,
            0.0,
            [1, 1],
            [[5.0, 7.0]],
        )


def test_source_report_flags_negative_strength() -> None:
    source = _source()
    with pytest.raises(ValueError):
        r2s.photon_source(
            source["mesh_shape"],
            source["lower_left_cm"],
            source["upper_right_cm"],
            source["voxel_volumes_cm3"],
            source["edges_ev"],
            source["cooling_s"],
            [[100.0, -1.0], [0.0, 300.0]],
        )
    draft = dict(source)
    draft["strengths_photons_per_s"] = [[100.0, -1.0], [0.0, 300.0]]
    report = r2s.source_report(draft)
    assert report["valid"] is False
    assert any(issue["code"] == "source-strength-negative" for issue in report["issues"])


def test_select_cooling_is_exact_and_order_preserving() -> None:
    assert r2s.select_cooling([0.0, 3600.0, 7200.0], [7200.0, 0.0]) == [2, 0]
    with pytest.raises(ValueError):
        r2s.select_cooling([0.0, 3600.0], [3600.0, 3600.0])
    with pytest.raises(ValueError):
        r2s.select_cooling([0.0, 3600.0], [3600.0000001])
    # Close but distinct offsets stay distinct.
    assert r2s.select_cooling([3600.0, 3600.0000000000005], [3600.0000000000005]) == [1]


# ---------------------------------------------------------------------------
# W3: OpenMC adapter (fake 0.16.0 protocol)


class _FakeTabular:
    def __init__(self, x: list[float], p: list[float], interpolation: str) -> None:
        assert interpolation == "histogram"
        self.x = list(x)
        self.p = list(p)
        self.interpolation = interpolation


class _FakeUniform:
    def __init__(self, a: float, b: float) -> None:
        self.a = a
        self.b = b


class _FakeIndependentSource:
    def __init__(self, *, energy: Any, strength: float, particle: str) -> None:
        self.energy = energy
        self.strength = strength
        self.particle = particle


class _FakeMesh:
    def __init__(self) -> None:
        self.dimension: tuple[int, int, int] | None = None
        self.lower_left: tuple[float, float, float] | None = None
        self.upper_right: tuple[float, float, float] | None = None

    @property
    def n_elements(self) -> int:
        assert self.dimension is not None
        return self.dimension[0] * self.dimension[1] * self.dimension[2]


class _FakeMeshSource:
    def __init__(self, *, mesh: _FakeMesh, sources: Any) -> None:
        self.mesh = mesh
        # OpenMC 0.16.0 flattens multidimensional arrays in Fortran order.
        flat = []
        assert mesh.dimension is not None
        nx, ny, nz = mesh.dimension
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    flat.append(sources[i][j][k])
        self.sources = flat
        self.strength = sum(src.strength for src in flat)


def _install_fake_openmc(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = ModuleType("openmc")
    fake.RegularMesh = _FakeMesh  # type: ignore[attr-defined]
    fake.IndependentSource = _FakeIndependentSource  # type: ignore[attr-defined]
    fake.MeshSource = _FakeMeshSource  # type: ignore[attr-defined]
    fake.__version__ = "0.16.0"  # type: ignore[attr-defined]
    stats = ModuleType("openmc.stats")
    stats.Tabular = _FakeTabular  # type: ignore[attr-defined]
    stats.Uniform = _FakeUniform  # type: ignore[attr-defined]
    fake.stats = stats  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "openmc", fake)
    monkeypatch.setitem(sys.modules, "openmc.stats", stats)


def test_mesh_source_maps_2x3x2_cells_in_c_order(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_openmc(monkeypatch)
    strengths = [[float(100 + i), 1.0] for i in range(12)]
    source = _source(
        strengths,
        mesh_shape=[2, 3, 2],
        volumes=[1.0] * 12,
        edges=[0.0, 1.0, 2.0],
    )
    converted = omc.mesh_source(source)
    # Fake MeshSource stores native OpenMC order (x fastest, then y, then z).
    for i in range(2):
        for j in range(3):
            for k in range(2):
                c_order = (i * 3 + j) * 2 + k
                native = i + 2 * (j + 3 * k)
                assert converted.sources[native].strength == pytest.approx(sum(strengths[c_order]))
    assert converted.mesh.dimension == (2, 3, 2)


def test_mesh_source_energy_pdf_is_per_ev(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_openmc(monkeypatch)
    source = _source(
        [[100.0, 0.0], [0.0, 300.0]],
        edges=[0.0, 2.0, 10.0],
        volumes=[1.0, 1.0],
    )
    converted = omc.mesh_source(source)
    first = converted.sources[0].energy
    assert isinstance(first, _FakeTabular)
    # Group strength / (voxel total * width): 100/(100*2), 0.
    assert first.p == pytest.approx([0.5, 0.0])
    second = converted.sources[1].energy
    assert second.p == pytest.approx([0.0, 300.0 / (300.0 * 8.0)])


def test_mesh_source_zero_voxel_gets_zero_strength_uniform(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_openmc(monkeypatch)
    source = _source([[100.0, 0.0], [0.0, 0.0]])
    converted = omc.mesh_source(source)
    assert converted.sources[0].strength == pytest.approx(100.0)
    assert converted.sources[1].strength == 0.0
    assert isinstance(converted.sources[1].energy, _FakeUniform)


def test_mesh_source_rejects_nonregular_volumes(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_openmc(monkeypatch)
    source = _source(volumes=[1.0, 2.0])
    with pytest.raises(ValueError, match="regular-mesh cell volume"):
        omc.mesh_source(source)


# ---------------------------------------------------------------------------
# W4b: reproducibility manifest


def test_manifest_roundtrip_hash_and_tamper(tmp_path: Path) -> None:
    source = _source()
    input_file = tmp_path / "flux.csv"
    input_file.write_text("raw\n", encoding="utf-8")
    output_file = tmp_path / "source.xml"
    output_file.write_text("<settings/>\n", encoding="utf-8")
    manifest = omc.source_manifest(
        source=source,
        input_files=[input_file],
        output_files=[output_file],
        external_solver="openmc",
        external_version="0.16.0",
        evaluated_library="caller-library-v1",
    )
    path = tmp_path / "manifest.json"
    omc.write_manifest(manifest, path)
    first = path.read_bytes()
    omc.write_manifest(omc.read_manifest(path), path)
    assert path.read_bytes() == first
    parsed = omc.read_manifest(path)
    assert parsed["source"]["energy_grid_sha256"] == manifest["source"]["energy_grid_sha256"]
    assert parsed["source"]["voxel_order"] == "c-order-(i,j,k)"
    assert parsed["external"]["solver_version_source"] == "caller-supplied"
    assert parsed["evaluated_library"]["hash"] is None
    check = omc.verify_manifest(path)
    assert check["ok"] is True and check["checked"] == 2
    output_file.write_text("<tampered/>\n", encoding="utf-8")
    tampered = omc.verify_manifest(path)
    assert tampered["ok"] is False
    assert tampered["mismatches"][0]["path"].endswith("source.xml")


def test_manifest_missing_output_and_unknown_schema(tmp_path: Path) -> None:
    source = _source()
    manifest = omc.source_manifest(
        source=source,
        output_files=[tmp_path / "not-written.xml"],
    )
    check = omc.verify_manifest(manifest, base_dir=tmp_path)
    assert check["ok"] is False
    assert check["missing"][0]["path"].endswith("not-written.xml")
    with pytest.raises(ValueError, match="unknown manifest schema"):
        omc.write_manifest({"schema": "other"}, tmp_path / "bad.json")
    with pytest.raises(ValueError, match="not JSON-serializable"):
        bad = dict(manifest)
        bad["source"] = dict(bad["source"])
        bad["source"]["total_photons_per_s"] = float("nan")
        omc.write_manifest(bad, tmp_path / "bad-nan.json")


# ---------------------------------------------------------------------------
# W5: cooling-time batch export


def _writer(source: dict[str, Any], path: Path) -> None:
    path.write_text(
        json.dumps({"cooling_s": source["cooling_s"], "total": source["total_photons_per_s"]}),
        encoding="utf-8",
    )


def test_cooling_batch_exports_distinct_times_and_manifest(tmp_path: Path) -> None:
    zero = _source(cooling_s=0.0)
    hour = _source([[7.0, 11.0], [13.0, 17.0]], cooling_s=3600.0)
    batch = omc.export_cooling_batch([hour, zero], [0.0, 3600.0], tmp_path, writer=_writer)
    assert batch["manifest"]["status"] == "complete"
    assert [entry["time_s"] for entry in batch["outputs"]] == [0.0, 3600.0]
    assert [entry["total_photons_per_s"] for entry in batch["outputs"]] == pytest.approx(
        [400.0, 48.0]
    )
    for entry in batch["outputs"]:
        assert (tmp_path / entry["path"]).is_file()
    manifest = omc.read_manifest(tmp_path / "batch-manifest.json")
    assert manifest["schema"] == "nucleide-source-batch-v1"
    assert manifest["units"] == "photons/s"
    assert manifest["voxel_order"] == "c-order-(i,j,k)"
    assert manifest["outputs"][0]["validation"]["valid"] is True
    assert len(manifest["outputs"][0]["energy_grid_sha256"]) == 64
    assert omc.verify_manifest(manifest, base_dir=tmp_path)["ok"] is True


def test_cooling_batch_rejects_duplicates_missing_and_existing(tmp_path: Path) -> None:
    zero = _source(cooling_s=0.0)
    hour = _source(cooling_s=3600.0)
    with pytest.raises(ValueError):
        omc.export_cooling_batch([zero, hour], [3600.0, 3600.0], tmp_path, writer=_writer)
    with pytest.raises(ValueError):
        omc.export_cooling_batch([zero], [1.0], tmp_path, writer=_writer)
    omc.export_cooling_batch([zero], [0.0], tmp_path, writer=_writer)
    with pytest.raises(FileExistsError):
        omc.export_cooling_batch([zero], [0.0], tmp_path, writer=_writer)
    # Distinct names for close-but-distinct offsets.
    close_a = _source(cooling_s=3600.0)
    close_b = _source(cooling_s=3600.0000000000005)
    batch = omc.export_cooling_batch(
        [close_a, close_b],
        [3600.0, 3600.0000000000005],
        tmp_path / "close",
        writer=_writer,
    )
    assert len({entry["path"] for entry in batch["outputs"]}) == 2


def test_cooling_batch_zero_policy_and_partial_failure(tmp_path: Path) -> None:
    zero = _source([[0.0, 0.0], [0.0, 0.0]], cooling_s=0.0)
    live = _source(cooling_s=3600.0)
    with pytest.raises(ValueError, match="zero total"):
        omc.export_cooling_batch([zero], [0.0], tmp_path, writer=_writer)
    skipped = omc.export_cooling_batch(
        [zero, live], [0.0, 3600.0], tmp_path, writer=_writer, skip_zero=True
    )
    assert skipped["skipped"] == [{"time_s": 0.0, "reason": "all-zero-source"}]
    assert [entry["time_s"] for entry in skipped["outputs"]] == [3600.0]

    def failing_writer(source: dict[str, Any], path: Path) -> None:
        raise OSError("simulated disk failure")

    with pytest.raises(omc.BatchExportError) as caught:
        omc.export_cooling_batch([live], [3600.0], tmp_path / "fail", writer=failing_writer)
    assert caught.value.partial["status"] == "partial"
    assert caught.value.partial["outputs"] == []
    assert not (tmp_path / "fail" / "batch-manifest.json").exists()
