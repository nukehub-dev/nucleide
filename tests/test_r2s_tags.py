"""Python-side tests for R2S per-voxel photon-source tags (synthetic only).

The `.photonSrc` text below is hand-built in the documented
`<nuclide> <time> <strengths...>` row format; all expected values are
closed-form mapping math (no transport, no MOAB, no HDF5).
"""

import pytest

import nucleide.r2s as r2s

PHOTON_TEXT = """\
mn-56 shutdown 6.0 2.0
co-60 shutdown 1.0 3.0
mn-56 1 h 5.0 1.0
TOTAL shutdown 7.0 5.0
"""


def test_tag_copies_zone_totals() -> None:
    tags = r2s.r2s_tag_zone_strength([400.0, 50.0], [0, 0, 1])
    assert tags["n_zones"] == 2
    assert tags["zone_of_voxel"] == [0, 0, 1]
    assert tags["source_strength"] == pytest.approx([400.0, 400.0, 50.0])
    assert tags["decay_time_s"] == pytest.approx([0.0, 0.0, 0.0])
    assert tags["total"] == pytest.approx(850.0)


def test_tag_splits_conservatively() -> None:
    tags = r2s.r2s_tag_zone_strength([400.0, 50.0], [0, 0, 1], split=True)
    assert tags["source_strength"] == pytest.approx([200.0, 200.0, 50.0])
    assert tags["total"] == pytest.approx(450.0)


def test_tag_rejects_dangling_zone() -> None:
    with pytest.raises(ValueError):
        r2s.r2s_tag_zone_strength([400.0], [0, 7])


def test_photon_group_sums_select_and_add() -> None:
    out = r2s.r2s_photon_group_sums(PHOTON_TEXT, ["mn-56", "co-60"], 0.0)
    assert [g["nuclide"] for g in out["groups"]] == ["mn-56", "co-60"]
    assert out["sums"] == pytest.approx([7.0, 5.0])
    assert out["total"] == pytest.approx(12.0)


def test_photon_group_sums_time_and_total_rows() -> None:
    out = r2s.r2s_photon_group_sums(PHOTON_TEXT, ["mn-56"], 3600.0)
    assert len(out["groups"]) == 1
    assert out["sums"] == pytest.approx([5.0, 1.0])
    out = r2s.r2s_photon_group_sums(PHOTON_TEXT, ["ghost"], 0.0)
    assert out["groups"] == []
    assert out["sums"] == []
    assert out["total"] == pytest.approx(0.0)


def test_volume_weighted_integrated_tags() -> None:
    tags = r2s.tag_zone_strength(
        [120.0, 40.0], [0, 0, 1], split=True, voxel_volumes_cm3=[1.0, 3.0, 2.0]
    )
    assert tags["source_strength"] == pytest.approx([30.0, 90.0, 40.0])
    assert tags["total"] == pytest.approx(160.0)
    assert tags["decay_time_s"] == [0.0, 0.0, 0.0]


@pytest.mark.parametrize(  # type: ignore[untyped-decorator]
    "volumes", [[1.0], [0.0, 1.0], [-1.0, 1.0], [float("inf"), 1.0]]
)
def test_volume_tag_rejects_invalid_volumes(volumes: list[float]) -> None:
    with pytest.raises(ValueError):
        r2s.tag_zone_strength([12.0], [0, 0], split=True, voxel_volumes_cm3=volumes)


def test_volumes_require_split_and_mapped_positive_sources() -> None:
    with pytest.raises(ValueError, match="split"):
        r2s.tag_zone_strength([12.0], [0], voxel_volumes_cm3=[2.0])
    with pytest.raises(ValueError):
        r2s.tag_zone_strength([12.0, 3.0], [0], split=True, voxel_volumes_cm3=[2.0])


def test_nonunit_photon_density_readback() -> None:
    text = "TOTAL 600 s 8 4\nTOTAL 3600 s 2 4\n"
    source = r2s.photon_group_sums(text, ["TOTAL"], 3600.0)
    assert source["sums"] == [2.0, 4.0]
    assert source["total"] == 6.0  # photons/(s cm³)
    tags = r2s.tag_zone_strength(
        [source["total"] * 2.0], [0, 0], split=True, voxel_volumes_cm3=[0.5, 1.5]
    )
    assert tags["source_strength"] == [3.0, 9.0]
    assert tags["total"] == 12.0  # photons/s


def test_volume_total_retains_many_small_contributions() -> None:
    import math

    count = 100_000
    tiny = 2.0**-54
    volumes = [1.0] + [tiny] * count
    tags = r2s.tag_zone_strength([1.0], [0] * len(volumes), split=True, voxel_volumes_cm3=volumes)
    expected_volume = 1.0 + count * tiny
    assert tags["source_strength"][0] == 1.0 / expected_volume
    assert math.isclose(math.fsum(tags["source_strength"]), 1.0, rel_tol=1e-15, abs_tol=0.0)
    assert tags["total"] == math.fsum(tags["source_strength"])
