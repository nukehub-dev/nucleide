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
