"""Golden tests for the cccc-io readers/writer against synthetic fixtures."""

from pathlib import Path
from typing import Any

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "cccc"


def read_fixture(name: str) -> str:
    return (FIX / name).read_text()


def sample_deck() -> dict[str, Any]:
    return {
        "title": "synthetic slab",
        "dim": 1,
        "zones": [
            {"id": 1, "material": "fuel", "isotxs_labels": ["U235"], "density": 10.0},
            {"id": 2, "material": "blanket", "isotxs_labels": ["PU239"], "density": 5.0},
        ],
        "source": "isotropic",
    }


class TestIsotxs:
    def test_fixture_two_by_three(self) -> None:
        lib = nucleide.cccc.isotxs_parse(read_fixture("isotxs_sample"))
        assert [n["label"] for n in lib["nuclides"]] == ["U235", "PU239"]
        u5 = lib["nuclides"][0]
        assert u5["zaid"] == "92235"
        assert u5["groups"] == 3
        assert u5["total_xs"] == [1.1, 2.2, 3.3]
        pu9 = lib["nuclides"][1]
        assert pu9["total_xs"] == pytest.approx([4.4, 5.5, 6.6])

    def test_alias(self) -> None:
        assert nucleide.cccc.cccc_parse_isotxs is nucleide.cccc.isotxs_parse

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.cccc.isotxs_parse("  \n")

    def test_group_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.cccc.isotxs_parse("ISOTXS 3\nNUCLIDE U235 92235 2\n1.0 2.0\n")


class TestRtflux:
    def test_fixture_two_points_three_groups(self) -> None:
        flux = nucleide.cccc.rtflux_parse(read_fixture("rtflux_sample"))
        assert flux["kind"] == "RTFLUX"
        assert flux["groups"] == 3
        assert flux["per_point"] == 3
        assert flux["npoints"] == 2
        assert flux["values"] == pytest.approx([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        assert flux["total"] == pytest.approx(21.0)

    def test_kind_variants(self) -> None:
        for kind, keyword in [("rtflux", "RTFLUX"), ("atflux", "ATFLUX"), ("rzflux", "RZFLUX")]:
            flux = nucleide.cccc.rtflux_parse(f"{keyword} 1 2\n1.0 2.0\n", kind=kind)
            assert flux["kind"] == keyword
            assert flux["npoints"] == 1

    def test_alias(self) -> None:
        assert nucleide.cccc.cccc_parse_rtflux is nucleide.cccc.rtflux_parse

    def test_kind_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.cccc.rtflux_parse("ATFLUX 1 2\n1.0 2.0\n", kind="rtflux")

    def test_bad_kind_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.cccc.rtflux_parse("RTFLUX 1 2\n1.0 2.0\n", kind="nope")


class TestPartisn:
    def test_render(self) -> None:
        text = nucleide.cccc.partisn_render(sample_deck())
        assert "synthetic slab" in text
        assert "DIM 1" in text
        assert "ZONE 1" in text
        assert "MATERIAL fuel" in text
        assert "U235" in text
        assert "SOURCE isotropic" in text
        assert text.endswith("END\n")

    def test_render_without_source(self) -> None:
        deck = sample_deck()
        deck["source"] = None
        text = nucleide.cccc.partisn_render(deck)
        assert "SOURCE" not in text
        assert text.endswith("END\n")

    def test_validate(self) -> None:
        nucleide.cccc.partisn_validate(sample_deck(), read_fixture("isotxs_sample"))

    def test_validate_aliases(self) -> None:
        assert nucleide.cccc.cccc_render_partisn is nucleide.cccc.partisn_render
        assert nucleide.cccc.cccc_validate_partisn is nucleide.cccc.partisn_validate

    def test_dangling_label_raises(self) -> None:
        deck = sample_deck()
        zones = deck["zones"]
        assert isinstance(zones, list)
        zones[0]["isotxs_labels"].append("U238")
        with pytest.raises(ValueError, match="U238"):
            nucleide.cccc.partisn_validate(deck, read_fixture("isotxs_sample"))

    def test_bad_dim_raises(self) -> None:
        deck = sample_deck()
        deck["dim"] = 5
        with pytest.raises(ValueError):
            nucleide.cccc.partisn_validate(deck, read_fixture("isotxs_sample"))
