"""Integration tests for the R2S workflow builder against ALARA fixtures."""

from pathlib import Path
from typing import Any

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "alara"


def read_fixture(*parts: str) -> str:
    return (FIX.joinpath(*parts)).read_text()


class TestFromDeck:
    def test_sample2_workflow(self) -> None:
        workflow = nucleide.r2s.r2s_from_deck(read_fixture("decks", "sample2"))
        assert workflow["steps"] == [
            {"zone": "inner_zone", "flux": "flux_1"},
            {"zone": "outer_zone", "flux": "flux_1"},
        ]
        assert len(workflow["cooling_s"]) == 4
        assert workflow["cooling_s"][0] == pytest.approx(86_400.0)
        assert workflow["top_schedule"] == "1_year"

    def test_aliases(self) -> None:
        assert nucleide.r2s.from_deck is nucleide.r2s.r2s_from_deck
        assert nucleide.r2s.validate is nucleide.r2s.r2s_validate
        assert nucleide.r2s.expand is nucleide.r2s.r2s_expand
        assert nucleide.r2s.assemble is nucleide.r2s.r2s_assemble

    def test_invalid_deck_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.r2s.r2s_from_deck("bogus 1\n")


class TestValidate:
    def test_sample2_validates(self) -> None:
        deck_text = read_fixture("decks", "sample2")
        workflow = nucleide.r2s.r2s_from_deck(deck_text)
        nucleide.r2s.r2s_validate(workflow, deck_text)

    def test_unknown_zone_raises(self) -> None:
        deck_text = read_fixture("decks", "sample2")
        workflow = {
            "steps": [{"zone": "ghost_zone", "flux": "flux_1"}],
            "cooling_s": [86_400.0],
            "top_schedule": "1_year",
        }
        with pytest.raises(ValueError, match="ghost_zone"):
            nucleide.r2s.r2s_validate(workflow, deck_text)

    def test_unknown_flux_raises(self) -> None:
        deck_text = read_fixture("decks", "sample2")
        workflow = {
            "steps": [{"zone": "inner_zone", "flux": "ghost_flux"}],
            "cooling_s": [86_400.0],
            "top_schedule": "1_year",
        }
        with pytest.raises(ValueError, match="ghost_flux"):
            nucleide.r2s.r2s_validate(workflow, deck_text)


class TestExpand:
    def test_sample2_steady_state(self) -> None:
        steps = nucleide.r2s.r2s_expand(read_fixture("decks", "sample2"))
        assert len(steps) == 1
        assert steps[0]["flux"] == "flux_1"
        assert steps[0]["duration_s"] == pytest.approx(365.25 * 86_400.0)
        assert steps[0]["is_cooling"] is False

    def test_sample3_pulsed(self) -> None:
        steps = nucleide.r2s.r2s_expand(read_fixture("decks", "sample3"), top="total")
        assert len(steps) == 19
        assert steps[0]["duration_s"] == pytest.approx(1.0)
        assert steps[1]["is_cooling"] is True

    def test_unknown_top_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.r2s.r2s_expand(read_fixture("decks", "sample2"), top="missing")


class TestFromSnapshot:
    @staticmethod
    def snapshot(**overrides: Any) -> dict[str, Any]:
        base: dict[str, Any] = {
            "zones": [
                {
                    "id": "B0120-005",
                    "volume_cm3": 1200.0,
                    "zbottom_cm": 0.0,
                    "ztop_cm": 10.0,
                    "material": "fuel",
                    "xs_type": "A",
                    "temperature_C": 600.0,
                    "composition": {"U235": 1.0e-3, "nU238": 2.0e-2},
                },
                {
                    "id": "B0120-006",
                    "volume_cm3": 800.0,
                    "composition": {"PU239": 5.0e-4},
                },
            ],
            "flux_defs": [{"name": "snap_flux", "file": "data/flux", "scale": 1.0}],
            "cooling_s": [86_400.0],
        }
        base.update(overrides)
        return base

    def test_broadcast_workflow_and_decks(self) -> None:
        bundle = nucleide.r2s.r2s_from_snapshot(self.snapshot())
        workflow = bundle["workflow"]
        assert workflow["steps"] == [
            {"zone": "B0120-005", "flux": "snap_flux"},
            {"zone": "B0120-006", "flux": "snap_flux"},
        ]
        assert workflow["cooling_s"] == pytest.approx([86_400.0])
        assert workflow["top_schedule"] == "snap_schedule"
        assert "mixture mix_B0120-005" in bundle["deck"]
        assert len(bundle["decks"]) == 2
        # Template deck validates and expands through the existing machinery.
        nucleide.r2s.r2s_validate(workflow, bundle["deck"])
        steps = nucleide.r2s.r2s_expand(bundle["deck"])
        assert len(steps) == 1
        assert steps[0]["duration_s"] == pytest.approx(86_400.0)

    def test_alias(self) -> None:
        assert nucleide.r2s.from_snapshot is nucleide.r2s.r2s_from_snapshot

    def test_void_zone_skipped(self) -> None:
        snap = self.snapshot()
        snap["zones"][1]["composition"] = {}
        bundle = nucleide.r2s.r2s_from_snapshot(snap)
        assert bundle["workflow"]["steps"] == [{"zone": "B0120-005", "flux": "snap_flux"}]
        assert len(bundle["decks"]) == 1

    def test_name_matched_fluxes(self) -> None:
        snap = self.snapshot(
            flux_defs=[
                {"name": "B0120-005", "file": "data/f1", "scale": 1.0},
                {"name": "B0120-006", "file": "data/f2", "scale": 2.0},
            ]
        )
        bundle = nucleide.r2s.r2s_from_snapshot(snap)
        assert bundle["workflow"]["steps"] == [
            {"zone": "B0120-005", "flux": "B0120-005"},
            {"zone": "B0120-006", "flux": "B0120-006"},
        ]

    def test_empty_zones_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.r2s.r2s_from_snapshot(self.snapshot(zones=[]))

    def test_unknown_composition_key_raises(self) -> None:
        snap = self.snapshot()
        snap["zones"][0]["composition"] = {"DUMP1": 1.0e-3}
        with pytest.raises(ValueError, match="DUMP1"):
            nucleide.r2s.r2s_from_snapshot(snap)

    def test_empty_cooling_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.r2s.r2s_from_snapshot(self.snapshot(cooling_s=[]))


class TestAssemble:
    def test_sample2_inner_zone(self) -> None:
        output_text = read_fixture("output", "sample2.out")
        source = nucleide.r2s.r2s_assemble(output_text, "sample2", "inner_zone", 4)
        assert source["zone"] == "inner_zone"
        assert len(source["groups"]) == 4
        # Uniform split: every group carries one quarter of the total.
        assert source["groups"][0] == pytest.approx(source["total"] / 4.0)
        assert all(g == pytest.approx(source["groups"][0]) for g in source["groups"])
        assert source["total"] > 0.0

    def test_missing_zone_is_zero(self) -> None:
        output_text = read_fixture("output", "sample2.out")
        source = nucleide.r2s.r2s_assemble(output_text, "sample2", "missing", 3)
        assert source["groups"] == pytest.approx([0.0, 0.0, 0.0])
        assert source["total"] == pytest.approx(0.0)

    def test_zero_groups(self) -> None:
        output_text = read_fixture("output", "sample2.out")
        source = nucleide.r2s.r2s_assemble(output_text, "sample2", "inner_zone", 0)
        assert source["groups"] == []
        assert source["total"] == pytest.approx(0.0)
