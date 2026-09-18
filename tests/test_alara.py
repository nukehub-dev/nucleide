"""Golden tests for the alara-io readers against vendored fixtures."""

from pathlib import Path
from typing import Any

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "alara"

ROW_KEYS = {
    "time_s",
    "time_label",
    "nuclide",
    "half_life_s",
    "run_lbl",
    "block",
    "block_name",
    "block_num",
    "variable",
    "var_unit",
    "value",
}


def read_fixture(*parts: str) -> str:
    return (FIX.joinpath(*parts)).read_text()


class TestDeck:
    def test_sample2_blocks(self) -> None:
        deck = nucleide.alara.alara_parse_deck(read_fixture("decks", "sample2"))
        assert deck["block_kinds"] == [
            "geometry",
            "dimension",
            "mat_loading",
            "material_lib",
            "element_lib",
            "mixture",
            "mixture",
            "flux",
            "schedule",
            "pulsehistory",
            "dump_file",
            "data_library",
            "cooling",
            "output",
            "truncation",
        ]
        assert deck["geometry"] == "rectangular"
        assert [m["name"] for m in deck["mixtures"]] == ["inner_mix", "outer_mix"]
        assert deck["fluxes"][0]["name"] == "flux_1"
        assert deck["fluxes"][0]["file"] == "data/fluxin1"
        assert deck["fluxes"][0]["scale"] == pytest.approx(1.0)
        assert deck["fluxes"][0]["skip"] == 1
        assert deck["fluxes"][0]["format"] == "default"
        assert len(deck["cooling_times_s"]) == 4
        assert deck["cooling_times_s"][0] == pytest.approx(86_400.0)
        assert deck["cooling_times_s"][1] == pytest.approx(100.0 * 86_400.0)
        assert deck["outputs"][0]["resolution"] == "interval"
        assert deck["truncation"] == pytest.approx(1e-7)

    def test_sample3_blocks(self) -> None:
        deck = nucleide.alara.alara_parse_deck(read_fixture("decks", "sample3"))
        assert deck["geometry"] == "cylindrical"
        assert [m["name"] for m in deck["mixtures"]] == ["inner_mix", "outer_mix"]
        assert deck["fluxes"][0]["file"] == "data/fluxin2"
        assert deck["fluxes"][0]["scale"] == pytest.approx(1e6)
        assert deck["fluxes"][0]["skip"] == 0
        assert len(deck["cooling_times_s"]) == 4
        assert deck["cooling_times_s"][0] == pytest.approx(1.0)
        assert deck["outputs"][0]["resolution"] == "zone"
        assert deck["truncation"] == pytest.approx(1e-8)
        assert "schedule" in deck["block_kinds"]
        assert "pulsehistory" in deck["block_kinds"]

    def test_invalid_deck_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_parse_deck("bogus 1\n")
        with pytest.raises(ValueError):
            nucleide.alara.alara_parse_deck("   \n")


class TestFlux:
    def test_fluxin2_three_by_175(self) -> None:
        flux = nucleide.alara.alara_parse_flux(read_fixture("flux", "fluxin2"), "fluxin2")
        assert flux["name"] == "fluxin2"
        assert flux["groups_per_interval"] == 175
        assert flux["num_intervals"] == 3
        assert len(flux["intervals"]) == 3
        assert all(len(iv) == 175 for iv in flux["intervals"])
        assert len(flux["totals"]) == 3
        assert flux["total"] == pytest.approx(sum(flux["totals"]))
        assert flux["total"] > 0.0

    def test_fluxin_zeros(self) -> None:
        flux = nucleide.alara.alara_parse_flux(read_fixture("flux", "fluxin_zeros"), "zeros")
        assert flux["groups_per_interval"] == 175
        assert flux["num_intervals"] == 1
        assert flux["total"] == pytest.approx(0.0)

    def test_invalid_flux_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_parse_flux("1.0\nnope\n", "bad")


class TestOutput:
    def test_sample2_rows(self) -> None:
        rows = nucleide.alara.alara_parse_output(read_fixture("output", "sample2.out"), "sample2")
        assert len(rows) == 2838
        assert all(set(r) == ROW_KEYS for r in rows)
        variables = {r["variable"] for r in rows}
        assert variables == {"Number Density", "Specific Activity"}
        assert {r["block"] for r in rows} == {"Interval"}
        assert {(r["block_num"], r["block_name"]) for r in rows} == {
            (1, "inner_zone"),
            (2, "outer_zone"),
        }
        times = sorted({r["time_s"] for r in rows})
        assert times == [
            -1.0,
            0.0,
            86_400.0,
            8_640_000.0,
            31_536_000.0,
            3_153_600_000.0,
        ]
        assert all(r["run_lbl"] == "sample2" for r in rows)
        totals = [r for r in rows if r["nuclide"] == "total"]
        assert len(totals) == 24
        h1 = [
            r
            for r in rows
            if r["nuclide"] == "h-1"
            and r["variable"] == "Number Density"
            and r["block_num"] == 1
            and r["time_s"] == -1.0
        ]
        assert len(h1) == 1
        assert h1[0]["value"] == pytest.approx(6.6354e21)
        assert h1[0]["half_life_s"] == pytest.approx(-1.0)
        assert h1[0]["time_label"] == "pre-irrad"
        assert h1[0]["var_unit"] == "atoms/cm3"

    def test_empty_output_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_parse_output("  \n", "empty")


class TestSchedule:
    def test_expand_sample2_steady_state(self) -> None:
        steps = nucleide.alara.alara_expand_schedule(read_fixture("decks", "sample2"))
        assert len(steps) == 1
        assert steps[0]["flux"] == "flux_1"
        assert steps[0]["duration_s"] == pytest.approx(365.25 * 86_400.0)
        assert steps[0]["is_cooling"] is False

    def test_expand_sample3_pulsed(self) -> None:
        steps = nucleide.alara.alara_expand_schedule(read_fixture("decks", "sample3"), top="total")
        assert len(steps) == 19
        assert steps[0]["duration_s"] == pytest.approx(1.0)
        assert steps[1]["is_cooling"] is True
        assert steps[1]["duration_s"] == pytest.approx(5.0)
        total = sum(s["duration_s"] for s in steps)
        assert total == pytest.approx(10.0 * 1.0 + 9.0 * 5.0)

    def test_expand_unknown_top_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_expand_schedule(read_fixture("decks", "sample2"), top="missing")


def s1_entry(
    nuclide: str, activity_bq: float, irt: int, alpha_frac: float | None = None
) -> dict[str, Any]:
    entry: dict[str, Any] = {"nuclide": nuclide, "activity_bq": activity_bq, "irt": irt}
    if alpha_frac is not None:
        entry["alpha_frac"] = alpha_frac
    return entry


def s2_entry(
    nuclide: str, activity_bq: float, e_alpha: float, e_beta: float, e_gamma: float
) -> dict[str, Any]:
    return {
        "nuclide": nuclide,
        "activity_bq": activity_bq,
        "e_alpha_ev": e_alpha,
        "e_beta_ev": e_beta,
        "e_gamma_ev": e_gamma,
    }


class TestSubletS1:
    def test_hand_vectors_at_exact_equality(self) -> None:
        out = nucleide.alara.alara_total_activity([s1_entry("Co60", 10.0, 1)])
        assert out["total_bq"] == 10.0
        assert out["beta_bq"] == 10.0
        assert out["alpha_bq"] == 0.0
        assert out["gamma_bq"] == 0.0
        out = nucleide.alara.alara_total_activity(
            [
                s1_entry("Co60", 10.0, 1),
                s1_entry("Po210", 5.0, 4),
                s1_entry("Tc99_m1", 4.0, 3),
            ]
        )
        assert out["total_bq"] == 19.0
        assert out["alpha_bq"] == 5.0
        assert out["beta_bq"] == 10.0
        assert out["gamma_bq"] == 4.0

    def test_split_irts_and_conservation(self) -> None:
        out = nucleide.alara.alara_total_activity([s1_entry("U235", 100.0, 12, 0.25)])
        assert out["alpha_bq"] == 25.0
        assert out["beta_bq"] == 75.0
        out = nucleide.alara.alara_total_activity([s1_entry("Pu240", 10.0, 15, 0.5)])
        assert out["alpha_bq"] == 5.0
        assert out["gamma_bq"] == 5.0
        out = nucleide.alara.alara_total_activity(
            [
                s1_entry("Co60", 10.0, 1),
                s1_entry("Po210", 5.0, 4),
                s1_entry("U235", 100.0, 12, 0.25),
                s1_entry("Pu240", 10.0, 15, 0.5),
            ]
        )
        assert out["total_bq"] == out["alpha_bq"] + out["beta_bq"] + out["gamma_bq"]

    def test_ex_tritium(self) -> None:
        out = nucleide.alara.alara_total_activity(
            [s1_entry("Co60", 10.0, 1), s1_entry("H3", 5.0, 1)]
        )
        assert out["total_bq"] == 15.0
        assert out["ex_tritium_bq"] == 10.0

    def test_malformed_is_loud(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("Co60", -1.0, 1)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("Co60", 1.0, 10)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("U235", 1.0, 12)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("Co60", 1.0, 1, 0.5)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("U235", 1.0, 12, 1.5)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_total_activity([s1_entry("Xx999", 1.0, 1)])


class TestSubletS2:
    def test_hand_vectors_at_exact_equality(self) -> None:
        c1 = 1.602176634e-22
        out = nucleide.alara.alara_decay_heat([s2_entry("Co60", 1e10, 0.0, 1e6, 1e6)])
        assert out["beta_kw"] == 1e10 * 1e6 * c1
        assert out["gamma_kw"] == 1e10 * 1e6 * c1
        assert out["alpha_kw"] == 0.0
        assert out["total_kw"] == out["alpha_kw"] + out["beta_kw"] + out["gamma_kw"]
        out = nucleide.alara.alara_decay_heat(
            [
                s2_entry("Po210", 2e10, 5e6, 0.0, 0.0),
                s2_entry("Co60", 1e10, 0.0, 1e6, 2e6),
            ]
        )
        assert out["alpha_kw"] == 2e10 * 5e6 * c1
        assert out["beta_kw"] == 1e10 * 1e6 * c1
        assert out["gamma_kw"] == 1e10 * 2e6 * c1

    def test_ex_tritium(self) -> None:
        c1 = 1.602176634e-22
        out = nucleide.alara.alara_decay_heat(
            [
                s2_entry("Co60", 10.0, 0.0, 1e6, 0.0),
                s2_entry("H3", 5.0, 0.0, 6e3, 0.0),
            ]
        )
        assert out["total_kw"] == 10.0 * 1e6 * c1 + 5.0 * 6e3 * c1
        assert out["ex_tritium_kw"] == 10.0 * 1e6 * c1

    def test_malformed_is_loud(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_decay_heat([s2_entry("Co60", -1.0, 0.0, 1e6, 0.0)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_decay_heat([s2_entry("Co60", 1.0, 0.0, -2.0, 0.0)])
        with pytest.raises(ValueError):
            nucleide.alara.alara_decay_heat([s2_entry("Xx999", 1.0, 0.0, 1e6, 0.0)])
