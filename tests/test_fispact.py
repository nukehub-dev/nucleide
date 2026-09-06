"""Golden tests for the fispact-io parser against the synthetic fixture."""

from pathlib import Path

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "fispact"

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


def read_fixture(name: str) -> str:
    return (FIX / name).read_text()


class TestParseOutput:
    def test_fixture_rows(self) -> None:
        rows = nucleide.fispact.fispact_parse_output(read_fixture("inventory.fis"), "synth")
        # 3 cooling steps x 4 nuclide rows (incl. `total`) x 3 variables.
        assert len(rows) == 36
        assert all(set(r) == ROW_KEYS for r in rows)
        assert all(r["run_lbl"] == "synth" for r in rows)
        assert {r["block"] for r in rows} == {"Material"}
        assert all(r["block_name"] == "inventory" for r in rows)
        assert all(r["block_num"] == -1 for r in rows)
        assert all(r["half_life_s"] == -1.0 for r in rows)
        assert {r["variable"] for r in rows} == {
            "Number Density",
            "Specific Activity",
            "Total Decay Heat",
        }
        times = sorted({r["time_s"] for r in rows})
        assert times == pytest.approx([0.0, 86_400.0, 31_536_000.0])

    def test_values_and_units(self) -> None:
        rows = nucleide.fispact.fispact_parse_output(read_fixture("inventory.fis"), "synth")
        h3_atoms = [
            r
            for r in rows
            if r["nuclide"] == "h-3"
            and r["variable"] == "Number Density"
            and r["time_s"] == pytest.approx(86_400.0)
        ]
        assert len(h3_atoms) == 1
        assert h3_atoms[0]["value"] == pytest.approx(9.99e19)
        assert h3_atoms[0]["var_unit"] == "atoms"
        h3_act = [
            r
            for r in rows
            if r["nuclide"] == "h-3"
            and r["variable"] == "Specific Activity"
            and r["time_s"] == pytest.approx(86_400.0)
        ][0]
        assert h3_act["value"] == pytest.approx(1.199e9)
        assert h3_act["var_unit"] == "Bq"
        totals = [r for r in rows if r["nuclide"] == "total"]
        assert len(totals) == 9

    def test_alias(self) -> None:
        assert nucleide.fispact.parse_output is nucleide.fispact.fispact_parse_output

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.fispact.fispact_parse_output("  \n", "empty")

    def test_bad_nuclide_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.fispact.fispact_parse_output("TIME = 0.0 SECS\nXx999 1.0 2.0 3.0\n", "r")

    def test_row_outside_step_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.fispact.fispact_parse_output("h-3 1.0 2.0 3.0\n", "r")
