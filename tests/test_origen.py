"""Golden tests for the scoped ORIGEN TAPE readers against synthetic fixtures."""

from pathlib import Path

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "origen"


def read_fixture(name: str) -> str:
    return (FIX / name).read_text()


class TestTape5:
    def test_sample(self) -> None:
        tape = nucleide.origen.origen_parse_tape5(read_fixture("tape5_sample"))
        assert tape["titles"] == ["SYNTHETIC PWR PIN - TAPE5 SAMPLE", "CASE 1 - BASE DEPLETION"]
        assert tape["irradiation_steps"] == [
            {"flux": pytest.approx(3.0e13), "days": pytest.approx(100.0)},
            {"flux": pytest.approx(0.0), "days": pytest.approx(30.0)},
        ]
        assert [m["name"] for m in tape["materials"]] == ["fuel", "clad"]
        assert tape["materials"][0]["entries"] == [
            {"nuclide": "U235", "grams": pytest.approx(10.5)},
            {"nuclide": "U238", "grams": pytest.approx(1000.0)},
            {"nuclide": "Pu239", "grams": pytest.approx(0.25)},
        ]
        assert tape["materials"][1]["entries"] == [
            {"nuclide": "Zr90", "grams": pytest.approx(50.0)}
        ]

    def test_alias(self) -> None:
        assert nucleide.origen.tape5_parse is nucleide.origen.origen_parse_tape5

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.origen.origen_parse_tape5("  \n")

    def test_bad_flux_raises(self) -> None:
        with pytest.raises(ValueError, match="FLUX"):
            nucleide.origen.origen_parse_tape5("TITLE\nFLUX= hot DAYS= 10.0\n")


class TestTape6:
    def test_sample(self) -> None:
        tape = nucleide.origen.origen_parse_tape6(read_fixture("tape6_sample"))
        assert len(tape["records"]) == 4
        assert tape["records"][0] == {
            "nuclide": "U235",
            "grams": pytest.approx(10.2),
            "activity_bq": pytest.approx(8.16e5),
        }
        assert tape["records"][3]["nuclide"] == "Cs137"
        expected = 8.16e5 + 1.23e4 + 5.63e8 + 3.84e10
        assert tape["total_activity"] == pytest.approx(expected)

    def test_alias(self) -> None:
        assert nucleide.origen.tape6_parse is nucleide.origen.origen_parse_tape6

    def test_unknown_nuclide_raises(self) -> None:
        with pytest.raises(ValueError, match="Xx999"):
            nucleide.origen.origen_parse_tape6("U235 10.0 8.16e5\nXx999 1.0 2.0\n")

    def test_bad_grams_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.origen.origen_parse_tape6("U235 lots 8.16e5\n")


class TestTape9:
    def test_sample(self) -> None:
        entries = nucleide.origen.origen_parse_tape9(read_fixture("tape9_sample"))
        assert len(entries) == 5
        assert entries[0] == {"nuclide": "U235", "decay_const": pytest.approx(3.1209e-17)}
        assert entries[4] == {"nuclide": "Co60", "decay_const": pytest.approx(4.1674e-09)}
        by_name = {e["nuclide"]: e["decay_const"] for e in entries}
        assert by_name["Cs137"] == pytest.approx(7.3217e-10)

    def test_alias(self) -> None:
        assert nucleide.origen.tape9_parse is nucleide.origen.origen_parse_tape9

    def test_bad_lambda_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.origen.origen_parse_tape9("U235 3.1209e-17\nU238 fast\n")

    def test_negative_lambda_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.origen.origen_parse_tape9("U235 -1.0e-9\n")
