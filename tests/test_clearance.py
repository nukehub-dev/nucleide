"""Clearance / waste-classification analytics tests.

Covers the `alara-io` clearance kernel facade (clearance index + sum-of-
fractions over caller-supplied or EU-default limit tables) and the
`fispact-io` clearance-block reader against the synthetic fixture in the real
FISPACT-II wide-table grammar.
"""

from pathlib import Path

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures" / "fispact"

TOY_LIMITS = {"Co60": 10.0, "H3": 5.0, "Fe55": 2.0}


class TestEuTable:
    def test_spot_values_from_directive(self) -> None:
        table = nucleide.alara.alara_clearance_eu_table()
        assert len(table) > 200
        # Annex VII Table A transcription spot checks (Bq/g).
        assert table["H-3"] == 100.0
        assert table["C-14"] == 1.0
        assert table["Co-60"] == 0.1
        assert table["Cs-137"] == 0.1
        assert table["Sr-90"] == 1.0
        assert table["Pu-239"] == 0.1
        # Part 2 parents + K-40.
        assert table["U-238"] == 1.0
        assert table["Th-232"] == 1.0
        assert table["K-40"] == 10.0
        assert all(limit > 0.0 for limit in table.values())


class TestEsCsnTables:
    """Spanish CSN conditional NORM landfill tables (Tablas 1-3).

    Spots are hand-read from CSN/PDT/AICD/TGE/2503/02 (TGE/VAR/2025/1),
    flattened per the documented chain-expansion rule (last claiming key
    wins, so subchain levels surface on their members).
    """

    def test_entry_counts_and_positive_limits(self) -> None:
        for landfill in ("inert", "non_hazardous", "hazardous"):
            table = nucleide.alara.alara_clearance_es_table(landfill, "rocks")
            assert len(table) == 40
            assert all(limit > 0.0 for limit in table.values())

    def test_inert_rocks_spots(self) -> None:
        table = nucleide.alara.alara_clearance_es_table("inert", "rocks")
        assert table["U-238"] == 10.0  # U-nat
        assert table["Pa-234m"] == 10.0  # U-nat
        assert table["Th-230"] == 10.0
        assert table["Ra-226"] == 10.0  # Ra-226+
        assert table["Bi-214"] == 10.0  # Ra-226+
        assert table["Pb-210"] == 10.0  # Pb-210+
        assert table["Po-210"] == 5.0
        assert table["Ac-227"] == 5.0  # Ac-227+
        assert table["Tl-207"] == 5.0  # Ac-227+
        assert table["Th-232"] == 5.0
        assert table["Ra-228"] == 10.0  # Ra-228+
        assert table["Th-228"] == 5.0  # Th-228+
        assert table["K-40"] == 10.0

    def test_material_columns_differ(self) -> None:
        inert_rocks = nucleide.alara.alara_clearance_es_table("inert", "rocks")
        inert_gas = nucleide.alara.alara_clearance_es_table("inert", "oil_gas")
        assert inert_gas["U-238"] == 500.0
        assert inert_gas["Ra-226"] == 50.0
        assert inert_gas["Pb-210"] == 100.0
        assert inert_gas["Th-232"] == 100.0
        assert inert_gas["K-40"] == 100.0
        assert inert_rocks["U-238"] != inert_gas["U-238"]

    def test_landfill_types_differ(self) -> None:
        rocks = {
            landfill: nucleide.alara.alara_clearance_es_table(landfill, "rocks")
            for landfill in ("inert", "non_hazardous", "hazardous")
        }
        assert rocks["inert"]["Po-210"] == 5.0
        assert rocks["non_hazardous"]["Po-210"] == 10.0
        assert rocks["hazardous"]["Po-210"] == 100.0
        assert rocks["hazardous"]["K-40"] == 50.0
        assert rocks["hazardous"]["U-238"] == 10.0
        haz_gas = nucleide.alara.alara_clearance_es_table("hazardous", "oil_gas")
        assert haz_gas["U-238"] == 500.0
        assert haz_gas["Po-210"] == 500.0
        assert haz_gas["K-40"] == 500.0

    def test_chain_expansion_parent_value_on_every_member(self) -> None:
        # Tabla 4 subchains whose members no later key claims keep the parent
        # value on every member (inert landfill, rocks column).
        table = nucleide.alara.alara_clearance_es_table("inert", "rocks")
        for member in ("Ra-226", "Rn-222", "Po-218", "Pb-214", "Bi-214", "Po-214"):
            assert table[member] == 10.0  # Ra-226+
        for member in ("Pb-210", "Bi-210"):
            assert table[member] == 10.0  # Pb-210+
        for member in ("Ac-227", "Th-227", "Fr-223", "Ra-223", "Rn-219", "Po-215"):
            assert table[member] == 5.0  # Ac-227+
        for member in ("Th-228", "Ra-224", "Rn-220", "Po-216", "Pb-212", "Bi-212"):
            assert table[member] == 5.0  # Th-228+
        assert table["Tl-208"] == 5.0
        assert table["Po-212"] == 5.0
        assert table["K-40"] == 10.0

    def test_hand_computed_screening_vectors_exact(self) -> None:
        inert_rocks = nucleide.alara.alara_clearance_es_table("inert", "rocks")
        assert (
            nucleide.alara.alara_clearance_index(
                {"U-238": 5.0, "Th-230": 5.0, "Po-210": 2.5, "K-40": 5.0}, inert_rocks
            )
            == 2.0
        )
        non_haz_sands = nucleide.alara.alara_clearance_es_table("non_hazardous", "sands")
        assert (
            nucleide.alara.alara_clearance_index(
                {"Pb-210": 5.0, "Bi-210": 5.0, "K-40": 5.0, "Ra-228": 10.0},
                non_haz_sands,
            )
            == 2.5
        )
        haz_gas = nucleide.alara.alara_clearance_es_table("hazardous", "oil_gas")
        assert (
            nucleide.alara.alara_clearance_index(
                {"U-238": 250.0, "Ra-226": 25.0, "Po-210": 250.0, "K-40": 250.0},
                haz_gas,
            )
            == 2.0
        )

    def test_bad_selector_is_loud(self) -> None:
        with pytest.raises(ValueError, match="landfill"):
            nucleide.alara.alara_clearance_es_table("bogus", "rocks")  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="material"):
            nucleide.alara.alara_clearance_es_table("inert", "lava")  # type: ignore[arg-type]


class TestClearanceIndex:
    def test_hand_computed_vectors_exact(self) -> None:
        assert nucleide.alara.alara_clearance_index({"Co60": 10.0}, TOY_LIMITS) == 1.0
        assert nucleide.alara.alara_clearance_index({"Co60": 5.0, "H-3": 2.5}, TOY_LIMITS) == 1.0
        inv = {"Co60": 10.0, "H3": 5.0, "Fe-55": 4.0}
        assert nucleide.alara.alara_clearance_index(inv, TOY_LIMITS) == 4.0
        assert nucleide.alara.alara_clearance_index({}, TOY_LIMITS) == 0.0
        assert nucleide.alara.alara_clearance_index({"Fe55": 0.0}, TOY_LIMITS) == 0.0

    def test_eu_default_table_matches_units(self) -> None:
        # 0.05 Bq/g of Co-60 against a 0.1 Bq/g limit -> CI == 0.5.
        assert nucleide.alara.alara_clearance_index({"Co-60": 0.05}) == pytest.approx(0.5)

    def test_missing_limit_is_loud(self) -> None:
        with pytest.raises(ValueError, match="Mn54"):
            nucleide.alara.alara_clearance_index({"Mn54": 1.0}, TOY_LIMITS)

    def test_bad_activity_is_loud(self) -> None:
        for bad in (-1.0, float("nan"), float("inf")):
            with pytest.raises(ValueError, match="Co60"):
                nucleide.alara.alara_clearance_index({"Co60": bad}, TOY_LIMITS)

    def test_bad_limit_is_loud(self) -> None:
        with pytest.raises(ValueError):
            nucleide.alara.alara_clearance_index({"Co60": 1.0}, {"Co60": 0.0})
        with pytest.raises(ValueError):
            nucleide.alara.alara_clearance_index({"Co60": 1.0}, {"Co60": -2.0})

    def test_bad_name_is_loud(self) -> None:
        with pytest.raises(ValueError, match="Xx999"):
            nucleide.alara.alara_clearance_index({"Xx999": 1.0}, TOY_LIMITS)


class TestSumOfFractions:
    def test_boundary_probes_both_sides(self) -> None:
        at = nucleide.alara.alara_sum_of_fractions({"Co60": 10.0}, TOY_LIMITS)
        assert at["sum"] == 1.0
        assert at["class"] == "satisfied"
        assert at["max_fraction"] == 1.0
        assert at["max_nuclide"] == "Co-60"

        above = nucleide.alara.alara_sum_of_fractions({"Co60": 10.000000000000002}, TOY_LIMITS)
        assert above["sum"] > 1.0
        assert above["class"] == "exceeded"

        below = nucleide.alara.alara_sum_of_fractions({"Co60": 9.999999999999998}, TOY_LIMITS)
        assert below["sum"] < 1.0
        assert below["class"] == "satisfied"

    def test_dominant_nuclide(self) -> None:
        out = nucleide.alara.alara_sum_of_fractions(
            {"Co60": 6.0, "H3": 1.0, "Fe55": 0.2}, TOY_LIMITS
        )
        assert out["class"] == "satisfied"
        assert out["sum"] == pytest.approx(0.9)
        assert out["max_nuclide"] == "Co-60"
        assert out["max_fraction"] == pytest.approx(0.6)


class TestFispactClearanceReader:
    def read_fixture(self) -> str:
        return (FIX / "clearance.out").read_text()

    def test_fixture_rows(self) -> None:
        rows = nucleide.fispact.fispact_parse_clearance(self.read_fixture())
        assert len(rows) == 7
        keys = {
            "interval",
            "time_s",
            "time_label",
            "cooling",
            "nuclide",
            "flags",
            "activity_bq",
            "clearance_index",
            "half_life_s",
        }
        assert all(set(r) == keys for r in rows)

        step1 = [r for r in rows if r["interval"] == 1]
        assert len(step1) == 4
        assert not any(r["cooling"] for r in step1)
        assert step1[0]["nuclide"] == "V-55"
        assert step1[0]["flags"] == ">"
        assert step1[0]["half_life_s"] == -1.0
        assert step1[1]["flags"] == "#>"

        co60 = step1[2]
        assert co60["nuclide"] == "Co-60"
        assert co60["activity_bq"] == pytest.approx(4.0e6)
        assert co60["clearance_index"] == pytest.approx(1.0e6)

        rb86m = step1[3]
        assert rb86m["nuclide"] == "Rb-86m"
        assert rb86m["flags"] == "&"

        step2 = [r for r in rows if r["interval"] == 2]
        assert len(step2) == 3
        assert all(r["cooling"] for r in step2)
        assert step2[0]["time_s"] == pytest.approx(86_400.0)
        assert step2[0]["time_label"] == "COOLING TIME IS 8.6400E+04 SECS OR 1.0000E+00 DAYS"

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.fispact.fispact_parse_clearance("  \n")

    def test_malformed_row_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.fispact.fispact_parse_clearance(
                "1 * * * TIME INTERVAL 1 * * * TIME IS 0.0 SECS\n"
                " NUCLIDE ATOMS GRAMS Bq b a g DOSE ING INH CLEARANCE Bq/A2 HALF LIFE\n"
                " kW kW kW Sv/hr DOSE(Sv) DOSE(Sv) INDEX Ratio seconds\n"
                "Co 60 > 1.0 1.0 notafloat 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n"
            )

    def test_end_to_end_screening(self) -> None:
        # Parse the clearance block, take the shutdown-step inventory, and
        # screen it against a caller-supplied limit table (zero-activity
        # nuclides still need table entries).
        rows = nucleide.fispact.fispact_parse_clearance(self.read_fixture())
        step1 = {r["nuclide"]: r["activity_bq"] for r in rows if r["interval"] == 1}
        limits = {"V-55": 1.0, "Fe-56": 1.0, "Co-60": 4.0e6, "Rb-86m": 8.0e9}
        out = nucleide.alara.alara_sum_of_fractions(step1, limits)
        assert out["sum"] == pytest.approx(2.0)
        assert out["class"] == "exceeded"
        # Both contributors sit at fraction 1.0; the tie keeps the first.
        assert out["max_nuclide"] == "Co-60"
        assert out["max_fraction"] == pytest.approx(1.0)
