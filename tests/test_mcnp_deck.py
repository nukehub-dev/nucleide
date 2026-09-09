"""Golden tests for full-deck parse/edit/write round-trip (synthetic fixtures)."""

from pathlib import Path

import pytest

import nucleide

FIXTURES = Path(__file__).parent.parent / "fixtures" / "mcnp" / "inp"


def _deck(name: str) -> "nucleide.mcnp.DeckProblem":
    return nucleide.mcnp.read_deck(str(FIXTURES / name))


class TestDeckRoundTrip:
    def test_byte_identical(self) -> None:
        for name in ["deck_minimal.txt", "deck_macro.txt", "deck_params.txt"]:
            path = FIXTURES / name
            deck = nucleide.mcnp.read_deck(str(path))
            assert deck.dumps() == path.read_text(), name

    def test_parse_deck_text_matches_file(self) -> None:
        path = FIXTURES / "deck_minimal.txt"
        assert nucleide.mcnp.parse_deck(path.read_text()).dumps() == path.read_text()


class TestDeckModel:
    def test_minimal_contents(self) -> None:
        deck = _deck("deck_minimal.txt")
        assert deck.title == "Minimal pin-cell deck"
        assert [c["num"] for c in deck.cells] == ["1", "2", "3"]
        assert deck.cells[0]["mat"] == "1"
        assert deck.cells[0]["dens"] == "-10"  # Rust shortest-render floats
        assert deck.cells[0]["geom"] == "-1"
        assert deck.cells[2]["dens"] == ""  # void cell
        assert [s["kind"] for s in deck.surfs] == ["SO", "SO", "PX"]
        assert deck.material_numbers == [1, 2]
        assert "MODE" in deck.data_names

    def test_macro_kinds_and_reflecting(self) -> None:
        deck = _deck("deck_macro.txt")
        kinds = {s["kind"] for s in deck.surfs}
        assert {"SO", "RPP", "RCC"} <= kinds
        assert deck.surfs[0]["reflecting"] == "true"  # Rust bool rendering

    def test_setters(self) -> None:
        deck = _deck("deck_minimal.txt")
        deck.set_cell_density(1, -7.0)
        assert deck.cells[0]["dens"] == "-7"
        assert "1 1 -7 -1 imp:n=1" in deck.dumps()
        deck.set_cell_material(2, 1)
        assert deck.cells[1]["mat"] == "1"
        with pytest.raises(ValueError, match="no cell"):
            deck.set_cell_density(99, -1.0)
        with pytest.raises(ValueError):
            deck.set_cell_density(3, -1.0)  # void cell

    def test_structural_errors(self) -> None:
        with pytest.raises(ValueError):
            nucleide.mcnp.parse_deck("only a message")
        with pytest.raises(ValueError):
            nucleide.mcnp.parse_deck("msg\ntitle\n1 1 -1.0 -1\n")


class TestDeckL3Semantics:
    def test_l3_fixture_round_trip_and_views(self) -> None:
        path = FIXTURES / "deck_l3.txt"
        deck = nucleide.mcnp.read_deck(str(path))
        assert deck.dumps() == path.read_text()

        assert deck.mode == {"particles": "N P"}
        transforms = deck.transforms
        assert [t["number"] for t in transforms] == ["1", "2"]
        assert transforms[0]["displacement"] == "0 0 5"
        assert transforms[0]["rotation"] == ""
        assert transforms[0]["in_degrees"] == "false"
        assert transforms[1]["in_degrees"] == "true"
        assert transforms[1]["main_to_aux"] == "true"

        universes = {u["number"]: u for u in deck.universes}
        assert set(universes) == {"0", "1", "2"}
        assert universes["1"]["cells"] == "1 2"
        assert universes["0"]["cells"] == "4"

        assert deck.lattices == [{"cell": "3", "lattice": "1"}]
        fills = deck.fills
        assert len(fills) == 1
        assert fills[0]["cell"] == "3"
        assert fills[0]["kind"] == "matrix"
        assert fills[0]["min_index"] == "0 0 0"
        assert fills[0]["max_index"] == "1 0 0"
        assert fills[0]["universes"] == "1 1"

        assert {v["particle"] for v in deck.importances} == {"N"}
        assert deck.volumes == [{"cell": "4", "volume": "100"}]

        tallies = deck.tallies
        assert len(tallies) == 1
        assert tallies[0]["number"] == "4"
        assert tallies[0]["type"] == "4"
        assert tallies[0]["particles"] == "N"
        assert tallies[0]["entries"] == "1 2"
        assert tallies[0]["fm"] == "1.0 1 101"
        assert tallies[0]["e_bins"] == "0.01 0.1 1.0"

        deck.validate()
        assert deck.validation_notes() == []

    def test_validate_catches_link_errors(self) -> None:
        deck = _deck("deck_l3.txt")
        deck.set_cell_fill(4, 9)
        with pytest.raises(ValueError, match="missing universe 9"):
            deck.validate()

    def test_validate_catches_duplicates_and_bad_cards(self) -> None:
        dup = "msg\ntitle\n1 0 -1\n1 0 -2\n\n1 so 1\n2 so 2\n\n"
        with pytest.raises(ValueError, match="duplicate cell number 1"):
            nucleide.mcnp.parse_deck(dup).validate()
        dup_mode = "msg\ntitle\n1 0 -1\n\n1 so 1\n\nmode n\nmode p\n"
        with pytest.raises(ValueError, match="duplicate MODE card"):
            nucleide.mcnp.parse_deck(dup_mode).validate()
        bad_lat = "msg\ntitle\n1 0 -1 lat=3\n\n1 so 1\n\n"
        with pytest.raises(ValueError, match="LAT must be 1 or 2"):
            nucleide.mcnp.parse_deck(bad_lat).validate()
        lat_no_fill = "msg\ntitle\n1 0 -1 lat=1\n\n1 so 1\n\n"
        with pytest.raises(ValueError, match="LAT but no FILL"):
            nucleide.mcnp.parse_deck(lat_no_fill).validate()
        dangling = "msg\ntitle\n1 5 -1.0 -1\n\n1 so 1\n\n"
        with pytest.raises(ValueError, match="missing material 5"):
            nucleide.mcnp.parse_deck(dangling).validate()

    def test_mode_mismatch_is_a_note(self) -> None:
        deck = _deck("deck_l3.txt")
        deck.set_mode(["P"])
        deck.validate()  # notes are not errors
        notes = deck.validation_notes()
        assert len(notes) == 5
        assert any("IMP:N" in note for note in notes)
        assert any(note.startswith("F4:N") for note in notes)

    def test_l3_setters(self) -> None:
        deck = _deck("deck_l3.txt")
        deck.set_cell_universe(4, 5, True)
        assert "u=-5" in deck.dumps()
        deck.validate()
        deck.set_cell_lattice(4, 2)
        with pytest.raises(ValueError, match="LAT but no FILL"):
            deck.validate()
        deck.set_cell_fill(4, 0)
        deck.validate()
        deck.set_cell_lattice(4, None)
        deck.set_cell_universe(4, 0, False)
        deck.validate()
        with pytest.raises(ValueError, match="LAT must be 1 or 2"):
            deck.set_cell_lattice(4, 3)
        with pytest.raises(ValueError, match="unknown particle"):
            deck.set_mode(["Q2"])

    def test_negative_universe_flag_and_data_jumps(self) -> None:
        text = "msg\ntitle\n1 0 -1 u=-1\n2 0 -2\n\n1 so 1\n2 so 2\n\nU J 2\n"
        deck = nucleide.mcnp.parse_deck(text)
        universes = {u["number"]: u for u in deck.universes}
        assert set(universes) == {"0", "1", "2"}
        assert universes["1"]["not_truncated"] == "1"
        assert universes["2"]["cells"] == "2"
        deck.validate()

    def test_periodic_surface_link(self) -> None:
        deck = _deck("deck_l3.txt")
        surfs = {s["num"]: s for s in deck.surfs}
        assert surfs["4"]["periodic"] == "3"
        assert surfs["3"]["transform"] == "1"
        assert "4 -3 pz -20" in deck.dumps()
        deck.validate()
        bad = "msg\ntitle\n1 0 -1\n\n1 -9 pz 0\n\n"
        with pytest.raises(ValueError, match="missing periodic surface 9"):
            nucleide.mcnp.parse_deck(bad).validate()
