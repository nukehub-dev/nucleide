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
