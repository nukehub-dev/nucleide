"""Python-side tests for single-material code emission.

Exercises `nucleide.emit` on a synthetic uranium-metal composition: all five
cards are produced, the drift table is lossless, and FLUKA's missing O16
isotope entry shows up as reported drift on a UO2-style mix.
"""

from nucleide import emit


def test_emit_cards_five_codes() -> None:
    cards = emit.emit_cards({"U235": 5.0, "U238": 95.0}, "umetal", density=19.1)
    assert set(cards) == {"MCNP", "Serpent", "FLUKA", "ALARA", "PARTISN"}
    assert cards["MCNP"].startswith("m1 ")
    assert "92235.80c" in cards["MCNP"]
    assert cards["Serpent"].startswith("mat umetal -19.1\n")
    assert "92235.03c -0.05" in cards["Serpent"]
    assert "COMPOUND" in cards["FLUKA"]
    assert "235-U" in cards["FLUKA"]
    assert cards["ALARA"].startswith("mixture umetal\n")
    assert cards["PARTISN"].startswith("TITLE umetal")


def test_emit_drift_table_lossless() -> None:
    rows = emit.emit_drift_table({"U235": 5.0, "U238": 95.0}, "umetal", density=19.1)
    assert [r["code"] for r in rows] == ["MCNP", "Serpent", "FLUKA", "ALARA", "PARTISN"]
    for row in rows:
        assert row["mass_in"] == 100.0
        assert row["mass_out"] == 100.0
        assert row["rel_drift"] == 0.0
        assert row["dropped"] == []
    reparsed = {r["code"] for r in rows if r["reparsed"]}
    assert reparsed == {"MCNP", "ALARA"}


def test_emit_drift_table_reports_fluka_drop() -> None:
    rows = emit.emit_drift_table({"U235": 1.0, "O16": 2.0}, "uox", density=10.0)
    fluka = next(r for r in rows if r["code"] == "FLUKA")
    assert fluka["mass_out"] == 1.0
    assert fluka["rel_drift"] == 2.0 / 3.0
    assert [d["nuclide"] for d in fluka["dropped"]] == ["O16"]
    assert all("no-fluka-name" in d["reason"] for d in fluka["dropped"])


def test_emit_needs_density_for_serpent() -> None:
    import pytest

    with pytest.raises(ValueError, match="density"):
        emit.emit_cards({"U235": 1.0}, "umetal")
