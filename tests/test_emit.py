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


def test_emit_armi_matches_gnds_cards_and_drift() -> None:
    armi = {"nU235": 5.0, "nU238": 95.0}
    gnds = {"U235": 5.0, "U238": 95.0}
    assert emit.emit_armi_cards(armi, "umetal", density=19.1) == emit.emit_cards(
        gnds, "umetal", density=19.1
    )
    armi_rows = emit.emit_armi_drift_table(armi, "umetal", density=19.1)
    gnds_rows = emit.emit_drift_table(gnds, "umetal", density=19.1)
    assert armi_rows == gnds_rows
    assert armi_rows[0]["code"] == "MCNP"
    assert all(r["rel_drift"] == 0.0 for r in armi_rows)


def test_emit_armi_isomer_mcc3_and_zaid_keys() -> None:
    cards = emit.emit_armi_cards({"nAm242m": 2.0, "U235": 3.0}, "mix", density=10.0)
    assert "95242.80c" in cards["MCNP"]
    assert "92235.80c" in cards["MCNP"]
    assert cards["Serpent"].startswith("mat mix -10\n")
    assert "95242.03c" in cards["Serpent"]
    mcc3 = emit.emit_armi_cards({"U-2355": 5.0, "U235_7": 1.0, "PU2397": 2.0}, "mix", density=10.0)
    assert mcc3 == emit.emit_cards({"U235": 6.0, "Pu239": 2.0}, "mix", density=10.0)
    zaid = emit.emit_armi_cards({"92235": 5.0, "92238": 95.0}, "umetal", density=19.1)
    assert zaid == emit.emit_cards({"U235": 5.0, "U238": 95.0}, "umetal", density=19.1)


def test_emit_armi_rejects_bad_keys() -> None:
    import pytest

    for bad in ["nXx999", "", "DUMP1", "DUMP2", "LFP38"]:
        with pytest.raises(ValueError):
            emit.emit_armi_cards({bad: 1.0}, "mix", density=1.0)


def test_emit_armi_rejects_elemental_keys() -> None:
    import pytest

    for el in ["ZR", "FE", "NA", "O", "C", "W"]:
        with pytest.raises(ValueError, match="expand first"):
            emit.emit_armi_cards({el: 1.0}, "mix", density=1.0)


def test_emit_armi_bare_am242_needs_explicit_state() -> None:
    import pytest

    with pytest.raises(ValueError, match="AM242M"):
        emit.emit_armi_cards({"AM242": 1.0}, "mix", density=1.0)
    ground = emit.emit_armi_cards({"AM242G": 1.0, "U235": 3.0}, "mix", density=1.0)
    assert ground == emit.emit_cards({"Am242": 1.0, "U235": 3.0}, "mix", density=1.0)
    meta = emit.emit_armi_cards({"AM242M": 1.0, "U235": 3.0}, "mix", density=1.0)
    assert meta == emit.emit_cards({"Am242_m1": 1.0, "U235": 3.0}, "mix", density=1.0)
