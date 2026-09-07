"""Tests for ARMI dialect labels and composition checks (no ARMI needed)."""

import pytest

import nucleide


class TestArmiDialects:
    def test_db_and_zaid_forms(self) -> None:
        assert nucleide.nuclei.armi_to_nucid("nU235").nucid == 922_350_000
        assert nucleide.nuclei.armi_to_nucid("U235").nucid == 922_350_000
        assert nucleide.nuclei.armi_to_nucid("92235").nucid == 922_350_000
        assert nucleide.nuclei.armi_to_nucid("92238").nucid == 922_380_000

    def test_round_trip_labels(self) -> None:
        u235 = nucleide.nuclei.Nuclide("U235")
        assert nucleide.nuclei.nucid_to_armi(u235) == "nU235"
        assert nucleide.nuclei.mcc3_to_nucid("U-2355").nucid == 922_350_000

    def test_isomers_and_errors(self) -> None:
        assert nucleide.nuclei.armi_to_nucid("nAm242m").nucid == 952_420_001
        with pytest.raises(ValueError):
            nucleide.nuclei.armi_to_nucid("nXx999")
        with pytest.raises(ValueError):
            nucleide.nuclei.armi_to_nucid("")
        with pytest.raises(ValueError):
            nucleide.nuclei.mcc3_to_nucid("U235_8")


class TestCompositionChecks:
    def test_clean_composition(self) -> None:
        assert nucleide.material.check_labels({"U235": 0.04, "U238": 0.96}) == []
        assert nucleide.material.audit_material({"U235": 0.04, "U238": 0.96}) == []

    def test_collisions_reported(self) -> None:
        hits = nucleide.material.check_labels({"Am242_m1": 0.5, "Am242_m2": 0.5})
        assert any(h["width"] == 6 and set(h["members"]) == {"Am242_m1", "Am242_m2"} for h in hits)

    def test_audit_finds_problems(self) -> None:
        # Absolute masses normalize, so fractions always sum; the sum check
        # fires only on non-positive totals.
        assert nucleide.material.audit_material({})[0]["kind"] == "FractionsDontSum"
        issues = nucleide.material.audit_material({"U235": -1.0})
        assert any(i["kind"] == "NegativeMass" for i in issues)
