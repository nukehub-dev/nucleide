"""Golden tests for data accessors, input parsing, enrichment, materials."""

import math
from pathlib import Path
from typing import Any

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures"


class TestDataAccessors:
    def test_half_life_anchors(self) -> None:
        i135 = nucleide.half_life("I135")
        u238 = nucleide.half_life("U238")
        assert i135 == pytest.approx(2.36520e4, rel=1e-4)
        assert u238 == pytest.approx(1.41e17, rel=1e-2)  # 4.47 Gyr in seconds
        assert nucleide.decay_constant("I135") is not None
        lam = nucleide.decay_constant("I135")
        hl = nucleide.half_life("I135")
        assert lam is not None and hl is not None
        assert abs(lam * hl - math.log(2)) < 1e-12

    def test_q_values(self) -> None:
        hcap = nucleide.q_value_capture("H1")
        assert hcap == pytest.approx(2.2246, abs=1e-3)
        u5 = nucleide.q_value_capture("U235")
        assert u5 == pytest.approx(6.5455, abs=1e-3)
        qa = nucleide.q_value_alpha("U238")
        assert qa == pytest.approx(4.2699, abs=1e-3)


class TestInpParsing:
    def test_materials(self) -> None:
        mats = nucleide.read_inp(str(FIX / "mcnp" / "inp" / "mcnp_inp.txt"))
        assert len(mats) >= 1
        m = mats[0]
        assert m["number"] >= 1
        assert len(m["fractions"]) > 0

    def test_comments_fixture(self) -> None:
        mats = nucleide.read_inp(str(FIX / "mcnp" / "inp" / "mcnp_inp_comments.txt"))
        # commented-out M cards must not produce entries
        joined = " ".join(str(m["number"]) for m in mats)
        assert joined.count("4") <= 1


class TestEnrichment:
    def test_default_cascade_solves(self) -> None:
        c = nucleide.Cascade.default_uranium()
        c.solve()
        # Natural feed 0.72% -> product must exceed tails assay
        assert 0.0072 < c.x_prod_j < 1.0
        assert 0.0 < c.x_tail_j < 0.0072
        assert c.swu_per_feed > 0
        swu_pp = c.separative_work_per_product()
        assert swu_pp != 0.0

    def test_deterministic(self) -> None:
        a = nucleide.Cascade.default_uranium()
        a.solve()
        b = nucleide.Cascade.default_uranium()
        b.solve()
        assert a.x_prod_j == b.x_prod_j


class TestMaterials:
    def test_formula_water(self) -> None:
        comp = nucleide.from_formula("H2O")
        h_total = sum(v for k, v in comp.items() if k.startswith("H"))
        o_total = sum(v for k, v in comp.items() if k.startswith("O"))
        assert h_total == pytest.approx(2.0 / 3.0, abs=1e-9)
        assert o_total == pytest.approx(1.0 / 3.0, abs=1e-9)
        # natural hydrogen carries deuterium
        assert "H2" in comp and comp["H2"] > 0

    def test_formula_nested_parens(self) -> None:
        comp = nucleide.from_formula("Ca(OH)2")
        assert sum(comp.values()) == pytest.approx(1.0, abs=1e-9)
        assert any(k.startswith("Ca") for k in comp)

    def test_formula_errors(self) -> None:
        with pytest.raises(ValueError):
            nucleide.from_formula("Xx2O")
        with pytest.raises(ValueError):
            nucleide.from_formula("(H2O")

    def test_activity_cs137(self) -> None:
        grams = 1e-3  # 1 mg of Cs137
        out = nucleide.activity({"Cs137": grams})
        total = sum(v for k, v in out.items() if k != "specific")
        specific = out["specific"]
        # Cs137 specific activity ~3.2 TBq/g
        assert 1e12 < specific < 1e13
        assert total == pytest.approx(specific * grams, rel=1e-9)


class TestMaterialsCompendium:
    COMPENDIUM = FIX / "data" / "MaterialsCompendium.json"

    def test_load_and_lookup(self) -> None:
        lib = nucleide.MaterialsCompendium.load(str(self.COMPENDIUM))
        assert len(lib) == 411
        air = lib.get("air (DRY, near sea level)")
        assert air is not None
        assert air["mat_num"] == 4
        fractions: dict[Any, Any] = air["fractions"]  # keyed by integer ZAID
        assert fractions[7014] == pytest.approx(0.752316, abs=1e-5)

    def test_as_material_named_fractions(self) -> None:
        lib = nucleide.MaterialsCompendium.load(str(self.COMPENDIUM))
        bone = lib.get("Bone Equivalent Plastic, B-110", as_material=True)
        assert bone is not None
        fractions: dict[str, Any] = bone["fractions"]
        assert "H1" in fractions

    def test_missing_returns_none(self) -> None:
        lib = nucleide.MaterialsCompendium.load(str(self.COMPENDIUM))
        assert lib.get("Unobtainium-999") is None
