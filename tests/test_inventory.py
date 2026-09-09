"""Tests for the unit-aware decay inventory facade (synthetic chains only)."""

import math
from pathlib import Path

import pytest

import nucleide

T_HALF_CO60 = 1.663442e8
CHAIN_CO60 = f"""<depletion_chain>
  <nuclide name="Co60" half_life="{T_HALF_CO60!r}"/>
  <nuclide name="Ni60"/>
</depletion_chain>
"""

CHAIN_ABC = f"""<depletion_chain>
  <nuclide name="A" half_life="{math.log(2) / 1e-6!r}">
    <decay type="beta" target="B" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="B" half_life="{math.log(2) / 1e-5!r}">
    <decay type="beta" target="C" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="C"/>
</depletion_chain>
"""


def _chain(tmp_path: Path, text: str, name: str = "chain.xml") -> str:
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_units_round_trip(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_CO60, "co.xml"))
    inv = nucleide.depletion.Inventory(chain, {"Co60": 1.0}, units="Ci")
    assert inv.activities("Bq")["Co60"] == pytest.approx(3.7e10, rel=1e-9)
    grams = inv.masses("g")["Co60"]
    back = nucleide.depletion.Inventory(chain, {"Co60": grams}, units="g")
    assert back.activities("Ci")["Co60"] == pytest.approx(1.0, rel=1e-9)
    assert inv.moles("mol")["Co60"] == pytest.approx(grams / 59.93, rel=1e-3)
    assert inv.mass_fractions() == {"Co60": 1.0}
    assert inv.half_lives_readable()["Co60"].endswith("y")


def test_decay_matches_bateman(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_ABC, "abc.xml"))
    inv = nucleide.depletion.Inventory(chain, {"A": 1e15})
    out = inv.decay(1e5, time_unit="s")
    lam = 1e-6
    assert out.numbers()["A"] == pytest.approx(1e15 * math.exp(-lam * 1e5), rel=1e-8)
    out_y = inv.decay(1e5 / 31_557_600, time_unit="y")
    assert out_y.numbers()["A"] == pytest.approx(out.numbers()["A"], rel=1e-12)


def test_arithmetic_and_csv(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_ABC, "abc.xml"))
    a = nucleide.depletion.Inventory(chain, {"A": 1e15})
    b = nucleide.depletion.Inventory(chain, {"B": 2e14})
    assert (a.add(b)).numbers() == {"A": 1e15, "B": 2e14}
    assert a.mul(2.0).numbers()["A"] == 2e15
    assert a.sub(a).numbers() == {"A": 0.0}
    text = a.to_csv()
    assert nucleide.depletion.Inventory.from_csv(chain, text).numbers() == a.numbers()
    with pytest.raises(ValueError):
        nucleide.depletion.Inventory.from_csv(chain, text + "Nope,1.0\n")


def test_cumulative_and_progeny(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_ABC, "abc.xml"))
    cum = nucleide.depletion.cumulative_decays(chain, {"A": 1e15}, 1e5)
    # Every decayed A atom decayed once (A has no production channel here).
    lam = 1e-6
    assert cum["A"] == pytest.approx(1e15 * (1 - math.exp(-lam * 1e5)), rel=1e-9)
    assert nucleide.depletion.progeny(chain, "A") == [("B", 1.0, "beta")]
    assert nucleide.depletion.branching_fraction(chain, "A", "B") == 1.0
    assert nucleide.depletion.decay_mode(chain, "A", "B") == "beta"
    assert nucleide.depletion.branching_fraction(chain, "A", "C") is None
    edges = nucleide.depletion.chain_edges(chain)
    assert ("A", "B", 1.0, "beta") in edges
    assert ("B", "C", 1.0, "beta") in edges


def test_decay_with_method_matches_cram(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_ABC, "abc.xml"))
    inv = nucleide.depletion.Inventory(chain, {"A": 1e15})
    ref = inv.decay(1e5, time_unit="s")
    for method in ("bateman", "bateman_hp"):
        out = inv.decay(1e5, time_unit="s", method=method)
        for nuc in ("A", "B", "C"):
            assert out.numbers()[nuc] == pytest.approx(ref.numbers()[nuc], rel=1e-8)


def test_error_paths(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_chain(tmp_path, CHAIN_ABC, "abc.xml"))
    # Atom counts store verbatim; chain membership is enforced at solve time.
    stray = nucleide.depletion.Inventory(chain, {"Nope": 1.0})
    with pytest.raises(ValueError, match="unknown nuclide"):
        stray.decay(1.0)
    with pytest.raises(ValueError, match="[Uu]nit"):
        nucleide.depletion.Inventory(chain, {"A": 1.0}, units="furlongs")
    with pytest.raises(ValueError, match="[Uu]nit"):
        nucleide.depletion.Inventory(chain, {"A": 1.0}).decay(1.0, time_unit="eons")
