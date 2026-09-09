"""Python-side tests for the Bateman/HP decay fast path (`method=`).

Mirrors the Rust oracles in `crates/depletion/src/bateman.rs` through the
public API: `deplete` / `deplete_series` / `DepletionSystem.solve` /
`Inventory.decay` with `method="bateman"` / `"bateman_hp"`.
"""

import math
from pathlib import Path

import pytest

import nucleide

CHAIN = Path(__file__).parent.parent / "fixtures" / "depletion" / "chain_simple.xml"
CHAIN_ABC = Path(__file__).parent.parent / "fixtures" / "depletion" / "chain_abc.xml"

LAM_A = 1e-6
LAM_B = 1e-5


def bateman_abc(n0: float, t: float) -> dict[str, float]:
    na = n0 * math.exp(-LAM_A * t)
    nb = n0 * LAM_A / (LAM_B - LAM_A) * (math.exp(-LAM_A * t) - math.exp(-LAM_B * t))
    return {"A": na, "B": nb, "C": n0 - na - nb}


def _chain_from_xml(tmp_path: Path, text: str, name: str):  # type: ignore[no-untyped-def]
    path = tmp_path / name
    path.write_text(text)
    return nucleide.depletion.read_chain(str(path))


def test_bateman_matches_analytic_i135() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    t_half_i = 2.36520e4
    lam = math.log(2) / t_half_i
    dt = 1e5
    for method in ("bateman", "bateman_hp"):
        out = nucleide.depletion.deplete(chain, {"I135": 1e16}, dt, method=method)
        assert out["I135"] == pytest.approx(1e16 * math.exp(-lam * dt), rel=1e-8)
        assert out["Xe135"] > 0


def test_bateman_matches_analytic_abc() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN_ABC))
    n0 = {"A": 1e15}
    dt = 1e5
    for method in ("bateman", "bateman_hp"):
        out = nucleide.depletion.deplete(chain, n0, dt, method=method)
        want = bateman_abc(1e15, dt)
        assert out["A"] == pytest.approx(want["A"], rel=1e-8)
        assert out["B"] == pytest.approx(want["B"], rel=1e-7)
        assert out["C"] == pytest.approx(want["C"], rel=1e-9)
        total = out["A"] + out["B"] + out["C"]
        assert total == pytest.approx(1e15, rel=1e-8)


def test_bateman_vs_cram48_series_band() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN_ABC))
    n0 = {"A": 1e15}
    dts = [2e4] * 5
    ref_rows = nucleide.depletion.deplete_series(chain, n0, dts, order=48)["atoms"]
    for method in ("bateman", "bateman_hp"):
        out = nucleide.depletion.deplete_series(chain, n0, dts, method=method)
        assert out["times"] == pytest.approx([2e4, 4e4, 6e4, 8e4, 1e5])
        for k, row in enumerate(out["atoms"]):
            want = bateman_abc(1e15, (k + 1) * 2e4)
            for nuc in ("A", "B", "C"):
                assert row[nuc] == pytest.approx(want[nuc], rel=1e-6), f"{method} node {k}"
                assert row[nuc] == pytest.approx(ref_rows[k][nuc], rel=1e-6)
            total = row["A"] + row["B"] + row["C"]
            assert total == pytest.approx(1e15, rel=1e-8)


def test_repeated_calls_reuse_cache_bit_identical() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN_ABC))
    n0 = {"A": 1e14, "B": 5e13, "C": 1e10}
    for method in ("bateman", "bateman_hp"):
        a = nucleide.depletion.deplete(chain, n0, 2e4, method=method)
        b = nucleide.depletion.deplete(chain, n0, 2e4, method=method)
        assert a == b


def test_wide_spread_chain_hp_beats_or_matches(tmp_path: Path) -> None:
    # Es-254-style 10-order spread: ~1 y parent, 1 us daughter (synthetic).
    lam_p = math.log(2) / 3.15576e7
    lam_d = math.log(2) / 1.0e-6
    text = f"""<depletion_chain>
  <nuclide name="P" half_life="{math.log(2) / lam_p!r}">
    <decay type="beta" target="D" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="D" half_life="{math.log(2) / lam_d!r}">
    <decay type="beta" target="S" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="S"/>
</depletion_chain>
"""
    chain = _chain_from_xml(tmp_path, text, "wide.xml")
    n0 = {"P": 1e15}
    dt = 1e6
    exp_p = math.exp(-lam_p * dt)
    exp_d = math.exp(-lam_d * dt)
    want_p = 1e15 * exp_p
    want_d = 1e15 * lam_p / (lam_d - lam_p) * (exp_p - exp_d)
    want_s = 1e15 - want_p - want_d
    std = nucleide.depletion.deplete(chain, n0, dt, method="bateman")
    hp = nucleide.depletion.deplete(chain, n0, dt, method="bateman_hp")

    def err(row: dict[str, float]) -> float:
        return max(
            abs(row[n] - w) / max(abs(w), 1e-30)
            for n, w in (("P", want_p), ("D", want_d), ("S", want_s))
        )

    assert err(hp) <= err(std)
    for row in (std, hp):
        assert row["P"] == pytest.approx(want_p, rel=1e-8)
        assert row["D"] == pytest.approx(want_d, rel=1e-8)
        assert row["S"] == pytest.approx(want_s, rel=1e-8)
        assert sum(row.values()) == pytest.approx(1e15, rel=1e-8)


def test_zero_dt_returns_input() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN_ABC))
    n0 = {"A": 1e14, "B": 5e13, "C": 1e10}
    for method in ("bateman", "bateman_hp"):
        assert nucleide.depletion.deplete(chain, n0, 0.0, method=method) == n0


def test_d1_near_equal_falls_back_to_cram(tmp_path: Path) -> None:
    lam = 1e-6
    text = f"""<depletion_chain>
  <nuclide name="A" half_life="{math.log(2) / lam!r}">
    <decay type="beta" target="B" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="B" half_life="{math.log(2) / (lam * (1.0 + 1e-13))!r}"/>
</depletion_chain>
"""
    chain = _chain_from_xml(tmp_path, text, "d1.xml")
    n0 = {"A": 1e12}
    for method in ("bateman", "bateman_hp"):
        got = nucleide.depletion.deplete(chain, n0, 1e5, method=method)
        want = nucleide.depletion.deplete(chain, n0, 1e5, order=48)
        assert got == want


def test_d2_cyclic_falls_back_to_cram(tmp_path: Path) -> None:
    text = f"""<depletion_chain>
  <nuclide name="A" half_life="{math.log(2) / LAM_A!r}">
    <decay type="beta" target="B" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="B" half_life="{math.log(2) / LAM_B!r}">
    <decay type="beta" target="A" branching_ratio="1.0"/>
  </nuclide>
</depletion_chain>
"""
    chain = _chain_from_xml(tmp_path, text, "d2.xml")
    n0 = {"A": 1e12}
    for method in ("bateman", "bateman_hp"):
        got = nucleide.depletion.deplete(chain, n0, 1e5, method=method)
        want = nucleide.depletion.deplete(chain, n0, 1e5, order=48)
        assert got == want


def test_d3_stable_limit_forms(tmp_path: Path) -> None:
    text = """<depletion_chain>
  <nuclide name="S1"/>
  <nuclide name="S2"/>
</depletion_chain>
"""
    chain = _chain_from_xml(tmp_path, text, "d3.xml")
    n0 = {"S1": 3.0e11, "S2": 2.0e11}
    for method in ("bateman", "bateman_hp"):
        assert nucleide.depletion.deplete(chain, n0, 1e5, method=method) == n0


def test_d4_reactions_fall_back_to_cram(tmp_path: Path) -> None:
    text = f"""<depletion_chain>
  <nuclide name="A" half_life="{math.log(2) / LAM_A!r}">
    <decay type="beta" target="B" branching_ratio="1.0"/>
    <reaction type="(n,gamma)" target="B" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="B" half_life="{math.log(2) / LAM_B!r}">
    <decay type="beta" target="C" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="C"/>
</depletion_chain>
"""
    chain = _chain_from_xml(tmp_path, text, "d4.xml")
    n0 = {"A": 1e12}
    rates = {"A:(n,gamma)": 1e-7}
    for method in ("bateman", "bateman_hp"):
        got = nucleide.depletion.deplete(chain, n0, 1e5, rates=rates, method=method)
        want = nucleide.depletion.deplete(chain, n0, 1e5, rates=rates, order=48)
        assert got == want


def test_bad_method_rejected() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    with pytest.raises(ValueError, match="unsupported method"):
        nucleide.depletion.deplete(chain, {"I135": 1.0}, 1.0, method="cram32")


def test_order_still_selects_cram_when_method_default() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    n0 = {"I135": 1e15, "U235": 1e20}
    a = nucleide.depletion.deplete(chain, n0, 5e4, order=16)
    b = nucleide.depletion.deplete(chain, n0, 5e4, order=16, method="cram48")
    assert a == b


def test_solve_and_inventory_with_method(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN_ABC))
    sys = nucleide.depletion.build_depletion_system(chain, {})
    want = bateman_abc(1e15, 1e5)
    got = sys.solve({"A": 1e15}, 1e5, method="bateman_hp")
    assert got["B"] == pytest.approx(want["B"], rel=1e-7)
    vec = sys.solve_vec([1e15, 0.0, 0.0], 1e5, method="bateman")
    assert vec[1] == pytest.approx(want["B"], rel=1e-7)
    inv = nucleide.depletion.Inventory(chain, {"A": 1e15})
    out = inv.decay(1e5, method="bateman")
    assert out.numbers()["A"] == pytest.approx(want["A"], rel=1e-8)
