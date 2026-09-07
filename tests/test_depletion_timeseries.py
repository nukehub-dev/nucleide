"""Python-side tests for depletion time-series behavior.

Exercises multi-step integration through the public `deplete` API against
the analytic Bateman solution on the synthetic A -> B -> C fixture chain,
plus the activity (A = lambda*N) and decay-heat (H = lambda*N*E) formulas
documented in `nucleide-depletion`'s `integrate` module.
"""

import math
from pathlib import Path

import pytest

import nucleide

CHAIN = Path(__file__).parent.parent / "fixtures" / "depletion" / "chain_abc.xml"

LAM_A = 1e-6
LAM_B = 1e-5
MEV_TO_J = 1.602e-13
E_A_MEV = 1.0
E_B_MEV = 2.0


def bateman(n0: float, t: float) -> dict[str, float]:
    na = n0 * math.exp(-LAM_A * t)
    nb = n0 * LAM_A / (LAM_B - LAM_A) * (math.exp(-LAM_A * t) - math.exp(-LAM_B * t))
    return {"A": na, "B": nb, "C": n0 - na - nb}


def run_series(n0: dict[str, float], dt: float, count: int) -> list[dict[str, float]]:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    state = {nuc: n0.get(nuc, 0.0) for nuc in chain.nuclides}
    rows = [dict(state)]
    for _ in range(count):
        state = nucleide.depletion.deplete(chain, state, dt, order=48)
        rows.append(dict(state))
    return rows


def test_multistep_matches_bateman() -> None:
    n0 = {"A": 1e15}
    dt, count = 2e4, 5
    rows = run_series(n0, dt, count)
    assert len(rows) == count + 1
    for k, row in enumerate(rows):
        want = bateman(n0["A"], k * dt)
        for nuc in ("A", "B", "C"):
            assert row[nuc] == pytest.approx(want[nuc], rel=1e-6), f"node {k} {nuc}"


def test_atoms_conserved() -> None:
    rows = run_series({"A": 1e15}, 2e4, 5)
    for k, row in enumerate(rows):
        total = row["A"] + row["B"] + row["C"]
        assert total == pytest.approx(1e15, rel=1e-8), f"node {k}"


def test_activity_and_decay_heat_formulas() -> None:
    rows = run_series({"A": 1e15}, 2e4, 5)
    for k, row in enumerate(rows):
        t = k * 2e4
        want = bateman(1e15, t)
        assert LAM_A * row["A"] == pytest.approx(LAM_A * want["A"], rel=1e-6)
        assert LAM_B * row["B"] == pytest.approx(LAM_B * want["B"], rel=1e-6)
        heat_a = LAM_A * row["A"] * E_A_MEV * MEV_TO_J
        heat_b = LAM_B * row["B"] * E_B_MEV * MEV_TO_J
        assert heat_a == pytest.approx(LAM_A * want["A"] * E_A_MEV * MEV_TO_J, rel=1e-6)
        assert heat_b == pytest.approx(LAM_B * want["B"] * E_B_MEV * MEV_TO_J, rel=1e-6)
        # Stable C carries no activity or heat.
        assert row["C"] == pytest.approx(want["C"], rel=1e-6)


def test_invalid_dt_rejected() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    with pytest.raises(ValueError):
        nucleide.depletion.deplete(chain, {"A": 1.0}, 0.0, order=48)
    with pytest.raises(ValueError):
        nucleide.depletion.deplete(chain, {"A": 1.0}, -1.0, order=48)
