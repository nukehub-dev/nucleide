"""Python-side tests for the CRAM depletion core."""

import math
from pathlib import Path

import pytest

import nucleide

CHAIN = Path(__file__).parent.parent / "fixtures" / "depletion" / "chain_simple.xml"


def test_chain_parse() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    assert len(chain.nuclides) == 9
    assert "I135" in chain.nuclides
    assert chain.index_of("Xe135") is not None
    assert chain.index_of("Nope135") is None


def test_cram_decay_matches_analytic() -> None:
    # I135 -> Xe135 with branching ratio 1; half lives from the chain.
    # Build a minimal two-isotope check against the analytic Bateman result.
    chain = nucleide.depletion.read_chain(str(CHAIN))
    t_half_i = 2.36520e4
    lam = math.log(2) / t_half_i
    dt = 1e5
    n0 = {"I135": 1e16}
    out = nucleide.depletion.deplete(chain, n0, dt, order=48)
    expected = 1e16 * math.exp(-lam * dt)
    assert out["I135"] == pytest.approx(expected, rel=1e-8)
    # Xe135 gains from decay but also decays itself (half life 32904 s)
    assert out["Xe135"] > 0


def test_order16_and_48_agree() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    n0 = {"I135": 1e15, "U235": 1e20}
    a = nucleide.depletion.deplete(chain, n0, 5e4, order=48)
    b = nucleide.depletion.deplete(chain, n0, 5e4, order=16)
    for k in a:
        assert a[k] == pytest.approx(b[k], rel=1e-4)


def test_bad_order_rejected() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    with pytest.raises(ValueError):
        nucleide.depletion.deplete(chain, {"I135": 1.0}, 1.0, order=32)


def test_unknown_rate_nuclide_rejected() -> None:
    chain = nucleide.depletion.read_chain(str(CHAIN))
    with pytest.raises(ValueError):
        nucleide.depletion.deplete(chain, {"I135": 1.0}, 1.0, rates={"Zzz9:(n,gamma)": 1e-3})
