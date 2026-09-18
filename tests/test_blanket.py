"""Python-side tests for the TBR / blanket bookkeeping core.

Stellaris Point-A anchors (raw TBR 1.1070, 1.074 after the 3% ECRH-port
haircut, energy multiplication 1.20, burn 416.6 g/day at Point A) gate the
arithmetic. The burn constant is transcribed again, independently, below —
the same two-spelling discipline as the damage Ballabio checks. All inputs
are synthetic hand-picked values; nothing reads files and nothing solves
transport.
"""

import math

import pytest

import nucleide.blanket as blk

RAW_TBR = 1.1070
PENALIZED_TBR = 1.07379  # 1.1070 * 0.97
POINT_A_FUSION_MW = 2714.82  # power implied by 416.6 g/day through (F1)


def _burn_constant_py() -> float:
    """Independent transcription of the (F1) burn constant.

    One triton (3.0160492 u) per 17.6 MeV of D-T fusion energy, scaled to
    g/day per MW: ``seconds/day * W/MW * kg-per-reaction * g/kg /
    J-per-reaction``. Uses ``math.pow`` for the powers of ten rather than
    ``**`` — this repo's mypy types float**float as Any — and groups the
    factors differently from the Rust const expression.
    """
    seconds_per_day = 24.0 * 60.0 * 60.0
    watts_per_mw = math.pow(10.0, 6.0)
    mev_per_reaction = 17.6
    joules_per_mev = math.pow(10.0, 6.0) * 1.602_176_634e-19
    grams_per_reaction = 3.0160492 * 1.660_539_066_60e-27 * math.pow(10.0, 3.0)
    reactions_per_mw_s = watts_per_mw / (mev_per_reaction * joules_per_mev)
    return seconds_per_day * reactions_per_mw_s * grams_per_reaction


def test_burn_constant_matches_independent_transcription() -> None:
    assert blk.tritium_burn_g_per_day(1.0) == pytest.approx(_burn_constant_py(), rel=1e-12)
    assert blk.tritium_burn_g_per_day(1000.0) == pytest.approx(153.45399532, rel=1e-9)


def test_stellaris_raw_tbr_gate() -> None:
    assert blk.tbr_from_tallies(RAW_TBR, 1.0) == pytest.approx(1.1070, rel=0, abs=1e-12)
    assert blk.tbr_from_tallies(2.0 * RAW_TBR, 2.0) == pytest.approx(1.1070, rel=0, abs=1e-12)


def test_stellaris_port_penalty_gate() -> None:
    got = blk.apply_port_penalty(RAW_TBR, [0.03])
    assert got == pytest.approx(PENALIZED_TBR, rel=0, abs=1e-12)
    assert got == pytest.approx(1.074, rel=0, abs=5e-4)
    assert blk.apply_port_penalty(RAW_TBR, []) == pytest.approx(RAW_TBR, rel=0, abs=0.0)
    assert blk.apply_port_penalty(1.0, [0.03, 0.02]) == pytest.approx(0.97 * 0.98, rel=1e-12)


def test_stellaris_power_multiplication_gate() -> None:
    assert blk.energy_multiplication(3600.0, 3000.0) == pytest.approx(1.20, rel=1e-12)
    assert blk.blanket_power(3000.0, 1.20) == pytest.approx(3600.0, rel=1e-12)


def test_stellaris_burn_gate() -> None:
    assert blk.tritium_burn_g_per_day(POINT_A_FUSION_MW) == pytest.approx(416.6, rel=0, abs=1e-2)


def test_stellaris_margin_gates() -> None:
    assert blk.breeding_margin(PENALIZED_TBR) == pytest.approx(0.07379, rel=0, abs=1e-12)
    assert blk.meets_requirement(PENALIZED_TBR, 1.0) is True
    assert blk.meets_requirement(PENALIZED_TBR, 1.05) is True
    assert blk.meets_requirement(PENALIZED_TBR, 1.10) is False
    net = blk.net_surplus_g_per_day(PENALIZED_TBR, POINT_A_FUSION_MW)
    assert net == pytest.approx((PENALIZED_TBR - 1.0) * 416.6, rel=1e-4)
    assert net == pytest.approx(30.74, rel=0, abs=5e-2)
    # Sub-breakeven is a reported deficit, never clamped.
    assert blk.net_surplus_g_per_day(0.9, 1000.0) < 0.0


def test_malformed_inputs_are_loud() -> None:
    with pytest.raises(ValueError, match="source_neutrons"):
        blk.tbr_from_tallies(1.0, 0.0)
    with pytest.raises(ValueError, match="tritons_bred"):
        blk.tbr_from_tallies(-1.0, 1.0)
    with pytest.raises(ValueError, match="port fraction"):
        blk.apply_port_penalty(1.0, [1.0])
    with pytest.raises(ValueError, match="port fraction"):
        blk.apply_port_penalty(1.0, [-0.01])
    with pytest.raises(ValueError, match="effective_tbr"):
        blk.breeding_margin(-0.5)
    with pytest.raises(ValueError, match="required_tbr"):
        blk.meets_requirement(1.1, 0.0)
    with pytest.raises(ValueError, match="fusion_power_mw"):
        blk.energy_multiplication(1.0, -3.0)
    with pytest.raises(ValueError, match="fusion_power_mw"):
        blk.blanket_power(0.0, 1.2)
    with pytest.raises(ValueError, match="multiplication"):
        blk.blanket_power(100.0, -1.2)
    with pytest.raises(ValueError, match="fusion_power_mw"):
        blk.tritium_burn_g_per_day(0.0)
    with pytest.raises(ValueError, match="effective_tbr"):
        blk.net_surplus_g_per_day(-1.0, 100.0)
