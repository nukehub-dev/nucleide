"""Python-side tests for the fuel-cycle micro-adds (separate/blend/CUSUM/SWU)."""

import math

import pytest

import nucleide.enrichment as enr
import nucleide.material as mat

FEED = {
    "U235": 10.0,
    "U238": 90.0,
    "Pu239": 1.0,
    "Pu240": 2.0,
    "Am241": 3.0,
    "Am242": 2.8,
}
EFFS = {"U235": 0.7, "U238": 0.7, "Pu239": 0.4, "Pu240": 0.4, "Am241": 0.4}


def test_separate_splits_and_conserves_mass() -> None:
    product, tails = mat.separate_material(FEED, EFFS)
    assert product == pytest.approx(
        {"U235": 7.0, "U238": 63.0, "Pu239": 0.4, "Pu240": 0.8, "Am241": 1.2}
    )
    assert tails == pytest.approx(
        {"U235": 3.0, "U238": 27.0, "Pu239": 0.6, "Pu240": 1.2, "Am241": 1.8, "Am242": 2.8}
    )
    assert sum(product.values()) == pytest.approx(72.4)
    assert sum(tails.values()) == pytest.approx(36.4)
    for nuc, mass in FEED.items():
        assert product.get(nuc, 0.0) + tails.get(nuc, 0.0) == pytest.approx(mass)


def test_separate_rejects_bad_efficiencies() -> None:
    with pytest.raises(ValueError):
        mat.separate_material(FEED, {"U235": 1.5})
    with pytest.raises(ValueError):
        mat.separate_material(FEED, {"U235": float("nan")})


def test_blend_normalizes_ratios() -> None:
    out = mat.blend_material([({"U235": 1.0, "Pu239": 1.0}, 1.0), ({"U238": 1.0}, 2.0)])
    assert out == pytest.approx({"U235": 1.0 / 3.0, "Pu239": 1.0 / 3.0, "U238": 2.0 / 3.0})
    assert sum(out.values()) == pytest.approx(4.0 / 3.0)


def test_blend_rejects_degenerate_recipes() -> None:
    with pytest.raises(ValueError):
        mat.blend_material([])
    with pytest.raises(ValueError):
        mat.blend_material([({"U235": 1.0}, 0.0)])
    with pytest.raises(ValueError):
        mat.blend_material([({"U235": 1.0}, -1.0)])


def test_cusum_step_sequence() -> None:
    cusum = mat.Cusum()
    for _ in range(10):
        assert cusum.update(1.0) is False
    assert cusum.statistic() == pytest.approx(0.0)
    assert cusum.update(2.0) is False
    assert cusum.statistic() == pytest.approx(0.7583352368020273)
    assert cusum.update(2.0) is False
    assert cusum.update(2.0) is True
    assert cusum.status() is True
    assert cusum.count() == 13
    cusum.reset()
    assert cusum.status() is False
    assert cusum.statistic() == pytest.approx(0.0)
    assert cusum.count() == 0


def test_cusum_rejects_bad_tuning() -> None:
    with pytest.raises(ValueError):
        mat.Cusum(ref_shift_k=-0.1)
    with pytest.raises(ValueError):
        mat.Cusum(alarm_h=0.0)


def test_swu_closed_form() -> None:
    def value(x: float) -> float:
        return (2.0 * x - 1.0) * math.log(x / (1.0 - x))

    xf, xp, xt = 0.0072, 0.05, 0.002
    assert enr.value_func(xf) == pytest.approx(value(xf), rel=1e-12)
    assert enr.value_func(xp) == pytest.approx(value(xp), rel=1e-12)
    assert enr.value_func(xt) == pytest.approx(value(xt), rel=1e-12)
    feed, tails = 92.3076923076923, 82.3076923076923
    expected = 10.0 * value(xp) + tails * value(xt) - feed * value(xf)
    assert expected == pytest.approx(87.59916188589864, rel=1e-12)
    assert 10.0 * enr.swu_per_prod(xf, xp, xt) == pytest.approx(expected, rel=1e-12)
    assert feed * enr.swu_per_feed(xf, xp, xt) == pytest.approx(expected, rel=1e-12)
    assert tails * enr.swu_per_tail(xf, xp, xt) == pytest.approx(expected, rel=1e-12)
