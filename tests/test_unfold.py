"""Python-side tests for the unfolding core (round-trip gates + input errors)."""

import math

import pytest

import nucleide.unfold as unfold

N_GROUPS = 6


def _synthetic_response(n_det: int, n_groups: int, width: float) -> list[list[float]]:
    """Deterministic synthetic response matrix (hand-built, never evaluated data)."""
    rows = []
    for i in range(n_det):
        center = i * (n_groups - 1) / max(n_det - 1, 1)
        rows.append([math.exp(-(((j - center) / width) ** 2)) + 1e-3 for j in range(n_groups)])
    return rows


def _log_midpoints(n_groups: int, lo: float, hi: float) -> list[float]:
    step = math.log(hi / lo) / n_groups
    return [math.exp(math.log(lo) + (j + 0.5) * step) for j in range(n_groups)]


def _maxwellian_plus_inv_e(midpoints: list[float], kt: float) -> list[float]:
    return [e * math.exp(-e / kt) + 1e-6 / e for e in midpoints]


def test_forward_fold_matches_hand_matvec() -> None:
    response = _synthetic_response(3, 4, 1.5)
    spectrum = [1.0, 2.0, 3.0, 4.0]
    rates = unfold.forward_fold(response, spectrum)
    for i, row in enumerate(response):
        assert rates[i] == pytest.approx(sum(r * p for r, p in zip(row, spectrum, strict=True)))


def test_exact_guess_is_a_fixed_point() -> None:
    response = _synthetic_response(4, 12, 2.0)
    midpoints = _log_midpoints(12, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    out = unfold.sandii(response, rates, truth)
    assert out["iterations"] == 1
    assert out["max_rel_change"] == 0.0
    assert out["spectrum"] == truth
    assert out["rates"] == pytest.approx(rates, rel=1e-12)
    assert out["rate_factors"] == pytest.approx([1.0] * len(rates), rel=1e-12)


def test_recovers_determined_system_from_biased_guess() -> None:
    response = _synthetic_response(N_GROUPS, N_GROUPS, 0.9)
    midpoints = _log_midpoints(N_GROUPS, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [3.0 * p for p in truth]
    out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=10_000)
    assert out["spectrum"] == pytest.approx(truth, rel=1e-6)
    assert out["rate_factors"] == pytest.approx([1.0] * N_GROUPS, rel=1e-9)
    assert out["iterations"] >= 1


def test_underdetermined_round_trip_reproduces_rates() -> None:
    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    # Non-uniform bias: part lies in the response null space and must keep
    # the guess; the recoverable part must still converge.
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]

    def max_rel_err(a: list[float], b: list[float]) -> float:
        return max(abs(x - y) / y for x, y in zip(a, b, strict=True))

    out = unfold.sandii(response, rates, guess, tolerance=1e-9, max_iterations=50_000)
    assert out["rate_factors"] == pytest.approx([1.0] * 6, rel=1e-6)
    assert max_rel_err(out["spectrum"], truth) < max_rel_err(guess, truth) / 2.0
    assert all(p > 0.0 for p in out["spectrum"])


def test_irdff_ii_analytical_benchmark_shapes_recover() -> None:
    # IRDFF-II (Trkov et al., Nucl. Data Sheets 163 (2020) 1, arXiv:1909.03336)
    # names analytical benchmark-field shapes; three are reused here as
    # published facts with citation. The library's tabulated group spectra
    # are IAEA-copyright data files and are deliberately NOT used.
    n_groups = 6
    midpoints = _log_midpoints(n_groups, 1e-9, 20.0)
    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K, eV -> MeV
    shapes = {
        "thermal": [e * math.exp(-e / thermal_kt) + 1e-30 for e in midpoints],
        "1/E": [1.0 / e + 1e-30 for e in midpoints],
        "25keV": [math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in midpoints],
    }
    for shape in shapes.values():
        response = _synthetic_response(n_groups, n_groups, 0.8)
        rates = unfold.forward_fold(response, shape)
        guess = [5.0 * p for p in shape]
        out = unfold.sandii(response, rates, guess, tolerance=1e-11, max_iterations=100_000)
        assert out["spectrum"] == pytest.approx(shape, rel=1e-5)


def test_non_convergence_is_a_hard_fail() -> None:
    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.0 * p for p in truth]
    with pytest.raises(ValueError, match="did not converge"):
        unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=1)


def test_validation_errors_name_the_cause() -> None:
    response = _synthetic_response(2, 3, 1.0)
    rates = [1.0, 2.0]
    guess = [1.0, 1.0, 1.0]
    with pytest.raises(ValueError, match="shape mismatch"):
        unfold.sandii(response, [1.0], guess)
    with pytest.raises(ValueError, match="non-finite"):
        unfold.sandii(response, [1.0, float("nan")], guess)
    with pytest.raises(ValueError, match="negative entry"):
        unfold.sandii([[1.0, -1.0, 1.0], [1.0, 1.0, 1.0]], rates, guess)
    with pytest.raises(ValueError, match="strictly positive"):
        unfold.sandii(response, rates, [1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="tolerance"):
        unfold.sandii(response, rates, guess, tolerance=0.0)
    with pytest.raises(ValueError, match="max_iterations"):
        unfold.sandii(response, rates, guess, max_iterations=0)
    with pytest.raises(ValueError, match="unreachable"):
        unfold.sandii([[0.0, 0.0, 0.0]], [1.0], guess)
    with pytest.raises(ValueError, match="shape mismatch"):
        unfold.forward_fold(response, [1.0])


def test_staysl_exact_guess_is_a_fixed_point() -> None:
    response = _synthetic_response(4, 12, 2.0)
    midpoints = _log_midpoints(12, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    out = unfold.staysl(response, rates, sigmas, truth)
    assert out["iterations"] == 1
    assert out["max_rel_change"] == 0.0
    assert out["spectrum"] == truth
    assert out["rate_factors"] == pytest.approx([1.0] * len(rates), rel=1e-12)


def test_staysl_recovers_determined_system_from_biased_guess() -> None:
    response = _synthetic_response(N_GROUPS, N_GROUPS, 0.9)
    midpoints = _log_midpoints(N_GROUPS, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [3.0 * p for p in truth]
    out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-12)
    assert out["spectrum"] == pytest.approx(truth, rel=1e-6)
    assert out["rate_factors"] == pytest.approx([1.0] * N_GROUPS, rel=1e-6)


def test_staysl_caller_sigmas_are_the_weights() -> None:
    # Two redundant readings of the same group, one precise (2.0 +/- 0.01) and
    # one sloppy (2.2 +/- 1.0): the weighted fit must follow the precise
    # detector, not the average.
    response = [[1.0, 0.0], [1.0, 0.0]]
    rates = [2.0, 2.2]
    sigmas = [0.01, 1.0]
    guess = [1.5, 7.0]
    out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-12, damping=1e-8)
    assert out["spectrum"][0] == pytest.approx(2.0, abs=0.02)
    # The unseen second group keeps the guess exactly.
    assert out["spectrum"][1] == 7.0


def test_staysl_underdetermined_run_reproduces_the_rates() -> None:
    # Fewer detectors than groups: the least-squares fixed point interpolates
    # the measured rates. Unlike SAND-II's multiplicative spreading it does
    # NOT guarantee a spectrum closer to the truth (weakly constrained group
    # directions take their values from the fit, not the prior), so the gate
    # pins rate reproduction only.
    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]
    out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-6)
    assert out["rate_factors"] == pytest.approx([1.0] * 6, rel=1e-4)
    assert all(math.isfinite(p) for p in out["spectrum"])


def test_staysl_irdff_ii_analytical_benchmark_shapes_recover() -> None:
    # Same published shapes as the SAND-II probe (Trkov et al., NDS 163
    # (2020) 1, arXiv:1909.03336), with group bounds bracketing each shape's
    # support so every group sees a significant flux.
    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K, eV -> MeV
    shapes = [
        [e * math.exp(-e / thermal_kt) + 1e-30 for e in _log_midpoints(18, 1e-9, 1e-3)],
        [1.0 / e + 1e-30 for e in _log_midpoints(18, 1e-6, 10.0)],
        [math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in _log_midpoints(18, 1e-3, 0.3)],
    ]
    for shape in shapes:
        response = _synthetic_response(18, 18, 0.8)
        rates = unfold.forward_fold(response, shape)
        sigmas = [1.0] * len(rates)
        guess = [5.0 * p for p in shape]
        out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-10)
        assert out["spectrum"] == pytest.approx(shape, rel=1e-5)


def test_staysl_sandii_fixed_point_is_a_least_squares_fixed_point() -> None:
    # Both methods consume the same fold: once SAND-II has converged (its
    # rate factors pinned to 1), the anchored least-squares cycle about the
    # SAND-II spectrum is a no-op — one confirming solve.
    response = _synthetic_response(6, 6, 1.0)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.5 * p for p in truth]
    sandii_out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=100_000)
    sigmas = [1.0] * len(rates)
    out = unfold.staysl(
        response,
        rates,
        sigmas,
        sandii_out["spectrum"],
        tolerance=1e-4,
        damping=1e-6,
    )
    assert out["iterations"] == 1
    assert out["max_rel_change"] < 1e-4
    assert out["spectrum"] == pytest.approx(sandii_out["spectrum"], rel=1e-5)


def test_staysl_zero_measurement_pulls_toward_zero_and_clips() -> None:
    response = [[1.0, 1.0]]
    rates = [0.0]
    sigmas = [1.0]
    guess = [2.0, 3.0]
    out = unfold.staysl(response, rates, sigmas, guess)
    assert out["spectrum"][0] == 0.0
    assert 0.0 < out["spectrum"][1] < 1e-2
    assert out["rates"][0] == pytest.approx(0.0, abs=1e-2)
    assert out["rate_factors"] == [0.0]


def test_staysl_non_convergence_is_a_hard_fail() -> None:
    response = [[1.0, 0.0, 0.0]]
    rates = [4.0]
    sigmas = [1.0]
    guess = [2.0, 7.0, 3.0]
    with pytest.raises(ValueError, match="did not converge"):
        unfold.staysl(response, rates, sigmas, guess, max_iterations=1)


def test_staysl_validation_errors_name_the_cause() -> None:
    response = _synthetic_response(2, 3, 1.0)
    rates = [1.0, 2.0]
    sigmas = [1.0, 1.0]
    guess = [1.0, 1.0, 1.0]
    with pytest.raises(ValueError, match="shape mismatch"):
        unfold.staysl(response, rates, [1.0], guess)
    with pytest.raises(ValueError, match="non-finite sigma"):
        unfold.staysl(response, rates, [1.0, float("nan")], guess)
    with pytest.raises(ValueError, match="sigma must be > 0"):
        unfold.staysl(response, rates, [1.0, 0.0], guess)
    with pytest.raises(ValueError, match="strictly positive"):
        unfold.staysl(response, rates, sigmas, [1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="tolerance"):
        unfold.staysl(response, rates, sigmas, guess, tolerance=0.0)
    with pytest.raises(ValueError, match="max_iterations"):
        unfold.staysl(response, rates, sigmas, guess, max_iterations=0)
    with pytest.raises(ValueError, match="damping"):
        unfold.staysl(response, rates, sigmas, guess, damping=0.0)
    with pytest.raises(ValueError, match="unreachable"):
        unfold.staysl([[0.0, 0.0, 0.0]], [1.0], [1.0], guess)


def test_gravel_exact_guess_is_a_fixed_point() -> None:
    response = _synthetic_response(4, 12, 2.0)
    midpoints = _log_midpoints(12, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    out = unfold.gravel(response, rates, sigmas, truth)
    assert out["iterations"] == 1
    assert out["max_rel_change"] == 0.0
    assert out["spectrum"] == truth
    assert out["rate_factors"] == pytest.approx([1.0] * len(rates), rel=1e-12)


def test_gravel_recovers_determined_system_from_biased_guess() -> None:
    response = _synthetic_response(N_GROUPS, N_GROUPS, 0.9)
    midpoints = _log_midpoints(N_GROUPS, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [3.0 * p for p in truth]
    out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=10_000)
    assert out["spectrum"] == pytest.approx(truth, rel=1e-6)
    assert out["rate_factors"] == pytest.approx([1.0] * N_GROUPS, rel=1e-9)


def test_gravel_caller_sigmas_are_the_weights() -> None:
    # Two redundant readings of the same group, one precise (2.0 +/- 0.01) and
    # one sloppy (2.2 +/- 1.0): the N^2/sigma^2 factors (40000 vs 4.84) must
    # follow the precise detector, not the average.
    response = [[1.0, 0.0], [1.0, 0.0]]
    rates = [2.0, 2.2]
    sigmas = [0.01, 1.0]
    guess = [1.5, 7.0]
    out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-12)
    assert out["spectrum"][0] == pytest.approx(2.0, abs=0.02)
    # The unseen second group keeps the guess exactly.
    assert out["spectrum"][1] == 7.0


def test_gravel_underdetermined_round_trip_reproduces_rates() -> None:
    # Fewer detectors than groups with counting-statistics sigmas
    # (sigma^2 = N): the chi-square-weighted fixed point reproduces the
    # rates and improves on the guess.
    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [math.sqrt(r) for r in rates]
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]

    def max_rel_err(a: list[float], b: list[float]) -> float:
        return max(abs(x - y) / y for x, y in zip(a, b, strict=True))

    out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-9, max_iterations=50_000)
    assert out["rate_factors"] == pytest.approx([1.0] * 6, rel=1e-6)
    assert max_rel_err(out["spectrum"], truth) < max_rel_err(guess, truth)
    assert all(p > 0.0 for p in out["spectrum"])


def test_gravel_irdff_ii_analytical_benchmark_shapes_recover() -> None:
    # Same published shapes as the SAND-II probe (Trkov et al., NDS 163
    # (2020) 1, arXiv:1909.03336). The library's tabulated group spectra
    # are IAEA-copyright data files and are deliberately NOT used.
    n_groups = 6
    midpoints = _log_midpoints(n_groups, 1e-9, 20.0)
    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K, eV -> MeV
    shapes = {
        "thermal": [e * math.exp(-e / thermal_kt) + 1e-30 for e in midpoints],
        "1/E": [1.0 / e + 1e-30 for e in midpoints],
        "25keV": [math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in midpoints],
    }
    for shape in shapes.values():
        response = _synthetic_response(n_groups, n_groups, 0.8)
        rates = unfold.forward_fold(response, shape)
        sigmas = [1.0] * len(rates)
        guess = [5.0 * p for p in shape]
        out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-11, max_iterations=100_000)
        assert out["spectrum"] == pytest.approx(shape, rel=1e-5)


def test_gravel_sandii_fixed_points_agree() -> None:
    # Both methods consume the same fold: at a fixed point every rate ratio
    # is 1, so the weights are irrelevant and each method holds the other's
    # fixed point after one confirming no-op adjustment.
    response = _synthetic_response(6, 6, 1.0)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [2.5 * p for p in truth]
    sandii_out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=100_000)
    out = unfold.gravel(response, rates, sigmas, sandii_out["spectrum"], tolerance=1e-4)
    assert out["iterations"] == 1
    assert out["max_rel_change"] < 1e-4
    assert out["spectrum"] == pytest.approx(sandii_out["spectrum"], rel=1e-5)
    gravel_out = unfold.gravel(
        response, rates, sigmas, guess, tolerance=1e-12, max_iterations=100_000
    )
    back = unfold.sandii(response, rates, gravel_out["spectrum"], tolerance=1e-4)
    assert back["iterations"] == 1
    assert back["max_rel_change"] < 1e-4
    assert back["spectrum"] == pytest.approx(gravel_out["spectrum"], rel=1e-5)


def test_gravel_zero_measurement_carries_no_weight() -> None:
    # Pinned divergence from SAND-II: the N^2/sigma^2 factor of a zero
    # measurement is zero, so the detector is skipped rather than pinning
    # its groups to zero — the guess survives untouched.
    out = unfold.gravel([[1.0, 1.0]], [0.0], [1.0], [2.0, 3.0])
    assert out["spectrum"] == [2.0, 3.0]
    assert out["iterations"] == 1
    assert out["max_rel_change"] == 0.0


def test_gravel_non_convergence_is_a_hard_fail() -> None:
    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [2.0 * p for p in truth]
    with pytest.raises(ValueError, match="did not converge"):
        unfold.gravel(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=1)


def test_gravel_validation_errors_name_the_cause() -> None:
    response = _synthetic_response(2, 3, 1.0)
    rates = [1.0, 2.0]
    sigmas = [1.0, 1.0]
    guess = [1.0, 1.0, 1.0]
    with pytest.raises(ValueError, match="shape mismatch"):
        unfold.gravel(response, rates, [1.0], guess)
    with pytest.raises(ValueError, match="non-finite sigma"):
        unfold.gravel(response, rates, [1.0, float("nan")], guess)
    with pytest.raises(ValueError, match="sigma must be > 0"):
        unfold.gravel(response, rates, [1.0, 0.0], guess)
    with pytest.raises(ValueError, match="strictly positive"):
        unfold.gravel(response, rates, sigmas, [1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="tolerance"):
        unfold.gravel(response, rates, sigmas, guess, tolerance=0.0)
    with pytest.raises(ValueError, match="max_iterations"):
        unfold.gravel(response, rates, sigmas, guess, max_iterations=0)
    with pytest.raises(ValueError, match="unreachable"):
        unfold.gravel([[0.0, 0.0, 0.0]], [1.0], [1.0], guess)
