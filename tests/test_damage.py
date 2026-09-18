"""Python-side tests for the damage/gas metric core (analytic gates + errors).

The constants of the Lindhard partition (Robinson fit) and the NRT/arc-dpa
piecewise forms are transcribed again, independently, in this file — the
same two-spelling discipline as the plasma-source Ballabio checks. Inputs
are synthetic hand-picked values throughout; the SPECTER report spots live
in the validation harness, not here — except the vendored-Table-VII checks
at the end, which pin the committed transcription (same stance as the EU
Annex VII table tests), not the report.
"""

import math

import pytest

import nucleide.damage as dmg

ED_EV = 40.0  # synthetic Fe-like average threshold displacement energy
KAPPA = 0.8  # NRT displacement efficiency (NRT 1975)


def _lindhard_partition_py(t_ev: float, z: float, a: float) -> float:
    """Independent transcription: Robinson fit, self-recoil closed forms.

    Self-recoil closed forms E_L = 30.724·Z²·(2·Z^{2/3})^{1/2}·2 and
    k_L = 0.0793·2^{3/4}·Z^{2/3}/A^{1/2} (Lindhard 1963; Robinson 1970).
    ``math.pow`` rather than ``**`` — this repo's mypy types float**float as Any.
    """
    e_l = 30.724 * z * z * math.pow(2.0 * math.pow(z, 2.0 / 3.0), 0.5) * 2.0
    k_l = 0.0793 * math.pow(2.0, 0.75) * math.pow(z, 2.0 / 3.0) / math.pow(a, 0.5)
    eps = t_ev / e_l
    g = eps + 3.4008 * math.pow(eps, 1.0 / 6.0) + 0.40244 * math.pow(eps, 0.75)
    return 1.0 / (1.0 + k_l * g)


def _saturation_limit_py(z: float, a: float) -> float:
    """E_L/k_L rearranged term-by-term (shares no intermediate with above)."""
    return (
        (30.724 / 0.0793)
        * math.pow(z, 1.0 / 3.0)
        * math.pow(z, 0.5)
        * math.pow(2.0 * math.pow(z, 2.0 / 3.0), 1.25)
        * math.pow(a, 1.5)
        / (math.pow(a, 0.5) * math.pow(2.0 * a, 0.5))
    )


def test_lindhard_partition_matches_independent_transcription() -> None:
    # Fe-56 self-recoil across 12 decades of PKA energy.
    for t in [1e-3, 1.0, 1e2, 1e4, 1e5, 1e6, 1e8, 1e10, 1e12]:
        got = dmg.lindhard_partition(t, "Fe56", "Fe56")
        want = _lindhard_partition_py(t, 26.0, 56.0)
        assert got == pytest.approx(want, rel=1e-12)
    # int nucid spelling agrees with the name spelling.
    fe56_nucid = 260560000
    assert dmg.lindhard_partition(1e5, fe56_nucid, fe56_nucid) == pytest.approx(
        dmg.lindhard_partition(1e5, "Fe56", "Fe56"), rel=0, abs=0.0
    )
    # At a fixed PKA energy a heavier lattice keeps the larger damage
    # fraction (its reduced energy ε is smaller: E_L scales ~Z^{7/3}).
    assert dmg.lindhard_partition(1e5, "W184", "W184") > dmg.lindhard_partition(1e5, "C12", "C12")


def test_partition_limits() -> None:
    # P -> 1 as T -> 0 (logarithmically: g(ε) ~ ε^{1/6}).
    assert 0.99 < dmg.lindhard_partition(1e-6, "Fe56", "Fe56") < 1.0
    # Monotone decrease.
    prev = 1.0
    for t in [1e0, 1e2, 1e4, 1e6, 1e8, 1e10]:
        p = dmg.lindhard_partition(t, "Fe56", "Fe56")
        assert 0.0 < p < prev
        prev = p
    # Damage energy saturates at E_L/k_L (asymptote, within 1% at 1e12 eV).
    t_dam = dmg.damage_energy(1e12, "Fe56", "Fe56")
    assert t_dam == pytest.approx(_saturation_limit_py(26.0, 56.0), rel=1e-2)
    # Never above the PKA energy.
    for t in [0.5, 100.0, 1e5, 1e8, 1e11]:
        assert dmg.damage_energy(t, "Fe56", "Fe56") <= t


def test_nrt_displacements_piecewise() -> None:
    # Below threshold / plateau / high branch against the transcribed form.
    assert dmg.nrt_displacements(ED_EV * 0.5, ED_EV, "Fe56") == 0.0
    assert dmg.nrt_displacements(ED_EV, ED_EV, "Fe56") == 1.0
    assert dmg.nrt_displacements(2.0 * ED_EV, ED_EV, "Fe56") == 1.0
    branch = 2.0 * ED_EV / KAPPA
    assert dmg.nrt_displacements(branch * (1.0 - 1e-12), ED_EV, "Fe56") == 1.0
    t = 1.0e5
    t_dam = dmg.damage_energy(t, "Fe56", "Fe56")
    assert dmg.nrt_displacements(t, ED_EV, "Fe56") == pytest.approx(
        KAPPA * t_dam / (2.0 * ED_EV), rel=1e-15
    )


def test_arc_efficiency_construction() -> None:
    b_arc, c_arc = -0.55, 0.3
    junction = 2.0 * ED_EV / KAPPA
    # ξ(2E_d/κ) = 1 exactly (continuity construction, Nordlund 2018 Eq. 7).
    assert dmg.arc_efficiency(junction, ED_EV, b_arc, c_arc) == 1.0
    # Saturation to c_arc and monotone decrease.
    assert dmg.arc_efficiency(1e9, ED_EV, b_arc, c_arc) == pytest.approx(c_arc, abs=1e-3)
    xs = [1.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    xi = [dmg.arc_efficiency(t, ED_EV, b_arc, c_arc) for t in xs]
    assert all(xi[i + 1] < xi[i] for i in range(len(xi) - 1))
    # High branch: NRT × ξ(T_dam); plateau: identical to NRT.
    t = 1.0e5
    t_dam = dmg.damage_energy(t, "Fe56", "Fe56")
    want = KAPPA * t_dam * dmg.arc_efficiency(t_dam, ED_EV, b_arc, c_arc) / (2.0 * ED_EV)
    assert dmg.arc_displacements(t, ED_EV, "Fe56", b_arc, c_arc) == pytest.approx(want, rel=1e-15)
    assert dmg.arc_displacements(ED_EV * 1.5, ED_EV, "Fe56", b_arc, c_arc) == 1.0
    assert dmg.arc_displacements(ED_EV * 0.5, ED_EV, "Fe56", b_arc, c_arc) == 0.0


def test_folds_match_hand_sums() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    flux = [1.0e12, 2.0e12, 4.0e12]
    resp = [100.0, 200.0, 50.0]
    hand = sum(f * r for f, r in zip(flux, resp, strict=True))
    assert dmg.nrt_dpa(flux, resp, bounds, 2.0) == pytest.approx(2.0e-24 * hand, rel=1e-15)
    assert dmg.arc_dpa(flux, resp, bounds, 2.0) == pytest.approx(2.0e-24 * hand, rel=1e-15)
    assert dmg.gas_appm(flux, resp, bounds, 2.0) == pytest.approx(2.0e-18 * hand, rel=1e-15)
    # Piecewise-constant semantics: bounds width never enters the sum.
    wide = [0.0, 5.0, 10.0, 20.0]
    assert dmg.nrt_dpa(flux, resp, wide, 2.0) == pytest.approx(2.0e-24 * hand, rel=1e-15)


def test_zero_flux_groups_and_zero_dpa() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    resp = [1.0e3, 1.5e2, 7.0e3]
    live = [0.0, 2.0e12, 0.0]
    assert dmg.nrt_dpa(live, resp, bounds, 1.0) == pytest.approx(1.0e-24 * 2.0e12 * 1.5e2)
    dead = [0.0, 0.0, 0.0]
    assert dmg.nrt_dpa(dead, resp, bounds, 1.0) == 0.0
    with pytest.raises(ValueError, match="zero dpa"):
        dmg.he_dpa_ratio(dead, resp, resp, bounds, 1.0)
    with pytest.raises(ValueError, match="zero dpa"):
        dmg.he_dpa_ratio(live, resp, [0.0, 0.0, 0.0], bounds, 1.0)


def test_he_dpa_ratio_uses_one_flux_fold() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    flux = [1.0e13, 1.0e12, 1.0e11]
    he = [0.2, 0.4, 0.6]
    dpa_xs = [50.0, 100.0, 200.0]
    ratio = dmg.he_dpa_ratio(flux, he, dpa_xs, bounds, 3.0)
    assert ratio == pytest.approx(
        dmg.gas_appm(flux, he, bounds, 3.0) / dmg.nrt_dpa(flux, dpa_xs, bounds, 3.0), rel=1e-15
    )


def test_fold_uq_gate_and_determinism() -> None:
    flux = [1.0e13, 3.0e12]
    resp = [80.0, 160.0]
    bounds = [0.0, 0.5, 20.0]
    mean = [0.0, 0.0, 0.0, 0.0]
    cov = [
        [0.0004, 0.0, 0.0, 0.0],
        [0.0, 0.0009, 0.0, 0.0],
        [0.0, 0.0, 0.0025, 0.0],
        [0.0, 0.0, 0.0, 0.0001],
    ]
    for metric in ("nrt_dpa", "arc_dpa", "gas_appm"):
        out = dmg.fold_uq(metric, flux, resp, bounds, 2.0, mean, cov, 20_000, 20260915, 5.0)
        assert out["passed"] is True
        assert out["metric"] == metric
        assert out["n"] == 20_000
        assert out["expected"] == pytest.approx(out["nominal"], rel=1e-15)
        mean_se = out["analytic_std"] / math.sqrt(20_000)
        assert abs(out["mean"] - out["expected"]) <= 5.0 * mean_se
        std_se = out["analytic_std"] / math.sqrt(2.0 * (20_000 - 1))
        assert abs(out["std"] - out["analytic_std"]) <= 5.0 * std_se
        again = dmg.fold_uq(metric, flux, resp, bounds, 2.0, mean, cov, 20_000, 20260915, 5.0)
        assert again == out
    # The ratio UQ is a loud named-open.
    with pytest.raises(ValueError, match="not yet supported"):
        dmg.fold_uq("he_dpa_ratio", flux, resp, bounds, 1.0, mean, cov, 16, 1, 5.0)
    # Unknown metric names are rejected loudly too.
    with pytest.raises(ValueError, match="unknown fold metric"):
        dmg.fold_uq("dpa", flux, resp, bounds, 1.0, mean, cov, 16, 1, 5.0)


def test_error_paths_name_their_cause() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    flux = [1.0e12, 2.0e12, 4.0e12]
    resp = [100.0, 200.0, 50.0]
    with pytest.raises(ValueError, match="response.*expected 3, got 2"):
        dmg.nrt_dpa(flux, [1.0, 2.0], bounds, 1.0)
    with pytest.raises(ValueError, match="bounds.*expected 4, got 2"):
        dmg.nrt_dpa(flux, resp, [0.0, 1.0], 1.0)
    with pytest.raises(ValueError, match="seconds.*strictly positive"):
        dmg.nrt_dpa(flux, resp, bounds, 0.0)
    with pytest.raises(ValueError, match="flux.*non-negative"):
        dmg.nrt_dpa([-1.0, 1.0, 1.0], resp, bounds, 1.0)
    with pytest.raises(ValueError, match="response.*non-negative"):
        dmg.nrt_dpa(flux, [1.0, -1.0, 1.0], bounds, 1.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        dmg.nrt_dpa(flux, resp, [0.0, 0.5, 0.5, 1.0], 1.0)
    with pytest.raises(ValueError, match="non-finite"):
        dmg.nrt_dpa([float("nan"), 1.0, 1.0], resp, bounds, 1.0)
    with pytest.raises(ValueError, match="pka_energy_ev"):
        dmg.nrt_displacements(-1.0, ED_EV, "Fe56")
    with pytest.raises(ValueError, match="ed_ev"):
        dmg.nrt_displacements(100.0, 0.0, "Fe56")
    with pytest.raises(ValueError, match="b_arc"):
        dmg.arc_efficiency(100.0, ED_EV, 0.5, 0.3)
    with pytest.raises(ValueError, match="c_arc"):
        dmg.arc_efficiency(100.0, ED_EV, -0.5, 1.2)
    with pytest.raises(ValueError):
        dmg.nrt_displacements(100.0, ED_EV, 999_999_999)


def test_specter_table_row_count_and_spectra() -> None:
    # The committed Table VII transcription pins 24 element rows; the
    # spectrum vocabulary is the seven report columns.
    assert dmg.specter_spectra() == [
        "thermal",
        "fission",
        "14mev",
        "hfir",
        "ebr2",
        "fftf",
        "fusion",
    ]
    for spectrum in dmg.specter_spectra():
        table = dmg.specter_table(spectrum)
        assert len(table) == 24
        assert len(dmg.specter_damage_energy(spectrum)) == 24


def test_specter_transcription_spots() -> None:
    # Hand-read from the report scan (Table VII p. -31-, Table II p. -14-).
    hfir = dmg.specter_table("hfir")
    assert hfir["Fe"] == pytest.approx(19.1 * 0.8 / (2.0 * 0.040), rel=1e-15)
    assert hfir["Ti"] == pytest.approx(21.9 * 0.8 / (2.0 * 0.040), rel=1e-15)
    assert hfir["Cu"] == pytest.approx(18.6 * 0.8 / (2.0 * 0.040), rel=1e-15)
    energy = dmg.specter_damage_energy("hfir")
    assert energy["Fe"] == 19.1
    fusion = dmg.specter_damage_energy("fusion")
    assert fusion["Ni"] == 109.4
    assert fusion["W"] == 68.0
    assert dmg.specter_damage_energy("14mev")["Ni"] == 300.0
    assert dmg.specter_damage_energy("thermal")["Co"] == 13.38
    assert dmg.specter_damage_energy("fission")["V"] == 101.0
    # Table II E_d spots (eV).
    for element, expected in [
        ("Be", 31.0),
        ("Al", 27.0),
        ("Fe", 40.0),
        ("Mo", 60.0),
        ("Ag", 60.0),
        ("Ta", 53.0),
        ("W", 90.0),
        ("Au", 30.0),
    ]:
        assert dmg.specter_ed(element) == expected


def test_specter_hfir_reproduces_oracle_spots_at_table_print_precision() -> None:
    # Table VII prints 3 significant figures, so the vendored HFIR column
    # reproduces the validation harness's 5-digit Table VI spots
    # (Fe 191.18, Ti 218.73, Cu 186.49 barns) at 5e-3 relative — the
    # finest agreement the transcribed inputs allow.
    hfir = dmg.specter_table("hfir")
    for element, printed in [("Fe", 191.18), ("Ti", 218.73), ("Cu", 186.49)]:
        assert hfir[element] == pytest.approx(printed, rel=5e-3)
    # Same through the one-group fluence fold the oracle runs.
    fluence, bounds = [4.78373e22], [0.0, 20.0]
    for element, printed_dpa in [("Fe", 9.1455), ("Ti", 10.464), ("Cu", 8.9212)]:
        got = dmg.nrt_dpa(fluence, [hfir[element]], bounds, 1.0)
        assert got == pytest.approx(printed_dpa, rel=5e-3)


def test_specter_fold_equals_caller_slice_fold() -> None:
    # The table is just another caller input: identical floats fold
    # identically, matching the hand product.
    xs = dmg.specter_table("fusion")["Fe"]
    assert dmg.nrt_dpa([3.0e14], [xs], [0.0, 20.0], 1.0) == pytest.approx(
        1.0e-24 * 3.0e14 * xs, rel=1e-15
    )


def test_specter_unknown_keys_are_loud_errors() -> None:
    with pytest.raises(ValueError, match="unknown SPECTER table spectrum"):
        dmg.specter_table("pwr")
    with pytest.raises(ValueError, match="unknown SPECTER table element"):
        dmg.specter_ed("U")
    with pytest.raises(ValueError, match="unknown SPECTER table spectrum"):
        dmg.specter_damage_energy("HFIR-CTR32")


def test_coil_fast_flux_threshold_pins_groups() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    flux = [1.0e12, 2.0e12, 4.0e12]
    assert dmg.coil_fast_flux(flux, bounds, 0.0) == pytest.approx(7.0e12, rel=1e-15)
    assert dmg.coil_fast_flux(flux, bounds, 0.1) == pytest.approx(6.0e12, rel=1e-15)
    assert dmg.coil_fast_flux(flux, bounds, 1.0) == pytest.approx(4.0e12, rel=1e-15)
    assert dmg.coil_fast_flux(flux, bounds, 20.0) == 0.0
    # Threshold cutting through a group includes the whole group.
    assert dmg.coil_fast_flux(flux, bounds, 0.5) == pytest.approx(6.0e12, rel=1e-15)
    rate = dmg.coil_fast_flux(flux, bounds, 0.1)
    assert dmg.coil_fast_fluence(flux, bounds, 0.1, 2.0) == pytest.approx(rate * 2.0, rel=1e-15)


def test_coil_stellaris_anchor_is_ten_fpy_class() -> None:
    # Caller-supplied design numbers (gates, not defaults): 99th-percentile
    # coil fast flux 9.5e13 m^-2 s^-1 against fast-fluence limits
    # 3e22 / 1.5e23 m^-2 -> ~10 FPY, limited by the first channel.
    life = dmg.coil_lifetime([3.0e22, 1.5e23], [9.5e13, 9.5e13])
    assert life["limiting"] == 0
    assert life["seconds"] == pytest.approx(3.0e22 / 9.5e13, rel=1e-9)
    assert life["seconds"] / dmg.FPY_SECONDS == pytest.approx(10.0, rel=1e-2)
    fluence = dmg.coil_fast_fluence([9.5e13], [0.0, 20.0], 0.1, 10.0 * dmg.FPY_SECONDS)
    assert fluence == pytest.approx(3.0e22, rel=1e-3)
    cu = dmg.coil_lifetime([1.5e23], [9.5e13])
    assert cu["seconds"] / dmg.FPY_SECONDS == pytest.approx(50.0, rel=1e-2)


def test_coil_weakest_link_and_remaining() -> None:
    assert dmg.coil_lifetime([10.0, 10.0], [1.0, 2.0]) == {"seconds": 5.0, "limiting": 1}
    # Zero rates never breach and are skipped; all zero is infinite life.
    assert dmg.coil_lifetime([10.0, 10.0], [0.0, 2.0]) == {"seconds": 5.0, "limiting": 1}
    idle = dmg.coil_lifetime([10.0, 10.0], [0.0, 0.0])
    assert idle["seconds"] == math.inf and idle["limiting"] is None
    assert dmg.coil_lifetime([0.0, 10.0], [1.0, 1.0]) == {"seconds": 0.0, "limiting": 0}
    # dpa spectral weighting stays in the landed fold: the NRT fold at one
    # second is the rate the life kernel consumes unchanged.
    rate = dmg.nrt_dpa([1.0e12, 2.0e12], [100.0, 200.0], [0.0, 0.5, 20.0], 1.0)
    assert rate == pytest.approx(1.0e-24 * 5.0e14, rel=1e-15)
    assert dmg.coil_lifetime([rate], [rate]) == {"seconds": 1.0, "limiting": 0}
    # Remaining life from accumulated damage, clamped at spent limits.
    assert dmg.coil_remaining([10.0, 10.0], [4.0, 9.0], [2.0, 1.0]) == {
        "seconds": 1.0,
        "limiting": 1,
    }
    assert dmg.coil_remaining([10.0, 10.0], [11.0, 5.0], [1.0, 1.0]) == {
        "seconds": 0.0,
        "limiting": 0,
    }
    # History accumulation is rate x time per interval.
    assert dmg.coil_accumulate([1.0e13, 2.0e13], [10.0, 5.0]) == pytest.approx(2.0e14, rel=1e-15)


def test_coil_error_paths_name_their_cause() -> None:
    bounds = [0.0, 0.1, 1.0, 20.0]
    flux = [1.0e12, 2.0e12, 4.0e12]
    with pytest.raises(ValueError, match="bounds.*expected 4, got 2"):
        dmg.coil_fast_flux(flux, [0.0, 1.0], 0.1)
    with pytest.raises(ValueError, match="threshold_mev"):
        dmg.coil_fast_flux(flux, bounds, -0.1)
    with pytest.raises(ValueError, match="seconds.*strictly positive"):
        dmg.coil_fast_fluence(flux, bounds, 0.1, 0.0)
    with pytest.raises(ValueError, match="flux.*non-negative"):
        dmg.coil_fast_flux([-1.0, 1.0, 1.0], bounds, 0.1)
    with pytest.raises(ValueError, match="strictly increasing"):
        dmg.coil_fast_flux(flux, [0.0, 0.5, 0.5, 1.0], 0.1)
    with pytest.raises(ValueError, match="durations.*expected 2, got 1"):
        dmg.coil_accumulate([1.0, 2.0], [1.0])
    with pytest.raises(ValueError, match="rates.*expected 2, got 1"):
        dmg.coil_lifetime([1.0, 2.0], [1.0])
    with pytest.raises(ValueError, match="rates.*non-negative"):
        dmg.coil_lifetime([1.0], [-1.0])
    with pytest.raises(ValueError, match="accumulated.*expected 1, got 2"):
        dmg.coil_remaining([1.0], [0.5, 0.5], [1.0])
    with pytest.raises(ValueError, match="accumulated.*non-negative"):
        dmg.coil_remaining([1.0], [-0.5], [1.0])
