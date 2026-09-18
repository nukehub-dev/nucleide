"""Python-side tests for the ECRH accessibility kernel (hand vectors + errors)."""

import pytest

import nucleide.plasma_source as ps

# Hand-computed vectors from the SI constants (e, m_e, eps0):
# f_ce/B = 27.99248987233304 GHz/T, O1 coef = 1.2404426061150442e16 m^-3/GHz^2.
FREQ_GHZ = 170.0
B_COLD_T = 6.073057479892957
B_REL_T = 6.310751020354947  # 20 keV Maxwell-Juttner shift
N_O1 = 3.584879131672478e20
N_X1_AT_5T = 6.334175791171868e19


def test_ecrh_scalars_match_hand_vectors() -> None:
    out = ps.ecrh_scalars(FREQ_GHZ, field_t=5.0, electron_temperature_kev=20.0)
    assert out["f_ce_per_t_ghz"] == pytest.approx(27.99248987233304, rel=1e-9)
    assert out["f_ce_ghz"] == pytest.approx(139.9624493616652, rel=1e-9)
    assert out["b_cold_t"] == pytest.approx(B_COLD_T, rel=1e-9)
    assert out["b_rel_t"] == pytest.approx(B_REL_T, rel=1e-9)
    assert out["n_o1_m3"] == pytest.approx(N_O1, rel=1e-9)
    assert out["n_x1_m3"] == pytest.approx(N_X1_AT_5T, rel=1e-9)


def test_ecrh_scalars_second_harmonic_and_evanescent() -> None:
    out = ps.ecrh_scalars(FREQ_GHZ, harmonic=2)
    assert out["b_cold_t"] == pytest.approx(B_COLD_T / 2.0, rel=1e-12)
    assert out["b_rel_t"] is None
    assert out["f_ce_ghz"] is None
    # 100 GHz at 5 T sits below f_ce: X1 evanescent, no cut-off.
    evanescent = ps.ecrh_scalars(100.0, field_t=5.0)
    assert evanescent["n_x1_m3"] is None


def test_ecrh_accessibility_beamline_positions() -> None:
    s = [0.0, 1.0, 2.0]
    b = [7.0, 6.0, 5.0]
    ne = [1.0e20, 3.0e20, 5.0e20]
    out = ps.ecrh_accessibility(s, b, ne, FREQ_GHZ, electron_temperature_kev=20.0)
    assert out["b_res_t"] == pytest.approx(B_REL_T, rel=1e-9)
    assert out["resonance_m"] == pytest.approx([0.689248979645053], rel=1e-9)
    assert out["n_o1_m3"] == pytest.approx(N_O1, rel=1e-9)
    # O1 cut-off at 170 GHz (3.585e20) sits between the last two nodes.
    assert out["o1_cutoff_m"] == pytest.approx([1.292439565836239], rel=1e-9)
    # Cold (no-temperature) call locates the unshifted resonance instead.
    cold = ps.ecrh_accessibility(s, b, ne, FREQ_GHZ)
    assert cold["b_res_t"] == pytest.approx(B_COLD_T, rel=1e-9)
    assert cold["resonance_m"] == pytest.approx([0.9269425201070431], rel=1e-9)


def test_ecrh_malformed_inputs_are_loud() -> None:
    with pytest.raises(ValueError):
        ps.ecrh_scalars(0.0)
    with pytest.raises(ValueError):
        ps.ecrh_scalars(FREQ_GHZ, harmonic=0)
    with pytest.raises(ValueError):
        ps.ecrh_scalars(FREQ_GHZ, electron_temperature_kev=-1.0)
    with pytest.raises(ValueError):
        ps.ecrh_accessibility([], [], [], FREQ_GHZ)
    with pytest.raises(ValueError):
        ps.ecrh_accessibility([0.0, 1.0], [5.0], [1e19, 1e20], FREQ_GHZ)
    with pytest.raises(ValueError):
        ps.ecrh_accessibility([1.0, 0.0], [5.0, 5.0], [1e19, 1e20], FREQ_GHZ)
    with pytest.raises(ValueError):
        ps.ecrh_accessibility([0.0, 1.0], [5.0, 5.0], [1e19, -1.0], FREQ_GHZ)
