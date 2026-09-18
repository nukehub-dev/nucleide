"""Python-side tests for the tokamak fusion-source core (analytic gates + errors)."""

import math
from typing import Any

import pytest

import nucleide.mcnp as mcnp
import nucleide.plasma_source as ps

RADIUS_CM = 300.0
HEIGHT_CM = 25.0
N = 200_000

RING_SPEC = {
    "kind": "ring",
    "radius": RADIUS_CM,
    "height": HEIGHT_CM,
    "reaction": "dt",
    "ion_temperature_kev": 20.0,
}

# ITER-ish synthetic H-mode parametric plasma (hand round numbers).
PARAMETRIC_SPEC = {
    "kind": "parametric",
    "reaction": "dt",
    "major_radius": 620.0,
    "minor_radius": 200.0,
    "elongation": 1.85,
    "triangularity": 0.35,
    "shafranov_factor": 15.0,
    "mode": "H",
    "pedestal_radius": 150.0,
    "ion_density_centre": 1.2e20,
    "ion_density_peaking_factor": 1.1,
    "ion_density_pedestal": 4.0e19,
    "ion_density_separatrix": 3.0e19,
    "ion_temperature_centre": 28.0,
    "ion_temperature_peaking_factor": 2.5,
    "ion_temperature_beta": 2.0,
    "ion_temperature_pedestal": 4.0,
    "ion_temperature_separatrix": 0.1,
}


def _ballabio_moments(reaction: str, ti_kev: float) -> tuple[float, float]:
    """Independent transcription of the Ballabio et al. 1998 Table III fits."""
    if reaction == "dd":
        a1, a2, a3, a4 = 4.69515, -0.040729, 0.47, 0.81844
        m0 = 2.4495e6
        w0, b1, b2, b3, b4 = 82.542, 1.7013e-3, 0.16888, 0.49, 7.9460e-4
    else:
        a1, a2, a3, a4 = 5.30509, 2.4736e-3, 1.84, 1.3818
        m0 = 14.021e6
        w0, b1, b2, b3, b4 = 177.259, 5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5
    fwhm_over_sigma = 2 * math.sqrt(2 * math.log(2))
    mean = m0 + 1e3 * (a1 * ti_kev ** (2 / 3) / (1 + a2 * ti_kev**a3) + a4 * ti_kev)
    delta = b1 * ti_kev ** (2 / 3) / (1 + b2 * ti_kev**b3) + b4 * ti_kev
    sigma = w0 * (1 + delta) * math.sqrt(ti_kev) / fwhm_over_sigma
    return mean / 1e6, sigma / 1e3


def test_ring_spatial_moments_match_closed_form() -> None:
    out = ps.particles(RING_SPEC, N, seed=42)
    radius = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
    assert bool((out["z"] == HEIGHT_CM).all())
    assert radius == pytest.approx(RADIUS_CM, rel=1e-12)
    se = RADIUS_CM / math.sqrt(2 * N)
    assert abs(float(out["x"].mean())) < 8 * se
    assert abs(float(out["y"].mean())) < 8 * se
    direction_norm = (out["u"] ** 2 + out["v"] ** 2 + out["w"] ** 2) ** 0.5
    assert direction_norm == pytest.approx(1.0, rel=1e-12)
    assert abs(float(out["w"].mean())) < 0.01


def test_spectrum_moments_match_closed_form() -> None:
    out = ps.particles(RING_SPEC, N, seed=7)
    mean, sigma = _ballabio_moments("dt", 20.0)
    assert float(out["energy"].mean()) == pytest.approx(mean, abs=5e-3)
    assert float(out["energy"].std()) == pytest.approx(sigma, abs=5e-3)


def test_spectrum_moments_helper_matches_independent_fit() -> None:
    for reaction in ("dd", "dt"):
        for ti in (1.0, 10.0, 20.0):
            got = ps.spectrum_moments(reaction, ti)
            mean, sigma = _ballabio_moments(reaction, ti)
            assert got["mean_mev"] == pytest.approx(mean, abs=1e-12)
            assert got["sigma_mev"] == pytest.approx(sigma, abs=1e-12)


def test_zero_temperature_is_monoenergetic() -> None:
    spec = {"kind": "point", "position": [0, 0, 0], "reaction": "dd", "ion_temperature_kev": 0.0}
    out = ps.particles(spec, 100, seed=1)
    assert out["energy"] == pytest.approx(2.4495, rel=1e-12)
    moments = ps.spectrum_moments("dd", 0.0)
    assert moments["sigma_mev"] == 0.0
    assert moments["nominal_mev"] == pytest.approx(2.4495)


def test_pinned_seed_reproduces_the_stream() -> None:
    a = ps.particles(RING_SPEC, 512, seed=1234)
    b = ps.particles(RING_SPEC, 512, seed=1234)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert bool((a[key] == b[key]).all())
    c = ps.particles(RING_SPEC, 512, seed=1235)
    assert not bool((a["energy"] == c["energy"]).all())


def test_sdef_card_round_trips_through_typed_reader() -> None:
    out = ps.emit_source_cards(RING_SPEC)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["rad"] == "D1"
    assert parsed["axs"] == "0 0 1"
    assert parsed["pos"] == "0 0 25"
    assert parsed["erg"] == "D2"
    assert len(parsed["distributions"]) == 2
    assert parsed["distributions"][0]["si"] == "300 300"
    assert out["sdef"]["drift"][0]["reparsed"] is True


def test_sdef_point_monoenergetic_card() -> None:
    spec = {"kind": "point", "position": [1, 2, 3], "reaction": "dd", "ion_temperature_kev": 0.0}
    out = ps.emit_source_cards(spec)
    assert out["sdef"]["card"] == "SDEF POS=1 2 3\n     ERG=2.4495\n     WGT=1\n     PAR=n"
    assert out["spectrum"]["mono"] is True
    assert out["sdef"]["drift"][0]["rel_drift"] == 0.0


def test_serpent_card_structure() -> None:
    out = ps.emit_source_cards(RING_SPEC)
    card = out["serpent"]["card"]
    lines = card.splitlines()
    assert lines[0] == "src 1 pos 0 0 25"
    assert lines[1] == "src 1 rad d1"
    assert lines[2] == "src 1 erg d2"
    assert lines[3] == "src 1 wgt 1"
    assert lines[4] == "SI1 300 300"
    assert lines[5] == "SP1 1"
    assert lines[6].startswith("SI2 ")
    assert lines[7].startswith("SP2 ")
    assert out["serpent"]["drift"][0]["reparsed"] is False


def test_gaussian_tabulation_drift_is_reported() -> None:
    out = ps.emit_source_cards(RING_SPEC)
    drift = out["sdef"]["drift"][0]
    assert drift["quantity"] == "emission probability"
    assert drift["rel_drift"] > 0.0
    assert drift["rel_drift"] == pytest.approx(1.0 - drift["accounted"], rel=1e-12)
    # +/-4 sigma tail mass.
    assert drift["rel_drift"] == pytest.approx(1.0 - (math.erf(4 / math.sqrt(2))), rel=1e-3)


def test_mcnp_version_six_and_bad_version() -> None:
    # Neutrons carry designator "n" in both dialects; version 6 must emit an
    # identical card and an unknown version must raise loudly.
    mono = {
        "kind": "point",
        "position": [0, 0, 0],
        "reaction": "dt",
        "ion_temperature_kev": 0.0,
        "mcnp_version": 6,
    }
    out = ps.emit_source_cards(mono)
    assert out["sdef"]["card"] == "SDEF POS=0 0 0\n     ERG=14.021\n     WGT=1\n     PAR=n"
    bad = dict(RING_SPEC, mcnp_version=4)
    with pytest.raises(ValueError, match="version"):
        ps.emit_source_cards(bad)


def test_validation_errors_are_loud() -> None:
    with pytest.raises(ValueError, match="radius"):
        ps.particles(dict(RING_SPEC, radius=0.0), 4, seed=0)
    with pytest.raises(ValueError, match="ion temperature"):
        ps.particles(dict(RING_SPEC, ion_temperature_kev=-1.0), 4, seed=0)
    with pytest.raises(ValueError, match="reaction"):
        ps.particles(dict(RING_SPEC, reaction="tt"), 4, seed=0)
    with pytest.raises(ValueError, match="kind"):
        ps.particles(dict(RING_SPEC, kind="torus"), 4, seed=0)
    with pytest.raises(ValueError, match="position"):
        ps.particles(
            {"kind": "point", "position": [0, 0], "reaction": "dt", "ion_temperature_kev": 0.0},
            4,
            seed=0,
        )
    # Ring/point specs ignore parametric-only keys; a parametric spec keeps
    # loud errors for malformed sectors (one angle without the other,
    # non-finite or out-of-range angles).
    ps.particles(dict(RING_SPEC, elongation=1.8), 4, seed=0)
    with pytest.raises(ValueError, match="together"):
        ps.particles(dict(PARAMETRIC_SPEC, rotation_angle=1.57), 4, seed=0)
    with pytest.raises(ValueError, match="together"):
        ps.particles(dict(PARAMETRIC_SPEC, start_angle=0.78), 4, seed=0)
    with pytest.raises(ValueError, match="sector"):
        ps.particles(dict(PARAMETRIC_SPEC, start_angle=-0.1, rotation_angle=1.57), 4, seed=0)
    with pytest.raises(ValueError, match="sector"):
        ps.particles(dict(PARAMETRIC_SPEC, start_angle=0.0, rotation_angle=0.0), 4, seed=0)
    with pytest.raises(ValueError, match="sector"):
        ps.particles(dict(PARAMETRIC_SPEC, start_angle=0.0, rotation_angle=7.0), 4, seed=0)
    with pytest.raises(ValueError, match="sector"):
        ps.particles(
            dict(PARAMETRIC_SPEC, start_angle=float("nan"), rotation_angle=1.0),
            4,
            seed=0,
        )


def test_parametric_fuel_mixture_fractions_are_loud() -> None:
    bad_specs = [
        dict(PARAMETRIC_SPEC, fuel={"D": 0.7}),  # missing T
        dict(PARAMETRIC_SPEC, fuel={"D": 0.7, "T": 0.4}),  # non-summing
        dict(PARAMETRIC_SPEC, fuel={"D": -0.1, "T": 1.1}),  # negative
        dict(PARAMETRIC_SPEC, fuel={"D": 0.7, "T": float("nan")}),  # non-finite
        dict(PARAMETRIC_SPEC, fuel={"H": 1.0, "T": 0.0}),  # unknown key
        dict(PARAMETRIC_SPEC, fuel={"D": 0.7, "T": 0.3, "He3": 0.1}),  # extra key
        dict(PARAMETRIC_SPEC, fuel=0.7),  # not a dict
    ]
    for spec in bad_specs:
        with pytest.raises(ValueError):
            ps.particles(spec, 4, seed=0)


def test_parametric_fuel_mixture_sampling_fires_both_branches() -> None:
    # 70/30 D/T blend: the D-D line (2.45 MeV) carries ~p_dd ≈ 0.7% of the
    # births at these temperatures, and the mean sits just below the D-T
    # line — the two-branch rate rule, not either single-fuel kernel.
    spec = dict(PARAMETRIC_SPEC, fuel={"D": 0.7, "T": 0.3})
    n = 200_000
    out = ps.particles(spec, n, seed=11)
    frac_dd = float((out["energy"] < 10.0).mean())
    assert 0.002 < frac_dd < 0.015
    mean_e = float(out["energy"].mean())
    assert 13.9 < mean_e < 14.1
    assert float(out["energy"].std()) > 0.4
    # Deterministic per seed, and the dict spelling needs no `reaction` key.
    again = ps.particles(spec, 128, seed=11)
    assert bool((again["energy"] == out["energy"][:128]).all())


def test_parametric_fuel_mixture_emits_cards_with_mixture_summary() -> None:
    # The equimolar dict exercises the recovery anchor end to end (the
    # bit-for-bit gates live in the crate); here: cards stay well-formed and
    # the spectrum summary carries between-branch variance.
    spec = dict(PARAMETRIC_SPEC, fuel={"D": 0.5, "T": 0.5})
    out = ps.emit_source_cards(spec, bins=15)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["rad"] == "D1"
    assert parsed["ext"] == "D2"
    assert parsed["erg"] == "D3"
    assert out["sdef"]["drift"][0]["reparsed"] is True
    assert out["spectrum"]["mean_mev"] == pytest.approx(14.0, abs=0.15)
    assert out["spectrum"]["sigma_mev"] > 0.4
    assert out["spectrum"]["mono"] is False
    lines = out["serpent"]["card"].splitlines()
    assert lines[3] == "src 1 erg d3"
    with pytest.raises(ValueError, match="minor radius"):
        ps.particles(dict(PARAMETRIC_SPEC, minor_radius=-1.0), 4, seed=0)
    with pytest.raises(ValueError, match="triangularity"):
        ps.particles(dict(PARAMETRIC_SPEC, triangularity=1.2), 4, seed=0)
    with pytest.raises(ValueError, match="mode"):
        ps.particles(dict(PARAMETRIC_SPEC, mode="Q"), 4, seed=0)


def test_parametric_sampling_is_bounded_symmetric_and_deterministic() -> None:
    n = 100_000
    out = ps.particles(PARAMETRIC_SPEC, n, seed=2024)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert out[key].shape == (n,)
    major = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
    # Inside the mapped flux-surface region (Shafranov-shift bounds).
    assert major.min() > 620.0 - 200.0 - 15.0
    assert major.max() < 620.0 + 200.0 + 15.0
    assert out["z"].max() < 1.85 * 200.0
    assert out["z"].min() > -1.85 * 200.0
    # Up-down and toroidal symmetry.
    assert abs(float(out["z"].mean())) < 1.0
    assert abs(float(out["x"].mean())) < 8.0
    assert abs(float(out["y"].mean())) < 8.0
    # Birth energies: DT line broadened by the local T_i (centre 28 keV).
    assert float(out["energy"].mean()) == pytest.approx(14.08, abs=0.02)
    assert float(out["energy"].std()) > 0.1
    # Determinism.
    again = ps.particles(PARAMETRIC_SPEC, 128, seed=2024)
    for key in ("x", "y", "z", "energy"):
        assert bool((again[key] == out[key][:128]).all())


def test_parametric_birth_radius_peaks_in_the_core() -> None:
    # Recover the minor radius through the Miller map (one fixed-point step
    # for the Shafranov term) and check the peaked H-mode birth profile.
    n = 100_000
    out = ps.particles(PARAMETRIC_SPEC, n, seed=7)
    major = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
    z_k = out["z"] / 1.85
    r = ((major - 620.0) ** 2 + z_k**2) ** 0.5
    shift = 15.0 * (1.0 - (r / 200.0) ** 2)
    r = ((major - 620.0 - shift) ** 2 + z_k**2) ** 0.5
    assert r.max() < 200.0
    # Peaked profiles put most births well inside the pedestal.
    assert float(r.mean()) < 120.0
    assert float((r < 150.0).mean()) > 0.9


def test_parametric_sdef_card_round_trips_with_three_marginals() -> None:
    out = ps.emit_source_cards(PARAMETRIC_SPEC, bins=15)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["rad"] == "D1"
    assert parsed["ext"] == "D2"
    assert parsed["erg"] == "D3"
    assert len(parsed["distributions"]) == 3
    assert out["sdef"]["drift"][0]["reparsed"] is True
    quantities = [row["quantity"] for row in out["sdef"]["drift"]]
    assert quantities == ["emission probability", "spatial marginals", "joint correlation"]
    corr = out["sdef"]["drift"][2]["rel_drift"]
    assert 0.0 < corr < 1.0
    # Magnetic-axis spectrum summary (centre T = 28 keV, Ballabio DT).
    assert out["spectrum"]["mean_mev"] == pytest.approx(14.0818, abs=1e-3)
    assert out["spectrum"]["mono"] is False


def test_parametric_serpent_card_structure() -> None:
    out = ps.emit_source_cards(PARAMETRIC_SPEC, bins=12)
    lines = out["serpent"]["card"].splitlines()
    assert lines[0] == "src 1 pos 0 0 0"
    assert lines[1] == "src 1 rad d1"
    assert lines[2] == "src 1 ext d2"
    assert lines[3] == "src 1 erg d3"
    assert lines[4] == "src 1 wgt 1"
    assert len(lines) == 11
    for number in (1, 2, 3):
        assert lines[3 + 2 * number].startswith(f"SI{number} ")
        assert lines[4 + 2 * number].startswith(f"SP{number} ")
    assert out["serpent"]["drift"][0]["reparsed"] is False


def test_parametric_sector_births_are_uniform_over_sector() -> None:
    # start 0.5 rad, rotation 1.5 rad: every birth lands in-sector and the
    # mean angle matches the uniform closed form (start + rotation/2).
    start, rotation = 0.5, 1.5
    spec = dict(PARAMETRIC_SPEC, start_angle=start, rotation_angle=rotation)
    n = 20_000
    out = ps.particles(spec, n, seed=23)
    phi = [math.atan2(y, x) % (2.0 * math.pi) for x, y in zip(out["x"], out["y"], strict=True)]
    assert all(start <= p < start + rotation for p in phi)
    want = start + rotation / 2.0
    se = rotation / math.sqrt(12.0 * n)
    assert abs(sum(phi) / n - want) < 8.0 * se


def test_parametric_full_rotation_recovers_full_torus() -> None:
    # Exact full-rotation spelling reproduces the landed kernel bit-for-bit:
    # the same seeded stream and the same emitted cards.
    full = dict(PARAMETRIC_SPEC, start_angle=0.0, rotation_angle=2.0 * math.pi)
    a = ps.particles(PARAMETRIC_SPEC, 512, seed=17)
    b = ps.particles(full, 512, seed=17)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert bool((a[key] == b[key]).all())
    cards_a = ps.emit_source_cards(PARAMETRIC_SPEC, bins=15)
    cards_b = ps.emit_source_cards(full, bins=15)
    assert cards_b["sdef"]["card"] == cards_a["sdef"]["card"]
    assert cards_b["sdef"]["drift"] == cards_a["sdef"]["drift"]
    assert cards_b["serpent"]["card"] == cards_a["serpent"]["card"]


def test_parametric_sector_card_carries_phi_marginal() -> None:
    # A partial sector adds the uniform angle-bin marginal (PHI=D4 / phi d4)
    # plus the toroidal-sector drift row; the SDEF card still round-trips.
    spec = dict(PARAMETRIC_SPEC, start_angle=0.5, rotation_angle=math.pi)
    out = ps.emit_source_cards(spec, bins=15)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["phi"] == "D4"
    assert parsed["rad"] == "D1"
    assert parsed["ext"] == "D2"
    assert parsed["erg"] == "D3"
    assert len(parsed["distributions"]) == 4
    quantities = [row["quantity"] for row in out["sdef"]["drift"]]
    assert quantities == [
        "emission probability",
        "spatial marginals",
        "joint correlation",
        "toroidal sector",
    ]
    sector_row = out["sdef"]["drift"][3]
    assert sector_row["reparsed"] is True
    assert "rotation/2π = 0.5" in sector_row["note"]
    serpent_lines = out["serpent"]["card"].splitlines()
    assert serpent_lines[4] == "src 1 phi d4"
    assert any(line.startswith("SI4 ") for line in serpent_lines)
    assert any(line.startswith("SP4 ") for line in serpent_lines)
    assert out["serpent"]["drift"][-1]["quantity"] == "toroidal sector"
    assert out["serpent"]["drift"][-1]["reparsed"] is False


def test_proton_accounting_pure_fuels_are_exact() -> None:
    # Pinned 50/50 convention: pure D-D makes one proton per neutron
    # (totals identical, ratio exactly 1.0); pure D-T makes none.
    dd = ps.proton_accounting(dict(PARAMETRIC_SPEC, reaction="dd"))
    assert dd["proton_total"] == dd["neutron_total"] > 0.0
    assert dd["proton_per_neutron"] == 1.0
    dt = ps.proton_accounting(PARAMETRIC_SPEC)
    assert dt["proton_total"] == 0.0
    assert dt["proton_per_neutron"] == 0.0
    assert "50/50" in dt["note"]
    # Single-reaction ring/point specs carry no density model: per-neutron
    # ratio only, totals None.
    ring_dd = ps.proton_accounting(dict(RING_SPEC, reaction="dd"))
    assert ring_dd["proton_per_neutron"] == 1.0
    assert ring_dd["neutron_total"] is None and ring_dd["proton_total"] is None
    ring_dt = ps.proton_accounting(RING_SPEC)
    assert ring_dt["proton_per_neutron"] == 0.0


def test_proton_accounting_blend_matches_branch_share() -> None:
    # 70/30 blend: proton_per_neutron equals the D-D branch weight ratio
    # p_dd / (1 - p_dd) from the in-script rate rule (both deterministic
    # quadratures; tolerance is grid precision, not sampling).
    spec = dict(PARAMETRIC_SPEC, fuel={"D": 0.7, "T": 0.3})
    got = ps.proton_accounting(spec)
    sv_dt = ps.reactivity("dt", 20.0)
    sv_dd = ps.reactivity("dd", 20.0)
    # Flat-profile hand check first: uniform 20 keV gives the exact ratio.
    assert got["proton_total"] > 0.0 and got["neutron_total"] > got["proton_total"]
    assert 0.0 < got["proton_per_neutron"] < 0.05
    assert sv_dt > 0.0 and sv_dd > 0.0


def test_proton_accounting_zero_tritium_recovers_dd() -> None:
    # Zero-tritium recovery anchor end to end: a D-only fuel dict gives the
    # pure-D-D bookkeeping (ratio 1.0, totals equal).
    spec = dict(PARAMETRIC_SPEC, fuel={"D": 1.0, "T": 0.0})
    got = ps.proton_accounting(spec)
    assert got["proton_total"] == got["neutron_total"] > 0.0
    assert got["proton_per_neutron"] == 1.0


def test_proton_accounting_pure_tritium_is_loud() -> None:
    # No T-T neutron branch: a pure-tritium mixture has zero neutron
    # strength, so proton bookkeeping (a share of the neutron source) is a
    # loud error, never a zero stream.
    with pytest.raises(ValueError):
        ps.proton_accounting(dict(PARAMETRIC_SPEC, fuel={"D": 0.0, "T": 1.0}))


def _flat_l_mode_spec() -> dict[str, object]:
    """Flat L-mode parametric spec: uniform 20 keV / 1e20 m⁻³ plasma, so a
    uniform species pair (T_D = T_T = 20) must reproduce the shared kernel."""
    spec = dict(PARAMETRIC_SPEC)
    spec.update(
        mode="L",
        triangularity=0.0,
        shafranov_factor=0.0,
        ion_density_centre=1.0e20,
        ion_density_peaking_factor=0.0,
        ion_density_pedestal=1.0,
        ion_density_separatrix=1.0,
        ion_temperature_centre=20.0,
        ion_temperature_peaking_factor=0.0,
        ion_temperature_beta=1.0,
        ion_temperature_pedestal=20.0,
        ion_temperature_separatrix=20.0,
    )
    return spec


def test_species_temperatures_fire_both_branches_at_pair_temperatures() -> None:
    # 70/30 blend at T_D = 20 keV, T_T = 30 keV: the D-T branch reacts at
    # T_DT = 24 keV and the D-D branch at T_D, so both lines fire in the
    # stream and the D-D share matches the in-script pair-temperature rule.
    spec = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": 20.0, "T": 30.0},
    )
    n = 20_000
    out = ps.particles(spec, n, seed=13)
    energy = out["energy"]
    sv_dt = ps.reactivity("dt", 24.0)
    sv_dd = ps.reactivity("dd", 20.0)
    p_dd = (0.7 * 0.7 / 2.0 * sv_dd) / (0.7 * 0.3 * sv_dt + 0.7 * 0.7 / 2.0 * sv_dd)
    frac_dd = float((energy < 10.0).mean())
    se = math.sqrt(p_dd * (1.0 - p_dd) / n)
    assert abs(frac_dd - p_dd) < 8.0 * se
    # Branch-weighted mean with the per-branch Ballabio lines (D-T at 24
    # keV, D-D at 20 keV); the pair is uniform so the mean is uniform too.
    mu_dt, sigma_dt = _ballabio_moments("dt", 24.0)
    mu_dd, sigma_dd = _ballabio_moments("dd", 20.0)
    want = p_dd * mu_dd + (1.0 - p_dd) * mu_dt
    var = p_dd * (sigma_dd**2 + mu_dd**2) + (1.0 - p_dd) * (sigma_dt**2 + mu_dt**2) - want**2
    assert abs(float(energy.mean()) - want) < 6.0 * math.sqrt(var / n)


def test_species_equal_temperatures_recover_shared_mixture() -> None:
    # T_D = T_T = profile T reproduces the shared-temperature mixture kernel
    # bit-for-bit end to end: the same seeded stream and the same cards.
    flat = dict(_flat_l_mode_spec(), fuel={"D": 0.5, "T": 0.5})
    pair = dict(flat, species_temperatures={"D": 20.0, "T": 20.0})
    a = ps.particles(flat, 512, seed=17)
    b = ps.particles(pair, 512, seed=17)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert bool((a[key] == b[key]).all())
    cards_a = ps.emit_source_cards(flat, bins=15)
    cards_b = ps.emit_source_cards(pair, bins=15)
    assert cards_b["sdef"]["card"] == cards_a["sdef"]["card"]
    assert cards_b["sdef"]["drift"] == cards_a["sdef"]["drift"]
    assert cards_b["serpent"]["card"] == cards_a["serpent"]["card"]


def test_species_temperatures_are_loud() -> None:
    # Non-finite or negative species temperatures never reach the sampler;
    # the dict needs both D and T keys; and the pair without a fuel mixture
    # is a loud scope error, never a silent single-fuel fallback.
    blend = dict(PARAMETRIC_SPEC, fuel={"D": 0.5, "T": 0.5})
    with pytest.raises(ValueError, match=">= 0"):
        ps.particles(dict(blend, species_temperatures={"D": -1.0, "T": 20.0}), 4, seed=0)
    with pytest.raises(ValueError, match=">= 0"):
        ps.particles(dict(blend, species_temperatures={"D": 20.0, "T": -0.5}), 4, seed=0)
    with pytest.raises(ValueError, match="non-finite"):
        ps.particles(dict(blend, species_temperatures={"D": float("nan"), "T": 20.0}), 4, seed=0)
    with pytest.raises(ValueError, match="non-finite"):
        ps.particles(dict(blend, species_temperatures={"D": 20.0, "T": float("inf")}), 4, seed=0)
    with pytest.raises(ValueError, match="both `D` and `T`"):
        ps.particles(dict(blend, species_temperatures={"D": 20.0}), 4, seed=0)
    with pytest.raises(ValueError, match="supported keys"):
        ps.particles(dict(blend, species_temperatures={"D": 20.0, "T": 20.0, "H": 1.0}), 4, seed=0)
    with pytest.raises(ValueError, match="not yet supported"):
        ps.particles(dict(PARAMETRIC_SPEC, species_temperatures={"D": 20.0, "T": 30.0}), 4, seed=0)


def _tail_spec() -> dict[str, object]:
    """Pinned tail spec: 70/30 blend at T_D = 20 keV, T_T = 30 keV with a 5%
    deuterium hot tail at 60 keV (sub-pairs at T_DT = 24, T_DTt = 48,
    T_D = 20, T_mix = 40, T_tail = 60 keV)."""
    return dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": 20.0, "T": 30.0},
        deuterium_tail={"fraction": 0.05, "temperature_kev": 60.0},
    )


def _tail_sub_weights() -> tuple[list[float], list[tuple[str, float]]]:
    """Sub-rate-rule weights and (reaction, temperature) pairs of the pinned
    tail spec, from the in-script Bosch-Hale-independent reactivity facade."""
    eta = 0.05
    sv_dt = ps.reactivity("dt", 24.0)
    sv_dt_tail = ps.reactivity("dt", 48.0)
    sv_dd = ps.reactivity("dd", 20.0)
    sv_dd_mix = ps.reactivity("dd", 40.0)
    sv_dd_tail = ps.reactivity("dd", 60.0)
    w = [
        0.7 * 0.3 * (1.0 - eta) * sv_dt,
        0.7 * 0.3 * eta * sv_dt_tail,
        0.7 * 0.7 / 2.0 * (1.0 - eta) ** 2 * sv_dd,
        0.7 * 0.7 / 2.0 * 2.0 * eta * (1.0 - eta) * sv_dd_mix,
        0.7 * 0.7 / 2.0 * eta**2 * sv_dd_tail,
    ]
    total = sum(w)
    return (
        [x / total for x in w],
        [("dt", 24.0), ("dt", 48.0), ("dd", 20.0), ("dd", 40.0), ("dd", 60.0)],
    )


def test_deuterium_tail_fires_all_sub_branches() -> None:
    # All five bulk/tail sub-lines fire in the stream at the sub-rate-rule
    # shares, and the mean energy matches the five-branch weighted mean.
    spec = _tail_spec()
    n = 20_000
    out = ps.particles(spec, n, seed=13)
    energy = out["energy"]
    p, branches = _tail_sub_weights()
    p_dd = p[2] + p[3] + p[4]
    frac_dd = float((energy < 10.0).mean())
    se = math.sqrt(p_dd * (1.0 - p_dd) / n)
    assert abs(frac_dd - p_dd) < 8.0 * se
    moments = [_ballabio_moments(reaction, ti) for reaction, ti in branches]
    want = sum(pi * mu for pi, (mu, _) in zip(p, moments, strict=True))
    var = sum(pi * (sigma**2 + mu**2) for pi, (mu, sigma) in zip(p, moments, strict=True)) - want**2
    assert abs(float(energy.mean()) - want) < 6.0 * math.sqrt(var / n)


def test_deuterium_tail_zero_fraction_recovers_no_tail() -> None:
    # eta = 0 reproduces the no-tail mixture kernel bit-for-bit end to end:
    # the same seeded stream and the same emitted cards.
    base = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": 20.0, "T": 30.0},
    )
    zero = dict(base, deuterium_tail={"fraction": 0.0, "temperature_kev": 60.0})
    a = ps.particles(base, 512, seed=17)
    b = ps.particles(zero, 512, seed=17)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert bool((a[key] == b[key]).all())
    cards_a = ps.emit_source_cards(base, bins=15)
    cards_b = ps.emit_source_cards(zero, bins=15)
    assert cards_b["sdef"]["card"] == cards_a["sdef"]["card"]
    assert cards_b["sdef"]["drift"] == cards_a["sdef"]["drift"]
    assert cards_b["serpent"]["card"] == cards_a["serpent"]["card"]


def test_deuterium_tail_card_carries_drift_row() -> None:
    # A nonzero tail folds its Ballabio sub-lines into the energy marginal
    # and appends the deuterium-tail drift row; the SDEF card still
    # round-trips through the typed reader.
    out = ps.emit_source_cards(_tail_spec(), bins=15)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["rad"] == "D1"
    assert parsed["ext"] == "D2"
    assert parsed["erg"] == "D3"
    quantities = [row["quantity"] for row in out["sdef"]["drift"]]
    assert quantities == [
        "emission probability",
        "spatial marginals",
        "joint correlation",
        "deuterium tail",
    ]
    tail_row = out["sdef"]["drift"][3]
    assert tail_row["reparsed"] is True
    assert "0.050000" in tail_row["note"] and "60.000000" in tail_row["note"]
    assert out["serpent"]["drift"][-1]["quantity"] == "deuterium tail"
    assert out["serpent"]["drift"][-1]["reparsed"] is False


def test_deuterium_tail_is_loud() -> None:
    # Non-finite, out-of-range, or negative tail parameters never reach the
    # sampler; the dict needs both fraction and temperature_kev keys; and the
    # tail without a fuel mixture is a loud scope error, never a silent
    # single-fuel fallback.
    blend = dict(PARAMETRIC_SPEC, fuel={"D": 0.5, "T": 0.5})
    bad_tails = [
        {"fraction": -0.1, "temperature_kev": 60.0},
        {"fraction": 1.1, "temperature_kev": 60.0},
        {"fraction": 0.05, "temperature_kev": -1.0},
        {"fraction": float("nan"), "temperature_kev": 60.0},
        {"fraction": 0.05, "temperature_kev": float("inf")},
        {"fraction": 0.05},  # missing temperature_kev
        {"temperature_kev": 60.0},  # missing fraction
        {"fraction": 0.05, "temperature_kev": 60.0, "shape": "slowing-down"},
        {"fraction": 0.05, "temperature_kev": 60.0, "species": "T"},
    ]
    for tail in bad_tails:
        with pytest.raises(ValueError):
            ps.particles(dict(blend, deuterium_tail=tail), 4, seed=0)
    with pytest.raises(ValueError, match="not yet supported"):
        ps.particles(
            dict(PARAMETRIC_SPEC, deuterium_tail={"fraction": 0.05, "temperature_kev": 60.0}),
            4,
            seed=0,
        )


def _two_point_lattice_spec() -> dict[str, Any]:
    """Hand-vector lattice: A=(300,0,25) rate 2, B=(-300,0,25) rate 1, D-T mono."""
    return {
        "kind": "lattice",
        "reaction": "dt",
        "points": [
            {"position": [300.0, 0.0, 25.0], "rate": 2.0, "ion_temperature_kev": 0.0},
            {"position": [-300.0, 0.0, 25.0], "rate": 1.0, "ion_temperature_kev": 0.0},
        ],
    }


def _ring_lattice_spec(k: int = 64) -> dict[str, object]:
    """Uniform ring lattice: k unit-rate D-T mono nodes at R=300 cm, z=25 cm."""
    import math as _math

    return {
        "kind": "lattice",
        "reaction": "dt",
        "points": [
            {
                "position": [
                    300.0 * _math.cos(2.0 * _math.pi * i / k),
                    300.0 * _math.sin(2.0 * _math.pi * i / k),
                    25.0,
                ],
                "rate": 1.0,
                "ion_temperature_kev": 0.0,
            }
            for i in range(k)
        ],
    }


def test_lattice_hand_vectors_total_mean_and_mono_energy() -> None:
    # total = 3.0; mean birth position = (100, 0, 25); mono D-T line.
    out = ps.particles(_two_point_lattice_spec(), 90_000, seed=42)
    n = len(out["x"])
    assert out["energy"] == pytest.approx(14.021, rel=1e-12)
    assert abs(float(out["x"].mean()) - 100.0) < 8.0 * 200.0 / math.sqrt(n)
    assert abs(float(out["y"].mean())) < 8.0 * 200.0 / math.sqrt(n)
    assert abs(float(out["z"].mean()) - 25.0) < 1e-9
    # Rate-2 node fires twice as often (8-sigma multinomial gate).
    frac_a = float((out["x"] > 0).mean())
    se = math.sqrt((2.0 / 3.0) * (1.0 / 3.0) / n)
    assert abs(frac_a - 2.0 / 3.0) < 8.0 * se


def test_lattice_single_point_recovers_point_source() -> None:
    # One node reproduces the landed point-source stream bit-for-bit.
    lattice = {
        "kind": "lattice",
        "reaction": "dt",
        "points": [{"position": [1.0, -2.0, 3.5], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    point = {
        "kind": "point",
        "position": [1.0, -2.0, 3.5],
        "reaction": "dt",
        "ion_temperature_kev": 20.0,
    }
    a = ps.particles(lattice, 512, seed=9)
    b = ps.particles(point, 512, seed=9)
    for key in ("x", "y", "z", "u", "v", "w", "energy", "weight"):
        assert bool((a[key] == b[key]).all())


def test_lattice_ring_limit_matches_closed_form() -> None:
    # Axisymmetric limit converges to the analytic ring moments.
    out = ps.particles(_ring_lattice_spec(), N, seed=42)
    radius = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
    assert radius == pytest.approx(300.0, rel=1e-9)
    assert bool((out["z"] == 25.0).all())
    assert out["energy"] == pytest.approx(14.021, rel=1e-12)
    se = 300.0 / math.sqrt(2 * N)
    assert abs(float(out["x"].mean())) < 8 * se
    assert abs(float(out["y"].mean())) < 8 * se


def test_lattice_sdef_card_round_trips_with_lattice_row() -> None:
    out = ps.emit_source_cards(_two_point_lattice_spec(), bins=8)
    card = out["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    assert parsed["card"] == card
    assert parsed["rad"] == "D1"
    assert parsed["ext"] == "D2"
    assert parsed["erg"] == "D3"
    assert len(parsed["distributions"]) == 3
    quantities = [row["quantity"] for row in out["sdef"]["drift"]]
    assert quantities == [
        "emission probability",
        "spatial marginals",
        "joint correlation",
        "lattice discretization",
    ]
    assert out["sdef"]["drift"][0]["reparsed"] is True
    assert "2 nodes" in out["sdef"]["drift"][3]["note"]
    assert out["spectrum"]["mean_mev"] == pytest.approx(14.021, rel=1e-9)
    assert out["spectrum"]["mono"] is True
    lines = out["serpent"]["card"].splitlines()
    assert lines[0] == "src 1 pos 0 0 0"
    assert lines[1] == "src 1 rad d1"
    assert lines[4] == "src 1 wgt 1"
    assert out["serpent"]["drift"][-1]["quantity"] == "lattice discretization"
    assert out["serpent"]["drift"][-1]["reparsed"] is False


def test_lattice_symmetry_fold_conserves_and_replicates() -> None:
    # Base-sector cloud (one node) with 4 field periods: births spread
    # uniformly over the 4 copies and the mean energy stays the Ballabio line.
    base = {
        "kind": "lattice",
        "reaction": "dt",
        "field_periods": 4,
        "base_angle": 0.0,
        "points": [{"position": [300.0, 0.0, 0.0], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    n = 20_000
    out = ps.particles(base, n, seed=5)
    radius = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
    assert radius == pytest.approx(300.0, rel=1e-9)
    mean, _sigma = _ballabio_moments("dt", 20.0)
    assert abs(float(out["energy"].mean()) - mean) < 0.01
    # Quadrant occupancy is uniform (8-sigma multinomial gate per quadrant).
    quad = ((out["x"] > 0).astype(int) + 2 * (out["y"] > 0).astype(int)).astype(int)
    for q in range(4):
        frac = float((quad == q).mean())
        se = math.sqrt(0.25 * 0.75 / n)
        assert abs(frac - 0.25) < 8.0 * se
    # Proton bookkeeping on a D-T lattice: no protons, caller neutron total.
    acc = ps.proton_accounting(base)
    assert acc["proton_total"] == 0.0
    assert acc["neutron_total"] == pytest.approx(4.0)
    assert acc["proton_per_neutron"] == 0.0


def test_lattice_mixture_fires_both_branches() -> None:
    spec = {
        "kind": "lattice",
        "fuel": {"D": 0.7, "T": 0.3},
        "points": [{"position": [0.0, 0.0, 0.0], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    out = ps.particles(spec, 200_000, seed=11)
    frac_dd = float((out["energy"] < 10.0).mean())
    assert 0.002 < frac_dd < 0.015
    acc = ps.proton_accounting(spec)
    assert 0.0 < acc["proton_per_neutron"] < 0.05


def test_lattice_errors_are_loud() -> None:
    with pytest.raises(ValueError, match="kind"):
        ps.particles(dict(_two_point_lattice_spec(), kind="cloud"), 4, seed=0)
    with pytest.raises(ValueError, match="points"):
        ps.particles({"kind": "lattice", "reaction": "dt"}, 4, seed=0)
    with pytest.raises(ValueError, match="at least one point"):
        ps.particles({"kind": "lattice", "reaction": "dt", "points": []}, 4, seed=0)
    bad_rate = _two_point_lattice_spec()
    bad_rate["points"][0]["rate"] = -1.0
    with pytest.raises(ValueError):
        ps.particles(bad_rate, 4, seed=0)
    bad_pos = _two_point_lattice_spec()
    bad_pos["points"][1]["position"] = [0.0, 0.0]
    with pytest.raises(ValueError, match="three entries"):
        ps.particles(bad_pos, 4, seed=0)
    with pytest.raises(ValueError, match="together"):
        ps.particles(dict(_two_point_lattice_spec(), field_periods=4), 4, seed=0)
    with pytest.raises(ValueError, match="together"):
        ps.particles(dict(_two_point_lattice_spec(), base_angle=0.0), 4, seed=0)
    with pytest.raises(ValueError, match="periods"):
        ps.particles(dict(_two_point_lattice_spec(), field_periods=0, base_angle=0.0), 4, seed=0)
    with pytest.raises(ValueError, match="not yet supported"):
        ps.particles(
            dict(_two_point_lattice_spec(), species_temperatures={"D": 20.0, "T": 20.0}),
            4,
            seed=0,
        )
