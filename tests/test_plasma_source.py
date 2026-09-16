"""Python-side tests for the tokamak fusion-source core (analytic gates + errors)."""

import math

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
    # the loud boundary for toroidal sectors and for bad fuel fractions.
    ps.particles(dict(RING_SPEC, elongation=1.8), 4, seed=0)
    with pytest.raises(ValueError, match="not yet supported"):
        ps.particles(dict(PARAMETRIC_SPEC, rotation_angle=1.57), 4, seed=0)
    with pytest.raises(ValueError, match="not yet supported"):
        ps.particles(dict(PARAMETRIC_SPEC, start_angle=0.78), 4, seed=0)


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
