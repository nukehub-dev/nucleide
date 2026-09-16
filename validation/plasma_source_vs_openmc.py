"""Tokamak fusion-source cross-check (`nucleide.plasma_source` vs analytic gates +
openmc-plasma-source).

Two parts:

1. Analytic gates (always run): P1-P4 on the synthetic ring spec
   (closed-form ring/spectrum moments, sampler determinism, card emission),
   P5-P7 on a parametric Miller-geometry plasma (flat-profile Jacobian
   closed forms at zero triangularity, Bosch–Hale reactivity transcription,
   sampled birth moments vs independent fine quadrature of S·R·|J|), and P8
   on a pinned 70/30 D/T fuel blend (mixture-rule sampled moments and the
   D-D branch share vs in-script mixture quadrature).
2. openmc-plasma-source cross-check (container oracle): O1-O3 on the
   ring/point legs (Ballabio helpers, D-D energy moments, ring geometry),
   O4-O7 on the parametric leg (Miller map, L/H profiles, reactivity vs
   NeSST, sampled moments vs quadrature built from the upstream functions),
   and O8 on the pinned 70/30 blend (mixture moments vs quadrature over the
   upstream map/profiles with the NeSST D-T and D-D reactivities). The
   upstream package is an optional oracle dependency: if it cannot be
   imported, the oracle checks are reported as SKIP with their reason (never
   silently).

Units: Nucleide works in cm/MeV/keV; openmc-plasma-source in m/eV — the
oracle converts at the boundary.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from common import Report, fmt, rel_diff

import nucleide.mcnp as mcnp
import nucleide.plasma_source as ps

FAILURES = 0

# Tokamak ring: R = 3 m at z = 0.5 m, D-T, T_i = 20 keV (ITER-ish numbers,
# hand-picked round values).
RADIUS_CM = 300.0
HEIGHT_CM = 50.0
TI_KEV = 20.0
N = 200_000

RING_SPEC = {
    "kind": "ring",
    "radius": RADIUS_CM,
    "height": HEIGHT_CM,
    "reaction": "dt",
    "ion_temperature_kev": TI_KEV,
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


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _ballabio_moments(reaction: str, ti_kev: float) -> tuple[float, float]:
    """Independent in-script transcription of the Table III fits (keV in, MeV out)."""
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


def _reactivity_m3_s(reaction: str, ti_kev: float) -> float:
    """Independent in-script transcription of the Bosch–Hale/Atzeni fits."""
    if reaction == "dt":
        c1, a = 643.41e-22, 6.6610
        num = -1.0675e-4 * ti_kev**3 + 4.6064e-3 * ti_kev**2 + 1.5136e-2 * ti_kev
        den = 1.366e-5 * ti_kev**3 + 1.35e-2 * ti_kev**2 + 7.5189e-2 * ti_kev + 1.0
    else:
        c1, a = 3.5741e-22, 6.2696
        num = 5.8577e-3 * ti_kev
        den = -2.964e-6 * ti_kev**2 + 7.6822e-3 * ti_kev + 1.0
    xi = a * ti_kev ** (-1 / 3)
    eta = 1.0 - num / den
    return c1 * eta ** (-5 / 6) * xi**2 * math.exp(-3 * eta ** (1 / 3) * xi)


def _map_rz(spec: dict, a: float, theta: float) -> tuple[float, float]:
    """Miller forward map (Fausser spelling), independent transcription."""
    r0 = spec["major_radius"]
    esh = spec["shafranov_factor"]
    shift = esh * (1.0 - (a / spec["minor_radius"]) ** 2)
    big_r = r0 + a * math.cos(theta + spec["triangularity"] * math.sin(theta)) + shift
    z = spec["elongation"] * a * math.sin(theta)
    return big_r, z


def _volume_element(spec: dict, a: float, theta: float) -> float:
    """R(a,theta) * |J(a,theta)| from the analytic derivatives."""
    am = spec["minor_radius"]
    esh = spec["shafranov_factor"]
    delta = spec["triangularity"]
    kappa = spec["elongation"]
    th_m = theta + delta * math.sin(theta)
    dthm = 1.0 + delta * math.cos(theta)
    drda = math.cos(th_m) - 2.0 * esh * a / am**2
    drdt = -a * math.sin(th_m) * dthm
    dzda = kappa * math.sin(theta)
    dzdt = kappa * a * math.cos(theta)
    big_r, _ = _map_rz(spec, a, theta)
    return big_r * abs(drda * dzdt - drdt * dzda)


def _profile_density(spec: dict, r: float, mode: str | None = None) -> float:
    """Fausser density profile [m^-3], independent transcription."""
    am = spec["minor_radius"]
    rped = spec["pedestal_radius"]
    if (mode or spec["mode"]) == "L":
        return (
            spec["ion_density_centre"] * (1.0 - (r / am) ** 2) ** spec["ion_density_peaking_factor"]
        )
    if r < rped:
        core = (1.0 - (r / rped) ** 2) ** spec["ion_density_peaking_factor"]
        return (
            spec["ion_density_pedestal"]
            + (spec["ion_density_centre"] - spec["ion_density_pedestal"]) * core
        )
    return spec["ion_density_separatrix"] + (
        spec["ion_density_pedestal"] - spec["ion_density_separatrix"]
    ) * (am - r) / (am - rped)


def _profile_temperature(spec: dict, r: float, mode: str | None = None) -> float:
    """Fausser temperature profile [keV], independent transcription."""
    am = spec["minor_radius"]
    rped = spec["pedestal_radius"]
    if (mode or spec["mode"]) == "L":
        return (
            spec["ion_temperature_centre"]
            * (1.0 - (r / am) ** 2) ** spec["ion_temperature_peaking_factor"]
        )
    if r < rped:
        arg = 1.0 - (r / rped) ** spec["ion_temperature_beta"]
        core = arg ** spec["ion_temperature_peaking_factor"]
        return (
            spec["ion_temperature_pedestal"]
            + (spec["ion_temperature_centre"] - spec["ion_temperature_pedestal"]) * core
        )
    return spec["ion_temperature_separatrix"] + (
        spec["ion_temperature_pedestal"] - spec["ion_temperature_separatrix"]
    ) * (am - r) / (am - rped)


def _strength(spec: dict, r: float) -> float:
    """Relative neutron source density [1/s/cm^3 scale] at minor radius r."""
    n_cm3 = _profile_density(spec, r) * 1e-6
    sv = _reactivity_m3_s(spec["reaction"], _profile_temperature(spec, r)) * 1e6
    factor = 0.25 if spec["reaction"] == "dt" else 0.5
    return factor * n_cm3 * n_cm3 * sv


def _recover_minor_radius(spec: dict, major: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Invert the Miller map for the sampled minor radius: fixed point on the
    Shafranov term, converged (contraction factor ~ 2*esh*r/a^2 <= 0.075)."""
    r = ((major - spec["major_radius"]) ** 2 + (z / spec["elongation"]) ** 2) ** 0.5
    for _ in range(6):
        shift = spec["shafranov_factor"] * (1.0 - (r / spec["minor_radius"]) ** 2)
        r = ((major - spec["major_radius"] - shift) ** 2 + (z / spec["elongation"]) ** 2) ** 0.5
    return r


def _quadrature_moments(spec: dict, n_r: int = 400, n_t: int = 400) -> dict[str, float]:
    """Fine-quadrature reference for strength-weighted birth moments."""
    am = spec["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        s = _strength(spec, r)
        mu, _ = _ballabio_moments(spec["reaction"], _profile_temperature(spec, r))
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = s * _volume_element(spec, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * mu
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den}


def _quadrature_moments_mixture(
    spec: dict, f_d: float, f_t: float, n_r: int = 400, n_t: int = 400
) -> dict[str, float]:
    """Fine-quadrature reference under the pinned mixture rule
    S = n^2·[f_D·f_T·<sv>_DT + (f_D^2/2)·<sv>_DD]: strength-weighted <r^2>,
    branch-weighted mean birth energy, and the total D-D branch share."""
    am = spec["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = num_dd = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        n_cm3 = _profile_density(spec, r) * 1e-6
        ti = _profile_temperature(spec, r)
        sv_dt = _reactivity_m3_s("dt", ti) * 1e6
        sv_dd = _reactivity_m3_s("dd", ti) * 1e6
        w_dt = f_d * f_t * sv_dt
        w_dd = f_d * f_d / 2.0 * sv_dd
        mu_dt, _ = _ballabio_moments("dt", ti)
        mu_dd, _ = _ballabio_moments("dd", ti)
        p_dd = w_dd / (w_dt + w_dd) if w_dt + w_dd > 0.0 else 0.0
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = n_cm3 * n_cm3 * (w_dt + w_dd) * _volume_element(spec, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * (p_dd * mu_dd + (1.0 - p_dd) * mu_dt)
            num_dd += w * p_dd
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den, "p_dd": num_dd / den}


def analytic_gates() -> tuple[list[list[str]], list[str]]:
    """P1-P4 on the synthetic ring spec; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    out = ps.particles(RING_SPEC, N, seed=42)
    radius = np.hypot(out["x"], out["y"])
    notes.append(
        f"P1 ring sampling: n={N}, radius min/max {radius.min():.6f}/{radius.max():.6f} cm."
    )
    ok = bool(np.all(out["z"] == HEIGHT_CM)) and bool(
        np.all(np.abs(radius - RADIUS_CM) < 1e-9 * RADIUS_CM)
    )
    rows.append(["P1 ring radius/height", "-", "exact", _check(ok, "P1 radius/height")])
    se = RADIUS_CM / math.sqrt(2 * N)
    for label, value in (
        ("P1 <x>", float(np.mean(out["x"]))),
        ("P1 <y>", float(np.mean(out["y"]))),
    ):
        err = abs(value) / se
        rows.append([label, fmt(value), "< 8 sigma", _check(err < 8.0, label)])
    mean_w = float(np.mean(out["w"]))
    rows.append(
        [
            "P1 direction isotropy <w>",
            fmt(mean_w),
            "< 1e-2",
            _check(abs(mean_w) < 1e-2, "P1 isotropy"),
        ]
    )

    for reaction in ("dd", "dt"):
        for ti in (1.0, 10.0, 20.0):
            got = ps.spectrum_moments(reaction, ti)
            mean, sigma = _ballabio_moments(reaction, ti)
            worst = max(rel_diff(got["mean_mev"], mean), rel_diff(got["sigma_mev"], sigma))
            rows.append(
                [
                    f"P2 {reaction.upper()} Ti={ti:g} keV moments",
                    fmt(worst),
                    "< 1e-12",
                    _check(worst < 1e-12, f"P2 {reaction} {ti}"),
                ]
            )
    notes.append("P2 spectrum moments vs in-script Ballabio 1998 Table III transcription.")

    a = ps.particles(RING_SPEC, 512, seed=1234)
    b = ps.particles(RING_SPEC, 512, seed=1234)
    c = ps.particles(RING_SPEC, 512, seed=1235)
    ok = all(bool(np.array_equal(a[k], b[k])) for k in a) and not np.array_equal(
        a["energy"], c["energy"]
    )
    rows.append(["P3 pinned-seed determinism", "-", "identical", _check(ok, "P3 determinism")])

    emitted = ps.emit_source_cards(RING_SPEC)
    card = emitted["sdef"]["card"]
    parsed = mcnp.parse_sdef(card)
    ok = parsed["card"] == card and parsed["rad"] == "D1" and parsed["erg"] == "D2"
    rows.append(
        ["P4 SDEF reader round trip", "-", "byte-identical", _check(ok, "P4 SDEF round trip")]
    )
    drift = emitted["sdef"]["drift"][0]
    tail = 1.0 - math.erf(4.0 / math.sqrt(2.0))
    rows.append(
        [
            "P4 tabulation tail drift",
            fmt(drift["rel_drift"]),
            "~ +/-4 sigma tail",
            _check(rel_diff(drift["rel_drift"], tail) < 1e-3, "P4 tail drift"),
        ]
    )
    serpent = emitted["serpent"]["card"].splitlines()
    ok = (
        serpent[1] == "src 1 rad d1"
        and serpent[4] == f"SI1 {RADIUS_CM:g} {RADIUS_CM:g}"
        and emitted["serpent"]["drift"][0]["reparsed"] is False
    )
    rows.append(["P4 Serpent src structure", "-", "shape", _check(ok, "P4 Serpent structure")])
    notes.append("P4 Serpent drift rows are analytic by design (no reader in the workspace).")
    return rows, notes


def parametric_gates() -> tuple[list[list[str]], list[str]]:
    """P5-P7 on the parametric plasma spec; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    # P5: flat-profile, zero-triangularity closed forms.
    flat = dict(
        PARAMETRIC_SPEC,
        mode="L",
        elongation=1.85,
        triangularity=0.0,
        shafranov_factor=0.0,
        ion_density_peaking_factor=0.0,
        ion_temperature_peaking_factor=0.0,
        ion_temperature_beta=1.0,
        ion_temperature_centre=20.0,
        ion_temperature_pedestal=20.0,
        ion_temperature_separatrix=20.0,
    )
    out = ps.particles(flat, N, seed=11)
    r0 = flat["major_radius"]
    a_minor = flat["minor_radius"]
    kappa = flat["elongation"]
    major = np.hypot(out["x"], out["y"])
    r = ((major - r0) ** 2 + (out["z"] / kappa) ** 2) ** 0.5
    r2 = float(np.mean(r**2))
    tol = 8.0 * a_minor / math.sqrt(2 * N) * a_minor
    rows.append(
        [
            "P5 flat <r^2> = a^2/2",
            fmt(rel_diff(r2, a_minor * a_minor / 2.0)),
            "< 5e-2",
            _check((r2 - a_minor * a_minor / 2.0) < tol, "P5 r^2"),
        ]
    )
    z_mean = float(np.mean(out["z"]))
    rows.append(
        [
            "P5 <z> = 0",
            fmt(z_mean),
            "< 8 sigma",
            _check(abs(z_mean) < 8.0 * a_minor / math.sqrt(2 * N), "P5 z"),
        ]
    )
    want_r = r0 + a_minor * a_minor / (4.0 * r0)
    err_r = rel_diff(float(np.mean(major)), want_r)
    rows.append(["P5 <R> = R0 + a^2/4R0", fmt(err_r), "< 5e-2", _check(err_r < 5e-2, "P5 R")])
    notes.append(
        "P5 flat-profile closed forms at delta=0: <r^2> = a^2/2, <z> = 0, and "
        "the outboard volume bias <R> = R0 + a^2/(4 R0) (volume element R|J|)."
    )

    # P6: reactivity transcription (Bosch-Hale/Atzeni).
    worst = 0.0
    for reaction in ("dt", "dd"):
        for ti in (1.0, 10.0, 20.0, 40.0):
            err = rel_diff(ps.reactivity(reaction, ti), _reactivity_m3_s(reaction, ti))
            worst = max(worst, err)
    rows.append(["P6 reactivity transcription", fmt(worst), "< 1e-12", _check(worst < 1e-12, "P6")])
    notes.append("P6 reactivity <sigma v> vs in-script Bosch-Hale/Atzeni transcription.")

    # P7: H-mode sampled moments vs independent fine quadrature.
    out = ps.particles(PARAMETRIC_SPEC, N, seed=5)
    quad = _quadrature_moments(PARAMETRIC_SPEC)
    e_mean = float(np.mean(out["energy"]))
    err_e = rel_diff(e_mean, quad["mean_e"])
    rows.append(["P7 mean birth energy", fmt(err_e), "< 1e-2", _check(err_e < 1e-2, "P7 E")])
    major = np.hypot(out["x"], out["y"])
    r = _recover_minor_radius(PARAMETRIC_SPEC, major, out["z"])
    err_r2 = rel_diff(float(np.mean(r**2)), quad["<r^2>"])
    rows.append(["P7 birth <r^2>", fmt(err_r2), "< 1e-2", _check(err_r2 < 1e-2, "P7 r2")])
    notes.append(
        f"P7 H-mode sampled moments vs fine quadrature of S·R·|J| "
        f"(in-script Fausser/Bosch-Hale transcription), n={N}."
    )

    # P8: pinned 70/30 D/T blend — sampled moments vs the in-script mixture
    # quadrature, plus the D-D branch share of the stream (the upstream
    # `fuel`-dict spelling; always run).
    f_d, f_t = 0.7, 0.3
    blend = dict(PARAMETRIC_SPEC, fuel={"D": f_d, "T": f_t})
    out = ps.particles(blend, N, seed=101)
    quad = _quadrature_moments_mixture(PARAMETRIC_SPEC, f_d, f_t)
    err_e = rel_diff(float(np.mean(out["energy"])), quad["mean_e"])
    rows.append(
        ["P8 mixture mean birth energy", fmt(err_e), "< 1e-2", _check(err_e < 1e-2, "P8 E")]
    )
    major = np.hypot(out["x"], out["y"])
    r = _recover_minor_radius(PARAMETRIC_SPEC, major, out["z"])
    err_r2 = rel_diff(float(np.mean(r**2)), quad["<r^2>"])
    rows.append(["P8 mixture birth <r^2>", fmt(err_r2), "< 1e-2", _check(err_r2 < 1e-2, "P8 r2")])
    frac_dd = float(np.mean(out["energy"] < 10.0))
    se_dd = math.sqrt(quad["p_dd"] * (1.0 - quad["p_dd"]) / N)
    rows.append(
        [
            "P8 D-D branch share",
            fmt(frac_dd - quad["p_dd"]),
            "< 8 sigma",
            _check(abs(frac_dd - quad["p_dd"]) < 8.0 * se_dd, "P8 share"),
        ]
    )
    notes.append(
        f"P8 70/30 D/T blend (fuel dict) vs in-script mixture quadrature "
        f"(S = n^2·[f_D·f_T·<sv>_DT + (f_D^2/2)·<sv>_DD], Eriksson/DRESS rate "
        f"rule), n={N}; the D-D branch share gate proves both Ballabio lines "
        f"fire in the sampled stream."
    )
    return rows, notes


def oracle_check_openmc_plasma_source() -> tuple[list[list[str]], list[str], bool]:
    """O1-O7 vs the upstream MIT-licensed package. Returns (rows, notes, skipped)."""
    # NeSST (< 1.2) still imports scipy.integrate.cumtrapz, removed in scipy
    # 1.14; the modern spelling is a drop-in for its usage. Shim the alias in
    # this oracle process only (never in the shipped library).
    try:
        import scipy.integrate
        from scipy.integrate import cumulative_trapezoid

        scipy.integrate.cumtrapz = cumulative_trapezoid  # type: ignore[attr-defined]
    except ImportError:
        pass
    try:
        from openmc_plasma_source import fusion_point_source, fusion_ring_source
        from openmc_plasma_source.fuel_types import (
            neutron_energy_mean,
            neutron_energy_std_dev,
        )
        from openmc_plasma_source.tokamak_source import (
            tokamak_convert_a_alpha_to_R_Z,
            tokamak_ion_density,
            tokamak_ion_temperature,
        )
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Oracle check (openmc-plasma-source) SKIPPED: {exc}"
        print(note)
        return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []

    # O1: coefficient transcription vs the upstream published-fit helpers.
    worst = 0.0
    for reaction, upstream in (("dd", "DD"), ("dt", "DT")):
        for ti_kev in (1.0, 10.0, 20.0):
            got = ps.spectrum_moments(reaction, ti_kev)
            mean_ev = neutron_energy_mean(ion_temperature=ti_kev * 1e3, reaction=upstream)
            std_ev = neutron_energy_std_dev(ion_temperature=ti_kev * 1e3, reaction=upstream)
            err = max(
                rel_diff(got["mean_mev"], mean_ev / 1e6),
                rel_diff(got["sigma_mev"], std_ev / 1e6),
            )
            worst = max(worst, err)
    rows.append(
        ["O1 Ballabio fit transcription", fmt(worst), "< 1e-12", _check(worst < 1e-12, "O1 fits")]
    )
    notes.append(
        "O1 compares the Rust coefficient transcription against the upstream "
        "package's Table III helpers at three ion temperatures."
    )

    # O2: D-D energy moments on sampled particles (upstream fuel={"D": 1.0} is
    # the pure D-D Gaussian; mixtures are out of scope here).
    ti_ev = TI_KEV * 1e3
    (upstream,) = fusion_ring_source(
        radius=RADIUS_CM / 100.0,
        z_placement=HEIGHT_CM / 100.0,
        temperature=ti_ev,
        fuel={"D": 1.0},
    )
    rng = np.random.default_rng(2026)
    try:
        sampled = upstream.energy.sample(N, seed=1)
        # OpenMC >= 0.13 returns (samples, weights); older versions return
        # the bare sample array.
        if isinstance(sampled, tuple):
            sampled = sampled[0]
        upstream_e = np.asarray(sampled) / 1e6  # eV -> MeV
    except (AttributeError, TypeError):
        mean_ev = neutron_energy_mean(ion_temperature=ti_ev, reaction="DD")
        std_ev = neutron_energy_std_dev(ion_temperature=ti_ev, reaction="DD")
        upstream_e = rng.normal(mean_ev, std_ev, N) / 1e6
    ours = ps.particles(dict(RING_SPEC, reaction="dd"), N, seed=42)
    err_mean = rel_diff(float(np.mean(ours["energy"])), float(np.mean(upstream_e)))
    err_std = rel_diff(float(np.std(ours["energy"])), float(np.std(upstream_e)))
    rows.append(["O2 DD energy mean", fmt(err_mean), "< 1e-2", _check(err_mean < 1e-2, "O2 mean")])
    rows.append(["O2 DD energy std", fmt(err_std), "< 1e-2", _check(err_std < 1e-2, "O2 std")])
    notes.append(
        f"O2 D-D sampled moments at T_i={TI_KEV:g} keV, n={N} per side "
        "(upstream space object sampled analytically where OpenMC exposes no sampler)."
    )

    # O3: ring geometry moments vs the upstream CylindricalIndependent space
    # (geometry is fuel-independent; positions compared in cm).
    (upstream_dt,) = fusion_ring_source(
        radius=RADIUS_CM / 100.0,
        z_placement=HEIGHT_CM / 100.0,
        temperature=ti_ev,
        fuel={"D": 0.5, "T": 0.5},
    )
    space = upstream_dt.space
    r_ref = float(np.asarray(space.r.x)[0]) * 100.0  # m -> cm
    z_ref = float(np.asarray(space.z.x)[0]) * 100.0
    phi_a, phi_b = float(space.phi.a), float(space.phi.b)
    phi = rng.uniform(phi_a, phi_b, N)
    up_x = r_ref * np.cos(phi)
    up_y = r_ref * np.sin(phi)
    err_r = max(
        rel_diff(
            float(np.mean(np.hypot(ours["x"], ours["y"]))), float(np.mean(np.hypot(up_x, up_y)))
        ),
        rel_diff(r_ref, RADIUS_CM),
    )
    rows.append(
        ["O3 ring radius moments", fmt(err_r), "< 1e-12", _check(err_r < 1e-12, "O3 radius")]
    )
    rows.append(
        [
            "O3 ring height moments",
            fmt(rel_diff(float(np.mean(ours["z"])), z_ref)),
            "< 1e-12",
            _check(rel_diff(float(np.mean(ours["z"])), z_ref) < 1e-12, "O3 height"),
        ]
    )
    # Second moment <x^2> = R^2/2 for a uniform azimuth; comparing the two
    # independent samplings catches a biased/shrunk azimuth distribution
    # (a first-moment comparison of two noisy ~0 means is meaningless).
    err_phi = rel_diff(float(np.mean(ours["x"] ** 2)), float(np.mean(up_x**2)))
    rows.append(["O3 azimuth <x^2>", fmt(err_phi), "< 1e-2", _check(err_phi < 1e-2, "O3 azimuth")])
    notes.append(
        "O3 geometry cross-check uses the upstream ring's own CylindricalIndependent "
        "parameters (r delta, z delta, uniform azimuth) against Nucleide's sampled ring."
    )

    # O3b: point source geometry.
    (up_point,) = fusion_point_source(
        coordinate=(1.0, -2.0, 0.5), temperature=ti_ev, fuel={"D": 1.0}
    )
    px, py, pz = (float(v) * 100.0 for v in up_point.space.xyz)
    ours_p = ps.particles(
        {
            "kind": "point",
            "position": [px, py, pz],
            "reaction": "dd",
            "ion_temperature_kev": TI_KEV,
        },
        128,
        seed=5,
    )
    ok = bool(np.all(ours_p["x"] == px)) and bool(np.all(ours_p["z"] == pz))
    rows.append(["O3 point position", "-", "exact", _check(ok, "O3 point")])

    # O4: Miller map transcription vs the upstream forward map (exact).
    g = PARAMETRIC_SPEC
    worst = 0.0
    for a_cm in (0.1, 50.0, 120.0, 199.9):
        for theta in (0.1, 1.1, 2.6, 4.2, 5.9):
            r_up, z_up = tokamak_convert_a_alpha_to_R_Z(
                a=a_cm,
                alpha=theta,
                shafranov_factor=g["shafranov_factor"],
                minor_radius=g["minor_radius"],
                major_radius=g["major_radius"],
                triangularity=g["triangularity"],
                elongation=g["elongation"],
            )
            r_ref, z_ref = _map_rz(g, a_cm, theta)
            worst = max(
                worst, rel_diff(float(r_ref), float(r_up)), rel_diff(float(z_ref), float(z_up))
            )
    rows.append(
        ["O4 Miller map transcription", fmt(worst), "< 1e-12", _check(worst < 1e-12, "O4 map")]
    )
    notes.append("O4 forward map R(a,θ), Z(a,θ) (Shafranov shift included) at a 4×5 grid.")

    # O5: profile transcription vs the upstream L/H formulas (exact).
    worst_n = worst_t = 0.0
    for mode in ("L", "H"):
        for r_cm in (0.0, 30.0, 100.0, 149.9, 150.1, 180.0, 200.0):
            n_ref = _profile_density(g, r_cm, mode=mode)
            t_ref = _profile_temperature(g, r_cm, mode=mode)
            n_up = float(
                tokamak_ion_density(
                    mode=mode,
                    ion_density_centre=g["ion_density_centre"],
                    ion_density_peaking_factor=g["ion_density_peaking_factor"],
                    ion_density_pedestal=g["ion_density_pedestal"],
                    minor_radius=g["minor_radius"],
                    pedestal_radius=g["pedestal_radius"],
                    ion_density_separatrix=g["ion_density_separatrix"],
                    r=r_cm,
                )
            )
            t_up = (
                float(
                    tokamak_ion_temperature(
                        r=r_cm,
                        mode=mode,
                        pedestal_radius=g["pedestal_radius"],
                        ion_temperature_pedestal=g["ion_temperature_pedestal"],
                        ion_temperature_centre=g["ion_temperature_centre"],
                        ion_temperature_beta=g["ion_temperature_beta"],
                        ion_temperature_peaking_factor=g["ion_temperature_peaking_factor"],
                        ion_temperature_separatrix=g["ion_temperature_separatrix"],
                        minor_radius=g["minor_radius"],
                    )
                )
                / 1e3
            )  # eV -> keV
            worst_n = max(worst_n, rel_diff(n_ref, n_up))
            worst_t = max(worst_t, rel_diff(t_ref, t_up))
    rows.append(["O5 density profile", fmt(worst_n), "< 1e-12", _check(worst_n < 1e-12, "O5 n")])
    rows.append(
        ["O5 temperature profile", fmt(worst_t), "< 1e-12", _check(worst_t < 1e-12, "O5 T")]
    )
    notes.append(
        "O5 Fausser L/H density and temperature profiles at seven radii (keV/eV converted)."
    )

    # O6: reactivity vs NeSST (published-coefficient precision ~5e-9).
    from NeSST.spectral_model import reac_DD, reac_DT

    worst = 0.0
    for reaction, up in (("dt", reac_DT), ("dd", reac_DD)):
        for ti_kev in (1.0, 10.0, 20.0, 40.0):
            err = rel_diff(ps.reactivity(reaction, ti_kev), float(up(ti_kev * 1e3)))
            worst = max(worst, err)
    rows.append(
        ["O6 reactivity vs NeSST", fmt(worst), "< 1e-6", _check(worst < 1e-6, "O6 reactivity")]
    )
    notes.append(
        "O6 Bosch-Hale reactivity vs NeSST (tolerance 1e-6 reflects the "
        "published coefficients' precision, not the implementation)."
    )

    # O7: sampled parametric moments vs quadrature built from the upstream
    # map/profile functions with the NeSST D-T reactivity (their public
    # tokamak_source only offers fuel mixtures, so the end-to-end reference
    # integrates the same published formulas over their own map).
    ours = ps.particles(PARAMETRIC_SPEC, N, seed=77)
    quad = _quadrature_upstream(g, tokamak_ion_density, tokamak_ion_temperature, reac_DT)
    err_e = rel_diff(float(np.mean(ours["energy"])), quad["mean_e"])
    rows.append(["O7 mean birth energy", fmt(err_e), "< 1e-2", _check(err_e < 1e-2, "O7 E")])
    major = np.hypot(ours["x"], ours["y"])
    r = _recover_minor_radius(g, major, ours["z"])
    err_r2 = rel_diff(float(np.mean(r**2)), quad["<r^2>"])
    rows.append(["O7 birth <r^2>", fmt(err_r2), "< 1e-2", _check(err_r2 < 1e-2, "O7 r2")])
    notes.append(
        f"O7 end-to-end parametric moments vs quadrature over the upstream "
        f"map/profiles with the NeSST D-T reactivity, n={N}."
    )

    # O8: pinned 70/30 D/T blend vs quadrature over the upstream
    # map/profiles with the NeSST D-T and D-D reactivities under the pinned
    # mixture rule. Upstream's tokamak_source additionally models T-T for a
    # D+T blend; Nucleide has no T-T branch, so the reference integrates the
    # two shared branches only (the O7 precedent).
    f_d, f_t = 0.7, 0.3
    blend = dict(PARAMETRIC_SPEC, fuel={"D": f_d, "T": f_t})
    ours = ps.particles(blend, N, seed=78)
    quad = _quadrature_upstream_mixture(
        g, tokamak_ion_density, tokamak_ion_temperature, reac_DT, reac_DD, f_d, f_t
    )
    err_e = rel_diff(float(np.mean(ours["energy"])), quad["mean_e"])
    rows.append(
        ["O8 mixture mean birth energy", fmt(err_e), "< 1e-2", _check(err_e < 1e-2, "O8 E")]
    )
    major = np.hypot(ours["x"], ours["y"])
    r = _recover_minor_radius(g, major, ours["z"])
    err_r2 = rel_diff(float(np.mean(r**2)), quad["<r^2>"])
    rows.append(["O8 mixture birth <r^2>", fmt(err_r2), "< 1e-2", _check(err_r2 < 1e-2, "O8 r2")])
    notes.append(
        f"O8 70/30 D/T blend vs quadrature over the upstream map/profiles "
        f"with NeSST reac_DT + reac_DD under the mixture rule, n={N} "
        f"(upstream T-T branch excluded — no T-T reaction in Nucleide)."
    )
    return rows, notes, False


def _quadrature_upstream_mixture(
    g: dict,
    density_fn,
    temperature_fn,
    reac_dt,
    reac_dd,
    f_d: float,
    f_t: float,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Quadrature of n^2·[f_D·f_T·reac_DT + (f_D^2/2)·reac_DD]·R·|J| using the
    upstream profile functions and NeSST reactivities; the volume element
    uses the in-script Miller transcription (its agreement with the upstream
    map is gate O4)."""
    am = g["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        n_m3 = float(
            density_fn(
                mode=g["mode"],
                ion_density_centre=g["ion_density_centre"],
                ion_density_peaking_factor=g["ion_density_peaking_factor"],
                ion_density_pedestal=g["ion_density_pedestal"],
                minor_radius=am,
                pedestal_radius=g["pedestal_radius"],
                ion_density_separatrix=g["ion_density_separatrix"],
                r=r,
            )
        )
        t_kev = (
            float(
                temperature_fn(
                    r=r,
                    mode=g["mode"],
                    pedestal_radius=g["pedestal_radius"],
                    ion_temperature_pedestal=g["ion_temperature_pedestal"],
                    ion_temperature_centre=g["ion_temperature_centre"],
                    ion_temperature_beta=g["ion_temperature_beta"],
                    ion_temperature_peaking_factor=g["ion_temperature_peaking_factor"],
                    ion_temperature_separatrix=g["ion_temperature_separatrix"],
                    minor_radius=am,
                )
            )
            / 1e3
        )
        sv_dt = float(reac_dt(t_kev * 1e3)) * 1e6  # m^3/s -> cm^3/s
        sv_dd = float(reac_dd(t_kev * 1e3)) * 1e6
        w_dt = f_d * f_t * sv_dt
        w_dd = f_d * f_d / 2.0 * sv_dd
        mu_dt, _ = _ballabio_moments("dt", t_kev)
        mu_dd, _ = _ballabio_moments("dd", t_kev)
        p_dd = w_dd / (w_dt + w_dd) if w_dt + w_dd > 0.0 else 0.0
        strength = n_m3 * n_m3 * 1e-12 * (w_dt + w_dd)  # (n·1e-6)^2 = n^2·1e-12
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = strength * _volume_element(g, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * (p_dd * mu_dd + (1.0 - p_dd) * mu_dt)
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den}


def _quadrature_upstream(
    g: dict,
    density_fn,
    temperature_fn,
    reac_dt,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Quadrature of n^2/4 * <sigma v>_DT * R * |J| using the upstream profile
    functions and NeSST reactivity; the volume element uses the in-script
    Miller transcription (its agreement with the upstream map is gate O4)."""
    am = g["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        n_m3 = float(
            density_fn(
                mode=g["mode"],
                ion_density_centre=g["ion_density_centre"],
                ion_density_peaking_factor=g["ion_density_peaking_factor"],
                ion_density_pedestal=g["ion_density_pedestal"],
                minor_radius=am,
                pedestal_radius=g["pedestal_radius"],
                ion_density_separatrix=g["ion_density_separatrix"],
                r=r,
            )
        )
        t_kev = (
            float(
                temperature_fn(
                    r=r,
                    mode=g["mode"],
                    pedestal_radius=g["pedestal_radius"],
                    ion_temperature_pedestal=g["ion_temperature_pedestal"],
                    ion_temperature_centre=g["ion_temperature_centre"],
                    ion_temperature_beta=g["ion_temperature_beta"],
                    ion_temperature_peaking_factor=g["ion_temperature_peaking_factor"],
                    ion_temperature_separatrix=g["ion_temperature_separatrix"],
                    minor_radius=am,
                )
            )
            / 1e3
        )
        sv = float(reac_dt(t_kev * 1e3)) * 1e6  # m^3/s -> cm^3/s
        strength = 0.25 * (n_m3 * 1e-6) ** 2 * sv
        mu, _ = _ballabio_moments("dt", t_kev)
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = strength * _volume_element(g, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * mu
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den}


def main() -> int:
    report = Report("plasma_source", "Tokamak fusion sources (`plasma_source_vs_openmc.py`)")
    report.prose(
        "Two-part oracle for `nucleide.plasma_source`: analytic gates P1-P4 on "
        "the synthetic ring spec (closed-form ring/spectrum moments, sampler "
        "determinism, card emission), P5-P7 on a parametric Miller-geometry "
        "plasma (flat-profile Jacobian closed forms, Bosch-Hale reactivity "
        "transcription, sampled moments vs fine quadrature), and P8 on a "
        "pinned 70/30 D/T fuel blend (mixture-rule moments and the D-D "
        "branch share vs in-script mixture quadrature) — always run — plus "
        "container-only cross-checks O1-O8 against the upstream "
        "openmc-plasma-source package and NeSST (Ballabio helpers, sampled "
        "ring/point sources, Miller map, Fausser profiles, reactivities, "
        "end-to-end parametric moments, and the pinned 70/30 blend)."
    )
    rows1, notes1 = analytic_gates()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows1)
    rows2, notes2 = parametric_gates()
    for note in notes2:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows2)
    rows3, notes3, skipped = oracle_check_openmc_plasma_source()
    for note in notes3:
        report.prose(note)
    if skipped:
        skip_rows = [
            [f"{gate}", "SKIP (openmc-plasma-source unavailable)"]
            for gate in ("O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8")
        ]
        report.table(["Gate", "Status"], skip_rows)
    else:
        report.table(["Gate", "Rel err", "Tol", "Status"], rows3)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} plasma-source check(s) failed", file=sys.stderr)
        return 1
    print("All plasma-source checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
