"""Tokamak fusion-source cross-check (`nucleide.plasma_source` vs analytic gates +
openmc-plasma-source).

Two parts:

1. Analytic gates (always run): P1-P4 on the synthetic ring spec
   (closed-form ring/spectrum moments, sampler determinism, card emission),
   P5-P7 on a parametric Miller-geometry plasma (flat-profile Jacobian
   closed forms at zero triangularity, Bosch–Hale reactivity transcription,
   sampled birth moments vs independent fine quadrature of S·R·|J|), and P8
   on a pinned 70/30 D/T fuel blend (mixture-rule sampled moments and the
   D-D branch share vs in-script mixture quadrature), plus P9 on the
   D(d,p)T proton bookkeeping (pure-fuel exact ratios, the 70/30 share vs
   quadrature, the no-T recovery anchor, and the loud pure-T error), plus
   P10 on toroidal sectors (in-sector uniform births, bit-for-bit
   full-rotation recovery, the PHI=D4 card marginal with its drift row),
   plus P11 on per-species ion temperatures (pair-temperature sampled
   moments and branch share vs quadrature, bit-for-bit equal-T recovery),
   plus P12 on the arbitrary-3D birth-rate lattice (hand-vector totals and
   means, single-point bit-for-bit recovery of the point source, ring-limit
   moment closure, symmetry-fold replication and totals, card round trip
   with the lattice drift row, mixture branch fire, loud errors),
   plus P13 on the deuterium hot-tail fraction (sub-rate-rule moments and
   branch share vs quadrature, the card drift tail-neutron share, bit-for-bit
   zero-fraction recovery, loud errors).
2. openmc-plasma-source cross-check (container oracle): O1-O3 on the
   ring/point legs (Ballabio helpers, D-D energy moments, ring geometry),
   O4-O7 on the parametric leg (Miller map, L/H profiles, reactivity vs
   NeSST, sampled moments vs quadrature built from the upstream functions),
   O8 on the pinned 70/30 blend (mixture moments vs quadrature over the
   upstream map/profiles with the NeSST D-T and D-D reactivities), and O9
   on a tritium-rich 10/90 blend (two-branch kernel vs the three-branch
   upstream model: the deferred T-T branch quantified, never ignored).
   Protons have no oracle — the upstream package models neutrons only —
   so D(d,p)T accounting rests on the always-run P9 gates — and O10 on a
   partial toroidal sector ((r, z, E) moments vs the upstream quadrature
   plus the uniform angle marginal) — and O11 on a 70/30 blend at distinct
   species temperatures (moments vs the upstream quadrature with the NeSST
   reactivities evaluated at the pair temperatures) — and O12 on the lattice
    ring limit (dense monoenergetic ring cloud vs the upstream ring geometry
    and Ballabio helpers, where those models apply) — and O13 on the
    deuterium hot tail (pinned tail blend vs the upstream quadrature with
    the NeSST reactivities evaluated at the five bulk/tail sub-pair
    temperatures). The
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


def _quadrature_moments_mixture_species(
    spec: dict, f_d: float, f_t: float, t_d: float, t_t: float, n_r: int = 400, n_t: int = 400
) -> dict[str, float]:
    """Fine-quadrature reference under the per-species-temperature mixture rule
    S = n^2·[f_D·f_T·<sv>_DT(T_DT) + (f_D^2/2)·<sv>_DD(T_D)] with the
    mass-weighted T_DT = T_D + (2/5)·(T_T − T_D): strength-weighted <r^2>,
    branch-weighted mean birth energy, and the total D-D branch share."""
    t_dt = t_d + 0.4 * (t_t - t_d)
    am = spec["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = num_dd = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        n_cm3 = _profile_density(spec, r) * 1e-6
        sv_dt = _reactivity_m3_s("dt", t_dt) * 1e6
        sv_dd = _reactivity_m3_s("dd", t_d) * 1e6
        w_dt = f_d * f_t * sv_dt
        w_dd = f_d * f_d / 2.0 * sv_dd
        mu_dt, _ = _ballabio_moments("dt", t_dt)
        mu_dd, _ = _ballabio_moments("dd", t_d)
        p_dd = w_dd / (w_dt + w_dd) if w_dt + w_dd > 0.0 else 0.0
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = n_cm3 * n_cm3 * (w_dt + w_dd) * _volume_element(spec, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * (p_dd * mu_dd + (1.0 - p_dd) * mu_dt)
            num_dd += w * p_dd
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den, "p_dd": num_dd / den}


def _quadrature_moments_mixture_tail(
    spec: dict,
    f_d: float,
    f_t: float,
    t_d: float,
    t_t: float,
    eta: float,
    t_tail: float,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Fine-quadrature reference under the deuterium hot-tail sub-rate rule
    (module-rustdoc mixture v3): the bulk/tail sub-pairs react at T_DT, T_DTt
    (D-T) and T_D, T_mix, T_tail (D-D) with the sub-rate weights. Returns the
    strength-weighted <r^2>, the sub-branch-weighted mean birth energy, the
    total D-D branch share, and the tail neutron share (the eta-scaled
    sub-rate over the total)."""
    t_dt = t_d + 0.4 * (t_t - t_d)
    t_dtt = t_tail + 0.4 * (t_t - t_tail)
    t_mix = 0.5 * (t_d + t_tail)
    am = spec["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = num_dd = num_tail = 0.0
    for i in range(n_r):
        r = (i + 0.5) * dr
        n_cm3 = _profile_density(spec, r) * 1e-6
        sv_dt = _reactivity_m3_s("dt", t_dt) * 1e6
        sv_dt_tail = _reactivity_m3_s("dt", t_dtt) * 1e6
        sv_dd = _reactivity_m3_s("dd", t_d) * 1e6
        sv_dd_mix = _reactivity_m3_s("dd", t_mix) * 1e6
        sv_dd_tail = _reactivity_m3_s("dd", t_tail) * 1e6
        w = [
            f_d * f_t * (1.0 - eta) * sv_dt,
            f_d * f_t * eta * sv_dt_tail,
            f_d * f_d / 2.0 * (1.0 - eta) ** 2 * sv_dd,
            f_d * f_d / 2.0 * 2.0 * eta * (1.0 - eta) * sv_dd_mix,
            f_d * f_d / 2.0 * eta**2 * sv_dd_tail,
        ]
        mus = [
            _ballabio_moments("dt", t_dt)[0],
            _ballabio_moments("dt", t_dtt)[0],
            _ballabio_moments("dd", t_d)[0],
            _ballabio_moments("dd", t_mix)[0],
            _ballabio_moments("dd", t_tail)[0],
        ]
        total = sum(w)
        p = [x / total for x in w] if total > 0.0 else [1.0, 0.0, 0.0, 0.0, 0.0]
        mean_e = sum(pi * mu for pi, mu in zip(p, mus, strict=True))
        p_dd = p[2] + p[3] + p[4]
        p_tail = (w[1] + w[3] + w[4]) / total if total > 0.0 else 0.0
        for j in range(n_t):
            theta = (j + 0.5) * dt
            cell = n_cm3 * n_cm3 * total * _volume_element(spec, r, theta) * dr * dt
            den += cell
            num_r2 += cell * r * r
            num_e += cell * mean_e
            num_dd += cell * p_dd
            num_tail += cell * p_tail
    return {
        "<r^2>": num_r2 / den,
        "mean_e": num_e / den,
        "p_dd": num_dd / den,
        "p_tail": num_tail / den,
    }


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

    # P9: D(d,p)T proton accounting + zero-tritium recovery (always run).
    # Proton bookkeeping is analytic (no sampling noise): pure D-D gives one
    # proton per neutron exactly, D-T-only gives none, the 70/30 ratio
    # matches the in-script mixture quadrature share, a D-only fuel dict
    # recovers the pure-D-D bookkeeping exactly, and pure tritium (no T-T
    # neutron branch) is a loud error, never a zero stream.
    dd_acc = ps.proton_accounting(dict(PARAMETRIC_SPEC, reaction="dd"))
    rows.append(
        [
            "P9 pure D-D proton/neutron",
            fmt(dd_acc["proton_per_neutron"] - 1.0),
            "exact",
            _check(
                dd_acc["proton_per_neutron"] == 1.0
                and dd_acc["proton_total"] == dd_acc["neutron_total"],
                "P9 pure DD",
            ),
        ]
    )
    dt_acc = ps.proton_accounting(PARAMETRIC_SPEC)
    rows.append(
        [
            "P9 D-T-only protons",
            fmt(dt_acc["proton_total"]),
            "exact 0",
            _check(
                dt_acc["proton_total"] == 0.0 and dt_acc["proton_per_neutron"] == 0.0,
                "P9 DT",
            ),
        ]
    )
    got_acc = ps.proton_accounting(blend)
    want_share = quad["p_dd"] / (1.0 - quad["p_dd"])
    err_p = rel_diff(got_acc["proton_per_neutron"], want_share)
    rows.append(["P9 70/30 proton share", fmt(err_p), "< 1e-2", _check(err_p < 1e-2, "P9 share")])
    d_only = ps.proton_accounting(dict(PARAMETRIC_SPEC, fuel={"D": 1.0, "T": 0.0}))
    rows.append(
        [
            "P9 no-T recovery",
            fmt(d_only["proton_per_neutron"] - 1.0),
            "exact",
            _check(
                d_only["proton_per_neutron"] == 1.0
                and d_only["proton_total"] == d_only["neutron_total"],
                "P9 no-T",
            ),
        ]
    )
    try:
        ps.proton_accounting(dict(PARAMETRIC_SPEC, fuel={"D": 0.0, "T": 1.0}))
        loud_ok = False
    except Exception:  # noqa: BLE001 — any loud error satisfies the gate
        loud_ok = True
    rows.append(["P9 pure-T loud", "-", "raises", _check(loud_ok, "P9 pure-T")])
    notes.append(
        "P9 D(d,p)T proton bookkeeping (pinned 50/50 with the D-D neutron "
        "branch; sampler and cards stay neutron-only): exact pure-fuel "
        "ratios, the 70/30 share vs the P8 quadrature (tolerance is the "
        "256-vs-400 grid precision on the H-mode pedestal kink — exactness "
        "is pinned by the crate hand vectors), the no-T recovery anchor, "
        "and the loud pure-T error."
    )

    # P10: toroidal sectors (always run). Births are uniform over
    # [start, start + rotation) (mean-angle gate at 8 sigma, every birth
    # in-sector), the exact full-rotation spelling reproduces the landed
    # kernel bit-for-bit (same seeded stream, same emitted cards), and a
    # partial sector adds the PHI=D4 angle-bin marginal plus the
    # toroidal-sector drift row while the SDEF card still round-trips.
    start, rotation = 0.5, 1.5
    sector = dict(PARAMETRIC_SPEC, start_angle=start, rotation_angle=rotation)
    out = ps.particles(sector, N, seed=23)
    phi = np.arctan2(out["y"], out["x"]) % (2.0 * math.pi)
    in_sector = bool(((phi >= start) & (phi < start + rotation)).all())
    rows.append(["P10 sector births in-sector", "-", "all", _check(in_sector, "P10 in-sector")])
    want_phi = start + rotation / 2.0
    se_phi = rotation / math.sqrt(12.0 * N)
    mean_phi = float(np.mean(phi))
    rows.append(
        [
            "P10 sector mean angle",
            fmt(mean_phi - want_phi),
            "< 8 sigma",
            _check(abs(mean_phi - want_phi) < 8.0 * se_phi, "P10 angle"),
        ]
    )
    full = dict(PARAMETRIC_SPEC, start_angle=0.0, rotation_angle=2.0 * math.pi)
    a = ps.particles(PARAMETRIC_SPEC, 512, seed=17)
    b = ps.particles(full, 512, seed=17)
    recovered = all(bool((a[key] == b[key]).all()) for key in ("x", "y", "z", "energy"))
    rows.append(["P10 full-rotation stream", "-", "exact", _check(recovered, "P10 stream")])
    cards_a = ps.emit_source_cards(PARAMETRIC_SPEC, bins=15)
    cards_b = ps.emit_source_cards(full, bins=15)
    same_cards = (
        cards_b["sdef"]["card"] == cards_a["sdef"]["card"]
        and cards_b["serpent"]["card"] == cards_a["serpent"]["card"]
    )
    rows.append(["P10 full-rotation cards", "-", "exact", _check(same_cards, "P10 cards")])
    sec_cards = ps.emit_source_cards(sector, bins=15)
    sec_parsed = mcnp.parse_sdef(sec_cards["sdef"]["card"])
    phi_ok = (
        sec_parsed["card"] == sec_cards["sdef"]["card"]
        and sec_parsed["phi"] == "D4"
        and len(sec_parsed["distributions"]) == 4
    )
    rows.append(["P10 sector PHI marginal", "-", "round-trip", _check(phi_ok, "P10 PHI")])
    drift_quantities = [row["quantity"] for row in sec_cards["sdef"]["drift"]]
    drift_ok = drift_quantities[-1] == "toroidal sector" and bool(
        sec_cards["sdef"]["drift"][-1]["reparsed"]
    )
    rows.append(["P10 sector drift row", "-", "present", _check(drift_ok, "P10 drift")])
    notes.append(
        f"P10 toroidal sector (start {start} rad, rotation {rotation} rad): "
        f"uniform birth angles over the sector (n={N}), bit-for-bit "
        f"full-rotation recovery of the landed kernel (stream and cards), "
        f"and the PHI=D4 angle-bin marginal with its drift row."
    )

    # P11: per-species ion temperatures (always run). A 70/30 blend at
    # T_D = 20 keV, T_T = 30 keV reacts the D-T branch at T_DT = 24 keV and
    # the D-D branch at T_D: sampled mean energy, radius moment, and D-D
    # branch share sit on the pair-temperature quadrature. On a flat
    # uniform-20 keV L-mode plasma the equal pair (20, 20) reproduces the
    # shared-temperature mixture kernel bit-for-bit (stream and cards).
    t_d11, t_t11 = 20.0, 30.0
    species11 = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": t_d11, "T": t_t11},
    )
    ours11 = ps.particles(species11, N, seed=81)
    quad11 = _quadrature_moments_mixture_species(PARAMETRIC_SPEC, 0.7, 0.3, t_d11, t_t11)
    err_e11 = rel_diff(float(np.mean(ours11["energy"])), quad11["mean_e"])
    rows.append(
        ["P11 species mean birth energy", fmt(err_e11), "< 1e-2", _check(err_e11 < 1e-2, "P11 E")]
    )
    major11 = np.hypot(ours11["x"], ours11["y"])
    r11 = _recover_minor_radius(PARAMETRIC_SPEC, major11, ours11["z"])
    err_r11 = rel_diff(float(np.mean(r11**2)), quad11["<r^2>"])
    rows.append(
        ["P11 species birth <r^2>", fmt(err_r11), "< 1e-2", _check(err_r11 < 1e-2, "P11 r2")]
    )
    frac_dd11 = float(np.mean(ours11["energy"] < 10.0))
    se_dd11 = math.sqrt(quad11["p_dd"] * (1.0 - quad11["p_dd"]) / N)
    rows.append(
        [
            "P11 species D-D branch share",
            fmt(frac_dd11 - quad11["p_dd"]),
            "< 8 sigma",
            _check(abs(frac_dd11 - quad11["p_dd"]) < 8.0 * se_dd11, "P11 share"),
        ]
    )
    flat11 = dict(
        PARAMETRIC_SPEC,
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
        fuel={"D": 0.5, "T": 0.5},
    )
    pair11 = dict(flat11, species_temperatures={"D": 20.0, "T": 20.0})
    a11 = ps.particles(flat11, 512, seed=17)
    b11 = ps.particles(pair11, 512, seed=17)
    recovered11 = all(bool((a11[key] == b11[key]).all()) for key in ("x", "y", "z", "energy"))
    rows.append(["P11 equal-T stream recovery", "-", "exact", _check(recovered11, "P11 stream")])
    cards_a11 = ps.emit_source_cards(flat11, bins=15)
    cards_b11 = ps.emit_source_cards(pair11, bins=15)
    same_cards11 = (
        cards_b11["sdef"]["card"] == cards_a11["sdef"]["card"]
        and cards_b11["serpent"]["card"] == cards_a11["serpent"]["card"]
    )
    rows.append(["P11 equal-T card recovery", "-", "exact", _check(same_cards11, "P11 cards")])
    notes.append(
        f"P11 per-species ion temperatures (T_D={t_d11:g} keV, T_T={t_t11:g} keV, "
        f"70/30 blend, n={N}): D-T reacts at T_DT={t_d11 + 0.4 * (t_t11 - t_d11):g} keV, "
        f"D-D at T_D (sampled moments and branch share vs the pair-temperature "
        f"quadrature); the uniform equal pair on a flat 20 keV plasma recovers "
        f"the shared kernel bit-for-bit (stream and cards)."
    )
    return rows, notes


LATTICE_TWO_POINT = {
    "kind": "lattice",
    "reaction": "dt",
    "points": [
        {"position": [300.0, 0.0, 25.0], "rate": 2.0, "ion_temperature_kev": 0.0},
        {"position": [-300.0, 0.0, 25.0], "rate": 1.0, "ion_temperature_kev": 0.0},
    ],
}


def _ring_lattice_spec(k: int = 64) -> dict:
    """Uniform ring lattice: k unit-rate D-T mono nodes at R=300 cm, z=25 cm."""
    return {
        "kind": "lattice",
        "reaction": "dt",
        "points": [
            {
                "position": [
                    300.0 * math.cos(2.0 * math.pi * i / k),
                    300.0 * math.sin(2.0 * math.pi * i / k),
                    25.0,
                ],
                "rate": 1.0,
                "ion_temperature_kev": 0.0,
            }
            for i in range(k)
        ],
    }


def lattice_gates() -> tuple[list[list[str]], list[str]]:
    """P12 on the arbitrary-3D birth-rate lattice; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    # P12a: hand-vector totals and means on the two-point cloud (D-T mono).
    out = ps.particles(LATTICE_TWO_POINT, 90_000, seed=42)
    n12 = len(out["x"])
    ok_e = bool((out["energy"] == 14.021).all())
    rows.append(["P12 hand-vector mono energy", "-", "exact", _check(ok_e, "P12 energy")])
    se_x = 200.0 / math.sqrt(n12)
    mean_x = float(np.mean(out["x"]))
    rows.append(
        [
            "P12 hand-vector mean x",
            fmt(mean_x - 100.0),
            "< 8 sigma",
            _check(abs(mean_x - 100.0) < 8.0 * se_x, "P12 mean x"),
        ]
    )
    frac_a = float(np.mean(out["x"] > 0))
    se_a = math.sqrt((2.0 / 3.0) * (1.0 / 3.0) / n12)
    rows.append(
        [
            "P12 rate-2 node share",
            fmt(frac_a - 2.0 / 3.0),
            "< 8 sigma",
            _check(abs(frac_a - 2.0 / 3.0) < 8.0 * se_a, "P12 share"),
        ]
    )

    # P12b: single-point lattice reproduces the point source bit-for-bit.
    lattice_1 = {
        "kind": "lattice",
        "reaction": "dt",
        "points": [{"position": [1.0, -2.0, 3.5], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    point_1 = {
        "kind": "point",
        "position": [1.0, -2.0, 3.5],
        "reaction": "dt",
        "ion_temperature_kev": 20.0,
    }
    a = ps.particles(lattice_1, 512, seed=9)
    b = ps.particles(point_1, 512, seed=9)
    recovered = all(bool((a[key] == b[key]).all()) for key in ("x", "y", "z", "energy"))
    rows.append(["P12 single-point recovery", "-", "exact", _check(recovered, "P12 recovery")])

    # P12c: ring-limit moment closure (axisymmetric convergence anchor).
    out = ps.particles(_ring_lattice_spec(), N, seed=42)
    radius = np.hypot(out["x"], out["y"])
    ok = bool(np.all(out["z"] == 25.0)) and bool(np.all(np.abs(radius - 300.0) < 1e-9 * 300.0))
    rows.append(["P12 ring radius/height", "-", "exact", _check(ok, "P12 ring")])
    se = 300.0 / math.sqrt(2 * N)
    mean_x = float(np.mean(out["x"]))
    rows.append(
        [
            "P12 ring <x>",
            fmt(mean_x),
            "< 8 sigma",
            _check(abs(mean_x) < 8.0 * se, "P12 ring x"),
        ]
    )

    # P12d: symmetry-fold replication — one base node with 4 field periods
    # spreads uniformly over the 4 copies; totals carry the period factor.
    base4 = {
        "kind": "lattice",
        "reaction": "dt",
        "field_periods": 4,
        "base_angle": 0.0,
        "points": [{"position": [300.0, 0.0, 0.0], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    m = 20_000
    out = ps.particles(base4, m, seed=5)
    radius = np.hypot(out["x"], out["y"])
    ok = bool(np.all(np.abs(radius - 300.0) < 1e-9 * 300.0))
    rows.append(["P12 fold radius", "-", "exact", _check(ok, "P12 fold r")])
    quad = ((out["x"] > 0).astype(int) + 2 * (out["y"] > 0).astype(int)).astype(int)
    quad_ok = True
    for q in range(4):
        frac = float(np.mean(quad == q))
        if abs(frac - 0.25) >= 8.0 * math.sqrt(0.25 * 0.75 / m):
            quad_ok = False
    rows.append(["P12 fold uniform copies", "-", "uniform", _check(quad_ok, "P12 copies")])
    acc = ps.proton_accounting(base4)
    totals_ok = acc["proton_total"] == 0.0 and abs(acc["neutron_total"] - 4.0) < 1e-12
    rows.append(["P12 fold totals", "-", "4x base", _check(totals_ok, "P12 totals")])

    # P12e: card round trip with the lattice drift row; mixture branch fire.
    cards = ps.emit_source_cards(LATTICE_TWO_POINT, bins=8)
    sec_parsed = mcnp.parse_sdef(cards["sdef"]["card"])
    card_ok = (
        sec_parsed["card"] == cards["sdef"]["card"]
        and sec_parsed["rad"] == "D1"
        and sec_parsed["ext"] == "D2"
        and sec_parsed["erg"] == "D3"
        and len(sec_parsed["distributions"]) == 3
    )
    rows.append(["P12 lattice card round trip", "-", "byte-identical", _check(card_ok, "P12 card")])
    quantities = [row["quantity"] for row in cards["sdef"]["drift"]]
    drift_ok = quantities == [
        "emission probability",
        "spatial marginals",
        "joint correlation",
        "lattice discretization",
    ]
    rows.append(["P12 lattice drift rows", "-", "four rows", _check(drift_ok, "P12 drift")])
    blend = {
        "kind": "lattice",
        "fuel": {"D": 0.7, "T": 0.3},
        "points": [{"position": [0.0, 0.0, 0.0], "rate": 1.0, "ion_temperature_kev": 20.0}],
    }
    out = ps.particles(blend, N, seed=11)
    frac_dd = float(np.mean(out["energy"] < 10.0))
    rows.append(
        [
            "P12 mixture D-D branch fires",
            fmt(frac_dd),
            "0.002-0.015",
            _check(0.002 < frac_dd < 0.015, "P12 branches"),
        ]
    )
    notes.append(
        f"P12 arbitrary-3D birth-rate lattice: two-point hand vectors (total 3, "
        f"mean x 100 cm, mono 14.021 MeV), single-point bit-for-bit recovery of "
        f"the landed point source, ring-limit moment closure (n={N}), 4-period "
        f"symmetry-fold replication with period-scaled totals, card round trip "
        f"with the lattice drift row, and 70/30 mixture branch fire."
    )
    return rows, notes


def tail_gates() -> tuple[list[list[str]], list[str]]:
    """P13 on the deuterium hot-tail fraction; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    # P13a: pinned tail spec — 70/30 blend at T_D = 20 keV, T_T = 30 keV with
    # a 5% deuterium hot tail at 60 keV. Sampled mean energy, radius moment,
    # and D-D branch share sit on the sub-rate-rule quadrature.
    t_d13, t_t13, eta13, t_tail13 = 20.0, 30.0, 0.05, 60.0
    tail13 = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": t_d13, "T": t_t13},
        deuterium_tail={"fraction": eta13, "temperature_kev": t_tail13},
    )
    ours13 = ps.particles(tail13, N, seed=84)
    quad13 = _quadrature_moments_mixture_tail(
        PARAMETRIC_SPEC, 0.7, 0.3, t_d13, t_t13, eta13, t_tail13
    )
    err_e13 = rel_diff(float(np.mean(ours13["energy"])), quad13["mean_e"])
    rows.append(
        ["P13 tail mean birth energy", fmt(err_e13), "< 1e-2", _check(err_e13 < 1e-2, "P13 E")]
    )
    major13 = np.hypot(ours13["x"], ours13["y"])
    r13 = _recover_minor_radius(PARAMETRIC_SPEC, major13, ours13["z"])
    err_r13 = rel_diff(float(np.mean(r13**2)), quad13["<r^2>"])
    rows.append(["P13 tail birth <r^2>", fmt(err_r13), "< 1e-2", _check(err_r13 < 1e-2, "P13 r2")])
    frac_dd13 = float(np.mean(ours13["energy"] < 10.0))
    se_dd13 = math.sqrt(quad13["p_dd"] * (1.0 - quad13["p_dd"]) / N)
    rows.append(
        [
            "P13 tail D-D branch share",
            fmt(frac_dd13 - quad13["p_dd"]),
            "< 8 sigma",
            _check(abs(frac_dd13 - quad13["p_dd"]) < 8.0 * se_dd13, "P13 share"),
        ]
    )

    # P13b: the card drift row's tail neutron share reproduces the
    # independent quadrature share (cross-checks total_tail_strength; the
    # note prints 6 decimals, so the 1e-6 window covers the rounding).
    cards13 = ps.emit_source_cards(tail13, bins=15)
    tail_row13 = next(
        row for row in cards13["sdef"]["drift"] if row["quantity"] == "deuterium tail"
    )
    share13 = float(tail_row13["note"].split("volume-integrated source = ")[1])
    rows.append(
        [
            "P13 tail neutron share",
            fmt(share13 - quad13["p_tail"]),
            "< 1e-6",
            _check(abs(share13 - quad13["p_tail"]) < 1e-6, "P13 tail share"),
        ]
    )
    drift_ok13 = [row["quantity"] for row in cards13["sdef"]["drift"]] == [
        "emission probability",
        "spatial marginals",
        "joint correlation",
        "deuterium tail",
    ] and bool(tail_row13["reparsed"])
    rows.append(["P13 tail drift row", "-", "present", _check(drift_ok13, "P13 drift")])

    # P13c: eta = 0 reproduces the no-tail mixture kernel bit-for-bit.
    base13 = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": t_d13, "T": t_t13},
    )
    zero13 = dict(base13, deuterium_tail={"fraction": 0.0, "temperature_kev": t_tail13})
    a13 = ps.particles(base13, 512, seed=17)
    b13 = ps.particles(zero13, 512, seed=17)
    recovered13 = all(bool((a13[key] == b13[key]).all()) for key in ("x", "y", "z", "energy"))
    rows.append(["P13 zero-tail stream recovery", "-", "exact", _check(recovered13, "P13 stream")])
    cards_a13 = ps.emit_source_cards(base13, bins=15)
    cards_b13 = ps.emit_source_cards(zero13, bins=15)
    same_cards13 = (
        cards_b13["sdef"]["card"] == cards_a13["sdef"]["card"]
        and cards_b13["serpent"]["card"] == cards_a13["serpent"]["card"]
    )
    rows.append(["P13 zero-tail card recovery", "-", "exact", _check(same_cards13, "P13 cards")])

    # P13d: tail parameters are loud — out-of-range fraction, negative tail
    # temperature, and a tail without a fuel mixture never reach the sampler.
    loud13 = True
    for bad in (
        {"fraction": 1.5, "temperature_kev": 60.0},
        {"fraction": 0.05, "temperature_kev": -1.0},
    ):
        try:
            ps.particles(dict(base13, deuterium_tail=bad), 4, seed=0)
            loud13 = False
        except ValueError:
            pass
    try:
        ps.particles(
            dict(PARAMETRIC_SPEC, deuterium_tail={"fraction": 0.05, "temperature_kev": 60.0}),
            4,
            seed=0,
        )
        loud13 = False
    except ValueError:
        pass
    rows.append(["P13 tail loud errors", "-", "raise", _check(loud13, "P13 loud")])
    notes.append(
        f"P13 deuterium hot tail (70/30 blend, T_D={t_d13:g} keV, T_T={t_t13:g} keV, "
        f"eta={eta13:g} at T_tail={t_tail13:g} keV, n={N}): sampled moments and "
        f"D-D branch share vs the sub-rate-rule quadrature, the card drift "
        f"tail-neutron share vs the same quadrature, bit-for-bit zero-fraction "
        f"recovery of the no-tail kernel (stream and cards), and loud tail errors."
    )
    return rows, notes


def oracle_check_openmc_plasma_source() -> tuple[list[list[str]], list[str], bool]:
    """O1-O13 vs the upstream MIT-licensed package. Returns (rows, notes, skipped)."""
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
            r_map, z_map = _map_rz(g, a_cm, theta)
            worst = max(
                worst, rel_diff(float(r_map), float(r_up)), rel_diff(float(z_map), float(z_up))
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
    from NeSST.spectral_model import TT_model, reac_DD, reac_DT, reac_TT

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

    # O9: tritium-rich (10/90) blend — the landed two-branch kernel vs the
    # three-branch upstream model. Nucleide carries no T-T branch (no
    # publishable closed form exists; the oracle's vendored Hale table +
    # continuum files stay out of scope), so the sampled stream must sit on
    # the two-branch quadrature while the three-branch reference pulls away
    # by exactly the T-T share — the deferred branch quantified, not
    # ignored. Protons have no oracle (upstream models neutrons only); P9
    # pins the D(d,p)T bookkeeping analytically.
    f_d9, f_t9 = 0.1, 0.9
    blend9 = dict(PARAMETRIC_SPEC, fuel={"D": f_d9, "T": f_t9})
    ours9 = ps.particles(blend9, N, seed=79)
    mean9 = float(np.mean(ours9["energy"]))
    quad2 = _quadrature_upstream_mixture(
        g, tokamak_ion_density, tokamak_ion_temperature, reac_DT, reac_DD, f_d9, f_t9
    )
    quad3 = _quadrature_upstream_mixture_tt(
        g,
        tokamak_ion_density,
        tokamak_ion_temperature,
        reac_DT,
        reac_DD,
        reac_TT,
        TT_model,
        f_d9,
        f_t9,
    )
    err_2 = rel_diff(mean9, quad2["mean_e"])
    rows.append(
        ["O9 T-rich two-branch mean", fmt(err_2), "< 1e-2", _check(err_2 < 1e-2, "O9 2-branch")]
    )
    predicted = abs(quad3["mean_e"] - quad2["mean_e"])
    measured = abs(mean9 - quad3["mean_e"])
    rows.append(
        ["O9 T-T pull significant", fmt(measured), "> 0.1 MeV", _check(measured > 0.1, "O9 pull")]
    )
    err_tt = abs(measured - predicted) / predicted
    rows.append(
        [
            "O9 T-T pull at predicted scale",
            fmt(err_tt),
            "< 2e-2",
            _check(err_tt < 2e-2, "O9 scale"),
        ]
    )
    notes.append(
        f"O9 10/90 D/T blend, n={N}: sampled mean birth energy sits on the "
        f"two-branch upstream quadrature (no T-T on either side) while the "
        f"three-branch reference (NeSST reac_TT x2 multiplicity + Brune "
        f"continuum mean ~4.6 MeV) pulls away by the T-T-predicted "
        f"{predicted:.3f} MeV (TT neutron share {quad3['tt_share']:.4f}); "
        f"measured pull {measured:.3f} MeV."
    )

    # O10: toroidal sector vs the upstream quadrature. A partial sector must
    # leave the (r, z, E) physics untouched — the sampled energy and radius
    # moments sit on the same O7 two-branch quadrature built from the
    # upstream map/profiles — while birth angles are uniform over the sector
    # (the upstream full-torus phi draw has no sector spelling, so the angle
    # marginal itself rests on the always-run P10 gates).
    start10, rotation10 = 0.5, 1.5
    sector10 = dict(PARAMETRIC_SPEC, start_angle=start10, rotation_angle=rotation10)
    ours10 = ps.particles(sector10, N, seed=80)
    quad10 = _quadrature_upstream(g, tokamak_ion_density, tokamak_ion_temperature, reac_DT)
    err_e10 = rel_diff(float(np.mean(ours10["energy"])), quad10["mean_e"])
    rows.append(
        ["O10 sector mean birth energy", fmt(err_e10), "< 1e-2", _check(err_e10 < 1e-2, "O10 E")]
    )
    major10 = np.hypot(ours10["x"], ours10["y"])
    r10 = _recover_minor_radius(g, major10, ours10["z"])
    err_r10 = rel_diff(float(np.mean(r10**2)), quad10["<r^2>"])
    rows.append(
        ["O10 sector birth <r^2>", fmt(err_r10), "< 1e-2", _check(err_r10 < 1e-2, "O10 r2")]
    )
    phi10 = np.arctan2(ours10["y"], ours10["x"]) % (2.0 * math.pi)
    mean10 = float(np.mean(phi10))
    want10 = start10 + rotation10 / 2.0
    se10 = rotation10 / math.sqrt(12.0 * N)
    rows.append(
        [
            "O10 sector mean angle",
            fmt(mean10 - want10),
            "< 8 sigma",
            _check(abs(mean10 - want10) < 8.0 * se10, "O10 angle"),
        ]
    )
    notes.append(
        f"O10 partial toroidal sector (start {start10} rad, rotation {rotation10} rad), "
        f"n={N}: (r, z, E) moments sit on the O7 upstream quadrature (the sector "
        f"only remaps the toroidal draw) and birth angles are uniform over the sector."
    )

    # O11: per-species ion temperatures vs the upstream quadrature. A 70/30
    # blend at T_D = 20 keV, T_T = 30 keV reacts D-T at T_DT = 24 keV and D-D
    # at T_D through the scalar NeSST reactivity functions (the upstream
    # package has no distinct-temperature spelling — the pair temperatures
    # enter there, so the cross-check is the independent NeSST transcription
    # plus the upstream profiles and map).
    t_d11, t_t11 = 20.0, 30.0
    blend11 = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": t_d11, "T": t_t11},
    )
    ours11 = ps.particles(blend11, N, seed=82)
    quad11 = _quadrature_upstream_mixture_species(
        g, tokamak_ion_density, reac_DT, reac_DD, 0.7, 0.3, t_d11, t_t11
    )
    err_e11 = rel_diff(float(np.mean(ours11["energy"])), quad11["mean_e"])
    rows.append(
        ["O11 species mean birth energy", fmt(err_e11), "< 1e-2", _check(err_e11 < 1e-2, "O11 E")]
    )
    major11 = np.hypot(ours11["x"], ours11["y"])
    r11 = _recover_minor_radius(g, major11, ours11["z"])
    err_r11 = rel_diff(float(np.mean(r11**2)), quad11["<r^2>"])
    rows.append(
        ["O11 species birth <r^2>", fmt(err_r11), "< 1e-2", _check(err_r11 < 1e-2, "O11 r2")]
    )
    notes.append(
        f"O11 70/30 D/T blend at T_D={t_d11:g} keV, T_T={t_t11:g} keV "
        f"(D-T at T_DT={t_d11 + 0.4 * (t_t11 - t_d11):g} keV, D-D at T_D), n={N}: "
        f"sampled moments vs quadrature over the upstream map/profiles with "
        f"NeSST reac_DT + reac_DD evaluated at the pair temperatures."
    )
    # O12: lattice ring limit vs the upstream models where they apply. A dense
    # monoenergetic ring cloud (360 unit-rate D-T nodes at R = 3 m, z = 0.5 m)
    # sits on the upstream ring geometry (r delta at R, uniform azimuth via
    # <x^2> = R^2/2) and the Ballabio mean helper (Ti = 0 line).
    ring_lat = _ring_lattice_spec(360)
    ours12 = ps.particles(ring_lat, N, seed=83)
    err_r12 = rel_diff(float(np.mean(np.hypot(ours12["x"], ours12["y"]))), r_ref)
    rows.append(
        ["O12 lattice ring radius", fmt(err_r12), "< 1e-9", _check(err_r12 < 1e-9, "O12 radius")]
    )
    err_phi12 = rel_diff(float(np.mean(ours12["x"] ** 2)), float(np.mean(up_x**2)))
    rows.append(
        ["O12 lattice azimuth <x^2>", fmt(err_phi12), "< 1e-2", _check(err_phi12 < 1e-2, "O12 phi")]
    )
    mean_ev12 = neutron_energy_mean(ion_temperature=0.0, reaction="DT")
    err_e12 = rel_diff(float(np.mean(ours12["energy"])), mean_ev12 / 1e6)
    rows.append(
        ["O12 lattice mono energy", fmt(err_e12), "< 1e-12", _check(err_e12 < 1e-12, "O12 E")]
    )
    notes.append(
        f"O12 lattice ring limit (360-node monoenergetic cloud, n={N}): radius "
        f"and azimuth moments vs the upstream ring CylindricalIndependent "
        f"parameters, mean birth energy vs the upstream Ballabio mean helper."
    )
    # O13: deuterium hot tail vs the upstream quadrature. A 70/30 blend at
    # T_D = 20 keV, T_T = 30 keV with a 5% tail at 60 keV reacts the five
    # bulk/tail sub-pairs through the scalar NeSST reactivity functions at
    # the sub-pair temperatures (the upstream package has no tail spelling —
    # the O11 stance: the cross-check is the independent NeSST transcription
    # plus the upstream profiles and map).
    t_d13, t_t13, eta13, t_tail13 = 20.0, 30.0, 0.05, 60.0
    blend13 = dict(
        PARAMETRIC_SPEC,
        fuel={"D": 0.7, "T": 0.3},
        species_temperatures={"D": t_d13, "T": t_t13},
        deuterium_tail={"fraction": eta13, "temperature_kev": t_tail13},
    )
    ours13 = ps.particles(blend13, N, seed=85)
    quad13 = _quadrature_upstream_mixture_tail(
        g, tokamak_ion_density, reac_DT, reac_DD, 0.7, 0.3, t_d13, t_t13, eta13, t_tail13
    )
    err_e13 = rel_diff(float(np.mean(ours13["energy"])), quad13["mean_e"])
    rows.append(
        ["O13 tail mean birth energy", fmt(err_e13), "< 1e-2", _check(err_e13 < 1e-2, "O13 E")]
    )
    major13 = np.hypot(ours13["x"], ours13["y"])
    r13 = _recover_minor_radius(g, major13, ours13["z"])
    err_r13 = rel_diff(float(np.mean(r13**2)), quad13["<r^2>"])
    rows.append(["O13 tail birth <r^2>", fmt(err_r13), "< 1e-2", _check(err_r13 < 1e-2, "O13 r2")])
    notes.append(
        f"O13 70/30 D/T blend with a 5% deuterium hot tail at T_tail={t_tail13:g} keV "
        f"(T_D={t_d13:g} keV, T_T={t_t13:g} keV, n={N}): sampled moments vs "
        f"quadrature over the upstream map/profiles with NeSST reac_DT + "
        f"reac_DD evaluated at the five bulk/tail sub-pair temperatures."
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


def _quadrature_upstream_mixture_species(
    g: dict,
    density_fn,
    reac_dt,
    reac_dd,
    f_d: float,
    f_t: float,
    t_d: float,
    t_t: float,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Quadrature of n^2·[f_D·f_T·reac_DT(T_DT) + (f_D^2/2)·reac_DD(T_D)]·R·|J|
    using the upstream density profile and NeSST reactivities at the uniform
    pair temperatures (T_DT = T_D + (2/5)·(T_T − T_D)); the volume element
    uses the in-script Miller transcription (gate O4). The upstream package
    has no distinct-temperature spelling, so the pair temperatures enter
    through its scalar reactivity functions — the cross-check is the
    independent NeSST transcription plus the upstream profiles."""
    t_dt = t_d + 0.4 * (t_t - t_d)
    am = g["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = 0.0
    sv_dt = float(reac_dt(t_dt * 1e3)) * 1e6  # m^3/s -> cm^3/s
    sv_dd = float(reac_dd(t_d * 1e3)) * 1e6
    w_dt = f_d * f_t * sv_dt
    w_dd = f_d * f_d / 2.0 * sv_dd
    mu_dt, _ = _ballabio_moments("dt", t_dt)
    mu_dd, _ = _ballabio_moments("dd", t_d)
    p_dd = w_dd / (w_dt + w_dd) if w_dt + w_dd > 0.0 else 0.0
    mean_e = p_dd * mu_dd + (1.0 - p_dd) * mu_dt
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
        strength = n_m3 * n_m3 * 1e-12 * (w_dt + w_dd)  # (n·1e-6)^2 = n^2·1e-12
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = strength * _volume_element(g, r, theta) * dr * dt
            den += w
            num_r2 += w * r * r
            num_e += w * mean_e
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den}


def _quadrature_upstream_mixture_tail(
    g: dict,
    density_fn,
    reac_dt,
    reac_dd,
    f_d: float,
    f_t: float,
    t_d: float,
    t_t: float,
    eta: float,
    t_tail: float,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Quadrature of the deuterium hot-tail sub-rate rule using the upstream
    density profile and NeSST reactivities at the sub-pair temperatures
    (T_DT, T_DTt for D-T; T_D, T_mix, T_tail for D-D). The upstream package
    has no tail spelling — like O11, the sub-pair temperatures enter through
    its scalar reactivity functions, so the cross-check is the independent
    NeSST transcription plus the upstream profiles. Returns the
    strength-weighted <r^2> and the sub-branch-weighted mean birth energy."""
    t_dt = t_d + 0.4 * (t_t - t_d)
    t_dtt = t_tail + 0.4 * (t_t - t_tail)
    t_mix = 0.5 * (t_d + t_tail)
    am = g["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_r2 = den = num_e = 0.0
    sv_dt = float(reac_dt(t_dt * 1e3)) * 1e6  # m^3/s -> cm^3/s
    sv_dt_tail = float(reac_dt(t_dtt * 1e3)) * 1e6
    sv_dd = float(reac_dd(t_d * 1e3)) * 1e6
    sv_dd_mix = float(reac_dd(t_mix * 1e3)) * 1e6
    sv_dd_tail = float(reac_dd(t_tail * 1e3)) * 1e6
    w = [
        f_d * f_t * (1.0 - eta) * sv_dt,
        f_d * f_t * eta * sv_dt_tail,
        f_d * f_d / 2.0 * (1.0 - eta) ** 2 * sv_dd,
        f_d * f_d / 2.0 * 2.0 * eta * (1.0 - eta) * sv_dd_mix,
        f_d * f_d / 2.0 * eta**2 * sv_dd_tail,
    ]
    total = sum(w)
    mus = [
        _ballabio_moments("dt", t_dt)[0],
        _ballabio_moments("dt", t_dtt)[0],
        _ballabio_moments("dd", t_d)[0],
        _ballabio_moments("dd", t_mix)[0],
        _ballabio_moments("dd", t_tail)[0],
    ]
    mean_e = sum(x / total * mu for x, mu in zip(w, mus, strict=True))
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
        strength = n_m3 * n_m3 * 1e-12 * total
        for j in range(n_t):
            theta = (j + 0.5) * dt
            cell = strength * _volume_element(g, r, theta) * dr * dt
            den += cell
            num_r2 += cell * r * r
            num_e += cell * mean_e
    return {"<r^2>": num_r2 / den, "mean_e": num_e / den}


def _quadrature_upstream_mixture_tt(
    g: dict,
    density_fn,
    temperature_fn,
    reac_dt,
    reac_dd,
    reac_tt,
    tt_model,
    f_d: float,
    f_t: float,
    n_r: int = 400,
    n_t: int = 400,
) -> dict[str, float]:
    """Three-branch quadrature adding the oracle's T-T model to
    `_quadrature_upstream_mixture`: the T-T neutron weight is
    `f_T^2·reac_TT` (the oracle's `n_T^2/2·reac_TT` fuel density with the
    x2 neutron multiplicity folded in — the pinned normalization), and the
    T-T mean birth energy is the oracle continuum mean (Brune model,
    ~4.6 MeV, nearly T-independent) evaluated per cell where the branch
    burns. Returns the three-branch mean plus the T-T neutron share."""
    import numpy as _np

    e_grid = _np.linspace(1e3, 12e6, 200)  # eV (oracle continuum grid)
    am = g["minor_radius"]
    dr = am / n_r
    dt = 2.0 * math.pi / n_t
    num_e = den = num_tt = 0.0
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
        sv_tt = float(reac_tt(t_kev * 1e3)) * 1e6
        w_dt = f_d * f_t * sv_dt
        w_dd = f_d * f_d / 2.0 * sv_dd
        w_tt = f_t * f_t * sv_tt
        mu_dt, _ = _ballabio_moments("dt", t_kev)
        mu_dd, _ = _ballabio_moments("dd", t_kev)
        if w_tt > 0.0:
            spec = _np.asarray(tt_model.spec(e_grid, max(t_kev, 0.5) * 1e3, "Brune"))
            mu_tt = float(_np.trapezoid(e_grid * spec, e_grid) / _np.trapezoid(spec, e_grid))
            mu_tt /= 1e6  # eV -> MeV
        else:
            mu_tt = mu_dt  # weightless branch; value unused
        tot = w_dt + w_dd + w_tt
        mean_e = (w_dt * mu_dt + w_dd * mu_dd + w_tt * mu_tt) / tot if tot > 0.0 else mu_dt
        strength = n_m3 * n_m3 * 1e-12 * tot
        for j in range(n_t):
            theta = (j + 0.5) * dt
            w = strength * _volume_element(g, r, theta) * dr * dt
            den += w
            num_e += w * mean_e
            num_tt += w * (w_tt / tot if tot > 0.0 else 0.0)
    return {"mean_e": num_e / den, "tt_share": num_tt / den}


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
        "transcription, sampled moments vs fine quadrature), P8 on a "
        "pinned 70/30 D/T fuel blend (mixture-rule moments and the D-D "
        "branch share vs in-script mixture quadrature), and P9 on the "
        "D(d,p)T proton bookkeeping (exact pure-fuel ratios, the 70/30 "
        "share, the no-T recovery anchor, the loud pure-T error), and P10 "
        "on toroidal sectors (in-sector uniform births, bit-for-bit "
        "full-rotation recovery, the PHI=D4 marginal with its drift row), and P11 "
        "on per-species ion temperatures (pair-temperature moments and branch "
        "share, bit-for-bit equal-T recovery), and P12 "
        "on the arbitrary-3D birth-rate lattice (hand-vector totals and "
        "means, single-point recovery, ring-limit closure, symmetry-fold "
        "replication, card round trip, mixture branch fire), and P13 "
        "on the deuterium hot-tail fraction (sub-rate-rule moments and "
        "branch share, the drift tail-neutron share, zero-fraction "
        "recovery) — "
        "always run — plus container-only cross-checks O1-O13 against the upstream "
        "openmc-plasma-source package and NeSST (Ballabio helpers, sampled "
        "ring/point sources, Miller map, Fausser profiles, reactivities, "
        "end-to-end parametric moments, the pinned 70/30 blend, the "
        "tritium-rich 10/90 T-T divergence probe, the partial-sector "
        "moment cross-check, the distinct-temperature blend probe, "
        "the lattice ring-limit probe, and the hot-tail blend probe; protons have no oracle — "
        "upstream models neutrons only)."
    )
    rows1, notes1 = analytic_gates()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows1)
    rows2, notes2 = parametric_gates()
    for note in notes2:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows2)
    rows2b, notes2b = lattice_gates()
    for note in notes2b:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows2b)
    rows2c, notes2c = tail_gates()
    for note in notes2c:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows2c)
    rows3, notes3, skipped = oracle_check_openmc_plasma_source()
    for note in notes3:
        report.prose(note)
    if skipped:
        skip_rows = [
            [f"{gate}", "SKIP (openmc-plasma-source unavailable)"]
            for gate in (
                "O1",
                "O2",
                "O3",
                "O4",
                "O5",
                "O6",
                "O7",
                "O8",
                "O9",
                "O10",
                "O11",
                "O12",
                "O13",
            )
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
