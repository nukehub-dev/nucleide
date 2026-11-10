"""1D tritium-kernel parity probe (`nucleide.tritium` vs FESTIM 2.1, black-box).

Two parts:

1. Always-run Nucleide legs: VV Val-1da (the published analytic
   permeation breakthrough, ``J(L,t)/J_ss = 1 + 2*sum_{n>=1} (-1)^n
   exp(-D n^2 pi^2 t / L^2)`` with ``t_lag = L^2/6D``) and M1
   (manufactured-solution order check of the steady diffusion operator
   with an explicit source). These run everywhere and gate the Nucleide
   side of the probe.
2. Black-box parity legs F1-F4 against FESTIM 2.1 (Apache-2.0,
   ``festim==2.1``, v2.x API): F1 trap-free Dirichlet-Dirichlet, F2 one
   McNabb-Foster trap, F3 Sieverts ends, F4 two-layer stack. Each leg
   is computed independently on both sides (Nucleide fine/coarse plus
   FESTIM fine/coarse when the oracle runs) and compared on J_ss,
   J(t), t_lag, and the steady c(x) profile. Recombination ends are
   excluded by design: the Nucleide exact-face law versus the FESTIM
   surface-flux formulation is an uncharacterized mismatch zone.

Probe environment (never the main validation container, never the Rust
workspace): FESTIM 2.1 needs DOLFINx, so the probe runs in the separate
probe-only image ``dolfinx/dolfinx:stable`` (see ``validation/README.md``,
section "Tritium/FESTIM parity probe (probe-only image)"). The image has
no Rust toolchain, so run ``maturin develop`` first to write
``python/nucleide/_internal.abi3.so`` into the mounted checkout::

    maturin develop
    docker run --rm -v "$PWD:/work" -w /work dolfinx/dolfinx:stable \\
        bash -c "pip install festim==2.1 && python validation/tritium_vs_festim.py"

Outside that image ``import festim`` fails (or carries the wrong version)
and every F leg is reported as SKIP with its reason — never silently,
never as a pass. FESTIM is used black-box only: this script drives the
documented v2.1 surface (``HydrogenTransportProblem`` /
``HydrogenTransportProblemDiscontinuous``, ``Mesh1D``,
``Material(D_0, E_D)``, ``Species`` + ``Reaction`` for McNabb-Foster,
``FixedConcentrationBC`` / ``SievertsBC``, ``Settings`` + ``Stepsize``,
``SurfaceFlux`` / ``Profile1DExport``) and reads only the exported flux
and profile numbers. No FESTIM code is copied, vendored, or referenced.

All probe parameters are hand-picked synthetic values (SI units,
isothermal: E_D = E_S = E_k = E_p = 0, so the Arrhenius factors are 1 and
the run temperature is immaterial). Gates (discretization-dominated,
10^2-10^4x looser than the algebraic crate gates): J_ss 2%, J(t) 3% at
gated times (t >= 0.2*t_lag and above a 1e-3*J_ss floor), t_lag 5%, steady
c(x) 2% max-norm. Each leg must first self-converge on both sides to
< 1/3 of its gate; a leg that has not self-converged is INCONCLUSIVE,
never PASS.
"""

from __future__ import annotations

import importlib.metadata
import inspect
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Probe-image bootstrap: inside the probe-only dolfinx image the only
# nucleide on disk is this repo's ``python/`` tree (no maturin-develop
# editable install), so add it to sys.path (harmless when already present).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))

import common
import numpy as np
from common import Report, fmt, rel_diff

import nucleide.tritium as tri

_REAL_ENVIRONMENT = common.environment


def _tolerant_environment() -> dict[str, str]:
    """Best-effort environment metadata when distributions are missing.

    Follow the ``activation_vs_refs.py`` precedent: the canonical container
    always has nucleide/PyNE/OpenMC installed, so the real collector
    succeeds there unchanged. The probe image has no PyNE/OpenMC/nucleide
    distributions, so fall back to "not installed" markers instead of
    raising (same convention as ``common._dist_version``).
    """
    try:
        return _REAL_ENVIRONMENT()
    except Exception:
        try:
            nucleide = importlib.metadata.version("nucleide")
        except importlib.metadata.PackageNotFoundError:
            nucleide = "not installed"
        try:
            pyne = importlib.metadata.version("pyne")
        except importlib.metadata.PackageNotFoundError:
            pyne = "not installed"
        try:
            openmc = importlib.metadata.version("openmc")
        except importlib.metadata.PackageNotFoundError:
            openmc = "not installed"
        return {
            "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "nucleide": nucleide,
            "pyne": pyne,
            "openmc": openmc,
        }


common.environment = _tolerant_environment

FAILURES = 0

# Synthetic probe parameters (SI; hand-picked round values).
L = 1.0e-3  # slab thickness [m]
D = 1.0e-9  # diffusivity [m^2/s]
C0 = 1.0  # inlet concentration [mol/m^3] (F1/F2/F4)
CL = 0.0  # outlet concentration [mol/m^3]
TLAG = L * L / (6.0 * D)  # 166.666... s

# Trap leg (F2): K = k0/p0 = 1.0 m^3/mol, K*N = 2, Oriani D_eff = D/3.
TRAP_K0 = 0.5  # trapping rate [1/(mol/m^3)/s]
TRAP_P0 = 0.5  # detrapping rate [1/s]
TRAP_N = 2.0  # site density [mol/m^3]
TLAG_TRAP = L * L / (6.0 * tri.oriani(D, TRAP_K0 / TRAP_P0, TRAP_N))  # 500 s

# Sieverts leg (F3): c = K_S*sqrt(p) -> 8.0 / 2.0 mol/m^3.
SIEVERT_K = 2.0  # [mol/m^3/Pa^0.5]
P_IN = 16.0  # [Pa]
P_OUT = 1.0  # [Pa]

# Two-layer leg (F4): second stack line of the layered gate vectors.
LAYERS_2 = [
    {"thickness": 5e-4, "cells": 0, "D": 1e-9, "solubility": 2.0},
    {"thickness": 5e-4, "cells": 0, "D": 5e-10, "solubility": 0.5},
]

# G2-grid time outputs. The trap leg extends past 1000 s (t_lag_eff=500 s).
# G2-grid outputs are used by the Val-1da leg; parity legs use LEG_TIMES.
G2_TIMES = [1.0, 3.1622776601683795, 10.0, 31.622776601683793, 100.0, 316.22776601683796, 1000.0]
T_END_TRAP = 2500.0


# Log-spaced output grids for the parity legs (dense late-time sampling
# for the cumulative-flux t_lag fit; stepping cost is set by dt_max).
def _log_grid(t_end: float, n: int = 25) -> list[float]:
    grid = [float(t) for t in np.logspace(0.0, math.log10(t_end), n)]
    grid[-1] = t_end
    return grid


LEG_TIMES: dict[str, list[float]] = {
    "F1": _log_grid(1000.0),
    "F2": _log_grid(T_END_TRAP),
    "F3": _log_grid(1000.0),
    "F4": _log_grid(1000.0),
}

# Nucleide fine/coarse resolutions for the self-convergence legs. The
# Sieverts leg (IC discontinuous at both faces) needs dt <= ~0.05 s for
# the early boundary layers; the cL=0 legs stay converged at dt=0.1 s.
CELLS_FINE = 400
CELLS_COARSE = 200
DT_PER_LEG = {
    "F1": (0.1, 0.2),
    "F2": (0.1, 0.2),
    "F3": (0.025, 0.05),
    "F4": (0.1, 0.2),
}

# FESTIM fine/coarse mesh sizes (vertices) and fixed stepsizes [s]. The early
# breakthrough front (first gated output at t ~= 42 s, J/J_ss ~= 1.5%) is
# step-dominated: dt = 1 s over-smooths it by ~34% vs the analytic series
# while dt = 0.1 s errs by ~3.4%; fine/coarse self-change at 0.025/0.05 s is
# 6.6e-3 (< 1/3 of the 3% J(t) gate). The F4 interface front is worse: its
# Nitsche-coupled interface needs half the step again (self-change
# 1.05e-2 -> 1.4e-3 at 0.0125/0.025 s), so F4 overrides the pair. Calibrated
# from the first probe runs; gates unchanged.
MESH_FINE = 200
MESH_COARSE = 100
DT_FESTIM_FINE = 0.025
DT_FESTIM_COARSE = 0.05
DT_FESTIM_F4 = (0.0125, 0.025)


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


class ProbeSkip(Exception):
    """Raised when a probe leg cannot run (reason recorded, never a pass)."""


def _analytical_breakthrough(t: float) -> float:
    """VV Val-1da analytic outlet flux ratio ``J(L,t)/J_ss`` (trap-free)."""
    if t <= 0.0:
        return 0.0
    total = 1.0
    n = 1
    while n <= 5000:
        term = 2.0 * (-1.0) ** n * math.exp(-n * n * math.pi**2 * D * t / (L * L))
        total += term
        if n > 5 and abs(term) < 1e-15:
            break
        n += 1
    return total


def _tlag_from_flux(times: list[float], flux: list[float], jss: float) -> tuple[float | None, str]:
    """Permeation lag from the cumulative-flux asymptote intercept.

    Q(t) = int_0^t J dt asymptotes to J_ss*(t - t_lag); the lag is the
    negative Q-intercept of a line fit over the asymptotic regime where
    the flux sits above 90% of J_ss (falling back to the 70-95% window).
    Returns (t_lag, note); t_lag is None when the window has too few
    points (recorded, not fit).
    """
    if jss <= 0.0 or len(times) < 2:
        return None, "no positive J_ss"
    ts = [0.0, *times]
    js = [0.0, *flux]
    q = [0.0]
    for i in range(1, len(ts)):
        q.append(q[-1] + 0.5 * (js[i - 1] + js[i]) * (ts[i] - ts[i - 1]))
    window = [(t, qv) for t, j, qv in zip(ts, js, q, strict=True) if j >= 0.90 * jss]
    if len(window) < 3:
        window = [
            (t, qv) for t, j, qv in zip(ts, js, q, strict=True) if 0.70 * jss <= j <= 0.95 * jss
        ]
    if len(window) < 2:
        return None, f"only {len(window)} point(s) in the asymptotic window"
    tw = np.array([p[0] for p in window])
    qw = np.array([p[1] for p in window])
    slope, intercept = np.polyfit(tw, qw, deg=1)
    if slope <= 0.0:
        return None, "non-positive asymptote slope"
    return float(-intercept / slope), f"fit over {len(window)} points"


def _max_norm_profile(
    ref_x: list[float],
    ref_c: list[float],
    other_x: list[float],
    other_c: list[float],
    splits: list[float] | None = None,
) -> float:
    """Max-norm profile mismatch, normalised to the inlet concentration.

    Comparison is restricted to the overlap of the two grids (no endpoint
    clamping) and, per ``splits``, to each continuous segment (never
    interpolated across a solubility jump at an interface). Reference points
    sitting exactly on a split are skipped: there the field is
    double-valued by construction, so its single-sided limits are compared
    at the neighbouring nodes instead.
    """
    if len(ref_x) < 2 or len(other_x) < 2:
        return float("nan")
    split_points = set(splits or [])
    bounds = [float(ref_x[0]), *(splits or []), float(ref_x[-1])]
    worst = 0.0
    saw_point = False
    for lo, hi in zip(bounds, bounds[1:], strict=False):
        ref_seg = [(x, c) for x, c in zip(ref_x, ref_c, strict=True) if lo <= x <= hi]
        other_seg = [(x, c) for x, c in zip(other_x, other_c, strict=True) if lo <= x <= hi]
        if len(ref_seg) < 1 or len(other_seg) < 2:
            continue
        saw_point = True
        o_x = np.asarray([p[0] for p in other_seg])
        o_c = np.asarray([p[1] for p in other_seg])
        for x, c in ref_seg:
            if x in split_points:
                continue
            if o_x[0] <= x <= o_x[-1]:
                other_at = float(np.interp(x, o_x, o_c))
                worst = max(worst, abs(float(c) - other_at) / C0)
    return worst if saw_point else float("nan")


def _gated_j_comparison(
    times: list[float],
    ref_j: list[float],
    other_j: list[float],
    jss_ref: float,
    tlag_ref: float,
) -> tuple[int, float]:
    """Return (n_gated, worst rel err) under the J(t) gate mask."""
    if jss_ref <= 0.0:
        return 0, float("nan")
    floor = 1.0e-3 * jss_ref
    worst = 0.0
    n_gated = 0
    for t, rj, oj in zip(times, ref_j, other_j, strict=True):
        if t < 0.2 * tlag_ref:
            continue
        if rj < floor:
            continue
        n_gated += 1
        worst = max(worst, abs(rj - oj) / max(abs(rj), floor))
    return n_gated, worst


def _val1da_nucleide() -> tuple[list[list[str]], list[str]]:
    """Always-run VV Val-1da leg: Nucleide vs the analytic series."""
    rows: list[list[str]] = []
    notes: list[str] = []
    ref = [_analytical_breakthrough(t) for t in G2_TIMES]
    got = tri.breakthrough(D, L, G2_TIMES)
    # Gate on the ratio series directly (J/J_ss weights cancel J_ss=1e-6).
    worst = max(
        (abs(g - r) / max(abs(r), 1.0e-3))
        for g, r, t in zip(got, ref, G2_TIMES, strict=True)
        if t >= 0.2 * TLAG and r >= 1.0e-3
    )
    n_gated = sum(1 for r, t in zip(ref, G2_TIMES, strict=True) if t >= 0.2 * TLAG and r >= 1.0e-3)
    rows.append(
        ["V1 J/J_ss vs Val-1da series", fmt(worst), "< 3%", _check(worst < 0.03, "V1 J(t)")]
    )
    jss = D * C0 / L
    steady = tri.steady(
        L, 400, D, {"kind": "dirichlet", "value": C0}, {"kind": "dirichlet", "value": CL}
    )
    err_jss = rel_diff(steady["flux_right"], jss)
    rows.append(["V1 J_ss vs D*c0/L", fmt(err_jss), "< 2%", _check(err_jss < 0.02, "V1 Jss")])
    err_lag = rel_diff(tri.time_lag(L, D), TLAG)
    rows.append(["V1 t_lag vs L^2/6D", fmt(err_lag), "< 5%", _check(err_lag < 0.05, "V1 tlag")])
    notes.append(
        f"V1 VV Val-1da (analytic permeation, {n_gated} gated of {len(G2_TIMES)} outputs): "
        f"worst J/J_ss deviation {worst:.3e} on the published series definition."
    )
    return rows, notes


def _mms_nucleide() -> tuple[list[list[str]], list[str]]:
    """Always-run M1: manufactured-solution order on the steady operator."""
    rows: list[list[str]] = []
    notes: list[str] = []
    amp_b = 0.2
    wave = math.pi / L
    errors: list[float] = []
    grids = [50, 100, 200, 400]
    for cells in grids:
        dx = L / cells
        source = [D * amp_b * wave * wave * math.sin(wave * ((i + 0.5) * dx)) for i in range(cells)]
        out = tri.steady(
            L,
            cells,
            D,
            {"kind": "dirichlet", "value": 1.0},
            {"kind": "dirichlet", "value": 0.0},
            source=source,
        )
        err = max(
            abs(c - ((1.0 - x / L) + amp_b * math.sin(wave * x)))
            for x, c in zip(out["centres"], out["mobile"], strict=True)
        )
        errors.append(err)
    orders = [
        math.log2(errors[i] / errors[i + 1]) if errors[i + 1] > 0.0 else float("nan")
        for i in range(len(errors) - 1)
    ]
    finest_order = orders[-1]
    rows.append(
        [
            "M1 MMS order (200->400)",
            f"{finest_order:.3f}",
            "1.5-2.5",
            _check(1.5 <= finest_order <= 2.5, "M1 order"),
        ]
    )
    rows.append(
        [
            "M1 MMS max error (400 cells)",
            fmt(errors[-1]),
            "< 5e-6",
            _check(errors[-1] < 5.0e-6, "M1 error"),
        ]
    )
    notes.append(
        f"M1 manufactured steady solution c* = (1-x/L) + {amp_b}*sin(pi x/L) with "
        f"source S = D*B*(pi/L)^2*sin(pi x/L); max-norm errors "
        + ", ".join(f"{e:.2e}" for e in errors)
        + f" on {grids} cells (orders {', '.join(f'{o:.2f}' for o in orders)})."
    )
    return rows, notes


def _layer_specs(cells: int) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for layer in LAYERS_2:
        specs.append({**layer, "cells": cells})
    return specs


def _nucleide_leg(leg: str, cells: int, dt_max: float) -> dict[str, Any]:
    """One Nucleide solve pair (steady + transient) for a numbered leg."""
    right: dict[str, Any]
    out: dict[str, Any]
    if leg == "F1":
        left = {"kind": "dirichlet", "value": C0}
        right = {"kind": "dirichlet", "value": CL}
        times = LEG_TIMES["F1"]
        traps = None
        steady = tri.steady(L, cells, D, left, right)
        transient = tri.transient(
            L, cells, D, left, right, times, dt_max=dt_max, rtol=1e-10, atol=1e-14
        )
        out = {
            "times": times,
            "x": steady["centres"],
            "c": steady["mobile"],
            "jss": steady["flux_right"],
            "flux_t": transient["flux_right"],
        }
    elif leg == "F2":
        left = {"kind": "dirichlet", "value": C0}
        right = {"kind": "dirichlet", "value": CL}
        times = LEG_TIMES["F2"]
        traps = [{"k0": TRAP_K0, "p0": TRAP_P0, "site_density": TRAP_N}]
        steady = tri.steady(L, cells, D, left, right, traps=traps)
        transient = tri.transient(
            L, cells, D, left, right, times, traps=traps, dt_max=dt_max, rtol=1e-10, atol=1e-14
        )
        out = {
            "times": times,
            "x": steady["centres"],
            "c": steady["mobile"],
            "jss": steady["flux_right"],
            "flux_t": transient["flux_right"],
        }
    elif leg == "F3":
        left = {"kind": "sieverts", "solubility": SIEVERT_K, "pressure": P_IN}
        right = {"kind": "sieverts", "solubility": SIEVERT_K, "pressure": P_OUT}
        times = LEG_TIMES["F3"]
        steady = tri.steady(L, cells, D, left, right)
        transient = tri.transient(
            L, cells, D, left, right, times, dt_max=dt_max, rtol=1e-10, atol=1e-14
        )
        out = {
            "times": times,
            "x": steady["centres"],
            "c": steady["mobile"],
            "jss": steady["flux_right"],
            "flux_t": transient["flux_right"],
        }
    elif leg == "F4":
        left = {"kind": "dirichlet", "value": C0}
        right = {"kind": "dirichlet", "value": CL}
        times = LEG_TIMES["F4"]
        layers = _layer_specs(cells // 2)
        steady = tri.steady_layers(layers, left, right)
        transient = tri.transient_layers(
            layers, left, right, times, dt_max=dt_max, rtol=1e-10, atol=1e-14
        )
        out = {
            "times": times,
            "x": steady["centres"],
            "c": steady["mobile"],
            "jss": steady["flux_right"],
            "flux_t": transient["flux_right"],
        }
    else:
        raise ValueError(f"unknown leg {leg}")
    tlag, note = _tlag_from_flux(out["times"], out["flux_t"], out["jss"])
    out["tlag"] = tlag
    out["tlag_note"] = note
    return out


def _self_convergence_row(leg: str, quantity: str, change: float, gate: float) -> list[str]:
    return [
        f"{leg} {quantity} (Nucleide fine-coarse)",
        fmt(change),
        f"< {gate / 3.0:.2e}",
        _check(change < gate / 3.0, f"{leg} self-conv {quantity}"),
    ]


def _nucleide_self_convergence() -> tuple[list[list[str]], list[str], dict[str, dict]]:
    """Run Nucleide fine/coarse for F1-F4; returns (rows, notes, fine legs)."""
    rows: list[list[str]] = []
    notes: list[str] = []
    fine_legs: dict[str, dict] = {}
    for leg in ("F1", "F2", "F3", "F4"):
        dt_fine, dt_coarse = DT_PER_LEG[leg]
        fine = _nucleide_leg(leg, CELLS_FINE, dt_fine)
        coarse = _nucleide_leg(leg, CELLS_COARSE, dt_coarse)
        fine_legs[leg] = fine
        d_jss = rel_diff(fine["jss"], coarse["jss"])
        rows.append(_self_convergence_row(leg, "J_ss", d_jss, 0.02))
        tlag_ref = (
            fine["tlag"] if fine["tlag"] is not None else (TLAG_TRAP if leg == "F2" else TLAG)
        )
        _, d_jt = _gated_j_comparison(
            fine["times"], fine["flux_t"], coarse["flux_t"], fine["jss"], tlag_ref
        )
        rows.append(_self_convergence_row(leg, "J(t)", d_jt, 0.03))
        if fine["tlag"] is not None and coarse["tlag"] is not None:
            d_lag = rel_diff(fine["tlag"], coarse["tlag"])
        else:
            d_lag = 0.0
        rows.append(_self_convergence_row(leg, "t_lag", d_lag, 0.05))
        d_c = _max_norm_profile(
            fine["x"], fine["c"], coarse["x"], coarse["c"], splits=[5.0e-4] if leg == "F4" else None
        )
        rows.append(_self_convergence_row(leg, "c(x)", d_c, 0.02))
        notes.append(
            f"{leg} Nucleide reference: J_ss={fine['jss']:.6e} mol/m^2/s, "
            f"t_lag={fine['tlag']} ({fine['tlag_note']})."
        )
    return rows, notes, fine_legs


def _require_festim() -> Any:
    """Import FESTIM 2.1 or raise ProbeSkip with the loud reason."""
    try:
        import festim
    except Exception as exc:  # noqa: BLE001 — absence is a recorded skip, not a crash
        raise ProbeSkip(f"FESTIM unavailable: {exc}") from exc
    version = getattr(festim, "__version__", None)
    if version is not None and not str(version).startswith("2.1"):
        raise ProbeSkip(f"FESTIM version {version} != 2.1 (probe is pinned to 2.1)")
    return festim


def _festim_run(festim: Any, model: Any, flux_export: Any, profile_exports: list) -> dict:
    """Initialise and run one FESTIM model; extract flux/profile numbers."""
    model.show_progress_bar = False
    model.initialise()
    model.run()
    n = len(flux_export.data)
    if n < 2:
        raise ProbeSkip("FESTIM flux export has fewer than 2 timesteps")
    if len(model.timesteps) >= 2:
        m = min(len(model.timesteps), n)
        times = [float(t) for t in model.timesteps[:m]]
    else:
        # HydrogenTransportProblemDiscontinuous (FESTIM 2.1) does not
        # populate ``model.timesteps`` on the discontinuous path; the export
        # is sampled every fixed step, so the probe's uniform-stepsize
        # cadence gives the export times exactly.
        dt_fixed = float(model.settings.stepsize.initial_value)
        times = [i * dt_fixed for i in range(n)]
        m = n
    # Signed flux: FESTIM's SurfaceFlux export at the right face is already
    # positive in the +x direction for F1/F2/F4. abs() was destructive for
    # F3, whose physical flux is negative for ~120 s (the interior starts at
    # zero, below the 2 mol/m^3 outlet Sieverts value, and fills up from the
    # left); flipping its sign corrupts the cumulative-flux t_lag fit.
    flux = [float(v) for v in flux_export.data[:m]]
    xs: list[float] = []
    cs: list[float] = []
    for profile in profile_exports:
        if profile.x is None or not profile.data:
            raise ProbeSkip("FESTIM profile export is empty")
        for xv, cv in zip(profile.x, profile.data[-1], strict=True):
            xs.append(float(xv))
            cs.append(float(cv))
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in order]
    cs = [cs[i] for i in order]
    return {"times": times, "flux": flux, "x": xs, "c": cs}


def _festim_slab(
    festim: Any, vertices: Any, leg: str, times_grid: list[float], final_time: float, dt: float
) -> dict:
    """Build, run, and reduce one single-material FESTIM leg."""
    model = festim.HydrogenTransportProblem()
    model.mesh = festim.Mesh1D(vertices=vertices)
    material = festim.Material(D_0=D, E_D=0.0)
    volume = festim.VolumeSubdomain1D(id=1, borders=[0.0, L], material=material)
    left = festim.SurfaceSubdomain1D(id=1, x=0.0)
    right = festim.SurfaceSubdomain1D(id=2, x=L)
    model.subdomains = [volume, left, right]
    tritium = festim.Species("T")
    model.species = [tritium]
    model.temperature = 500.0
    if leg == "F3":
        # FESTIM 2.1 SievertsBC asserts its value function returns a UFL
        # expression; with all-float S_0/E_S/pressure the isothermal
        # (E_S = 0) law collapses to a plain float ("<class 'float'>" assert).
        # The documented spelling is S_0 as float-or-fem.Constant
        # (SievertsBC docstring), so pass the synthetic constant as a
        # mesh-bound dolfinx Constant. Physics unchanged.
        from dolfinx import fem

        sievert_k = fem.Constant(model.mesh.mesh, SIEVERT_K)
        model.boundary_conditions = [
            festim.SievertsBC(
                subdomain=left, S_0=sievert_k, E_S=0.0, pressure=P_IN, species=tritium
            ),
            festim.SievertsBC(
                subdomain=right, S_0=sievert_k, E_S=0.0, pressure=P_OUT, species=tritium
            ),
        ]
    else:
        model.boundary_conditions = [
            festim.FixedConcentrationBC(subdomain=left, value=C0, species=tritium),
            festim.FixedConcentrationBC(subdomain=right, value=CL, species=tritium),
        ]
    if leg == "F2":
        empty = festim.Species("empty_traps", mobile=False)
        trapped = festim.Species("trapped_T", mobile=False)
        model.species = [tritium, empty, trapped]
        model.reactions = [
            festim.Reaction(
                reactant=[tritium, empty],
                product=[trapped],
                k_0=TRAP_K0,
                E_k=0.0,
                p_0=TRAP_P0,
                E_p=0.0,
                volume=volume,
            )
        ]
        model.initial_conditions = [
            festim.InitialConcentration(value=TRAP_N, species=empty, volume=volume)
        ]
    stepsize = festim.Stepsize(initial_value=dt, max_stepsize=dt)
    model.settings = festim.Settings(
        atol=1e-10, rtol=1e-10, transient=True, final_time=final_time, stepsize=stepsize
    )
    flux_export = festim.SurfaceFlux(field=tritium, surface=right)
    profile_export = festim.Profile1DExport(field=tritium, subdomain=volume)
    model.exports = [flux_export, profile_export]
    result = _festim_run(festim, model, flux_export, [profile_export])
    flux_on_grid = [
        float(np.interp(t, np.asarray(result["times"]), np.asarray(result["flux"])))
        for t in times_grid
    ]
    jss = float(result["flux"][-1])
    tlag, note = _tlag_from_flux(times_grid, flux_on_grid, jss)
    result = {
        "times": times_grid,
        "flux_t": flux_on_grid,
        "x": result["x"],
        "c": result["c"],
        "jss": jss,
        "tlag": tlag,
        "tlag_note": note,
    }
    return result


_scifem_shim_installed = False


def _compat_scifem_transfer_shim() -> None:
    """Adapt FESTIM 2.1's meshtag transfer call to the scifem 0.27 signature.

    FESTIM 2.1 calls ``scifem.mesh.transfer_meshtags_to_submesh(tag, submesh,
    v_map, cell_map)`` positionally; scifem 0.27 made ``vertex_to_parent`` /
    ``cell_to_parent`` keyword-only, so the discontinuous-leg initialise
    raises ``TypeError: takes 2 positional arguments but 4 were given``.
    Rebind the module-global festim imports to the keyword form of the same
    call (identical arguments, semantics); no-op when the positional form
    still exists (detected by signature). Only the F4 leg hits this path.
    """
    global _scifem_shim_installed
    if _scifem_shim_installed:
        return
    _scifem_shim_installed = True
    try:
        import festim.subdomain.volume_subdomain as vsd

        original = vsd.transfer_meshtags_to_submesh
        keyword_only = {
            p.name
            for p in inspect.signature(original).parameters.values()
            if p.kind == inspect.Parameter.KEYWORD_ONLY
        }
        if "vertex_to_parent" not in keyword_only:
            return

        def _keyword_transfer(entity_tag, submesh, vertex_to_parent=None, cell_to_parent=None):
            return original(
                entity_tag,
                submesh,
                vertex_to_parent=vertex_to_parent,
                cell_to_parent=cell_to_parent,
            )

        vsd.transfer_meshtags_to_submesh = _keyword_transfer
    except Exception:  # noqa: BLE001 — shim failure surfaces as the leg's own error
        pass


def _festim_two_layer(festim: Any, n_per_layer: int, times_grid: list[float], dt: float) -> dict:
    """Two-layer discontinuous FESTIM leg (Sieverts interface law).

    Species carry the new 2.1 ``subdomains`` attribute (the discontinuous
    problem rejects a species without it), and the internal interface is
    represented by ``Interface`` alone: a separate mid ``SurfaceSubdomain1D``
    in ``subdomains`` is connected to both volumes and is rejected.
    """
    _compat_scifem_transfer_shim()
    l1, l2 = 5e-4, 5e-4
    total = l1 + l2
    verts = np.concatenate(
        [
            np.linspace(0.0, l1, n_per_layer),
            np.linspace(l1, total, n_per_layer)[1:],
        ]
    )
    mat1 = festim.Material(D_0=1e-9, E_D=0.0, K_S_0=2.0, E_K_S=0.0, solubility_law="sievert")
    mat2 = festim.Material(D_0=5e-10, E_D=0.0, K_S_0=0.5, E_K_S=0.0, solubility_law="sievert")
    vol1 = festim.VolumeSubdomain1D(id=1, borders=[0.0, l1], material=mat1)
    vol2 = festim.VolumeSubdomain1D(id=2, borders=[l1, total], material=mat2)
    left = festim.SurfaceSubdomain1D(id=1, x=0.0)
    right = festim.SurfaceSubdomain1D(id=2, x=total)
    interface = festim.Interface(id=3, subdomains=[vol1, vol2])
    try:
        model = festim.HydrogenTransportProblemDiscontinuous(interfaces=[interface])
    except TypeError:
        model = festim.HydrogenTransportProblemDiscontinuous()
        model.interfaces = [interface]
    model.mesh = festim.Mesh1D(vertices=verts)
    model.subdomains = [vol1, vol2, left, right]
    tritium = festim.Species("T", subdomains=[vol1, vol2])
    model.species = [tritium]
    model.temperature = 500.0
    model.boundary_conditions = [
        festim.FixedConcentrationBC(subdomain=left, value=C0, species=tritium),
        festim.FixedConcentrationBC(subdomain=right, value=CL, species=tritium),
    ]
    stepsize = festim.Stepsize(initial_value=dt, max_stepsize=dt)
    model.settings = festim.Settings(
        atol=1e-10, rtol=1e-10, transient=True, final_time=1000.0, stepsize=stepsize
    )
    flux_export = festim.SurfaceFlux(field=tritium, surface=right)
    prof1 = festim.Profile1DExport(field=tritium, subdomain=vol1)
    prof2 = festim.Profile1DExport(field=tritium, subdomain=vol2)
    model.exports = [flux_export, prof1, prof2]
    result = _festim_run(festim, model, flux_export, [prof1, prof2])
    flux_on_grid = [
        float(np.interp(t, np.asarray(result["times"]), np.asarray(result["flux"])))
        for t in times_grid
    ]
    jss = float(result["flux"][-1])
    tlag, note = _tlag_from_flux(times_grid, flux_on_grid, jss)
    return {
        "times": times_grid,
        "flux_t": flux_on_grid,
        "x": result["x"],
        "c": result["c"],
        "jss": jss,
        "tlag": tlag,
        "tlag_note": note,
    }


def _festim_leg(festim: Any, leg: str, n_nodes: int, dt: float) -> dict:
    """Run one FESTIM leg at one resolution."""
    if leg == "F4":
        return _festim_two_layer(festim, n_nodes // 2, LEG_TIMES["F4"], dt)
    if leg == "F2":
        return _festim_slab(
            festim, np.linspace(0.0, L, n_nodes), leg, LEG_TIMES["F2"], T_END_TRAP, dt
        )
    return _festim_slab(festim, np.linspace(0.0, L, n_nodes), leg, LEG_TIMES[leg], 1000.0, dt)


def _cross_gate_rows(
    leg: str, fine_nucleide: dict, festim_fine: dict
) -> tuple[list[list[str]], list[str]]:
    """Cross-code gates for one leg plus its FESTIM self-convergence note."""
    rows: list[list[str]] = []
    notes: list[str] = []
    d_jss = rel_diff(fine_nucleide["jss"], festim_fine["jss"])
    rows.append([f"{leg} J_ss parity", fmt(d_jss), "< 2%", _check(d_jss < 0.02, f"{leg} Jss")])
    tlag_ref = fine_nucleide["tlag"] if fine_nucleide["tlag"] is not None else TLAG
    _, d_jt = _gated_j_comparison(
        fine_nucleide["times"],
        fine_nucleide["flux_t"],
        festim_fine["flux_t"],
        fine_nucleide["jss"],
        tlag_ref,
    )
    rows.append([f"{leg} J(t) parity", fmt(d_jt), "< 3%", _check(d_jt < 0.03, f"{leg} J(t)")])
    if fine_nucleide["tlag"] is not None and festim_fine["tlag"] is not None:
        d_lag = rel_diff(fine_nucleide["tlag"], festim_fine["tlag"])
        status = _check(d_lag < 0.05, f"{leg} tlag")
    else:
        d_lag = float("nan")
        status = "SKIP"
        notes.append(
            f"{leg} t_lag parity SKIPPED: "
            f"Nucleide t_lag={fine_nucleide['tlag']} "
            f"({fine_nucleide['tlag_note']}), "
            f"FESTIM t_lag={festim_fine['tlag']} ({festim_fine['tlag_note']})."
        )
    rows.append([f"{leg} t_lag parity", f"{d_lag:.6e}", "< 5%", status])
    d_c = _max_norm_profile(
        fine_nucleide["x"],
        fine_nucleide["c"],
        festim_fine["x"],
        festim_fine["c"],
        splits=[5.0e-4] if leg == "F4" else None,
    )
    rows.append(
        [f"{leg} c(x) parity (max-norm)", fmt(d_c), "< 2%", _check(d_c < 0.02, f"{leg} c(x)")]
    )
    notes.append(
        f"{leg} parity: Nucleide J_ss={fine_nucleide['jss']:.6e} vs "
        f"FESTIM J_ss={festim_fine['jss']:.6e}; "
        f"FESTIM t_lag={festim_fine['tlag']} ({festim_fine['tlag_note']})."
    )
    return rows, notes


def probe_legs(fine_legs: dict[str, dict]) -> tuple[list[list[str]], list[str], list[list[str]]]:
    """Run the FESTIM probe legs; SKIP loudly when FESTIM is unavailable.

    Returns (cross-gate rows, prose notes, FESTIM self-convergence rows).
    """
    rows: list[list[str]] = []
    notes: list[str] = []
    self_rows: list[list[str]] = []
    try:
        festim = _require_festim()
    except ProbeSkip as exc:
        note = f"FESTIM probe legs SKIPPED: {exc}"
        print(note)
        for leg in ("F1", "F2", "F3", "F4"):
            rows.append([f"{leg} probe legs", "FESTIM unavailable", "SKIP"])
        notes.append(note)
        return rows, notes, self_rows
    version = getattr(festim, "__version__", "unknown")
    notes.append(f"FESTIM probe running: version {version} (pinned 2.1, black-box).")
    for leg in ("F1", "F2", "F3", "F4"):
        dt_fine, dt_coarse = DT_FESTIM_F4 if leg == "F4" else (DT_FESTIM_FINE, DT_FESTIM_COARSE)
        try:
            fine_festim = _festim_leg(festim, leg, MESH_FINE, dt_fine)
            coarse_festim = _festim_leg(festim, leg, MESH_COARSE, dt_coarse)
        except ProbeSkip as exc:
            note = f"{leg} probe leg SKIPPED: {exc}"
            print(note)
            notes.append(note)
            rows.append([f"{leg} probe leg", "SKIP", "SKIP"])
            continue
        except Exception as exc:  # noqa: BLE001 — probe failure is recorded, never a pass
            note = f"{leg} probe leg SKIPPED: FESTIM run failed: {exc}"
            print(note)
            notes.append(note)
            rows.append([f"{leg} probe leg", "SKIP", "SKIP"])
            continue
        d_jss = rel_diff(fine_festim["jss"], coarse_festim["jss"])
        self_rows.append(
            [
                f"{leg} J_ss (FESTIM fine-coarse)",
                fmt(d_jss),
                "< 6.67e-3",
                _check(d_jss < 0.02 / 3.0, f"{leg} festim self-conv Jss"),
            ]
        )
        tlag_ref = fine_festim["tlag"] if fine_festim["tlag"] is not None else TLAG
        _, d_jt = _gated_j_comparison(
            fine_festim["times"],
            fine_festim["flux_t"],
            coarse_festim["flux_t"],
            fine_festim["jss"],
            tlag_ref,
        )
        self_rows.append(
            [
                f"{leg} J(t) (FESTIM fine-coarse)",
                fmt(d_jt),
                "< 1.00e-2",
                _check(d_jt < 0.03 / 3.0, f"{leg} festim self-conv Jt"),
            ]
        )
        if fine_festim["tlag"] is not None and coarse_festim["tlag"] is not None:
            d_lag = rel_diff(fine_festim["tlag"], coarse_festim["tlag"])
        else:
            d_lag = 0.0
        self_rows.append(
            [
                f"{leg} t_lag (FESTIM fine-coarse)",
                fmt(d_lag),
                "< 1.67e-2",
                _check(d_lag < 0.05 / 3.0, f"{leg} festim self-conv tlag"),
            ]
        )
        d_c = _max_norm_profile(
            fine_festim["x"],
            fine_festim["c"],
            coarse_festim["x"],
            coarse_festim["c"],
            splits=[5.0e-4] if leg == "F4" else None,
        )
        self_rows.append(
            [
                f"{leg} c(x) (FESTIM fine-coarse)",
                fmt(d_c),
                "< 6.67e-3",
                _check(d_c < 0.02 / 3.0, f"{leg} festim self-conv c"),
            ]
        )
        leg_rows, leg_notes = _cross_gate_rows(leg, fine_legs[leg], fine_festim)
        rows.extend(leg_rows)
        notes.extend(leg_notes)
    return rows, notes, self_rows


def main() -> int:
    report = Report(
        "tritium_festim", "Tritium 1D kernel vs FESTIM 2.1 probe (`tritium_vs_festim.py`)"
    )
    report.prose(
        "Black-box parity probe between the `nucleide.tritium` 1D "
        "diffusion-trapping kernel and FESTIM 2.1 (Apache-2.0, "
        "`festim==2.1`, v2.x API), driven in the probe-only "
        "`dolfinx/dolfinx:stable` image (docker; never the main validation "
        "container, never the Rust workspace). Progression F1 trap-free "
        "Dirichlet-Dirichlet, F2 one McNabb-Foster trap, F3 Sieverts ends, "
        "F4 two-layer stack; recombination ends are excluded. Always-run "
        "legs VV Val-1da (analytic permeation) and M1 (manufactured-solution "
        "order) gate the Nucleide side everywhere; every FESTIM leg SKIPs "
        "loudly with its reason outside the probe image."
    )
    val_rows, val_notes = _val1da_nucleide()
    report.table(["Gate", "Value", "Expected", "Status"], val_rows)
    for note in val_notes:
        report.prose(note)
    mms_rows, mms_notes = _mms_nucleide()
    report.table(["Gate", "Value", "Expected", "Status"], mms_rows)
    for note in mms_notes:
        report.prose(note)

    self_rows, self_notes, fine_legs = _nucleide_self_convergence()
    probe_rows, probe_notes, festim_self_rows = probe_legs(fine_legs)
    self_rows.extend(festim_self_rows)
    report.table(["Self-convergence (< gate/3)", "Change", "Threshold", "Status"], self_rows)
    for note in self_notes:
        report.prose(note)
    report.table(["Cross-code gate", "Value", "Gate", "Status"], probe_rows)
    for note in probe_notes:
        report.prose(note)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} tritium/FESTIM probe check(s) failed", file=sys.stderr)
        return 1
    print("All tritium/FESTIM probe checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
