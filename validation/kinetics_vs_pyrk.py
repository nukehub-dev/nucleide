"""Point-kinetics cross-check (`nucleide-kinetics` vs analytic gates + PyRK).

Two tiers:

1. Analytic gates (always run): O1 initial-rate identity, O2 prompt-jump
   plateau, O3 1-group closed form, O4 6-group stable-period tail, plus
   invariants (zero-reactivity stationarity, positivity). Inputs come from
   the hand-built synthetic `fixtures/kinetics/` files.
2. PyRK cross-check (O5): the same ramp transient solved with the upstream
   PyRK neutronics block (`dpdt`/`dzetadt`, BSD-3) driven by its own
   `dopri5` loop, against `nucleide.kinetics` on identical runtime-read
   precursor data. PyRK is an optional oracle dependency: if it cannot be
   imported, tier 2 is reported as SKIP with its reason (never silently).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.kinetics as kin

FAILURES = 0
REPO_ROOT = Path(__file__).resolve().parent.parent
FIX = REPO_ROOT / "fixtures" / "kinetics"

# Six-group synthetic set shared by the ramp-table fixture test (hand-picked
# round numbers, matching fixtures/kinetics/README.md provenance).
SYN_BETAS = [0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032]
SYN_LAMBDAS = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
SYN_LAMBDA = 1e-5


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _load(name: str) -> dict:
    return json.loads((FIX / name).read_text())


def _ramp_table() -> tuple[list[float], list[float]]:
    times: list[float] = []
    values: list[float] = []
    for k, line in enumerate((FIX / "ramp_table.csv").read_text().splitlines()):
        if k == 0:
            continue
        a, b = line.split(",")
        times.append(float(a))
        values.append(float(b))
    return times, values


def tier1() -> tuple[list[list[str]], list[str], list[list[str]], str]:
    """Analytic gates O1-O4 + invariants.

    Returns (gate rows, prose notes, figure-series rows, prompt-jump level).
    The figure rows tabulate the O3 1-group step transient (``t since step``,
    ``Nucleide n``, ``Analytic n``) for ``make_figures.py``; no gates live
    in those tables.
    """
    rows: list[list[str]] = []
    notes: list[str] = []
    step = _load("step_oracle.json")
    p1 = step["params"]
    s1 = step["step"]

    rate = kin.initial_rate(
        p1["betas"],
        p1["lambdas"],
        p1["Lambda"],
        {"kind": "constant", "rho": s1["rho_final"]},
        step["n0"],
    )
    err = rel_diff(rate, step["initial_rate"])
    notes.append(f"O1 initial rate: got {rate:.6e}, want {step['initial_rate']:.6e}.")
    rows.append(["O1 initial rate", fmt(err), "< 1e-12", _check(err < 1e-12, "O1 initial rate")])

    pj = kin.prompt_jump(1.0, s1["rho_init"], s1["rho_final"], sum(p1["betas"]))
    err = rel_diff(pj, step["prompt_jump_factor"])
    notes.append(f"O2 prompt-jump factor: got {pj:.6e}.")
    rows.append(
        ["O2 jump factor (algebraic)", fmt(err), "< 1e-12", _check(err < 1e-12, "O2 jump factor")]
    )
    out = kin.solve(
        SYN_BETAS,
        SYN_LAMBDAS,
        SYN_LAMBDA,
        {"kind": "step", "t_step": 1.0, "rho_init": 0.0, "rho_final": 0.003},
        [1.05],
        1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    plateau = kin.prompt_jump(1.0, 0.0, 0.003, sum(SYN_BETAS))
    err = rel_diff(out["n"][0], plateau)
    notes.append(f"O2 plateau: solver {out['n'][0]:.6e} vs PJA {plateau:.6e}.")
    rows.append(["O2 prompt plateau", fmt(err), "< 2e-2", _check(err < 2e-2, "O2 prompt plateau")])

    ts = [s1["t_step"] + t for t in step["analytic_since_step"]["t"]]
    out = kin.solve(
        p1["betas"],
        p1["lambdas"],
        p1["Lambda"],
        {"kind": "step", **s1},
        ts,
        step["n0"],
        rtol=1e-10,
        atol=1e-14,
    )
    want = step["analytic_since_step"]["n"]
    worst = max(rel_diff(n, w) for n, w in zip(out["n"], want, strict=True))
    # Figure-source rows use repr (round-trip precision): fmt's 6 decimals
    # would collapse the ~1e-8 Nucleide-vs-analytic gap the figure checks.
    series_rows = [
        [repr(t), repr(n), repr(w)]
        for t, n, w in zip(step["analytic_since_step"]["t"], out["n"], want, strict=True)
    ]
    pj_level = repr(step["n0"] * step["prompt_jump_factor"])
    notes.append(f"O3 1-group closed form: worst rel err {worst:.3e} over {len(ts)} nodes.")
    rows.append(
        ["O3 1-group transient", fmt(worst), "< 1e-6", _check(worst < 1e-6, "O3 1-group transient")]
    )

    check = _load("inhour_check.json")
    p6 = check["params"]
    period = kin.stable_period(p6["betas"], p6["lambdas"], p6["Lambda"], check["rho"])
    err = rel_diff(period, check["stable_period_s"])
    notes.append(f"O4 stable period: got {period:.6e} s.")
    rows.append(["O4 stable period", fmt(err), "< 1e-9", _check(err < 1e-9, "O4 stable period")])
    out = kin.solve(
        p6["betas"],
        p6["lambdas"],
        p6["Lambda"],
        {"kind": "constant", "rho": check["rho"]},
        [50.0, 100.0],
        1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    slope = math.log(out["n"][1] / out["n"][0]) / 50.0
    err = rel_diff(slope, 1.0 / period)
    notes.append(f"O4 tail slope: got {slope:.6e} 1/s vs inhour root {1.0 / period:.6e} 1/s.")
    rows.append(["O4 tail slope", fmt(err), "< 2e-2", _check(err < 2e-2, "O4 tail slope")])

    out = kin.solve(
        SYN_BETAS,
        SYN_LAMBDAS,
        SYN_LAMBDA,
        {"kind": "constant", "rho": 0.0},
        [1.0, 10.0, 100.0],
        2.5,
    )
    worst = max(abs(n - 2.5) / 2.5 for n in out["n"])
    rows.append(
        [
            "Invariant: rho=0 stationary",
            fmt(worst),
            "< 1e-9",
            _check(worst < 1e-9, "rho=0 stationary"),
        ]
    )
    times, values = _ramp_table()
    out = kin.solve(
        SYN_BETAS,
        SYN_LAMBDAS,
        SYN_LAMBDA,
        {"kind": "polyline", "times": times, "values": values},
        [t for t in times if t > 0],
        1.0,
    )
    ok = all(n > 0.0 for n in out["n"]) and out["n"][-1] > 1.0
    notes.append(f"Ramp-table transient: n(10 s) = {out['n'][-1]:.6e}, all positive.")
    rows.append(
        [
            "Invariant: ramp positivity/growth",
            f"n_end={out['n'][-1]:.6e}",
            "> 1",
            _check(ok, "ramp positivity"),
        ]
    )
    return rows, notes, series_rows, pj_level


def tier2_pyrk() -> tuple[list[list[str]], list[str], bool]:
    """PyRK cross-check (O5). Returns (rows, notes, skipped)."""
    try:
        import numpy as np
        from pyrk import neutronics
        from pyrk import reactivity_insertion as ri
        from pyrk.timer import Timer
        from pyrk.utilities.ur import units
        from scipy.integrate import ode
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Tier 2 (PyRK cross-check) SKIPPED: {exc}"
        print(note)
        return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []
    dt, tf = 0.005, 5.0
    ti = Timer(t0=0.0 * units.seconds, tf=tf * units.seconds, dt=dt * units.seconds)
    rho_ext = ri.RampReactivityInsertion(
        timer=ti,
        t_start=1.0 * units.seconds,
        t_end=2.0 * units.seconds,
        rho_init=0.0 * units.delta_k,
        rho_rise=0.003 * units.delta_k,
        rho_final=0.003 * units.delta_k,
    )
    ne = neutronics.Neutronics(
        iso="u235", e="thermal", n_precursors=6, n_decay=0, timer=ti, rho_ext=rho_ext
    )
    # Precursor data are read from the oracle package at runtime (never
    # vendored): identical tables feed both codes below.
    betas = list(ne._pd.betas())
    lambdas = list(ne._pd.lambdas())
    lam_gen = float(ne._pd.Lambda())
    notes.append(f"PyRK u235/thermal runtime data: beta={sum(betas):.6e}, Lambda={lam_gen:.3e} s.")
    n_pg = 6
    y0 = np.zeros(1 + n_pg)
    y0[0] = 1.0
    for j in range(n_pg):
        y0[1 + j] = y0[0] * betas[j] / (lambdas[j] * lam_gen)

    def f_n(t: float, y: np.ndarray, neu: object, timer: object) -> np.ndarray:
        assert isinstance(neu, neutronics.Neutronics)
        assert isinstance(timer, Timer)
        idx = timer.t_idx(t * units.seconds)
        f = np.zeros(1 + n_pg)
        f[0] = neu.dpdt(idx, [], y[0], y[1:])
        for j in range(n_pg):
            f[1 + j] = neu.dzetadt(t, y[0], y[1 + j], j)
        return f

    solver = ode(f_n).set_integrator("dopri5")
    solver.set_initial_value(y0, 0.0)
    solver.set_f_params(ne, ti)
    probe_idx = [int(2.0 / dt), int(3.5 / dt), int(tf / dt)]
    pyrk_n: dict[float, float] = {}
    for k in range(1, ti.timesteps()):
        t_end = ti.t(k).magnitude
        solver.integrate(t_end)
        if k in probe_idx:
            pyrk_n[t_end] = float(solver.y[0])
    notes.append("PyRK driver loop mirrors upstream driver.solve (dopri5, per-step).")

    out = kin.solve(
        betas,
        lambdas,
        lam_gen,
        {
            "kind": "ramp",
            "t_start": 1.0,
            "t_end": 2.0,
            "rho_init": 0.0,
            "rho_rise": 0.003,
            "rho_final": 0.003,
        },
        sorted(pyrk_n),
        1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    worst = 0.0
    for t, n in zip(out["times"], out["n"], strict=True):
        err = rel_diff(n, pyrk_n[t])
        worst = max(worst, err)
        rows.append(
            [
                f"O5 ramp t={t:.2f} s",
                fmt(pyrk_n[t]),
                fmt(n),
                fmt(err),
                _check(err < 1e-3, f"O5 ramp t={t}"),
            ]
        )
    notes.append(f"O5 ramp vs PyRK: worst rel err {worst:.3e} at {len(pyrk_n)} probes.")
    return rows, notes, False


def main() -> int:
    report = Report("kinetics", "Point kinetics (`kinetics_vs_pyrk.py`)")
    report.prose(
        "Two-tier oracle for `nucleide.kinetics`: analytic gates O1-O4 plus "
        "invariants on synthetic fixtures (tier 1, always run), and a ramp "
        "cross-check against the upstream PyRK neutronics block on "
        "runtime-read precursor data (tier 2 / O5)."
    )
    rows1, notes1, series_rows, pj_level = tier1()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Rel err", "Tol", "Status"], rows1)
    report.heading("Step-transient series (figure source)", level=3)
    report.table(["t since step (s)", "Nucleide n", "Analytic n"], series_rows)
    report.table(["Quantity", "Value"], [["Prompt-jump level (n0 = 1)", pj_level]])
    rows2, notes2, skipped = tier2_pyrk()
    for note in notes2:
        report.prose(note)
    if skipped:
        report.table(["Gate", "Status"], [["O5 ramp vs PyRK", "SKIP (pyrk unavailable)"]])
    else:
        report.table(["Gate", "PyRK n", "Nucleide n", "Rel err", "Status"], rows2)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} kinetics check(s) failed", file=sys.stderr)
        return 1
    print("All kinetics checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
