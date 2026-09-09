"""Compare Nucleide's CRAM depletion against OpenMC CRAM-48 and analytic Bateman."""

from __future__ import annotations

import importlib.metadata
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from common import Report, fmt

import nucleide

try:
    import openmc.deplete

    HAS_OPENMC = True
    OPENMC_SKIP = ""
except Exception as exc:
    openmc = None  # type: ignore[no-redef]
    HAS_OPENMC = False
    OPENMC_SKIP = f"OpenMC oracle skipped: cannot import openmc.deplete ({exc})."

REPO_ROOT = Path(__file__).resolve().parent.parent
CHAIN_NI = REPO_ROOT / "fixtures" / "depletion" / "chain_ni.xml"
CHAIN_SIMPLE = REPO_ROOT / "fixtures" / "depletion" / "chain_simple.xml"

# Targets for <decay> entries in chain_ni.xml that lack a `target` attribute.
# These must be added for Nucleide's parser; we add them identically to the
# patch used in crates/depletion/benches/depletion_bench.rs.
DECAY_DAUGHTERS = {
    "Fe55": "Mn55",
    "Fe59": "Co59",
    "Ni57": "Co57",
    "Ni59": "Co59",
    "Ni63": "Cu63",
    "Ni65": "Cu65",
}


def patched_chain_ni_path() -> str:
    """Return a path to a patched chain_ni.xml usable by both Nucleide and OpenMC."""
    text = CHAIN_NI.read_text()

    current = ""
    out_lines: list[str] = []
    for line in text.splitlines():
        start = line.find('<nuclide name="')
        if start != -1:
            sub = line[start + 15 :]
            end = sub.find('"')
            if end != -1:
                current = sub[:end]
        if "<decay" in line and "target=" not in line:
            daughter = DECAY_DAUGHTERS.get(current, "Nothing")
            line = line.replace("/>", f' target="{daughter}" />')
        out_lines.append(line)

    extra: list[str] = []
    seen: set[str] = set()
    for daughter in DECAY_DAUGHTERS.values():
        if daughter != "Nothing" and f'name="{daughter}"' not in text and daughter not in seen:
            extra.append(f'  <nuclide name="{daughter}" reactions="0"/>')
            seen.add(daughter)

    out = "\n".join(out_lines)
    if extra:
        out = out.replace("</depletion_chain>", "\n".join(extra) + "\n</depletion_chain>")

    with tempfile.NamedTemporaryFile("w", suffix=".xml", delete=False) as f:
        f.write(out)
        return f.name


def nucleide_rates(chain: nucleide.depletion.Chain) -> dict[str, float]:
    """Build one-group reaction rates for chain_ni.xml in Nucleide's format."""
    rates: dict[str, float] = {}
    for nuc in chain.nuclides:
        # Reaction kinds are discovered from the chain XML in OpenMC below.
        # We apply a uniform small rate so both solvers see the same matrix.
        rates[f"{nuc}:(n,gamma)"] = 1.0e-7
        rates[f"{nuc}:(n,2n)"] = 1.0e-8
        rates[f"{nuc}:(n,p)"] = 1.0e-8
        rates[f"{nuc}:(n,a)"] = 1.0e-8
    return rates


def openmc_rates(chain: openmc.deplete.Chain) -> openmc.deplete.ReactionRates:
    """Build matching one-group reaction rates for OpenMC's form_matrix."""
    reactions = sorted({r.type for nuc in chain.nuclides for r in nuc.reactions})
    nuc_names = [nuc.name for nuc in chain.nuclides]
    rates = openmc.deplete.ReactionRates(["0"], nuc_names, reactions)
    for nuc in chain.nuclides:
        i = nuc_names.index(nuc.name)
        for r in nuc.reactions:
            j = reactions.index(r.type)
            if r.type == "(n,gamma)":
                val = 1.0e-7
            elif r.type in ("(n,2n)", "(n,p)", "(n,a)"):
                val = 1.0e-8
            else:
                val = 0.0
            rates[0, i, j] = val
    return rates


def run_chain_ni_vectors() -> tuple[dict[str, float], dict[str, float]]:
    """Deplete chain_ni.xml with Nucleide and OpenMC; return final density dicts."""
    chain_path = patched_chain_ni_path()
    nuc_chain = nucleide.depletion.read_chain(chain_path)
    om_chain = openmc.deplete.Chain.from_xml(chain_path)

    dt = 2.592e6  # 30 days in seconds
    n0 = {nuc: (1.0e24 if nuc == "Ni58" else 0.0) for nuc in nuc_chain.nuclides}

    nuc_rates = nucleide_rates(nuc_chain)
    nuc_out = nucleide.depletion.deplete(nuc_chain, n0, dt, rates=nuc_rates, order=48)

    om_rates = openmc_rates(om_chain)
    A = om_chain.form_matrix(om_rates[0])
    n0_vec = np.array([n0.get(nuc.name, 0.0) for nuc in om_chain.nuclides])
    om_out_vec = openmc.deplete.cram.CRAM48(A, n0_vec, dt)
    om_out = {nuc.name: float(om_out_vec[i]) for i, nuc in enumerate(om_chain.nuclides)}
    return dict(nuc_out), om_out


def run_chain_ni_comparison() -> dict[str, float]:
    """Deplete chain_ni.xml with Nucleide and OpenMC and return difference stats."""
    nuc_out, om_out = run_chain_ni_vectors()

    max_dens = max(max(abs(v) for v in nuc_out.values()), max(abs(v) for v in om_out.values()))
    tiny = 1.0e-12 * max_dens
    diffs: list[float] = []
    for nuc in om_out:
        a = nuc_out.get(nuc, 0.0)
        b = om_out[nuc]
        scale = max(abs(a), abs(b))
        if scale < tiny:
            # Numerical noise around zero: scale by the largest density.
            diffs.append(abs(a - b) / max_dens)
        else:
            diffs.append(abs(a - b) / scale)
    return {
        "max_rel_diff": max(diffs),
        "mean_rel_diff": sum(diffs) / len(diffs),
        "max_density": max_dens,
    }


def bateman_3_nuclide() -> dict[str, float]:
    """Compare Nucleide/OpenMC for I135 -> Xe135 -> Cs135 against Bateman."""
    chain = nucleide.depletion.read_chain(str(CHAIN_SIMPLE))
    om_chain = openmc.deplete.Chain.from_xml(str(CHAIN_SIMPLE))

    t_half_i = 2.36520e4
    t_half_xe = 3.29040e4
    lam_i = math.log(2) / t_half_i
    lam_xe = math.log(2) / t_half_xe
    dt = 1.0e5
    n0_i = 1.0e16
    n0 = {"I135": n0_i, "Xe135": 0.0, "Cs135": 0.0}

    # Nucleide
    nuc_out = nucleide.depletion.deplete(chain, n0, dt, order=48)

    # OpenMC: build a zero-reaction-rates matrix so only decay remains.
    reactions = sorted({r.type for nuc in om_chain.nuclides for r in nuc.reactions})
    om_rates = openmc.deplete.ReactionRates(
        ["0"], [nuc.name for nuc in om_chain.nuclides], reactions
    )
    n0_vec = np.array([n0.get(nuc.name, 0.0) for nuc in om_chain.nuclides])
    A = om_chain.form_matrix(om_rates[0])
    om_out_vec = openmc.deplete.cram.CRAM48(A, n0_vec, dt)
    om_out = {nuc.name: float(om_out_vec[i]) for i, nuc in enumerate(om_chain.nuclides)}

    # Analytic (closed-form Bateman) for A -> B -> C (C stable).
    exp_i = math.exp(-lam_i * dt)
    exp_xe = math.exp(-lam_xe * dt)
    anal_i = n0_i * exp_i
    anal_xe = n0_i * lam_i / (lam_xe - lam_i) * (exp_i - exp_xe)
    anal_cs = n0_i - anal_i - anal_xe

    results: dict[str, float] = {}
    for label, nuc_val, om_val, anal_val in [
        ("I135", nuc_out["I135"], om_out["I135"], anal_i),
        ("Xe135", nuc_out["Xe135"], om_out["Xe135"], anal_xe),
        ("Cs135", nuc_out["Cs135"], om_out["Cs135"], anal_cs),
    ]:
        rel_nuc = abs(nuc_val - anal_val) / max(abs(anal_val), 1.0)
        rel_omc = abs(om_val - anal_val) / max(abs(anal_val), 1.0)
        results[f"nucleide_vs_analytic_{label}"] = rel_nuc
        results[f"openmc_vs_analytic_{label}"] = rel_omc
    return results


# --- 0.3.0 Predictor time-series check (synthetic A -> B -> C, no fixtures) ---
#: Decay constants for the synthetic convergence chain [1/s].
ABC_LAM_A = 1.0e-6
ABC_LAM_B = 1.0e-5
#: Initial atoms of A (B and C start at zero) and the five step lengths [s].
ABC_N0 = 1.0e15
ABC_DTS = [2.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4]


def _write_abc_chain() -> str:
    """Write the synthetic A -> B -> C decay chain to a temp file (no fixtures)."""
    text = (
        "<depletion_chain>\n"
        f'  <nuclide name="A" half_life="{math.log(2) / ABC_LAM_A!r}">\n'
        '    <decay type="beta" target="B" branching_ratio="1.0"/>\n'
        "  </nuclide>\n"
        f'  <nuclide name="B" half_life="{math.log(2) / ABC_LAM_B!r}">\n'
        '    <decay type="beta" target="C" branching_ratio="1.0"/>\n'
        "  </nuclide>\n"
        '  <nuclide name="C" reactions="0"/>\n'
        "</depletion_chain>\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".xml", delete=False) as f:
        f.write(text)
        return f.name


def _abc_bateman(t: float) -> dict[str, float]:
    """Analytic Bateman solution for A -> B -> C (C stable) from pure A."""
    exp_a = math.exp(-ABC_LAM_A * t)
    exp_b = math.exp(-ABC_LAM_B * t)
    a = ABC_N0 * exp_a
    b = ABC_N0 * ABC_LAM_A / (ABC_LAM_B - ABC_LAM_A) * (exp_a - exp_b)
    return {"A": a, "B": b, "C": ABC_N0 - a - b}


def _heat_total(entry) -> float:
    """Normalize one `decay_heat` series entry (scalar or per-nuclide dict)."""
    if isinstance(entry, dict):
        return sum(float(v) for v in entry.values())
    return float(entry)


def run_bateman_crosscheck() -> dict:
    """Cross-check Bateman/HP fast path vs CRAM-48 on A -> B -> C.

    In-memory only (synthetic chain from `_write_abc_chain`, no fixtures):
    single-step `deplete` with `method="bateman"` / `"bateman_hp"` at
    ``dt = 1e5`` s vs CRAM-48 and vs the analytic Bateman solution, plus a
    Predictor `deplete_series` band check. Returns the worst relative
    differences and the atom-conservation errors.
    """
    path = _write_abc_chain()
    try:
        chain = nucleide.depletion.read_chain(path)
        n0 = {"A": ABC_N0, "B": 0.0, "C": 0.0}
        dt = 1.0e5
        ref = nucleide.depletion.deplete(chain, n0, dt, order=48)
        want = _abc_bateman(dt)
        series_ref = nucleide.depletion.deplete_series(
            chain, n0, list(ABC_DTS), integrator="predictor", order=48
        )
        out: dict[str, float] = {}
        for method in ("bateman", "bateman_hp"):
            got = nucleide.depletion.deplete(chain, n0, dt, method=method)
            series = nucleide.depletion.deplete_series(
                chain, n0, list(ABC_DTS), integrator="predictor", method=method
            )
            vs_cram = max(abs(got[n] - ref[n]) / max(abs(ref[n]), 1.0) for n in ("A", "B", "C"))
            vs_analytic = max(
                abs(got[n] - want[n]) / max(abs(want[n]), 1.0) for n in ("A", "B", "C")
            )
            total = sum(got.values())
            cons = abs(total - ABC_N0) / ABC_N0
            t = 0.0
            series_band = 0.0
            for step, (row, ref_row) in enumerate(
                zip(series["atoms"], series_ref["atoms"], strict=True)
            ):
                t += ABC_DTS[step]
                bat = _abc_bateman(t)
                series_band = max(
                    series_band,
                    max(abs(row[n] - bat[n]) / max(abs(bat[n]), 1.0) for n in "ABC"),
                    max(
                        abs(row[n] - ref_row[n]) / max(abs(ref_row[n]), 1.0)
                        for n in ("A", "B", "C")
                    ),
                )
            out[f"{method}_vs_cram"] = vs_cram
            out[f"{method}_vs_analytic"] = vs_analytic
            out[f"{method}_cons"] = cons
            out[f"{method}_series_band"] = series_band
        return out
    finally:
        os.unlink(path)


def run_predictor_series() -> dict:
    """Run Predictor `deplete_series` on A -> B -> C; check vs Bateman.

    Returns node table rows plus summary scalars for the activity/heat
    channel spot checks (heat resolves to 0.0: synthetic nuclides carry no
    decay-energy data, per the documented zero-heat behavior).
    """
    path = _write_abc_chain()
    try:
        chain = nucleide.depletion.read_chain(path)
        out = nucleide.depletion.deplete_series(
            chain,
            {"A": ABC_N0, "B": 0.0, "C": 0.0},
            list(ABC_DTS),
            integrator="predictor",
            order=48,
        )
    finally:
        os.unlink(path)

    node_rows: list[list[str]] = []
    node_diffs: list[float] = []
    cons_errs: list[float] = []
    act_self: list[float] = []
    act_bateman: list[float] = []
    act_c: list[float] = []
    t = 0.0
    for step, (atoms, act) in enumerate(zip(out["atoms"], out["activity"], strict=True)):
        t += ABC_DTS[step]
        want = _abc_bateman(t)
        diffs = [abs(atoms[n] - want[n]) / max(abs(want[n]), 1.0) for n in ("A", "B", "C")]
        node_diffs.append(max(diffs))
        total = sum(atoms.values())
        cons_errs.append(abs(total - ABC_N0) / ABC_N0)
        node_rows.append([str(step + 1), fmt(t), fmt(node_diffs[-1]), fmt(cons_errs[-1])])
        for n, lam in (("A", ABC_LAM_A), ("B", ABC_LAM_B)):
            denom = max(abs(lam * atoms[n]), 1.0e-30)
            act_self.append(abs(act[n] - lam * atoms[n]) / denom)
        act_bateman.append(abs(act["B"] - ABC_LAM_B * want["B"]) / max(ABC_LAM_B * want["B"], 1.0))
        act_c.append(abs(act["C"]))

    heat = [_heat_total(entry) for entry in out["decay_heat"]]
    return {
        "node_rows": node_rows,
        "node_max": max(node_diffs),
        "cons_max": max(cons_errs),
        "act_self": max(act_self),
        "act_bateman": max(act_bateman),
        "act_c": max(act_c),
        "heat": max(abs(v) for v in heat),
    }


def emit_report(report: Report) -> bool:
    """Write the JSON report; skip loudly when env metadata is missing.

    `common.environment()` needs installed PyNE/OpenMC distributions, so a
    bare checkout outside the container cannot write `results/*.json`. That
    write is skipped (never hand-written) while every check above still runs.
    """
    try:
        report.emit()
        return True
    except importlib.metadata.PackageNotFoundError as exc:
        print(
            "SKIPPED report write: validation environment metadata unavailable"
            f" outside the container ({exc}); checks above still ran."
        )
        return False


def main() -> int:
    report = Report("depletion", "Depletion (`depletion_vs_openmc.py`)")

    ni_stats = run_chain_ni_comparison() if HAS_OPENMC else None
    bate = bateman_3_nuclide() if HAS_OPENMC else None
    if not HAS_OPENMC:
        report.heading("CRAM-48 oracles (OpenMC)")
        print(OPENMC_SKIP)
        report.prose(
            f"SKIPPED: {OPENMC_SKIP} The chain_ni CRAM-48 and 3-nuclide Bateman"
            " comparisons require OpenMC; they rerun inside the validation container."
        )
    else:
        assert ni_stats is not None and bate is not None
        report.heading("CRAM-48 on `chain_ni.xml`")
        report.table(
            ["Metric", "Value"],
            [
                ["Max relative difference", fmt(ni_stats["max_rel_diff"])],
                ["Mean relative difference", fmt(ni_stats["mean_rel_diff"])],
            ],
        )
        report.prose(
            "Differences are relative to OpenMC's own CRAM-48 solver on the same"
            " patched chain file."
        )

        report.heading("3-nuclide analytic Bateman check")
        report.prose("Relative difference vs. analytic solution:")
        report.table(
            ["Nuclide", "Nucleide vs analytic", "OpenMC vs analytic"],
            [
                [
                    "I-135",
                    fmt(bate["nucleide_vs_analytic_I135"]),
                    fmt(bate["openmc_vs_analytic_I135"]),
                ],
                [
                    "Xe-135",
                    fmt(bate["nucleide_vs_analytic_Xe135"]),
                    fmt(bate["openmc_vs_analytic_Xe135"]),
                ],
                [
                    "Cs-135",
                    fmt(bate["nucleide_vs_analytic_Cs135"]),
                    fmt(bate["openmc_vs_analytic_Cs135"]),
                ],
            ],
        )

    series = run_predictor_series()
    bateman = run_bateman_crosscheck()
    report.heading("Predictor time-series convergence (synthetic A->B->C)")
    report.prose(
        "Predictor `deplete_series` (sequential CRAM-48 steps) on a synthetic"
        " A->B->C pure-decay chain built inline in a temp file (no fixtures:"
        " A: 1e-6 /s, B: 1e-5 /s, C: stable; five 2e4 s steps from 1e15 atoms"
        " of A) vs the analytic Bateman solution at every time node."
    )
    report.table(
        ["Step", "Time (s)", "Max rel diff vs Bateman", "Atom-conservation rel err"],
        series["node_rows"],
    )

    report.heading("Activity and decay-heat channels")
    report.prose(
        "Activity is A = λ·N per nuclide (C is stable, so exactly 0.0). The"
        " synthetic A/B/C nuclides carry no decay-energy data, so the heat"
        " channel resolves to 0.0 — the documented zero-heat behavior (the"
        " check asserts that resolution, not nonzero values)."
    )
    report.table(
        ["Check", "Value", "Tolerance"],
        [
            ["A/B activity self-consistency (max rel)", fmt(series["act_self"]), "1.0e-09"],
            ["B activity vs Bateman (max rel)", fmt(series["act_bateman"]), "1.0e-06"],
            ["C activity (max abs, Bq)", fmt(series["act_c"]), "0.0"],
            ["Decay heat (max abs, W)", fmt(series["heat"]), "0.0"],
        ],
    )

    print(
        f"predictor series: node_max={fmt(series['node_max'])}"
        f" cons_max={fmt(series['cons_max'])} act_self={fmt(series['act_self'])}"
        f" act_bateman={fmt(series['act_bateman'])} heat={fmt(series['heat'])}"
    )
    print(
        "bateman cross-check:"
        f" std_vs_cram={fmt(bateman['bateman_vs_cram'])}"
        f" hp_vs_cram={fmt(bateman['bateman_hp_vs_cram'])}"
        f" std_series={fmt(bateman['bateman_series_band'])}"
        f" hp_series={fmt(bateman['bateman_hp_series_band'])}"
    )

    report.heading("Bateman-vs-CRAM cross-check (synthetic A->B->C)")
    report.prose(
        "Analytic Bateman fast path (`method='bateman'` / `'bateman_hp'`) on"
        " the same synthetic A->B->C pure-decay chain built inline in a temp"
        " file (no fixtures): single-step `deplete` at dt = 1e5 s vs CRAM-48"
        " and vs the analytic solution, plus a Predictor `deplete_series`"
        " band check (max vs Bateman analytic and vs CRAM-48 at every node)."
    )
    report.table(
        ["Method", "Vs CRAM-48", "Vs analytic", "Conservation", "Series band"],
        [
            [
                "bateman",
                fmt(bateman["bateman_vs_cram"]),
                fmt(bateman["bateman_vs_analytic"]),
                fmt(bateman["bateman_cons"]),
                fmt(bateman["bateman_series_band"]),
            ],
            [
                "bateman_hp",
                fmt(bateman["bateman_hp_vs_cram"]),
                fmt(bateman["bateman_hp_vs_analytic"]),
                fmt(bateman["bateman_hp_cons"]),
                fmt(bateman["bateman_hp_series_band"]),
            ],
        ],
    )

    emit_report(report)

    # Exit nonzero if results are unexpectedly loose.
    if ni_stats is not None and ni_stats["max_rel_diff"] > 1.0e-8:
        print("FAIL: chain_ni CRAM diff larger than expected 1e-8", file=sys.stderr)
        return 1
    if bate is not None and any(v > 1.0e-8 for v in bate.values()):
        print("FAIL: analytic Bateman diff larger than expected 1e-8", file=sys.stderr)
        return 1
    if series["node_max"] > 1.0e-6:
        print("FAIL: Predictor series drift from Bateman beyond 1e-6", file=sys.stderr)
        return 1
    if series["cons_max"] > 1.0e-8:
        print("FAIL: Predictor series atom conservation beyond 1e-8", file=sys.stderr)
        return 1
    if series["act_self"] > 1.0e-9 or series["act_bateman"] > 1.0e-6:
        print("FAIL: Predictor series activity channel mismatch", file=sys.stderr)
        return 1
    if series["act_c"] != 0.0 or series["heat"] != 0.0:
        print("FAIL: stable-nuclide activity or zero-heat resolution broken", file=sys.stderr)
        return 1
    for key, value in bateman.items():
        if (key.endswith("_vs_cram") or key.endswith("_series_band")) and value > 1.0e-6:
            print(f"FAIL: Bateman cross-check {key} beyond 1e-6", file=sys.stderr)
            return 1
        if (key.endswith("_vs_analytic") or key.endswith("_cons")) and value > 1.0e-8:
            print(f"FAIL: Bateman cross-check {key} beyond 1e-8", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
