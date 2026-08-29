"""Compare Nucleide's CRAM depletion against OpenMC CRAM-48 and analytic Bateman."""

from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import openmc.deplete
from common import Report, fmt

import nucleide

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


def run_chain_ni_comparison() -> dict[str, float]:
    """Deplete chain_ni.xml with Nucleide and OpenMC and return difference stats."""
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


def main() -> int:
    report = Report("depletion", "Depletion (`depletion_vs_openmc.py`)")

    ni_stats = run_chain_ni_comparison()
    report.heading("CRAM-48 on `chain_ni.xml`")
    report.table(
        ["Metric", "Value"],
        [
            ["Max relative difference", fmt(ni_stats["max_rel_diff"])],
            ["Mean relative difference", fmt(ni_stats["mean_rel_diff"])],
        ],
    )
    report.prose(
        "Differences are relative to OpenMC's own CRAM-48 solver on the same patched\nchain file."
    )

    bate = bateman_3_nuclide()
    report.heading("3-nuclide analytic Bateman check")
    report.prose("Relative difference vs. analytic solution:")
    report.table(
        ["Nuclide", "Nucleide vs analytic", "OpenMC vs analytic"],
        [
            ["I-135", fmt(bate["nucleide_vs_analytic_I135"]), fmt(bate["openmc_vs_analytic_I135"])],
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

    report.emit()

    # Exit nonzero if results are unexpectedly loose.
    if ni_stats["max_rel_diff"] > 1.0e-8:
        print("FAIL: chain_ni CRAM diff larger than expected 1e-8", file=sys.stderr)
        return 1
    if any(v > 1.0e-8 for v in bate.values()):
        print("FAIL: analytic Bateman diff larger than expected 1e-8", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
