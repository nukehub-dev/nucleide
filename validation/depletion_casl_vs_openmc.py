"""Compare Nucleide's CRAM-48 depletion against OpenMC on the CASL/VERA chain.

The CASL simplified depletion chain (VERA Depletion Benchmark,
CASL-U-2015-1014-000; thermal-spectrum branching, 228 nuclides, fission
product yields) is downloaded once from OpenMC's pregenerated-chain archive
and cached under ``validation/.cache/`` (git-ignored).
"""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import openmc.deplete
from common import Report, fmt

import nucleide

SCRIPT_DIR = Path(__file__).resolve().parent
CACHE_DIR = SCRIPT_DIR / ".cache"
CHAIN_CASL = CACHE_DIR / "chain_casl_thermal.xml"
CHAIN_URL = "https://anl.box.com/shared/static/3nvnasacm2b56716oh5hyndxdyauh5gs.xml"

DT = 2.592e6  # 30 days in seconds

# Uniform one-group reaction rates [1/s] applied to both codes. The exercise
# compares solvers on identical matrices, not physics.
RATE_TABLE = {
    "(n,gamma)": 1.0e-7,
    "fission": 1.0e-8,
    "(n,2n)": 1.0e-8,
    "(n,3n)": 1.0e-9,
    "(n,4n)": 1.0e-9,
    "(n,p)": 1.0e-8,
    "(n,a)": 1.0e-8,
}

# Fresh-UO2-style initial inventory (atom counts; ~5 wt% U-235 enrichment).
INITIAL = {
    "U235": 5.0e22,
    "U238": 9.5e23,
    "O16": 2.0e24,
}

NOTABLE = ["U235", "U238", "Pu239", "Xe135", "I135", "Cs137", "Sm149"]


def ensure_chain() -> Path:
    """Return the cached CASL chain path, downloading it once if absent."""
    if CHAIN_CASL.exists():
        return CHAIN_CASL
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Downloading CASL chain from {CHAIN_URL} -> {CHAIN_CASL}")
    req = urllib.request.Request(CHAIN_URL, headers={"User-Agent": "nucleide-validation"})
    with urllib.request.urlopen(req) as resp, CHAIN_CASL.open("wb") as out:
        out.write(resp.read())
    return CHAIN_CASL


def nucleide_casl_rates(chain: nucleide.depletion.Chain) -> dict[str, float]:
    """Build one-group rates for the CASL chain in Nucleide's dict format."""
    rates: dict[str, float] = {}
    for nuc in chain.nuclides:
        for kind, value in RATE_TABLE.items():
            rates[f"{nuc}:{kind}"] = value
    return rates


def openmc_casl_rates(chain: openmc.deplete.Chain) -> openmc.deplete.ReactionRates:
    """Build matching one-group rates for OpenMC's form_matrix."""
    reactions = sorted({r.type for nuc in chain.nuclides for r in nuc.reactions})
    nuc_names = [nuc.name for nuc in chain.nuclides]
    rates = openmc.deplete.ReactionRates(["0"], nuc_names, reactions)
    for nuc in chain.nuclides:
        i = nuc_names.index(nuc.name)
        for r in nuc.reactions:
            j = reactions.index(r.type)
            rates[0, i, j] = RATE_TABLE.get(r.type, 0.0)
    return rates


def run_casl_vectors() -> tuple[dict[str, float], dict[str, float]]:
    """Deplete the CASL chain with both codes; return final density dicts."""
    chain_path = str(ensure_chain())
    nuc_chain = nucleide.depletion.read_chain(chain_path)
    om_chain = openmc.deplete.Chain.from_xml(chain_path)

    n0 = {nuc: INITIAL.get(nuc, 0.0) for nuc in nuc_chain.nuclides}

    nuc_out = nucleide.depletion.deplete(
        nuc_chain, n0, DT, rates=nucleide_casl_rates(nuc_chain), order=48
    )

    om_rates = openmc_casl_rates(om_chain)
    A = om_chain.form_matrix(om_rates[0])
    n0_vec = np.array([n0.get(nuc.name, 0.0) for nuc in om_chain.nuclides])
    om_out_vec = openmc.deplete.cram.CRAM48(A, n0_vec, DT)
    om_out = {nuc.name: float(om_out_vec[i]) for i, nuc in enumerate(om_chain.nuclides)}
    return dict(nuc_out), om_out


def vector_diff_stats(nuc_out: dict[str, float], om_out: dict[str, float]) -> dict[str, float]:
    """Relative-difference statistics between two density vectors."""
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


def main() -> int:
    report = Report("depletion_casl", "Depletion on the CASL chain (`depletion_casl_vs_openmc.py`)")

    nuc_out, om_out = run_casl_vectors()
    stats = vector_diff_stats(nuc_out, om_out)

    report.prose(
        "CASL simplified depletion chain (VERA Depletion Benchmark, CASL-U-2015-1014-000;"
        f"\n{len(om_out)} nuclides, thermal-spectrum branching, fission product yields),"
        " downloaded\n"
        f"from the OpenMC pregenerated-chain archive and cached at `.cache/{CHAIN_CASL.name}`.\n"
        "Fresh-UO2-style initial inventory (U-235/U-238/O-16); identical uniform one-group\n"
        "fission and capture rates on both sides; single 30-day CRAM-48 step. This compares\n"
        "the two solvers on identical matrices, not physics."
    )

    report.heading("Full density vector after one 30-day CRAM-48 step")
    report.table(
        ["Metric", "Value"],
        [
            ["Max relative difference", fmt(stats["max_rel_diff"])],
            ["Mean relative difference", fmt(stats["mean_rel_diff"])],
        ],
    )
    report.prose(
        "Relative differences are computed per nuclide where the larger of the two"
        "\ndensities exceeds 1e-12 x the maximum density; smaller densities are scaled by"
        "\nthe maximum density instead."
    )

    report.heading("Notable nuclides")
    max_dens = stats["max_density"]
    rows: list[list[str]] = []
    for nuc in NOTABLE:
        if nuc not in om_out:
            continue
        a = nuc_out.get(nuc, 0.0)
        b = om_out[nuc]
        scale = max(abs(a), abs(b))
        diff = abs(a - b) / (scale if scale >= 1.0e-12 * max_dens else max_dens)
        rows.append([nuc, fmt(a), fmt(b), fmt(diff)])
    report.table(["Nuclide", "Nucleide", "OpenMC", "Rel diff"], rows)

    report.emit()

    if stats["max_rel_diff"] > 1.0e-8:
        print("FAIL: CASL chain CRAM diff larger than expected 1e-8", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
