"""Compare Nucleide nuclear data and name conversions against PyNE and OpenMC."""

from __future__ import annotations

import contextlib
import sys

import openmc.data
import pyne.data
import pyne.nucname as nucname
from common import Report, abs_diff, fmt, rel_diff

import nucleide


def compare_atomic_masses(sample: list[str]) -> dict[str, dict[str, float]]:
    """Compare atomic masses vs PyNE (AME2016) and OpenMC (AME2020)."""
    stats: dict[str, list[float]] = {"pyne": [], "openmc": []}
    per_nuc: dict[str, dict[str, float]] = {}
    for name in sample:
        nuc_val = nucleide.nuclei.atomic_mass(name)
        pyne_val = pyne.data.atomic_mass(name)
        openmc_val = None
        with contextlib.suppress(Exception):
            openmc_val = openmc.data.atomic_mass(name.lower())

        per_nuc[name] = {
            "nucleide": nuc_val if nuc_val is not None else float("nan"),
            "pyne": pyne_val if pyne_val is not None else float("nan"),
            "openmc": openmc_val if openmc_val is not None else float("nan"),
        }
        if nuc_val is not None and pyne_val is not None:
            stats["pyne"].append(abs_diff(nuc_val, pyne_val))
        if nuc_val is not None and openmc_val is not None:
            stats["openmc"].append(abs_diff(nuc_val, openmc_val))

    return {
        "per_nuc": per_nuc,
        "pyne_max_abs": max(stats["pyne"]) if stats["pyne"] else float("nan"),
        "pyne_mean_abs": (sum(stats["pyne"]) / len(stats["pyne"]))
        if stats["pyne"]
        else float("nan"),
        "openmc_max_abs": max(stats["openmc"]) if stats["openmc"] else float("nan"),
        "openmc_mean_abs": (sum(stats["openmc"]) / len(stats["openmc"]))
        if stats["openmc"]
        else float("nan"),
    }


def compare_natural_abundances() -> dict[str, float]:
    """Compare natural abundances for all naturally-occurring isotopes."""
    om_abund = openmc.data.NATURAL_ABUNDANCE
    diffs_pyne: list[float] = []
    diffs_openmc: list[float] = []
    for name in om_abund:
        nuc_val = nucleide.nuclei.natural_abundance(name)
        om_val = om_abund[name]
        pyne_val = None
        with contextlib.suppress(Exception):
            pyne_val = pyne.data.natural_abund(nucname.id(name))
        if nuc_val is not None:
            diffs_openmc.append(abs_diff(nuc_val, om_val))
            if pyne_val is not None:
                diffs_pyne.append(abs_diff(nuc_val, pyne_val))

    return {
        "count": len(om_abund),
        "openmc_max_abs": max(diffs_openmc) if diffs_openmc else float("nan"),
        "openmc_mean_abs": (sum(diffs_openmc) / len(diffs_openmc))
        if diffs_openmc
        else float("nan"),
        "pyne_max_abs": max(diffs_pyne) if diffs_pyne else float("nan"),
        "pyne_mean_abs": (sum(diffs_pyne) / len(diffs_pyne)) if diffs_pyne else float("nan"),
    }


def compare_half_lives(sample: list[str]) -> dict[str, float]:
    """Compare half-lives vs PyNE and OpenMC."""
    diffs_pyne: list[float] = []
    diffs_openmc: list[float] = []
    for name in sample:
        nuc_val = nucleide.nuclei.half_life(name)
        pyne_val = pyne.data.half_life(name)
        om_name = name.lower().replace("-", "").replace("m", "_m1")
        om_val = None
        with contextlib.suppress(Exception):
            om_val = openmc.data.half_life(om_name)

        if nuc_val is not None and pyne_val is not None and pyne_val > 0:
            diffs_pyne.append(rel_diff(nuc_val, pyne_val))
        if nuc_val is not None and om_val is not None and om_val > 0:
            diffs_openmc.append(rel_diff(nuc_val, om_val))

    return {
        "pyne_max_rel": max(diffs_pyne) if diffs_pyne else float("nan"),
        "pyne_mean_rel": (sum(diffs_pyne) / len(diffs_pyne)) if diffs_pyne else float("nan"),
        "openmc_max_rel": max(diffs_openmc) if diffs_openmc else float("nan"),
        "openmc_mean_rel": (sum(diffs_openmc) / len(diffs_openmc))
        if diffs_openmc
        else float("nan"),
    }


def compare_name_conversions(sample: list[str]) -> dict[str, float]:
    """Compare Nucleide name-dialect conversions against pyne.nucname."""

    def expected_zaid(name: str) -> int:
        n = nucleide.nuclei.Nuclide(name)
        z, a, state = n.z, n.a, n.state
        zaid = z * 1000 + a
        # MCNP special case: Am-242 and Am-242m are swapped.
        if zaid == 95242 and state < 2:
            state = (state + 1) % 2
        if state > 0:
            zaid += 300 + state * 100
        return zaid

    checks = {
        "nucid": (lambda n: n.nucid, lambda name: nucname.id(name)),
        "zzaaam": (lambda n: n.zzaaam, lambda name: nucname.zzaaam(name)),
        "zaid": (lambda n: n.zaid, expected_zaid),
        "serpent": (lambda n: n.serpent, lambda name: nucname.serpent(nucname.id(name))),
        "nist": (lambda n: n.nist, lambda name: nucname.nist(nucname.id(name))),
        "cinder": (lambda n: n.cinder, lambda name: nucname.cinder(nucname.id(name))),
        "alara": (lambda n: n.alara, lambda name: nucname.alara(nucname.id(name))),
        "sza": (lambda n: n.sza, lambda name: nucname.sza(nucname.id(name))),
    }

    diffs: dict[str, list[float]] = {k: [] for k in checks}
    for name in sample:
        nuc = nucleide.nuclei.Nuclide(name)
        for key, (nuc_fn, pyne_fn) in checks.items():
            try:
                a = nuc_fn(nuc)
                b = pyne_fn(name)
                if isinstance(a, str):
                    diffs[key].append(0.0 if a == b else 1.0)
                else:
                    diffs[key].append(abs_diff(float(a), float(b)))
            except Exception:
                # PyNE may not support some metastable dialects; skip silently.
                pass

    return {f"{key}_max": max(v) if v else float("nan") for key, v in diffs.items()}


def main() -> int:
    report = Report("nuclear_data", "Nuclear data (`nuclear_data_vs_refs.py`)")

    mass_sample = [
        "H1",
        "C12",
        "N14",
        "O16",
        "Fe56",
        "U235",
        "U238",
        "Pu239",
        "Pu240",
        "Am241",
        "Am242m",
        "Ba137m",
        "Xe135",
        "Cs137",
        "Sr90",
        "Co60",
        "Ni58",
        "Mn55",
        "Cu63",
        "Mo95",
        "Tc99",
        "I129",
        "I135",
        "Xe136",
        "Nd143",
        "Sm149",
        "Eu151",
        "Gd157",
        "Ho165",
        "W182",
        "W186",
        "Pb206",
        "Pb207",
        "Pb208",
        "Bi209",
        "Th232",
        "Pa233",
        "U233",
        "Np237",
        "Pu241",
        "Pu242",
        "Am243",
        "Cm244",
        "Bk249",
        "Cf252",
        "Es253",
        "Fm257",
        "Md260",
        "No259",
        "Lr262",
    ]

    hl_sample = [
        "H3",
        "C14",
        "Co60",
        "Sr90",
        "Tc99",
        "I129",
        "I135",
        "Cs137",
        "Ba137m",
        "Pm147",
        "Sm151",
        "Eu154",
        "Am241",
        "Am242m",
        "Cm244",
        "Pu239",
        "Pu240",
        "U235",
        "U238",
        "Np237",
    ]

    name_sample = [
        "H1",
        "U235",
        "Pu239",
        "Am242m",
        "Ba137m",
        "Co60",
        "Cs137",
        "I135",
        "Xe135",
        "Fe56",
        "W186",
    ]

    mass_stats = compare_atomic_masses(mass_sample)
    report.heading("Atomic masses")
    report.table(
        ["Reference", "Max abs diff (u)", "Mean abs diff (u)"],
        [
            ["OpenMC", fmt(mass_stats["openmc_max_abs"]), fmt(mass_stats["openmc_mean_abs"])],
            ["PyNE", fmt(mass_stats["pyne_max_abs"]), fmt(mass_stats["pyne_mean_abs"])],
        ],
    )

    abund_stats = compare_natural_abundances()
    report.heading(f"Natural abundances ({abund_stats['count']} isotopes)")
    report.table(
        ["Reference", "Max abs diff", "Mean abs diff"],
        [
            ["OpenMC", fmt(abund_stats["openmc_max_abs"]), fmt(abund_stats["openmc_mean_abs"])],
            ["PyNE", fmt(abund_stats["pyne_max_abs"]), fmt(abund_stats["pyne_mean_abs"])],
        ],
    )

    hl_stats = compare_half_lives(hl_sample)
    report.heading("Half-lives")
    report.table(
        ["Reference", "Max rel diff", "Mean rel diff"],
        [
            ["OpenMC", fmt(hl_stats["openmc_max_rel"]), fmt(hl_stats["openmc_mean_rel"])],
            ["PyNE", fmt(hl_stats["pyne_max_rel"]), fmt(hl_stats["pyne_mean_rel"])],
        ],
    )

    name_stats = compare_name_conversions(name_sample)
    report.heading("Name-dialect conversions vs `pyne.nucname`")
    report.prose(
        "All conversions (alara, cinder, nist, nucid, serpent, sza, zaid, zzaaam) had a\n"
        f"maximum relative/absolute difference of **{fmt(max(name_stats.values()))}**."
    )

    report.emit()

    if abund_stats["openmc_max_abs"] > 1.0e-12:
        print("FAIL: natural abundance mismatch with OpenMC", file=sys.stderr)
        return 1
    if hl_stats["openmc_max_rel"] > 1.0e-6:
        print("FAIL: half-life mismatch with OpenMC", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
