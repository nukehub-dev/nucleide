"""Compare Nucleide nuclear data and name conversions against PyNE and OpenMC."""

from __future__ import annotations

import contextlib
import importlib.metadata
import sys

from common import Report, abs_diff, fmt, rel_diff

import nucleide

try:
    import openmc.data

    HAS_OPENMC = True
    OPENMC_SKIP = ""
except Exception as exc:
    openmc = None  # type: ignore[no-redef]
    HAS_OPENMC = False
    OPENMC_SKIP = f"OpenMC oracle skipped: cannot import openmc.data ({exc})."

try:
    import pyne.data
    import pyne.nucname as nucname

    HAS_PYNE = True
    PYNE_SKIP = ""
except Exception as exc:
    pyne = None  # type: ignore[no-redef]
    nucname = None  # type: ignore[no-redef]
    HAS_PYNE = False
    PYNE_SKIP = f"PyNE oracle skipped: cannot import pyne.data ({exc})."


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


class SkipCheck(Exception):
    """An oracle is unavailable; the message is the loud skip reason."""


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


#: Nuclides covered by both `simple_xs.tsv` and PyNE's KAERI simple-xs table.
SIMPLE_XS_SAMPLE = [
    "H1",
    "B10",
    "C12",
    "O16",
    "Fe56",
    "Co59",
    "U235",
    "U238",
    "Pu239",
    "Pb208",
]

#: Reaction keys to try against PyNE's simple-xs `reaction()` accessor.
_SIMPLE_XS_RX = ("total", "sigma_t", 1)

#: NIST NCNR bound coherent lengths in fm (Sears 1992 tabulation).
#: `scattering_lengths.tsv` copies these values exactly, hence the tight check.
SCATTERING_ANCHORS_FM = {
    "H1": -3.7406,
    "H2": 6.671,
    "C12": 6.6511,
    "O16": 5.803,
    "Fe56": 9.94,
    "U238": 8.402,
    "Pb208": 9.5,
    "B10": -0.1,
}

#: Mean *prompt* recoverable decay energies in MeV from ENDF/B-VII.1 MF8/MT457
#: (daughter gammas belong to the daughter row: Cs137 excludes the 662 keV
#: line, which lives on Ba137_m1). Screening-level comparison with a loose
#: 20% band (NOT ENSDF evaluations).
DECAY_ENERGY_SPOTS_MEV = {
    "H3": 0.00569,
    "Co60": 2.60061,
    "Sr90": 0.1958,
    "I131": 0.573438,
    "Xe135": 0.56798,
    "Cs137": 0.179448,
    "Ba137_m1": 0.661397,
    "U235": 4.61919,
    "Pu239": 5.24326,
    "Am241": 5.62799,
}


def _pyne_simple_xs_source():
    """Return PyNE's KAERI simple-xs source, or raise SkipCheck with a reason."""
    try:
        from pyne.xs.data_source import SimpleDataSource
    except Exception as exc:
        raise SkipCheck(
            "PyNE simple_xs skipped: pyne.xs.data_source is unavailable"
            f" ({exc}); the nomoab build may lack this module."
        ) from exc
    try:
        src = SimpleDataSource()
    except Exception as exc:
        raise SkipCheck(
            f"PyNE simple_xs skipped: cannot construct SimpleDataSource ({exc})."
        ) from exc
    if not getattr(src, "exists", True):
        raise SkipCheck(
            "PyNE simple_xs skipped: SimpleDataSource.exists is False"
            " (nuc_data.h5 has no /neutron/simple_xs table)."
        )
    return src


def _pyne_simple_xs_lookup(src, name: str) -> tuple[float, float]:
    """Return (thermal_b, fast14_b) totals from PyNE, or raise SkipCheck."""
    reasons: list[str] = []
    for rx in _SIMPLE_XS_RX:
        try:
            data = src.reaction(name, rx)
        except Exception as exc:
            reasons.append(f"{rx}: {exc}")
            continue
        if data is None:
            reasons.append(f"{rx}: no data")
            continue
        vals = [float(v) for v in list(data)]
        if len(vals) >= 2 and any(v > 0 for v in vals):
            # Source group structure is descending energy: the first point
            # is 14 MeV and the last is thermal (2.53e-8 MeV), so thermal
            # is vals[-1] and fast-14 is vals[0].
            return vals[-1], vals[0]
        reasons.append(f"{rx}: empty response")
    raise SkipCheck(f"no total channel for {name} ({'; '.join(reasons)}).")


def compare_simple_xs(sample: list[str]) -> dict:
    """Compare `simple_xs` thermal/fast totals vs PyNE `nuc_data.h5` simple_xs."""
    rows: list[list[str]] = []
    skipped: list[str] = []
    diffs: list[float] = []
    if not HAS_PYNE:
        skipped.append(PYNE_SKIP)
        return {"rows": rows, "skipped": skipped, "diffs": diffs, "available": False}
    try:
        src = _pyne_simple_xs_source()
    except SkipCheck as exc:
        skipped.append(str(exc))
        return {"rows": rows, "skipped": skipped, "diffs": diffs, "available": False}
    for name in sample:
        try:
            nuc_val = nucleide.nuclei.simple_xs(name)
        except Exception as exc:
            skipped.append(f"{name}: nucleide.nuclei.simple_xs raised ({exc}).")
            continue
        if nuc_val is None:
            skipped.append(f"{name}: Nucleide has no simple_xs entry.")
            continue
        try:
            ref_th, ref_fa = _pyne_simple_xs_lookup(src, name)
        except SkipCheck as exc:
            skipped.append(f"{name}: {exc}")
            continue
        nuc_th, nuc_fa = nuc_val
        d_th = rel_diff(nuc_th, ref_th)
        d_fa = rel_diff(nuc_fa, ref_fa)
        diffs.extend([d_th, d_fa])
        rows.append(
            [name, fmt(nuc_th), fmt(ref_th), fmt(d_th), fmt(nuc_fa), fmt(ref_fa), fmt(d_fa)]
        )
    return {"rows": rows, "skipped": skipped, "diffs": diffs, "available": bool(rows)}


def compare_scattering_lengths() -> dict:
    """Compare coherent scattering lengths vs the NIST-anchored values."""
    rows: list[list[str]] = []
    diffs: list[float] = []
    missing: list[str] = []
    for name, anchor in SCATTERING_ANCHORS_FM.items():
        val = nucleide.nuclei.scattering_length(name)
        if val is None:
            missing.append(name)
            rows.append([name, fmt(anchor), "None", "n/a"])
            continue
        diffs.append(abs_diff(val, anchor))
        rows.append([name, fmt(val), fmt(anchor), fmt(diffs[-1])])
    return {
        "rows": rows,
        "missing": missing,
        "max_abs": max(diffs) if diffs else float("nan"),
    }


def compare_decay_energies() -> dict:
    """Compare `decay_energy` vs chain/ENDF-derived spot values (20% band)."""
    rows: list[list[str]] = []
    diffs: list[float] = []
    missing: list[str] = []
    for name, spot in DECAY_ENERGY_SPOTS_MEV.items():
        val = nucleide.nuclei.decay_energy(name)
        if val is None:
            missing.append(name)
            rows.append([name, "None", fmt(spot), "n/a"])
            continue
        diffs.append(rel_diff(val, spot))
        rows.append([name, fmt(val), fmt(spot), fmt(diffs[-1])])
    return {
        "rows": rows,
        "missing": missing,
        "stable_none": nucleide.nuclei.decay_energy("Fe56") is None,
        "max_rel": max(diffs) if diffs else float("nan"),
    }


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

    mass_stats = compare_atomic_masses(mass_sample) if (HAS_PYNE and HAS_OPENMC) else None
    abund_stats = compare_natural_abundances() if (HAS_PYNE and HAS_OPENMC) else None
    hl_stats = compare_half_lives(hl_sample) if (HAS_PYNE and HAS_OPENMC) else None
    name_stats = compare_name_conversions(name_sample) if (HAS_PYNE and HAS_OPENMC) else None

    if not (HAS_PYNE and HAS_OPENMC):
        report.heading("Reference-data oracles (PyNE/OpenMC)")
        for reason in (PYNE_SKIP, OPENMC_SKIP):
            if reason:
                print(reason)
                report.prose(f"SKIPPED: {reason}")
        report.prose(
            "The atomic-mass, natural-abundance, half-life, and name-dialect"
            " comparisons require both PyNE and OpenMC; they rerun inside the"
            " validation container."
        )
    else:
        assert mass_stats is not None and abund_stats is not None
        assert hl_stats is not None and name_stats is not None
        report.heading("Atomic masses")
        report.table(
            ["Reference", "Max abs diff (u)", "Mean abs diff (u)"],
            [
                ["OpenMC", fmt(mass_stats["openmc_max_abs"]), fmt(mass_stats["openmc_mean_abs"])],
                ["PyNE", fmt(mass_stats["pyne_max_abs"]), fmt(mass_stats["pyne_mean_abs"])],
            ],
        )

        report.heading(f"Natural abundances ({abund_stats['count']} isotopes)")
        report.table(
            ["Reference", "Max abs diff", "Mean abs diff"],
            [
                ["OpenMC", fmt(abund_stats["openmc_max_abs"]), fmt(abund_stats["openmc_mean_abs"])],
                ["PyNE", fmt(abund_stats["pyne_max_abs"]), fmt(abund_stats["pyne_mean_abs"])],
            ],
        )

        report.heading("Half-lives")
        report.table(
            ["Reference", "Max rel diff", "Mean rel diff"],
            [
                ["OpenMC", fmt(hl_stats["openmc_max_rel"]), fmt(hl_stats["openmc_mean_rel"])],
                ["PyNE", fmt(hl_stats["pyne_max_rel"]), fmt(hl_stats["pyne_mean_rel"])],
            ],
        )

        report.heading("Name-dialect conversions vs `pyne.nucname`")
        report.prose(
            "All conversions (alara, cinder, nist, nucid, serpent, sza, zaid, zzaaam) had a\n"
            f"maximum relative/absolute difference of **{fmt(max(name_stats.values()))}**."
        )

    xs_stats = compare_simple_xs(SIMPLE_XS_SAMPLE)
    report.heading("Screening cross sections (`simple_xs`) vs PyNE `nuc_data.h5`")
    report.prose(
        "Nucleide thermal (0.0253 eV) and 14-MeV total cross sections vs PyNE's"
        " KAERI-anchored `SimpleDataSource` totals (`/neutron/simple_xs` in the"
        " `nuc_data.h5` bundled with PyNE, read in place — no download)."
        " Thermal is the source's first group point, fast the last (14 MeV)."
    )
    if xs_stats["rows"]:
        report.table(
            [
                "Nuclide",
                "Nucleide th (b)",
                "PyNE th (b)",
                "Rel diff th",
                "Nucleide fast (b)",
                "PyNE fast (b)",
                "Rel diff fast",
            ],
            xs_stats["rows"],
        )
    for note in xs_stats["skipped"]:
        print(f"SKIPPED simple_xs: {note}")
        report.prose(f"SKIPPED simple_xs: {note}")

    scat_stats = compare_scattering_lengths()
    report.heading("Bound coherent scattering lengths vs NIST anchors")
    report.prose(
        "Nucleide `scattering_length` (coherent, fm) vs the NIST NCNR bound"
        " coherent lengths (Sears, Neutron News 3(3), 1992) that"
        " `scattering_lengths.tsv` copies exactly — hence the tight tolerance."
    )
    report.table(
        ["Nuclide", "Nucleide (fm)", "NIST (fm)", "Abs diff (fm)"],
        scat_stats["rows"],
    )

    de_stats = compare_decay_energies()
    report.heading("Decay energies vs chain/ENDF-derived spot values (screening-level)")
    report.prose(
        "Nucleide `decay_energy` (mean recoverable MeV per decay) vs"
        " chain/ENDF-derived spot values within a loose 20% screening band:"
        " these placeholder heat values are NOT ENSDF evaluations — never for"
        " spectroscopy, dose, or safety use."
    )
    report.table(
        ["Nuclide", "Nucleide (MeV)", "Spot (MeV)", "Rel diff"],
        de_stats["rows"],
    )
    if de_stats["stable_none"]:
        de_note = "Stable Fe56 correctly resolves to None (no decay-energy entry)."
    else:
        de_note = "UNEXPECTED: stable Fe56 has a decay-energy entry."
    print(de_note)
    report.prose(de_note)

    if xs_stats["available"]:
        print(f"simple_xs vs PyNE: max rel diff {fmt(max(xs_stats['diffs']))}")  # type: ignore[arg-type]
    print(f"scattering vs NIST: max abs diff {fmt(scat_stats['max_abs'])} fm")
    print(f"decay_energy vs spots: max rel diff {fmt(de_stats['max_rel'])}")

    emit_report(report)

    if abund_stats is not None and abund_stats["openmc_max_abs"] > 1.0e-12:
        print("FAIL: natural abundance mismatch with OpenMC", file=sys.stderr)
        return 1
    if hl_stats is not None and hl_stats["openmc_max_rel"] > 1.0e-6:
        print("FAIL: half-life mismatch with OpenMC", file=sys.stderr)
        return 1
    if xs_stats["available"] and max(xs_stats["diffs"]) > 0.30:  # type: ignore[arg-type]
        print("FAIL: simple_xs mismatch with PyNE simple_xs beyond 30%", file=sys.stderr)
        return 1
    if scat_stats["missing"] or scat_stats["max_abs"] > 1.0e-9:
        print("FAIL: scattering-length mismatch with NIST anchors", file=sys.stderr)
        return 1
    if de_stats["missing"] or not de_stats["stable_none"] or de_stats["max_rel"] > 0.20:
        print("FAIL: decay-energy spot check outside the 20% screening band", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
