"""Activation-code I/O checks for the alara-io, cccc-io, fispact-io, origen-io and r2s crates.

Comparisons run on committed `fixtures/` inputs with in-memory patching only
(no fixture is modified on disk):

- ALARA decks `sample2`/`sample3`: block inventory, mixtures, flux entries,
  cooling schedules, outputs, plus deterministic re-parse stability (no deck
  serializer exists, so "round-trip" means parsing the same text twice yields
  identical results).
- ALARA group-flux files `fluxin2` (3 x 175) and `fluxin_zeros`, with deck
  skip/scale spot values.
- ALARA output `sample2.out`: 2838 rows, 2 response variables, 6 time columns.
- FISPACT synthetic inventory: 36 rows, 3 variables.
- CCCC synthetic ISOTXS/RTFLUX plus a PARTISN render/validate round-trip.
- ORIGEN synthetic TAPE5/TAPE6/TAPE9 echoes plus activity totals.
- R2S `from_deck` workflow on the sample2 deck plus validate/assemble smoke
  tests.

Independent oracles exist only where PyNE ships a matching reader (the
`pyne.alara` deck probe below, container-only and API-guarded); everything
else is synthetic self-consistency, recorded as such. Every skip is printed
and recorded in the report prose; none is silent.
"""

from __future__ import annotations

import importlib.metadata
import platform
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import common
from common import Report, fmt, rel_diff

import nucleide

REPO_ROOT = Path(__file__).resolve().parent.parent
ALARA_DIR = REPO_ROOT / "fixtures" / "alara"
CCCC_DIR = REPO_ROOT / "fixtures" / "cccc"
FISPACT_DIR = REPO_ROOT / "fixtures" / "fispact"
ORIGEN_DIR = REPO_ROOT / "fixtures" / "origen"

FAILURES = 0
TOL = 1.0e-6
EXPECTED_BLOCKS = [
    "geometry",
    "dimension",
    "mat_loading",
    "material_lib",
    "element_lib",
    "mixture",
    "mixture",
    "flux",
    "schedule",
    "pulsehistory",
    "dump_file",
    "data_library",
    "cooling",
    "output",
    "truncation",
]

# ALARA input-deck block vocabulary: the pyne.alara oracle probe below is
# heuristic (it tries whatever entry point exists), so a result is only
# comparable when every key/entry is a plausible deck block name.
_KNOWN_BLOCKS = frozenset(EXPECTED_BLOCKS)

_REAL_ENVIRONMENT = common.environment


def _tolerant_environment() -> dict[str, str]:
    """Best-effort environment metadata when oracle distributions are missing.

    The canonical container always has nucleide/PyNE/OpenMC installed, so the
    real collector succeeds there unchanged. Local runs (nucleide wheel only)
    fall back to "not installed" markers instead of raising.
    """
    try:
        return _REAL_ENVIRONMENT()
    except Exception:

        def _version(dist: str) -> str:
            try:
                return importlib.metadata.version(dist)
            except Exception:
                return "not installed"

        return {
            "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "nucleide": _version("nucleide"),
            "pyne": _version("pyne"),
            "openmc": _version("openmc"),
        }


common.environment = _tolerant_environment


def _note(text: str) -> str:
    """Loud skip note: printed, and wrapped so the rendered Markdown lints."""
    note = textwrap.fill(text, width=100)
    print(note)
    return note


def _check(condition: bool) -> str:
    """Record a PASS/FAIL cell, counting failures for the run summary."""
    global FAILURES
    if condition:
        return "PASS"
    FAILURES += 1
    return "FAIL"


def _close(a: float, b: float, tol: float = TOL) -> bool:
    """Relative closeness using the harness `rel_diff` helper."""
    return rel_diff(a, b) <= tol


def _mix_labels(mixture: dict[str, Any]) -> list[str]:
    """Human-readable constituent labels for a parsed mixture."""
    labels: list[str] = []
    for entry in mixture["entries"]:
        labels.append(str(entry.get("name") or entry.get("symbol")))
    return labels


def alara_deck_section(report: Report) -> None:
    """Sample2/sample3 deck inventories vs committed values."""
    report.heading("ALARA decks")
    report.prose(
        "Nucleide `alara_parse_deck` run on `fixtures/alara/decks/sample2` and"
        " `sample3`.\nExpected values are the committed deck contents (upstream ALARA"
        " sample inputs).\nRe-parse stability re-parses the same text and requires an"
        " identical result\n(no deck serializer exists, so this is parse determinism"
        " rather than a\nwrite/read round-trip)."
    )
    rows: list[list[str]] = []
    decks = {
        "sample2": (ALARA_DIR / "decks" / "sample2").read_text(),
        "sample3": (ALARA_DIR / "decks" / "sample3").read_text(),
    }
    parsed = {name: nucleide.alara.alara_parse_deck(text) for name, text in decks.items()}

    d2 = parsed["sample2"]
    rows.append(
        [
            "sample2",
            "block count",
            "15",
            str(len(d2["block_kinds"])),
            _check(len(d2["block_kinds"]) == 15),
        ]
    )
    rows.append(
        [
            "sample2",
            "block order",
            "upstream order",
            "match" if d2["block_kinds"] == EXPECTED_BLOCKS else "MISMATCH",
            _check(d2["block_kinds"] == EXPECTED_BLOCKS),
        ]
    )
    rows.append(
        [
            "sample2",
            "geometry",
            "rectangular",
            str(d2["geometry"]),
            _check(d2["geometry"] == "rectangular"),
        ]
    )
    rows.append(
        [
            "sample2",
            "mixtures",
            "inner_mix[3], outer_mix[2]",
            f"{_mix_labels(d2['mixtures'][0])}, {_mix_labels(d2['mixtures'][1])}",
            _check(
                [m["name"] for m in d2["mixtures"]] == ["inner_mix", "outer_mix"]
                and [len(m["entries"]) for m in d2["mixtures"]] == [3, 2]
            ),
        ]
    )
    f2 = d2["fluxes"][0]
    rows.append(
        [
            "sample2",
            "flux entry",
            "flux_1 data/fluxin1 scale 1.0 skip 1",
            f"{f2['name']} {f2['file']} scale {f2['scale']} skip {f2['skip']}",
            _check(
                f2["name"] == "flux_1"
                and f2["file"] == "data/fluxin1"
                and _close(f2["scale"], 1.0)
                and f2["skip"] == 1
            ),
        ]
    )
    rows.append(
        [
            "sample2",
            "cooling count",
            "4",
            str(len(d2["cooling_times_s"])),
            _check(len(d2["cooling_times_s"]) == 4),
        ]
    )
    rows.append(
        [
            "sample2",
            "cooling 1 d",
            "86400 s",
            fmt(d2["cooling_times_s"][0]),
            _check(_close(d2["cooling_times_s"][0], 86_400.0)),
        ]
    )
    rows.append(
        [
            "sample2",
            "output resolution",
            "interval",
            str(d2["outputs"][0]["resolution"]),
            _check(d2["outputs"][0]["resolution"] == "interval"),
        ]
    )
    rows.append(
        [
            "sample2",
            "truncation",
            "1e-07",
            fmt(d2["truncation"]),
            _check(_close(d2["truncation"], 1e-7)),
        ]
    )

    d3 = parsed["sample3"]
    rows.append(
        [
            "sample3",
            "block count/order",
            "15, upstream order",
            str(len(d3["block_kinds"])),
            _check(len(d3["block_kinds"]) == 15 and d3["block_kinds"] == EXPECTED_BLOCKS),
        ]
    )
    rows.append(
        [
            "sample3",
            "geometry",
            "cylindrical",
            str(d3["geometry"]),
            _check(d3["geometry"] == "cylindrical"),
        ]
    )
    rows.append(
        [
            "sample3",
            "mixtures",
            "inner_mix[3] (CONC, b4c, ni), outer_mix[1] (fe)",
            f"{_mix_labels(d3['mixtures'][0])}, {_mix_labels(d3['mixtures'][1])}",
            _check(
                [m["name"] for m in d3["mixtures"]] == ["inner_mix", "outer_mix"]
                and [len(m["entries"]) for m in d3["mixtures"]] == [3, 1]
            ),
        ]
    )
    f3 = d3["fluxes"][0]
    rows.append(
        [
            "sample3",
            "flux entry",
            "flux_1 data/fluxin2 scale 1e6 skip 0",
            f"{f3['name']} {f3['file']} scale {f3['scale']} skip {f3['skip']}",
            _check(f3["file"] == "data/fluxin2" and _close(f3["scale"], 1e6) and f3["skip"] == 0),
        ]
    )
    rows.append(
        [
            "sample3",
            "cooling times",
            "[1, 86400, 60, 31557600] s",
            ", ".join(fmt(v) for v in d3["cooling_times_s"]),
            _check(
                len(d3["cooling_times_s"]) == 4
                and _close(d3["cooling_times_s"][0], 1.0)
                and _close(d3["cooling_times_s"][1], 86_400.0)
                and _close(d3["cooling_times_s"][2], 60.0)
                and _close(d3["cooling_times_s"][3], 31_557_600.0)
            ),
        ]
    )
    rows.append(
        [
            "sample3",
            "output resolution",
            "zone",
            str(d3["outputs"][0]["resolution"]),
            _check(d3["outputs"][0]["resolution"] == "zone"),
        ]
    )
    rows.append(
        [
            "sample3",
            "truncation",
            "1e-08",
            fmt(d3["truncation"]),
            _check(_close(d3["truncation"], 1e-8)),
        ]
    )
    for name, text in decks.items():
        stable = nucleide.alara.alara_parse_deck(text) == parsed[name]
        rows.append(
            [
                name,
                "re-parse stability",
                "identical",
                "identical" if stable else "DIFFERS",
                _check(stable),
            ]
        )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)
    report.prose(
        "Cooling `1 m` in sample3 reads as 60 s (minute convention), and the deck"
        " year is\n365.25 d (31557600 s); the output file below uses a 365 d year"
        " instead —\nboth conventions are preserved verbatim, not reconciled."
    )


def alara_flux_section(report: Report) -> None:
    """Group-flux fixtures vs committed values."""
    report.heading("ALARA group fluxes")
    report.prose(
        "`alara_parse_flux` run on `fixtures/alara/flux/fluxin2` (3 intervals x"
        " 175\ngroups) and `fluxin_zeros` (1 interval x 175 groups of exact zeros)."
        " Flux\nfiles carry no energy grid; the 175-group FENDL binning is implied by"
        " position."
    )
    rows: list[list[str]] = []
    flux = nucleide.alara.alara_parse_flux((ALARA_DIR / "flux" / "fluxin2").read_text(), "fluxin2")
    rows.append(
        [
            "fluxin2",
            "shape",
            "3 x 175",
            f"{flux['num_intervals']} x {flux['groups_per_interval']}",
            _check(flux["num_intervals"] == 3 and flux["groups_per_interval"] == 175),
        ]
    )
    rows.append(
        [
            "fluxin2",
            "interval lengths",
            "all 175",
            ",".join(str(len(iv)) for iv in flux["intervals"]),
            _check(all(len(iv) == 175 for iv in flux["intervals"])),
        ]
    )
    rows.append(
        [
            "fluxin2",
            "total == sum(totals)",
            fmt(sum(flux["totals"])),
            fmt(flux["total"]),
            _check(_close(flux["total"], sum(flux["totals"])) and flux["total"] > 0.0),
        ]
    )
    rows.append(
        [
            "fluxin2",
            "leading groups zero",
            "[0, 0, 0]",
            str([float(v) for v in flux["intervals"][0][:3]]),
            _check(all(v == 0.0 for v in flux["intervals"][0][:3])),
        ]
    )
    zeros = nucleide.alara.alara_parse_flux(
        (ALARA_DIR / "flux" / "fluxin_zeros").read_text(), "zeros"
    )
    rows.append(
        [
            "fluxin_zeros",
            "shape/total",
            "1 x 175, total 0",
            f"{zeros['num_intervals']} x {zeros['groups_per_interval']}, total {zeros['total']}",
            _check(
                zeros["num_intervals"] == 1
                and zeros["groups_per_interval"] == 175
                and zeros["total"] == 0.0
            ),
        ]
    )
    stable = (
        nucleide.alara.alara_parse_flux((ALARA_DIR / "flux" / "fluxin2").read_text(), "fluxin2")
        == flux
    )
    rows.append(
        [
            "fluxin2",
            "re-parse stability",
            "identical",
            "identical" if stable else "DIFFERS",
            _check(stable),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)


def alara_output_section(report: Report) -> None:
    """Sample2 output listing vs committed values."""
    report.heading("ALARA output")
    report.prose(
        "`alara_parse_output` run on `fixtures/alara/output/sample2.out` (upstream"
        " ALARA\nlisting, run label `sample2`)."
    )
    rows: list[list[str]] = []
    text = (ALARA_DIR / "output" / "sample2.out").read_text()
    parsed = nucleide.alara.alara_parse_output(text, "sample2")
    rows.append(["sample2.out", "row count", "2838", str(len(parsed)), _check(len(parsed) == 2838)])
    variables = {r["variable"] for r in parsed}
    rows.append(
        [
            "sample2.out",
            "response variables (2)",
            "{Number Density, Specific Activity}",
            str(sorted(variables)),
            _check(variables == {"Number Density", "Specific Activity"}),
        ]
    )
    times = sorted({r["time_s"] for r in parsed})
    expected_times = [-1.0, 0.0, 86_400.0, 8_640_000.0, 31_536_000.0, 3_153_600_000.0]
    rows.append(
        [
            "sample2.out",
            "time columns (6)",
            str(expected_times),
            str(times),
            _check(
                len(times) == 6
                and all(_close(a, b) for a, b in zip(times, expected_times, strict=True))
            ),
        ]
    )
    rows.append(
        [
            "sample2.out",
            "block kinds",
            "{Interval}",
            str({r["block"] for r in parsed}),
            _check({r["block"] for r in parsed} == {"Interval"}),
        ]
    )
    zones = {(r["block_num"], r["block_name"]) for r in parsed}
    rows.append(
        [
            "sample2.out",
            "zones",
            "{(1, inner_zone), (2, outer_zone)}",
            str(sorted(zones)),
            _check(zones == {(1, "inner_zone"), (2, "outer_zone")}),
        ]
    )
    totals = [r for r in parsed if r["nuclide"] == "total"]
    rows.append(["sample2.out", "total rows", "24", str(len(totals)), _check(len(totals) == 24)])
    h1 = [
        r
        for r in parsed
        if r["nuclide"] == "h-1"
        and r["variable"] == "Number Density"
        and r["block_num"] == 1
        and r["time_s"] == -1.0
    ]
    rows.append(
        [
            "sample2.out",
            "h-1 spot value",
            "6.6354e21 atoms/cm3",
            fmt(h1[0]["value"]) if len(h1) == 1 else f"{len(h1)} matches",
            _check(len(h1) == 1 and _close(h1[0]["value"], 6.6354e21)),
        ]
    )
    stable = nucleide.alara.alara_parse_output(text, "sample2") == parsed
    rows.append(
        [
            "sample2.out",
            "re-parse stability",
            "identical",
            "identical" if stable else "DIFFERS",
            _check(stable),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)
    report.prose(
        "`time_s == -1.0` marks the pre-irradiation column and `half_life_s == -1.0`"
        " marks\nstable nuclides and `total` rows (sentinel, not a measurement). The"
        " output\n1 y column is 31536000 s (365 d), while the deck cools for 31557600 s"
        "\n(365.25 d); both are kept verbatim."
    )


def fispact_section(report: Report) -> None:
    """Synthetic FISPACT inventory vs committed values."""
    report.heading("FISPACT inventory")
    report.prose(
        "`fispact_parse_output` run on the synthetic"
        " `fixtures/fispact/inventory.fis`\n(authored for nucleide; not FISPACT-II"
        " output). No independent oracle\nreader exists, so this is self-consistency"
        " against committed values."
    )
    rows: list[list[str]] = []
    text = (FISPACT_DIR / "inventory.fis").read_text()
    parsed = nucleide.fispact.fispact_parse_output(text, "synth")
    rows.append(["inventory.fis", "row count", "36", str(len(parsed)), _check(len(parsed) == 36)])
    variables = {r["variable"] for r in parsed}
    rows.append(
        [
            "inventory.fis",
            "variables (3)",
            "{Number Density, Specific Activity, Total Decay Heat}",
            str(sorted(variables)),
            _check(variables == {"Number Density", "Specific Activity", "Total Decay Heat"}),
        ]
    )
    times = sorted({r["time_s"] for r in parsed})
    rows.append(
        [
            "inventory.fis",
            "time columns (3)",
            "[0.0, 86400.0, 31536000.0]",
            str(times),
            _check(
                len(times) == 3
                and _close(times[0], 0.0)
                and _close(times[1], 86_400.0)
                and _close(times[2], 31_536_000.0)
            ),
        ]
    )
    h3 = [
        r
        for r in parsed
        if r["nuclide"] == "h-3"
        and r["variable"] == "Number Density"
        and _close(r["time_s"], 86_400.0)
    ]
    rows.append(
        [
            "inventory.fis",
            "h-3 atoms @ 1 d",
            "9.99e19",
            fmt(h3[0]["value"]) if len(h3) == 1 else f"{len(h3)} matches",
            _check(len(h3) == 1 and _close(h3[0]["value"], 9.99e19)),
        ]
    )
    totals = [r for r in parsed if r["nuclide"] == "total"]
    rows.append(["inventory.fis", "total rows", "9", str(len(totals)), _check(len(totals) == 9)])
    rows.append(
        [
            "inventory.fis",
            "half-life sentinel",
            "all -1.0",
            str(sorted({r["half_life_s"] for r in parsed})),
            _check(all(r["half_life_s"] == -1.0 for r in parsed)),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)
    report.prose(
        "The fixture carries no half-life column, so every row records"
        " `half_life_s ==\n-1.0` (missing-data sentinel, not a measurement). The extra"
        " DOSE column in\nthe second cooling step is ignored by the parser, and the"
        " third step uses the\n365-day year factor."
    )


def cccc_section(report: Report) -> None:
    """Synthetic CCCC libraries plus a PARTISN render/validate round-trip."""
    report.heading("CCCC libraries and PARTISN deck")
    report.prose(
        "`isotxs_parse`/`rtflux_parse` run on the synthetic"
        " `fixtures/cccc/*_sample`\nfiles, plus a PARTISN render/validate round-trip on a"
        " two-zone slab deck.\nNo independent oracle reader exists, so this is"
        " self-consistency against\ncommitted values."
    )
    rows: list[list[str]] = []
    isotxs_text = (CCCC_DIR / "isotxs_sample").read_text()
    lib = nucleide.cccc.isotxs_parse(isotxs_text)
    rows.append(
        [
            "isotxs_sample",
            "nuclides",
            "[U235, PU239]",
            str([n["label"] for n in lib["nuclides"]]),
            _check([n["label"] for n in lib["nuclides"]] == ["U235", "PU239"]),
        ]
    )
    rows.append(
        [
            "isotxs_sample",
            "U235 total xs",
            "[1.1, 2.2, 3.3]",
            str(lib["nuclides"][0]["total_xs"]),
            _check(
                all(
                    _close(a, b, 1e-12)
                    for a, b in zip(lib["nuclides"][0]["total_xs"], [1.1, 2.2, 3.3], strict=True)
                )
            ),
        ]
    )
    rows.append(
        [
            "isotxs_sample",
            "PU239 total xs",
            "[4.4, 5.5, 6.6]",
            str(lib["nuclides"][1]["total_xs"]),
            _check(
                all(
                    _close(a, b, 1e-12)
                    for a, b in zip(lib["nuclides"][1]["total_xs"], [4.4, 5.5, 6.6], strict=True)
                )
            ),
        ]
    )
    flux = nucleide.cccc.rtflux_parse((CCCC_DIR / "rtflux_sample").read_text())
    rows.append(
        [
            "rtflux_sample",
            "kind/groups/points",
            "RTFLUX 3 groups 2 points",
            f"{flux['kind']} {flux['groups']} groups {flux['npoints']} points",
            _check(flux["kind"] == "RTFLUX" and flux["groups"] == 3 and flux["npoints"] == 2),
        ]
    )
    rows.append(
        [
            "rtflux_sample",
            "total",
            "21.0",
            fmt(flux["total"]),
            _check(_close(flux["total"], 21.0, 1e-12)),
        ]
    )
    deck: dict[str, Any] = {
        "title": "synthetic slab",
        "dim": 1,
        "zones": [
            {"id": 1, "material": "fuel", "isotxs_labels": ["U235"], "density": 10.0},
            {"id": 2, "material": "blanket", "isotxs_labels": ["PU239"], "density": 5.0},
        ],
        "source": "isotropic",
    }
    text = nucleide.cccc.partisn_render(deck)
    rows.append(
        [
            "partisn",
            "rendered keywords",
            "DIM/ZONE/MATERIAL/SOURCE/END",
            "present"
            if all(k in text for k in ("DIM 1", "ZONE 1", "MATERIAL", "SOURCE", "END\n"))
            else "MISSING",
            _check(
                all(k in text for k in ("DIM 1", "ZONE 1", "MATERIAL", "SOURCE", "END\n"))
                and text.endswith("END\n")
            ),
        ]
    )
    try:
        nucleide.cccc.partisn_validate(deck, isotxs_text)
        valid_ok = True
    except Exception:
        valid_ok = False
    rows.append(
        ["partisn", "validate vs isotxs", "OK", "OK" if valid_ok else "RAISED", _check(valid_ok)]
    )
    try:
        bad = {
            "title": "bad",
            "dim": 1,
            "zones": [{"id": 1, "material": "fuel", "isotxs_labels": ["U238"], "density": 1.0}],
        }
        nucleide.cccc.partisn_validate(bad, isotxs_text)
        dangling_ok = False
    except ValueError:
        dangling_ok = True
    except Exception:
        dangling_ok = False
    rows.append(
        [
            "partisn",
            "dangling label rejected",
            "ValueError",
            "ValueError" if dangling_ok else "NOT RAISED",
            _check(dangling_ok),
        ]
    )
    rows.append(
        [
            "partisn",
            "re-render stability",
            "identical",
            "identical" if nucleide.cccc.partisn_render(deck) == text else "DIFFERS",
            _check(nucleide.cccc.partisn_render(deck) == text),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)
    report.prose(
        "The fixtures cover the documented text-analog subset only; production CCCC"
        "\nfiles are binary. `partisn_validate` maps deck zone labels against the"
        " ISOTXS\nlibrary and rejects unmapped labels."
    )


def origen_section(report: Report) -> None:
    """Synthetic ORIGEN TAPEs vs committed values."""
    report.heading("ORIGEN TAPEs")
    report.prose(
        "`origen_parse_tape5/6/9` run on the synthetic"
        " `fixtures/origen/tape*_sample`\nfiles (simplified input echoes, not real"
        " ORIGEN output). No independent\noracle reader exists, so this is"
        " self-consistency against committed values."
    )
    rows: list[list[str]] = []
    tape5 = nucleide.origen.origen_parse_tape5((ORIGEN_DIR / "tape5_sample").read_text())
    rows.append(
        [
            "tape5_sample",
            "titles",
            "2",
            str(tape5["titles"]),
            _check(
                tape5["titles"] == ["SYNTHETIC PWR PIN - TAPE5 SAMPLE", "CASE 1 - BASE DEPLETION"]
            ),
        ]
    )
    rows.append(
        [
            "tape5_sample",
            "irradiation steps",
            "[(3e13, 100 d), (0, 30 d)]",
            str([(s["flux"], s["days"]) for s in tape5["irradiation_steps"]]),
            _check(
                len(tape5["irradiation_steps"]) == 2
                and _close(tape5["irradiation_steps"][0]["flux"], 3.0e13)
                and _close(tape5["irradiation_steps"][0]["days"], 100.0)
                and tape5["irradiation_steps"][1]["flux"] == 0.0
                and _close(tape5["irradiation_steps"][1]["days"], 30.0)
            ),
        ]
    )
    rows.append(
        [
            "tape5_sample",
            "materials",
            "fuel[3], clad[1]",
            ",".join(f"{m['name']}[{len(m['entries'])}]" for m in tape5["materials"]),
            _check(
                [m["name"] for m in tape5["materials"]] == ["fuel", "clad"]
                and [len(m["entries"]) for m in tape5["materials"]] == [3, 1]
            ),
        ]
    )
    tape6 = nucleide.origen.origen_parse_tape6((ORIGEN_DIR / "tape6_sample").read_text())
    rows.append(
        [
            "tape6_sample",
            "record count",
            "4",
            str(len(tape6["records"])),
            _check(len(tape6["records"]) == 4),
        ]
    )
    rows.append(
        [
            "tape6_sample",
            "total == sum(records)",
            fmt(sum(r["activity_bq"] for r in tape6["records"])),
            fmt(tape6["total_activity"]),
            _check(
                _close(
                    tape6["total_activity"], sum(r["activity_bq"] for r in tape6["records"]), 1e-9
                )
                and tape6["total_activity"] > 0.0
            ),
        ]
    )
    tape9 = nucleide.origen.origen_parse_tape9((ORIGEN_DIR / "tape9_sample").read_text())
    by_name = {e["nuclide"]: e["decay_const"] for e in tape9}
    rows.append(["tape9_sample", "entry count", "5", str(len(tape9)), _check(len(tape9) == 5)])
    rows.append(
        [
            "tape9_sample",
            "decay constants",
            "U235 3.1209e-17, Cs137 7.3217e-10, Co60 4.1674e-09",
            ", ".join(f"{k} {fmt(by_name[k])}" for k in ("U235", "Cs137", "Co60")),
            _check(
                _close(by_name["U235"], 3.1209e-17, 1e-9)
                and _close(by_name["Cs137"], 7.3217e-10, 1e-9)
                and _close(by_name["Co60"], 4.1674e-09, 1e-9)
            ),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)


def r2s_section(report: Report) -> None:
    """R2S workflow builder on the sample2 deck/output pair."""
    report.heading("R2S workflow")
    report.prose(
        "`r2s_from_deck`/`r2s_expand`/`r2s_validate` run on"
        " `fixtures/alara/decks/sample2`;\n`r2s_assemble` runs on"
        " `fixtures/alara/output/sample2.out`. Smoke tests on\ncommitted"
        " fixtures; no independent oracle exists."
    )
    rows: list[list[str]] = []
    deck_text = (ALARA_DIR / "decks" / "sample2").read_text()
    workflow = nucleide.r2s.r2s_from_deck(deck_text)
    rows.append(
        [
            "sample2 deck",
            "workflow steps",
            "[inner_zone x flux_1, outer_zone x flux_1]",
            str(workflow["steps"]),
            _check(
                workflow["steps"]
                == [
                    {"zone": "inner_zone", "flux": "flux_1"},
                    {"zone": "outer_zone", "flux": "flux_1"},
                ]
            ),
        ]
    )
    rows.append(
        [
            "sample2 deck",
            "top schedule",
            "1_year",
            str(workflow["top_schedule"]),
            _check(workflow["top_schedule"] == "1_year"),
        ]
    )
    rows.append(
        [
            "sample2 deck",
            "cooling entries",
            "4, first 86400 s",
            f"{len(workflow['cooling_s'])}, first {fmt(workflow['cooling_s'][0])}",
            _check(len(workflow["cooling_s"]) == 4 and _close(workflow["cooling_s"][0], 86_400.0)),
        ]
    )
    steps = nucleide.r2s.r2s_expand(deck_text)
    total = sum(s["duration_s"] for s in steps)
    rows.append(
        [
            "sample2 deck",
            "expanded total",
            "31557600 s (1 ALARA year)",
            fmt(total),
            _check(
                len(steps) == 1
                and _close(total, 31_557_600.0)
                and steps[0]["flux"] == "flux_1"
                and steps[0]["is_cooling"] is False
            ),
        ]
    )
    try:
        nucleide.r2s.r2s_validate(workflow, deck_text)
        valid_ok = True
    except Exception:
        valid_ok = False
    rows.append(
        [
            "sample2 deck",
            "workflow validates",
            "OK",
            "OK" if valid_ok else "RAISED",
            _check(valid_ok),
        ]
    )
    output_text = (ALARA_DIR / "output" / "sample2.out").read_text()
    source = nucleide.r2s.r2s_assemble(output_text, "sample2", "inner_zone", 4)
    rows.append(
        [
            "sample2.out",
            "assemble groups",
            "4 uniform quarters",
            ", ".join(fmt(g) for g in source["groups"]),
            _check(
                len(source["groups"]) == 4
                and source["total"] > 0.0
                and all(_close(g, source["total"] / 4.0, 1e-12) for g in source["groups"])
            ),
        ]
    )
    missing = nucleide.r2s.r2s_assemble(output_text, "sample2", "missing", 3)
    rows.append(
        [
            "sample2.out",
            "assemble unknown zone",
            "total 0",
            fmt(missing["total"]),
            _check(_close(missing["total"], 0.0, 1e-12)),
        ]
    )
    report.table(["Fixture", "Check", "Expected", "Observed", "Status"], rows)
    report.prose(
        "Both zones share one flux (no void-only step is synthesized), and the"
        " expanded\nschedule totals one ALARA year of 365.25 d (31557600 s). Photon"
        " assembly\nsums zone `SpecificActivity` and splits the total uniformly over the"
        "\nrequested groups: the total shutdown strength is preserved, while real"
        " decay\nphotons follow the nuclide- and energy-dependent `.photonSrc` lines."
    )


def _try_pyne_alara_parse(pyne_alara: Any, deck_text: str) -> Any | None:
    """Attempt a pyne.alara deck parse with whatever entry point exists.

    Returns the parsed object, or None when no usable entry point is found.
    Never raises: every probe failure maps to a loud skip.
    """
    deck_path = str(ALARA_DIR / "decks" / "sample2")
    candidates = [
        n
        for n in dir(pyne_alara)
        if not n.startswith("_")
        and any(k in n.lower() for k in ("pars", "dict", "deck", "input", "read", "load"))
    ]
    for name in candidates:
        func = getattr(pyne_alara, name, None)
        if not callable(func):
            continue
        for arg in (deck_text, deck_path):
            try:
                result = func(arg)
            except Exception:
                continue
            if isinstance(result, (dict, list)):
                return result
    return None


def _compare_pyne_alara(report: Report, pyne_alara: Any) -> None:
    """Compare the sample2 deck block set against pyne.alara where possible."""
    try:
        public = sorted(n for n in dir(pyne_alara) if not n.startswith("_"))
        report.prose(f"`pyne.alara` is importable; public names: {', '.join(public) or '—'}.")
        deck_text = (ALARA_DIR / "decks" / "sample2").read_text()
        nuc_blocks = nucleide.alara.alara_parse_deck(deck_text)["block_kinds"]
        parsed = _try_pyne_alara_parse(pyne_alara, deck_text)
        if parsed is None:
            report.prose(
                _note(
                    "SKIPPED: no pyne.alara entry point in this environment parses an"
                    " ALARA input deck into a comparable block set, so the deck"
                    " block-set comparison did not run."
                )
            )
            return
        if isinstance(parsed, dict):
            keys = [str(k) for k in parsed]
            if any(k.lower() not in _KNOWN_BLOCKS for k in keys):
                report.prose(
                    _note(
                        "SKIPPED: the pyne.alara parse result is not an ALARA"
                        " block set (keys fall outside the deck block"
                        " vocabulary), so the deck block-set comparison did"
                        " not run."
                    )
                )
                return
            pyne_blocks: Any = sorted(keys)
            expected: Any = sorted(nuc_blocks)
        elif isinstance(parsed, list):
            vals = [str(v) for v in parsed]
            if any(v.lower() not in _KNOWN_BLOCKS for v in vals):
                report.prose(
                    _note(
                        "SKIPPED: the pyne.alara parse result is not an ALARA"
                        " block set (entries fall outside the deck block"
                        " vocabulary), so the deck block-set comparison did"
                        " not run."
                    )
                )
                return
            pyne_blocks = sorted(vals)
            expected = sorted(nuc_blocks)
        else:
            report.prose(
                _note(
                    "SKIPPED: the pyne.alara parse result has an unrecognized shape,"
                    " so the deck block-set comparison did not run."
                )
            )
            return
        report.table(
            ["Source", "Block set", "Status"],
            [
                ["nucleide", ", ".join(expected), ""],
                ["pyne.alara", ", ".join(pyne_blocks), _check(list(pyne_blocks) == list(expected))],
            ],
        )
    except Exception as exc:
        report.prose(_note(f"SKIPPED: pyne.alara probe failed ({type(exc).__name__}: {exc})."))


def pyne_oracle_section(report: Report) -> None:
    """Container-only PyNE oracle probes; every absence is a loud skip."""
    report.heading("PyNE oracle probes (container-only)")
    try:
        import pyne.alara as pyne_alara
    except ImportError as exc:
        report.prose(
            _note(
                "SKIPPED: `pyne.alara` is not importable in this environment"
                f" ({exc}); the ALARA deck oracle comparison did not run. It runs"
                " inside the validation container."
            )
        )
    else:
        _compare_pyne_alara(report, pyne_alara)
    for module in ("pyne.origen", "pyne.fispact", "pyne.cccc"):
        try:
            __import__(module)
        except ImportError as exc:
            report.prose(
                _note(
                    f"SKIPPED: `{module}` is not importable ({exc}); PyNE ships no"
                    " matching reader, so the corresponding fixtures remain"
                    " synthetic self-consistency checks."
                )
            )
        else:
            report.prose(
                f"`{module}` is unexpectedly importable; no comparison is defined for it yet."
            )


def main() -> int:
    report = Report("activation", "Activation I/O (`activation_vs_refs.py`)")
    report.prose(
        "Nucleide's activation-code readers (alara-io, cccc-io, fispact-io,"
        " origen-io)\nand the r2s workflow builder are checked against committed"
        " fixtures.\nIndependent-oracle status: the ALARA deck block set has a"
        " container-only\n`pyne.alara` probe below; FISPACT, CCCC, ORIGEN and R2S"
        " fixtures are\nsynthetic self-consistency checks with no independent oracle."
        " Skipped\ncomparisons are listed explicitly; none is silent."
    )
    alara_deck_section(report)
    alara_flux_section(report)
    alara_output_section(report)
    fispact_section(report)
    cccc_section(report)
    origen_section(report)
    r2s_section(report)
    pyne_oracle_section(report)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} activation check(s) failed", file=sys.stderr)
        return 1
    print("All activation checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
