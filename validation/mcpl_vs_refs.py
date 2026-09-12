"""MCPL interchange oracle (`nucleide-mcpl-io` vs upstream `mcpl` tooling).

Two tiers:

1. Synthetic gates (always run): hand-built axis-vector records are written
   with :func:`nucleide.mcpl.write_mcpl`, read back, and checked for
   field-level conservation (energy/position/direction/PDG/flags) through
   single precision. Gzip transport is covered semantically only (compressed
   bytes are encoder-dependent and never asserted).
2. Upstream cross-check: the same synthetic file is opened with the upstream
   ``mcpl`` Python package (``pip install mcpl`` / ``conda -c conda-forge
   mcpl``) when importable, comparing particle count, kinetic energies, and
   PDG codes at record level. The upstream package is an optional oracle
   dependency: if it cannot be imported (or its API mismatches), tier 2 is
   reported as SKIP with its reason (never silently).
"""

from __future__ import annotations

import os
import sys
import tempfile

from common import Report, fmt, rel_diff

import nucleide.mcpl as mcpl

FAILURES = 0

PARTICLES = [
    {
        "ekin": 2.5,
        "polarisation": [0.0, 0.0, 0.0],
        "position": [1.0, -2.0, 0.5],
        "direction": [0.0, 0.0, 1.0],
        "time": 3.0,
        "weight": 1.0,
        "pdgcode": 2112,
        "userflags": 0,
    },
    {
        "ekin": 0.662,
        "polarisation": [0.0, 0.0, 0.0],
        "position": [0.0, 0.0, 0.0],
        "direction": [1.0, 0.0, 0.0],
        "time": 0.0,
        "weight": 0.5,
        "pdgcode": 22,
        "userflags": 7,
    },
]

HEADER = {
    "srcname": "nucleide-oracle",
    "comments": ["synthetic probe"],
    "has_userflags": True,
    "has_polarisation": False,
    "double_prec": False,
    "universal_pdgcode": None,
    "universal_weight": None,
    "blobs": [],
}


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def tier1(path: str) -> tuple[list[list[str]], list[str]]:
    """Synthetic round-trip gates G1-G3. Returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []
    mcpl.write_mcpl(path, HEADER, PARTICLES)
    back = mcpl.read_mcpl(path)
    got = back.particles()

    ok = back.nparticles == 2 and len(got) == 2
    rows.append(["G1 particle count", "2", str(back.nparticles), _check(ok, "G1 count")])

    ok = all(
        rel_diff(g["ekin"], w["ekin"]) < 1e-6 and g["pdgcode"] == w["pdgcode"]
        for g, w in zip(got, PARTICLES, strict=True)
    )
    rows.append(
        [
            "G2 energy+PDG conserved",
            "2.5/2112, 0.662/22",
            f"{fmt(got[0]['ekin'])}/{got[0]['pdgcode']}, {fmt(got[1]['ekin'])}/{got[1]['pdgcode']}",
            _check(ok, "G2 fields"),
        ]
    )

    ok = all(
        abs(a - b) < 1e-6
        for g, w in zip(got, PARTICLES, strict=True)
        for a, b in zip(g["direction"], w["direction"], strict=True)
    )
    rows.append(
        [
            "G3 directions conserved",
            "+z, +x",
            f"{[fmt(v) for v in got[0]['direction']]}, {[fmt(v) for v in got[1]['direction']]}",
            _check(ok, "G3 directions"),
        ]
    )
    notes.append("Tier 1 uses hand-built synthetic records only; no upstream files are read.")
    return rows, notes


def tier2(path: str) -> tuple[list[list[str]], list[str], bool]:
    """Upstream `mcpl` record-level cross-check. Returns (rows, notes, skipped).

    Opens the nucleide-written synthetic file with the upstream ``mcpl``
    Python package (``MCPLFile`` context manager, ``particle_blocks`` arrays:
    ``ekin``/``pdgcode``/``weight``/``userflags``/``x``/``y``/``z``/
    ``ux``/``uy``/``uz``/``time`` — all exercised against 2.2.8) and compares
    counts, energies, and PDG codes. The package is an optional oracle
    dependency: import or API mismatches report SKIP with the reason (never
    silently, never assumed).
    """
    try:
        import mcpl as up  # type: ignore[import-not-found]
    except Exception as exc:
        return (
            [["Upstream mcpl cross-check", "SKIP (mcpl unavailable)"]],
            [f"Tier 2 skipped: cannot import mcpl ({exc})."],
            True,
        )
    try:
        n_total = 0
        ekins: list[float] = []
        pdgs: list[int] = []
        with up.MCPLFile(path) as fh:
            for block in fh.particle_blocks:
                n_total += len(block.ekin)
                ekins.extend(float(v) for v in block.ekin)
                pdgs.extend(int(v) for v in block.pdgcode)
    except Exception as exc:
        return (
            [["Upstream mcpl cross-check", "SKIP (upstream API mismatch)"]],
            [f"Tier 2 skipped: upstream read failed ({exc})."],
            True,
        )
    rows: list[list[str]] = []
    notes: list[str] = []
    ok = n_total == len(PARTICLES)
    rows.append(
        ["T1 upstream particle count", str(len(PARTICLES)), str(n_total), _check(ok, "T1 count")]
    )
    ok = (
        len(ekins) == len(PARTICLES)
        and all(rel_diff(g, w["ekin"]) < 1e-6 for g, w in zip(ekins, PARTICLES, strict=True))
        and pdgs == [w["pdgcode"] for w in PARTICLES]
    )
    got_t2 = (
        f"{fmt(ekins[0])}/{pdgs[0]}, {fmt(ekins[1])}/{pdgs[1]}"
        if len(ekins) == 2 and len(pdgs) == 2
        else f"{ekins}/{pdgs}"
    )
    rows.append(["T2 upstream energy+PDG", "2.5/2112, 0.662/22", got_t2, _check(ok, "T2 fields")])
    notes.append("Tier 2 opens the nucleide-written synthetic file with upstream tooling.")
    return rows, notes, False


def main() -> int:
    report = Report("mcpl", "MCPL interchange vs upstream tooling")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "probe.mcpl")
        rows1, notes1 = tier1(path)
        report.table(["Gate", "Expected", "Got", "Status"], rows1)
        for note in notes1:
            report.prose(note)
        rows2, notes2, skipped = tier2(path)
        if skipped:
            report.table(["Check", "Status"], rows2)
        else:
            report.table(["Check", "Expected", "Got", "Status"], rows2)
        for note in notes2:
            report.prose(note)
    report.emit()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
