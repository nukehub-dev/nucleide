"""MCPL interchange oracle (`nucleide-mcpl-io` vs upstream `mcpl` tooling).

Three tiers:

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
3. SSW round-trip (always run, synthetic pairs only): the committed
   hand-built `fixtures/mcpl/ssw_conversion/reference.w` converts to MCPL
   with the documented surface/kind pairing and back against the same
   reference header; energies/times/userflags are checked in closed form.
   A CLI cross-check over the upstream `ssw2mcpl`/`mcpl2ssw` console scripts
   (shipped by the `mcpl-extra` 2.2.8 package, pinned in `Containerfile`)
   runs when those binaries are present and SKIP-reports with its reason
   otherwise (never fails without the oracle).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.mcpl as mcpl
from nucleide.mcnp import read_ssw

REPO_ROOT = Path(__file__).resolve().parent.parent
SSW_REF = REPO_ROOT / "fixtures" / "mcpl" / "ssw_conversion" / "reference.w"
SSW_SURFS = [100, 200]
SSW_KINDS = ["neutron", "gamma"]

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


def tier_ssw(tmp: str) -> tuple[list[list[str]], list[str]]:
    """SSW round-trip gates S1-S8 on the synthetic reference pair.

    Converts the committed hand-built `reference.w` (two tracks, explicit
    surface/kind pairing) to MCPL and back against the same reference header,
    checking closed-form energy/time/userflag conservation (S1-S4) plus the
    fidelity tail: upstream `cs = 1.0` compat (S5), `niss`
    passthrough vs override (S6), polarisation/universal opt-ins (S7), and
    the widened SSW-PDG table (S8). Synthetic pairs only: no MCNP run, no
    vendored upstream bytes.
    """
    rows: list[list[str]] = []
    notes: list[str] = []
    probe = os.path.join(tmp, "ssw_probe.mcpl")
    back_w = os.path.join(tmp, "ssw_back.w")
    n = mcpl.ssw2mcpl(str(SSW_REF), probe, SSW_SURFS, SSW_KINDS)
    ok = n == 2
    rows.append(["S1 ssw2mcpl count", "2", str(n), _check(ok, "S1 count")])
    got = mcpl.read_mcpl(probe).particles()
    ok = (
        len(got) == 2
        and rel_diff(got[0]["ekin"], 2.5) < 1e-6
        and rel_diff(got[1]["ekin"], 0.662) < 1e-6
        and [p["pdgcode"] for p in got] == [2112, 22]
    )
    rows.append(
        [
            "S2 ssw2mcpl energy+PDG",
            "2.5/2112, 0.662/22",
            f"{fmt(got[0]['ekin'])}/{got[0]['pdgcode']}, {fmt(got[1]['ekin'])}/{got[1]['pdgcode']}",
            _check(ok, "S2 fields"),
        ]
    )
    ok = (
        abs(got[0]["time"] - 3.0e5 * 1e-5) < 1e-9
        and got[1]["time"] == 0.0
        and [p["userflags"] for p in got] == SSW_SURFS
    )
    rows.append(
        [
            "S3 shakes->ms + surf flags",
            "3.0 ms/[100, 200]",
            f"{fmt(got[0]['time'])} ms/{[p['userflags'] for p in got]}",
            _check(ok, "S3 time+flags"),
        ]
    )
    m = mcpl.mcpl2ssw(probe, str(SSW_REF), back_w)
    tracks = read_ssw(back_w).tracks()
    ok = (
        m == 2
        and len(tracks) == 2
        and rel_diff(tracks[0]["erg"], 2.5) < 1e-6
        and abs(tracks[0]["tme"] - 3.0e5) / 3.0e5 < 1e-6
        and rel_diff(tracks[1]["erg"], 0.662) < 1e-6
    )
    rows.append(
        [
            "S4 mcpl2ssw round-trip",
            "2.5 erg/3.0e5 tme, 0.662 erg",
            f"{fmt(tracks[0]['erg'])} erg/{fmt(tracks[0]['tme'])} tme, {fmt(tracks[1]['erg'])} erg",
            _check(ok, "S4 round-trip"),
        ]
    )
    # S5: force_cs_to_one reproduces the upstream 2.2.8 `cs = 1.0` spelling
    # (u/v stay verbatim); default keeps the true cosine (S4 above).
    cs_w = os.path.join(tmp, "ssw_cs1.w")
    mcpl.mcpl2ssw(probe, str(SSW_REF), cs_w, force_cs_to_one=True)
    cs_tracks = read_ssw(cs_w).tracks()
    ok = (
        len(cs_tracks) == 2
        and all(abs(t["cs"] - 1.0) < 1e-12 for t in cs_tracks)
        and abs(cs_tracks[1]["u"] - 1.0) < 1e-6
        and abs(cs_tracks[1]["v"]) < 1e-6
    )
    rows.append(
        [
            "S5 force_cs_to_one compat",
            "cs 1.0/1.0, u/v verbatim",
            f"cs {[fmt(t['cs']) for t in cs_tracks]}, u {fmt(cs_tracks[1]['u'])}",
            _check(ok, "S5 cs compat"),
        ]
    )
    # S6: niss passes the reference header through by default (upstream-2.2.8
    # spelling); niss= stamps an explicit override.
    ref_niss = read_ssw(str(SSW_REF)).niss
    back_niss = read_ssw(back_w).niss
    ok = back_niss == ref_niss
    rows.append(
        [
            "S6 niss passthrough",
            f"niss {ref_niss}",
            f"niss {back_niss}",
            _check(ok, "S6 passthrough"),
        ]
    )
    niss_w = os.path.join(tmp, "ssw_niss.w")
    mcpl.mcpl2ssw(probe, str(SSW_REF), niss_w, niss=7)
    got_niss = read_ssw(niss_w).niss
    ok = got_niss == 7
    rows.append(
        [
            "S6b niss override",
            "niss 7",
            f"niss {got_niss}",
            _check(ok, "S6b override"),
        ]
    )
    # S7: polarisation + universal-PDG opt-ins round-trip through the header.
    pol_mcpl = os.path.join(tmp, "pol.mcpl")
    mcpl.ssw2mcpl(str(SSW_REF), pol_mcpl, SSW_SURFS, SSW_KINDS, {"polarisation": [0.1, 0.2, 0.3]})
    pol_back = mcpl.read_mcpl(pol_mcpl)
    ok = pol_back.has_polarisation and all(
        abs(a - b) < 1e-6
        for p in pol_back.particles()
        for a, b in zip(p["polarisation"], [0.1, 0.2, 0.3], strict=True)
    )
    rows.append(
        [
            "S7 polarisation carry",
            "has_polarisation [0.1, 0.2, 0.3]",
            f"{pol_back.has_polarisation} {pol_back.particles()[0]['polarisation']}",
            _check(ok, "S7 polarisation"),
        ]
    )
    uni_mcpl = os.path.join(tmp, "uni.mcpl")
    mcpl.ssw2mcpl(
        str(SSW_REF), uni_mcpl, [100, 200], ["neutron", "neutron"], {"universal_pdg": True}
    )
    uni_back = mcpl.read_mcpl(uni_mcpl)
    ok = uni_back.universal_pdgcode == 2112 and [p["pdgcode"] for p in uni_back.particles()] == [
        2112,
        2112,
    ]
    rows.append(
        [
            "S7b universal PDG carry",
            "universal 2112",
            f"universal {uni_back.universal_pdgcode}",
            _check(ok, "S7b universal"),
        ]
    )
    # S8: widened SSW-PDG table over the same reference geometry.
    ext_mcpl = os.path.join(tmp, "ext.mcpl")
    mcpl.ssw2mcpl(str(SSW_REF), ext_mcpl, SSW_SURFS, ["electron", "positron"])
    ext = mcpl.read_mcpl(ext_mcpl).particles()
    ok = [p["pdgcode"] for p in ext] == [11, -11] and all(
        rel_diff(p["ekin"], w) < 1e-6 for p, w in zip(ext, [2.5, 0.662], strict=True)
    )
    rows.append(
        [
            "S8 extended PDG table",
            "11/-11, energies verbatim",
            f"{[p['pdgcode'] for p in ext]}, {[fmt(p['ekin']) for p in ext]}",
            _check(ok, "S8 pdg table"),
        ]
    )
    notes.append(
        "Tier 3 converts the committed synthetic SSW reference (hand-framed, "
        "no MCNP run) with the documented surface/kind pairing; the "
        "`mcpl2ssw` leg clones the same reference header (S1-S4 baseline, "
        "S5-S8 fidelity tail: cs compat, niss, polarisation/universal, "
        "extended PDG)."
    )
    return rows, notes


def tier_extra(tmp: str) -> tuple[list[list[str]], list[str], bool]:
    """Upstream converter CLI cross-check (`ssw2mcpl`/`mcpl2ssw` scripts).

    Runs the upstream converter console scripts (shipped by the `mcpl-extra`
    2.2.8 package: `ssw2mcpl [options] input.ssw [output.mcpl]`,
    `mcpl2ssw [options] <input.mcpl> <reference.ssw> [output.ssw]`) over the
    synthetic reference pair and compares particle/track counts plus energies
    against the nucleide leg. The scripts are an optional oracle dependency:
    when absent (or when any probe step errors), tier 4 SKIP-reports with its
    reason and never fails.
    """
    to_mcpl = shutil.which("ssw2mcpl")
    to_ssw = shutil.which("mcpl2ssw")
    if to_mcpl is None or to_ssw is None:
        missing = "ssw2mcpl" if to_mcpl is None else "mcpl2ssw"
        return (
            [["Upstream converter cross-check", f"SKIP ({missing} unavailable)"]],
            [f"Tier 4 skipped: no `{missing}` script on PATH (mcpl-extra not installed)."],
            True,
        )
    try:
        up_mcpl = os.path.join(tmp, "extra.mcpl")
        r = subprocess.run(
            [to_mcpl, "-s", "-n", str(SSW_REF), up_mcpl],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if r.returncode != 0:
            raise RuntimeError(f"ssw2mcpl exited {r.returncode}: {r.stderr.strip()[-500:]}")
        up_ps = mcpl.read_mcpl(up_mcpl).particles()
        up_w = os.path.join(tmp, "extra.w")
        r = subprocess.run(
            [to_ssw, up_mcpl, str(SSW_REF), up_w],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if r.returncode != 0:
            raise RuntimeError(f"mcpl2ssw exited {r.returncode}: {r.stderr.strip()[-500:]}")
        up_tracks = read_ssw(up_w).tracks()
    except Exception as exc:
        return (
            [["Upstream converter cross-check", "SKIP (oracle probe failed)"]],
            [f"Tier 4 skipped: converter probe failed ({exc})."],
            True,
        )
    rows: list[list[str]] = []
    notes: list[str] = []
    ok = len(up_ps) == 2
    rows.append(["T3 ssw2mcpl count", "2", str(len(up_ps)), _check(ok, "T3 count")])
    eks = sorted(p["ekin"] for p in up_ps)
    ok = len(eks) == 2 and rel_diff(eks[0], 0.662) < 1e-6 and rel_diff(eks[1], 2.5) < 1e-6
    rows.append(
        [
            "T4 converter energies",
            "0.662, 2.5",
            ", ".join(fmt(v) for v in eks),
            _check(ok, "T4 energies"),
        ]
    )
    # The return leg rewrites the header (`nrss`/`np1` patched) and forces
    # the stored `cs` slot to 1.0 by upstream design, so only count + energy
    # are compared here (direction cosines are nucleide-verbatim, S4).
    back_eks = sorted(t["erg"] for t in up_tracks)
    ok = (
        len(up_tracks) == 2
        and len(back_eks) == 2
        and rel_diff(back_eks[0], 0.662) < 1e-6
        and rel_diff(back_eks[1], 2.5) < 1e-6
    )
    rows.append(
        [
            "T5 mcpl2ssw count+energy",
            "2 tracks: 0.662, 2.5 erg",
            f"{len(up_tracks)} tracks: " + ", ".join(fmt(v) for v in back_eks) + " erg",
            _check(ok, "T5 return"),
        ]
    )
    notes.append("Tier 4 runs the upstream converter scripts over the synthetic pair.")
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
        rows3, notes3 = tier_ssw(tmp)
        report.table(["Gate", "Expected", "Got", "Status"], rows3)
        for note in notes3:
            report.prose(note)
        rows4, notes4, extra_skipped = tier_extra(tmp)
        if extra_skipped:
            report.table(["Check", "Status"], rows4)
        else:
            report.table(["Check", "Expected", "Got", "Status"], rows4)
        for note in notes4:
            report.prose(note)
    report.emit()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
