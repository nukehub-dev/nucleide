"""Damage/gas metric cross-check (`nucleide.damage` vs analytic gates + SPECTER).

Two parts:

1. Analytic gates (always run): G1 NRT piecewise + Lindhard partition against
   an independent in-script transcription of the published constants
   (Robinson fit; NRT 1975), G2 the arc-dpa efficiency construction
   (Nordlund et al. 2018 Eq. (7): junction value exactly 1, saturation at
   ``c_arc``, monotone decrease, high branch = NRT x xi), G3 spectral folds
   against hand sums (piecewise-constant per group; appm = dpa x 1e6), G4
   invariants (zero-flux groups contribute exactly 0; the He/dpa ratio at
   zero dpa is a loud named error, never ``inf``; malformed folds name their
   cause), G5 the UQ hook (pinned seed, k-standard-error gates against the
   exact expectation and the first-order propagation, determinism).
2. SPECTER oracle check (container): the SPECTER manual (Greenwood &
   Smither, ANL/FPP/TM-197, January 1985 — a US-government work, public
   domain) is downloaded once from OSTI into ``validation/.cache/`` (SHA-256
   pinned) and converted to text with ``pdftotext``; each transcribed spot
   below is verified to occur in the report text (OCR-tolerant matching —
   the scan's embedded text layer confuses 0/O, 5/S, 1/l), then folded
   through ``nucleide.damage`` and compared against the report's printed
   result at report print precision (5 significant digits; tolerance
   1e-4 relative ~ 2 half-ulp of the print). SPECTER is an oracle only:
   no table is vendored. Outside the container (no download, no
   ``pdftotext``) the check is a recorded SKIP, never a silent pass. SPECTER
   stays the oracle for this script; the crate's opt-in Table VII fallback
   (`nucleide.damage.specter_table`) is a separate caller input the folds
   never consult implicitly.

Units: fluence in n/cm², cross sections in barns (gas spots in millibarns as
printed), energies in MeV bounds / eV damage functions, exposure in seconds.
"""

from __future__ import annotations

import hashlib
import math
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.damage as dmg

FAILURES = 0
REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_DIR = Path(__file__).resolve().parent
CACHE_DIR = SCRIPT_DIR / ".cache"

# --- SPECTER provenance (Greenwood & Smither, ANL/FPP/TM-197, Jan 1985) ---
SPECTER_URL = "https://www.osti.gov/servlets/purl/6022143"
SPECTER_SHA256 = "489bca1482ffb33c160329283c0bbb2245841d2808458749ff431caf8d89ce44"
SPECTER_PDF = CACHE_DIR / "specter_anl_fpp_tm_197.pdf"
# HFIR-CTR32 irradiation printed in the report's Tables V/VI: total fluence
# 4.78373E+22 n/cm^2 (Table VI header, "TOTAL FLUENCE").
SPECTER_FLUENCE = 4.78373e22
SPECTER_BOUNDS = [0.0, 20.0]  # one-group fluence fold; bounds documented only

# Report print precision: 5 significant digits in every transcribed value.
PRINT_TOL = 1e-4

ED_EV = 40.0  # synthetic Fe-like threshold displacement energy (gates)
KAPPA = 0.8


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _lindhard_partition_py(t_ev: float, z: float, a: float) -> float:
    """Independent transcription: Robinson fit, self-recoil closed forms."""
    e_l = 30.724 * z * z * (2.0 * z ** (2.0 / 3.0)) ** 0.5 * 2.0
    k_l = 0.0793 * 2.0**0.75 * z ** (2.0 / 3.0) / a**0.5
    eps = t_ev / e_l
    g = eps + 3.4008 * eps ** (1.0 / 6.0) + 0.40244 * eps**0.75
    return 1.0 / (1.0 + k_l * g)


def analytic_gates() -> tuple[list[list[str]], list[str]]:
    """G1-G5 on synthetic inputs; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    # G1: NRT piecewise + Lindhard partition vs the independent transcription.
    worst = 0.0
    for t in [1e-3, 1.0, 1e2, 1e4, 1e5, 1e6, 1e8, 1e10, 1e12]:
        got = dmg.lindhard_partition(t, "Fe56", "Fe56")
        want = _lindhard_partition_py(t, 26.0, 56.0)
        worst = max(worst, rel_diff(got, want))
    rows.append(
        [
            "G1 Lindhard partition transcription",
            fmt(worst),
            "< 1e-12",
            _check(worst < 1e-12, "G1 partition"),
        ]
    )
    t_dam = dmg.damage_energy(1e5, "Fe56", "Fe56")
    nrt = dmg.nrt_displacements(1e5, ED_EV, "Fe56")
    err = rel_diff(nrt, KAPPA * t_dam / (2.0 * ED_EV))
    rows.append(["G1 NRT high branch", fmt(err), "< 1e-15", _check(err < 1e-15, "G1 NRT branch")])
    ok = (
        dmg.nrt_displacements(ED_EV * 0.5, ED_EV, "Fe56") == 0.0
        and dmg.nrt_displacements(ED_EV, ED_EV, "Fe56") == 1.0
    )
    rows.append(["G1 NRT piecewise low/plateau", "-", "exact", _check(ok, "G1 NRT piecewise")])
    notes.append(
        "G1: Fe-56 synthetic target, E_d = 40 eV; partition and branch "
        "recomputed from published constants."
    )

    # G2: arc-dpa construction.
    b_arc, c_arc = -0.55, 0.3
    junction = 2.0 * ED_EV / KAPPA
    xi_junction = dmg.arc_efficiency(junction, ED_EV, b_arc, c_arc)
    rows.append(
        [
            "G2 xi(2E_d/0.8) = 1",
            fmt(xi_junction),
            "exact",
            _check(xi_junction == 1.0, "G2 junction"),
        ]
    )
    err = rel_diff(dmg.arc_efficiency(1e9, ED_EV, b_arc, c_arc), c_arc)
    rows.append(["G2 xi saturation", fmt(err), "< 1e-3", _check(err < 1e-3, "G2 saturation")])
    xs = [1.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    xi = [dmg.arc_efficiency(t, ED_EV, b_arc, c_arc) for t in xs]
    ok = all(xi[i + 1] < xi[i] for i in range(len(xi) - 1))
    rows.append(["G2 xi monotone (b<0)", "-", "decreasing", _check(ok, "G2 monotone")])
    t_dam = dmg.damage_energy(1e5, "Fe56", "Fe56")
    want = KAPPA * t_dam * dmg.arc_efficiency(t_dam, ED_EV, b_arc, c_arc) / (2.0 * ED_EV)
    err = rel_diff(dmg.arc_displacements(1e5, ED_EV, "Fe56", b_arc, c_arc), want)
    rows.append(
        ["G2 arc high branch = NRT x xi", fmt(err), "< 1e-15", _check(err < 1e-15, "G2 arc branch")]
    )
    notes.append(
        "G2: synthetic arc constants (b_arc=-0.55, c_arc=0.3); the fitted "
        "material table is caller data by design."
    )

    # G3: folds against hand sums.
    bounds3 = [0.0, 0.1, 1.0, 20.0]
    flux3 = [1.0e12, 2.0e12, 4.0e12]
    resp3 = [100.0, 200.0, 50.0]
    hand = sum(f * r for f, r in zip(flux3, resp3, strict=True))
    err = rel_diff(dmg.nrt_dpa(flux3, resp3, bounds3, 2.0), 2.0e-24 * hand)
    rows.append(["G3 nrt_dpa hand sum", fmt(err), "< 1e-15", _check(err < 1e-15, "G3 nrt_dpa")])
    err = rel_diff(dmg.gas_appm(flux3, resp3, bounds3, 2.0), 2.0e-18 * hand)
    rows.append(["G3 gas_appm hand sum", fmt(err), "< 1e-15", _check(err < 1e-15, "G3 gas_appm")])
    wide = [0.0, 5.0, 10.0, 20.0]
    ok = dmg.nrt_dpa(flux3, resp3, wide, 2.0) == dmg.nrt_dpa(flux3, resp3, bounds3, 2.0)
    rows.append(
        ["G3 piecewise-constant (bounds width-free)", "-", "identical", _check(ok, "G3 width-free")]
    )

    # G4: invariants and loud errors. "Exact" means exact equality with the
    # IEEE hand product of the live groups — 1e-24 * 4e14 rounds one ulp
    # below the 4e-10 literal, so a decimal-literal target is unachievable.
    live = [0.0, 2.0e12, 0.0]
    got = dmg.nrt_dpa(live, resp3, bounds3, 1.0)
    want = 1.0e-24 * (2.0e12 * resp3[1])
    rows.append(
        ["G4 zero-flux groups", fmt(got), "exact hand product", _check(got == want, "G4 zero flux")]
    )
    try:
        dmg.he_dpa_ratio([0.0, 0.0, 0.0], resp3, resp3, bounds3, 1.0)
        ok = False
    except ValueError as exc:
        ok = "zero dpa" in str(exc)
    rows.append(["G4 He/dpa at zero dpa", "-", "loud named error", _check(ok, "G4 ZeroDpa")])
    for label, call in [
        ("G4 negative flux", lambda: dmg.nrt_dpa([-1.0, 1.0, 1.0], resp3, bounds3, 1.0)),
        ("G4 non-monotonic bounds", lambda: dmg.nrt_dpa(flux3, resp3, [0.0, 0.5, 0.5, 1.0], 1.0)),
    ]:
        try:
            call()
            ok = False
        except ValueError:
            ok = True
        rows.append([label, "-", "ValueError", _check(ok, label)])

    # G5: UQ hook — pinned seed, k-SE gates, determinism.
    flux2 = [1.0e13, 3.0e12]
    resp2 = [80.0, 160.0]
    bounds2 = [0.0, 0.5, 20.0]
    mean4 = [0.0, 0.0, 0.0, 0.0]
    cov4 = [
        [0.0004, 0.0, 0.0, 0.0],
        [0.0, 0.0009, 0.0, 0.0],
        [0.0, 0.0, 0.0025, 0.0],
        [0.0, 0.0, 0.0, 0.0001],
    ]
    out = dmg.fold_uq("nrt_dpa", flux2, resp2, bounds2, 2.0, mean4, cov4, 20_000, 20260915, 5.0)
    mean_se = out["analytic_std"] / math.sqrt(20_000)
    std_se = out["analytic_std"] / math.sqrt(2.0 * (20_000 - 1))
    ok = (
        out["passed"]
        and abs(out["mean"] - out["expected"]) <= 5.0 * mean_se
        and abs(out["std"] - out["analytic_std"]) <= 5.0 * std_se
        and out["expected"] == out["nominal"]
    )
    rows.append(["G5 UQ k-SE gate (n=20000)", "-", "passed", _check(ok, "G5 k-SE gate")])
    again = dmg.fold_uq("nrt_dpa", flux2, resp2, bounds2, 2.0, mean4, cov4, 20_000, 20260915, 5.0)
    rows.append(
        ["G5 pinned-seed determinism", "-", "identical", _check(again == out, "G5 determinism")]
    )
    try:
        dmg.fold_uq("he_dpa_ratio", flux2, resp2, bounds2, 1.0, mean4, cov4, 16, 1, 5.0)
        ok = False
    except ValueError as exc:
        ok = "not yet supported" in str(exc)
    rows.append(["G5 ratio UQ named-open", "-", "NotYetSupported", _check(ok, "G5 ratio open")])
    notes.append(
        "G5: seeded MVN over the caller [flux, response] block (linalg "
        "engine); exact bilinear expectation, first-order propagated std."
    )
    return rows, notes


# --------------------------------------------------------------------------
# SPECTER oracle section (container): transcribed spots from ANL/FPP/TM-197.
# --------------------------------------------------------------------------

# (label, response value in barns as printed, printed metric, provenance needles)
# Needles are sign-free digit/E strings (see _normalize).
DPA_SPOTS = [
    ("Fe dpa", 1.9118e2, 9.1455e0, ["19118E02", "91455E00"]),
    ("Ti dpa", 2.1873e2, 1.0464e1, ["21873E02", "10464E01"]),
    ("Cu dpa", 1.8649e2, 8.9212e0, ["18649E02", "89212E00"]),
]

# (label, sigma in mb as printed, printed appm, provenance needles)
GAS_SPOTS = [
    ("C12 He appm", 4.5151e-1, 2.1599e1, ["45151E01", "21599E01"]),
    ("Li7 He appm", 2.8447e1, 1.3608e3, ["28447E01", "13608E03"]),
    ("B10 H appm", 1.7095e0, 8.1780e1, ["17095E00", "81780E01"]),
    ("N14 H appm", 6.1271e2, 2.9310e4, ["61271E02", "29310E04"]),
]

# Fe He/dpa ratio from Table VI: helium appm 3.1944E+00 and DPA 9.1455E+00,
# with the spectrum-averaged He production cross section 6.6776E-05 barns
# printed in the same Fe row (the two prints satisfy appm = Phi*sigma*1e-18
# exactly, so the sigma transcription is corroborated in-report).
FE_HE_SPOT = ("Fe He/dpa", 6.6776e-5, 1.9118e2, 3.1944e0 / 9.1455e0, ["66776E05", "31944E00"])

_OCR_TABLE = str.maketrans(
    {"O": "0", "o": "0", "S": "5", "s": "5", "l": "1", "I": "1", "B": "8", "G": "6"}
)


def _normalize(text: str) -> str:
    """OCR-tolerant normalization: fix scan errors first (0/O, 5/S, 1/l,
    8/B, 6/G), then keep only digits and the exponent marker so line breaks,
    table spacing, and +/- sign confusion cannot split or corrupt a value.
    The check is a presence anchor only — the numeric gates do the judging.
    """
    fixed = text.translate(_OCR_TABLE)
    return "".join(c for c in fixed if c.isdigit() or c == "E")


def _fetch_specter() -> Path:
    """Return the cached SPECTER PDF, downloading once if absent."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if not SPECTER_PDF.exists():
        with urllib.request.urlopen(SPECTER_URL, timeout=120) as resp:
            SPECTER_PDF.write_bytes(resp.read())
    digest = hashlib.sha256(SPECTER_PDF.read_bytes()).hexdigest()
    if digest != SPECTER_SHA256:
        raise RuntimeError(
            f"SPECTER PDF SHA-256 mismatch: {digest} (expected {SPECTER_SHA256}); "
            "delete the cache file to re-fetch"
        )
    return SPECTER_PDF


def oracle_check_specter() -> tuple[list[list[str]], list[str], bool]:
    """SPECTER spot folds vs the report's printed results. Returns (rows, notes, skipped)."""
    try:
        pdf = _fetch_specter()
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Oracle check (SPECTER) SKIPPED: {exc}"
        print(note)
        return [], [note], True
    pdftotext = shutil.which("pdftotext")
    if pdftotext is None:
        note = "Oracle check (SPECTER) SKIPPED: pdftotext (poppler-utils) not installed"
        print(note)
        return [], [note], True
    text = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"], capture_output=True, text=True, check=True
    ).stdout
    blob = _normalize(text)
    for needle in ["478373E22"]:  # the HFIR-CTR32 fluence header
        if needle not in blob:
            note = (
                f"Oracle check (SPECTER) SKIPPED: report text does not contain the "
                f"HFIR-CTR32 fluence ({needle}); unexpected OCR layer"
            )
            print(note)
            return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []
    fluence = [SPECTER_FLUENCE]

    def _transcription_ok(needles: list[str]) -> bool:
        return all(n in blob for n in needles)

    for label, xs_barns, printed, needles in DPA_SPOTS:
        got = dmg.nrt_dpa(fluence, [xs_barns], SPECTER_BOUNDS, 1.0)
        err = rel_diff(got, printed)
        ok_text = _transcription_ok(needles)
        ok = err < PRINT_TOL and ok_text
        if not ok and not ok_text:
            print(
                f"FAIL: {label} transcription not found in report text: {needles}", file=sys.stderr
            )
        rows.append(
            [f"S1 {label}", fmt(err), f"< {PRINT_TOL:g} (print)", _check(ok, f"S1 {label}")]
        )

    for label, sigma_mb, printed_appm, needles in GAS_SPOTS:
        got = dmg.gas_appm(fluence, [sigma_mb * 1e-3], SPECTER_BOUNDS, 1.0)
        err = rel_diff(got, printed_appm)
        ok_text = _transcription_ok(needles)
        ok = err < PRINT_TOL and ok_text
        if not ok and not ok_text:
            print(
                f"FAIL: {label} transcription not found in report text: {needles}", file=sys.stderr
            )
        rows.append(
            [f"S2 {label}", fmt(err), f"< {PRINT_TOL:g} (print)", _check(ok, f"S2 {label}")]
        )

    label, he_xs_b, dpa_xs_b, printed_ratio, needles = FE_HE_SPOT
    got = dmg.he_dpa_ratio(fluence, [he_xs_b], [dpa_xs_b], SPECTER_BOUNDS, 1.0)
    err = rel_diff(got, printed_ratio)
    ok_text = _transcription_ok(needles)
    ok = err < PRINT_TOL and ok_text
    if not ok and not ok_text:
        print(f"FAIL: {label} transcription not found in report text: {needles}", file=sys.stderr)
    rows.append([f"S3 {label}", fmt(err), f"< {PRINT_TOL:g} (print)", _check(ok, f"S3 {label}")])

    notes.append(
        "Spots transcribed from ANL/FPP/TM-197 (US-gov PD): Table VI (HFIR-CTR32, "
        "fluence 4.78373E+22) spectrum-averaged dpa cross sections and DPA for Fe/Ti/Cu; "
        "Table V spectral-averaged gas production cross sections (mb) and GAS(APPM). "
        "One-group fluence folds (seconds=1); bounds are documentation-only. Tolerance "
        "1e-4 relative = 2 half-ulp of the 5-digit print, the finest agreement "
        "transcribed inputs allow. SPECTER stays this script's oracle; the "
        "crate's opt-in Table VII fallback is a separate caller input."
    )
    notes.append(
        "Fe He/dpa uses the He production cross section printed in the Fe row "
        "(6.6776E-05 barns), which the report's own helium appm (3.1944E+00) "
        "corroborates via appm = fluence x sigma x 1e-18."
    )
    return rows, notes, False


def main() -> int:
    report = Report("damage", "Damage/gas metrics (`damage_vs_specter.py`)")
    report.prose(
        "Two-part oracle for `nucleide.damage`: analytic gates G1-G5 on synthetic "
        "inputs (always run — NRT/arc closed forms, fold conventions, invariants, "
        "and the pinned-seed MVN UQ gate), and a SPECTER cross-check "
        "(ANL/FPP/TM-197, US-gov PD) folding report-transcribed spots at print "
        "precision (container only; loud SKIP outside)."
    )
    rows1, notes1 = analytic_gates()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows1)
    rows2, notes2, skipped = oracle_check_specter()
    for note in notes2:
        report.prose(note)
    if skipped:
        report.table(
            ["Gate", "Status"],
            [
                [gate, "SKIP (SPECTER unavailable)"]
                for gate in ("S1 dpa spots", "S2 gas spots", "S3 He/dpa")
            ],
        )
    else:
        report.table(["Gate", "Rel err", "Tol", "Status"], rows2)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} damage check(s) failed", file=sys.stderr)
        return 1
    print("All damage checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
