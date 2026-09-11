"""Spectroscopy cross-check (`nucleide-spectroscopy` vs hand values + PyNE).

Two tiers:

1. Synthetic gates (always run): E1 rectangular / E2 five-point smoothing,
   E3 background / E4 gross / E5 net counts, E6 energy bins, E7 efficiency
   (including the 6-coefficient golden), E8 X-ray algebra, and the
   dollar/plain `.spe` fixture parse with cross-format counts equality.
   Inputs are the hand-built synthetic `fixtures/spectroscopy/` files.
2. PyNE cross-check: the same vectors and fixture files driven through the
   upstream ``pyne.spectanalysis`` / ``pyne.gammaspec`` readers on identical
   runtime-read inputs. PyNE is an optional oracle dependency: if it cannot
   be imported, tier 2 is reported as SKIP with its reason (never silently).
   The X-ray algebra has no container check — the upstream routine needs its
   HDF5 atomic table and no atomic values are vendored here — recorded below.
"""

from __future__ import annotations

import sys
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.spectroscopy as sp

FAILURES = 0
REPO_ROOT = Path(__file__).resolve().parent.parent
FIX = REPO_ROOT / "fixtures" / "spectroscopy"

COUNTS = [2.0, 5.0, 1.0, 6.0, 3.0, 8.0, 4.0]
CHANNELS = [float(c) for c in range(len(COUNTS))]
EFF_COEFF = [
    -2.81861504261204,
    -0.727352820018942,
    -0.0395798886481904,
    -0.0592305254664096,
    0.023772637347443,
    0.0325306475072671,
]
ATOMIC = {
    "k_shell_fluor": 0.9,
    "l_shell_fluor": 0.4,
    "prob": 0.8,
    "kb_to_ka": 0.2,
    "ka2_to_ka1": 0.5,
    "ka1_en_kev": 10.0,
    "ka2_en_kev": 20.0,
    "kb_en_kev": 30.0,
    "l_en_kev": 40.0,
}


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _worst_rel(got: list[float], want: list[float]) -> float:
    return max(rel_diff(g, w) for g, w in zip(got, want, strict=True))


def tier1() -> tuple[list[list[str]], list[str], list[list[str]], str]:
    """Synthetic gates E1-E8.

    Returns (gate rows, prose notes, overlay rows, background level). The
    overlay rows tabulate the per-channel smoothing overlay (``Channel``,
    ``Raw counts``, ``Rect-smoothed (m=5)``, ``Five-point smoothed``) for
    ``make_figures.py``; no gates live in those tables.
    """
    rows: list[list[str]] = []
    notes: list[str] = []

    rect = sp.rect_smooth(COUNTS, 5)
    err = _worst_rel(rect, [2.0, 5.0, 3.4, 4.6, 4.4, 8.0, 4.0])
    notes.append(f"E1 rect m=5: worst rel err {err:.3e} (hand values).")
    rows.append(["E1 rect smooth", fmt(err), "< 1e-12", _check(err < 1e-12, "E1")])

    five = sp.five_point_smooth(COUNTS)
    err = _worst_rel(five, [2.0, 5.0, 30.0 / 9.0, 39.0 / 9.0, 42.0 / 9.0, 8.0, 4.0])
    notes.append(f"E2 five-point: worst rel err {err:.3e} (hand values).")
    rows.append(["E2 five-point smooth", fmt(err), "< 1e-12", _check(err < 1e-12, "E2")])

    overlay_rows = [
        [str(c), repr(r), repr(s), repr(f)]
        for c, r, s, f in zip(range(len(COUNTS)), COUNTS, rect, five, strict=True)
    ]
    bg_level = repr(sp.calc_bg(COUNTS, CHANNELS, 2, 5, 1))

    err = rel_diff(sp.calc_bg(COUNTS, CHANNELS, 2, 5, 1), 76.0 / 6.0)
    rows.append(["E3 background", fmt(err), "< 1e-12", _check(err < 1e-12, "E3")])
    err = rel_diff(sp.gross_count(COUNTS, CHANNELS, 2, 5), 10.0)
    rows.append(["E4 gross count", fmt(err), "< 1e-12", _check(err < 1e-12, "E4")])
    err = rel_diff(sp.net_counts(COUNTS, CHANNELS, 2, 5, 1), -16.0 / 6.0)
    rows.append(["E5 net counts", fmt(err), "< 1e-12", _check(err < 1e-12, "E5")])

    got = sp.energy_bins([0.0, 1.0, 2.0], [1.5, 2.0, 0.5])
    err = _worst_rel(got, [1.5, 4.0, 7.5])
    rows.append(["E6 energy bins", fmt(err), "< 1e-12", _check(err < 1e-12, "E6")])

    err = rel_diff(sp.detector_efficiency(1.0, EFF_COEFF, 1), 0.059688551591347033)
    notes.append("E7 golden: 6-coefficient fit-1 efficiency at 1 MeV.")
    rows.append(["E7 efficiency golden", fmt(err), "< 1e-12", _check(err < 1e-12, "E7")])

    lines = sp.xray_lines(ATOMIC, k_conv=2.0, l_conv=3.0)
    err = _worst_rel([i for _, i in lines], [1.0, 0.5, 0.3, 1.84])
    notes.append("E8 X-ray: k=2/l=3 hand intensities (Ka1, Ka2, Kb, L).")
    rows.append(["E8 xray lines", fmt(err), "< 1e-12", _check(err < 1e-12, "E8")])

    dollar = sp.read_dollar_spe(str(FIX / "dollar_min.spe"))
    plain = sp.read_spe(str(FIX / "plain_min.spe"))
    ok = dollar["counts"] == plain["counts"] and len(dollar["counts"]) == 8
    notes.append(
        f"SPE fixtures: cross-format counts equality over {len(dollar['counts'])} channels."
    )
    rows.append(["SPE cross-format counts", "n=8", "equal", _check(ok, "SPE counts")])
    return rows, notes, overlay_rows, bg_level


def tier2_pyne() -> tuple[list[list[str]], list[str], bool]:
    """PyNE cross-check. Returns (rows, notes, skipped)."""
    try:
        from pyne import gammaspec
        from pyne import spectanalysis as sa
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Tier 2 (PyNE cross-check) SKIPPED: {exc}"
        print(note)
        return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []

    upstream = sa.PhSpectrum()
    upstream.channels = list(range(len(COUNTS)))
    upstream.counts = list(COUNTS)

    def smooth_case(label: str, ours: list[float], theirs: list[float]) -> None:
        err = _worst_rel(list(ours), list(theirs))
        rows.append([label, fmt(err), "< 1e-12", _check(err < 1e-12, label)])

    smooth_case(
        "rect m=5 vs PyNE",
        sp.rect_smooth(COUNTS, 5),
        list(sa.rect_smooth(upstream, 5).counts),
    )
    smooth_case(
        "five-point vs PyNE",
        sp.five_point_smooth(COUNTS),
        list(sa.five_point_smooth(upstream).counts),
    )

    for label, ours, theirs in [
        ("bg vs PyNE", sp.calc_bg(COUNTS, CHANNELS, 2, 5, 1), sa.calc_bg(upstream, 2, 5, 1)),
        ("gross vs PyNE", sp.gross_count(COUNTS, CHANNELS, 2, 5), sa.gross_count(upstream, 2, 5)),
        (
            "net vs PyNE",
            sp.net_counts(COUNTS, CHANNELS, 2, 5, 1),
            sa.net_counts(upstream, 2, 5, 1),
        ),
    ]:
        err = rel_diff(ours, float(theirs))
        rows.append([label, fmt(err), "< 1e-12", _check(err < 1e-12, label)])

    g = gammaspec.GammaSpectrum(calib_e_fit=[1.5, 2.0, 0.5])
    g.channels = [0.0, 1.0, 2.0]
    g.calc_ebins()
    err = _worst_rel(sp.energy_bins([0.0, 1.0, 2.0], [1.5, 2.0, 0.5]), list(g.ebin))
    rows.append(["ebins vs PyNE", fmt(err), "< 1e-12", _check(err < 1e-12, "ebins vs PyNE")])

    err = rel_diff(
        sp.detector_efficiency(1.0, EFF_COEFF, 1),
        float(gammaspec.calc_e_eff(1, EFF_COEFF, 1)),
    )
    rows.append(["efficiency vs PyNE", fmt(err), "< 1e-12", _check(err < 1e-12, "eff vs PyNE")])

    for label, path, reader, ours_fn in [
        ("dollar", "dollar_min.spe", gammaspec.read_dollar_spe_file, sp.read_dollar_spe),
        ("plain", "plain_min.spe", gammaspec.read_spe_file, sp.read_spe),
    ]:
        ref = reader(str(FIX / path))
        got_spe = ours_fn(str(FIX / path))
        err = _worst_rel(list(got_spe["counts"]), [float(c) for c in ref.counts])
        rows.append([f"{label} counts vs PyNE", fmt(err), "< 1e-12", _check(err < 1e-12, label)])
        err = _worst_rel(list(got_spe["ebin"]), [float(e) for e in ref.ebin])
        rows.append([f"{label} ebin vs PyNE", fmt(err), "< 1e-12", _check(err < 1e-12, label)])

    notes.append(
        "X-ray algebra has no container check: the upstream routine reads "
        "its HDF5 atomic table at run time and no atomic values are "
        "vendored here, so E8 is pinned by the synthetic hand values above."
    )
    return rows, notes, False


def main() -> int:
    report = Report("spectroscopy", "Spectroscopy (`spectroscopy_vs_pyne.py`)")
    report.prose(
        "Two-tier oracle for `nucleide.spectroscopy`: synthetic E1-E8 gates "
        "on hand-built fixtures (tier 1, always run), and a cross-check "
        "against the upstream `pyne.spectanalysis` / `pyne.gammaspec` "
        "routines on identical runtime inputs (tier 2)."
    )
    rows1, notes1, overlay_rows, bg_level = tier1()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Rel err", "Tol", "Status"], rows1)
    report.heading("Smoothing overlay (figure source)", level=3)
    report.table(
        ["Channel", "Raw counts", "Rect-smoothed (m=5)", "Five-point smoothed"],
        overlay_rows,
    )
    report.table(["Quantity", "Value"], [["Background level (E3, channels 2..5)", bg_level]])
    rows2, notes2, skipped = tier2_pyne()
    for note in notes2:
        report.prose(note)
    if skipped:
        report.table(["Gate", "Status"], [["PyNE cross-check", "SKIP (pyne unavailable)"]])
    else:
        report.table(["Gate", "Rel err", "Tol", "Status"], rows2)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} spectroscopy check(s) failed", file=sys.stderr)
        return 1
    print("All spectroscopy checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
