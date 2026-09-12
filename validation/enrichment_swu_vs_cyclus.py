"""SWU cross-check (`nucleide-enrichment` vs closed form + cyclus toolkit).

Two tiers:

1. Closed-form gates (always run): the Dirac separation potential
   ``V(x) = (2x - 1) ln(x / (1 - x))`` and the SWU balance
   ``SWU = P*V(xP) + T*V(xT) - F*V(xF)`` on the cyclus assay ladder
   (feed 0.0072, product 0.05, tails 0.002, 10 units of product), with
   hand-recomputed expected values quoted below.
2. cyclus toolkit cross-check: the same quantities via
   ``cyclus.toolkit.enrichment`` (``ValueFunc``/``SwuRequired``) when the
   cyclus Python package is importable. cyclus is an optional oracle
   dependency: if it cannot be imported, tier 2 is reported as SKIP with
   its reason (never silently).
"""

from __future__ import annotations

import math
import sys

from common import Report, fmt, rel_diff

import nucleide.enrichment as enr

FAILURES = 0

# Cyclus toolkit assay ladder (feed, product, tails) and product quantity.
XF, XP, XT, P = 0.0072, 0.05, 0.002, 10.0
# Hand-recomputed closed-form values: V(x) = (2x-1)*ln(x/(1-x));
# F = P*(xP-xT)/(xF-xT) = 92.3076923076923; T = F-P = 82.3076923076923;
# SWU = P*V(xP)+T*V(xT)-F*V(xF) = 87.59916188589864.
V_FEED = 4.855507353675083
V_PROD = 2.6499950812497963
V_TAIL = 6.187755671368513
FEED_QTY = 92.3076923076923
TAILS_QTY = 82.3076923076923
SWU_TOTAL = 87.59916188589864


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _value(x: float) -> float:
    """Independent closed-form separation potential (no nucleide import)."""
    return (2.0 * x - 1.0) * math.log(x / (1.0 - x))


def tier1() -> tuple[list[list[str]], list[str]]:
    """Closed-form gates G1-G3. Returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    for label, x, want in (("feed", XF, V_FEED), ("product", XP, V_PROD), ("tails", XT, V_TAIL)):
        got = enr.value_func(x)
        err = rel_diff(got, want)
        notes.append(f"G1 V({label}={x}): nucleide {got:.6e}, closed form {want:.6e}.")
        rows.append([f"G1 value_func ({label})", fmt(err), "< 1e-12", _check(err < 1e-12, label)])

    swu_feed = FEED_QTY * enr.swu_per_feed(XF, XP, XT)
    swu_prod = P * enr.swu_per_prod(XF, XP, XT)
    swu_tail = TAILS_QTY * enr.swu_per_tail(XF, XP, XT)
    for label, got in (("feed", swu_feed), ("product", swu_prod), ("tails", swu_tail)):
        err = rel_diff(got, SWU_TOTAL)
        notes.append(f"G2 SWU via {label}: nucleide {got:.6e}, closed form {SWU_TOTAL:.6e}.")
        rows.append(
            [f"G2 SWU total (via {label})", fmt(err), "< 1e-12", _check(err < 1e-12, label)]
        )

    spread = max(swu_feed, swu_prod, swu_tail) - min(swu_feed, swu_prod, swu_tail)
    err = spread / max(abs(SWU_TOTAL), 1.0e-30)
    notes.append(f"G3 three-stream spread: {spread:.3e} SWU across feed/product/tails views.")
    rows.append(["G3 stream agreement", fmt(err), "< 1e-12", _check(err < 1e-12, "G3 spread")])
    return rows, notes


def tier2_cyclus() -> tuple[list[list[str]], list[str], bool]:
    """cyclus toolkit cross-check. Returns (rows, notes, skipped)."""
    try:
        import cyclus.toolkit.enrichment as cy_enr
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Tier 2 (cyclus toolkit cross-check) SKIPPED: {exc}"
        print(note)
        return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []
    for label, x in (("feed", XF), ("product", XP), ("tails", XT)):
        got = float(cy_enr.ValueFunc(x))
        want = enr.value_func(x)
        err = rel_diff(got, want)
        rows.append(
            [
                f"cyclus ValueFunc ({label})",
                fmt(want),
                fmt(got),
                fmt(err),
                _check(err < 1e-12, f"cyclus {label}"),
            ]
        )
    assays = cy_enr.Assays(XF, XP, XT)
    got = float(cy_enr.SwuRequired(P, assays))
    err = rel_diff(got, SWU_TOTAL)
    notes.append(f"cyclus SwuRequired({P}, assays): {got:.6e} vs closed form {SWU_TOTAL:.6e}.")
    rows.append(
        [
            "cyclus SwuRequired",
            fmt(SWU_TOTAL),
            fmt(got),
            fmt(err),
            _check(err < 1e-9, "cyclus SwuRequired"),
        ]
    )
    return rows, notes, False


def main() -> int:
    report = Report("enrichment_swu", "Enrichment SWU (`enrichment_swu_vs_cyclus.py`)")
    report.prose(
        "Two-tier oracle for the `nucleide-enrichment` separative-work helpers: "
        "closed-form gates G1-G3 on the cyclus assay ladder (tier 1, always run), "
        "and a live cross-check against `cyclus.toolkit.enrichment` "
        "(`ValueFunc`/`SwuRequired`) when the cyclus package is importable "
        "(tier 2)."
    )
    rows1, notes1 = tier1()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Rel err", "Tol", "Status"], rows1)
    rows2, notes2, skipped = tier2_cyclus()
    for note in notes2:
        report.prose(note)
    if skipped:
        report.table(
            ["Gate", "Status"],
            [["cyclus ValueFunc/SwuRequired", "SKIP (cyclus unavailable)"]],
        )
    else:
        report.table(["Gate", "Nucleide", "cyclus", "Rel err", "Status"], rows2)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} enrichment-SWU check(s) failed", file=sys.stderr)
        return 1
    print("All enrichment-SWU checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
