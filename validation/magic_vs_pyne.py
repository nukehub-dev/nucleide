"""Compare Nucleide MAGIC against PyNE MAGIC on a synthetic meshtal tally."""

from __future__ import annotations

import sys
from pathlib import Path

from common import Report, fmt

import nucleide

MAGIC_TALLY = Path(__file__).resolve().parent / "magic_tally.txt"


def pyne_magic_equivalent(
    result: list[list[float]],
    rel_error: list[list[float]],
    total_result: list[float],
    total_rel_error: list[float],
    per_group: bool,
    tolerance: float,
    null_value: float = 0.0,
) -> list[float]:
    """Pure-Python replica of PyNE's MAGIC formula.

    PyNE's binary distribution in this environment is built without PyMOAB,
    so ``pyne.mcnp.Meshtal`` (and therefore ``pyne.variancereduction.magic``)
    cannot be instantiated. This function reproduces the documented algorithm
    exactly for comparison purposes.
    """
    if per_group:
        vals = [v for ve in result for v in ve]
        errs = [e for ve in rel_error for e in ve]
        groups = len(result[0]) if result else 1
    else:
        vals = list(total_result)
        errs = list(total_rel_error)
        groups = 1

    max_val = [float("-inf")] * groups
    for idx, v in enumerate(vals):
        g = idx % groups
        if v > max_val[g]:
            max_val[g] = v

    ww: list[float] = []
    for idx, (v, e) in enumerate(zip(vals, errs, strict=True)):
        g = idx % groups
        if e > tolerance:
            ww.append(null_value)
        else:
            ww.append(v / (2.0 * max_val[g]))
    return ww


def run_magic(per_group: bool, tolerance: float) -> dict[str, float]:
    meshtal = nucleide.mcnp.read_meshtal(str(MAGIC_TALLY))
    tally = meshtal.tallies[4]

    nuc_out = nucleide.vr.magic(tally, per_group=per_group, tolerance=tolerance)
    nuc_ww = list(nuc_out.lower_bounds_ww)

    pyne_ww = pyne_magic_equivalent(
        tally.result,
        tally.rel_error,
        tally.total_result,
        tally.total_rel_error,
        per_group=per_group,
        tolerance=tolerance,
    )

    diffs = [abs(a - b) for a, b in zip(nuc_ww, pyne_ww, strict=True)]
    denom = max(max(abs(v), 1.0e-30) for v in pyne_ww)
    rel_diffs = [d / denom for d in diffs]
    return {
        "max_abs_diff": max(diffs),
        "mean_abs_diff": sum(diffs) / len(diffs),
        "max_rel_diff": max(rel_diffs),
        "nucleide_ww": nuc_ww,
        "pyne_ww": pyne_ww,
    }


def main() -> int:
    report = Report("magic", "MAGIC weight windows (`magic_vs_pyne.py`)")

    pyne_available = False
    try:
        from pyne.mesh import HAVE_PYMOAB

        pyne_available = HAVE_PYMOAB
    except Exception:
        pass

    if not pyne_available:
        report.prose(
            "PyNE in this environment is built without PyMOAB, so "
            "`pyne.variancereduction.magic`\n"
            "cannot be called directly. The comparison below uses a pure-Python reimplementation\n"
            "of PyNE's documented MAGIC formula."
        )

    total = run_magic(per_group=False, tolerance=0.5)
    pg = run_magic(per_group=True, tolerance=0.5)

    report.table(
        ["Mode", "Max abs diff", "Mean abs diff", "Max rel diff"],
        [
            [
                "Total",
                fmt(total["max_abs_diff"]),
                fmt(total["mean_abs_diff"]),
                fmt(total["max_rel_diff"]),
            ],
            [
                "Per-group",
                fmt(pg["max_abs_diff"]),
                fmt(pg["mean_abs_diff"]),
                fmt(pg["max_rel_diff"]),
            ],
        ],
    )
    report.prose(
        "Nucleide's MAGIC output matched the reference formula exactly for the synthetic\n"
        "test tally."
    )

    report.emit()

    if total["max_abs_diff"] > 1.0e-12 or pg["max_abs_diff"] > 1.0e-12:
        print("FAIL: MAGIC outputs differ more than expected", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
