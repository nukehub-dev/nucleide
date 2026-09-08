"""Self-consistency checks for single-material code emission (`nucleide-emit`).

No external code emits one material to MCNP/Serpent/FLUKA/ALARA/PARTISN cards
with a mass-drift report, so there is no independent oracle: uranium metal
must round-trip losslessly on all five dialects, the known FLUKA O16 gap must
appear as reported drift (not silent loss), and the MCNP/ALARA texts must
carry the re-parse-verified flag.
"""

from __future__ import annotations

import sys

from common import Report, fmt

import nucleide.emit as emit

FAILURES = 0


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def main() -> int:
    report = Report("emit", "Emission drift (`emit_vs_self.py`)")
    report.prose(
        "Self-consistency oracle for `nucleide.emit` (no external code offers"
        " this comparison). Uranium metal must emit losslessly on all five"
        " dialects; the FLUKA O16 gap is asserted as reported drift."
    )

    rows: list[list[str]] = []
    table = emit.emit_drift_table({"U235": 5.0, "U238": 95.0}, "umetal", density=19.1)
    for row in table:
        rows.append(
            [
                row["code"],
                fmt(row["mass_out"]),
                fmt(row["rel_drift"]),
                str(len(row["dropped"])),
                str(row["reparsed"]),
                _check(
                    abs(row["mass_out"] - 100.0) < 1e-9
                    and abs(row["rel_drift"]) < 1e-9
                    and row["dropped"] == [],
                    f"{row['code']} lossless uranium metal",
                ),
            ]
        )
    report.table(["Code", "Mass out (g)", "Rel drift", "Dropped", "Reparsed", "Status"], rows)
    reparsed = {r["code"] for r in table if r["reparsed"]}
    report.prose(
        "Re-parse verified dialects: " + ", ".join(sorted(reparsed)) + "."
        if reparsed
        else "No dialect re-parse verified."
    )
    if reparsed != {"MCNP", "ALARA"}:
        _check(False, f"reparsed dialects == {{MCNP, ALARA}}, got {sorted(reparsed)}")

    water = emit.emit_drift_table({"U235": 1.0, "O16": 2.0}, "uox", density=10.0)
    fluka = next(r for r in water if r["code"] == "FLUKA")
    report.prose(
        f"Water-like mix: FLUKA accounts {fmt(fluka['mass_out'])} of 3.0 g with"
        f" {len(fluka['dropped'])} dropped nuclide(s)."
    )
    _check(
        abs(fluka["mass_out"] - 1.0) < 1e-12
        and [d["nuclide"] for d in fluka["dropped"]] == ["O16"],
        "FLUKA O16 gap reported as drift",
    )

    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} emission check(s) failed", file=sys.stderr)
        return 1
    print("All emission checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
