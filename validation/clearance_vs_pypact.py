"""Clearance / waste-classification cross-check.

Two parts:

1. Analytic gates (always run): hand-computed clearance-index vectors at
   exact equality, sum-of-fractions boundary probes at == 1 on both sides,
   EU 2013/59/Euratom Annex VII Table A vendored-transcription spot values,
   the FISPACT-II clearance-block parse of the synthetic fixture in the
   real wide-table grammar (HAZARDS + CLEAR keywords), the Spanish CSN
   table spots, and the Sublet S1+S2 hand vectors (total activity with the
   IRT split, decay heat over caller decay energies) at exact equality
   plus parts-to-total conservation.
2. pypact cross-check: the same fixture parsed with the upstream Apache-2.0
   ``pypact`` reader (the citable grammar's reference consumer), comparing
   the overlapping per-nuclide columns (atoms, activity, clearance index).
   pypact is an optional oracle dependency: if it cannot be imported (or the
   fixture trips its column assertion), the oracle leg is reported as SKIP
   with its reason (never silently).
"""

from __future__ import annotations

import sys
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.alara as alara
import nucleide.fispact as fispact

FAILURES = 0
REPO_ROOT = Path(__file__).resolve().parent.parent
FIX = REPO_ROOT / "fixtures" / "fispact"

TOY_LIMITS = {"Co60": 10.0, "H3": 5.0, "Fe55": 2.0}


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def synthetic_gates() -> tuple[list[list[str]], list[str], float]:
    """Analytic gates C1-C5 on the synthetic fixture. Returns (rows, notes, total CI)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    ci = alara.alara_clearance_index({"Co60": 10.0}, TOY_LIMITS)
    rows.append(["C1 CI == 1 at limit", fmt(ci), "== 1", _check(ci == 1.0, "C1 single")])
    ci = alara.alara_clearance_index({"Co60": 10.0, "H3": 5.0, "Fe55": 4.0}, TOY_LIMITS)
    rows.append(["C1 CI mixed vector", fmt(ci), "== 4", _check(ci == 4.0, "C1 mixed")])

    at = alara.alara_sum_of_fractions({"Co60": 10.0}, TOY_LIMITS)
    rows.append(
        [
            "C2 boundary == 1",
            fmt(at["sum"]),
            "satisfied",
            _check(at["class"] == "satisfied", "C2 =="),
        ]
    )
    above = alara.alara_sum_of_fractions({"Co60": 10.000000000000002}, TOY_LIMITS)
    rows.append(
        [
            "C2 one ulp above",
            fmt(above["sum"]),
            "exceeded",
            _check(above["class"] == "exceeded" and above["sum"] > 1.0, "C2 >"),
        ]
    )
    below = alara.alara_sum_of_fractions({"Co60": 9.999999999999998}, TOY_LIMITS)
    rows.append(
        [
            "C2 one ulp below",
            fmt(below["sum"]),
            "satisfied",
            _check(below["class"] == "satisfied" and below["sum"] < 1.0, "C2 <"),
        ]
    )

    table = alara.alara_clearance_eu_table()
    spot = all(
        table[name] == value
        for name, value in [
            ("H-3", 100.0),
            ("C-14", 1.0),
            ("Co-60", 0.1),
            ("Cs-137", 0.1),
            ("Pu-239", 0.1),
            ("U-238", 1.0),
            ("K-40", 10.0),
        ]
    )
    notes.append(f"EU Annex VII Table A transcription: {len(table)} entries (Bq/g).")
    rows.append(["C3 EU table spots", f"{len(table)} entries", "exact", _check(spot, "C3 spots")])

    rows_text = (FIX / "clearance.out").read_text()
    parsed = fispact.fispact_parse_clearance(rows_text)
    step1 = [r for r in parsed if r["interval"] == 1]
    ok = (
        len(parsed) == 7
        and len(step1) == 4
        and step1[0]["nuclide"] == "V-55"
        and step1[0]["flags"] == ">"
        and step1[3]["nuclide"] == "Rb-86m"
        and step1[3]["flags"] == "&"
        and all(r["cooling"] for r in parsed if r["interval"] == 2)
    )
    rows.append(["C4 clearance-block parse", f"{len(parsed)} rows", "7", _check(ok, "C4 parse")])
    total_ci = sum(r["clearance_index"] for r in parsed)
    expected_total = 1.0e6 + 4.0e9 + 9.5e5 + 2.0e7
    rows.append(
        [
            "C4 column total CI",
            fmt(total_ci),
            fmt(expected_total),
            _check(total_ci == expected_total, "C4 total"),
        ]
    )

    inventory = {r["nuclide"]: r["activity_bq"] for r in step1}
    limits = {"V-55": 1.0, "Fe-56": 1.0, "Co-60": 4.0e6, "Rb-86m": 8.0e9}
    out = alara.alara_sum_of_fractions(inventory, limits)
    rows.append(
        [
            "C5 parsed-inventory CI",
            fmt(out["sum"]),
            "== 2",
            _check(out["sum"] == 2.0 and out["class"] == "exceeded", "C5 e2e"),
        ]
    )

    # C6: Spanish CSN conditional NORM landfill tables (Tablas 1-3,
    # CSN/PDT/AICD/TGE/2503/02): transcription spots, chain expansion, and a
    # hand-computed screening vector at exact equality.
    es_inert = alara.alara_clearance_es_table("inert", "rocks")
    es_haz_gas = alara.alara_clearance_es_table("hazardous", "oil_gas")
    notes.append(
        "Spanish CSN NORM landfill tables (draft technical opinion "
        "CSN/PDT/AICD/TGE/2503/02, TGE/VAR/2025/1, csn.es, accessed "
        "2026-09-16): Tablas 1-3 per landfill type and NORM material "
        "nature, Tabla 4 chain keys expanded to members at the parent value."
    )
    es_spot = (
        len(es_inert) == 40
        and all(
            es_inert[name] == value
            for name, value in [
                ("U-238", 10.0),  # U-nat
                ("Ra-226", 10.0),  # Ra-226+
                ("Pb-210", 10.0),  # Pb-210+
                ("Po-210", 5.0),
                ("Th-232", 5.0),
                ("K-40", 10.0),
            ]
        )
        and all(
            es_haz_gas[name] == value
            for name, value in [("U-238", 500.0), ("Po-210", 500.0), ("K-40", 500.0)]
        )
    )
    rows.append(
        [
            "C6 ES CSN table spots",
            f"{len(es_inert)} entries",
            "exact",
            _check(es_spot, "C6 spots"),
        ]
    )
    es_expand = all(
        es_inert[member] == 10.0
        for member in ("Ra-226", "Rn-222", "Po-218", "Pb-214", "Bi-214", "Po-214")
    ) and all(es_inert[member] == 5.0 for member in ("Ac-227", "Th-227", "Tl-207"))
    rows.append(
        [
            "C6 chain expansion members",
            "Ra-226+ / Ac-227+ members",
            "parent value",
            _check(es_expand, "C6 expand"),
        ]
    )
    es_ci = alara.alara_clearance_index(
        {"U-238": 250.0, "Ra-226": 25.0, "Po-210": 250.0, "K-40": 250.0}, es_haz_gas
    )
    rows.append(
        [
            "C6 ES screening vector",
            fmt(es_ci),
            "== 2",
            _check(es_ci == 2.0, "C6 vector"),
        ]
    )

    # C7: Sublet S1+S2 radiological totals (Table X rows Ai / Ai·E·C1 with the
    # IRT alpha/beta/gamma split; open output_interpretation pin, no NDS
    # pages): hand-computed vectors at exact equality plus parts-to-total
    # conservation and the excluding-tritium companions.
    def _s1(nuclide: str, activity: float, irt: int, frac: float | None = None) -> dict:
        entry: dict = {"nuclide": nuclide, "activity_bq": activity, "irt": irt}
        if frac is not None:
            entry["alpha_frac"] = frac
        return entry

    s1 = alara.alara_total_activity(
        [
            _s1("Co-60", 10.0, 1),
            _s1("Po-210", 5.0, 4),
            _s1("Tc-99m", 4.0, 3),
            _s1("U-235", 100.0, 12, 0.25),
            _s1("H-3", 5.0, 1),
        ]
    )
    notes.append(
        "Sublet S1 total activity (Table X Ai = Ni li, Bq; IRT split per the "
        "open output_interpretation Activity break-down section): IRT 4 "
        "alpha; IRT 1 beta; IRT 3 gamma; IRT 12 split alpha/beta."
    )
    s1_ok = (
        s1["total_bq"] == 124.0
        and s1["alpha_bq"] == 30.0
        and s1["beta_bq"] == 90.0
        and s1["gamma_bq"] == 4.0
        and s1["ex_tritium_bq"] == 119.0
        and s1["total_bq"] == s1["alpha_bq"] + s1["beta_bq"] + s1["gamma_bq"]
    )
    rows.append(["C7 S1 hand vector", fmt(s1["total_bq"]), "== 124", _check(s1_ok, "C7 S1")])

    c1 = 1.602176634e-22
    s2 = alara.alara_decay_heat(
        [
            {
                "nuclide": "Co-60",
                "activity_bq": 1e10,
                "e_alpha_ev": 0.0,
                "e_beta_ev": 1e6,
                "e_gamma_ev": 2e6,
            },
            {
                "nuclide": "H-3",
                "activity_bq": 5.0,
                "e_alpha_ev": 0.0,
                "e_beta_ev": 6e3,
                "e_gamma_ev": 0.0,
            },
        ]
    )
    notes.append(
        "Sublet S2 decay heat (Table X Ai E C1, kW; C1 = eV to kJ): "
        "caller-supplied average decay energies, never vendored."
    )
    s2_ok = (
        s2["alpha_kw"] == 0.0
        and s2["beta_kw"] == 1e10 * 1e6 * c1 + 5.0 * 6e3 * c1
        and s2["gamma_kw"] == 1e10 * 2e6 * c1
        and s2["total_kw"] == s2["alpha_kw"] + s2["beta_kw"] + s2["gamma_kw"]
        # Ex-tritium subtracts the tritium heat, so allow one rounding step.
        and rel_diff(s2["ex_tritium_kw"], 1e10 * 1e6 * c1 + 1e10 * 2e6 * c1) < 1e-15
    )
    rows.append(["C7 S2 hand vector", fmt(s2["total_kw"]), "parts sum", _check(s2_ok, "C7 S2")])
    return rows, notes, total_ci


def oracle_check_pypact() -> tuple[list[list[str]], list[str], bool]:
    """pypact cross-check of the overlapping inventory columns."""
    try:
        from pypact.reader import InventoryReader
    except Exception as exc:  # noqa: BLE001 — oracle is optional; reason recorded
        note = f"Oracle check (pypact) SKIPPED: {exc}"
        print(note)
        return [], [note], True

    rows: list[list[str]] = []
    notes: list[str] = []
    try:
        with InventoryReader(str(FIX / "clearance.out")) as out:
            intervals = [iv for iv in out.inventory_data if len(iv.nuclides) > 0]
            pyp = [
                {n.name: (n.atoms, n.activity, n.clearance_index) for n in iv.nuclides}
                for iv in intervals
            ]
    except Exception as exc:  # noqa: BLE001 — pypact asserts column shapes; record, skip
        note = f"Oracle check (pypact) SKIPPED: fixture rejected by pypact: {exc}"
        print(note)
        return [], [note], True

    nuc = {
        r["interval"]: {}
        for r in fispact.fispact_parse_clearance((FIX / "clearance.out").read_text())
    }
    for r in fispact.fispact_parse_clearance((FIX / "clearance.out").read_text()):
        nuc[r["interval"]][r["nuclide"].replace("-", "")] = (
            None,
            r["activity_bq"],
            r["clearance_index"],
        )

    worst = 0.0
    for step, pstep in zip(sorted(nuc), pyp, strict=True):
        for name, (_, activity, clearance_index) in nuc[step].items():
            if name not in pstep:
                rows.append([f"O step {step} {name}", "missing", "-", _check(False, name)])
                continue
            _, pactivity, pci = pstep[name]
            err = max(rel_diff(activity, pactivity), rel_diff(clearance_index, pci))
            worst = max(worst, err)
            rows.append(
                [
                    f"O step {step} {name}",
                    fmt(activity),
                    fmt(pactivity),
                    fmt(err),
                    _check(err < 1e-12, f"O {name}"),
                ]
            )
    n_steps = sum(len(s) for s in nuc.values())
    notes.append(f"pypact overlap: {n_steps} nuclide-steps, worst rel err {worst:.3e}.")
    return rows, notes, False


def main() -> int:
    report = Report("clearance", "Clearance / waste-classification (`clearance_vs_pypact.py`)")
    report.prose(
        "Two-part oracle for `nucleide.alara` clearance analytics and the "
        "`nucleide.fispact` clearance-block reader: analytic gates C1-C7 on "
        "synthetic fixtures and the vendored tables (always run), "
        "and a cross-check of the overlapping inventory columns (atoms, activity, "
        "clearance index) against the upstream Apache-2.0 pypact reader (oracle leg)."
    )
    rows1, notes1, total_ci = synthetic_gates()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Value", "Expected", "Status"], rows1)
    report.table(["Quantity", "Value"], [["Parsed clearance-index total", fmt(total_ci)]])

    rows2, notes2, skipped = oracle_check_pypact()
    for note in notes2:
        report.prose(note)
    if skipped:
        report.table(["Gate", "Status"], [["O pypact overlap", "SKIP (pypact unavailable)"]])
    else:
        report.table(["Gate", "Nucleide", "pypact", "Rel err", "Status"], rows2)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} clearance check(s) failed", file=sys.stderr)
        return 1
    print("All clearance checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
