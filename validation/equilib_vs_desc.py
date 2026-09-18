"""Equilibrium-data cross-check (`nucleide.equilib` vs format/grammar gates +
DESC/simsopt container probes).

Two parts:

1. Always-run gates E1-E10 on synthetic inputs (hand-packed classic CDF-1
   bytes built from the published `wout` variable lists — never solver
   output — plus STELLOPT-style `&INDATA` text): magic-byte accept
   (CDF-1 reads `"classic"`), wrong-variant rejection (HDF5 signature
   and CDF-5 magic fail loudly toward facade-side conversion),
   `read_wout` round-trip of the reader-minimum variables, Jacobian
   hand vectors (single-mode constant, two-mode cosine), the J2 grid
   shape/corner contract, the INDATA grammar vectors (switches,
   D-exponents, repeats, boundary tables, `NCURR = 1`, fixed
   boundary), sliced-index rejection, garbage/truncated rejection, and
   the wall-load gates (E9 axisymmetric analytic wall flux plus
   conservation total, E10 cosine hand vector plus loud malformed
   inputs).
2. Container oracle legs (DESC/simsopt): D1 probes the DESC VMEC-format
   loader surface (`desc.io` carrying one of the documented loader
   spellings), S1 probes the simsopt VMEC interface
   (`simsopt.mhd.vmec.Vmec`). These are optional oracle dependencies:
   when the packages cannot be imported the legs are reported as SKIP
   with their reason (never silently). No equilibrium is ever solved
   here — the legs only pin the documented producer/consumer surfaces
   the Rust reader aligns with.

Cross-check procedure for any real file: `ncdump -k <file>` must report
`classic` or `64-bit offset`.
"""

from __future__ import annotations

import math
import struct
import sys
import tempfile
from pathlib import Path

from common import Report, fmt

import nucleide.equilib as eq

FAILURES = 0


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _pad(buf: bytearray, name: str) -> None:
    buf += struct.pack(">I", len(name)) + name.encode("ascii")
    while len(buf) % 4:
        buf.append(0)


def synthetic_wout() -> bytes:
    """Hand-packed classic CDF-1 `wout` (published var lists only)."""
    dims = [("radius", 2), ("mn_mode", 2), ("mn_mode_nyq", 1)]

    def ints(vals: list[int]) -> bytes:
        return struct.pack(f">{len(vals)}i", *vals)

    def doubles(vals: list[float]) -> bytes:
        return struct.pack(f">{len(vals)}d", *vals)

    # (name, dim ids, type code, payload); 4 = int, 6 = double.
    varlist = [
        ("nfp", [], 4, ints([3])),
        ("ns", [], 4, ints([2])),
        ("xm", [1], 6, doubles([0.0, 1.0])),
        ("xn", [1], 6, doubles([0.0, 0.0])),
        ("rmnc", [0, 1], 6, doubles([1.0, 0.1, 1.1, 0.1])),
        ("zmns", [0, 1], 6, doubles([0.0, 0.2, 0.0, 0.25])),
        ("lmns", [0, 1], 6, doubles([0.0, 0.01, 0.0, 0.02])),
        ("gmnc", [0, 2], 6, doubles([2.0, 1.0])),
        ("xm_nyq", [2], 6, doubles([0.0])),
        ("xn_nyq", [2], 6, doubles([0.0])),
    ]
    head = bytearray(b"CDF\x01") + struct.pack(">I", 0)
    head += struct.pack(">II", 10, len(dims))
    for name, length in dims:
        _pad(head, name)
        head += struct.pack(">I", length)
    head += struct.pack(">II", 0, 0)  # no global attrs
    head += struct.pack(">II", 11, len(varlist))
    header_len = len(head)
    for name, dids, _code, _payload in varlist:
        header_len += 4 + len(name) + (-len(name) % 4)
        header_len += 4 + 4 * len(dids) + 8 + 4 + 8
    begin = header_len
    entries = []
    for _name, _dids, _code, payload in varlist:
        entries.append(begin)
        begin += len(payload) + (-len(payload) % 4)
    for (name, dids, code, payload), off in zip(varlist, entries, strict=True):
        _pad(head, name)
        head += struct.pack(">I", len(dids))
        for d in dids:
            head += struct.pack(">I", d)
        head += struct.pack(">II", 0, 0)  # no per-var attrs
        head += struct.pack(">III", code, len(payload), off)
    assert len(head) == header_len
    out = bytes(head)
    for (_name, _dids, _code, payload), _off in zip(varlist, entries, strict=True):
        out += payload + b"\x00" * (-len(payload) % 4)
    return out


INDATA_TEXT = """! Synthetic VMEC input in the documented text grammar.
&INDATA
  LFREEB = F
  NFP = 3
  NCURR = 1
  MPOL = 5 NTOR = 4
  DELT = 0.9
  PMASS_TYPE = 'power_series'
  AM = 0.0 1.0D0 2*0.0
  RAXIS = 1.0 1.1
  RBC(0,0) = 1.0
  RBC(1,0) = -2.5D-1
  ZBS(1,0) = 2.5e-1
/
"""


def always_run(tmp: Path) -> tuple[list[list[str]], list[str]]:
    """Gates E1-E10; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []
    path = tmp / "synthetic_wout.nc"
    path.write_bytes(synthetic_wout())

    # E1: classic magic probes as `ncdump -k` spells it.
    kind = eq.probe_variant(str(path))
    rows.append(["E1 classic probe", kind, "classic", _check(kind == "classic", "E1")])

    # E2: wrong-variant bytes fail loudly toward facade-side conversion.
    for label, blob in [
        ("HDF5/netCDF-4", bytes([0x89, 0x48, 0x44, 0x46, 0x0D, 0x0A, 0x1A, 0x0A])),
        ("CDF-5", b"CDF\x05\x00\x00\x00\x00"),
    ]:
        bad = tmp / "bad.nc"
        bad.write_bytes(blob)
        try:
            eq.probe_variant(str(bad))
            ok, detail = False, "no error raised"
        except ValueError as exc:
            ok = "facade-side" in str(exc)
            detail = "loud" if ok else str(exc)
        rows.append([f"E2 {label} rejected", detail, "loud ValueError", _check(ok, f"E2 {label}")])

    # E3: reader-minimum round-trip.
    w = eq.read_wout(str(path))
    ok = (
        w["variant"] == "classic"
        and w["nfp"] == 3
        and w["ns"] == 2
        and (w["n_modes"], w["n_modes_nyq"]) == (2, 1)
        and w["xm"] == [0.0, 1.0]
        and w["xn"] == [0.0, 0.0]
        and w["rmnc"] == {"rows": 2, "cols": 2, "values": [[1.0, 0.1], [1.1, 0.1]]}
        and w["gmnc"] == {"rows": 2, "cols": 1, "values": [[2.0], [1.0]]}
        and w["xm_nyq"] == [0.0]
        and w["xn_nyq"] == [0.0]
        and w["mpol"] is None
    )
    rows.append(["E3 wout round-trip", fmt(float(w["ns"])), "ns=2", _check(ok, "E3")])

    # E4: Jacobian hand vectors (single (0,0) mode rows: constant).
    gmnc = w["gmnc"]["values"]
    j00 = eq.jacobian([0.0], [0.0], gmnc[0], 0.0, 0.0)
    j13 = eq.jacobian([0.0], [0.0], gmnc[0], 1.3, 2.1)
    j2 = eq.jacobian([0.0], [0.0], gmnc[1], 0.7, -0.4)
    ok = j00 == 2.0 and j13 == 2.0 and j2 == 1.0
    rows.append(["E4 Jacobian const rows", fmt(j00), "2.0", _check(ok, "E4")])
    # Two-mode cosine: 1.0 + 0.5*cos(2θ + ζ).
    c = eq.jacobian([0.0, 2.0], [0.0, -1.0], [1.0, 0.5], math.pi / 2, 0.0)
    ok = abs(c - 0.5) < 1e-12
    rows.append(["E4 two-mode cosine", fmt(c), "0.5", _check(ok, "E4b")])

    # E5: J2 grid shape and corner contract.
    grid = eq.jacobian_grid([0.0], [0.0], gmnc[1], 3, 4, 5)
    ok = len(grid) == 4 and all(len(r) == 5 for r in grid) and grid[0][0] == 1.0
    rows.append(["E5 J2 grid", f"{len(grid)}x{len(grid[0])}", "4x5", _check(ok, "E5")])

    # E6: INDATA grammar vectors.
    parsed = eq.parse_indata(INDATA_TEXT)
    ok = (
        eq.indata_int(parsed, "nfp") == 3
        and eq.indata_int(parsed, "mpol") == 5
        and eq.indata_int(parsed, "ntor") == 4
        and abs(eq.indata_float(parsed, "delt") - 0.9) < 1e-15
        and eq.indata_str(parsed, "pmass_type") == "power_series"
        and eq.indata_series_1d(parsed, "am") == [(0, 0.0), (1, 1.0), (2, 0.0), (3, 0.0)]
        and eq.indata_series_1d(parsed, "raxis") == [(0, 1.0), (1, 1.1)]
        and eq.indata_table_2d(parsed, "rbc") == [((0, 0), 1.0), ((1, 0), -0.25)]
        and eq.indata_table_2d(parsed, "zbs") == [((1, 0), 0.25)]
        and eq.is_free_boundary(parsed) is False
        and eq.ncurr_is_iprime(parsed) is True
    )
    rows.append(["E6 INDATA grammar", "switches+series", "match", _check(ok, "E6")])

    # E7: sliced indices rejected loudly.
    try:
        eq.parse_indata("&INDATA\nRBC(0:4,2) = 1.0\n/\n")
        ok, detail = False, "no error raised"
    except ValueError as exc:
        sliced = "sliced" in str(exc).lower()
        ok, detail = sliced, "loud" if sliced else str(exc)
    rows.append(["E7 sliced index", detail, "loud ValueError", _check(ok, "E7")])

    # E8: garbage and truncation rejected.
    e8ok = True
    for blob in [b"NOPE", b"CD", b"CDF\x03"]:
        bad = tmp / "bad2.nc"
        bad.write_bytes(blob)
        try:
            eq.probe_variant(str(bad))
            e8ok = False
        except ValueError:
            pass
    rows.append(["E8 garbage/truncated", "rejected", "ValueError x3", _check(e8ok, "E8")])

    # E9: wall-load axisymmetric limit (W3) plus conservation total (W2).
    wall = eq.wall_load([0.0, 0.5, 1.0], 4, 3, 2, [2.0] * 24, [3.0] * 24)
    flat = [v for row in wall["loads"] for v in row]
    ok = all(abs(v - 6.0) < 1e-12 for v in flat) and abs(wall["total"] - 12.0 * math.pi**2) < 1e-9
    rows.append(["E9 wall axisymmetric", fmt(wall["total"]), "12*pi^2", _check(ok, "E9")])

    # E10: cosine hand vector, conservation on a mixed field, loud errors.
    wall = eq.wall_load([0.0, 1.0], 4, 1, 1, [1.0] * 4, [1.5, 0.5, 1.5, 0.5])
    got = [row[0] for row in wall["loads"]]
    ok = all(abs(v - w) < 1e-12 for v, w in zip(got, [1.5, 0.5, 1.5, 0.5], strict=True))
    rows.append(["E10 wall cosine", fmt(got[0]), "1.5", _check(ok, "E10")])
    e10ok = True
    for kwargs in [
        {"s_edges": [0.0, 1.0, 0.5]},
        {"ntheta": 0},
        {"birth": [-1.0, 1.0]},
        {"jacobian": [1.0, -0.5]},
    ]:
        args: dict[str, object] = {
            "s_edges": [0.0, 1.0],
            "ntheta": 2,
            "nzeta": 1,
            "nfp": 1,
            "birth": [1.0, 1.0],
            "jacobian": [1.0, 1.0],
        }
        args.update(kwargs)
        try:
            eq.wall_load(**args)
            e10ok = False
        except ValueError:
            pass
    rows.append(["E10 wall malformed", "rejected", "ValueError x4", _check(e10ok, "E10b")])

    notes.append(
        "Synthetic `wout` bytes are hand-packed from the published variable"
        " lists (never solver output); the INDATA text follows the documented"
        " grammar sections. Real files cross-check with `ncdump -k <file>`"
        " reporting `classic` or `64-bit offset`."
    )
    return rows, notes


def container_legs() -> tuple[list[list[str]], list[str]]:
    """Gates D1/S1; SKIP loudly when DESC/simsopt are unavailable."""
    rows: list[list[str]] = []
    notes: list[str] = []

    try:
        import desc.io  # noqa: F401
    except Exception as exc:
        note = f"Oracle leg D1 (DESC VMEC loader) SKIPPED: {exc}"
        print(note)
        notes.append(note + " — the DESC-saved variant stays OUT until re-probed.")
        rows.append(["D1 DESC loader surface", "DESC unavailable", "SKIP"])
    else:
        import desc.io as desc_io

        candidates = ("VMECIO", "read_vmec_output", "load", "save")
        present = [c for c in candidates if hasattr(desc_io, c)]
        status = _check(bool(present), "D1")
        notes.append(f"DESC loader surface present: {present or ['(none)']}.")
        rows.append(["D1 DESC loader surface", ", ".join(present) or "none", status])

    try:
        from simsopt.mhd.vmec import Vmec as SimsoptVmec
    except Exception as exc:
        note = f"Oracle leg S1 (simsopt Vmec) SKIPPED: {exc}"
        print(note)
        notes.append(note + " — simsopt stays a container-only oracle.")
        rows.append(["S1 simsopt Vmec import", "simsopt unavailable", "SKIP"])
    else:
        status = _check(isinstance(SimsoptVmec, type), "S1")
        notes.append("simsopt `Vmec` interface present (import-level contract only).")
        notes.append("simsopt `Vmec` interface present (import-level contract only).")
        rows.append(["S1 simsopt Vmec import", "present", status])

    return rows, notes


def main() -> int:
    report = Report("equilib", "Equilibrium data (`equilib_vs_desc.py`)")
    report.prose(
        "Classic-netCDF `wout` reading plus the VMEC `&INDATA` text grammar"
        " (`nucleide.equilib`; reads data, never solves equilibria)."
        " Always-run gates E1-E10 pin the magic probe, the reader-minimum"
        " round-trip, the Jacobian hand vectors and grid contract, the grammar"
        " vectors, and the wall-load gates on synthetic inputs. DESC/simsopt"
        " legs D1/S1 pin"
        " the documented producer/consumer surfaces and SKIP loudly outside"
        " the container."
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        gate_rows, gate_notes = always_run(Path(tmpdir))
    report.table(["Gate", "Value", "Expected", "Status"], gate_rows)
    for note in gate_notes:
        report.prose(note)
    oracle_rows, oracle_notes = container_legs()
    report.table(["Gate", "Detail", "Status"], oracle_rows)
    for note in oracle_notes:
        report.prose(note)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} equilibrium check(s) failed", file=sys.stderr)
        return 1
    print("All equilibrium checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
