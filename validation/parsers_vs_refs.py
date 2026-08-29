"""Cross-validate Nucleide's file parsers against independent oracle readers.

Oracles (run inside the validation container):
- Serpent `*.m` files: serpentTools (pip-installed; see Containerfile).
- MCNP files: `pyne.mcnp` (Xsdir, SurfSrc, PtracReader; Wwinp/Meshtal need
  PyMOAB, which the nomoab PyNE build lacks — those probes skip loudly).
- FLUKA: `pyne.fluka.Usrbin` reads only binary USRBIN and needs PyMOAB; no
  working oracle exists for our ASCII `.lis` fixtures, so FLUKA is skipped.

Inputs are our own committed fixtures under `fixtures/`; no third-party files.
Every skip is printed and recorded in the report prose.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import numpy as np
from common import Report, fmt, rel_diff

import nucleide

REPO_ROOT = Path(__file__).resolve().parent.parent
SERPENT_DIR = REPO_ROOT / "fixtures" / "serpent"
MCNP_DIR = REPO_ROOT / "fixtures" / "mcnp"
FLUKA_DIR = REPO_ROOT / "fixtures" / "fluka"

# Worst relative difference over every compared numeric field.
WORST = 0.0


def _track(diff: float) -> float:
    global WORST
    WORST = max(WORST, diff)
    return diff


def _note(text: str) -> str:
    """Loud skip note: printed, and wrapped so the rendered Markdown lints."""
    note = textwrap.fill(text, width=100)
    print(note)
    return note


def arr_max_rel_diff(a, b) -> tuple[float, int]:
    """Max elementwise relative difference between two array-likes."""
    va = np.asarray(a, dtype=float).reshape(-1)
    vb = np.asarray(b, dtype=float).reshape(-1)
    if va.size != vb.size:
        raise ValueError(f"length mismatch: {va.size} vs {vb.size}")
    return max((_track(rel_diff(x, y)) for x, y in zip(va, vb, strict=True)), default=0.0), int(
        va.size
    )


def serpent_kind(path: Path) -> str | None:
    """Infer the Serpent file kind from its suffix."""
    name = path.name
    for suffix, kind in (("_res.m", "res"), ("_dep.m", "dep"), ("_det.m", "det")):
        if name.endswith(suffix):
            return kind
    return None


def compare_serpent_dep(path: Path, st_reader) -> list[list[str]]:
    """Compare depletion-file fields shared by Nucleide and serpentTools."""
    nuc = nucleide.serpent.read_serpent(str(path), "dep")
    rows: list[list[str]] = []

    # Scalar/vector metadata: ZAI, DAYS, BU.
    md = st_reader.metadata
    d, n = arr_max_rel_diff(nuc["ZAI"], md["zai"])
    rows.append(["ZAI", str(n), fmt(d)])
    d, n = arr_max_rel_diff(nuc["DAYS"], md["days"])
    rows.append(["DAYS", str(n), fmt(d)])
    d, n = arr_max_rel_diff(nuc["BU"], md["burnup"])
    rows.append(["BU (burnup)", str(n), fmt(d)])
    names_nuc = [str(x).strip() for x in nuc["NAMES"]]
    names_st = [str(x).strip() for x in md["names"]]
    n_mismatch = sum(1 for a, b in zip(names_nuc, names_st, strict=False) if a != b)
    n_mismatch += abs(len(names_nuc) - len(names_st))
    rows.append(["NAMES", str(len(names_nuc)), f"{n_mismatch} mismatches"])
    if n_mismatch:
        _track(1.0)

    # Per-material arrays present in both readers.
    fields = {"ADENS": "adens", "MDENS": "mdens", "VOLUME": "volume"}
    for mat_name, st_mat in sorted(st_reader.materials.items()):
        prefix = f"MAT_{mat_name}_" if mat_name != "total" else "TOT_"
        for nuc_key, st_key in sorted(fields.items()):
            full_key = prefix + nuc_key
            if full_key not in nuc or st_key not in st_mat.data:
                continue
            d, n = arr_max_rel_diff(nuc[full_key], st_mat.data[st_key])
            rows.append([full_key, str(n), fmt(d)])
    return rows


def compare_serpent_det(path: Path, st_reader) -> list[list[str]]:
    """Compare detector tallies, errors and grids shared by both readers."""
    nuc = nucleide.serpent.read_serpent(str(path), "det")
    rows: list[list[str]] = []
    for det_name, st_det in sorted(st_reader.detectors.items()):
        base = f"DET{det_name}"
        if base not in nuc:
            continue
        bins = np.asarray(nuc[base], dtype=float).reshape(-1, 13)
        d, n = arr_max_rel_diff(bins[:, 11], st_det.tallies)
        rows.append([f"{base} tallies", str(n), fmt(d)])
        d, n = arr_max_rel_diff(bins[:, 12], st_det.errors)
        rows.append([f"{base} rel errors", str(n), fmt(d)])
        for grid, st_grid in sorted(st_det.grids.items()):
            key = base + grid
            if key not in nuc:
                continue
            d, n = arr_max_rel_diff(nuc[key], st_grid)
            rows.append([f"{key} grid", str(n), fmt(d)])
    return rows


def serpent_section(report: Report) -> None:
    """Serpent fixtures vs serpentTools."""
    report.heading("Serpent vs serpentTools")
    try:
        import serpentTools
    except ImportError:
        note = (
            "SKIPPED: serpentTools is not installed in this environment; the Serpent "
            "oracle comparison did not run."
        )
        report.prose(_note(note))
        return

    st_reader_kind = {"res": "results", "dep": "dep", "det": "det"}
    report.prose(
        f"serpentTools {serpentTools.__version__} readers run on `fixtures/serpent/*.m`."
        "\nCompared fields — dep: `ZAI`, `DAYS`, `BU`/`burnup`, `NAMES`, per-material"
        " `ADENS`,\n`MDENS`, `VOLUME`; det: per-detector tally and relative-error columns plus"
        " shared\nbin grids (`E`, `T`, `X`, `Y`)."
    )
    rows: list[list[str]] = []
    skips: list[str] = []
    for path in sorted(SERPENT_DIR.glob("*.m")):
        kind = serpent_kind(path)
        if kind is None:
            continue
        try:
            reader = serpentTools.read(str(path), st_reader_kind[kind])
        except Exception as exc:
            note = (
                f"SKIPPED {path.name}: serpentTools {serpentTools.__version__} cannot read it "
                f"({type(exc).__name__}: {exc})"
            )
            skips.append(_note(note))
            continue
        if kind == "dep":
            for field, n, diff in compare_serpent_dep(path, reader):
                rows.append([path.name, field, n, diff])
        elif kind == "det":
            for field, n, diff in compare_serpent_det(path, reader):
                rows.append([path.name, field, n, diff])
        else:
            # res: nucleide side parses fine; compare nothing without an oracle
            # (serpentTools ResultsReader postprocessing fails on these files,
            # caught above).
            pass
    if rows:
        report.table(["File", "Field", "Values compared", "Max rel diff"], rows)
    for note in skips:
        report.prose(note)


def compare_xsdir() -> list[list[str]]:
    """Compare xsdir table entries against pyne.mcnp.Xsdir."""
    from pyne import mcnp

    path = MCNP_DIR / "xsdir" / "dummy_xsdir"
    nuc = nucleide.mcnp.read_xsdir(str(path))
    ref = mcnp.Xsdir(str(path))
    rows: list[list[str]] = []
    rows.append(["table count", str(len(nuc.tables)), fmt(abs(len(nuc.tables) - len(ref.tables)))])
    nuc_tables = {t.name: t for t in nuc.tables}
    for rt in ref.tables:
        nt = nuc_tables.get(rt.name)
        if nt is None:
            rows.append([rt.name, "—", "MISSING in Nucleide"])
            _track(1.0)
            continue
        str_ok = nt.filename == rt.filename and nt.ptable == rt.ptable
        int_ok = (
            nt.filetype == rt.filetype
            and nt.address == rt.address
            and nt.tablelength == rt.tablelength
        )
        d, _ = arr_max_rel_diff([nt.awr, nt.temperature or 0.0], [rt.awr, rt.temperature or 0.0])
        status = fmt(d) if (str_ok and int_ok) else "field mismatch"
        if not (str_ok and int_ok):
            _track(1.0)
        rows.append([f"{rt.name} (awr, temperature, name/dir/ints)", "2 + 5", status])
    return rows


def compare_surfsrc() -> tuple[list[list[str]], list[str]]:
    """Compare SSW surface-source files against pyne.mcnp.SurfSrc."""
    from pyne import mcnp

    rows: list[list[str]] = []
    skips: list[str] = []
    fields = ["nps", "wgt", "erg", "tme", "x", "y", "z", "u", "v", "w", "cs"]
    for path in sorted(MCNP_DIR.glob("ssw/*.w")):
        try:
            ref = mcnp.SurfSrc(str(path), "rb")
            ref.read_header()
            ref.read_tracklist()
        except Exception as exc:
            note = f"SKIPPED {path.name}: PyNE SurfSrc failed ({type(exc).__name__}: {exc})"
            skips.append(_note(note))
            continue
        nuc = nucleide.mcnp.read_ssw(str(path))
        header_ok = (
            nuc.kod.strip() == ref.kod.strip()
            and nuc.ver.strip() == ref.ver.strip()
            and nuc.np1 == ref.np1
            and nuc.nrss == ref.nrss
            and abs(nuc.ncrd) == abs(ref.ncrd)
            and nuc.njsw == ref.njsw
            and nuc.niss == ref.niss
        )
        nuc_tracks = nuc.tracks()
        worst = 0.0
        n_vals = 0
        for field in fields:
            a = [t[field] for t in nuc_tracks]
            b = [getattr(t, field) for t in ref.tracklist]
            d, n = arr_max_rel_diff(a, b)
            worst = max(worst, d)
            n_vals += n
        status = fmt(worst) if header_ok else "header mismatch"
        if not header_ok:
            _track(1.0)
        rows.append([f"{path.name} ({len(nuc_tracks)} tracks)", str(n_vals), status])
    return rows, skips


def compare_ptrac() -> tuple[list[list[str]], list[str]]:
    """Compare PTRAC headers against pyne.mcnp.PtracReader.

    PyNE ships no `PtracFile` class; `PtracReader` is its low-level binary
    helper, so only the header scalars it exposes (problem title, per-record
    variable counts) are compared.
    """
    from pyne.mcnp import PtracReader

    rows: list[list[str]] = []
    skips: list[str] = []
    for path in sorted(MCNP_DIR.glob("ptrac/*.ptrac")):
        try:
            ref = PtracReader(str(path))
        except Exception as exc:
            note = f"SKIPPED {path.name}: PyNE PtracReader failed ({type(exc).__name__}: {exc})"
            skips.append(_note(note))
            continue
        nuc = nucleide.mcnp.read_ptrac(str(path))
        ref_title = ref.problem_title
        if isinstance(ref_title, bytes):
            ref_title = ref_title.decode(errors="replace")
        title_ok = nuc.problem_title.strip() == ref_title.strip()
        nums_ok = dict(nuc.variable_nums) == dict(ref.variable_nums)
        n_bad = (not title_ok) + (not nums_ok)
        if n_bad:
            _track(1.0)
        rows.append([path.name, "problem_title + variable_nums", "OK" if not n_bad else "MISMATCH"])
    return rows, skips


def mcnp_section(report: Report) -> None:
    """MCNP fixtures vs pyne.mcnp."""
    report.heading("MCNP vs PyNE")
    report.prose(
        "PyNE 0.7.5 (`nomoab` build) `pyne.mcnp` readers run on the committed"
        "\n`fixtures/mcnp/` files. Compared fields — xsdir: per-table `name`, `awr`,"
        "\n`filename`, `filetype`, `address`, `tablelength`, `temperature`, `ptable`;"
        "\nssw: header (`kod`, `ver`, `np1`, `nrss`, `ncrd`, `njsw`, `niss`) and per-track"
        "\n`nps`, `wgt`, `erg`, `tme`, `x`, `y`, `z`, `u`, `v`, `w`, `cs` payloads;"
        "\nptrac: problem title and per-record variable counts (PyNE has no `PtracFile`"
        "\nclass; its low-level `PtracReader` exposes only headers)."
    )

    rows: list[list[str]] = []
    for name, fileds, status in compare_xsdir():
        rows.append(["xsdir", name, fileds, status])
    ssw_rows, ssw_skips = compare_surfsrc()
    for name, n_vals, status in ssw_rows:
        rows.append(["ssw", name, n_vals, status])
    ptrac_rows, ptrac_skips = compare_ptrac()
    for name, what, status in ptrac_rows:
        rows.append(["ptrac", name, what, status])
    report.table(["Format", "Item", "Values compared", "Max rel diff / status"], rows)

    skips = ssw_skips + ptrac_skips
    for what, reason in [
        ("wwinp", "`pyne.mcnp.Wwinp` requires PyMOAB to build its mesh (nomoab build)"),
        ("meshtal", "`pyne.mcnp.Meshtal` requires PyMOAB (nomoab build)"),
    ]:
        for path in sorted(MCNP_DIR.glob(f"{what}/*")):
            note = f"SKIPPED {path.name}: {reason}"
            skips.append(_note(note))
    for note in skips:
        report.prose(note)


def fluka_section(report: Report) -> None:
    """FLUKA fixtures: no working oracle."""
    report.heading("FLUKA vs PyNE")
    note = (
        "SKIPPED: no working FLUKA oracle exists in this environment. `pyne.fluka.Usrbin`"
        " requires PyMOAB (absent in the nomoab build) and reads only binary USRBIN output,"
        " while Nucleide's committed fixtures are ASCII `.lis` files"
        f" ({', '.join(p.name for p in sorted(FLUKA_DIR.glob('*.lis')))})."
    )
    report.prose(_note(note))


def main() -> int:
    report = Report("parsers", "Parser cross-validation (`parsers_vs_refs.py`)")
    report.prose(
        "Nucleide's readers are cross-checked against independent oracle readers on the"
        "\nsame committed fixture files. Skipped comparisons (missing or incapable oracle)"
        "\nare listed explicitly below; none are silent."
    )
    serpent_section(report)
    mcnp_section(report)
    fluka_section(report)
    report.emit()

    if WORST > 1.0e-6:
        print(f"FAIL: parser oracle mismatch up to {WORST:.3e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
