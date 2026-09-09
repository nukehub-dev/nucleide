#!/usr/bin/env python3
"""Regenerate the screening-level nuclear-data TSVs in ``crates/nuclei/src/data/``.

Each table is derived from a primary evaluated source — no hand-copied values:

- ``simple_xs.tsv``: thermal totals from the NIST NCNR bound cross
  sections converted to free-atom values (see below) and 14-MeV totals from
  ENDF/B-VII.1 MF3/MT1 (File 3 is complete above the resonance range).
  Rows need both sources, so isomers (no neutron tapes) are absent.
  Upstream: https://www.ncnr.nist.gov/resources/n-lengths/ and
  https://www.nndc.bnl.gov/endf/
- ``scattering_lengths.tsv``: NIST NCNR bound scattering lengths (Sears,
  Neutron News 3(3), 1992), parsed with the same row pattern PyNE's
  ``pyne/dbgen/scattering_lengths.py`` uses.  Upstream:
  https://www.ncnr.nist.gov/resources/n-lengths/
- ``decay_energy.tsv``: mean *prompt* recoverable energy per decay from
  ENDF/B-VII.1 decay tapes (MF8/MT457 average-energy summary:
  ``E = v0 + v2 + v4`` of the six-value record, i.e. the mean-energy
  components with the paired uncertainties skipped).  Prompt means the
  daughter's gammas belong to the daughter's row (e.g. the 662 keV line
  lives on Ba137_m1, not Cs137) — chain codes sum members.  Isomers with
  their own tapes (``*m1``/``*m2``) get their own rows.
- ``dose_factors.tsv``: external-air/soil, ingestion, and inhalation dose
  factors for 93 nuclides from the BSD-3 PyNE ``dbgen/dosefactors*.csv``
  tables (HNF-SD-WM-TI-707 Rev.1 / HNF-5636 App. O; GENII/EPA/DOE are 3
  parallel evaluations).  ``+D`` (plus daughters) folds into the parent.
  Air is EPA-only (GENII/DOE rows are ``-1`` sentinels, matching PyNE).

Usage::

    python3 scripts/gen-nuclear-data.py \
        --endf-neutrons /path/to/endf-b-vii.1/neutrons \
        --endf-decay /path/to/endf-b-vii.1/decay \
        --nist-html /path/to/scattering_lengths.html \
        --dose-air /tmp/dosefactors_external_air.csv \
        --dose-soil /tmp/dosefactors_external_soil.csv \
        --dose-ingest /tmp/dosefactors_ingest.csv \
        --dose-inhale /tmp/dosefactors_inhale.csv \
        --out crates/nuclei/src/data

All inputs are local checkouts (see help for download URLs); nothing is
fetched over the network.  Exit nonzero if any hard spot-check fails.
Pass only the ``--dose-*`` flags (plus ``--out``) for a dose-only run;
pass only the ENDF/NIST flags for the legacy tables.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys

THERMAL_EV = 0.0253
FAST_EV = 14.0e6

# Spot-checks: (table, name, lo, hi, unit).  Failure aborts the run.
SPOTS = [
    ("decay", "Co60", 2.55, 2.65, "MeV"),
    ("decay", "Cs137", 0.15, 0.22, "MeV"),  # prompt only; 662 keV is Ba137_m1's
    ("decay", "Ba137_m1", 0.60, 0.72, "MeV"),
    ("decay", "H3", 0.005, 0.0065, "MeV"),
    ("decay", "Sr90", 0.19, 0.21, "MeV"),
    ("simple_xs", "H1", 20.0, 21.5, "b thermal"),
    ("simple_xs", "U235", 680.0, 710.0, "b thermal"),
    ("scattering", "H1", -3.7407, -3.7405, "fm coherent"),
    ("scattering", "O16", 5.802, 5.804, "fm coherent"),
]

SYMBOLS = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)
Z_OF = {s: i + 1 for i, s in enumerate(SYMBOLS)}

FLOAT_RE = re.compile(r"([+-]?\d+\.\d+)([+-]\d+)")


def f11(raw: str) -> float:
    """Parse an 11-char ENDF floating field (exponent without ``e``)."""
    return float(FLOAT_RE.sub(r"\1e\2", raw.strip()))


def read_lines(path: str) -> list[str]:
    with open(path) as fh:
        return fh.read().splitlines()


def section_lines(lines: list[str], mf: int, mt: int) -> list[str]:
    """Return 66-char data bodies (no MAT/MF/MT tail) of one ENDF section."""
    out: list[str] = []
    in_section = False
    for line in lines:
        if len(line) < 75:
            continue
        try:
            lmf, lmt = int(line[70:72]), int(line[72:75])
        except ValueError:
            continue
        if (lmf, lmt) == (mf, mt):
            in_section = True
            out.append(line[:66])
        elif in_section:
            break
    return out


def fields(body: str) -> list[float]:
    """Split one 66-char ENDF body into its 11-char floats.

    Fixed columns (not whitespace) because high-precision evaluations omit
    the blank between adjacent fields (``1.164410e+2-1.271980e+2``).
    Trailing blank slots (short INT records) are skipped.
    """
    out = []
    for i in range(0, 66, 11):
        slot = body[i : i + 11]
        if slot.strip():
            out.append(f11(slot))
    return out


def chunks(fields: list[float], n: int):
    for i in range(0, len(fields), n):
        yield fields[i : i + n]


def tab1_interpolate(xs: list[str], energy: float) -> float:
    """Interpolate an MF3 TAB1 record at ``energy`` honoring INT laws 1-5."""
    head = fields(xs[1])
    nr, npts = int(head[4]), int(head[5])
    int_lines = xs[2 : 2 + (2 * nr + 5) // 6]
    ints: list[int] = []
    for chunk in chunks([v for line in int_lines for v in fields(line)], 2):
        ints.extend([int(chunk[0]), int(chunk[1])])
    nbt = ints[0::2]
    law = ints[1::2]
    pairs = [v for line in xs[2 + (2 * nr + 5) // 6 :] for v in fields(line)]
    grid = pairs[0::2]
    vals = pairs[1::2]
    assert len(grid) == npts, f"NP mismatch: {len(grid)} vs {npts}"
    if not grid[0] <= energy <= grid[-1]:
        raise ValueError(f"E={energy} outside grid [{grid[0]}, {grid[-1]}]")
    seg = next(i for i, b in enumerate(nbt) if energy <= grid[b - 1])
    kind = law[seg]
    i = next(j for j in range(len(grid) - 1) if grid[j] <= energy <= grid[j + 1])
    x0, x1, y0, y1 = grid[i], grid[i + 1], vals[i], vals[i + 1]
    if energy in (x0, x1) or y0 == y1:
        return y0 if energy == x0 else y1
    if kind == 1:
        return y0
    if kind == 2:
        return y0 + (y1 - y0) * (energy - x0) / (x1 - x0)
    if kind == 3:
        return y0 + (y1 - y0) * math.log(energy / x0) / math.log(x1 / x0)
    if kind == 4:
        return y0 * math.exp(math.log(y1 / y0) * (energy - x0) / (x1 - x0))
    if kind == 5:
        return y0 * (y1 / y0) ** (math.log(energy / x0) / math.log(x1 / x0))
    raise ValueError(f"unsupported INT law {kind}")


TAPE_RE = re.compile(r"^(?:dec|n)-(\d+)_([A-Za-z]+)_(\d+)(m\d+)?\.endf$")


def tape_id(fname: str) -> tuple[int, str, int, int] | None:
    """Return (Z, symbol, A, state) from a tape filename, or None."""
    m = TAPE_RE.match(os.path.basename(fname))
    if not m:
        return None
    _zdir, sym, a, msuf = m.groups()
    z = Z_OF.get(sym)
    if z is None or z != int(_zdir):
        return None
    return z, sym, int(a), int(msuf[1:]) if msuf else 0


def gnds(sym: str, a: int, state: int) -> str:
    name = f"{sym}{a}"
    return f"{name}_m{state}" if state else name


def gen_simple_xs(
    neutron_dir: str, nist: dict[str, dict[str, float]]
) -> tuple[dict[str, tuple[float, float]], list[str]]:
    """Thermal totals from NIST bound XS (free-atom corrected) + ENDF fast.

    ENDF File 3 alone cannot supply thermal totals when the evaluation has
    a resonance range: there File 3 holds background only (e.g. U235 MT1 is
    0 below ~77 eV) and the resonance contribution lives in File 2.  So:

    - tapes with no MF2 section: ENDF MF3/MT1 at 0.0253 eV (exact);
    - tapes with MF2: NIST NCNR bound XS converted to free-atom via
      ``xs*(A/(A+1))^2 + xs_a`` (rows without a NIST entry are skipped).

    The 14-MeV column always comes from ENDF/B-VII.1 MF3/MT1 (File 3 is
    complete there), interpolated per the tape's own INT laws.
    """
    rows: dict[str, tuple[float, float]] = {}
    log: list[str] = []
    n_endf = n_nist = 0
    for fname in sorted(os.listdir(neutron_dir)):
        ident = tape_id(fname)
        if ident is None or not fname.startswith("n-"):
            continue
        z, sym, a, state = ident
        if a == 0:
            log.append(f"skip elemental evaluation tape: {fname}")
            continue
        if state:
            log.append(f"skip isomer neutron tape (no isomer XS): {fname}")
            continue
        name = gnds(sym, a, 0)
        lines = read_lines(os.path.join(neutron_dir, fname))
        has_mf2 = any(len(line) >= 75 and line[70:72] == " 2" for line in lines)
        thermal: float | None = None
        if not has_mf2:
            sec = section_lines(lines, 3, 1)
            if len(sec) >= 3:
                try:
                    thermal = tab1_interpolate(sec, THERMAL_EV)
                    n_endf += 1
                except (ValueError, AssertionError, IndexError) as exc:
                    log.append(f"note {fname}: ENDF thermal failed ({exc}); trying NIST")
        if thermal is None:
            entry = nist.get(name)
            if entry is None or entry["xs"] is None:
                log.append(f"skip {fname}: resonance tape, no NIST row for {name}")
                continue
            thermal = entry["xs"] * (a / (a + 1.0)) ** 2 + (entry["xs_a"] or 0.0)
            n_nist += 1
        sec = section_lines(lines, 3, 1)
        if len(sec) < 3:
            log.append(f"skip {fname}: no MF3/MT1 section")
            continue
        try:
            za = int(fields(sec[0])[0])
            assert za == z * 1000 + a, f"ZA mismatch in {fname}"
            fa = tab1_interpolate(sec, FAST_EV)
        except (ValueError, AssertionError, IndexError) as exc:
            log.append(f"skip {fname}: {exc}")
            continue
        rows[name] = (thermal, fa)
    log.append(f"thermal sources: ENDF-direct={n_endf} NIST-derived={n_nist}")
    return rows, log


def gen_decay(decay_dir: str) -> tuple[dict[str, float], dict[str, str], list[str]]:
    rows: dict[str, float] = {}
    half: dict[str, str] = {}
    log: list[str] = []
    n_e5 = 0
    for fname in sorted(os.listdir(decay_dir)):
        ident = tape_id(fname)
        if ident is None or not fname.startswith("dec-"):
            continue
        z, sym, a, state = ident
        lines = read_lines(os.path.join(decay_dir, fname))
        sec = section_lines(lines, 8, 457)
        if len(sec) < 3:
            log.append(f"skip {fname}: no MF8/MT457 section")
            continue
        try:
            head = fields(sec[0])
            assert int(head[0]) == z * 1000 + a, f"ZA mismatch in {fname}"
            if int(head[2]) != state:
                log.append(f"note {fname}: LIS={int(head[2])} vs filename m{state}")
            tline = fields(sec[1])
            assert int(tline[4]) == 6, f"summary record N1 != 6 in {fname}"
            vals = fields(sec[2])[:6]
        except (AssertionError, IndexError, ValueError) as exc:
            log.append(f"skip {fname}: {exc}")
            continue
        if vals[4] != 0.0:
            n_e5 += 1
            log.append(f"note {fname}: third component E5={vals[4]:.6g} eV included")
        heat_mev = (vals[0] + vals[2] + vals[4]) / 1e6
        if heat_mev == 0.0:
            # All-zero summary: stable nuclide (or no heat data). Omitted so
            # that lookups return None, matching the half-life table's
            # "stable nuclides are simply absent" convention.
            continue
        rows[gnds(sym, a, state)] = heat_mev
        half[gnds(sym, a, state)] = f"{tline[0]:.6g} s"
    log.append(f"third-component (E5) nonzero in {n_e5} tapes")
    return rows, half, log


CELL = r"((?:[^<]|<i>|</i>)*?)"
NIST_ROW_RE = re.compile(
    r"<td>\s*(?P<iso>[A-Za-z\d]+)\s*"
    rf"<td>\s*(?P<conc>{CELL})\s*"
    rf"<td>\s*(?P<b_coherent>{CELL})\s*"
    rf"<td>\s*(?P<b_incoherent>{CELL})\s*"
    rf"<td>\s*(?P<xs_coherent>{CELL})\s*"
    rf"<td>\s*(?P<xs_incoherent>{CELL})\s*"
    rf"<td>\s*(?P<xs>{CELL})\s*"
    rf"<td>\s*(?P<xs_a>{CELL})\s*<tr>"
)

ISO_MASS_FIRST_RE = re.compile(r"^(\d+)([A-Za-z]+)(m\d+)?$")
ISO_SYM_FIRST_RE = re.compile(r"^([A-Za-z]+)(\d+)(m\d+)?$")


def split_iso(label: str) -> tuple[str, int, int] | None:
    """Split a NIST isotope label (``16O`` or ``O16``) into (sym, A, state)."""
    m = ISO_MASS_FIRST_RE.match(label.strip())
    if m:
        a, sym, msuf = m.groups()
        return sym, int(a), int(msuf[1:]) if msuf else 0
    m = ISO_SYM_FIRST_RE.match(label.strip())
    if m:
        sym, a, msuf = m.groups()
        return sym, int(a), int(msuf[1:]) if msuf else 0
    return None  # element rows (U, W, ...) carry no mass number


def monoisotopic_a(abundance_tsv: str) -> dict[str, int]:
    """Map element symbol to mass number for strictly monoisotopic elements.

    Derived from the vendored natural-abundance table: an element qualifies
    only if exactly one mass number occurs (isomers ignored).  This grounds
    natural-element NIST rows (which carry no A) without a hand-kept list.
    """
    masses: dict[str, set[int]] = {}
    with open(abundance_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            name = line.split()[0]
            m = re.match(r"([A-Za-z]+)(\d+)", name)
            if m is None:
                continue
            masses.setdefault(m.group(1), set()).add(int(m.group(2)))
    return {sym: next(iter(v)) for sym, v in masses.items() if len(v) == 1}


def nist_num(raw: str) -> float | complex | None:
    """Mirror PyNE's ``nist_num`` without ``eval``: strip "(...)" qualifiers."""
    text = re.sub(r"\(.*?\)", "", raw).replace("<i>i</i>", "j").strip()
    if text in ("---", ""):
        return None
    text = text.lstrip("+-") if text.startswith(("+/-", "-/+")) else text
    try:
        return complex(text) if "j" in text else float(text)
    except ValueError:
        return None


def gen_scattering(
    html_path: str,
) -> tuple[
    dict[str, tuple[float, float]],
    dict[str, dict[str, float | None]],
    dict[str, dict[str, float | None]],
    list[str],
]:
    """Parse the NIST table; also return bound ``xs``/``xs_a`` for thermal.

    Returns (isotope_rows, isotope_extra, elem100_extra, log) where
    ``elem100_extra`` holds natural-element rows tabulated at ~100%
    abundance — those elements are monoisotopic, so callers may attribute
    the row to the matching ENDF tape isotope (logged per row).
    """
    rows: dict[str, tuple[float, float]] = {}
    extra: dict[str, dict[str, float | None]] = {}
    elem100: dict[str, dict[str, float | None]] = {}
    log: list[str] = []
    n_imag = 0
    with open(html_path) as fh:
        html = fh.read()
    for m in NIST_ROW_RE.finditer(html):
        col = m.groupdict()
        label = col["iso"].strip()
        split = split_iso(label)
        b_coh = nist_num(col["b_coherent"])
        b_inc = nist_num(col["b_incoherent"])
        xs_inc = nist_num(col["xs_incoherent"])
        xs = nist_num(col["xs"])
        xs_a = nist_num(col["xs_a"])
        if b_coh is None:
            log.append(f"skip {label}: no coherent length")
            continue
        if isinstance(b_coh, complex):
            if b_coh.imag != 0.0:
                n_imag += 1
            b_coh = b_coh.real
        if isinstance(b_inc, complex):
            b_inc = b_inc.real
        if b_inc is None:
            xs_i = xs_inc.real if isinstance(xs_inc, complex) else xs_inc
            # No incoherent data anywhere: scattered none.
            b_inc = 0.0 if xs_i is None or xs_i <= 0.0 else 10.0 * math.sqrt(xs_i / (4.0 * math.pi))
        bundle = {
            "b_coh": float(b_coh),
            "b_inc": float(b_inc),
            "xs": _real_or_none(xs),
            "xs_a": _real_or_none(xs_a),
        }
        if split is None:
            conc = nist_num(col["conc"])
            if isinstance(conc, complex):
                conc = conc.real
            if conc is not None and abs(conc - 100.0) < 0.5 and label in Z_OF:
                elem100[label] = bundle
            continue
        sym, a, state = split
        if sym not in Z_OF:
            log.append(f"skip unmapped isotope label: {label}")
            continue
        name = gnds(sym, a, state)
        rows[name] = (bundle["b_coh"], bundle["b_inc"])
        extra[name] = bundle
    log.append(f"complex coherent lengths (real part kept): {n_imag} rows")
    return rows, extra, elem100, log


def _real_or_none(v: float | complex | None) -> float | None:
    if v is None:
        return None
    return v.real if isinstance(v, complex) else v


# Dose-factor spots: (GNDS name, pathway, source, expected DF).
# Values are the raw table factors (mrem/h per Ci/m^3 for air,
# mrem/h per Ci/m^2 for soil, mrem/pCi for ingest/inhale), copied from the
# upstream PyNE CSVs. Failure aborts the run.
DOSE_SPOTS = [
    ("Co60", "ingest", "EPA", 2.69e-05),
    ("Cs137", "inhale", "EPA", 3.19e-05),  # upstream Cs-137+D folds to Cs137
    ("H3", "air", "EPA", 4.41e-012),
    ("K40", "soil", "EPA", 4.33e02),
]

DOSE_PATHWAYS = ("air", "soil", "ingest", "inhale")
DOSE_SOURCES = ("EPA", "DOE", "GENII")


def pyne_to_gnds(raw: str) -> str:
    """Convert a PyNE dose-table nuclide (``H-3``, ``Cs-137+D``) to GNDS.

    ``+D`` (plus daughters) folds into the parent per ``dosefactors.py``
    ``read_row``; dashes are dropped and a trailing metastable ``m`` becomes
    ``_m1`` (all dose-table isomers are first isomers).
    """
    base = raw.strip()
    if base.endswith("+D"):
        base = base[:-2]
    base = base.replace("-", "")
    if base and base[-1] in ("m", "M") and len(base) >= 2 and base[-2].isdigit():
        base = base[:-1] + "_m1"
    m = re.match(r"^([A-Za-z]+)(.*)$", base)
    if m:
        sym, rest = m.groups()
        sym = sym[0].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()
        base = sym + rest
    return base


def _dose_csv_rows(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"{path}: missing header row")
    return list(reader.fieldnames), list(reader)


def gen_dose(
    air_path: str, soil_path: str, ingest_path: str, inhale_path: str
) -> tuple[dict[tuple[str, str, str], tuple[float, float | None, str]], list[str]]:
    """Build dose-factor rows from four local PyNE CSV copies.

    Returns ({(gnds, pathway, source): (factor, f1, lung_model)}, log) where
    ``f1`` is set only on ingest rows and ``lung_model`` only on inhale rows
    (empty string elsewhere).  Air is EPA-only: GENII/DOE air factors are
    ``-1`` sentinels, matching PyNE's ``grab_dose_factors`` (which stores
    ``-1`` for missing GENII/DOE air and ratio).
    """
    rows: dict[tuple[str, str, str], tuple[float, float | None, str]] = {}
    log: list[str] = []

    _, air = _dose_csv_rows(air_path)
    for r in air:
        name = pyne_to_gnds(r["Nuclide"])
        factor = float(r["Air Dose Rate Factor"])
        rows[(name, "air", "EPA")] = (factor, None, "")
        rows[(name, "air", "GENII")] = (-1.0, None, "")
        rows[(name, "air", "DOE")] = (-1.0, None, "")

    _, soil = _dose_csv_rows(soil_path)
    for r in soil:
        name = pyne_to_gnds(r["Nuclide"])
        for src in DOSE_SOURCES:
            rows[(name, "soil", src)] = (float(r[src]), None, "")

    _, ingest = _dose_csv_rows(ingest_path)
    for r in ingest:
        name = pyne_to_gnds(r["Nuclide"])
        f1 = float(r["f1"])
        for src in DOSE_SOURCES:
            rows[(name, "ingest", src)] = (float(r[src]), f1, "")

    _, inhale = _dose_csv_rows(inhale_path)
    for r in inhale:
        name = pyne_to_gnds(r["Nuclide"])
        lung = r["Lung Model"].strip()
        for src in DOSE_SOURCES:
            rows[(name, "inhale", src)] = (float(r[src]), None, lung)

    names = sorted({k[0] for k in rows})
    log.append(f"dose nuclides (folded +D): {len(names)}")
    log.append(f"dose rows (nuclide x pathway x source): {len(rows)}")
    # The air table carries Sb-125 where the other three carry Sb-125+D;
    # both fold to Sb125, so the four inputs must agree after folding.
    for label, got in (("air", air), ("soil", soil), ("ingest", ingest), ("inhale", inhale)):
        folded = sorted({pyne_to_gnds(r["Nuclide"]) for r in got})
        if folded != names:
            missing = sorted(set(names) - set(folded))
            extra = sorted(set(folded) - set(names))
            log.append(f"note {label}: nuclide set differs (missing={missing} extra={extra})")
    return rows, log


def write_dose_tsv(
    path: str,
    header: list[str],
    rows: dict[tuple[str, str, str], tuple[float, float | None, str]],
) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join("# " + h for h in header) + "\n")
        for name, pathway, source in sorted(rows):
            factor, f1, lung = rows[(name, pathway, source)]
            f1_txt = "" if f1 is None else f"{f1:.6g}"
            fh.write(f"{name}\t{pathway}\t{source}\t{factor:.6g}\t{f1_txt}\t{lung}\n")


def write_tsv(path: str, header: list[str], rows: dict, fmt) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join("# " + h for h in header) + "\n")
        for name in sorted(rows):
            fh.write(f"{name}\t" + fmt(rows[name]) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--endf-neutrons",
        required=False,
        default=None,
        help="ENDF/B-VII.1 neutron tapes dir (download: https://www.nndc.bnl.gov/endf/)",
    )
    ap.add_argument(
        "--endf-decay",
        required=False,
        default=None,
        help="ENDF/B-VII.1 decay tapes dir (same source)",
    )
    ap.add_argument(
        "--nist-html",
        required=False,
        default=None,
        help="NIST scattering-lengths page copy "
        "(download: https://www.ncnr.nist.gov/resources/n-lengths/; "
        "a copy ships in PyNE as pyne/dbgen/scattering_lengths.html)",
    )
    ap.add_argument(
        "--dose-air",
        required=False,
        default=None,
        help="Local copy of PyNE dosefactors_external_air.csv "
        "(download: https://raw.githubusercontent.com/pyne/pyne/develop/pyne/dbgen/dosefactors_external_air.csv)",
    )
    ap.add_argument(
        "--dose-soil",
        required=False,
        default=None,
        help="Local copy of PyNE dosefactors_external_soil.csv "
        "(download: https://raw.githubusercontent.com/pyne/pyne/develop/pyne/dbgen/dosefactors_external_soil.csv)",
    )
    ap.add_argument(
        "--dose-ingest",
        required=False,
        default=None,
        help="Local copy of PyNE dosefactors_ingest.csv "
        "(download: https://raw.githubusercontent.com/pyne/pyne/develop/pyne/dbgen/dosefactors_ingest.csv)",
    )
    ap.add_argument(
        "--dose-inhale",
        required=False,
        default=None,
        help="Local copy of PyNE dosefactors_inhale.csv "
        "(download: https://raw.githubusercontent.com/pyne/pyne/develop/pyne/dbgen/dosefactors_inhale.csv)",
    )
    ap.add_argument("--out", required=True, help="crates/nuclei/src/data dir")
    args = ap.parse_args()

    legacy = [args.endf_neutrons, args.endf_decay, args.nist_html]
    dose_args = [args.dose_air, args.dose_soil, args.dose_ingest, args.dose_inhale]
    if all(v is None for v in legacy) and all(v is None for v in dose_args):
        print("nothing to do: pass ENDF/NIST flags and/or --dose-* flags")
        return 2
    if any(v is None for v in legacy) and not all(v is None for v in legacy):
        print("legacy tables need --endf-neutrons, --endf-decay, and --nist-html together")
        return 2
    if any(v is None for v in dose_args) and not all(v is None for v in dose_args):
        print("dose table needs --dose-air, --dose-soil, --dose-ingest, --dose-inhale together")
        return 2

    if all(v is None for v in dose_args):
        dose_rows = None
    else:
        assert args.dose_air and args.dose_soil and args.dose_ingest and args.dose_inhale
        dose_rows, dose_log = gen_dose(
            args.dose_air, args.dose_soil, args.dose_ingest, args.dose_inhale
        )
        for name, pathway, source, expected in DOSE_SPOTS:
            got = dose_rows.get((name, pathway, source))
            if got is None:
                print(f"SPOT FAIL: {(name, pathway, source)} missing from dose")
                return 1
            if abs(got[0] - expected) / max(abs(expected), 1e-30) > 1e-6:
                print(f"SPOT FAIL: {(name, pathway, source)}={got[0]:.6g} != {expected:.6g}")
                return 1
            print(f"spot ok: dose       {name:10s} {pathway:7s} {source:5s} {got[0]:.6g}")
        print(f"rows: dose={len(dose_rows)}")
        for line in dose_log:
            print(f"note: {line}")
        write_dose_tsv(
            os.path.join(args.out, "dose_factors.tsv"),
            [
                "GNDS name\tpathway\tsource\tfactor\tf1\tlung_model",
                "GNDS name\tpathway(air|soil|ingest|inhale)\tsource(EPA|DOE|GENII)",
                "factor\tf1(to body fluids, ingest only)\tlung(D|W|Y|V|O, inhale only)",
                "Dose factors from the BSD-3 PyNE dbgen tables (no separate CSV license):",
                "dosefactors_external_air.csv (93 rows; Nuclide,Air Dose Rate Factor,Ratio;",
                "mrem/h per Ci/m^3; EPA-only), dosefactors_external_soil.csv (93 rows;",
                "Nuclide,GENII,EPA,DOE,...; mrem/h per Ci/m^2, 15cm slab),",
                "dosefactors_ingest.csv (93 rows; Nuclide,f1,GENII,EPA,DOE,...; mrem/pCi),",
                "dosefactors_inhale.csv (93 rows; Nuclide,Lung Model,GENII,EPA,DOE,...;",
                "mrem/pCi; lung D/W/Y/V/O).",
                "Proximate provenance: HNF-SD-WM-TI-707 Rev.1 (Dec 1999), App. O of",
                "HNF-5636 (2001); GENII/EPA/DOE are 3 parallel evaluations in that report.",
                "Semantics follow pyne/dbgen/dosefactors.py (260 lines): Nuclide keys like",
                "H-3 with +D = plus daughters fold into the parent (suffix stripped).",
                "Air is EPA-only: GENII/DOE air rows are -1 sentinels (PyNE convention);",
                "accessors treat negative factors as missing. Liability: not for safety",
                "decisions (upstream PyNE disclaimer).",
                "Regenerate: python3 scripts/gen-nuclear-data.py --help.",
            ],
            dose_rows,
        )
        if all(v is None for v in legacy):
            return 0

    assert args.endf_neutrons and args.endf_decay and args.nist_html
    xs_nist, xs_extra, elem100, sl_log = gen_scattering(args.nist_html)
    # Attribute ~100%-abundance element rows to the matching tape isotope,
    # but ONLY for provably monoisotopic elements (single mass number in the
    # vendored natural-abundance table).  This stops e.g. the Cs element row
    # from leaking onto artificial Cs134-137 while letting Be9/Na23/Al27 in.
    repo_data = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "crates", "nuclei", "src", "data"
    )
    mono = monoisotopic_a(os.path.join(repo_data, "natural_abundance.tsv"))
    tape_isotopes: set[str] = set()
    for fname in sorted(os.listdir(args.endf_neutrons)):
        ident = tape_id(fname)
        if ident is not None and fname.startswith("n-"):
            _, sym, a, state = ident
            if a != 0 and state == 0:
                tape_isotopes.add(gnds(sym, a, 0))
    for name in sorted(tape_isotopes - set(xs_nist)):
        sym = re.match(r"[A-Za-z]+", name).group(0)  # type: ignore[union-attr]
        want_a = int(re.match(r"[A-Za-z]+(\d+)", name).group(1))  # type: ignore[union-attr]
        if sym in elem100 and mono.get(sym) == want_a:
            bundle = elem100[sym]
            xs_nist[name] = (bundle["b_coh"], bundle["b_inc"])  # type: ignore[misc]
            xs_extra[name] = bundle
            sl_log.append(f"monoisotopic element row attributed: {sym} -> {name}")
    xs, xs_log = gen_simple_xs(args.endf_neutrons, xs_extra)
    de, _half, de_log = gen_decay(args.endf_decay)
    sl = xs_nist

    tables = {"simple_xs": xs, "decay": de, "scattering": sl}
    for kind, name, lo, hi, _unit in SPOTS:
        val = tables[kind].get(name)
        if val is None:
            print(f"SPOT FAIL: {name} missing from {kind}")
            return 1
        probe = val[0] if isinstance(val, tuple) else val
        if not lo <= probe <= hi:
            print(f"SPOT FAIL: {name}={probe} outside [{lo}, {hi}]")
            return 1
        print(f"spot ok: {kind:10s} {name:10s} {probe:.6g}")

    print(f"rows: simple_xs={len(xs)} decay={len(de)} scattering={len(sl)}")
    for line in xs_log + de_log + sl_log:
        print(f"note: {line}")

    write_tsv(
        os.path.join(args.out, "simple_xs.tsv"),
        [
            "GNDS name\tthermal_barn\tfast14mev_barn",
            "Thermal totals: NIST NCNR bound XS converted to free-atom via",
            "xs*(A/(A+1))^2 + xs_a. Fast totals: ENDF/B-VII.1 MF3/MT1 at 14 MeV.",
            "Screening-level only: use evaluated libraries for transport.",
            "Regenerate: python3 scripts/gen-nuclear-data.py --help.",
        ],
        xs,
        lambda v: f"{v[0]:.6g}\t{v[1]:.6g}",
    )
    write_tsv(
        os.path.join(args.out, "decay_energy.tsv"),
        [
            "GNDS name\tmev_per_decay",
            "Mean PROMPT recoverable energy per decay (MeV) from ENDF/B-VII.1",
            "decay tapes (MF8/MT457 summary components, uncertainties skipped).",
            "Daughter gammas belong to the daughter row: chain codes must sum",
            "members (e.g. Cs137 prompt + Ba137_m1 662 keV). Mean-field decay",
            "data evaluation; not for spectroscopy or safety calculations.",
            "Regenerate: python3 scripts/gen-nuclear-data.py --help.",
        ],
        de,
        lambda v: f"{v:.6g}",
    )
    write_tsv(
        os.path.join(args.out, "scattering_lengths.tsv"),
        [
            "GNDS name\tb_coherent_fm\tb_incoherent_fm",
            "Bound lengths (fm) from the NIST NCNR tabulation (Sears 1992).",
            "Complex coherent lengths kept by real part; incoherent derived",
            "via b = 10*sqrt(sigma_i/4*pi) where the table gives only sigma.",
            "Regenerate: python3 scripts/gen-nuclear-data.py --help.",
        ],
        sl,
        lambda v: f"{v[0]:.6g}\t{v[1]:.6g}",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
