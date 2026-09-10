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
- ``decay_branches.tsv``: per-branch daughters from the ENDF/B-VIII.0
  decay sublibrary (MF8/MT457 NDK records: RTYP decay-mode code, RFS
  daughter state flag, BR branching fraction).  One row per kept branch:
  ``parent_GNDS``, ``progeny_GNDS``, ``bf``, ``mode``.  Spontaneous-fission
  and fission-family branches are dropped (depletion matrices skip ``sf``
  gains); tapes with zero half-life or stable flags yield no rows, so
  effectively-stable entries (e.g. Te123) stay absent, matching the
  half-life table's stable-absent convention.
  Upstream: https://www.nndc.bnl.gov/endf-b8.0/
  (``zips/ENDF-B-VIII.0_decay.zip``).
- ``ame2020.tsv`` isomer rows: for every ENDF/B-VIII.0 isomer tape
  (``*m1``/``*m2``) with a positive File-1 MT451 ``ELIS`` excitation
  energy, one row keyed by the full state-bearing nucid:
  ``m = m_ground(AME2020) + ELIS[eV]/1e6/931.49410242``.  Ground rows are
  preserved verbatim; no NUBASE import (ENDF excitation energies only).
- ``half_life.tsv``: ``T1/2`` in seconds for every ENDF/B-VIII.0 decay tape
  with a usable MF8/MT457 NDK half-life (stable tapes flagged ``NST != 0``
  and zero-half-life evaluation dummies yield no rows — the same
  stable-absent convention as the branch table).  ``T1/2`` is the first
  value of the MT457 summary LIST record (ENDF-102 section 8.4.1:
  ``[MAT, 8,457/ T1/2, dT1/2, 0, 0, 2*NC, 0 / ...]``; units are seconds
  per the ENDF-102 quantity dictionary).  Values are emitted with
  ``repr`` (shortest round-trip) so the writer reproduces the committed
  table byte-identically (data rows).
  Upstream: https://www.nndc.bnl.gov/endf-b8.0/
  (``zips/ENDF-B-VIII.0_decay.zip``).

Per-table evaluation bases (kept distinct on purpose): ``decay_energy.tsv``
stays on ENDF/B-VII.1 while ``decay_branches.tsv``, ``half_life.tsv``, and
the ``ame2020.tsv`` isomer rows use ENDF/B-VIII.0.  Each TSV header records
its own basis.

Generator path for ``ame2020.tsv`` (so the ground rows are never orphaned):
ground-state rows are condensed once from the AME2020 mass table
(``mass.mas20``; Huang et al., Chinese Physics C 45, 030002/030003, 2021)
and carried verbatim by ``read_ame_grounds`` on every ``--endf-decay8``
run; ``gen_isomer_masses`` only appends isomer rows
(``m_ground + ELIS/931.49410242``) and never rewrites grounds.  There is
no AME2020 re-parse path in this script: to refresh grounds, replace the
verbatim block in ``ame2020.tsv`` from a new ``mass.mas20`` and re-run
``--endf-decay8`` to re-append isomers.

Usage::

    python3 scripts/gen-nuclear-data.py \
        --endf-neutrons /path/to/endf-b-vii.1/neutrons \
        --endf-decay /path/to/endf-b-vii.1/decay \
        --nist-html /path/to/scattering_lengths.html \
        --endf-decay8 /path/to/endf-b-viii.0/decay \
        --dose-air /tmp/dosefactors_external_air.csv \
        --dose-soil /tmp/dosefactors_external_soil.csv \
        --dose-ingest /tmp/dosefactors_ingest.csv \
        --dose-inhale /tmp/dosefactors_inhale.csv \
        --out crates/nuclei/src/data

All inputs are local checkouts (see help for download URLs); nothing is
fetched over the network.  Exit nonzero if any hard spot-check fails.
Pass only the ``--dose-*`` flags (plus ``--out``) for a dose-only run;
pass only the ENDF/NIST flags for the legacy tables; pass only
``--endf-decay8`` (plus ``--out``) for the VIII.0 branch/half-life/isomer tables.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys

THERMAL_EV = 0.0253  # NIST 2200 m/s thermal energy (NCNR n-lengths page).
FAST_EV = 14.0e6  # ENDF 14-MeV reference energy for fast totals (MF3/MT1).

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


# ---------------------------------------------------------------------------
# ENDF/B-VIII.0 decay branches + isomer masses.
#
# RTYP provenance: exact float codes confirmed against ENDF-102 section 8.4
# (2023 manual, BNL-224854-2023-INRE: the decay-mode table lists 0 gamma —
# not used in MT457 — 1 beta-, 2 EC/beta+, 3 IT, 4 alpha, 5 neutron, 6 SF,
# 7 proton, 10 unknown, plus the multi-particle construction rule where the
# digits apply in emission order, e.g. RTYP = 1.5 is beta- decay followed by
# neutron emission) and the MIT-licensed OpenMC implementation
# (openmc/data/decay.py ``_DECAY_MODES`` + ``get_decay_modes``: digit-wise
# parse of the RTYP float with 10.0 special-cased to unknown).  Mode-name
# strings below reuse OpenMC's tokens verbatim so depletion matrices keep
# their ``sf``-skip / ``alpha`` / ``p`` substring semantics.
# ---------------------------------------------------------------------------

#: MeV per u, matching ``MEV_PER_U`` in ``crates/nuclei/src/data.rs``.
MEV_PER_U = 931.494_102_42

# Single-digit RTYP codes: (mode token, dA, dZ). Digits apply in emission
# order to (Z, A); the NDK RFS flag sets the daughter's isomer state.
# Digits 0/8/9 (gamma / conversion-electron / x-ray) carry no nucleons and
# are no-ops for progeny. Digit 6 (SF) and code 10 (unknown origin) have no
# single-daughter mapping: those branches are dropped.
_RTYP_DIGITS = {
    1: ("beta-", 0, 1),
    2: ("ec/beta+", 0, -1),
    3: ("IT", 0, 0),
    4: ("alpha", -4, -2),
    5: ("n", -1, 0),
    7: ("p", -1, -1),
}


def rtyp_digits(value: float) -> list[int] | None:
    """Split an NDK RTYP float into emission digits (OpenMC ``get_decay_modes``).

    Returns None for unknown origin (10.0). Multi-particle codes expand
    digit-wise: 1.5 is beta- + neutron, 2.77 is EC + proton + proton, 1.1 is
    double-beta decay (two successive beta- steps).
    """
    if int(value) == 10:
        return None
    return [int(x) for x in str(value).strip("0").replace(".", "")]


def mt451_elis(lines: list[str]) -> tuple[float, int, int] | None:
    """Return (ELIS_eV, LIS, LISO) from a decay tape's File 1 MT451 record."""
    sec = section_lines(lines, 1, 451)
    if len(sec) < 2:
        return None
    try:
        rec = fields(sec[1])
        return rec[0], int(rec[2]), int(rec[3])
    except (IndexError, ValueError):
        return None


def ndk_modes(lines: list[str]) -> tuple[float, list[tuple[float, float, float]]] | None:
    """Return (half_life_s, [(RTYP, RFS, BR)]) from MF8/MT457, or None.

    Returns None for stable tapes (NST flag set), zero half-lives
    (evaluation dummies such as Te123/Ca46), missing sections, or tapes
    with no decay modes. Q values and uncertainties are not needed here.
    """
    sec = section_lines(lines, 8, 457)
    if len(sec) < 4:
        return None
    try:
        head = fields(sec[0])
        nst = int(head[4])
        if nst != 0:
            return None
        hl_head = fields(sec[1])
        t12 = hl_head[0]
        if not t12 > 0.0:
            return None
        n_avg = (int(hl_head[4]) + 5) // 6
        spin = fields(sec[2 + n_avg])
        ndk = int(spin[5])
        if ndk <= 0:
            return None
        vals: list[float] = []
        for body in sec[3 + n_avg :]:
            vals.extend(fields(body))
            if len(vals) >= 6 * ndk:
                break
        if len(vals) < 6 * ndk:
            return None
    except (IndexError, ValueError):
        return None
    modes = [(vals[6 * i], vals[6 * i + 1], vals[6 * i + 4]) for i in range(ndk)]
    return t12, modes


def gen_decay_branches(
    decay8_dir: str,
) -> tuple[list[tuple[str, str, float, str]], dict[str, float], list[str]]:
    """Build (parent, progeny, bf, mode) rows from ENDF/B-VIII.0 decay tapes.

    Progeny apply the RTYP digits in emission order (beta- ``Z+1``,
    EC/beta+ ``Z-1``, alpha ``Z-2/A-4``, IT unchanged, delayed neutrons and
    protons subtract the emitted nucleons); the mode token is the initial
    event per ENDF-102 procedure 8.4.2-3. SF/fission-family branches
    (any digit 6) and unknown-origin branches (10.0) are dropped. Also
    returns {parent: half_life_s} for spot-checks.
    """
    rows: list[tuple[str, str, float, str]] = []
    half: dict[str, float] = {}
    log: list[str] = []
    n_sf = n_unknown = n_zero_br = n_bad = n_tapes = 0
    kept_sum: dict[str, float] = {}
    for fname in sorted(os.listdir(decay8_dir)):
        ident = tape_id(fname)
        if ident is None or not fname.startswith("dec-"):
            continue
        z, sym, a, state = ident
        lines = read_lines(os.path.join(decay8_dir, fname))
        parsed = ndk_modes(lines)
        if parsed is None:
            continue
        t12, modes = parsed
        parent = gnds(sym, a, state)
        half[parent] = t12
        n_tapes += 1
        for rtyp, rfs, br in modes:
            digits = rtyp_digits(rtyp)
            if digits is None:
                n_unknown += 1
                log.append(f"drop {fname}: unknown-origin branch (RTYP=10)")
                continue
            if 6 in digits:
                n_sf += 1
                continue
            if any(d not in _RTYP_DIGITS and d not in (0, 8, 9) for d in digits):
                n_bad += 1
                log.append(f"drop {fname}: unmapped RTYP digits in {rtyp:.6g}")
                continue
            if br == 0.0:
                n_zero_br += 1
                continue
            first = next(d for d in digits if d in _RTYP_DIGITS)
            mode = _RTYP_DIGITS[first][0]
            zz, aa = z, a
            for d in digits:
                if d in _RTYP_DIGITS:
                    _, da, dz = _RTYP_DIGITS[d]
                    aa += da
                    zz += dz
            rfs_i = int(rfs)
            if not 1 <= zz <= 118 or not zz <= aa <= 999 or rfs_i < 0 or rfs_i > 9:
                n_bad += 1
                log.append(f"drop {fname}: invalid progeny Z={zz} A={aa} RFS={rfs_i}")
                continue
            progeny = gnds(SYMBOLS[zz - 1], aa, rfs_i)
            rows.append((parent, progeny, br, mode))
            kept_sum[parent] = kept_sum.get(parent, 0.0) + br
    for parent, total in sorted(kept_sum.items()):
        if abs(total - 1.0) > 0.01:
            log.append(f"note {parent}: kept branches sum to {total:.6g} (SF dropped)")
    log.append(
        f"branch tapes={n_tapes} rows={len(rows)} "
        f"dropped_sf={n_sf} dropped_unknown={n_unknown} "
        f"dropped_zero_br={n_zero_br} dropped_bad={n_bad}"
    )
    return rows, half, log


def gen_isomer_masses(
    decay8_dir: str, ground: dict[int, tuple[float, str]]
) -> tuple[dict[int, tuple[float, str]], list[str]]:
    """Build {full_nucid: (mass_u, unc_str)} for ENDF/B-VIII.0 isomer tapes.

    ``mass = m_ground(AME2020) + ELIS[eV]/1e6/MEV_PER_U`` with ELIS from
    the tape's File 1 MT451 record (LIS/LISO identify the level). ENDF
    excitation energies only; no NUBASE import. Tapes with unset (zero)
    ELIS still get a row carrying the ground-state mass, flagged in the log.
    """
    rows: dict[int, tuple[float, str]] = {}
    log: list[str] = []
    n_no_elis = n_no_ground = n_zero_elis = 0
    for fname in sorted(os.listdir(decay8_dir)):
        ident = tape_id(fname)
        if ident is None or not fname.startswith("dec-"):
            continue
        z, sym, a, state = ident
        if state == 0:
            continue
        lines = read_lines(os.path.join(decay8_dir, fname))
        info = mt451_elis(lines)
        if info is None:
            n_no_elis += 1
            log.append(f"skip {fname}: no MT451 ELIS record")
            continue
        elis_ev, _lis, liso = info
        if liso != state:
            log.append(f"note {fname}: filename m{state} vs MT451 LISO={liso}")
        if elis_ev <= 0.0:
            n_zero_elis += 1
            log.append(f"note {fname}: unset ELIS carries the ground-state mass")
        key = (z * 1000 + a) * 10_000
        entry = ground.get(key)
        if entry is None:
            n_no_ground += 1
            log.append(f"skip {fname}: no AME2020 ground row for {gnds(sym, a, 0)}")
            continue
        g_mass, g_unc = entry
        rows[key + state] = (g_mass + elis_ev / 1e6 / MEV_PER_U, g_unc)
    log.append(
        f"isomer rows={len(rows)} missing_elis={n_no_elis} "
        f"missing_ground={n_no_ground} zero_elis={n_zero_elis}"
    )
    return rows, log


def read_ame_grounds(path: str) -> tuple[list[str], dict[int, tuple[float, str]]]:
    """Return (verbatim ground lines, {ground_nucid: (mass, unc_str)}).

    State-bearing rows (if the file already carries isomer rows from a
    previous run) are excluded from ``ground`` but preserved by the caller.
    """
    lines: list[str] = []
    ground: dict[int, tuple[float, str]] = {}
    with open(path) as fh:
        for line in fh.read().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            nucid = int(cols[0])
            if nucid % 10 == 0:
                lines.append(line)
                ground[nucid] = (float(cols[1]), cols[2] if len(cols) > 2 else "")
    return lines, ground


def ame_isomer_lines(ground_lines: list[str], isomers: dict[int, tuple[float, str]]) -> list[str]:
    """Merge verbatim ground lines with isomer rows, sorted by nucid."""
    merged: dict[int, str] = {}
    for line in ground_lines:
        merged[int(line.split("\t")[0])] = line
    for nucid, (mass, unc) in isomers.items():
        merged[nucid] = f"{nucid}\t{mass:.9f}\t{unc}"
    return [merged[k] for k in sorted(merged)]


# ENDF/B-VIII.0 decay-sublibrary source pin for the branch/isomer/half-life
# tables.  Download URL + MD5 of the zip as served (verified 2026-09-10):
# the writer below reproduces the committed half_life.tsv data rows
# byte-identically from an unpack of this artifact (3,821 decay tapes,
# 3,561 with usable MF8/MT457 NDK half-lives).
DECAY8_ZIP_URL = "https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_decay.zip"
DECAY8_ZIP_MD5 = "aa80cd0a880d9d7e0905940b868370c3"

# Spot-checks on the VIII.0 branch/isomer tables. Failure aborts the run.
# All values are read off the ENDF/B-VIII.0 decay tapes themselves (see the
# per-branch notes); K-40 additionally pins the half-life table's K40 entry.
BRANCH_SPOTS = [
    # (parent, progeny, mode, expected_bf)
    ("K40", "Ca40", "beta-", 0.8914),
    ("K40", "Ar40", "ec/beta+", 0.1086),
    ("Es254", "Bk250", "alpha", 1.0),
    ("Ba137_m1", "Ba137", "IT", 1.0),
    ("He8", "Li8", "beta-", 0.84),
    ("He8", "Li7", "beta-", 0.16),
    ("Es254_m1", "Fm254", "beta-", 0.98),
    ("Es254_m1", "Es254", "IT", 0.0155),
]

# (parent, expected half-life in s, relative tolerance). K-40 pins the
# half-life table entry derived from the same VIII.0 checkout.
BRANCH_HALF_LIFE_SPOTS = [
    ("K40", 3.93839e16, 1e-6),
    ("Es254", 2.382048e7, 1e-6),
    ("Es254_m1", 141479.9, 1e-6),
]

# Parents that must NOT appear (zero-half-life evaluation dummies and the
# effectively-stable Te123, matching the half-life stable-absent rule).
BRANCH_ABSENT = ["Te123", "Ca46"]

# Isomer-mass spots: (GNDS isomer, ELIS in eV from its MT451 record).
ISOMER_ELIS_SPOTS_EV = [
    ("Ba137_m1", 661659.0),
    ("Te123_m1", 247470.0),
    ("Es254_m1", 80000.0),
]


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


def write_branch_tsv(path: str, header: list[str], rows: list[tuple[str, str, float, str]]) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join("# " + h for h in header) + "\n")
        for parent, progeny, bf, mode in sorted(rows):
            fh.write(f"{parent}\t{progeny}\t{bf:.6g}\t{mode}\n")


def run_decay8(out_dir: str, decay8_dir: str) -> int:
    """Generate ``decay_branches.tsv`` + ``half_life.tsv`` + extended ``ame2020.tsv``.

    Nonzero on spot failure.  The half-life writer emits ``repr`` values
    straight from ``ndk_modes`` (same stable-absent rule as the branch
    table), which reproduces the committed data rows byte-identically.
    """
    branches, branch_half, br_log = gen_decay_branches(decay8_dir)

    by_parent: dict[str, list[tuple[str, str, float, str]]] = {}
    for row in branches:
        by_parent.setdefault(row[0], []).append(row)

    for parent, progeny, mode, expected in BRANCH_SPOTS:
        cands = [r for r in by_parent.get(parent, []) if r[1] == progeny and r[3] == mode]
        if not cands:
            print(f"SPOT FAIL: {parent}->{progeny} [{mode}] missing from branches")
            return 1
        got = cands[0][2]
        if abs(got - expected) / max(abs(expected), 1e-30) > 1e-6:
            print(f"SPOT FAIL: {parent}->{progeny} [{mode}]={got:.6g} != {expected:.6g}")
            return 1
        print(f"spot ok: branch     {parent:10s} -> {progeny:10s} [{mode:8s}] {got:.6g}")

    for parent, expected, tol in BRANCH_HALF_LIFE_SPOTS:
        got = branch_half.get(parent)
        if got is None:
            print(f"SPOT FAIL: {parent} missing half-life in VIII.0 tapes")
            return 1
        if abs(got - expected) / expected > tol:
            print(f"SPOT FAIL: {parent} T1/2={got:.6g} != {expected:.6g}")
            return 1
        print(f"spot ok: halflife   {parent:10s} {got:.6g} s")

    for parent in BRANCH_ABSENT:
        if parent in by_parent:
            print(f"SPOT FAIL: {parent} must stay absent from branches (stable-absent rule)")
            return 1
        print(f"spot ok: absent     {parent:10s} (no branch rows)")

    # K-40 agreement: the two evaluated branches sum to unity.
    k_sum = sum(r[2] for r in by_parent.get("K40", []))
    if abs(k_sum - 1.0) > 1e-9:
        print(f"SPOT FAIL: K40 branches sum to {k_sum:.6g}, expected 1.0")
        return 1
    print(f"spot ok: agreement  K40 kept branches sum to {k_sum:.6g}")

    # Isomer masses from MT451 ELIS + vendored AME2020 grounds.
    ame_path = os.path.join(out_dir, "ame2020.tsv")
    ground_lines, ground = read_ame_grounds(ame_path)
    isomers, iso_log = gen_isomer_masses(decay8_dir, ground)
    for name, elis_ev in ISOMER_ELIS_SPOTS_EV:
        m = re.match(r"^([A-Za-z]+)(\d+)_m(\d+)$", name)
        assert m is not None
        key = (Z_OF[m.group(1)] * 1000 + int(m.group(2))) * 10_000 + int(m.group(3))
        got_mass, _unc = isomers.get(key, (None, ""))
        if got_mass is None:
            print(f"SPOT FAIL: {name} missing from isomer masses")
            return 1
        g_mass = ground[key - int(m.group(3))][0]
        expected_mass = g_mass + elis_ev / 1e6 / MEV_PER_U
        if abs(got_mass - expected_mass) > 1e-9:
            print(f"SPOT FAIL: {name} mass {got_mass:.9f} != {expected_mass:.9f}")
            return 1
        print(f"spot ok: isomer     {name:10s} {got_mass:.9f} u (ELIS {elis_ev:.6g} eV)")

    print(
        f"rows: branches={len(branches)} parents={len(by_parent)}"
        f" isomers={len(isomers)} halflives={len(branch_half)}"
    )
    for line in br_log + iso_log:
        print(f"note: {line}")

    pm_tape = "dec-061_Pm_137m1.endf"
    pm_present = os.path.exists(os.path.join(decay8_dir, pm_tape))
    pm_note = (
        f"{pm_tape} present: Pm137_m1 emitted from its ENDF tape."
        if pm_present
        else f"{pm_tape} absent from ENDF/B-VIII.0: no Pm137_m1 row (no fill-in)."
    )
    print(f"note: {pm_note}")

    write_branch_tsv(
        os.path.join(out_dir, "decay_branches.tsv"),
        [
            "parent_GNDS\tprogeny_GNDS\tbf\tmode",
            "Per-branch daughters from ENDF/B-VIII.0 decay tapes (MF8/MT457 NDK",
            "records: RTYP decay-mode code, RFS daughter state flag, BR branching",
            "fraction). RTYP digits apply in emission order (ENDF-102 8.4; OpenMC",
            "decay.py digit table): beta- Z+1, EC/beta+ Z-1, alpha Z-2/A-4, IT",
            "unchanged, delayed neutrons/protons subtract emitted nucleons; the",
            "mode token is the initial event. SF/fission-family branches (any",
            "digit 6) and unknown-origin branches (RTYP 10) are dropped, so kept",
            "branches of strong SF emitters sum to 1-BR(SF); depletion matrices",
            "skip sf gains. Tapes with zero half-life or stable flags yield no",
            "rows (Te123 stays absent: effectively stable in this evaluation).",
            "Basis note: this table and half_life.tsv use ENDF/B-VIII.0 while",
            "decay_energy.tsv stays on ENDF/B-VII.1 (different sublibrary vintages).",
            "Screening-level only: use evaluated libraries for transport.",
            "Regenerate: python3 scripts/gen-nuclear-data.py --endf-decay8 <dir> --out <dir>.",
        ],
        branches,
    )

    half_body = [f"{name}\t{branch_half[name]}" for name in sorted(branch_half)]
    with open(os.path.join(out_dir, "half_life.tsv"), "w") as fh:
        fh.write(
            "\n".join(
                [
                    "# GNDS name\thalf_life_seconds",
                    "# Half-lives (T1/2) in seconds from the ENDF/B-VIII.0 decay",
                    "# sublibrary (MF8/MT457): T1/2 is the first value of the",
                    "# summary LIST record (ENDF-102 section 8.4.1: [MAT, 8,457/",
                    "# T1/2, dT1/2, 0, 0, 2*NC, 0 / ...]; seconds per the",
                    "# ENDF-102 quantity dictionary). Stable tapes (NST != 0)",
                    "# and zero-half-life evaluation dummies (e.g. Te123, Ca46)",
                    "# yield no rows: stable nuclides are simply absent.",
                    "# Basis note: this table and decay_branches.tsv use",
                    "# ENDF/B-VIII.0 while decay_energy.tsv stays on ENDF/B-VII.1.",
                    f"# Source: {DECAY8_ZIP_URL}",
                    f"# Source zip MD5 (as served): {DECAY8_ZIP_MD5}",
                    "# Screening-level only: use evaluated libraries for",
                    "# transport; not for safety calculations (no warranty).",
                    "# Regenerate: python3 scripts/gen-nuclear-data.py --endf-decay8"
                    " <dir> --out <dir>.",
                ]
            )
            + "\n"
        )
        for line in half_body:
            fh.write(line + "\n")

    isomer_body = ame_isomer_lines(ground_lines, isomers)
    with open(ame_path, "w") as fh:
        fh.write(
            "\n".join(
                [
                    "# nucid\tmass_u\tuncertainty_u",
                    "# Ground states: AME2020 (Huang et al., Chinese Physics C 45,",
                    "# 030002/030003, 2021), condensed from the mass.mas20 table",
                    "# and carried verbatim by read_ame_grounds (this script has",
                    "# no AME2020 re-parse path: refresh grounds from mass.mas20,",
                    "# then re-run --endf-decay8 to re-append isomers).",
                    "# Isomer rows (state > 0, full nucid):",
                    "# m = m_ground(AME2020) + E*/931.49410242 with E* (eV) = ELIS",
                    "# from the ENDF/B-VIII.0 decay tape's File 1 MT451 record",
                    "# (LIS/LISO identify the level). ENDF excitation energies only;",
                    "# no NUBASE import. Isomer tapes with unset (zero) ELIS carry",
                    "# the ground-state mass. Uncertainty column carries the ground-state",
                    "# AME2020 uncertainty.",
                    f"# {pm_note}",
                    "# Regenerate: python3 scripts/gen-nuclear-data.py --endf-decay8"
                    " <dir> --out <dir>.",
                ]
            )
            + "\n"
        )
        for line in isomer_body:
            fh.write(line + "\n")
    return 0


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
        "--endf-decay8",
        required=False,
        default=None,
        help="ENDF/B-VIII.0 decay tapes dir for the branch/half-life/isomer tables "
        "(download: https://www.nndc.bnl.gov/endf-b8.0/ "
        "zips/ENDF-B-VIII.0_decay.zip). Ground rows of ame2020.tsv are "
        "carried verbatim (condensed from AME2020 mass.mas20, no re-parse "
        "path); only isomer rows are re-appended.",
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
    if (
        all(v is None for v in legacy)
        and all(v is None for v in dose_args)
        and args.endf_decay8 is None
    ):
        print("nothing to do: pass ENDF/NIST flags, --endf-decay8, and/or --dose-* flags")
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
                "Regenerate: python3 scripts/gen-nuclear-data.py --dose-air <air.csv>"
                " --dose-soil <soil.csv> --dose-ingest <ingest.csv>"
                " --dose-inhale <inhale.csv> --out <dir>.",
            ],
            dose_rows,
        )
        if all(v is None for v in legacy) and args.endf_decay8 is None:
            return 0

    if args.endf_decay8 is not None:
        rc = run_decay8(args.out, args.endf_decay8)
        if rc != 0:
            return rc
        if all(v is None for v in legacy) and all(v is None for v in dose_args):
            return 0

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
    # Legacy VII.1 half-life dict intentionally unused: the canonical
    # half-life table is the VIII.0 writer in run_decay8 (repr values from
    # ndk_modes, same stable-absent rule as the branch table).
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
            "Regenerate: python3 scripts/gen-nuclear-data.py --endf-neutrons <vii1-neutrons-dir>"
            " --endf-decay <vii1-decay-dir> --nist-html <scattering_lengths.html> --out <dir>.",
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
            "Regenerate: python3 scripts/gen-nuclear-data.py --endf-neutrons <vii1-neutrons-dir>"
            " --endf-decay <vii1-decay-dir> --nist-html <scattering_lengths.html> --out <dir>.",
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
            "Regenerate: python3 scripts/gen-nuclear-data.py --endf-neutrons <vii1-neutrons-dir>"
            " --endf-decay <vii1-decay-dir> --nist-html <scattering_lengths.html> --out <dir>.",
        ],
        sl,
        lambda v: f"{v[0]:.6g}\t{v[1]:.6g}",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
