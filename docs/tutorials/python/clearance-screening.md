---
title: Clearance screening
sidebar:
  order: 21
---

Clearance screening compares a radioactive inventory against per-nuclide
clearance levels. The clearance index `CI = sum_i A_i / CL_i` adds up each
nuclide's activity divided by its clearance level, and the sum-of-fractions
rule treats `CI <= 1` as satisfying the screening criterion (the boundary
included). Nucleide's `nucleide.alara` helpers run that arithmetic over plain
`{nuclide: activity}` dicts, and `nucleide.fispact` reads the clearance-bearing
inventory table from FISPACT-II output. This tutorial shows both; the
implementations live in `crates/alara-io` and `crates/fispact-io`.

## Reading the FISPACT-II clearance block

FISPACT-II prints a wide inventory table with a per-nuclide clearance column
when the input deck sets the `HAZARDS` and `CLEAR` keywords.
`fispact_parse_clearance` reads that table from the output text — one dict per
(nuclide, time step) row:

```python
from pathlib import Path

from nucleide.fispact import fispact_parse_clearance

text = Path("fixtures/fispact/clearance.out").read_text()
rows = fispact_parse_clearance(text)
print(len(rows))

step1 = [r for r in rows if r["interval"] == 1]
for r in step1:
    print(r["nuclide"], r["flags"], r["activity_bq"], r["clearance_index"], r["half_life_s"])
```

```text
7
V-55 > 0.0 0.0 -1.0
Fe-56 #> 0.0 0.0 -1.0
Co-60 > 4000000.0 1000000.0 166360000.0
Rb-86m & 8000000000.0 4000000000.0 61.0
```

Each row carries the 1-based `interval`, the step time in seconds (`time_s`)
with its verbatim `time_label`, a `cooling` flag (true for `COOLING TIME IS`
steps), the nuclide name in the shared dialect spelling (`Co-60`), the
verbatim flag characters, the `activity_bq` column, the half-life in seconds,
and the code's own `clearance_index` column — FISPACT-II's per-nuclide
`A_i / CL_i` against its internal clearance library, not against a caller
table. Stable nuclides are stored with `half_life_s == -1.0` (V-55 and Fe-56
above). The flags keep FISPACT-II's markers verbatim: `#` stable, `>` present
before irradiation, `&` gamma spectrum approximately calculated, `?`
convergence not reached.

## Screening against the default clearance table

`alara_clearance_index` computes the clearance index over an inventory dict.
With `limits=None` (the default) the comparison uses Nucleide's vendored
transcription of EU 2013/59/Euratom Annex VII Table A — 260
activity-concentration limits in Bq/g for solid material. Because that table
lists concentrations, the inventory must be massic activities (Bq/g), not
total activities (Bq):

```python
from nucleide.alara import (
    alara_clearance_eu_table,
    alara_clearance_index,
    alara_sum_of_fractions,
)

eu = alara_clearance_eu_table()
print(len(eu), eu["Co-60"], eu["H-3"])

inventory = {"Co-60": 0.05, "Cs-137": 0.2}  # Bq/g
print(alara_clearance_index(inventory))

out = alara_sum_of_fractions(inventory)
print(out["sum"], out["class"], out["max_nuclide"], out["max_fraction"])
```

```text
260 0.1 100.0
2.5
2.5 exceeded Cs-137 2.0
```

`alara_sum_of_fractions` returns the same sum plus the screening `class`
(`"satisfied"` when `sum <= 1`, otherwise `"exceeded"`) and the dominant
contributor — the nuclide with the largest single fraction `A_i / CL_i` (here
Cs-137 at 0.2 Bq/g against a 0.1 Bq/g limit). Nuclide keys accept any
shared-dialect spelling (`"Co60"`, `"Co-60"`, and `"co-60"` resolve to the
same nuclide); results report the canonical spelling.

## Spanish conditional tables

For Spanish conditional NORM landfill screening, `alara_clearance_es_table`
selects one of the three CSN Tables 1–3 explicitly by landfill type; the
caller picks the governing table and there is no cross-table logic:

```python
from nucleide.alara import alara_clearance_es_table, alara_sum_of_fractions

es = alara_clearance_es_table("non_hazardous", "slags")
print(len(es))

inventory = {"Ra-226": 0.2, "K-40": 5.0}  # Bq/g
print(alara_sum_of_fractions(inventory, es)["class"])
```

```text
40
satisfied
```

## The class boundary

The boundary sits exactly at `sum == 1` and belongs to the satisfied side: an
inventory sitting exactly at its combined limits classifies as satisfied, and
one float step above flips the class:

```python
toy_limits = {"Co60": 10.0, "H3": 5.0, "Fe55": 2.0}

at = alara_sum_of_fractions({"Co60": 10.0}, toy_limits)
print(at["sum"], at["class"])
above = alara_sum_of_fractions({"Co60": 10.000000000000002}, toy_limits)
print(above["sum"], above["class"])
```

```text
1.0 satisfied
1.0000000000000002 exceeded
```

## A custom limit table

Pass a `{nuclide: limit}` dict as `limits` to screen against site-specific or
national values instead of the EU table. Activities and limits must share one
unit basis — here both sides are total activities in Bq, so the comparison is
consistent even though the numbers look nothing like Bq/g concentrations:

```python
print(alara_clearance_index({"Co60": 10.0, "H3": 2.5}, toy_limits))
```

```text
1.5
```

Every inventory nuclide needs a table entry, and every activity must be finite
and non-negative; violations raise a clear error instead of being skipped:

```python
try:
    alara_clearance_index({"Mn54": 1.0}, toy_limits)
except ValueError as exc:
    print(exc)
```

```text
nuclide Mn54 has no clearance limit (table holds 3 entries)
```

## From a parsed inventory to a screening result

Tying the two halves together: take one time step from the parsed FISPACT-II
block and screen it against a caller-supplied limit table in the matching
basis (total Bq on both sides):

```python
step1_inv = {r["nuclide"]: r["activity_bq"] for r in rows if r["interval"] == 1}
limits = {"V-55": 1.0, "Fe-56": 1.0, "Co-60": 4.0e6, "Rb-86m": 8.0e9}
out = alara_sum_of_fractions(step1_inv, limits)
print(out["sum"], out["class"], out["max_nuclide"], out["max_fraction"])
```

```text
2.0 exceeded Co-60 1.0
```

## Screening arithmetic, not a compliance decision

These helpers are screening arithmetic, not a compliance decision. Real
clearance calls for the governing regulatory table for the material and
pathway, complete material bookkeeping, and the national transposition of the
underlying directive. The vendored EU table is a transcription of one
directive's annex for solid materials, and the Spanish tables transcribe the
CSN conditional NORM landfill levels; a sum below 1 marks an inventory as
passing this arithmetic check — not as cleared for release.

## See also

- `tests/test_clearance.py` for the boundary probes, the error cases, and an
  end-to-end parse-and-screen replay.
- [Sample data files](../../reference/fixtures.mdx) for `fixtures/fispact/`,
  including the clearance-block sample used above.
- [Python API reference](../../reference/python-api.mdx#nucleidealara) for
  the full `nucleide.alara` and `nucleide.fispact` signatures.
