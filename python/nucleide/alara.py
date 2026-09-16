"""ALARA activation-code interop (backed by the `nucleide-alara-io` crate).

Thin glue over ALARA input decks, group-flux files, activation-output
listings, photon sources, schedule expansion, and clearance / waste-
classification analytics (clearance index + sum-of-fractions screening).
The solver itself stays out of scope.

Schedule-expansion choice: :func:`alara_expand_schedule` takes deck text
(plus an optional top schedule name) instead of JSON schedule/history blobs,
so callers reuse the already-parsed deck blocks without a parallel schema.

Clearance screening choice: :func:`alara_clearance_index` and
:func:`alara_sum_of_fractions` take plain ``{nuclide: activity}`` dicts plus
an optional ``{nuclide: limit}`` dict (None = the vendored EU 2013/59/Euratom
Annex VII Table A default, Bq/g). Nuclide names accept any shared-dialect
spelling; activities and limits must share one unit basis. Screening
arithmetic only, never a compliance decision.

The vendored limit tables are selectable explicitly:
:func:`alara_clearance_eu_table` (EU Annex VII Table A) and
:func:`alara_clearance_es_table` (Spanish CSN conditional NORM landfill
levels, Tables 1-3 per landfill type and NORM material nature, with the
Tabla 4 chain keys expanded to per-member entries at the parent value). No
cross-table logic: the caller picks the governing table.
"""

from nucleide._internal import (
    alara_check_block,
    alara_clearance_es_table,
    alara_clearance_eu_table,
    alara_clearance_index,
    alara_expand_schedule,
    alara_flux_len,
    alara_flux_total,
    alara_output_total_activity,
    alara_output_totals,
    alara_parse_deck,
    alara_parse_flux,
    alara_parse_output,
    alara_photon_total_strength,
    alara_schedule_total_time,
    alara_sum_of_fractions,
    alara_validate_deck,
)

__all__ = [
    "alara_parse_deck",
    "alara_parse_flux",
    "alara_parse_output",
    "alara_expand_schedule",
    "alara_validate_deck",
    "alara_check_block",
    "alara_flux_total",
    "alara_flux_len",
    "alara_output_totals",
    "alara_output_total_activity",
    "alara_photon_total_strength",
    "alara_schedule_total_time",
    "alara_clearance_eu_table",
    "alara_clearance_es_table",
    "alara_clearance_index",
    "alara_sum_of_fractions",
]
