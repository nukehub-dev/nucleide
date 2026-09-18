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

Sublet S1/S2 radiological totals: :func:`alara_total_activity` sums
per-nuclide activities with the FISPACT-II IRT alpha/beta/gamma split
(plus the excluding-tritium companion), and :func:`alara_decay_heat` folds
caller-supplied average decay energies into per-class decay heat (kW, plus
the excluding-tritium companion). Pure arithmetic, never vendored data.

Sublet S4/S5 committed hazards: :func:`alara_ingestion_hazard` and
:func:`alara_inhalation_hazard` fold per-nuclide activities with the
caller-supplied 50-year committed ingestion/inhalation dose coefficients
(Sv/Bq, never vendored — the ICRP tables are copyrighted) into the
``TOTAL ... HAZARD FOR ALL MATERIALS`` dose (Sv, plus the
excluding-tritium companion). Pure arithmetic, never vendored data.

Sublet S3 gamma dose rate: :func:`alara_dose_slab` and
:func:`alara_dose_point` evaluate the FISPACT-II ``DOSE`` slab/point kernel
(Sv/h) over caller gamma groups (group yields plus air/mixture attenuation
coefficients — :func:`alara_dose_mixture_mu` folds elemental values with
mixture fractions; no table is vendored). Short point distances clamp to
0.3 m and report it loudly via the ``clamped`` flag, never silently.
Screening arithmetic only, never a compliance decision.

Sublet S6 transport ratio: :func:`alara_transport_ratio` folds per-nuclide
activities with the caller-supplied ``A2`` limits (TBq, never vendored) into
the dimensionless ``Total Bq/A2`` ratio plus the effective A2 it defines.

Sublet S7 IAEA clearance index: :func:`alara_iaea_clearance_index` folds
per-nuclide activities with the caller-supplied IAEA levels (Bq/kg, never
vendored — the IAEA tables are permission-gated) and the total mass into the
dimensionless clearance index with its ``<= 1`` screening class (boundary
included) and dominant contributor. Screening arithmetic only, never a
compliance decision.
"""

from nucleide._internal import (
    alara_check_block,
    alara_clearance_es_table,
    alara_clearance_eu_table,
    alara_clearance_index,
    alara_decay_heat,
    alara_dose_mixture_mu,
    alara_dose_point,
    alara_dose_slab,
    alara_expand_schedule,
    alara_flux_len,
    alara_flux_total,
    alara_iaea_clearance_index,
    alara_ingestion_hazard,
    alara_inhalation_hazard,
    alara_output_total_activity,
    alara_output_totals,
    alara_parse_deck,
    alara_parse_flux,
    alara_parse_output,
    alara_photon_total_strength,
    alara_schedule_total_time,
    alara_sum_of_fractions,
    alara_total_activity,
    alara_transport_ratio,
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
    "alara_total_activity",
    "alara_decay_heat",
    "alara_ingestion_hazard",
    "alara_inhalation_hazard",
    "alara_transport_ratio",
    "alara_iaea_clearance_index",
    "alara_dose_slab",
    "alara_dose_point",
    "alara_dose_mixture_mu",
]
