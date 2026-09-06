"""ALARA activation-code interop (backed by the `nucleide-alara-io` crate).

Thin glue over ALARA input decks, group-flux files, activation-output
listings, and schedule expansion. The solver itself stays out of scope.

Schedule-expansion choice: :func:`alara_expand_schedule` takes deck text
(plus an optional top schedule name) instead of JSON schedule/history blobs,
so callers reuse the already-parsed deck blocks without a parallel schema.
"""

from nucleide._internal import (
    alara_expand_schedule,
    alara_parse_deck,
    alara_parse_flux,
    alara_parse_output,
)

__all__ = [
    "alara_parse_deck",
    "alara_parse_flux",
    "alara_parse_output",
    "alara_expand_schedule",
]
