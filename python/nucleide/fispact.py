"""FISPACT-II output parser (backed by the `nucleide-fispact-io` crate).

Reads FISPACT-II inventory tables into ALARA-compatible response frames so
ALARA and FISPACT-II results share one analysis shape. No activation solving
is performed here.
"""

from nucleide._internal import fispact_is_output, fispact_parse_output

# Short alias for the single-parser module.
parse_output = fispact_parse_output

__all__ = ["fispact_parse_output", "parse_output", "fispact_is_output"]
