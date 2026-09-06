"""CCCC binary-standard readers and PARTISN deck writer (backed by the `nucleide-cccc-io` crate).

Thin glue over the documented text-analog subset: ISOTXS multigroup
libraries, RTFLUX/ATFLUX/RZFLUX flux files, and a minimal PARTISN input
writer with ISOTXS nuclide mapping. No transport solving is performed here.
"""

from nucleide._internal import (
    isotxs_parse,
    partisn_render,
    partisn_validate,
    rtflux_parse,
)

# Prefixed aliases mirroring the `alara_*` naming convention.
cccc_parse_isotxs = isotxs_parse
cccc_parse_rtflux = rtflux_parse
cccc_render_partisn = partisn_render
cccc_validate_partisn = partisn_validate

__all__ = [
    "isotxs_parse",
    "rtflux_parse",
    "partisn_render",
    "partisn_validate",
    "cccc_parse_isotxs",
    "cccc_parse_rtflux",
    "cccc_render_partisn",
    "cccc_validate_partisn",
]
