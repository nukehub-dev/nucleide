"""Scoped ORIGEN 2.2 TAPE readers (backed by the `origen-io` crate).

Covers the decay-data path (`TAPE9`-style decay constants), `TAPE5` input
echo, and `TAPE6` output inventories. Full ORIGEN burnup driving stays inside
ORIGEN; this module only reads its text interfaces.
"""

from nucleide._internal import (
    origen_parse_tape5,
    origen_parse_tape6,
    origen_parse_tape9,
)

# Short aliases for the scoped TAPE readers.
tape5_parse = origen_parse_tape5
tape6_parse = origen_parse_tape6
tape9_parse = origen_parse_tape9

__all__ = [
    "origen_parse_tape5",
    "origen_parse_tape6",
    "origen_parse_tape9",
    "tape5_parse",
    "tape6_parse",
    "tape9_parse",
]
