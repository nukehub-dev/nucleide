"""MCPL particle-list interchange (backed by the `nucleide-mcpl-io` crate).

Thin facade over the flat `_internal` extension: use :func:`read_mcpl` to
parse a file (`.gz` reads through gzip transparently) and :func:`write_mcpl`
to emit one from a header dict plus particle dicts. Particle units follow the
published format (kinetic energy in MeV, position in cm, time in ms).
Cross-tool byte compatibility beyond self-consistent round-trips is
oracle-gated (see `validation/mcpl_vs_refs.py`); SSW conversion is deferred.
"""

from nucleide._internal import McplFile, read_mcpl, write_mcpl

__all__ = [
    "McplFile",
    "read_mcpl",
    "write_mcpl",
]
