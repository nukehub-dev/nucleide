"""MCPL particle-list interchange (backed by the `nucleide-mcpl-io` crate).

Thin facade over the flat `_internal` extension: use :func:`read_mcpl` to
parse a file (`.gz` reads through gzip transparently) and :func:`write_mcpl`
to emit one from a header dict plus particle dicts. Particle units follow the
published format (kinetic energy in MeV, position in cm, time in ms).
SSW↔MCPL conversion (neutron/gamma-only v1) runs through :func:`ssw2mcpl`
and :func:`mcpl2ssw`: SSW tracks need one explicit surface id + particle
kind (`"neutron"`/`"gamma"`) per track, and `mcpl2ssw` clones a reference
SSW header with `nrss` patched. Other particle kinds are errors, never
silent skips.
Cross-tool byte compatibility beyond self-consistent round-trips is
oracle-gated (see `validation/mcpl_vs_refs.py`).
"""

from nucleide._internal import McplFile, mcpl2ssw, read_mcpl, ssw2mcpl, write_mcpl

__all__ = [
    "McplFile",
    "read_mcpl",
    "write_mcpl",
    "ssw2mcpl",
    "mcpl2ssw",
]
