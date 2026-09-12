"""Nucleide — a memory-safe toolkit for the nuclear-engineering glue layer.

Rust core (crates/*) exposed through PyO3; Python stays the user-facing API.
Functionality is grouped into domain submodules mirroring the Rust crates:
`nucleide.nuclei`, `nucleide.material`, `nucleide.mcnp`, `nucleide.serpent`,
`nucleide.fluka`, `nucleide.vr`, `nucleide.enrichment`,
`nucleide.depletion`, `nucleide.kinetics`, `nucleide.spectroscopy`,
`nucleide.alara`, `nucleide.cccc`,
`nucleide.fispact`, `nucleide.origen`, and `nucleide.r2s`.
`nucleide.mcpl` reads/writes MCPL particle lists.
`nucleide.data` fetches repo data files (compendium, sample chains) pinned
to the installed release.
"""

from nucleide import (
    alara,
    cccc,
    data,
    depletion,
    emit,
    enrichment,
    fispact,
    fluka,
    kinetics,
    material,
    mcnp,
    mcpl,
    nuclei,
    origen,
    r2s,
    serpent,
    spectroscopy,
    vr,
)
from nucleide._internal import version

__version__ = version()
__all__ = [
    "version",
    "alara",
    "nuclei",
    "material",
    "mcnp",
    "serpent",
    "fluka",
    "vr",
    "enrichment",
    "depletion",
    "kinetics",
    "spectroscopy",
    "cccc",
    "emit",
    "fispact",
    "origen",
    "r2s",
    "mcpl",
    "data",
]
