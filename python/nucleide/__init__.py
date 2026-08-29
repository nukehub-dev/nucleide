"""Nucleide — a memory-safe toolkit for the nuclear-engineering glue layer.

Rust core (crates/*) exposed through PyO3; Python stays the user-facing API.
Functionality is grouped into domain submodules mirroring the Rust crates:
`nucleide.nuclei`, `nucleide.material`, `nucleide.mcnp`, `nucleide.serpent`,
`nucleide.fluka`, `nucleide.vr`, `nucleide.enrichment`, and
`nucleide.depletion`.
"""

from nucleide import (
    depletion,
    enrichment,
    fluka,
    material,
    mcnp,
    nuclei,
    serpent,
    vr,
)
from nucleide._internal import version

__version__ = version()
__all__ = [
    "version",
    "nuclei",
    "material",
    "mcnp",
    "serpent",
    "fluka",
    "vr",
    "enrichment",
    "depletion",
]
