"""Nucleide — a memory-safe toolkit for the nuclear-engineering glue layer.

Rust core (crates/*) exposed through PyO3; Python stays the user-facing API.
Functionality is grouped into domain submodules mirroring the Rust crates:
`nucleide.nuclei`, `nucleide.material`, `nucleide.mcnp`, `nucleide.serpent`,
`nucleide.fluka`, `nucleide.vr`, `nucleide.enrichment`,
`nucleide.depletion`, `nucleide.kinetics`, `nucleide.spectroscopy`,
`nucleide.tritium`, `nucleide.alara`, `nucleide.cccc`,
`nucleide.fispact`, `nucleide.origen`, `nucleide.r2s`,
`nucleide.plasma_source`, and `nucleide.damage`.
`nucleide.blanket` books TBR and blanket power over caller transport
tallies (port penalties, energy multiplication, tritium burn and
fuel-cycle margin).
`nucleide.mcpl` reads/writes MCPL particle lists.
`nucleide.unfold` adjusts a guess neutron spectrum against measured
activation rates (SAND-II).
`nucleide.openmc` bridges OpenMC statepoint tallies to plain arrays
(caller-side OpenMC, never a hard dependency).
`nucleide.equilib` reads equilibrium data (classic-netCDF VMEC `wout`
files and `&INDATA` input text) with flux-surface Jacobian helpers —
reads data, never solves equilibria.
`nucleide.uq` is the seeded UQ-lite sampling kernel (caller-supplied blocks).
`nucleide.data` fetches repo data files (compendium, sample chains) pinned
to the installed release.
"""

from nucleide import (
    alara,
    blanket,
    cccc,
    damage,
    data,
    depletion,
    emit,
    enrichment,
    equilib,
    fispact,
    fluka,
    kinetics,
    material,
    mcnp,
    mcpl,
    nuclei,
    openmc,
    origen,
    plasma_source,
    r2s,
    serpent,
    spectroscopy,
    tritium,
    unfold,
    uq,
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
    "tritium",
    "cccc",
    "emit",
    "fispact",
    "origen",
    "r2s",
    "mcpl",
    "uq",
    "data",
    "plasma_source",
    "damage",
    "unfold",
    "blanket",
    "equilib",
    "openmc",
]
