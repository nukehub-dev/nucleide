"""Rigorous two-step (R2S) shutdown-dose-rate orchestration (backed by the `nucleide-r2s` crate).

Combines neutron flux meshes, ALARA activation decks, decay photon sources,
and mesh source sampling into one reproducible workflow description.
Transport and activation solving stay inside their respective codes; this
module only builds decks, links zones to spectra, and assembles photon
sources.

Photon-source approximation: :func:`r2s_assemble` sums shutdown
``SpecificActivity`` over a zone and splits the total uniformly over the
requested energy groups. The uniform split preserves only the total shutdown
strength; real decay photons follow the nuclide- and energy-dependent lines
in ALARA ``.photonSrc`` spectra.
"""

from nucleide._internal import (
    r2s_assemble,
    r2s_expand,
    r2s_from_deck,
    r2s_from_snapshot,
    r2s_validate,
)

# Short aliases for the workflow-builder summaries.
from_deck = r2s_from_deck
from_snapshot = r2s_from_snapshot
validate = r2s_validate
expand = r2s_expand
assemble = r2s_assemble

__all__ = [
    "r2s_from_deck",
    "r2s_from_snapshot",
    "r2s_validate",
    "r2s_expand",
    "r2s_assemble",
    "from_deck",
    "from_snapshot",
    "validate",
    "expand",
    "assemble",
]
