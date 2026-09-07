"""Depletion chains and CRAM solvers (backed by the `nucleide-depletion` crate)."""

from nucleide._internal import (
    Chain,
    DepletionSystem,
    Inventory,
    branching_fraction,
    build_depletion_system,
    chain_edges,
    cumulative_decays,
    decay_mode,
    deplete,
    deplete_series,
    progeny,
    read_chain,
)

__all__ = [
    "Chain",
    "DepletionSystem",
    "Inventory",
    "read_chain",
    "build_depletion_system",
    "deplete",
    "deplete_series",
    "cumulative_decays",
    "progeny",
    "branching_fraction",
    "decay_mode",
    "chain_edges",
]
