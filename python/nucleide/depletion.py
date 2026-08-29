"""Depletion chains and CRAM solvers (backed by the `depletion` crate)."""

from nucleide._internal import (
    Chain,
    DepletionSystem,
    build_depletion_system,
    deplete,
    read_chain,
)

__all__ = [
    "Chain",
    "DepletionSystem",
    "read_chain",
    "build_depletion_system",
    "deplete",
]
