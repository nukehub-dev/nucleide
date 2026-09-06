"""Variance-reduction tools (backed by the `nucleide-vr-tools` crate)."""

from nucleide._internal import (
    AliasTable,
    MagicOutput,
    MeshSourceSampler,
    magic,
)

__all__ = [
    "magic",
    "MagicOutput",
    "AliasTable",
    "MeshSourceSampler",
]
