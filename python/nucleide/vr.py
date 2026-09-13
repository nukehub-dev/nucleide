"""Variance-reduction tools (backed by the `nucleide-vr-tools` crate)."""

from nucleide._internal import (
    AliasTable,
    MagicOutput,
    MeshSourceSampler,
    magic,
    magic_with,
)

__all__ = [
    "magic",
    "magic_with",
    "MagicOutput",
    "AliasTable",
    "MeshSourceSampler",
]
