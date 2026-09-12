"""Enrichment cascade solving (backed by the `nucleide-enrichment` crate)."""

from nucleide._internal import (
    Cascade,
    enrichment_swu_per_feed,
    enrichment_swu_per_prod,
    enrichment_swu_per_tail,
    enrichment_value_func,
)

__all__ = [
    "Cascade",
    "value_func",
    "swu_per_feed",
    "swu_per_prod",
    "swu_per_tail",
]


def value_func(x: float) -> float:
    """Dirac separation potential ``V(x) = (2x - 1) ln(x / (1 - x))``."""
    return enrichment_value_func(x)


def swu_per_feed(x_feed: float, x_prod: float, x_tail: float) -> float:
    """SWU per unit mass of feed for assays ``x_feed``, ``x_prod``, ``x_tail``."""
    return enrichment_swu_per_feed(x_feed, x_prod, x_tail)


def swu_per_prod(x_feed: float, x_prod: float, x_tail: float) -> float:
    """SWU per unit mass of product for assays ``x_feed``, ``x_prod``, ``x_tail``."""
    return enrichment_swu_per_prod(x_feed, x_prod, x_tail)


def swu_per_tail(x_feed: float, x_prod: float, x_tail: float) -> float:
    """SWU per unit mass of tails for assays ``x_feed``, ``x_prod``, ``x_tail``."""
    return enrichment_swu_per_tail(x_feed, x_prod, x_tail)
