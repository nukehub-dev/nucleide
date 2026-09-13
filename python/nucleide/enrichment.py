"""Enrichment cascade solving (backed by the `nucleide-enrichment` crate)."""

from nucleide._internal import (
    Cascade,
    enrichment_alphastar_i,
    enrichment_feed_per_prod,
    enrichment_feed_per_tail,
    enrichment_prod_per_feed,
    enrichment_prod_per_tail,
    enrichment_swu_per_feed,
    enrichment_swu_per_prod,
    enrichment_swu_per_tail,
    enrichment_tail_per_feed,
    enrichment_tail_per_prod,
    enrichment_value_func,
)


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


def prod_per_feed(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Product-per-feed mass ratio for assays ``x_feed``, ``x_prod``, ``x_tail``."""
    return enrichment_prod_per_feed(x_feed, x_prod, x_tail)


def tail_per_feed(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Tails-per-feed mass ratio."""
    return enrichment_tail_per_feed(x_feed, x_prod, x_tail)


def tail_per_prod(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Tails-per-product mass ratio."""
    return enrichment_tail_per_prod(x_feed, x_prod, x_tail)


def feed_per_prod(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Feed-per-product mass ratio."""
    return enrichment_feed_per_prod(x_feed, x_prod, x_tail)


def feed_per_tail(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Feed-per-tails mass ratio."""
    return enrichment_feed_per_tail(x_feed, x_prod, x_tail)


def prod_per_tail(x_feed: float, x_prod: float, x_tail: float) -> float:
    """Product-per-tails mass ratio."""
    return enrichment_prod_per_tail(x_feed, x_prod, x_tail)


def alphastar_i(alpha: float, Mstar: float, M_i: float) -> float:
    """Stage separation factor for a component of mass ``M_i``."""
    return enrichment_alphastar_i(alpha, Mstar, M_i)


__all__ = [
    "Cascade",
    "value_func",
    "swu_per_feed",
    "swu_per_prod",
    "swu_per_tail",
    "prod_per_feed",
    "tail_per_feed",
    "tail_per_prod",
    "feed_per_prod",
    "feed_per_tail",
    "prod_per_tail",
    "alphastar_i",
]
