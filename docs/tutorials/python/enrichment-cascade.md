---
title: Enrichment Cascade
sidebar:
  order: 6
---

Nucleide solves multicomponent enrichment cascades with a numeric solver and
closed-form SWU helpers. The implementation lives in `crates/enrichment`. For
the MARC model and SWU derivation, see the
[Enrichment cascades theory](../../theory/enrichment.mdx) page.

## Default uranium cascade

```python
from nucleide.enrichment import Cascade

c = Cascade.default_uranium()
c.solve()
print(c.swu_per_feed, c.swu_per_prod)
```

## Custom cascade

Custom cascades are configured on the Rust side by specifying component
molecular weights, assays, and separation factors; the Python API exposes the
result objects.

## Standalone SWU and mass-ratio helpers

Without building a cascade, the closed-form helpers evaluate the Dirac
separation potential, SWU per unit feed/product/tails, and the six assay
mass ratios on any `(x_feed, x_prod, x_tail)` ladder:

```python
from nucleide.enrichment import (
    alphastar_i,
    feed_per_prod,
    prod_per_feed,
    swu_per_feed,
    swu_per_prod,
    tail_per_feed,
    value_func,
)

x_feed, x_prod, x_tail = 0.0072, 0.05, 0.002
print(value_func(x_prod))
print(swu_per_feed(x_feed, x_prod, x_tail))
print(swu_per_prod(x_feed, x_prod, x_tail))
print(prod_per_feed(x_feed, x_prod, x_tail))  # (xf - xt) / (xp - xt)
print(tail_per_feed(x_feed, x_prod, x_tail))  # 1 - prod_per_feed
print(feed_per_prod(x_feed, x_prod, x_tail))  # 1 / prod_per_feed
print(alphastar_i(1.05, 236.0, 235.0))  # stage factor for mass 235
```

`tail_per_prod`, `feed_per_tail`, and `prod_per_tail` complete the ratio set;
each is the quotient of the two corresponding per-feed ratios.

## See also

- [`crates/enrichment/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/enrichment/src/lib.rs)
  for cascade construction details.
- `tests/test_data_inp_enrichment.py` for integration examples.
- [Cross-code validation results](https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md)
  for these cascades, benchmarked against PyNE 0.7.5 `multicomponent`.
