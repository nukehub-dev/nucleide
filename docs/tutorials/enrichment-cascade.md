---
title: Enrichment Cascade
sidebar:
  order: 6
---

Nucleide solves multicomponent enrichment cascades with a numeric solver and
closed-form SWU helpers. The implementation lives in `crates/enrichment`. For
the MARC model and SWU derivation, see the
[Enrichment cascades theory](../theory/enrichment.mdx) page.

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

## See also

- [`crates/enrichment/src/lib.rs`](https://github.com/nukehub-dev/nucleide/blob/main/crates/enrichment/src/lib.rs)
  for cascade construction details.
- `tests/test_data_inp_enrichment.py` for integration examples.
