---
title: VR and MAGIC
sidebar:
  order: 14
---

Nucleide derives MCNP weight-window lower bounds from mesh tallies with the
MAGIC algorithm and samples birth voxels through Walker alias tables. The
implementation lives in `crates/vr-tools`; for the algorithm, see the
[Variance reduction theory](../../theory/variance-reduction.mdx) page.

## MAGIC lower bounds

`magic` runs on energy-integrated totals with default parameters; `magic_with`
selects the array (`"total"` or `"per_group"`) and the tuning parameters
(`tolerance`, `null_value`) explicitly:

```python
from nucleide.mcnp import read_meshtal
from nucleide.vr import magic, magic_with

meshtal = read_meshtal("mcnp_meshtal_single_meshtal.txt")
tally = meshtal.tallies[4]

out = magic(tally)
print(out.groups_per_ve, len(out.lower_bounds_ww))

grouped = magic_with(tally, selection="per_group", tolerance=0.5, null_value=0.0)
print(grouped.groups_per_ve, len(grouped.lower_bounds_ww))
```

Cells whose relative error exceeds `tolerance` receive `null_value` (0.0 by
default) instead of a scaled bound.

## Mesh source sampling

`MeshSourceSampler` builds a birth-voxel sampler over a tally's totals in
`"analog"`, `"uniform"`, or `"user"` mode (user mode takes one unnormalized
density per voxel in `user_pdf`):

```python
from nucleide.vr import AliasTable, MeshSourceSampler

analog = MeshSourceSampler(tally, "analog")
print(analog.mode(), analog.num_voxels(), analog.table_len())
print(analog.sample(0.3, 0.7))  # {"index", "i", "j", "k", "weight"}

uniform = MeshSourceSampler(tally, "uniform")
user = MeshSourceSampler(tally, "user", user_pdf=[1.0] * tally.num_ves())

table = AliasTable([0.9, 0.05, 0.03, 0.02])
print(table.sample(0.3, 0.7))
```

Analog birth weights are unity; biased modes reweight by the analog-to-bias
PDF ratio.

## See also

- `tests/test_serpent_fluka_vr.py` and `tests/test_facade_bundle.py`
  for runnable examples.
