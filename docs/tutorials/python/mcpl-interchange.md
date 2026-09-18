---
title: MCPL Particle Interchange
sidebar:
  order: 12
---

Nucleide reads and writes MCPL particle lists — the code-neutral interchange
format for particle tracks between Monte Carlo codes (see the
[MCPL format page](https://mctools.github.io/mcpl/format/) and Kittelmann et
al., *Computer Physics Communications* 218 (2017), pp. 17–42). This tutorial
covers the Python API; the Rust types live in the `mcpl-io` crate.

Particle units follow the published format: kinetic energy in MeV, position
in cm, time in ms.

```python
import os
import tempfile
from nucleide.mcpl import read_mcpl, write_mcpl

tmp = tempfile.mkdtemp(prefix="nucleide-mcpl-demo-")
demo_path = os.path.join(tmp, "demo.mcpl")

header = {
    "srcname": "demo",
    "comments": [],
    "has_userflags": False,
    "has_polarisation": False,
    "double_prec": False,
    "universal_pdgcode": None,
    "universal_weight": None,
    "blobs": [],
}
particles = [
    {
        "ekin": 2.5,
        "polarisation": [0.0, 0.0, 0.0],
        "position": [0.0, 0.0, 0.0],
        "direction": [0.0, 0.0, 1.0],
        "time": 0.0,
        "weight": 1.0,
        "pdgcode": 2112,
        "userflags": 0,
    }
]

write_mcpl(demo_path, header, particles)

mc = read_mcpl(demo_path)
print(mc.version, mc.nparticles, mc.srcname)
for p in mc.particles():
    print(p["ekin"], p["pdgcode"], p["direction"])
```

A `.gz` suffix reads and writes through gzip transparently (compressed bytes
are encoder-dependent and never asserted). Files with per-file strides —
universal PDG codes or weights, polarisation vectors, user flags, and double
precision — are selected with the matching header flags; `McplFile` exposes
each flag (`has_userflags`, `has_polarisation`, `double_prec`,
`universal_pdgcode`, `universal_weight`) plus the header `comments` and
`blobs` verbatim.

SSW↔MCPL conversion (neutron/gamma-only v1) runs through
`nucleide.mcpl.ssw2mcpl` and `nucleide.mcpl.mcpl2ssw`. The SSW format stores
no per-track surface id or particle kind, so `ssw2mcpl` takes one explicit
surface id plus `"neutron"`/`"gamma"` kind per track (lengths must match the
file's track count); energy maps verbatim in MeV, time maps shakes→ms, and
an options dict selects double precision, surf→userflags, gzip, a deck-embed
blob, `srcname`, and `comments`:

```python
import tempfile
from nucleide.mcpl import mcpl2ssw, ssw2mcpl

tmp = tempfile.mkdtemp(prefix="nucleide-mcpl-")
mcpl_path = f"{tmp}/surface.mcpl"
ref_w = "fixtures/mcpl/ssw_conversion/reference.w"
n = ssw2mcpl(ref_w, mcpl_path, [100, 200], ["neutron", "gamma"])
m = mcpl2ssw(mcpl_path, ref_w, f"{tmp}/surface.w")  # reference header cloned
```

`mcpl2ssw` clones the reference SSW header (code/version/deck passthrough,
counts patched to the converted tally) with surface ids from each particle's
`userflags` unless `surface=` overrides them. Anything outside neutrons
(PDG 2112) and gammas (PDG 22) is an error, never a silent skip; see the
`ssw` module docs for the mapping table and the list of cases that are not
yet supported (each raises a clear error).

## Merging, extracting, and statistics

Four utilities cover the day-to-day particle-list workflows (matching the
upstream `mcpltool` merge/extract/repair surface plus record statistics):

- `merge_mcpl(paths, out_path)` concatenates compatible files. The first
  file's header wins (`srcname`, `comments`, `blobs`) and a provenance
  comment is appended; header statistics comments (`stat:sum:...`) from the
  first file ride along verbatim and sums are never synthesized or updated.
  All inputs must agree on the header options except floating-point
  precision, which promotes to double when single- and double-precision
  files are mixed (the lossless direction); anything else raises a clear
  error.
- `extract_mcpl(src_path, out_path, options=None)` writes a subset with the
  source header preserved verbatim. The options dict selects either an index
  range (`{"start": 1, "stop": 3}` — the half-open interval `[start, stop)`,
  each bound optional) or a predicate (`{"predicate": f}` with `f` over one
  particle dict); omitting options copies the whole file.
- `mcpl_stats(path)` returns record counts, energy moments (sum/min/max/mean
  in MeV), the total weight, and a per-PDG-code histogram.
- `repair_mcpl(path)` fixes a file whose writer never updated the particle
  count (for example after an interrupted job): the count is recomputed
  from the complete records, a partially written trailing record is
  ignored, and the file is rewritten in place.

```python
from nucleide.mcpl import extract_mcpl, mcpl_stats, merge_mcpl

combined = f"{tmp}/combined.mcpl"
n = merge_mcpl([mcpl_path, mcpl_path], combined)
m = extract_mcpl(combined, f"{tmp}/neutrons.mcpl", {"predicate": lambda p: p["pdgcode"] == 2112})
stats = mcpl_stats(f"{tmp}/neutrons.mcpl")
print(stats["nparticles"], stats["ekin_mean"], stats["pdg_counts"])
```

Cross-tool byte compatibility beyond self-consistent round-trips is
cross-checked against the upstream `mcpl` package in the automated
[validation harness](../../development/validation.md); the check is skipped
when that package is not installed.
