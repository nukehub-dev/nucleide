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
from nucleide.mcpl import read_mcpl, write_mcpl

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

write_mcpl("demo.mcpl", header, particles)

mc = read_mcpl("demo.mcpl")
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
from nucleide.mcpl import mcpl2ssw, ssw2mcpl

n = ssw2mcpl("surface.w", "surface.mcpl", [100, 200], ["neutron", "gamma"])
m = mcpl2ssw("surface.mcpl", "surface.w", "back.w")  # reference header cloned
```

`mcpl2ssw` clones the reference SSW header (code/version/deck passthrough,
counts patched to the converted tally) with surface ids from each particle's
`userflags` unless `surface=` overrides them. Anything outside neutrons
(PDG 2112) and gammas (PDG 22) is an error, never a silent skip; see the
`ssw` module docs for the mapping table and the named-open list.

Cross-tool byte compatibility beyond self-consistent round-trips is
oracle-gated in `validation/mcpl_vs_refs.py` (loud SKIP when the upstream
`mcpl` package is absent).
