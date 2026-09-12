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

SSW↔MCPL conversion is deferred: this module only reads and writes MCPL.
Cross-tool byte compatibility beyond self-consistent round-trips is
oracle-gated in `validation/mcpl_vs_refs.py` (loud SKIP when the upstream
`mcpl` package is absent).
