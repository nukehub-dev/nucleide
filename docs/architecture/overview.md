# Architecture Overview

Nucleide is a Rust workspace with a Python facade. The architecture is
intentionally simple: each crate owns one capability area, and the Python
bindings are a thin re-export layer.

## Layers

```text
┌─────────────────────────────────────┐
│  Python facade (python/nucleide/)   │  <- typed stubs + re-exports
├─────────────────────────────────────┤
│  PyO3 bindings (bindings/python/)   │  <- nucleide._internal
├─────────────────────────────────────┤
│  Capability crates (crates/*)       │  <- parsers, materials, depletion, ...
├─────────────────────────────────────┤
│  Linear-algebra facade (linalg)     │  <- isolates numeric backend choice
└─────────────────────────────────────┘
```

## Design principles

- **Memory safety first.** Core logic is Rust; Python stays the user-facing API.
- **One crate per concern.** Parsers, materials, depletion, enrichment, and
  variance reduction each have their own crate and test suite.
- **Layering is enforced.** `bindings/python` depends on workspace crates;
  workspace crates never depend on bindings or on Python.
- **Parser parity.** Where a reader reproduces legacy output byte-for-byte,
  that behavior is intentional and protected by golden-byte fixtures.
- `enrichment` stays independent of `material` by design.

## Request/data flow

1. User calls `nucleide.read_*` or constructs a `Nuclide` from Python.
2. The pure-Python facade forwards to `nucleide._internal`.
3. The PyO3 extension calls into the relevant workspace crate.
4. The crate parses, computes, and returns a Rust value that PyO3 converts back
   to Python.

## See also

- [Crate responsibilities](crate-responsibilities.md) for crate-level details.
- [Crate overview](../reference/crate-overview.md) for a one-line summary of every crate.
