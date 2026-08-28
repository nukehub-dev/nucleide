# Nucleide

Nucleide is a modern Rust toolkit for nuclear-engineering workflow glue:
legacy transport-code I/O, nuclide identification, materials, CRAM depletion,
and enrichment analytics — exposed through a typed Python API.

The project is a fresh Rust implementation of capabilities pioneered by
[PyNE](https://github.com/pyne/pyne), focused on memory safety, fast builds,
and `pip install`-able wheels. Scope is intentionally narrow today and will
expand as more parsers and workflow pieces land.

## Why

Nuclear-engineering workflows spend most of their time moving data between
codes rather than solving transport itself. The established tooling for that
glue layer carries a heavy build chain (CMake + Fortran + Cython) and
hand-written parsers that are hard to extend and harder to embed. `Nucleide`
rebuilds the high-value subset in memory-safe Rust with one-command
`pip install` wheels, keeping Python as the user-facing API.

## Features

| Area | Capabilities |
| --- | --- |
| Nuclide core (`nuclei`) | Canonical nucid representation, particle registry, reaction-name registry (labels, MT mapping, hashes), name-dialect conversions (ZZAAAMM, ZAID/MCNP, Serpent, FLUKA, NIST, CINDER, ALARA, SZA), AME2020 masses, natural abundances, half-lives |
| Materials (`material`) | Compositions, mixing arithmetic, unit conversions, DOE/PNNL Materials Compendium loading, materials XML export |
| MCNP I/O (`mcnp-io`) | xsdir, meshtal, SSW/SURFSRC, PTRAC, WWINP, MCTAL readers; material extraction from input decks; mesh-to-geometry deck generation |
| Serpent I/O (`serpent-io`) | `_res.m`, `_dep.m`, `_det.m` readers producing structured records |
| FLUKA I/O (`fluka-io`) | USRBIN tally reader, material/compound card generation |
| Depletion (`depletion`) | CRAM (orders 16/48) matrix exponential, depletion-chain XML parsing |
| Enrichment (`enrichment`) | Multicomponent cascade solver (numeric + assignment), SWU closed-form helpers |
| Variance reduction (`vr-tools`) | MAGIC weight-window generation, mesh source sampling with alias tables |
| Python bindings | PyO3 extension module behind a typed pure-Python facade (`nucleide._internal`, `.pyi` stubs, `py.typed`) |

## Out of scope

Transport solvers, Fortran discrete-ordinates ports, ENSDF evaluators, MOAB-dependent
meshing, and GUIs. `Nucleide` complements transport codes; it does not replace them.

## Layout

```text
nucleide/
├── crates/
│   ├── nuclei/        # nuclide ids, naming conventions, physical data
│   ├── material/      # compositions, mixing, libraries, XML export
│   ├── mcnp-io/       # xsdir/meshtal/SSW/MCTAL/PTRAC/WWINP
│   ├── serpent-io/    # res/dep/det readers
│   ├── fluka-io/      # usrbin reader, material cards
│   ├── vr-tools/      # MAGIC weight windows, source sampling
│   ├── enrichment/    # cascades, SWU
│   ├── depletion/     # CRAM + chain files
│   └── linalg/        # isolation facade over the linear-algebra backend
├── bindings/python/   # PyO3 crate -> nucleide._internal
├── python/nucleide/   # typed pure-Python facade (maturin mixed layout)
├── fixtures/          # golden-byte test data
└── tests/             # Python-side tests
```

## Development

```bash
git clone https://github.com/nukehub-dev/nucleide.git
cd nucleide

# Rust side
cargo test                       # workspace unit tests
cargo clippy --all-targets -- -D warnings

# Python side (needs: rustup, pip install maturin)
pip install maturin pytest pytest-cov ruff mypy
maturin develop                  # build + install into current venv
pytest tests/

# Lint / type-check / format the Python surface
ruff format python tests
ruff check python tests
mypy                             # strict type-check against .pyi stubs

# Rust coverage (needs llvm-tools-preview component)
cargo llvm-cov --workspace       # or --lcov for CI upload
```

## Tooling

| Layer | Format | Lint | Types | Coverage |
| --- | --- | --- | --- | --- |
| Rust | rustfmt (`cargo fmt`) | clippy `-D warnings` | — | cargo-llvm-cov (CI) |
| Python | ruff format | ruff check | mypy `--strict` via `.pyi` stubs | pytest-cov |

Wheels are built with **maturin** (PyO3 mixed layout). One wheel serves all
Python >= 3.10 via abi3 — the same stack used by pydantic-core, polars, and ruff.

## Validation strategy

1. Parsers are validated against **golden-byte fixtures** in `fixtures/`;
   parser output must match recorded snapshots before any release.
2. Numeric kernels (CRAM, cascade solving) are checked against published
   analytic vectors and cross-code results on shared inputs.
3. Behavioral compatibility with legacy tool output is asserted wherever a
   fixture exists, so downstream workflows see identical data.

## Status

Pre-alpha. APIs may change without notice.

## Documentation

Additional tutorials, reference pages, and developer guides live in the
[`docs/`](docs/README.md) tree.

## Acknowledgments

Nucleide is a fresh Rust implementation of workflow-glue capabilities pioneered
by [PyNE](https://github.com/pyne/pyne) ("Python for Nuclear Engineering",
BSD-3-Clause). Some reference data and golden test fixtures — notably the
DOE/PNNL Materials Compendium — are vendored directly from PyNE; see
`fixtures/data/MaterialsCompendium.LICENSE` for its terms.

## License

[BSD-2-Clause](LICENSE).
