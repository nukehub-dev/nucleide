---
title: Python API
sidebar:
  order: 2
---

Nucleide exposes a typed pure-Python package at `python/nucleide/` that re-exports
a compiled PyO3 extension built from `bindings/python/`. The extension module is
named `nucleide._internal`; users import from the `nucleide` domain submodules,
which mirror the Rust workspace crates.

## Module layout

```python
import nucleide as nuc
```

The top level carries `version()`, `nucleide.__version__`, and the domain
submodules below (so `nuc.nuclei.Nuclide(...)` works directly). Each submodule
re-exports its symbols from `nucleide._internal`:

<!-- GEN:module-table:START -->

| Submodule | Backing crate | Contents |
| --- | --- | --- |
| `nucleide.nuclei` | `nuclei` | Nuclide identifiers, nuclear data, and reaction names |
| `nucleide.material` | `material` | Material composition, activation, and the PNNL compendium |
| `nucleide.mcnp` | `mcnp-io` | MCNP file readers and writers |
| `nucleide.serpent` | `serpent-io` | Serpent output readers |
| `nucleide.fluka` | `fluka-io` | FLUKA USRBIN readers |
| `nucleide.vr` | `vr-tools` | Variance-reduction tools |
| `nucleide.enrichment` | `enrichment` | Enrichment cascade solving |
| `nucleide.depletion` | `depletion` | Depletion chains and CRAM solvers |
| `nucleide.alara` | `alara-io` | ALARA activation-code interop |
| `nucleide.cccc` | `cccc-io` | CCCC binary-standard readers and PARTISN deck writer |
| `nucleide.fispact` | `fispact-io` | FISPACT-II output parser |
| `nucleide.origen` | `origen-io` | Scoped ORIGEN 2.2 TAPE readers |
| `nucleide.r2s` | `r2s` | Rigorous two-step (R2S) shutdown-dose-rate orchestration |
| `nucleide.data` | — (pure Python) | Download Nucleide data files pinned to a release tag, branch, or commit. |

<!-- GEN:module-table:END -->

## `nucleide.nuclei`

<!-- GEN:module-nuclei:START -->

- `Nuclide` — A nuclide identifier with naming-convention conversions.
  - `__init__(name: str) -> None`
  - `name: str` (property)
  - `nucid: int` (property)
  - `zzaaam: int` (property)
  - `z: int` (property)
  - `a: int` (property)
  - `state: int` (property)
  - `zaid: int` (property)
  - `zzllaaam: str` (property)
  - `serpent: str` (property)
  - `nist: str` (property)
  - `cinder: int` (property)
  - `alara: str` (property)
  - `sza: int` (property)
  - `mass: float | None` (property)
  - `abundance: float | None` (property)
  - `fluka() -> str`

- `Particle` — A particle species with cross-code name translations.
  - `__init__(spec: str | int) -> None`
  - `name: str` (property)
  - `describe: str` (property)
  - `mcnp() -> str | None`
  - `mcnp6() -> str | None`
  - `fluka() -> str | None`
  - `geant4() -> str | None`

- `from_zaid(zaid: int) -> Nuclide`
- `atomic_mass(key: int | str) -> float | None`
- `natural_abundance(key: int | str) -> float | None`
- `half_life(key: int | str) -> float | None`
- `decay_constant(key: int | str) -> float | None`
- `q_value_capture(key: int | str) -> float | None`
- `q_value_alpha(key: int | str) -> float | None`
- `rxname_id(name: str) -> int`
- `rxname_name(id: int) -> str | None`
- `rxname_mt(id: int) -> int`

<!-- GEN:module-nuclei:END -->

## `nucleide.material`

<!-- GEN:module-material:START -->

- `MaterialsCompendium` — PNNL/DOE Materials Compendium (411 named materials).
  - `load(path: str) -> MaterialsCompendium`
  - `__len__() -> int`
  - `names() -> list[str]`
  - `get(name: str, as_material: bool = False) -> dict[str, Any] | None`

- `from_formula(formula: str) -> dict[str, float]`
- `activity(comp: dict[str, float]) -> dict[str, float]`
- `to_xml(comp: dict[str, float], name: str, density: float, units: str = 'g/cm3') -> str`

<!-- GEN:module-material:END -->

## `nucleide.mcnp`

<!-- GEN:module-mcnp:START -->

- `Xsdir` — Parsed MCNP xsdir cross-section index.
  - `datapath: str | None` (property)
  - `awr: dict[int, float]` (property)
  - `tables: list[XsdirTable]` (property)
  - `find_table(name: str) -> list[XsdirTable]`
  - `nucs() -> list[int]`

- `XsdirTable` — One xsdir directory entry.
  - `name: str` (property)
  - `awr: float` (property)
  - `filename: str` (property)
  - `filetype: int` (property)
  - `address: int` (property)
  - `tablelength: int` (property)
  - `temperature: float | None` (property)
  - `ptable: bool` (property)
  - `zaid() -> str`
  - `to_serpent(directory: str = '') -> str`

- `Meshtal` — Parsed MCNP meshtal file.
  - `version: str` (property)
  - `ld: str` (property)
  - `title: str` (property)
  - `histories: int` (property)
  - `tallies: dict[int, MeshTally]` (property)

- `MeshTally` — One fmesh4 tally from a meshtal file.
  - `tally_number: int` (property)
  - `particle: str` (property)
  - `dose_response: bool` (property)
  - `x_bounds: list[float]` (property)
  - `y_bounds: list[float]` (property)
  - `z_bounds: list[float]` (property)
  - `e_bounds: list[float]` (property)
  - `result: list[list[float]]` (property)
  - `rel_error: list[list[float]]` (property)
  - `total_result: list[float]` (property)
  - `total_rel_error: list[float]` (property)
  - `dims() -> tuple[int, int, int]`
  - `num_ves() -> int`
  - `num_e_groups() -> int`
  - `cell(i: int, j: int, k: int) -> tuple[list[float], list[float]]`
  - `cell_total(i: int, j: int, k: int) -> tuple[float, float]`

- `Wwinp` — Parsed MCNP WWINP weight-window file (Cartesian).
  - `ni: int` (property)
  - `nr: int` (property)
  - `ne: list[int]` (property)
  - `nf: tuple[int, int, int]` (property)
  - `origin: tuple[float, float, float]` (property)
  - `nc: tuple[int, int, int]` (property)
  - `cm: list[list[float]]` (property)
  - `bounds: list[list[float]]` (property)
  - `e: list[list[float]]` (property)
  - `ww_row(particle: int, group: int) -> list[float]`
  - `ww_column(particle: int, ve: int) -> list[float]`

- `Mctal` — Parsed MCNP MCTAL kcode data (upstream subset).
  - `code_name: str` (property)
  - `comment: str` (property)
  - `n_histories: int` (property)
  - `n_cycles: int` (property)
  - `n_inactive: int` (property)
  - `vars_per_cycle: int` (property)
  - `k_col: list[float]` (property)
  - `k_abs: list[float]` (property)
  - `k_path: list[float]` (property)
  - `prompt_life_col: list[float]` (property)
  - `prompt_life_path: list[float]` (property)
  - `averages: list[dict[str, float]]` (property)

- `SurfSrc` — Parsed MCNP SSW surface-source file.
  - `kod: str` (property)
  - `ver: str` (property)
  - `np1: int` (property)
  - `nrss: int` (property)
  - `ncrd: int` (property)
  - `njsw: int` (property)
  - `niss: int` (property)
  - `print_header() -> str`
  - `tracks() -> list[dict[str, float]]`

- `PtracFile` — Parsed MCNP PTRAC event file.
  - `problem_title: str` (property)
  - `width_code: int` (property)
  - `variable_nums: dict[str, int]` (property)
  - `events() -> list[dict[str, float]]`

- `read_xsdir(path: str) -> Xsdir`
- `read_meshtal(path: str) -> Meshtal`
- `read_wwinp(path: str) -> Wwinp`
- `read_mctal(path: str) -> Mctal`
- `read_ssw(path: str) -> SurfSrc`
- `read_ptrac(path: str) -> PtracFile`
- `write_ssw(ssw: SurfSrc, path: str, tracks: list[dict[str, float]] | None = None) -> None`
- `read_inp(path: str) -> list[dict[str, Any]]`

- `mesh_to_geom`:

  ```python
  mesh_to_geom(x_bounds: Sequence[float], y_bounds: Sequence[float], z_bounds: Sequence[float], cell_materials: Sequence[tuple[str, float] | None], title_card: str) -> str
  ```

<!-- GEN:module-mcnp:END -->

## `nucleide.serpent`

<!-- GEN:module-serpent:START -->

- `read_serpent(path: str, kind: str) -> dict[str, Any]`

<!-- GEN:module-serpent:END -->

## `nucleide.fluka`

<!-- GEN:module-fluka:START -->

- `UsrbinTally` — One FLUKA USRBIN detector.
  - `name: str` (property)
  - `particle: str` (property)
  - `nx: int` (property)
  - `ny: int` (property)
  - `nz: int` (property)
  - `x_bounds: list[float]` (property)
  - `y_bounds: list[float]` (property)
  - `z_bounds: list[float]` (property)
  - `data: list[float]` (property)
  - `error: list[float]` (property)
  - `dims() -> tuple[int, int, int]`

- `read_usrbin(path: str) -> list[UsrbinTally]`

<!-- GEN:module-fluka:END -->

## `nucleide.vr`

<!-- GEN:module-vr:START -->

- `magic(tally: MeshTally, per_group: bool = False, tolerance: float = 0.5) -> MagicOutput`
- `MagicOutput` — MAGIC weight-window generation output.
  - `lower_bounds_ww: list[float]` (property)
  - `groups_per_ve: int` (property)
  - `scale_factors: list[float]` (property)
  - `e_upper_bounds: list[float]` (property)
  - `ww_tag_name: str` (property)

- `AliasTable` — Walker alias table for discrete sampling.
  - `__init__(pdf: list[float]) -> None`
  - `sample(r1: float, r2: float) -> int`
  - `pdf: list[float]` (property)
  - `__len__() -> int`

- `MeshSourceSampler` — Mesh source sampler over a meshtal tally (analog/uniform/user).
  - `__init__(tally: MeshTally, mode: str, user_pdf: list[float] | None = None) -> None`
  - `sample(r1: float, r2: float) -> dict[str, float]`

<!-- GEN:module-vr:END -->

## `nucleide.enrichment`

<!-- GEN:module-enrichment:START -->

- `Cascade` — Enrichment cascade with numeric multicomponent solving.
  - `default_uranium() -> Cascade`

  - `__init__`:

    ```python
    __init__(alpha: float, Mstar: float, j: int, k: int, N: float, M: float, x_feed_j: float, x_prod_j: float, x_tail_j: float, mat_feed: dict[str, float]) -> None
    ```

  - `solve(tolerance: float | None = None, max_iterations: int | None = None) -> None`
  - `solve_multicomponent(tolerance: float | None = None, max_iterations: int | None = None) -> None`
  - `alpha: float` (property)
  - `Mstar: float` (property)
  - `N: float` (property)
  - `M: float` (property)
  - `x_feed_j: float` (property)
  - `x_prod_j: float` (property)
  - `x_tail_j: float` (property)
  - `l_t_per_feed: float` (property)
  - `swu_per_feed: float` (property)
  - `swu_per_prod: float` (property)
  - `mat_feed: dict[str, float]` (property)
  - `mat_prod: dict[str, float]` (property)
  - `mat_tail: dict[str, float]` (property)
  - `separative_work_per_product() -> float`

<!-- GEN:module-enrichment:END -->

## `nucleide.depletion`

<!-- GEN:module-depletion:START -->

- `Chain` — Parsed depletion chain (XML format).
  - `nuclides: list[str]` (property)
  - `index_of(name: str) -> int | None`

- `DepletionSystem` — Pre-built depletion system for repeated CRAM solves.
  - `solve(n0: dict[str, float], dt: float, order: int = 48) -> dict[str, float]`
  - `solve_vec(n0: list[float], dt: float, order: int = 48) -> list[float]`

- `read_chain(path: str) -> Chain`
- `build_depletion_system(chain: Chain, rates: dict[str, float]) -> DepletionSystem`

- `deplete`:

  ```python
  deplete(chain: Chain, n0: dict[str, float], dt: float, rates: dict[str, float] | None = None, order: int = 48) -> dict[str, float]
  ```

<!-- GEN:module-depletion:END -->

## `nucleide.alara`

<!-- GEN:module-alara:START -->

- `alara_parse_deck(text: str) -> dict[str, Any]`
- `alara_parse_flux(text: str, name: str) -> dict[str, Any]`
- `alara_parse_output(text: str, run_lbl: str) -> list[dict[str, Any]]`
- `alara_expand_schedule(deck_text: str, top: str | None = None) -> list[dict[str, Any]]`

<!-- GEN:module-alara:END -->

## `nucleide.cccc`

<!-- GEN:module-cccc:START -->

- `isotxs_parse(text: str) -> dict[str, Any]`
- `rtflux_parse(text: str, kind: str = 'rtflux') -> dict[str, Any]`
- `partisn_render(deck: dict[str, Any]) -> str`
- `partisn_validate(deck: dict[str, Any], isotxs_text: str) -> None`
- `cccc_parse_isotxs(text: str) -> dict[str, Any]` (alias of `isotxs_parse`)
- `cccc_parse_rtflux(text: str, kind: str = 'rtflux') -> dict[str, Any]` (alias of `rtflux_parse`)
- `cccc_render_partisn(deck: dict[str, Any]) -> str` (alias of `partisn_render`)
- `cccc_validate_partisn(deck: dict[str, Any], isotxs_text: str) -> None` (alias of `partisn_validate`)

<!-- GEN:module-cccc:END -->

## `nucleide.fispact`

<!-- GEN:module-fispact:START -->

- `fispact_parse_output(text: str, run_lbl: str) -> list[dict[str, Any]]`
- `parse_output(text: str, run_lbl: str) -> list[dict[str, Any]]` (alias of `fispact_parse_output`)

<!-- GEN:module-fispact:END -->

## `nucleide.origen`

<!-- GEN:module-origen:START -->

- `origen_parse_tape5(text: str) -> dict[str, Any]`
- `origen_parse_tape6(text: str) -> dict[str, Any]`
- `origen_parse_tape9(text: str) -> list[dict[str, Any]]`
- `tape5_parse(text: str) -> dict[str, Any]` (alias of `origen_parse_tape5`)
- `tape6_parse(text: str) -> dict[str, Any]` (alias of `origen_parse_tape6`)
- `tape9_parse(text: str) -> list[dict[str, Any]]` (alias of `origen_parse_tape9`)

<!-- GEN:module-origen:END -->

## `nucleide.r2s`

<!-- GEN:module-r2s:START -->

- `r2s_from_deck(deck_text: str) -> dict[str, Any]`
- `r2s_validate(workflow: dict[str, Any], deck_text: str) -> None`
- `r2s_expand(deck_text: str, top: str | None = None) -> list[dict[str, Any]]`
- `r2s_assemble(output_text: str, run_lbl: str, zone: str, groups: int) -> dict[str, Any]`
- `from_deck(deck_text: str) -> dict[str, Any]` (alias of `r2s_from_deck`)
- `validate(workflow: dict[str, Any], deck_text: str) -> None` (alias of `r2s_validate`)
- `expand(deck_text: str, top: str | None = None) -> list[dict[str, Any]]` (alias of `r2s_expand`)
- `assemble(output_text: str, run_lbl: str, zone: str, groups: int) -> dict[str, Any]` (alias of `r2s_assemble`)

<!-- GEN:module-r2s:END -->

## `nucleide.data`

Pure-Python helpers (no backing crate) for downloading repo data files that
the wheel does not bundle, pinned to the installed release:

<!-- GEN:module-data:START -->

- `COMPENDIUM_PATH = 'fixtures/data/MaterialsCompendium.json'`
- `default_ref() -> str`
- `fetch(path: str, *, ref: str | None = None, dest: str | Path = '.') -> str`
- `fetch_compendium(*, ref: str | None = None, dest: str | Path = '.') -> str`

<!-- GEN:module-data:END -->

## Version

- `version()` → workspace version string
- `nucleide.__version__` — same value

The signatures above are generated from `python/nucleide/_internal.pyi` (the function source of truth)
and the `python/nucleide/*.py` facades; see those files for the authoritative definitions.
