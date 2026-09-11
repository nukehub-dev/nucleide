"""Minimal mypy shadow for NumPy (see `mypy_path` in `pyproject.toml`).

Site-packages NumPy (>= 2) bundles stubs using Python 3.12+ `type`-statement
syntax, which mypy cannot parse under this repo's `python_version = "3.10"`
target (fatal, before per-module overrides apply). This shadow resolves
`import numpy` first and declares only the names `python/nucleide/` needs, as
opaque `Any` (runtime behavior is unchanged; the wheel stays
dependency-free and NumPy-backed helpers require NumPy only at call time).
"""

from typing import Any

ndarray: Any
empty: Any
asarray: Any
