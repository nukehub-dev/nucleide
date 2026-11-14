"""FESTIM caches must never carry solver results across runtime identities."""

from __future__ import annotations

import importlib.util
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import pytest

_MODULE_PATH = Path(__file__).resolve().parents[1] / "validation" / "festim_cache.py"
_SPEC = importlib.util.spec_from_file_location("festim_cache_test_module", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
CACHE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(CACHE)

_PACKAGES = (
    "festim",
    "fenics-dolfinx",
    "petsc4py",
    "mpi4py",
    "fenics-basix",
    "fenics-ufl",
    "fenics-ffcx",
    "scifem",
    "io4dolfinx",
    "fenicsx-compat",
    "numpy",
)


class _Distribution:
    def __init__(self, name: str, version: str = "1.0") -> None:
        self.metadata = {"Name": name}
        self.version = version


class _Calculator:
    def __init__(self) -> None:
        self.calls = 0

    def __call__(self) -> dict[str, Any]:
        self.calls += 1
        return {"computed": self.calls}


def _environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("NUCLEIDE_FESTIM_IMAGE_ID", "sha256:" + "1" * 64)
    distributions = [_Distribution(name) for name in _PACKAGES]
    monkeypatch.setattr(CACHE.importlib.metadata, "distributions", lambda: distributions)
    maps = tmp_path / "maps"
    lines = []
    for name in ("libpetsc.so.3", "libdolfinx.so", "other_native.so"):
        library = tmp_path / name
        library.write_bytes(name.encode())
        lines.append(f"1000-2000 r-xp 0000 00:01 1 {library}\n")
    maps.write_text("".join(lines))
    return maps


def _runtime(maps: Path) -> tuple[dict[str, Any] | None, str]:
    return cast(tuple[dict[str, Any] | None, str], CACHE.runtime_identity(maps_path=maps))


def _run(cache_dir: Path, runtime: dict[str, Any] | None, calculate: _Calculator) -> dict[str, Any]:
    run = cast(Callable[..., dict[str, Any]], CACHE.cached_run)
    return run(cache_dir, {"source": "same-source", "cells": 200}, runtime, calculate, label="F5")


def test_identical_verified_runtime_reuses_calculation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    calculate = _Calculator()
    first = _run(tmp_path / "cache", runtime, calculate)
    same_runtime, _ = _runtime(maps)
    assert _run(tmp_path / "cache", same_runtime, calculate) == first
    assert calculate.calls == 1


def test_image_package_native_bytes_and_solver_settings_each_invalidate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    calculate = _Calculator()
    cache_dir = tmp_path / "cache"
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    monkeypatch.setenv("NUCLEIDE_FESTIM_IMAGE_ID", "sha256:" + "2" * 64)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 2
    changed = [_Distribution(name, "2.0" if name == "petsc4py" else "1.0") for name in _PACKAGES]
    monkeypatch.setattr(CACHE.importlib.metadata, "distributions", lambda: changed)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 3
    (tmp_path / "libpetsc.so.3").write_bytes(b"different native solver build, same version")
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 4
    monkeypatch.setenv("PETSC_OPTIONS", "-ksp_type gmres")
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 5


def test_missing_image_marker_bypasses_existing_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    cache_dir = tmp_path / "cache"
    calculate = _Calculator()
    _run(cache_dir, runtime, calculate)
    original_files = {path.name: path.read_bytes() for path in cache_dir.iterdir()}
    monkeypatch.delenv("NUCLEIDE_FESTIM_IMAGE_ID")
    unverified, reason = _runtime(maps)
    assert unverified is None and "IMAGE_ID" in reason
    assert _run(cache_dir, unverified, calculate) == {"computed": 2}
    assert _run(cache_dir, unverified, calculate) == {"computed": 3}
    assert {path.name: path.read_bytes() for path in cache_dir.iterdir()} == original_files


def test_unreadable_native_library_disables_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    calculate = _Calculator()
    cache_dir = tmp_path / "cache"
    _run(cache_dir, runtime, calculate)

    def denied(path: Path) -> str:
        raise PermissionError(f"cannot read {path}")

    monkeypatch.setattr(CACHE, "_file_digest", denied)
    unverified, reason = _runtime(maps)
    assert unverified is None and "cannot read" in reason
    _run(cache_dir, unverified, calculate)
    assert calculate.calls == 2


def test_missing_critical_native_or_package_disables_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    calculate = _Calculator()
    cache_dir = tmp_path / "cache"
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    _run(cache_dir, runtime, calculate)
    originals = {path.name: path.read_bytes() for path in cache_dir.iterdir()}
    maps.write_text(
        "".join(line for line in maps.read_text().splitlines(True) if "libpetsc" not in line)
    )
    runtime, reason = _runtime(maps)
    assert runtime is None and "PETSc" in reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 2
    maps = _environment(tmp_path, monkeypatch)
    packages = [_Distribution(name) for name in _PACKAGES if name != "petsc4py"]
    monkeypatch.setattr(CACHE.importlib.metadata, "distributions", lambda: packages)
    runtime, reason = _runtime(maps)
    assert runtime is None and "petsc4py" in reason
    _run(cache_dir, runtime, calculate)
    assert calculate.calls == 3
    assert {path.name: path.read_bytes() for path in cache_dir.iterdir()} == originals


def test_generated_jit_forms_do_not_change_runtime_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maps = _environment(tmp_path, monkeypatch)
    runtime, reason = _runtime(maps)
    assert runtime is not None, reason
    with maps.open("a") as stream:
        stream.write(f"2000-3000 r-xp 0000 00:01 2 {tmp_path / 'libffcx_newform.so'} (deleted)\n")
    after_jit, _ = _runtime(maps)
    assert after_jit == runtime
    with maps.open("a") as stream:
        stream.write(f"3000-4000 r-xp 0000 00:01 3 {tmp_path / 'other_native.so'} (deleted)\n")
    unverifiable, reason = _runtime(maps)
    assert unverifiable is None and "deleted" in reason
