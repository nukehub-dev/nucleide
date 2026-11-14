"""Conservative cache identity for the optional FESTIM validation environment."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import re
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

_REQUIRED_PACKAGES = {
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
}
_ENV_PREFIXES = (
    "PETSC",
    "DOLFINX",
    "FFCX",
    "UFL",
    "BASIX",
    "OMP",
    "OPENBLAS",
    "RAYON",
    "MKL",
    "BLIS",
    "VECLIB",
    "MPI",
    "OMPI",
    "MPICH",
    "I_MPI",
    "UCX",
    "FI_",
)
_ENV_NAMES = {
    "LD_LIBRARY_PATH",
    "LD_PRELOAD",
    "PYTHONPATH",
    "PYTHONHOME",
    "PATH",
    "CC",
    "CXX",
    "CFLAGS",
    "CXXFLAGS",
    "LDFLAGS",
}


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def runtime_identity(
    *, maps_path: Path = Path("/proc/self/maps")
) -> tuple[dict[str, Any] | None, str]:
    """Collect verifiable runtime identity, or disable persistent caching.

    The runner supplies its actual inspected immutable image ID. Import
    FESTIM before collection, and do not mutate the native environment
    during the probe. Generated libffcx_* JIT forms are excluded: the caller's
    source/parameter identity covers forms and the image identifies compilers.
    """
    image_id = os.environ.get("NUCLEIDE_FESTIM_IMAGE_ID", "")
    if re.fullmatch(r"sha256:[0-9a-f]{64}", image_id) is None:
        return None, "actual inspected NUCLEIDE_FESTIM_IMAGE_ID is unavailable or invalid"
    try:
        packages = []
        for distribution in importlib.metadata.distributions():
            name = distribution.metadata.get("Name")
            version = distribution.version
            if not name or not version:
                return None, "installed distribution has incomplete name/version metadata"
            packages.append([re.sub(r"[-_.]+", "-", name).lower(), version])
        packages.sort()
        missing = _REQUIRED_PACKAGES - {entry[0] for entry in packages}
        if missing:
            return (
                None,
                f"critical oracle distribution metadata missing: {', '.join(sorted(missing))}",
            )
        paths: set[Path] = set()
        for line in maps_path.read_text().splitlines():
            fields = line.split(maxsplit=5)
            if len(fields) != 6 or not Path(fields[5].removesuffix(" (deleted)")).is_absolute():
                continue
            mapped = fields[5]
            if re.search(r"\.so(?:\.[^ /]+)*(?: \(deleted\))?$", mapped) is None:
                continue
            path = Path(mapped.removesuffix(" (deleted)"))
            if path.name.startswith("libffcx_"):
                continue
            if mapped.endswith(" (deleted)"):
                return None, f"loaded native library has been deleted: {mapped}"
            paths.add(path.resolve(strict=True))
        names = [path.name.lower() for path in paths]
        if not any("libpetsc" in name for name in names) or not any(
            "libdolfinx" in name for name in names
        ):
            return None, "loaded PETSc and DOLFINx native libraries could not both be identified"
        native = [[str(path), _file_digest(path)] for path in sorted(paths)]
    except Exception as exc:  # noqa: BLE001 — unverified identity disables reuse
        return None, f"runtime identity could not be verified: {exc}"
    settings = {
        name: value
        for name, value in sorted(os.environ.items())
        if name in _ENV_NAMES or name.startswith(_ENV_PREFIXES)
    }
    return {
        "image_id": image_id,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version,
        "executable": sys.executable,
        "packages": packages,
        "settings": settings,
        "native_libraries": native,
    }, "runtime image, distributions, settings and loaded native libraries verified"


def cached_run(
    cache_dir: Path,
    identity: dict[str, Any],
    runtime: dict[str, Any] | None,
    run: Callable[[], dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    """Reuse reduced calculations only with verified runtime provenance."""
    if runtime is None:
        print(f"{label}: cache disabled: unverified runtime; starting", flush=True)
        return run()
    combined = {
        **identity,
        "cache_helper_sha256": _file_digest(Path(__file__)),
        "runtime": runtime,
    }
    key = hashlib.sha256(json.dumps(combined, sort_keys=True).encode()).hexdigest()
    path = cache_dir / f"{key}.json"
    if path.is_file():
        saved = json.loads(path.read_text())
        if saved["identity"] != combined:
            raise ValueError(f"{label}: cache identity mismatch: {path}")
        print(f"{label}: cached {path.name}", flush=True)
        result: dict[str, Any] = saved["result"]
        return result
    print(f"{label}: starting", flush=True)
    result = run()
    cache_dir.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(".pending")
    pending.write_text(json.dumps({"identity": combined, "result": result}, allow_nan=False) + "\n")
    pending.replace(path)
    print(f"{label}: completed, cached {path.name}", flush=True)
    return result
