"""Optional synthetic ACTINV 1.7.2 export/readback probe; no activation solve.

Run after rebuilding Nucleide with ``maturin develop``. Supply the verified
Linux release binary explicitly with ``--actinv /path/to/actinv``. No program
or evaluated data is downloaded. Run directories and JSON receipts are kept
under ``--output-dir`` (default: ``validation/.cache/actinv-handoff``).

The pinned exporter is a leaf CLI process: subprocess.run's 30-second timeout
kills and reaps that owned process. The caller must apply the documented
memory/CPU/task limits; this script does not claim to enforce them itself.
Receipts record readable cgroup v2 controls and explain unavailable values.
This exchange check is separate from the independent physics comparisons
auto-discovered by run_all.sh.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ACTINV_RELEASE = "1.7.2"
ACTINV_COMMIT = "d3dc4e463346c1a6d3dc426432048022932e14c5"
ACTINV_SHA256 = "7d52a7b234f3125cc48ce222b32a1c22235acb756d89adbff386260387cf28d1"
ROOT = Path(__file__).resolve().parent.parent


def resource_limits() -> dict[str, Any]:
    """Record current cgroup v2 controls without changing or assuming limits."""
    result: dict[str, Any] = {
        "kind": "cgroup_v2",
        "status": "unavailable",
        "enforced_by_probe": False,
        "limits": {},
        "note": "Current cgroup controls; inherited ancestor limits may be tighter.",
    }
    try:
        membership = Path("/proc/self/cgroup").read_text(encoding="utf-8")
    except OSError as error:
        result["reason"] = f"Cannot read /proc/self/cgroup: {error}"
        return result
    result["proc_self_cgroup"] = membership
    unified = [line[3:] for line in membership.splitlines() if line.startswith("0::")]
    if len(unified) != 1:
        result["reason"] = "No unique unified cgroup v2 membership in /proc/self/cgroup"
        return result
    result["cgroup_path"] = unified[0]
    mount = Path("/sys/fs/cgroup").resolve()
    directory = (mount / unified[0].lstrip("/")).resolve()
    if not directory.is_relative_to(mount):
        result["reason"] = "Cgroup membership does not resolve within /sys/fs/cgroup"
        return result
    result["directory"] = str(directory)
    unavailable = {}
    for name in ("memory.max", "memory.swap.max", "pids.max", "cpu.max"):
        try:
            result["limits"][name] = (directory / name).read_text(encoding="utf-8").strip()
        except OSError as error:
            unavailable[name] = str(error)
    if unavailable:
        result["unavailable"] = unavailable
        result["status"] = "partial" if result["limits"] else "unavailable"
    else:
        result["status"] = "recorded"
    return result


def sha256(path: Path) -> str:
    """Hash a file without loading the full binary into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    """Keep acceptance checks active even when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def check_close(actual: float, expected: float, label: str) -> None:
    require(
        math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=0.0),
        f"{label}: expected {expected}, got {actual}",
    )


def probe(binary_path: Path, run_dir: Path, receipt: dict[str, Any]) -> None:
    """Export a known integrated source and check real Nucleide readers/tags."""
    binary = binary_path.expanduser().resolve(strict=True)
    require(binary.is_file(), f"ACTINV path is not a file: {binary}")
    actual_hash = sha256(binary)
    receipt.update({"actinv_path": str(binary), "actinv_binary_sha256": actual_hash})
    require(actual_hash == ACTINV_SHA256, "ACTINV binary SHA-256 does not match pinned 1.7.2")

    # Use the installed package; do not silently prepend another source tree.
    import nucleide
    from nucleide import _internal, alara, r2s

    extension = Path(_internal.__file__).resolve()
    try:
        distribution_version = importlib.metadata.version("nucleide")
    except importlib.metadata.PackageNotFoundError:
        distribution_version = "not installed as a distribution"
    receipt["nucleide"] = {
        "runtime_version": _internal.version(),
        "distribution_version": distribution_version,
        "package_path": str(Path(nucleide.__file__).resolve()),
        "extension_path": str(extension),
        "extension_sha256": sha256(extension),
        "checkout_sources_sha256": {
            name: sha256(ROOT / name)
            for name in (
                "Cargo.toml",
                "bindings/python/src/lib.rs",
                "crates/alara-io/src/photon.rs",
                "crates/r2s/src/tags.rs",
            )
        },
    }
    records = [
        {"record": "header", "schema": "actinv-r2s-source-1", "step": 2},
        {
            "record": "cell",
            "id": "zone1",
            "bounds_cm": [[0, 2], [0, 1], [0, 1]],
            "volume_cm3": 2,
            "photons_s": 12,
            "step_t_s": 3900,
            "groups": [
                {"centroid_eV": 1e5, "photons_s": 4},
                {"centroid_eV": 1e6, "photons_s": 8},
            ],
        },
        {"record": "footer", "cell_count": 1},
    ]
    source = run_dir / "synthetic-source.ndjson"
    source.write_text("\n".join(json.dumps(row) for row in records) + "\n", encoding="utf-8")
    output = run_dir / "alara"
    command = [
        str(binary),
        "export-source",
        "alara",
        str(source),
        str(output),
        "--shutdown-t-s",
        "300",
    ]
    receipt.update(
        {
            "command": command,
            "timeout_s": 30,
            "process_model": "pinned leaf exporter; timeout kills and reaps the child",
            "source_sha256": sha256(source),
        }
    )
    child = subprocess.run(command, text=True, capture_output=True, timeout=30, check=False)
    receipt.update({"returncode": child.returncode, "stdout": child.stdout, "stderr": child.stderr})
    require(child.returncode == 0, f"ACTINV export failed: {child.stderr.strip()}")

    index_path = output / "actinv-alara-index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    require(index["input_sha256"] == receipt["source_sha256"], "export index input hash mismatch")
    require(index["units"] == "photons/s/cm3", "unexpected exported density units")
    require(index["group_centroids_eV"] == [1e5, 1e6], "unexpected centroid labels/order")
    require(index["step"] == 2, "unexpected exported step")
    check_close(index["step_t_s"], 3900.0, "absolute source time")
    check_close(index["shutdown_t_s"], 300.0, "shutdown time")
    check_close(index["cooling_s"], 3600.0, "cooling time")
    require(len(index["cells"]) == 1, "expected exactly one exported cell")
    cell = index["cells"][0]
    require(cell["id"] == "zone1", "unexpected exported zone")
    check_close(cell["volume_cm3"], 2.0, "zone volume")
    check_close(cell["photons_s"], 12.0, "declared integrated strength")
    photon_path = (output / cell["file"]).resolve(strict=True)
    require(
        photon_path.is_relative_to(output.resolve()), "export index file escapes output directory"
    )
    require(sha256(photon_path) == cell["sha256"], "exported photon file hash mismatch")
    photon_text = photon_path.read_text(encoding="utf-8")
    # A different cooling row must not contaminate the actual export selection.
    selection_text = "TOTAL 600 s 8 4\n" + photon_text
    groups = r2s.photon_group_sums(selection_text, ["TOTAL"], 3600.0)
    require(groups["sums"] == [2.0, 4.0], "unexpected selected photon-group densities")
    check_close(groups["total"], 6.0, "selected density")
    check_close(alara.alara_photon_total_strength(photon_text), 6.0, "photon reader density")
    integrated = groups["total"] * cell["volume_cm3"]
    check_close(integrated, 12.0, "density times zone volume")
    tags = r2s.tag_zone_strength([integrated], [0, 0], split=True, voxel_volumes_cm3=[0.5, 1.5])
    require(tags["source_strength"] == [3.0, 9.0], "unexpected integrated voxel strengths")
    check_close(tags["total"], 12.0, "tagged integrated strength")
    receipt.update(
        {
            "index": index,
            "index_sha256": sha256(index_path),
            "photon_file_sha256": sha256(photon_path),
            "density_groups": groups["sums"],
            "density_total_photons_s_cm3": groups["total"],
            "zone_volume_cm3": cell["volume_cm3"],
            "reconstructed_photons_s": integrated,
            "voxel_volumes_cm3": [0.5, 1.5],
            "voxel_strengths_photons_s": tags["source_strength"],
            "tagged_photons_s": tags["total"],
            "cooling_s": index["cooling_s"],
            "photon_centroids_are_bin_edges": False,
        }
    )


def main() -> int:
    """Run an explicitly configured probe and retain a success/failure receipt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actinv", required=True, type=Path, help="pinned 1.7.2 Linux binary")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "validation/.cache/actinv-handoff"
    )
    args = parser.parse_args()
    try:
        output_root = args.output_dir.expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        run_dir = Path(tempfile.mkdtemp(prefix="run-", dir=output_root))
    except OSError as error:
        print(f"FAIL: cannot create ACTINV probe output: {error}", file=sys.stderr)
        return 1
    receipt: dict[str, Any] = {
        "kind": "synthetic_export_readback_only",
        "status": "FAIL",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "actinv_release": ACTINV_RELEASE,
        "actinv_source_commit": ACTINV_COMMIT,
        "expected_actinv_binary_sha256": ACTINV_SHA256,
        "activation_solve_performed": False,
        "evaluated_data_downloaded": False,
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "run_directory": str(run_dir),
        "resource_limits": resource_limits(),
    }
    try:
        probe(args.actinv, run_dir, receipt)
        receipt["status"] = "PASS"
    except Exception as error:
        receipt["error"] = f"{type(error).__name__}: {error}"
        if isinstance(error, subprocess.TimeoutExpired):
            receipt["timeout_s"] = error.timeout
            for name, value in (("stdout", error.stdout), ("stderr", error.stderr)):
                receipt[name] = (
                    value.decode(errors="replace") if isinstance(value, bytes) else value
                )
    serialized = json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    receipt_path = run_dir / "receipt.json"
    receipt_path.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
    print(f"Receipt: {receipt_path}", file=sys.stderr)
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
