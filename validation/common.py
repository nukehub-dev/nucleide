"""Shared helpers and report writer for the validation harness."""

from __future__ import annotations

import importlib.metadata
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = Path(__file__).resolve().parent / "results"


def _pyne_version() -> str:
    """Return the true PyNE package version.

    ``pyne.__version__`` carries stale upstream metadata in the conda package,
    so prefer ``conda list pyne --json``. Fall back to importlib metadata when
    conda is unavailable.
    """
    try:
        proc = subprocess.run(
            ["conda", "list", "pyne", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        entries = json.loads(proc.stdout)
        if isinstance(entries, list) and entries:
            version = entries[0].get("version")
            if version:
                return str(version)
    except Exception:
        pass
    return importlib.metadata.version("pyne")


def environment() -> dict[str, str]:
    """Return environment metadata for the validation run."""
    return {
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "nucleide": importlib.metadata.version("nucleide"),
        "pyne": _pyne_version(),
        "openmc": importlib.metadata.version("openmc"),
    }


def rel_diff(a: float, b: float) -> float:
    """Relative difference robust to near-zero values."""
    denom = max(abs(a), abs(b), 1.0e-30)
    return abs(a - b) / denom


def abs_diff(a: float, b: float) -> float:
    """Absolute difference."""
    return abs(a - b)


def fmt(value: float) -> str:
    """Format a numeric value for table cells."""
    return f"{value:.6e}"


class Report:
    """Accumulates sections and writes a JSON report plus a plain-text summary."""

    def __init__(self, name: str, title: str) -> None:
        self.name = name
        self.title = title
        self.sections: list[dict[str, Any]] = []

    def heading(self, text: str, level: int = 3) -> None:
        """Add a heading section."""
        self.sections.append({"kind": "heading", "level": level, "text": text})

    def prose(self, text: str) -> None:
        """Add a prose paragraph."""
        self.sections.append({"kind": "prose", "text": text})

    def table(self, headers: list[str], rows: list[list[str]]) -> None:
        """Add a table section."""
        self.sections.append({"kind": "table", "headers": headers, "rows": rows})

    def _summary(self) -> str:
        """Return a human-readable plain-text summary of the report."""
        lines: list[str] = [f"=== {self.title} ==="]
        for section in self.sections:
            kind = section["kind"]
            if kind == "heading":
                level = section["level"]
                marker = "=" if level <= 3 else "-"
                lines.append(f"{marker * 3} {section['text']} {marker * 3}")
            elif kind == "prose":
                lines.append(section["text"])
            elif kind == "table":
                headers = section["headers"]
                rows = section["rows"]
                lines.append(" | ".join(headers))
                for row in rows:
                    lines.append(" | ".join(row))
        return "\n".join(lines)

    def emit(self) -> None:
        """Write the JSON report and print a plain-text summary."""
        data = {
            "schema": 1,
            "name": self.name,
            "title": self.title,
            "environment": environment(),
            "sections": self.sections,
        }
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        path = RESULTS_DIR / f"{self.name}.json"
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(f"Report written to {path}")
        print(self._summary())
