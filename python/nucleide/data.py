"""Download Nucleide data files pinned to a release tag, branch, or commit.

The wheel ships no data files; this module fetches them from the GitHub
repository on demand (e.g. the Materials Compendium JSON, sample depletion
chains), defaulting to the tag matching the installed version.
"""

import urllib.error
import urllib.request
from pathlib import Path

from nucleide._internal import version

_RAW_BASE = "https://raw.githubusercontent.com/nukehub-dev/nucleide"

COMPENDIUM_PATH = "fixtures/data/MaterialsCompendium.json"

__all__ = [
    "COMPENDIUM_PATH",
    "default_ref",
    "fetch",
    "fetch_compendium",
]


def default_ref() -> str:
    """Git ref matching the installed Nucleide version (e.g. ``"v0.1.0"``)."""
    return f"v{version()}"


def fetch(path: str, *, ref: str | None = None, dest: str | Path = ".") -> str:
    """Download a repository data file and return its local path as a string.

    ``path`` is repo-relative (e.g. ``"fixtures/depletion/chain_simple.xml"``);
    the file keeps its basename under ``dest``. ``ref`` is a tag, branch, or
    commit and defaults to the installed version's tag — pass ``"main"`` or a
    commit SHA for development installs whose tag does not exist yet.
    """
    ref = ref or default_ref()
    url = f"{_RAW_BASE}/{ref}/{path}"
    out = Path(dest) / Path(path).name
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url) as response:
            out.write_bytes(response.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"failed to download {url}: HTTP {exc.code}. "
            "For a development install, pass ref='main' or a commit SHA."
        ) from exc
    return str(out)


def fetch_compendium(*, ref: str | None = None, dest: str | Path = ".") -> str:
    """Download the DOE/PNNL Materials Compendium JSON and return its path."""
    return fetch(COMPENDIUM_PATH, ref=ref, dest=dest)
