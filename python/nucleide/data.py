"""Download Nucleide data files pinned to a release tag, branch, or commit.

The wheel ships no data files; this module fetches them from the GitHub
repository on demand (e.g. the Materials Compendium JSON, sample depletion
chains), defaulting to the tag matching the installed version.

It also fetches third-party datasets that cannot be vendored, pinned by
content hash rather than by git ref (currently the EPA FGR 15
external-dosimetry coefficient zip: users download directly from EPA, and
the pinned SHA-256 guards against silent revisions).
"""

import hashlib
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

from nucleide._internal import version

_RAW_BASE = "https://raw.githubusercontent.com/nukehub-dev/nucleide"

COMPENDIUM_PATH = "fixtures/data/MaterialsCompendium.json"

#: EPA Federal Guidance Report No. 15 coefficient zip (EPA 402-R-25-001,
#: July 2025; data file dated 2025-05-28). US government work, fetched
#: directly from EPA at runtime — nothing from this file is vendored.
FGR15_URL = "https://www.epa.gov/system/files/other-files/2025-07/fgr15_data_2025_05_28.zip"

#: SHA-256 of the published FGR 15 zip; a mismatch means EPA revised the
#: file after this Nucleide version was pinned.
FGR15_SHA256 = "71314b3f1d73c197da8b589e74c450f61befce48c66ace6ac1f6e0597560bd91"

#: IAEA IRDFF-II dosimetry cross sections, pre-grouped in the SAND-II
#: 725-group structure (release of January 2020, files updated 2020-12-09;
#: Trkov et al., Nucl. Data Sheets 163 (2020) 1). IAEA copyright —
#: permission-based reproduction — so nothing is vendored; users download
#: directly from the IAEA at runtime and the pinned SHA-256 guards against
#: silent revisions.
IRDFF_URL = "https://www-nds.iaea.org/IRDFF/IRDFF-II_g725.zip"

#: SHA-256 of the published IRDFF-II_g725.zip; a mismatch means the IAEA
#: revised the file after this Nucleide version was pinned.
IRDFF_SHA256 = "6ec2b33c0f67bed46d46be062a24ccedaa5ffea9bbba919958da4b1349f48c85"

#: Name of the data member inside the IRDFF-II_g725 zip.
IRDFF_MEMBER = "IRDFF-II.g725"

__all__ = [
    "COMPENDIUM_PATH",
    "FGR15_SHA256",
    "FGR15_URL",
    "IRDFF_MEMBER",
    "IRDFF_SHA256",
    "IRDFF_URL",
    "default_ref",
    "fetch",
    "fetch_compendium",
    "fetch_fgr15",
    "fetch_irdff",
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


def _default_cache_dir() -> Path:
    """Per-user platform cache shared by runtime-fetched datasets.

    Windows ``%LOCALAPPDATA%\\nucleide\\Cache``, macOS
    ``~/Library/Caches/nucleide``, other POSIX ``$XDG_CACHE_HOME/nucleide``
    or ``~/.cache/nucleide``.
    """
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA")
        root = Path(base) if base else Path.home() / "AppData" / "Local"
        return root / "nucleide" / "Cache"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "nucleide"
    base = os.environ.get("XDG_CACHE_HOME")
    return (Path(base) if base else Path.home() / ".cache") / "nucleide"


def _sha256(path: Path) -> str:
    """Hex SHA-256 of a file, chunked."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_fgr15(*, dest: str | Path | None = None, url: str | None = None) -> str:
    """Download the EPA FGR 15 coefficient zip (hash-pinned) and return its path.

    FGR 15 (EPA 402-R-25-001, July 2025) external-dosimetry tables are not
    vendored: this fetches the official EPA zip — by default into the
    per-user platform cache (``%LOCALAPPDATA%\\nucleide\\Cache`` on Windows,
    ``~/Library/Caches/nucleide`` on macOS, ``$XDG_CACHE_HOME/nucleide`` or
    ``~/.cache/nucleide`` elsewhere) — and verifies its SHA-256 against
    the pinned ``FGR15_SHA256``. A verified cached file is reused as-is (no
    re-download); a hash mismatch raises loudly naming both digests (EPA may
    have revised the file), and a download failure with no usable cache
    raises a clear error. Pass ``url=`` to fetch from a mirror or a local
    ``file://`` copy (tests); the pinned hash is always enforced.
    """
    url = url or FGR15_URL
    dest_dir = Path(dest) if dest is not None else _default_cache_dir()
    out = dest_dir / Path(url).name
    out.parent.mkdir(parents=True, exist_ok=True)

    def verify(path: Path) -> str:
        actual = _sha256(path)
        if actual != FGR15_SHA256:
            raise RuntimeError(
                f"sha256 mismatch for {path.name}: expected {FGR15_SHA256}, got {actual}. "
                "EPA may have revised the FGR 15 data file after this Nucleide version was "
                "pinned; update Nucleide (or pass url= pointing at a verified copy of the "
                "2025-05-28 release)."
            )
        return str(path)

    if out.exists():
        return verify(out)
    try:
        with urllib.request.urlopen(url) as response:
            out.write_bytes(response.read())
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError(
            f"failed to download the FGR 15 data zip from {url!r}: {exc}. "
            f"No usable cached copy exists at {out}."
        ) from exc
    try:
        return verify(out)
    except RuntimeError:
        out.unlink(missing_ok=True)
        raise


def fetch_irdff(*, dest: str | Path | None = None, url: str | None = None) -> str:
    """Download the IAEA IRDFF-II g725 zip (hash-pinned) and return its path.

    IRDFF-II (Trkov et al., Nucl. Data Sheets 163 (2020) 1) dosimetry cross
    sections are IAEA copyright and never vendored: this fetches the official
    IAEA zip — by default into the per-user platform cache (same layout as
    :func:`fetch_fgr15`) — and verifies its SHA-256 against the pinned
    ``IRDFF_SHA256``. A verified cached file is reused as-is (no
    re-download); a hash mismatch raises loudly naming both digests, and a
    download failure with no usable cache raises a clear error. Pass ``url=``
    to fetch from a mirror or a local ``file://`` copy (tests); the pinned
    hash is always enforced. Parse the ``IRDFF_MEMBER`` member text with
    ``nucleide._internal.parse_irdff_g725`` into caller-ready response rows
    for the ``unfold`` iterators.
    """
    url = url or IRDFF_URL
    dest_dir = Path(dest) if dest is not None else _default_cache_dir()
    out = dest_dir / Path(url).name
    out.parent.mkdir(parents=True, exist_ok=True)

    def verify(path: Path) -> str:
        actual = _sha256(path)
        if actual != IRDFF_SHA256:
            raise RuntimeError(
                f"sha256 mismatch for {path.name}: expected {IRDFF_SHA256}, got {actual}. "
                "IAEA may have revised the IRDFF-II g725 file after this Nucleide version was "
                "pinned; update Nucleide (or pass url= pointing at a verified copy of the "
                "January 2020 release)."
            )
        return str(path)

    if out.exists():
        return verify(out)
    try:
        with urllib.request.urlopen(url) as response:
            out.write_bytes(response.read())
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError(
            f"failed to download the IRDFF-II data zip from {url!r}: {exc}. "
            f"No usable cached copy exists at {out}."
        ) from exc
    try:
        return verify(out)
    except RuntimeError:
        out.unlink(missing_ok=True)
        raise
