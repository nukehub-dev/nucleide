"""Tests for `nucleide.data` (network access is mocked)."""

import io
import urllib.error
import urllib.request
from pathlib import Path

import pytest

import nucleide
from nucleide import data


class _FakeResponse:
    def __enter__(self) -> io.BytesIO:
        return io.BytesIO(b"{}")

    def __exit__(self, *args: object) -> None:
        return None


def test_default_ref_matches_version() -> None:
    assert data.default_ref() == f"v{nucleide.__version__}"


def test_fetch_downloads_and_returns_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    calls: dict[str, str] = {}

    def fake_urlopen(url: str) -> _FakeResponse:
        calls["url"] = url
        return _FakeResponse()

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    out = data.fetch("fixtures/data/MaterialsCompendium.json", ref="v9.9.9", dest=tmp_path)
    assert out == str(tmp_path / "MaterialsCompendium.json")
    assert Path(out).read_bytes() == b"{}"
    assert calls["url"] == (
        "https://raw.githubusercontent.com/nukehub-dev/nucleide/"
        "v9.9.9/fixtures/data/MaterialsCompendium.json"
    )


def test_fetch_defaults_to_version_tag(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    calls: dict[str, str] = {}

    def fake_urlopen(url: str) -> _FakeResponse:
        calls["url"] = url
        return _FakeResponse()

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    data.fetch_compendium(dest=tmp_path)
    assert f"/v{nucleide.__version__}/" in calls["url"]


def test_fetch_http_error_suggests_ref_override(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def fake_urlopen(url: str) -> _FakeResponse:
        raise urllib.error.HTTPError(url, 404, "Not Found", None, None)  # type: ignore[arg-type]

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(RuntimeError, match="ref='main'"):
        data.fetch("fixtures/data/MaterialsCompendium.json", ref="v0.0.0", dest=tmp_path)
