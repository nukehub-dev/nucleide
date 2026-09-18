"""Execute every runnable ```python block in docs/tutorials/python/*.md.

Each page's blocks run in document order inside one shared namespace, with
the repository root as the working directory (tutorials read committed
fixtures by repo-relative paths). A block preceded (whitespace only) by
``<!-- code-test: skip -->`` is skipped — reserved for intentionally
non-runnable snippets such as placeholder paths. This keeps the tutorials
doubling as verified examples: any stale snippet fails CI instead of
rotting silently.
"""

import re
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parent.parent
PAGES = sorted((ROOT / "docs" / "tutorials" / "python").glob("*.md"))

FENCE = re.compile(r"```python[^\S\n]*\n(.*?)```", re.DOTALL)
SKIP = re.compile(r"<!--\s*code-test:\s*skip\s*-->")


def _blocks(page: Path) -> list[tuple[int, str, bool]]:
    """(index, code, runnable) triples for one tutorial page."""
    # Explicit UTF-8: tutorial prose carries non-ASCII symbols (⋅, ⟨⟩,
    # subscripts) and Windows defaults to a locale codec that chokes on them.
    text = page.read_text(encoding="utf-8")
    out: list[tuple[int, str, bool]] = []
    for i, match in enumerate(FENCE.finditer(text)):
        gap = text[: match.start()]
        stripped = gap.rstrip()
        skip = False
        skip_match = None
        for probe in SKIP.finditer(stripped):
            skip_match = probe
        if skip_match is not None and not stripped[skip_match.end() :].strip():
            skip = True
        out.append((i, match.group(1), not skip))
    return out


def _page_id(page: Path) -> str:
    return page.stem


@pytest.mark.parametrize("page", PAGES, ids=_page_id)  # type: ignore[untyped-decorator]
def test_tutorial_code_runs(page: Path) -> None:
    blocks = _blocks(page)
    if not blocks:
        pytest.skip(f"{page.name} has no ```python blocks")
    namespace: dict[str, Any] = {}
    ran = 0
    for i, code, runnable in blocks:
        if not runnable:
            continue
        exec(compile(code, f"{page.name}#{i}", "exec"), namespace)
        ran += 1
    assert ran > 0, f"{page.name} has only skipped ```python blocks"
