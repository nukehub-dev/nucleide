"""Execute the notebooks/ workflows end to end, and keep them clean.

Every ``notebooks/*.ipynb`` runs top to bottom in a fresh kernel with the
repository root as the working directory (workflows read committed fixtures
by repo-relative paths and use fixed seeds, so execution is deterministic).
Committed notebooks must carry no outputs and no execution counts — results
are recomputed by the reader, never vendored in the file. This keeps the
diffs reviewable and the workflows honest: anything stale fails CI instead
of rotting silently.
"""

from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS = sorted((ROOT / "notebooks").glob("*.ipynb"))


def _notebook_ids() -> list[str]:
    return [p.stem for p in NOTEBOOKS]


def _assert_stripped(path: Path) -> None:
    nb = nbformat.read(path, as_version=4)
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        assert not cell.get("outputs"), f"{path.name} cell {i} carries outputs"
        assert cell.get("execution_count") is None, (
            f"{path.name} cell {i} carries an execution count"
        )


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=_notebook_ids())  # type: ignore[untyped-decorator]
def test_notebooks_are_stripped(notebook: Path) -> None:
    assert NOTEBOOKS, "no notebooks found under notebooks/"
    _assert_stripped(notebook)


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=_notebook_ids())  # type: ignore[untyped-decorator]
def test_notebooks_execute(notebook: Path) -> None:
    nb = nbformat.read(notebook, as_version=4)
    client = NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    )
    client.execute()
