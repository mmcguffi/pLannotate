"""Ensure every project Python module explains its purpose."""

import ast
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
PYTHON_FILES = [
    path
    for directory in ("plannotate", "tests", "tools")
    for path in sorted((PROJECT_ROOT / directory).rglob("*.py"))
]


@pytest.mark.parametrize(
    "path",
    PYTHON_FILES,
    ids=lambda path: str(path.relative_to(PROJECT_ROOT)),
)
def test_module_has_docstring(path: Path):
    module = ast.parse(path.read_text())

    assert ast.get_docstring(module), (
        f"{path.relative_to(PROJECT_ROOT)} has no docstring"
    )
