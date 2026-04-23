from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


def test_manual_jupyter_notebook_executes():
    try:
        import nbclient
        import nbformat
    except ImportError:
        pytest.fail(
            "Notebook execution test dependencies are missing. "
            "Install them with `pip install .[test]`.",
            pytrace=False,
        )

    notebook_path = Path(__file__).with_name("manual_jupyter_test.ipynb")
    notebook = nbformat.read(notebook_path, as_version=4)

    client = nbclient.NotebookClient(
        notebook,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_path.parent.parent)}},
        startup_timeout=180,
        timeout=900,
    )
    executed = client.execute()

    error_outputs = [
        output
        for cell in executed.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    assert error_outputs == []
