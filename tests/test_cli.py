from typer.testing import CliRunner

from plannotate.main import app


def test_cli_help_without_databases():
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "setupdb" in result.stdout
    assert "batch" in result.stdout


def test_batch_requires_input():
    result = CliRunner().invoke(app, ["batch"])

    assert result.exit_code == 2
