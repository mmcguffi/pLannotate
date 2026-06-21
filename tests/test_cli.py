from typer.testing import CliRunner

from plannotate import resources
from plannotate.main import app


def test_cli_help_without_databases():
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "setupdb" in result.stdout
    assert "batch" in result.stdout


def test_batch_requires_input():
    result = CliRunner().invoke(app, ["batch"])

    assert result.exit_code == 2


def test_setupdb_force_reinstalls_existing_databases(monkeypatch):
    downloads = []
    monkeypatch.setattr(resources, "databases_exist", lambda: True)
    monkeypatch.setattr(resources, "download_databases", lambda: downloads.append(True))

    result = CliRunner().invoke(app, ["setupdb", "--force"])

    assert result.exit_code == 0
    assert downloads == [True]


def test_databases_prints_installed_manifest(monkeypatch):
    monkeypatch.setattr(
        resources,
        "get_database_manifest",
        lambda: {"schema_version": 1, "databases": {"Rfam": {"version": "15.1"}}},
    )

    result = CliRunner().invoke(app, ["databases"])

    assert result.exit_code == 0
    assert '"Rfam"' in result.stdout
    assert '"15.1"' in result.stdout
