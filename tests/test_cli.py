"""Unit tests for command-line behavior."""

from typer.testing import CliRunner

from plannotate import __version__, _package_data
from plannotate import main as main_module
from plannotate.main import app


def test_cli_help_without_databases():
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "setupdb" in result.stdout
    assert "batch" in result.stdout
    assert "streamlit" in result.stdout


def test_cli_version(monkeypatch):
    monkeypatch.setattr(
        main_module._package_data,
        "get_database_manifest",
        lambda: {
            "bundle": "plannotate-databases-v2",
            "build_date": "2026-06-20",
            "databases": {"Rfam": {"version": "release 15.1"}},
        },
    )

    result = CliRunner().invoke(app, ["--version"])

    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "database path" in result.stdout
    assert "plannotate-databases-v2" in result.stdout
    assert "2026-06-20" in result.stdout
    # per-source versions are listed
    assert "Rfam" in result.stdout
    assert "release 15.1" in result.stdout


def test_cli_version_without_databases(monkeypatch):
    def _missing():
        raise FileNotFoundError

    monkeypatch.setattr(main_module._package_data, "get_database_manifest", _missing)

    result = CliRunner().invoke(app, ["--version"])

    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "not installed" in result.stdout


def test_streamlit_requires_the_server_extra(monkeypatch, caplog):
    monkeypatch.setattr(main_module, "_streamlit_available", lambda: False)

    result = CliRunner().invoke(app, ["streamlit"])

    assert result.exit_code == 1
    assert "plannotate[server]" in caplog.text


def test_streamlit_requires_databases(monkeypatch, caplog):
    monkeypatch.setattr(main_module, "_streamlit_available", lambda: True)
    monkeypatch.setattr(_package_data, "databases_exist", lambda: False)

    result = CliRunner().invoke(app, ["streamlit"])

    assert result.exit_code == 1
    assert "setupdb" in caplog.text


def test_batch_requires_input():
    result = CliRunner().invoke(app, ["batch"])

    assert result.exit_code == 2


def test_setupdb_force_reinstalls_existing_databases(monkeypatch):
    downloads = []
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    monkeypatch.setattr(
        _package_data, "download_databases", lambda: downloads.append(True)
    )

    result = CliRunner().invoke(app, ["setupdb", "--force"])

    assert result.exit_code == 0
    assert downloads == [True]


def test_databases_prints_installed_manifest(monkeypatch):
    monkeypatch.setattr(
        _package_data,
        "get_database_manifest",
        lambda: {"schema_version": 1, "databases": {"Rfam": {"version": "15.1"}}},
    )

    result = CliRunner().invoke(app, ["databases"])

    assert result.exit_code == 0
    assert '"Rfam"' in result.stdout
    assert '"15.1"' in result.stdout
