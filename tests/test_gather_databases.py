"""Tests for the public database-build API."""

from subprocess import CompletedProcess

import pytest

from plannotate import build_databases
from plannotate import gather_databases as gather_module


def test_build_databases_runs_packaged_workflow(monkeypatch, tmp_path):
    observed = {}

    def fake_run(command):
        observed["command"] = command
        return CompletedProcess(command, 0)

    monkeypatch.setattr(gather_module.shutil, "which", lambda _: "/bin/snakemake")
    monkeypatch.setattr(gather_module.subprocess, "run", fake_run)

    archive = build_databases(tmp_path / "database build", cores=3, dry_run=True)

    command = observed["command"]
    assert command[0] == "/bin/snakemake"
    assert command[command.index("--directory") + 1] == str(
        (tmp_path / "database build").resolve()
    )
    assert command[command.index("--cores") + 1] == "3"
    assert "--dry-run" in command
    assert (
        archive == (tmp_path / "database build").resolve() / gather_module.ARCHIVE_NAME
    )


def test_build_databases_requires_snakemake(monkeypatch, tmp_path):
    monkeypatch.setattr(gather_module.shutil, "which", lambda _: None)

    with pytest.raises(RuntimeError, match=r"pLannotate\[databases\]"):
        build_databases(tmp_path)


def test_build_databases_rejects_invalid_core_count(tmp_path):
    with pytest.raises(ValueError, match="cores must be at least 1"):
        build_databases(tmp_path, cores=0)
