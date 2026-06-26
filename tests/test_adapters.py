"""Tests for external annotation-tool integrations."""

import logging
from pathlib import Path
from subprocess import CompletedProcess

import pytest
from Bio import SeqIO

from plannotate import _package_data
from plannotate._tools import blast, common, diamond, infernal

TEST_DATA = Path(__file__).parent / "test_data"


def test_external_search_failure_includes_tool_diagnostic(monkeypatch):
    monkeypatch.setattr(
        common.subprocess,
        "run",
        lambda *args, **kwargs: CompletedProcess(args, 1, stderr="database missing"),
    )

    with pytest.raises(RuntimeError, match="cmscan.*database missing"):
        common.run_command("cmscan query.fa", "cmscan")


def test_blast_reports_execution_and_hit_count(monkeypatch, tmp_path, caplog):
    observed = {}
    caplog.set_level(logging.DEBUG, logger="plannotate._tools")

    def fake_run(arguments, **kwargs):
        observed["arguments"] = arguments
        return CompletedProcess(arguments, 0, "", "")

    monkeypatch.setattr(common.subprocess, "run", fake_run)
    database_path = tmp_path / "database files" / "custom"

    blast.search(
        "ACGT",
        {
            "db_loc": str(database_path),
            "parameters": "-evalue 1",
        },
    )

    arguments = observed["arguments"]
    assert arguments[arguments.index("-db") + 1] == str(database_path)
    assert "Starting BLAST search" in caplog.text
    assert "Executing blastn command" in caplog.text
    assert "BLAST found 0 candidate hits" in caplog.text


@pytest.mark.integration
@pytest.mark.parametrize(
    ("database_name", "adapter"),
    [
        ("snapgene", blast.search),
        ("swissprot", diamond.search),
        ("Rfam", infernal.search),
    ],
)
def test_external_search_adapters(database_name, adapter):
    databases = _package_data.get_yaml(_package_data.get_yaml_path())
    sequence = str(SeqIO.read(TEST_DATA / "pXampl3.fa", "fasta").seq)

    assert not adapter(sequence, databases[database_name]).empty
