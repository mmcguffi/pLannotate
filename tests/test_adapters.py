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


def test_read_table_keeps_text_columns_out_of_numeric_inference(tmp_path):
    # a gapless, fully identical alignment is reported as a bare match run, which
    # numeric inference would silently turn into an integer
    output = tmp_path / "results.tsv"
    output.write_text("query\t100\t300\nquery\t100\t12AG287\n")

    table = common.read_table(str(output), "qseqid pident btop", ("btop",))

    assert list(table["btop"]) == ["300", "12AG287"]
    assert list(table["pident"]) == [100, 100]


@pytest.mark.parametrize("adapter", [blast, diamond])
def test_search_adapters_request_the_alignment_traceback(adapter):
    assert "btop" in adapter.COLUMNS.split()
    assert "btop" in adapter.TEXT_COLUMNS


def test_diamond_subject_coordinates_are_nucleotide_equivalents(monkeypatch, tmp_path):
    output_row = "q0\t1\t90\tP12345\t100\t100\t" + "A" * 90
    output_row += "\t30\t2\t31\t90\t1e-20\t30\n"

    def fake_run(_arguments, **kwargs):
        # The adapter's temporary output path is the value following -o.
        arguments = _arguments
        Path(arguments[arguments.index("-o") + 1]).write_text(output_row)
        return CompletedProcess(arguments, 0, "", "")

    monkeypatch.setattr(common.subprocess, "run", fake_run)
    result = diamond.search(
        "A" * 90,
        {"db_loc": str(tmp_path / "protein"), "parameters": ""},
    )

    row = result.iloc[0]
    assert (row["sstart"], row["send"], row["slen"]) == (4, 93, 300)


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
