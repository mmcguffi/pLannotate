"""Integration tests for command-line outputs."""

from pathlib import Path

import pytest
from Bio import SeqIO
from typer.testing import CliRunner

from plannotate.main import app

pytestmark = pytest.mark.integration
TEST_DATA = Path(__file__).parent / "test_data"


def _run_batch(input_name, output, *options):
    result = CliRunner().invoke(
        app,
        [
            "batch",
            "--input",
            str(TEST_DATA / input_name),
            "--output",
            str(output),
            "--suffix",
            "",
            *options,
        ],
    )
    assert result.exit_code == 0, result.exception
    return output / Path(input_name).with_suffix("")


def test_batch_writes_annotated_genbank(tmp_path):
    output_stem = _run_batch("pXampl3.fa", tmp_path)

    record = SeqIO.read(output_stem.with_suffix(".gbk"), "genbank")
    assert len(record.features) > 15


def test_batch_handles_sequence_without_hits(tmp_path):
    output_stem = _run_batch("random_dna.fa", tmp_path, "--html")

    record = SeqIO.read(output_stem.with_suffix(".gbk"), "genbank")
    assert record.features == []
    assert output_stem.with_suffix(".html").is_file()


def test_batch_handles_missing_feature_metadata(tmp_path):
    output_stem = _run_batch("nan_feature.fa", tmp_path)

    record = SeqIO.read(output_stem.with_suffix(".gbk"), "genbank")
    assert len(record.features) == 2
