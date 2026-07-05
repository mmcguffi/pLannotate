"""Tests for DNA sequence and input-file validation."""

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from plannotate.validation import (
    InvalidSequenceError,
    validate_file,
    validate_records,
    validate_sequence,
)

TEST_DATA = Path(__file__).parent / "test_data"


def test_validate_records_accepts_multiple_entries(tmp_path):
    fasta = tmp_path / "multi.fa"
    fasta.write_text(">alpha\nACGTACGT\n>beta\nTTTTGGGG\n")

    records = validate_records(fasta, max_length=None)

    assert [record.id for record in records] == ["alpha", "beta"]


def test_validate_records_rejects_an_invalid_entry(tmp_path):
    fasta = tmp_path / "multi.fa"
    fasta.write_text(">alpha\nACGT\n>beta\nACGZ\n")

    with pytest.raises(InvalidSequenceError):
        validate_records(fasta, max_length=None)


def _record(identifier: str) -> SeqRecord:
    record = SeqRecord(Seq("ACGT"), id=identifier, name=identifier)
    record.annotations["molecule_type"] = "DNA"
    return record


def test_genbank_validation_preserves_features_and_accepts_uppercase_extension(
    tmp_path,
):
    record = _record("test")
    record.features.append(SeqFeature(FeatureLocation(0, 2), type="promoter"))
    path = tmp_path / "input.GBK"
    SeqIO.write(record, path, "genbank")

    validated = validate_file(str(path), path.suffix)

    assert len(validated.features) == 1
    assert validated.features[0].type == "promoter"


def test_genbank_validation_rejects_multiple_records(tmp_path):
    path = tmp_path / "multiple.gbk"
    SeqIO.write([_record("one"), _record("two")], path, "genbank")

    with pytest.raises(InvalidSequenceError, match="multiple entries"):
        validate_file(str(path), path.suffix)


def test_sequence_validation_rejects_empty_input():
    with pytest.raises(InvalidSequenceError, match="empty"):
        validate_sequence("")


def test_file_validation_infers_extension(tmp_path):
    path = tmp_path / "input.fasta"
    SeqIO.write(_record("test"), path, "fasta")

    assert validate_file(path).id == "test"


def test_sequence_validation_rejects_invalid_nucleotide():
    with pytest.raises(InvalidSequenceError, match="invalid characters"):
        validate_sequence("ACXT")


def test_name_and_extension_are_normalized():
    from plannotate.validation import get_name_ext

    assert get_name_ext("/a/long/path/test.FASTA") == ("test", ".fasta")


@pytest.mark.parametrize("extension", ["fasta", "fa", "fas", "fna"])
def test_supported_fasta_extensions(extension):
    record = validate_file(TEST_DATA / f"pAdDeltaF6.{extension}")

    assert len(record) == 15_420


def test_validation_rejects_unknown_extension():
    with pytest.raises(ValueError, match="must be a FASTA or GenBank file"):
        validate_file(TEST_DATA / "pAdDeltaF6.txt")
