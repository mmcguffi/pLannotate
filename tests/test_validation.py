import os.path as op

import pytest

from plannotate import resources


def test_valid_sequence_correct():
    resources.validate_sequence("ACTG")


def test_valid_sequence_incorrect_base():
    with pytest.raises(ValueError):
        resources.validate_sequence("ACTX")


def test_get_name_ext():
    name, ext = resources.get_name_ext("./a/long/path/test.fasta")
    assert name == "test"
    assert ext == ".fasta"


def test_validate_file_fa():
    df_path = op.join(__package__, "test_data", "pXampl3.fa")
    resources.validate_file(df_path, ".fasta")


def test_validate_file_gbk():
    df_path = op.join(__package__, "test_data", "pXampl3_detailed.gbk")
    resources.validate_file(df_path, ".gbk")


@pytest.mark.parametrize("ext", ["fasta", "fa", "fas", "fna"])
def test_validate_file_all_fasta_extensions(ext):
    input_file = f"tests/test_data/pAdDeltaF6.{ext}"
    name, ext = resources.get_name_ext(input_file)
    sequence = resources.validate_file(input_file, ext)
    assert len(sequence) == 15420


def test_validate_file_bad_extension():
    input_file = "tests/test_data/pAdDeltaF6.txt"
    name, ext = resources.get_name_ext(input_file)
    with pytest.raises(ValueError, match="must be a FASTA or GenBank file"):
        _ = resources.validate_file(input_file, ext)
