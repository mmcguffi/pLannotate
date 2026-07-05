"""Integration tests for --rotate (require external tools and databases)."""

from pathlib import Path

import pytest
from Bio import SeqIO
from typer.testing import CliRunner

from plannotate import Construct
from plannotate.main import app

pytestmark = pytest.mark.integration
TEST_DATA = Path(__file__).parent / "test_data"


def _read_fasta(name):
    return str(SeqIO.read(TEST_DATA / name, "fasta").seq)


def test_construct_rotates_pXampl3_ori_to_base_one():
    construct = Construct(seq=_read_fasta("pXampl3.fa"), rotate=True)

    assert construct.rotation is not None
    assert construct.rotation.fallback_used is False
    assert construct.rotation.ori_name == "ori"

    oris = [
        f
        for f in construct.features
        if f.feature_type == "rep_origin" and f.feature_name == "ori"
    ]
    assert oris, "expected a ColE1 'ori' annotation"
    ori = oris[0]
    assert ori.qstart == 0
    assert ori.sframe == 1


def test_rotation_preserves_sequence_as_a_true_rotation():
    original = _read_fasta("pXampl3.fa")
    construct = Construct(seq=original, rotate=True)
    rotated = str(construct.seq)
    assert len(rotated) == len(original)
    # forward-strand rotation: rotated must be a cyclic permutation of the input
    assert rotated in (original + original)


def test_cli_rotate_emits_ori_at_start(tmp_path):
    result = CliRunner().invoke(
        app,
        [
            "batch",
            "--input",
            str(TEST_DATA / "pXampl3.fa"),
            "--output",
            str(tmp_path),
            "--suffix",
            "",
            "--rotate",
        ],
    )
    assert result.exit_code == 0, result.exception

    record = SeqIO.read(tmp_path / "pXampl3.gbk", "genbank")
    ori_features = [
        f
        for f in record.features
        if f.type == "rep_origin" and f.qualifiers.get("label", [""])[0] == "ori"
    ]
    assert ori_features
    ori = ori_features[0]
    assert int(ori.location.start) == 0
    assert ori.location.strand == 1


def test_rotate_fallback_is_deterministic_without_ori(tmp_path):
    seq = _read_fasta("random_dna.fa")
    first = Construct(seq=seq, rotate=True)
    second = Construct(seq=seq, rotate=True)

    assert first.rotation is not None
    assert first.rotation.fallback_used is True
    assert str(first.seq) == str(second.seq)
