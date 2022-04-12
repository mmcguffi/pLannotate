from plannotate import resources
from plannotate import annotate
import pytest

RRNB = "tgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctact"


def test_annotate():
    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


def test_BLAST():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    hits = annotate.BLAST(RRNB, snapgene_db)
    assert not hits.empty


def test_nothing():
    assert True


def test_get_fasta():
    print(resources.get_example_fastas())
