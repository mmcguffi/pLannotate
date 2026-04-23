import warnings
from io import StringIO

import pytest
from Bio import BiopythonParserWarning, SeqIO

from plannotate import resources

pytestmark = pytest.mark.integration

with open("./tests/test_data/RRNB_fragment.txt") as f:
    RRNB = f.read()

SAM_RIBOSWITCH = (
    "TTTCTATCCAGAGAGGTGGAGGGACTGGCCCTATGAAACCTCGGCAACATTATTGTGCCAATTCCAGCAA"
    "GCGCTAGCTTGAAAGATAGGAA"
)
SAM_RIBOSWITCH_5P_PAD_1 = "A" + SAM_RIBOSWITCH
SAM_RIBOSWITCH_3P_PAD_1 = SAM_RIBOSWITCH + "A"
ISSUE_60_ORIGINAL_SEQUENCE = (
    "cccccacctggcgacaggtgcctctgcggccaaaagccacgtgtataagatacacctgcaaaggcggcaca"
    "accccagtgccacgttgtgagttggatagttgtggaaagagtcaaatggctctcctcaagcgtattcaaca"
    "aggggctgaaggatgcccagaaggtacc"
)


def get_hits_by_db(sequence, db, linear=True):
    from plannotate import annotate

    hits = annotate.annotate(sequence, linear=linear)
    return hits.loc[hits["db"] == db].reset_index(drop=True)


def assert_gbk_round_trips_without_parser_warning(
    gbk_text, expected_start, expected_end
):
    assert "0.." not in gbk_text

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gbk = SeqIO.read(StringIO(gbk_text), "genbank")

    assert not any(issubclass(w.category, BiopythonParserWarning) for w in caught)
    assert gbk.features[0].location is not None
    assert int(gbk.features[0].location.start) == expected_start
    assert int(gbk.features[0].location.end) == expected_end


def test_BLAST():
    from plannotate import annotate

    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    hits = annotate.BLAST(RRNB, snapgene_db)
    assert not hits.empty


def test_annotate():
    from plannotate import annotate

    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
    hits = annotate.annotate(RRNB, linear=True)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


@pytest.mark.parametrize("linear", [True, False])
def test_issue_60_rfam_left_edge_hit_round_trips_after_fix(linear):
    rfam_hits = get_hits_by_db(SAM_RIBOSWITCH, "Rfam", linear=linear)

    assert len(rfam_hits) == 1
    assert rfam_hits.loc[0, "qstart"] == 0
    assert rfam_hits.loc[0, "qend"] == 92
    assert rfam_hits.loc[0, "qseq"] == SAM_RIBOSWITCH

    gbk_text = resources.get_gbk(rfam_hits, SAM_RIBOSWITCH, is_linear=linear)
    assert_gbk_round_trips_without_parser_warning(gbk_text, 0, 92)


@pytest.mark.parametrize(
    ("sequence", "expected_start", "expected_end"),
    [
        (SAM_RIBOSWITCH_3P_PAD_1, 0, 92),
        (SAM_RIBOSWITCH_5P_PAD_1, 1, 93),
    ],
)
def test_issue_60_rfam_padding_behaves_by_padding_direction(
    sequence, expected_start, expected_end
):
    rfam_hits = get_hits_by_db(sequence, "Rfam")

    assert len(rfam_hits) == 1
    assert rfam_hits.loc[0, "qstart"] == expected_start
    assert rfam_hits.loc[0, "qend"] == expected_end

    gbk_text = resources.get_gbk(rfam_hits, sequence, is_linear=True)
    assert_gbk_round_trips_without_parser_warning(
        gbk_text, expected_start, expected_end
    )


def test_issue_60_original_sequence_snapgene_hit_is_zero_based_in_df_but_valid_in_gbk():
    snapgene_hits = get_hits_by_db(ISSUE_60_ORIGINAL_SEQUENCE, "snapgene")

    assert len(snapgene_hits) == 1
    assert snapgene_hits.loc[0, "qstart"] == 0
    assert snapgene_hits.loc[0, "qend"] == len(ISSUE_60_ORIGINAL_SEQUENCE)

    gbk_text = resources.get_gbk(
        snapgene_hits, ISSUE_60_ORIGINAL_SEQUENCE, is_linear=True
    )
    assert "1..170" in gbk_text

    assert_gbk_round_trips_without_parser_warning(
        gbk_text, 0, len(ISSUE_60_ORIGINAL_SEQUENCE)
    )
