import warnings
from io import StringIO

import pytest
from Bio import BiopythonParserWarning, SeqIO

from plannotate import annotate
from plannotate.models import Construct, df_to_features


pytestmark = pytest.mark.integration

SAM_RIBOSWITCH = (
    "TTTCTATCCAGAGAGGTGGAGGGACTGGCCCTATGAAACCTCGGCAACATTATTGTGCCAATTCCAGCAA"
    "GCGCTAGCTTGAAAGATAGGAA"
)
SAM_RIBOSWITCH_5P_PAD_1 = "A" + SAM_RIBOSWITCH
SAM_RIBOSWITCH_3P_PAD_1 = SAM_RIBOSWITCH + "A"


def _rfam_hits(sequence, linear):
    hits = annotate.annotate(sequence, linear=linear)
    return hits.loc[hits["db"] == "Rfam"].reset_index(drop=True)


def _to_genbank(sequence, hits, linear):
    construct = Construct(seq=sequence, linear=linear, _skip_annotation=True)
    construct.features = df_to_features(hits)
    return construct.to_genbank()


def _assert_round_trip(gbk_text, expected_start, expected_end):
    assert "0.." not in gbk_text
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gbk = SeqIO.read(StringIO(gbk_text), "genbank")
    assert not any(
        issubclass(warning.category, BiopythonParserWarning) for warning in caught
    )
    assert int(gbk.features[0].location.start) == expected_start
    assert int(gbk.features[0].location.end) == expected_end


@pytest.mark.parametrize("linear", [True, False])
def test_issue_60_rfam_left_edge_round_trip(linear):
    hits = _rfam_hits(SAM_RIBOSWITCH, linear)

    assert len(hits) == 1
    assert hits.loc[0, "qstart"] == 0
    assert hits.loc[0, "qend"] == 92
    assert hits.loc[0, "qseq"] == SAM_RIBOSWITCH
    _assert_round_trip(_to_genbank(SAM_RIBOSWITCH, hits, linear), 0, 92)


@pytest.mark.parametrize(
    ("sequence", "expected_start", "expected_end"),
    [
        (SAM_RIBOSWITCH_3P_PAD_1, 0, 92),
        (SAM_RIBOSWITCH_5P_PAD_1, 1, 93),
    ],
)
def test_issue_60_rfam_padding(sequence, expected_start, expected_end):
    hits = _rfam_hits(sequence, linear=True)

    assert len(hits) == 1
    assert hits.loc[0, "qstart"] == expected_start
    assert hits.loc[0, "qend"] == expected_end
    _assert_round_trip(
        _to_genbank(sequence, hits, linear=True), expected_start, expected_end
    )
