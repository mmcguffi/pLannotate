"""End-to-end annotation integration tests."""

import warnings
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import BiopythonParserWarning, SeqIO

from plannotate import annotate
from plannotate.models import Construct, df_to_features

pytestmark = pytest.mark.integration
TEST_DATA = Path(__file__).parent / "test_data"

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


@pytest.mark.parametrize("linear", [True, False])
def test_rrnb_terminator_annotation(linear):
    sequence = (TEST_DATA / "RRNB_fragment.txt").read_text()

    hits = annotate.annotate(sequence, linear=linear)

    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


def _serialized_features(construct):
    return pd.DataFrame(
        {
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "strand": feature.location.strand,
            "type": feature.type,
            "is_frag": feature.qualifiers["fragment"],
            "database": feature.qualifiers["database"],
            "name": feature.qualifiers["label"],
        }
        for feature in construct.to_seqrecord().features
    ).sort_values(["start", "end"], ignore_index=True)


@pytest.mark.parametrize(
    ("detailed", "ground_truth"),
    [
        (False, "RNAs_ground-truth.csv"),
        (True, "RNAs_ground-truth-detailed.csv"),
    ],
)
def test_rna_annotation_matches_ground_truth(detailed, ground_truth):
    sequence = SeqIO.read(TEST_DATA / "RNAs.fasta", "fasta").seq
    actual = _serialized_features(Construct(sequence, detailed=detailed))
    expected = pd.read_csv(TEST_DATA / ground_truth)

    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)


@pytest.mark.parametrize(
    ("linear", "ground_truth"),
    [
        (False, "origin_crossing_araC_ground_truth_circular.csv"),
        (True, "origin_crossing_araC_ground_truth_linear.csv"),
    ],
)
def test_origin_crossing_annotation(linear, ground_truth):
    sequence = SeqIO.read(TEST_DATA / "origin_crossing_araC.fa", "fasta").seq
    results = Construct(sequence, linear=linear).annotations_df
    key_columns = ["sseqid", "qstart", "qend", "sframe", "db", "fragment"]
    actual = results[key_columns].sort_values(["qstart", "qend"], ignore_index=True)
    expected = pd.read_csv(TEST_DATA / ground_truth).sort_values(
        ["qstart", "qend"], ignore_index=True
    )

    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    arac = results.loc[results["sseqid"] == "araC"]
    assert len(arac) == (2 if linear else 1)
    assert bool(arac["fragment"].all()) is linear
