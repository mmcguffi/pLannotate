"""Tests for the circular-coordinate contract in hit filtering.

These exercise ``filter_and_clean_hits`` directly so the doubled-query handling
is verified without needing the external search tools or databases. They rely on
``qlen`` already being the true sequence length (set when hits are collected),
which is what lets the windowed/doubled second copy wrap back correctly.
"""

import pandas as pd

from plannotate._filter import (
    PERFECT_MATCH_BONUS,
    calculate_hit_scores,
    filter_and_clean_hits,
)

SEQUENCE_LENGTH = 100


def _hit(sseqid, qstart, qend, *, kind=1, priority=3):
    """Build one raw hit row using 1-based coordinates, as the tools report them."""
    length = abs(qend - qstart) + 1
    return {
        "sseqid": sseqid,
        "qstart": qstart,
        "qend": qend,
        "qlen": SEQUENCE_LENGTH,
        "length": length,
        "slen": length,
        "pident": 100,
        "priority": priority,
        "evalue": 1e-10,
        "kind": kind,
    }


def test_second_copy_duplicate_collapses_onto_first_copy():
    # the same feature found in both copies of the doubled query must reduce to one
    hits = pd.DataFrame(
        [
            _hit("rna", 11, 30),  # first copy
            _hit("rna", 111, 130),  # second copy (+ SEQUENCE_LENGTH)
        ]
    )

    result = filter_and_clean_hits(hits, is_linear=False)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["qlen"] == SEQUENCE_LENGTH
    # surviving coordinates are the wrapped, 0-based first-copy interval
    assert (row["qstart"], row["qend"]) == (10, 29)


def test_origin_spanning_hit_wraps_to_a_split_interval():
    # a hit crossing the origin keeps qstart > qend so downstream builds a join()
    hits = pd.DataFrame([_hit("ori", 95, 110)])

    result = filter_and_clean_hits(hits, is_linear=False)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["qlen"] == SEQUENCE_LENGTH
    assert row["qstart"] == 94  # 0-based, stays on the real sequence
    assert row["qend"] == 9  # 110 wrapped back by SEQUENCE_LENGTH
    assert row["qstart"] > row["qend"]


def _scoring_row(**overrides):
    row = {
        "qstart": 1,
        "qend": 100,
        "sstart": 1,
        "send": 100,
        "length": 100,
        "slen": 100,
        "pident": 100.0,
        "priority": 1,
        "sequence_identity": True,
    }
    row.update(overrides)
    return row


def test_perfect_match_bonus_skips_a_confidence_based_identity():
    # A covariance-model hit reports the model's confidence in pident, not a count of
    # matching bases, and Infernal prints it to two decimals. Granting the exact-match
    # bonus at 1.00 would put a tenfold score step on a rounding boundary: the 0.99
    # hit below is no worse a match, and the two must stay within a factor of ten.
    hits = pd.DataFrame(
        [
            _scoring_row(sequence_identity=False, pident=100.0),
            _scoring_row(sequence_identity=False, pident=99.0),
            _scoring_row(sequence_identity=True, pident=100.0),
        ]
    )

    scored = calculate_hit_scores(hits)

    confident, nearly, sequence = scored["score"].tolist()
    assert confident / nearly < 2
    # a real full-length sequence match still earns it
    assert sequence == confident * PERFECT_MATCH_BONUS


def test_perfect_match_bonus_still_applies_without_the_column():
    # a frame built before the column existed keeps the historical behaviour
    hits = pd.DataFrame([_scoring_row()]).drop(columns="sequence_identity")

    assert calculate_hit_scores(hits).loc[0, "score"] == 100 * PERFECT_MATCH_BONUS
