import pandas as pd


def _hit(**overrides):
    row = {
        "sseqid": "keep",
        "qstart": 10,
        "qend": 30,
        "qlen": 100,
        "wstart": 12,
        "wend": 28,
        "kind": 1,
        "evalue": 0.0,
        "pi_permatch": 100.0,
    }
    row.update(overrides)
    return row


def test_clean_filters_known_bad_and_low_quality_hits():
    from plannotate.annotate import clean

    hits = pd.DataFrame(
        [
            _hit(sseqid="keep"),
            _hit(sseqid="P03851"),
            _hit(sseqid="high_evalue", evalue=1.0),
            _hit(sseqid="poor_match", pi_permatch=3.0),
        ]
    )

    cleaned = clean(hits)

    assert cleaned["sseqid"].to_list() == ["keep"]


def test_clean_keeps_first_same_kind_overlap():
    from plannotate.annotate import clean

    hits = pd.DataFrame(
        [
            _hit(sseqid="first", qstart=10, qend=40, wstart=15, wend=35, kind=1),
            _hit(sseqid="second", qstart=25, qend=50, wstart=30, wend=45, kind=1),
        ]
    )

    cleaned = clean(hits)

    assert cleaned["sseqid"].to_list() == ["first"]


def test_clean_allows_different_kind_overlap():
    from plannotate.annotate import clean

    hits = pd.DataFrame(
        [
            _hit(sseqid="cds", qstart=10, qend=40, wstart=15, wend=35, kind="CDS"),
            _hit(
                sseqid="promoter",
                qstart=25,
                qend=50,
                wstart=30,
                wend=45,
                kind="promoter",
            ),
        ]
    )

    cleaned = clean(hits)

    assert cleaned["sseqid"].to_list() == ["cds", "promoter"]


def test_clean_preserves_origin_crossing_coordinates_and_duplicates():
    from plannotate.annotate import clean

    hits = pd.DataFrame(
        [
            _hit(
                sseqid="origin_crossing",
                qstart=95,
                qend=105,
                wstart=97,
                wend=102,
            ),
        ]
    )

    cleaned = clean(hits)

    assert cleaned.loc[0, "qstart_dup"] == 95
    assert cleaned.loc[0, "qend_dup"] == 105
    assert cleaned.loc[0, "qstart"] == 95
    assert cleaned.loc[0, "qend"] == 5
    assert cleaned.loc[0, "wstart"] == 97
    assert cleaned.loc[0, "wend"] == 2
