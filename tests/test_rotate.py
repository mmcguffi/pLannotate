"""Unit tests for origin-of-replication rotation (no database required)."""

import random

import pandas as pd
import pytest
from Bio.Seq import Seq

from plannotate import _rotate


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _rot(seq: str, k: int) -> str:
    k %= len(seq)
    return seq[k:] + seq[:k]


# --- _rotate_to_anchor -------------------------------------------------------


def test_rotate_to_anchor_forward_places_ori_at_start():
    ori = "ACGTACGTAC"
    seq = "TTTTT" + ori + "GGGGG"
    n = len(seq)
    rotated, offset, flipped = _rotate._rotate_to_anchor(
        seq, qstart=5, qend=5 + len(ori), sframe=1, n=n
    )
    assert flipped is False
    assert offset == 5
    assert rotated.startswith(ori)
    assert sorted(rotated) == sorted(seq)  # a true rotation


def test_rotate_to_anchor_reverse_places_biological_ori_at_start():
    ori5 = "AAATTTGGGC"  # the ori read 5'->3'
    # embed its reverse complement on the forward strand
    seq = "CCCCC" + _rc(ori5) + "TATATA"
    n = len(seq)
    qstart = 5
    qend = qstart + len(ori5)
    rotated, offset, flipped = _rotate._rotate_to_anchor(
        seq, qstart=qstart, qend=qend, sframe=-1, n=n
    )
    assert flipped is True
    assert offset == (n - qend) % n
    assert rotated.startswith(ori5)


def test_rotate_to_anchor_handles_origin_spanning_forward():
    seq = "ACGTACGTAC"  # n = 10
    n = len(seq)
    # feature wraps: indices 8,9,0,1,2 -> exclusive qend=3 < qstart=8
    rotated, offset, _ = _rotate._rotate_to_anchor(seq, qstart=8, qend=3, sframe=1, n=n)
    assert offset == 8
    assert rotated[:5] == seq[8:] + seq[:3]


def test_rotate_to_anchor_round_trip_restores_original():
    seq = "ACGTTGCAACGGTTACG"
    n = len(seq)
    rotated, offset, flipped = _rotate._rotate_to_anchor(
        seq, qstart=4, qend=9, sframe=1, n=n
    )
    restored = _rot(rotated, -offset)
    assert restored == seq
    assert flipped is False


# --- _canonical_fallback -----------------------------------------------------


@pytest.mark.parametrize("seed", range(8))
def test_fallback_is_rotation_and_strand_invariant(seed):
    rng = random.Random(seed)
    n = rng.randint(20, 120)
    seq = "".join(rng.choice("ACGT") for _ in range(n))

    base = _rotate._canonical_fallback(seq)[0]
    rotated_input = _rotate._canonical_fallback(_rot(seq, rng.randint(1, n - 1)))[0]
    rc_input = _rotate._canonical_fallback(_rc(seq))[0]

    assert base == rotated_input
    assert base == rc_input
    # the canonical frame must itself be a rotation of one of the two strands
    assert base in (seq + seq) or base in (_rc(seq) + _rc(seq))


def test_fallback_is_deterministic():
    seq = "GATTACAGATTACAGGCCTTAAGGCCTT"
    assert _rotate._canonical_fallback(seq) == _rotate._canonical_fallback(seq)


def test_fallback_handles_palindrome():
    seq = "GAATTC" * 4  # reverse-complement palindrome
    first = _rotate._canonical_fallback(seq)
    second = _rotate._canonical_fallback(seq)
    assert first == second
    assert len(first[0]) == len(seq)


def test_fallback_handles_sequence_shorter_than_k():
    seq = "ACGTA"  # shorter than _FALLBACK_K
    rotated, _, _ = _rotate._canonical_fallback(seq)
    assert len(rotated) == len(seq)
    assert rotated in (seq + seq) or rotated in (_rc(seq) + _rc(seq))


# --- ranking + selection -----------------------------------------------------


def test_ranking_table_has_colE1_first_and_origin_types():
    ranking = _rotate._load_ori_ranking()
    bacterial = ranking.loc[ranking["type"] == "bacterial"]
    top = bacterial.sort_values("rank").iloc[0]
    assert top["name"] == "ori"
    assert top["rank"] == 1
    # non-bacterial origins carry their actual category
    assert ranking.loc[ranking["name"] == "f1 ori", "type"].iloc[0] == "phage"
    assert ranking.loc[ranking["name"] == "SV40 ori", "type"].iloc[0] == "viral"


def test_select_origin_prefers_higher_ranked_ori():
    ranking = _rotate._load_ori_ranking()
    candidates = pd.DataFrame(
        [
            {"name": "pSC101 ori", "qstart": 10, "qend": 200, "sframe": 1},
            {"name": "ori", "qstart": 4750, "qend": 5339, "sframe": 1},
        ]
    )
    chosen = _rotate._select_origin(candidates, ranking)
    assert chosen is not None
    assert chosen["name"] == "ori"


def test_select_origin_ignores_non_bacterial_origins():
    ranking = _rotate._load_ori_ranking()
    candidates = pd.DataFrame(
        [
            {"name": "f1 ori", "qstart": 10, "qend": 100, "sframe": 1},
            {"name": "SV40 ori", "qstart": 50, "qend": 150, "sframe": 1},
        ]
    )
    assert _rotate._select_origin(candidates, ranking) is None


def test_select_origin_breaks_ties_by_position():
    ranking = _rotate._load_ori_ranking()
    candidates = pd.DataFrame(
        [
            {"name": "ori", "qstart": 900, "qend": 1400, "sframe": 1},
            {"name": "ori", "qstart": 100, "qend": 600, "sframe": 1},
        ]
    )
    chosen = _rotate._select_origin(candidates, ranking)
    assert chosen is not None
    assert chosen["qstart"] == 100


# --- _detect_origins (monkeypatched search) ----------------------------------


def test_detect_origins_filters_type_and_fragments(monkeypatch):
    from plannotate import _filter, annotate

    # raw enriched hits returned by the snapgene search (mixed types)
    raw = pd.DataFrame(
        [
            {"sseqid": "ori_(1)", "type": "rep_origin", "name": "ori"},
            {"sseqid": "amp", "type": "CDS", "name": "AmpR"},
        ]
    )
    monkeypatch.setattr(
        annotate, "_collect_source_hits", lambda *a, **k: raw, raising=True
    )
    monkeypatch.setattr(
        _rotate._package_data,
        "get_yaml",
        lambda _: {"snapgene": {"method": "blastn"}},
        raising=True,
    )

    # stand-in filter: keep rows, attach the scoring columns _is_fragment needs
    def fake_filter(hits, is_linear=False):
        scored = hits.copy()
        scored["percmatch"] = [99.0]
        scored["pi_permatch"] = [99.0]
        scored["length"] = [589]
        scored["qstart"] = [4750]
        scored["qend"] = [5338]
        scored["sframe"] = [1]
        return scored

    monkeypatch.setattr(_filter, "filter_and_clean_hits", fake_filter, raising=True)

    result = _rotate._detect_origins("ACGT" * 100, yaml_file=None, cores=1)
    assert list(result["name"]) == ["ori"]
    assert bool(result["fragment"].iloc[0]) is False
    assert result["qend"].iloc[0] == 5339  # exclusive end (filter value + 1)


# --- rotate_to_origin orchestration ------------------------------------------


def test_rotate_to_origin_uses_fallback_when_no_ori(monkeypatch):
    monkeypatch.setattr(
        _rotate, "_detect_origins", lambda *a, **k: pd.DataFrame(), raising=True
    )
    seq = "ACGTACGTTTGGGCCCAAATTT"
    result = _rotate.rotate_to_origin(seq)
    assert result.fallback_used is True
    assert result.ori_name is None
    assert result.rotated_seq == _rotate._canonical_fallback(seq)[0]


def test_rotate_to_origin_accepts_seq_and_exposes_seq_output(monkeypatch):
    monkeypatch.setattr(
        _rotate, "_detect_origins", lambda *a, **k: pd.DataFrame(), raising=True
    )
    seq = Seq("ACGTACGTTTGGGCCCAAATTT")
    result = _rotate.rotate_to_origin(seq)
    assert isinstance(result.rotated_seq, str)
    assert isinstance(result.as_seq(), Seq)
    assert str(result.as_seq()) == result.rotated_seq
    # a Seq input yields the same frame as the equivalent str input
    assert result.rotated_seq == _rotate.rotate_to_origin(str(seq)).rotated_seq


def test_rotate_to_origin_rotates_to_detected_ori(monkeypatch):
    ori = "ACGTACGTAC"
    seq = "TTTTT" + ori + "GGGGG"
    candidates = pd.DataFrame(
        [{"name": "ori", "qstart": 5, "qend": 5 + len(ori), "sframe": 1}]
    )
    monkeypatch.setattr(
        _rotate, "_detect_origins", lambda *a, **k: candidates, raising=True
    )
    result = _rotate.rotate_to_origin(seq)
    assert result.fallback_used is False
    assert result.ori_name == "ori"
    assert result.rank == 1
    assert result.rotated_seq.startswith(ori)


# --- public API surface ------------------------------------------------------


def test_public_exports_are_available():
    import plannotate

    assert {"rotate_to_origin", "RotationResult"} <= set(plannotate.__all__)
    # the exported names are the same objects as the implementation module's
    assert plannotate.rotate_to_origin is _rotate.rotate_to_origin
    assert plannotate.RotationResult is _rotate.RotationResult


def test_public_rotate_to_origin_accepts_str_and_seq(monkeypatch):
    import plannotate

    monkeypatch.setattr(
        _rotate, "_detect_origins", lambda *a, **k: pd.DataFrame(), raising=True
    )
    raw = "ACGTACGTTTGGGCCCAAATTT"

    from_str = plannotate.rotate_to_origin(raw)
    from_seq = plannotate.rotate_to_origin(Seq(raw))

    assert isinstance(from_str, plannotate.RotationResult)
    assert from_str.rotated_seq == from_seq.rotated_seq
    assert isinstance(from_seq.as_seq(), Seq)
    assert str(from_seq.as_seq()) == from_seq.rotated_seq


def test_rotation_result_as_seq_round_trips():
    result = _rotate.RotationResult(
        rotated_seq="ACGTACGT",
        offset=0,
        flipped=False,
        ori_name=None,
        rank=None,
        fallback_used=True,
    )
    seq = result.as_seq()
    assert isinstance(seq, Seq)
    assert str(seq) == "ACGTACGT"
