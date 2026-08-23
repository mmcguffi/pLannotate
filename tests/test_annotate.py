"""Tests for candidate collection and annotation helpers."""

import os

import pandas as pd
import pytest

from plannotate import annotate


def test_cached_hits_are_independent_and_config_aware(monkeypatch, tmp_path):
    yaml_file = tmp_path / "databases.yml"
    yaml_file.write_text("first: {}")
    calls = []
    monkeypatch.setattr(annotate._package_data, "get_yaml", lambda _: {"first": {}})

    def fake_sources(*args):
        calls.append(args)
        return [pd.DataFrame({"db": ["first"]})]

    monkeypatch.setattr(annotate._concurrency, "run_sources", fake_sources)
    annotate._collect_hits_cached.cache_clear()

    first = annotate._collect_hits("ACGT", True, yaml_file, cores=1)
    first.loc[0, "db"] = "mutated"
    second = annotate._collect_hits("ACGT", True, yaml_file, cores=1)

    assert second.loc[0, "db"] == "first"
    assert len(calls) == 1

    previous_mtime = yaml_file.stat().st_mtime_ns
    yaml_file.write_text("first: {priority: 2}")
    if yaml_file.stat().st_mtime_ns == previous_mtime:
        os.utime(yaml_file, ns=(previous_mtime + 1, previous_mtime + 1))
    annotate._collect_hits("ACGT", True, yaml_file, cores=1)

    assert len(calls) == 2


def test_fast_mode_restricts_sources_to_cheap_subset(monkeypatch, tmp_path):
    yaml_file = tmp_path / "databases.yml"
    yaml_file.write_text("unused: {}")
    all_sources = {name: {} for name in ("snapgene", "fpbase", "swissprot", "Rfam")}
    monkeypatch.setattr(annotate._package_data, "get_yaml", lambda _: all_sources)

    scheduled = []

    def fake_sources(*args):
        sources = args[3]
        scheduled.append(set(sources))
        return [pd.DataFrame({"db": list(sources)})]

    monkeypatch.setattr(annotate._concurrency, "run_sources", fake_sources)
    annotate._collect_hits_cached.cache_clear()

    annotate._collect_hits("ACGT", True, yaml_file, cores=1, fast=True)
    assert scheduled[-1] == set(annotate.FAST_SOURCES)

    annotate._collect_hits("ACGT", True, yaml_file, cores=1, fast=False)
    assert scheduled[-1] == set(all_sources)


def _seam_fragment(**overrides):
    row: dict[str, object] = {column: 0 for column in annotate.ADAPTER_COLUMNS}
    row.update(
        sseqid="featA",
        sframe=1,
        evalue=1e-9,
        qseq="A",
        pident=100.0,
        slen=100,
        qlen=1000,
    )
    row.update(overrides)
    return row


def test_stitch_merges_adjacent_seam_fragments():
    # forward feature split by the seam: 3' fragment ends at qlen, 5' fragment
    # starts at base 1, and their subject coords are contiguous (60 -> 61).
    hits = pd.DataFrame(
        [
            _seam_fragment(
                qstart=941,
                qend=1000,
                sstart=1,
                send=60,
                length=60,
                qseq="RIGHT",
                btop="60",
            ),
            _seam_fragment(
                qstart=1,
                qend=40,
                sstart=61,
                send=100,
                length=40,
                qseq="LEFT",
                btop="40",
            ),
        ]
    )

    merged = annotate._stitch_seam_hits(hits)

    assert len(merged) == 1
    row = merged.iloc[0]
    assert row["qstart"] == 941
    assert row["qend"] == 1040  # qlen + left_end, wrapped downstream by _filter
    assert row["length"] == 100
    assert row["qseq"] == "RIGHTLEFT"
    assert row["sstart"] == 1
    assert row["send"] == 100
    # a fused hit reports no traceback: joining "60" and "40" would read as a single
    # 6040-base match run, not the 100 bases actually aligned
    assert row["btop"] == ""


def test_stitch_keeps_a_reverse_strand_subject_span_descending():
    # blastn reports a minus-strand hit descending (sstart > send), so the merged span
    # runs 100 -> 1. Pairing the two rows' starts against their ends instead would
    # report 60 -> 61 -- two bases at the join, for a 100-base match.
    hits = pd.DataFrame(
        [
            _seam_fragment(
                qstart=941, qend=1000, sstart=100, send=61, length=60, sframe=-1
            ),
            _seam_fragment(qstart=1, qend=40, sstart=60, send=1, length=40, sframe=-1),
        ]
    )

    merged = annotate._stitch_seam_hits(hits)

    assert len(merged) == 1
    row = merged.iloc[0]
    assert row["sstart"] == 100
    assert row["send"] == 1


def test_stitch_counts_an_overlapping_subject_region_once():
    # an aligner can extend both fragments a little past the seam, so the pair covers
    # subject 1-60 and 59-100: two units in common. Summing the fragment lengths would
    # report a 102-unit match against a 100-unit subject, pushing every downstream
    # match fraction above 100%.
    hits = pd.DataFrame(
        [
            _seam_fragment(qstart=941, qend=1000, sstart=1, send=60, length=60),
            _seam_fragment(qstart=1, qend=40, sstart=59, send=100, length=42),
        ]
    )

    merged = annotate._stitch_seam_hits(hits)

    assert len(merged) == 1
    row = merged.iloc[0]
    assert row["length"] == 100
    assert row["length"] <= row["slen"]
    assert (row["sstart"], row["send"]) == (1, 100)


def test_stitch_leaves_non_contiguous_fragments_alone():
    # both touch the termini but their subject ranges overlap heavily, so they are
    # two distinct features, not one seam-spanning feature.
    hits = pd.DataFrame(
        [
            _seam_fragment(qstart=941, qend=1000, sstart=1, send=60, length=60),
            _seam_fragment(qstart=1, qend=40, sstart=1, send=40, length=40),
        ]
    )

    merged = annotate._stitch_seam_hits(hits)

    assert len(merged) == 2


def test_stitch_preserves_three_residue_diamond_seam_tolerance():
    hits = pd.DataFrame(
        [
            _seam_fragment(
                qstart=941,
                qend=1000,
                sstart=1,
                send=180,
                length=60,
            ),
            _seam_fragment(
                qstart=1,
                qend=40,
                sstart=190,
                send=309,
                length=40,
            ),
        ]
    )

    # Subject coordinates are nucleotide-equivalents: 180 -> 190 leaves nine
    # positions, the same three-residue tolerance used before DIAMOND normalization.
    assert len(annotate._stitch_seam_hits(hits)) == 1
    too_far = hits.copy()
    too_far.at[1, "sstart"] = 193
    too_far.at[1, "send"] = 312
    assert len(annotate._stitch_seam_hits(too_far)) == 2


def test_build_search_queries_does_not_double_linear_or_fast():
    seqs = {"q0": "ACGT" * 25}  # length 100

    for kwargs in (
        {"is_linear": True, "fast": False},
        {"is_linear": False, "fast": True},
    ):
        queries, true_lens = annotate._build_search_queries(seqs, **kwargs)
        # linear and fast searches use a single copy; seam features are stitched later
        assert queries["q0"] == seqs["q0"]
        assert true_lens == {"q0": 100}


def test_missing_swissprot_description_has_no_priority_penalty():
    assert annotate._existence_level_priority(None) == 0


def test_load_feature_details_none_location_synthesizes_missing_columns():
    # a hand-configured blast/diamond source with details.location None: its hits
    # carry no name/type/blurb columns, so the details must be synthesized (name
    # defaults to the id) instead of raising KeyError.
    hits = pd.DataFrame({"sseqid": ["featA", "featA", "featB"]})
    config = {"details": {"location": None}, "priority": 1}

    details = annotate._load_feature_details(hits, "customdb", config)

    assert list(details.columns) == ["sseqid", "name", "type", "blurb"]
    # collapsed to one row per id so the downstream merge cannot fan out
    assert list(details["sseqid"]) == ["featA", "featB"]
    assert list(details["name"]) == ["featA", "featB"]  # defaulted from the id
    assert set(details["type"]) == {"misc_feature"}
    assert set(details["blurb"]) == {""}


def test_enrich_hits_recovers_underscore_id_rewritten_as_pdb(monkeypatch):
    hits = pd.DataFrame({"sseqid": ["pdb|12Pk|tag", "pdb|1ABC|A"]})

    def fake_details(candidate_hits, _source, _config):
        # The metadata bundle contains the original underscore id, while a normal
        # PDB record is keyed by accession alone.
        assert {"12Pk_tag", "1ABC"} <= set(candidate_hits["sseqid"])
        return pd.DataFrame(
            [
                {
                    "sseqid": "12Pk_tag",
                    "name": "12Pk tag",
                    "type": "CDS",
                    "blurb": "",
                },
                {
                    "sseqid": "1ABC",
                    "name": "structure",
                    "type": "misc_feature",
                    "blurb": "",
                },
            ]
        )

    monkeypatch.setattr(annotate, "_load_feature_details", fake_details)
    enriched = annotate._enrich_hits(
        hits,
        "snapgene",
        {"priority": 1, "details": {"location": "unused"}},
    )

    assert enriched["sseqid"].tolist() == ["12Pk_tag", "1ABC"]
    assert enriched["name"].tolist() == ["12Pk tag", "structure"]
    assert enriched["type"].tolist() == ["CDS", "misc_feature"]


def test_enrich_hits_preserves_inline_rfam_metadata():
    hits = pd.DataFrame(
        {
            "sseqid": ["RF00106"],
            "name": ["RNAI"],
            "type": ["ncRNA"],
            "blurb": ["ColE1 antisense regulator"],
        }
    )

    enriched = annotate._enrich_hits(
        hits,
        "Rfam",
        {"priority": 1, "details": {"location": None}},
    )

    row = enriched.iloc[0]
    assert row["name"] == "RNAI"
    assert row["type"] == "ncRNA"
    assert row["blurb"] == "ColE1 antisense regulator"


def test_circular_search_query_doubles_sequence():
    query = "ACGT" * 25  # length 100

    # circular queries are fully doubled so origin-spanning features are recovered
    assert annotate._circular_search_query(query) == query + query


def test_build_search_queries_doubles_circular_and_records_true_length():
    seqs = {"q0": "ACGT" * 25}  # length 100

    queries, true_lens = annotate._build_search_queries(
        seqs, is_linear=False, fast=False
    )

    # circular -> fully doubled (lossless)
    assert queries["q0"] == seqs["q0"] + seqs["q0"]
    # true length is tracked separately from the doubled search length
    assert true_lens == {"q0": 100}


def test_collect_source_hits_maps_true_length_per_query(monkeypatch):
    def fake_run(queries, config, threads):
        rows = []
        for query_id in queries:
            row = {column: 1 for column in annotate.ADAPTER_COLUMNS}
            row["qseqid"] = query_id
            rows.append(row)
        return pd.DataFrame(rows)

    monkeypatch.setattr(annotate, "run_tool", fake_run)
    monkeypatch.setattr(annotate, "_enrich_hits", lambda hits, name, config: hits)

    result = annotate._collect_source_hits(
        {"q0": "ACGT", "q1": "ACGTACGT"},
        "Rfam",
        {},
        is_linear=False,
        true_lens={"q0": 100, "q1": 250},
    )

    # each hit's qlen is its own query's real length, not the searched length
    assert dict(zip(result["qseqid"], result["qlen"], strict=True)) == {
        "q0": 100,
        "q1": 250,
    }


def test_source_adapter_reports_missing_candidate_columns(monkeypatch):
    monkeypatch.setattr(
        annotate,
        "run_tool",
        lambda *args, **kwargs: pd.DataFrame({"sseqid": ["feature"]}),
    )

    with pytest.raises(ValueError, match="Source 'caller'.*qstart"):
        annotate._collect_source_hits(
            {"q0": "ACGT"},
            "caller",
            {},
            is_linear=True,
            true_lens={"q0": 4},
        )


def test_finalize_annotations_returns_only_canonical_columns():
    # a non-empty finalize must not leak the internal qseqid routing id; its columns
    # must match the canonical schema (and the empty path) exactly, in order.
    raw = pd.DataFrame(
        {
            "qseqid": ["q0"],
            "sseqid": ["feat"],
            "qstart": [1],
            "qend": [50],
            "sstart": [1],
            "send": [50],
            "sframe": [1],
            "evalue": [1e-9],
            "qseq": ["A" * 50],
            "length": [50],
            "slen": [50],
            "pident": [100.0],
            "qlen": [100],
            "db": ["snapgene"],
            "name": ["feat"],
            "blurb": [""],
            "type": ["CDS"],
            "priority": [1],
            "btop": ["50"],
            "structure": [""],
        }
    )

    finalized = annotate._finalize_annotations(raw, is_detailed=False, is_linear=True)

    assert list(finalized.columns) == annotate.ANNOTATION_COLUMNS
    assert list(annotate._empty_annotations().columns) == annotate.ANNOTATION_COLUMNS


def test_finalize_missing_type_uses_one_overlap_kind():
    rows = []
    for length, end in ((100, 100), (80, 80)):
        rows.append(
            {
                "qseqid": "q0",
                "sseqid": "unknown",
                "qstart": 1,
                "qend": end,
                "sstart": 1,
                "send": end,
                "sframe": 1,
                "evalue": 1e-20,
                "qseq": "A" * length,
                "length": length,
                "slen": 100,
                "pident": 100.0,
                "qlen": 200,
                "db": "custom",
                "name": None,
                "blurb": None,
                "type": None,
                "priority": 1,
                "btop": str(length),
                "structure": "",
            }
        )

    finalized = annotate._finalize_annotations(
        pd.DataFrame(rows), is_detailed=True, is_linear=True
    )

    assert len(finalized) == 1
    assert finalized.iloc[0]["type"] == "misc_feature"


def test_finalize_can_keep_raw_nested_fragments(monkeypatch):
    raw = pd.DataFrame(
        {
            "qseqid": ["q0"],
            "sseqid": ["feat"],
            "qstart": [1],
            "qend": [50],
            "sstart": [1],
            "send": [50],
            "sframe": [1],
            "evalue": [1e-20],
            "qseq": ["A" * 50],
            "length": [50],
            "slen": [100],
            "pident": [100.0],
            "qlen": [200],
            "db": ["snapgene"],
            "name": ["feat"],
            "blurb": [""],
            "type": ["promoter"],
            "priority": [1],
            "btop": ["50"],
            "structure": [""],
        }
    )
    called = []
    monkeypatch.setattr(
        annotate._nested,
        "suppress_nested_fragments",
        lambda hits: called.append(True) or hits,
    )

    annotate._finalize_annotations(
        raw,
        is_detailed=True,
        is_linear=True,
        apply_nested_policy=False,
    )

    assert called == []


def test_annotate_batch_splits_results_by_caller_key(monkeypatch):
    # the batch search keys each hit by an internal id (q0, q1, ... in input order)
    combined = pd.DataFrame(
        {
            "qseqid": ["q0", "q1", "q1"],
            "sseqid": ["a-feat", "b-feat-1", "b-feat-2"],
        }
    )
    monkeypatch.setattr(annotate, "_collect_hits_batch", lambda *_a, **_k: combined)
    # finalize is covered elsewhere; here just echo each sequence's slice of hits
    monkeypatch.setattr(
        annotate,
        "_finalize_annotations",
        lambda hits, is_detailed, is_linear, **_kwargs: hits.reset_index(drop=True),
    )

    results = annotate.annotate_batch({"alpha": "ACGT", "beta": "ACGTACGT"})

    # caller keys and input order are preserved, hits routed to the right sequence
    assert list(results) == ["alpha", "beta"]
    assert results["alpha"]["sseqid"].tolist() == ["a-feat"]
    assert results["beta"]["sseqid"].tolist() == ["b-feat-1", "b-feat-2"]


def test_annotate_batch_gives_empty_frame_to_unmatched_sequence(monkeypatch):
    combined = pd.DataFrame({"qseqid": ["q0"], "sseqid": ["only-alpha"]})
    monkeypatch.setattr(annotate, "_collect_hits_batch", lambda *_a, **_k: combined)
    monkeypatch.setattr(
        annotate,
        "_finalize_annotations",
        lambda hits, is_detailed, is_linear, **_kwargs: hits.reset_index(drop=True),
    )

    results = annotate.annotate_batch({"alpha": "ACGT", "beta": "TTTT"})

    assert results["alpha"]["sseqid"].tolist() == ["only-alpha"]
    assert results["beta"].empty
