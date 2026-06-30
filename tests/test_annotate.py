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
    row = {column: 0 for column in annotate.ADAPTER_COLUMNS}
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
                qstart=941, qend=1000, sstart=1, send=60, length=60, qseq="RIGHT"
            ),
            _seam_fragment(
                qstart=1, qend=40, sstart=61, send=100, length=40, qseq="LEFT"
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


def test_fast_mode_does_not_double_circular_query(monkeypatch):
    captured = {}

    def fake_run(search_query, config, threads):
        captured["query"] = search_query
        return pd.DataFrame({column: [1] for column in annotate.ADAPTER_COLUMNS})

    monkeypatch.setattr(annotate, "run_tool", fake_run)
    monkeypatch.setattr(annotate, "_enrich_hits", lambda hits, name, config: hits)
    monkeypatch.setattr(annotate, "_stitch_seam_hits", lambda hits: hits)

    query = "ACGT" * 25  # length 100
    annotate._collect_source_hits(query, "snapgene", {}, is_linear=False, fast=True)

    # fast mode searches a single copy; seam features are stitched, not doubled
    assert captured["query"] == query


def test_missing_swissprot_description_has_no_priority_penalty():
    assert annotate._existence_level_priority(None) == 0


def test_circular_search_query_doubles_sequence():
    query = "ACGT" * 25  # length 100

    # circular queries are fully doubled so origin-spanning features are recovered
    assert annotate._circular_search_query(query) == query + query


def test_collect_source_hits_doubles_query_and_records_true_length(monkeypatch):
    captured = {}

    def fake_run(search_query, config, threads):
        captured["query"] = search_query
        return pd.DataFrame({column: [1] for column in annotate.ADAPTER_COLUMNS})

    monkeypatch.setattr(annotate, "run_tool", fake_run)
    monkeypatch.setattr(annotate, "_enrich_hits", lambda hits, name, config: hits)

    query = "ACGT" * 25  # length 100
    result = annotate._collect_source_hits(query, "Rfam", {}, is_linear=False)

    # circular -> fully doubled (lossless)
    assert captured["query"] == query + query
    # qlen reflects the real sequence length, not the doubled search length
    assert result["qlen"].tolist() == [100]


def test_collect_source_hits_does_not_double_linear_query(monkeypatch):
    captured = {}

    def fake_run(search_query, config, threads):
        captured["query"] = search_query
        return pd.DataFrame({column: [1] for column in annotate.ADAPTER_COLUMNS})

    monkeypatch.setattr(annotate, "run_tool", fake_run)
    monkeypatch.setattr(annotate, "_enrich_hits", lambda hits, name, config: hits)

    query = "ACGT" * 25  # length 100
    annotate._collect_source_hits(query, "Rfam", {}, is_linear=True)

    assert captured["query"] == query


def test_source_adapter_reports_missing_candidate_columns(monkeypatch):
    monkeypatch.setattr(
        annotate,
        "run_tool",
        lambda *args, **kwargs: pd.DataFrame({"sseqid": ["feature"]}),
    )

    with pytest.raises(ValueError, match="Source 'caller'.*qstart"):
        annotate._collect_source_hits(
            "ACGT",
            "caller",
            {},
            is_linear=True,
        )
