"""Tests for core allocation and source concurrency."""

import re

import pandas as pd
import pytest
from typer.testing import CliRunner

from plannotate import _concurrency, annotate
from plannotate.main import app


def test_sources_are_ordered_and_receive_thread_allocations():
    sources = {
        "first": {"method": "blastn"},
        "second": {"method": "infernal"},
        "third": {"method": "diamond"},
    }
    observed = []

    def fake_source(query, source_name, source_config, is_linear, threads):
        observed.append((query, source_name, source_config, is_linear, threads))
        if source_name == "second":
            return pd.DataFrame()
        return pd.DataFrame({"db": [source_name]})

    results = _concurrency.run_sources(
        fake_source,
        "ACGT",
        True,
        sources,
        5,
    )

    assert results[0].loc[0, "db"] == "first"
    assert results[1].empty
    assert results[2].loc[0, "db"] == "third"
    assert sorted(observed, key=lambda call: call[1]) == [
        ("ACGT", "first", sources["first"], True, 2),
        ("ACGT", "second", sources["second"], True, 2),
        ("ACGT", "third", sources["third"], True, 1),
    ]


@pytest.mark.parametrize(
    ("cores", "expected"),
    [
        (1, [1, 1, 1, 1]),
        (4, [1, 1, 1, 1]),
        (5, [2, 1, 1, 1]),
        (6, [2, 2, 1, 1]),
        (8, [2, 2, 2, 2]),
        (9, [3, 2, 2, 2]),
    ],
)
def test_thread_allocation_follows_configuration_order(cores, expected):
    sources = {
        "Rfam": {"method": "infernal"},
        "fpbase": {"method": "diamond"},
        "swissprot": {"method": "diamond"},
        "snapgene": {"method": "blastn"},
    }

    assert _concurrency.allocate_threads(sources, cores) == expected


@pytest.mark.parametrize(
    ("cores", "expected"),
    [
        (4, [1, 1, 1, 1]),
        (5, [2, 1, 1, 1]),
        (6, [2, 1, 2, 1]),
        (8, [3, 1, 3, 1]),
        (10, [4, 1, 3, 2]),
    ],
)
def test_thread_allocation_feeds_the_bottleneck_source(cores, expected):
    # spare cores go to whichever source is projected to finish last, so the cheap
    # source never collects threads it cannot use
    sources = {
        "Rfam": {"method": "infernal", "cost": 1.8},
        "fpbase": {"method": "diamond", "cost": 0.05},
        "swissprot": {"method": "diamond", "cost": 1.3},
        "snapgene": {"method": "blastn", "cost": 0.6},
    }

    assert _concurrency.allocate_threads(sources, cores) == expected


def test_thread_allocation_rejects_invalid_core_count():
    with pytest.raises(ValueError, match="cores must be at least 1"):
        _concurrency.allocate_threads({}, 0)


@pytest.mark.parametrize(
    ("parameters", "options", "canonical", "expected"),
    [
        (
            "-evalue 1 -num_threads 8",
            ("-num_threads",),
            "-num_threads",
            "-evalue 1 -num_threads 3",
        ),
        ("--fast --threads=8", ("--threads", "-p"), "--threads", "--fast --threads 3"),
        ("--cut_ga -p 8", ("--threads", "-p"), "--threads", "--cut_ga --threads 3"),
        ("--rfam --cpu=8", ("--cpu",), "--cpu", "--rfam --cpu 3"),
    ],
)
def test_allocated_thread_count_overrides_configured_parameters(
    parameters, options, canonical, expected
):
    assert (
        _concurrency.parameters_with_threads(parameters, options, canonical, 3)
        == expected
    )


def test_thread_option_requires_a_value():
    with pytest.raises(ValueError, match="requires a value"):
        _concurrency.parameters_with_threads(
            "--fast --threads", ("--threads",), "--threads", 2
        )


def test_annotate_passes_core_limit_to_collection(monkeypatch):
    observed = {}

    def fake_collection(sequence, linear, yaml_file, cores):
        observed["cores"] = cores
        return pd.DataFrame()

    monkeypatch.setattr(annotate, "_collect_hits", fake_collection)

    annotate.annotate("ACGT", cores=3)

    assert observed["cores"] == 3


def test_batch_help_documents_core_limit():
    result = CliRunner().invoke(app, ["batch", "--help"])
    help_text = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", result.stdout)

    assert result.exit_code == 0
    assert "--cores" in help_text
    assert "-j" in help_text
