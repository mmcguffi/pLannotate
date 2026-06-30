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
        (1, (1, [1, 1, 1, 1])),
        (4, (4, [1, 1, 1, 1])),
        (5, (4, [2, 1, 1, 1])),
        (6, (4, [2, 2, 1, 1])),
        (8, (4, [2, 2, 2, 2])),
        (9, (4, [3, 2, 2, 2])),
    ],
)
def test_equal_cost_sources_get_one_balanced_thread_each(cores, expected):
    # with equal costs and enough cores, every source runs in its own lane and
    # spare threads spread evenly -- the original balanced behaviour
    sources = {
        "Rfam": {"method": "infernal"},
        "fpbase": {"method": "diamond"},
        "swissprot": {"method": "diamond"},
        "snapgene": {"method": "blastn"},
    }

    assert _concurrency.plan_concurrency(sources, cores) == expected


@pytest.mark.parametrize(
    ("cores", "expected"),
    [
        # at/below the source count, fewer lanes let the scaling bottleneck (Rfam)
        # take extra threads instead of being pinned to one
        (4, (3, [2, 1, 1, 1])),
        (5, (3, [2, 1, 2, 1])),
        (6, (3, [3, 1, 2, 1])),
        (8, (4, [3, 1, 3, 1])),
        (10, (3, [5, 1, 3, 2])),
    ],
)
def test_lane_planner_feeds_the_scaling_bottleneck(cores, expected):
    sources = {
        "Rfam": {"method": "infernal", "cost": 1.8},
        "fpbase": {"method": "diamond", "cost": 0.05},
        "swissprot": {"method": "diamond", "cost": 1.3},
        "snapgene": {"method": "blastn", "cost": 0.6},
    }

    assert _concurrency.plan_concurrency(sources, cores) == expected


def test_lane_plan_never_exceeds_the_core_budget():
    # the heaviest `workers` allocations may run at once; their sum must fit cores
    sources = {
        "Rfam": {"method": "infernal", "cost": 1.8},
        "fpbase": {"method": "diamond", "cost": 0.05},
        "swissprot": {"method": "diamond", "cost": 1.3},
        "snapgene": {"method": "blastn", "cost": 0.6},
    }
    for cores in range(1, 13):
        workers, threads = _concurrency.plan_concurrency(sources, cores)
        peak = sum(sorted(threads, reverse=True)[:workers])
        assert peak <= cores
        assert min(threads) >= 1


def test_thread_allocation_rejects_invalid_core_count():
    with pytest.raises(ValueError, match="cores must be at least 1"):
        _concurrency.plan_concurrency({}, 0)


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

    def fake_collection(sequence, linear, yaml_file, cores, fast=False):
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
