import json
import subprocess
from pathlib import Path

import re

import pandas as pd
import pytest
from typer.testing import CliRunner

from plannotate import _annotation_worker, annotate, search
from plannotate.main import app


def test_worker_runs_selected_database(monkeypatch, tmp_path):
    query_file = tmp_path / "query.txt"
    query_file.write_text("ACGT")
    yaml_file = tmp_path / "databases.yml"
    yaml_file.touch()
    output_file = tmp_path / "result.pkl"
    databases = {
        "first": {"method": "blastn"},
        "second": {"method": "infernal"},
    }
    monkeypatch.setattr(_annotation_worker.rsc, "get_yaml", lambda _: databases)

    observed = {}

    def fake_search(
        query, database_name, database_config, yaml_path, is_linear, threads
    ):
        observed.update(
            query=query,
            database_name=database_name,
            database_config=database_config,
            yaml_path=yaml_path,
            is_linear=is_linear,
            threads=threads,
        )
        return pd.DataFrame({"db": [database_name]})

    monkeypatch.setattr(_annotation_worker, "_search_single_database", fake_search)

    _annotation_worker.run_database_search(
        query_file=query_file,
        yaml_file=yaml_file,
        database_index=1,
        output_file=output_file,
        linear=True,
        threads=3,
    )

    assert observed == {
        "query": "ACGT",
        "database_name": "second",
        "database_config": databases["second"],
        "yaml_path": yaml_file,
        "is_linear": True,
        "threads": 3,
    }
    assert pd.read_pickle(output_file).to_dict("records") == [{"db": "second"}]


def test_workflow_passes_cores_and_collects_results(monkeypatch, tmp_path):
    yaml_file = tmp_path / "databases.yml"
    yaml_file.touch()
    monkeypatch.setattr(search.shutil, "which", lambda _: "/usr/bin/snakemake")
    monkeypatch.setattr(
        search.rsc,
        "get_yaml",
        lambda _: {
            "rfam": {"method": "infernal"},
            "diamond": {"method": "diamond"},
            "blast": {"method": "blastn"},
        },
    )
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        observed["kwargs"] = kwargs
        config_path = Path(command[command.index("--configfile") + 1])
        config = json.loads(config_path.read_text())
        observed["config"] = config
        work_dir = Path(command[command.index("--directory") + 1])
        results_dir = work_dir / "results"
        results_dir.mkdir()
        pd.DataFrame({"db": ["first"]}).to_pickle(results_dir / "0.pkl")
        pd.DataFrame().to_pickle(results_dir / "1.pkl")
        pd.DataFrame({"db": ["third"]}).to_pickle(results_dir / "2.pkl")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(search.subprocess, "run", fake_run)

    results = search._run_annotation_workflow(
        query_sequence="ACGT",
        is_linear=True,
        yaml_file=yaml_file,
        cores=2,
    )

    assert observed["command"][observed["command"].index("--cores") + 1] == "2"
    assert observed["config"]["database_count"] == 3
    assert observed["config"]["linear"] is True
    assert observed["config"]["thread_allocations"] == [1, 1, 1]
    assert [result.loc[0, "db"] for result in results] == ["first", "third"]
    assert observed["kwargs"]["capture_output"] is True


def test_workflow_reports_scheduler_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(search.shutil, "which", lambda _: "/usr/bin/snakemake")
    monkeypatch.setattr(
        search.rsc,
        "get_yaml",
        lambda _: {"rfam": {"method": "infernal"}},
    )
    monkeypatch.setattr(
        search.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 1, "", "database job failed"
        ),
    )

    with pytest.raises(RuntimeError, match="database job failed"):
        search._run_annotation_workflow(
            query_sequence="ACGT",
            is_linear=False,
            yaml_file=tmp_path / "databases.yml",
            cores=1,
        )


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
def test_thread_allocation_prioritizes_rfam_then_diamond(cores, expected):
    databases = {
        "Rfam": {"method": "infernal"},
        "fpbase": {"method": "diamond"},
        "swissprot": {"method": "diamond"},
        "snapgene": {"method": "blastn"},
    }

    assert search._allocate_search_threads(databases, cores) == expected


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
def test_scheduler_thread_count_overrides_database_parameters(
    parameters, options, canonical, expected
):
    assert (
        search._parameters_with_threads(parameters, options, canonical, 3) == expected
    )


def test_annotate_passes_core_limit_to_search(monkeypatch):
    observed = {}

    def fake_search(sequence, linear, yaml_file, cores):
        observed["cores"] = cores
        return pd.DataFrame()

    monkeypatch.setattr(annotate, "search_all_databases", fake_search)

    annotate.annotate("ACGT", cores=3)

    assert observed["cores"] == 3


def test_batch_help_documents_core_limit():
    result = CliRunner().invoke(app, ["batch", "--help"])
    help_text = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", result.stdout)

    assert result.exit_code == 0
    assert "--cores" in help_text
    assert "-j" in help_text
