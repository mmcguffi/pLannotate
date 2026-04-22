import os
import os.path as op
import subprocess
from pathlib import Path

import pytest

from plannotate import resources


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


def test_get_image():
    name = "icon.png"
    path = ("plannotate", "data", "images", name)
    assert resources.get_image(name) == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_template():
    name = "blurb.html"
    path = ("plannotate", "data", "templates", name)
    assert resources.get_template(name) == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_fasta():
    # TODO: this is a hack -- fix it
    current_dir = Path((os.path.abspath(__file__))).parent.parent

    fasta_loc = str(current_dir / "plannotate" / "data" / "fastas")
    assert resources.get_example_fastas() == fasta_loc


def test_get_yaml_path():
    path = ("plannotate", "data", "data", "databases.yml")
    assert resources.get_yaml_path() == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_yaml():
    yaml = resources.get_yaml(resources.get_yaml_path())

    assert isinstance(yaml, dict)
    assert len(yaml) > 0

    first_key = list(yaml.keys())[0]
    expected_fields = set(
        ("version", "method", "location", "priority", "parameters", "details", "db_loc")
    )
    assert set(yaml[first_key].keys()) == expected_fields


@pytest.mark.integration
def test_databases_exist():
    assert resources.databases_exist() is True


def test_download_databases_runs_expected_commands(monkeypatch):
    calls = []

    def fake_run(command, check=False):
        calls.append((command, check))
        return subprocess.CompletedProcess(command, 0)

    def fake_exists(path):
        return path.endswith("BLAST_dbs.tar.gz")

    monkeypatch.setattr(resources.subprocess, "run", fake_run)
    monkeypatch.setattr(resources.os.path, "exists", fake_exists)

    resources.download_databases()

    archive_path = f"{resources.ROOT_DIR}/data/BLAST_dbs.tar.gz"
    assert calls == [
        (
            [
                "curl",
                "-L",
                "-o",
                archive_path,
                "https://github.com/mmcguffi/pLannotate/releases/download/v1.2.0/BLAST_dbs.tar.gz",
            ],
            True,
        ),
        (
            ["tar", "-xzf", archive_path, "-C", f"{resources.ROOT_DIR}/data/"],
            True,
        ),
        (["rm", archive_path], False),
    ]


def test_download_databases_exits_if_archive_missing(monkeypatch):
    def fake_run(command, check=False):
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(resources.subprocess, "run", fake_run)
    monkeypatch.setattr(resources.os.path, "exists", lambda path: False)

    with pytest.raises(SystemExit):
        resources.download_databases()
