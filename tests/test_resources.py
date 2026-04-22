import os
import os.path as op
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
