"""Tests for SQLite feature-description lookup."""

import sqlite3
from contextlib import closing
from pathlib import Path

import pandas as pd

from plannotate._sqlite import (
    get_descriptions_db_path,
    load_descriptions_from_sqlite,
)


def test_descriptions_db_path_is_per_source_by_default():
    config = {
        "db_loc": "/data/diamond_dbs/swissprot",
        "details": {"location": "Default"},
    }

    assert get_descriptions_db_path("swissprot", config) == Path(
        "/data/diamond_dbs/swissprot.db"
    )


def test_descriptions_db_path_honors_explicit_db_override():
    config = {
        "db_loc": "/data/diamond_dbs/swissprot",
        "details": {"location": "/custom/my_descriptions.db"},
    }

    assert get_descriptions_db_path("swissprot", config) == Path(
        "/custom/my_descriptions.db"
    )


def test_load_descriptions_binds_identifiers_with_apostrophes(tmp_path):
    database_dir = tmp_path / "BLAST_dbs"
    database_dir.mkdir()
    database_path = database_dir / "snapgene.db"
    expected = pd.DataFrame(
        [
            {
                "sseqid": "P_element_5'_end",
                "name": "P element 5' end",
                "type": "misc_feature",
                "blurb": "terminal sequence",
            }
        ]
    )
    with closing(sqlite3.connect(database_path)) as connection:
        expected.to_sql("snapgene", connection, index=False)

    actual = load_descriptions_from_sqlite(
        "snapgene",
        {"P_element_5'_end"},
        {"db_loc": str(database_dir / "snapgene")},
    )

    pd.testing.assert_frame_equal(actual, expected)
