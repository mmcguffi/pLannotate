import sqlite3
from contextlib import closing

import pandas as pd

from plannotate.sqlite_utils import (
    check_sqlite_database,
    load_descriptions_from_sqlite,
    query_description_by_sseqid,
)


def test_load_descriptions_binds_identifiers_with_apostrophes(tmp_path):
    database_dir = tmp_path / "BLAST_dbs"
    database_dir.mkdir()
    database_path = database_dir / "descriptions.db"
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
        tmp_path,
        {"P_element_5'_end"},
        {"method": "blastn"},
    )

    pd.testing.assert_frame_equal(actual, expected)
    assert query_description_by_sseqid(
        "snapgene",
        tmp_path,
        "P_element_5'_end",
        {"method": "blastn"},
    ) == {
        "name": "P element 5' end",
        "type": "misc_feature",
        "blurb": "terminal sequence",
    }
    assert check_sqlite_database("snapgene", tmp_path, {"method": "blastn"})
