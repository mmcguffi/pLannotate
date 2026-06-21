import sqlite3
from io import StringIO

from plannotate.gather_databases.fpbase.gather_fpbase import (
    write_fasta_output,
    write_tsv_output,
)
from plannotate.gather_databases.scripts.combine_tsv import combine_files
from plannotate.gather_databases.scripts.create_database_manifest import create_manifest
from plannotate.gather_databases.scripts.package_database_bundle import DATABASE_FILES
from plannotate.gather_databases.scripts.csv_to_sqlite import create_sqlite_from_csv


def test_database_creator_normalizes_named_columns(tmp_path):
    source = tmp_path / "snapgene.csv"
    source.write_text(
        "sseqid,Feature,Type,Description\n"
        "P_element_5'_end,P element 5' end,misc_feature,terminal sequence\n"
    )
    destination = tmp_path / "descriptions.db"

    assert create_sqlite_from_csv(
        source,
        destination,
        "snapgene",
        column_renames={
            "Feature": "name",
            "Type": "type",
            "Description": "blurb",
        },
    )

    with sqlite3.connect(destination) as connection:
        row = connection.execute(
            "SELECT sseqid, name, type, blurb FROM snapgene"
        ).fetchone()
    assert row == (
        "P_element_5'_end",
        "P element 5' end",
        "misc_feature",
        "terminal sequence",
    )


def test_fpbase_metadata_id_matches_fasta_header():
    proteins = [
        {
            "id": "1",
            "name": "Example FP",
            "slug": "example-fp",
            "seq": "MTEST",
        }
    ]
    metadata = StringIO()
    fasta = StringIO()

    write_tsv_output(proteins, metadata)
    write_fasta_output(proteins, fasta)

    assert metadata.getvalue().split("\t", 1)[0] == "example-fp"
    assert fasta.getvalue().startswith(">example-fp\n")


def test_combine_files_preserves_order_and_inputs(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    output = tmp_path / "combined.tsv"
    first.write_text("one\n")
    second.write_text("two\n")

    combine_files([first, second], output)

    assert output.read_text() == "one\ntwo\n"
    assert first.read_text() == "one\n"
    assert second.read_text() == "two\n"


def test_database_manifest_records_versions_and_payload_checksums(tmp_path):
    for relative_path in DATABASE_FILES:
        path = tmp_path / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(relative_path.encode())
    output = tmp_path / "database-manifest.json"

    manifest = create_manifest(
        tmp_path,
        output,
        {"Rfam": "release 15.1", "swissprot": "release 2026_02"},
        build_date="2026-06-20",
    )

    assert manifest["build_date"] == "2026-06-20"
    assert manifest["databases"]["Rfam"]["version"] == "release 15.1"
    assert set(manifest["files"]) == set(DATABASE_FILES)
