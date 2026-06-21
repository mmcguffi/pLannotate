import io
import json
import tarfile
from pathlib import Path

import pytest

from plannotate import resources


def test_default_database_locations_follow_search_method():
    databases = resources.get_yaml(resources.get_yaml_path())

    assert Path(databases["snapgene"]["db_loc"]).parent.name == "BLAST_dbs"
    assert Path(databases["fpbase"]["db_loc"]).parent.name == "diamond_dbs"
    assert Path(databases["swissprot"]["db_loc"]).parent.name == "diamond_dbs"
    assert "infernal_dbs" in databases["Rfam"]["db_loc"]


def test_database_asset_url_tracks_package_version(monkeypatch):
    monkeypatch.delenv("PLANNOTATE_DATABASE_URL", raising=False)

    assert resources._database_asset_url().endswith(
        "/v2.0.0/plannotate-databases-v2.tar.gz"
    )


def test_database_asset_url_can_be_overridden(monkeypatch):
    monkeypatch.setenv("PLANNOTATE_DATABASE_URL", "https://example.test/db.tar.gz")

    assert resources._database_asset_url() == "https://example.test/db.tar.gz"


def test_validate_database_tree_reports_missing_files(tmp_path):
    with pytest.raises(ValueError, match="snapgene.nhr"):
        resources._validate_database_tree(tmp_path)


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "malicious.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("../outside.txt")
        payload = b"unsafe"
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))

    with pytest.raises(ValueError, match="Unsafe path"):
        resources._safe_extract(archive, tmp_path / "output")


def test_safe_extract_rejects_links(tmp_path):
    archive = tmp_path / "link.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("BLAST_dbs/link")
        info.type = tarfile.SYMTYPE
        info.linkname = "../../outside.txt"
        tar.addfile(info)

    with pytest.raises(ValueError, match="Links are not allowed"):
        resources._safe_extract(archive, tmp_path / "output")


def test_download_databases_installs_verified_bundle(tmp_path, monkeypatch):
    source = tmp_path / "source"
    payload_paths = set(resources.REQUIRED_DATABASE_FILES) - {
        resources.DATABASE_MANIFEST_NAME
    }
    for relative_path in payload_paths:
        database_file = source / relative_path
        database_file.parent.mkdir(parents=True, exist_ok=True)
        database_file.write_bytes(relative_path.encode())
    manifest = {
        "schema_version": 1,
        "bundle": "plannotate-databases-v2",
        "build_date": "2026-06-20",
        "databases": {},
        "files": {
            path: {"sha256": resources._sha256(source / path)} for path in payload_paths
        },
    }
    (source / resources.DATABASE_MANIFEST_NAME).write_text(json.dumps(manifest))

    archive = tmp_path / resources.DATABASE_ASSET_NAME
    with tarfile.open(archive, "w:gz") as tar:
        for relative_path in resources.REQUIRED_DATABASE_FILES:
            tar.add(source / relative_path, arcname=relative_path)

    destination = tmp_path / "installed"
    monkeypatch.setattr(resources, "get_data_directory", lambda: destination)
    monkeypatch.setenv("PLANNOTATE_DATABASE_URL", archive.as_uri())
    monkeypatch.setenv("PLANNOTATE_DATABASE_SHA256", resources._sha256(archive))

    resources.download_databases()

    resources._validate_database_tree(destination)
    resources._validate_database_manifest(destination)
    assert resources.databases_exist()
    assert resources.get_database_manifest()["schema_version"] == 1
