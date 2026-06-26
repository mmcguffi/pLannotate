"""Tests for packaged assets, configuration, and database installation."""

import io
import json
import tarfile
from pathlib import Path

import pandas as pd
import pytest

from plannotate import _package_data
from plannotate._tools import methods


def test_default_database_locations_follow_search_method():
    databases = _package_data.get_yaml(_package_data.get_yaml_path())

    assert Path(databases["snapgene"]["db_loc"]).parent.name == "BLAST_dbs"
    assert Path(databases["fpbase"]["db_loc"]).parent.name == "diamond_dbs"
    assert Path(databases["swissprot"]["db_loc"]).parent.name == "diamond_dbs"
    assert Path(databases["Rfam"]["db_loc"]).parent.name == "infernal_dbs"
    assert Path(databases["Rfam"]["clanin_loc"]).name == "Rfam.clanin"


@pytest.mark.parametrize(
    ("resource", "parts"),
    [
        (_package_data.get_image("icon.png"), ("data", "images", "icon.png")),
        (
            _package_data.get_template("blurb.html"),
            ("data", "templates", "blurb.html"),
        ),
        (_package_data.get_example_fastas(), ("data", "fastas")),
        (_package_data.get_yaml_path(), ("data", "data", "databases.yml")),
    ],
)
def test_packaged_resource_paths(resource, parts):
    assert resource.parts[-len(parts) :] == parts


@pytest.mark.integration
def test_downloaded_databases_are_complete():
    assert _package_data.databases_exist()


def test_database_config_rejects_invalid_priority(tmp_path):
    config = tmp_path / "databases.yml"
    config.write_text(
        "custom:\n"
        "  method: blastn\n"
        "  location: /tmp/database\n"
        "  priority: 0\n"
        "  details:\n"
        "    location: None\n"
    )

    with pytest.raises(ValueError, match="priority must be a positive integer"):
        _package_data.get_yaml(config)


def test_database_config_accepts_parameter_strings(tmp_path):
    config = tmp_path / "databases.yml"
    config.write_text(
        "custom:\n"
        "  method: blastn\n"
        "  location: /tmp/database\n"
        "  priority: 1\n"
        "  parameters: '-word_size 12'\n"
        "  details:\n"
        "    location: None\n"
    )

    assert _package_data.get_yaml(config)["custom"]["parameters"] == "-word_size 12"


def test_database_config_reports_missing_fields(tmp_path):
    config = tmp_path / "databases.yml"
    config.write_text("custom: {method: blastn}\n")

    with pytest.raises(ValueError, match="missing.*details"):
        _package_data.get_yaml(config)


def test_database_free_source_needs_no_location(tmp_path, monkeypatch):
    monkeypatch.setitem(
        methods.SEARCHERS,
        "caller",
        lambda *args: pd.DataFrame(),
    )
    monkeypatch.setitem(methods.DATABASE_DIRECTORIES, "caller", None)
    config = tmp_path / "sources.yml"
    config.write_text(
        "genes:\n"
        "  method: caller\n"
        "  location: null\n"
        "  priority: 1\n"
        "  details:\n"
        "    location: null\n"
    )

    source = _package_data.get_yaml(config)["genes"]

    assert "db_loc" not in source


def test_database_asset_url_tracks_package_version(monkeypatch):
    monkeypatch.delenv("PLANNOTATE_DATABASE_URL", raising=False)

    assert _package_data._database_asset_url().endswith(
        "/v2.0.0/plannotate-databases-v2.tar.gz"
    )


def test_database_asset_url_can_be_overridden(monkeypatch):
    monkeypatch.setenv("PLANNOTATE_DATABASE_URL", "https://example.test/db.tar.gz")

    assert _package_data._database_asset_url() == "https://example.test/db.tar.gz"


def test_validate_database_tree_reports_missing_files(tmp_path):
    with pytest.raises(ValueError, match="snapgene.nhr"):
        _package_data._validate_database_tree(tmp_path)


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "malicious.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("../outside.txt")
        payload = b"unsafe"
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))

    with pytest.raises(ValueError, match="Unsafe path"):
        _package_data._safe_extract(archive, tmp_path / "output")


def test_safe_extract_rejects_links(tmp_path):
    archive = tmp_path / "link.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("BLAST_dbs/link")
        info.type = tarfile.SYMTYPE
        info.linkname = "../../outside.txt"
        tar.addfile(info)

    with pytest.raises(ValueError, match="Links are not allowed"):
        _package_data._safe_extract(archive, tmp_path / "output")


def test_safe_extract_rejects_special_files(tmp_path):
    archive = tmp_path / "device.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("BLAST_dbs/device")
        info.type = tarfile.CHRTYPE
        tar.addfile(info)

    with pytest.raises(ValueError, match="Special files are not allowed"):
        _package_data._safe_extract(archive, tmp_path / "output")


def test_download_databases_installs_verified_bundle(tmp_path, monkeypatch):
    source = tmp_path / "source"
    payload_paths = set(_package_data.REQUIRED_DATABASE_FILES) - {
        _package_data.DATABASE_MANIFEST_NAME
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
            path: {"sha256": _package_data._sha256(source / path)}
            for path in payload_paths
        },
    }
    (source / _package_data.DATABASE_MANIFEST_NAME).write_text(json.dumps(manifest))

    archive = tmp_path / _package_data.DATABASE_ASSET_NAME
    with tarfile.open(archive, "w:gz") as tar:
        for relative_path in _package_data.REQUIRED_DATABASE_FILES:
            tar.add(source / relative_path, arcname=relative_path)

    destination = tmp_path / "installed"
    monkeypatch.setattr(_package_data, "get_data_directory", lambda: destination)
    monkeypatch.setenv("PLANNOTATE_DATABASE_URL", archive.as_uri())
    monkeypatch.setenv("PLANNOTATE_DATABASE_SHA256", _package_data._sha256(archive))

    _package_data.download_databases()

    _package_data._validate_database_tree(destination)
    _package_data._validate_database_manifest(destination)
    assert _package_data.databases_exist()
    assert _package_data.get_database_manifest()["schema_version"] == 1
