"""Tests for the database content fingerprint used to gate database versions."""

import importlib.util
import io
import json
import tarfile
from pathlib import Path

_SCRIPT = Path("plannotate/gather_databases/scripts/database_fingerprint.py")
_spec = importlib.util.spec_from_file_location("database_fingerprint", _SCRIPT)
assert _spec is not None and _spec.loader is not None
fingerprint = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fingerprint)


def _manifest(files, build_date="2026-06-28"):
    return {
        "schema_version": 1,
        "bundle": "plannotate-databases-v2",
        "database_version": 2,
        "build_date": build_date,
        "databases": {"Rfam": {"version": "14"}},
        "files": {name: {"sha256": digest} for name, digest in files.items()},
    }


def test_fingerprint_ignores_build_date():
    first = _manifest({"a": "aaa", "b": "bbb"}, build_date="2026-01-01")
    second = _manifest({"a": "aaa", "b": "bbb"}, build_date="2030-12-31")

    assert fingerprint.fingerprint_from_manifest(
        first
    ) == fingerprint.fingerprint_from_manifest(second)


def test_fingerprint_changes_with_content():
    base = _manifest({"a": "aaa", "b": "bbb"})
    changed = _manifest({"a": "aaa", "b": "ccc"})

    assert fingerprint.fingerprint_from_manifest(
        base
    ) != fingerprint.fingerprint_from_manifest(changed)


def test_fingerprint_independent_of_file_order():
    forward = _manifest({"a": "aaa", "b": "bbb"})
    reordered = {"files": {"b": {"sha256": "bbb"}, "a": {"sha256": "aaa"}}}

    assert fingerprint.fingerprint_from_manifest(
        forward
    ) == fingerprint.fingerprint_from_manifest(reordered)


def test_load_manifest_from_tarball_matches_json(tmp_path):
    manifest = _manifest({"a": "aaa", "b": "bbb"})
    manifest_path = tmp_path / "database-manifest.json"
    manifest_path.write_text(json.dumps(manifest))

    archive = tmp_path / "bundle.tar.gz"
    payload = json.dumps(manifest).encode()
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("database-manifest.json")
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))

    from_json = fingerprint.fingerprint_from_manifest(
        fingerprint.load_manifest(manifest_path)
    )
    from_tarball = fingerprint.fingerprint_from_manifest(
        fingerprint.load_manifest(archive)
    )

    assert from_json == from_tarball
