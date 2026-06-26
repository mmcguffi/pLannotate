"""Tests for deterministic database-bundle packaging."""

import subprocess
import sys
import tarfile
from pathlib import Path

from plannotate._package_data import REQUIRED_DATABASE_FILES


def test_database_bundle_is_complete_and_reproducible(tmp_path):
    source = tmp_path / "source"
    for relative_path in REQUIRED_DATABASE_FILES:
        database_file = source / relative_path
        database_file.parent.mkdir(parents=True, exist_ok=True)
        database_file.write_bytes(f"fixture:{relative_path}".encode())

    script = Path("plannotate/gather_databases/scripts/package_database_bundle.py")
    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"
    for output in (first, second):
        subprocess.run(
            [
                sys.executable,
                str(script),
                "--source",
                str(source),
                "--output",
                str(output),
            ],
            check=True,
        )

    assert first.read_bytes() == second.read_bytes()
    assert first.with_name(f"{first.name}.sha256").is_file()
    with tarfile.open(first, "r:gz") as archive:
        assert tuple(archive.getnames()) == REQUIRED_DATABASE_FILES
