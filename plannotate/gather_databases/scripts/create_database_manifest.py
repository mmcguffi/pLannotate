#!/usr/bin/env python3
"""Create the versioned manifest installed with a pLannotate database bundle."""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

if __package__:
    from . import package_database_bundle as bundle
else:
    import package_database_bundle as bundle  # type: ignore[import-not-found,no-redef]

DATABASE_FILES = bundle.DATABASE_FILES


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _build_date() -> str:
    source_date_epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if source_date_epoch:
        instant = datetime.fromtimestamp(int(source_date_epoch), tz=timezone.utc)
    else:
        instant = datetime.now(timezone.utc)
    return instant.date().isoformat()


def create_manifest(
    source: Path,
    output: Path,
    versions: dict[str, str],
    build_date: str | None = None,
) -> dict:
    """Create a manifest with source versions and payload checksums."""
    missing = [path for path in DATABASE_FILES if not (source / path).is_file()]
    if missing:
        raise ValueError("Cannot manifest incomplete databases: " + ", ".join(missing))

    manifest = {
        "schema_version": 1,
        "bundle": "plannotate-databases-v2",
        "build_date": build_date or _build_date(),
        "databases": {
            name: {"version": version.strip()}
            for name, version in sorted(versions.items())
        },
        "files": {
            relative_path: {"sha256": _sha256(source / relative_path)}
            for relative_path in DATABASE_FILES
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rfam-version-file", type=Path, required=True)
    parser.add_argument("--snapgene-version-file", type=Path, required=True)
    parser.add_argument("--fpbase-version-file", type=Path, required=True)
    parser.add_argument("--swissprot-version-file", type=Path, required=True)
    parser.add_argument("--build-date")
    args = parser.parse_args()
    versions = {
        "Rfam": args.rfam_version_file.read_text(),
        "snapgene": args.snapgene_version_file.read_text(),
        "fpbase": args.fpbase_version_file.read_text(),
        "swissprot": args.swissprot_version_file.read_text(),
    }
    create_manifest(args.source, args.output, versions, args.build_date)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
