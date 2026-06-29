#!/usr/bin/env python3
"""Compute a content fingerprint for a pLannotate database bundle.

The fingerprint identifies the database *content* independently of when it was
built. It hashes the per-file checksums recorded in ``database-manifest.json``
and deliberately ignores volatile metadata such as ``build_date``, so two builds
of identical databases produce the same fingerprint. The build helper uses it to
decide whether a freshly built bundle differs from the published one.

Accepts either the manifest JSON itself or a bundle ``.tar.gz`` (the manifest is
read from inside the archive). Prints the fingerprint to stdout.
"""

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

MANIFEST_FILE = "database-manifest.json"


def fingerprint_from_manifest(manifest: dict) -> str:
    """Return a stable sha256 over the manifest's per-file checksums."""
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("Manifest does not contain a 'files' checksum map")
    lines = []
    for relative_path in sorted(files):
        entry = files[relative_path]
        if not isinstance(entry, dict) or "sha256" not in entry:
            raise ValueError(f"Invalid manifest entry for {relative_path}")
        lines.append(f"{relative_path}:{entry['sha256']}")
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_manifest(path: Path) -> dict:
    if path.suffix == ".json":
        return json.loads(path.read_text())
    # treat anything else as a bundle archive and read the manifest member
    with tarfile.open(path, "r:*") as tar:
        try:
            member = tar.getmember(MANIFEST_FILE)
        except KeyError as exc:
            raise ValueError(f"{path} does not contain {MANIFEST_FILE}") from exc
        handle = tar.extractfile(member)
        if handle is None:
            raise ValueError(f"Could not read {MANIFEST_FILE} from {path}")
        return json.loads(handle.read())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        type=Path,
        help="database-manifest.json or a plannotate-databases bundle .tar.gz",
    )
    args = parser.parse_args()
    print(fingerprint_from_manifest(load_manifest(args.path)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
