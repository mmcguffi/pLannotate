#!/usr/bin/env python3
"""Decide whether a freshly built database bundle needs a new version.

Run this after rebuilding the databases. It compares the candidate bundle's
content fingerprint (which ignores ``build_date``) against the currently published
bundle and reports one of two outcomes:

* unchanged -- the content matches the published bundle, so keep the current
  ``DATABASE_VERSION`` and reuse/republish the existing bundle (exit code 0);
* changed -- the content differs, so bump ``DATABASE_VERSION`` and publish the new
  bundle (exit code 10).

The published bundle defaults to ``PLANNOTATE_DATABASE_URL`` or the canonical
release asset for the current version; pass ``--published`` to override with a URL
or local path.
"""

import argparse
import os
import shutil
import tempfile
from pathlib import Path
from urllib.request import urlopen

if __package__:
    from . import database_fingerprint as fp
else:
    import database_fingerprint as fp  # type: ignore[import-not-found,no-redef]

CHANGED_EXIT_CODE = 10


def _current_identity() -> tuple[int, str]:
    try:
        from plannotate._package_data import DATABASE_ASSET_NAME, DATABASE_VERSION
    except ImportError:
        DATABASE_VERSION = 2
        DATABASE_ASSET_NAME = "plannotate-databases-v2.tar.gz"
    return DATABASE_VERSION, DATABASE_ASSET_NAME


def _default_published_url(asset_name: str) -> str:
    return os.environ.get(
        "PLANNOTATE_DATABASE_URL",
        (
            "https://github.com/mmcguffi/pLannotate/releases/latest/download/"
            f"{asset_name}"
        ),
    )


def _fingerprint_of(location: str, workdir: Path) -> str:
    """Fingerprint a bundle given as a local path or a URL."""
    if "://" in location and not location.startswith("file://"):
        local = workdir / "published.tar.gz"
        with urlopen(location, timeout=60) as response, local.open("wb") as out:
            shutil.copyfileobj(response, out)
    else:
        local = Path(location.removeprefix("file://"))
    return fp.fingerprint_from_manifest(fp.load_manifest(local))


def main() -> int:
    current_version, asset_name = _current_identity()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "candidate",
        type=Path,
        help="freshly built bundle .tar.gz (or its database-manifest.json)",
    )
    parser.add_argument(
        "--published",
        default=_default_published_url(asset_name),
        help="URL or path to the currently published bundle",
    )
    args = parser.parse_args()

    candidate_fp = fp.fingerprint_from_manifest(fp.load_manifest(args.candidate))
    with tempfile.TemporaryDirectory(prefix="plannotate-dbcheck-") as tmp:
        published_fp = _fingerprint_of(args.published, Path(tmp))

    print(f"current DATABASE_VERSION: {current_version}")
    print(f"candidate fingerprint:    {candidate_fp}")
    print(f"published fingerprint:    {published_fp}")

    if candidate_fp == published_fp:
        print(
            f"unchanged: database content matches the published bundle; keep "
            f"DATABASE_VERSION = {current_version}."
        )
        return 0

    next_version = current_version + 1
    print(
        f"changed: database content differs. Bump DATABASE_VERSION to "
        f"{next_version} in plannotate/_package_data.py (asset becomes "
        f"plannotate-databases-v{next_version}.tar.gz) and publish the new bundle "
        f"(update PLANNOTATE_DATABASE_URL and PLANNOTATE_DATABASE_SHA256)."
    )
    return CHANGED_EXIT_CODE


if __name__ == "__main__":
    raise SystemExit(main())
