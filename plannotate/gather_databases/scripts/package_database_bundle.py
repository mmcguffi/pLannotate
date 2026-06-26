#!/usr/bin/env python3
"""Create a deterministic pLannotate 2.x database release bundle."""

import argparse
import gzip
import hashlib
import tarfile
from pathlib import Path

DATABASE_FILES = (
    "BLAST_dbs/snapgene.nhr",
    "BLAST_dbs/snapgene.nin",
    "BLAST_dbs/snapgene.nsq",
    "BLAST_dbs/snapgene.ndb",
    "BLAST_dbs/snapgene.nog",
    "BLAST_dbs/snapgene.nos",
    "BLAST_dbs/snapgene.not",
    "BLAST_dbs/snapgene.ntf",
    "BLAST_dbs/snapgene.nto",
    "BLAST_dbs/descriptions.db",
    "diamond_dbs/fpbase.dmnd",
    "diamond_dbs/swissprot.dmnd",
    "diamond_dbs/descriptions.db",
    "infernal_dbs/Rfam.cm",
    "infernal_dbs/Rfam.clanin",
    "infernal_dbs/Rfam.cm.i1f",
    "infernal_dbs/Rfam.cm.i1i",
    "infernal_dbs/Rfam.cm.i1m",
    "infernal_dbs/Rfam.cm.i1p",
)
MANIFEST_FILE = "database-manifest.json"
REQUIRED_FILES = DATABASE_FILES + (MANIFEST_FILE,)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_tar_info(info: tarfile.TarInfo) -> tarfile.TarInfo:
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    return info


def package_database_bundle(source: Path, output: Path) -> str:
    missing = [path for path in REQUIRED_FILES if not (source / path).is_file()]
    if missing:
        raise ValueError(
            "Cannot package incomplete database tree; missing: " + ", ".join(missing)
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w") as tar:
                for relative_path in REQUIRED_FILES:
                    source_path = source / relative_path
                    tar.add(
                        source_path,
                        arcname=relative_path,
                        recursive=False,
                        filter=_normalize_tar_info,
                    )

    checksum = _sha256(output)
    output.with_name(f"{output.name}.sha256").write_text(f"{checksum}  {output.name}\n")
    return checksum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path("gathered_data"))
    parser.add_argument(
        "--output", type=Path, default=Path("plannotate-databases-v2.tar.gz")
    )
    args = parser.parse_args()
    checksum = package_database_bundle(args.source, args.output)
    print(f"Created {args.output} ({checksum})")


if __name__ == "__main__":
    main()
