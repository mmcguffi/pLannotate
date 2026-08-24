#!/usr/bin/env python3
"""Audit accession pins and source geometry in packaged curation tables.

The curated tables in ``plannotate/data/data`` are keyed on ``(db, sseqid)``: every
claim is frozen to the exact packaged source records a human checked. That precision
is what stops a Swiss-Prot catalase from inheriting a chloramphenicol-resistance
claim, but it couples the tables to one database bundle. When ``plannotate setupdb``
installs a new bundle the pins can drift in two directions:

* a pinned accession disappears -- the row goes silently dead, annotating nothing;
* a new record appears under a curated name -- it is silently *not* annotated.

Both fail closed, so no test catches either. This script makes the drift visible.

``EXCLUSIONS`` records the records that were examined and deliberately left unpinned,
with the reason. Anything unpinned that is *not* listed there is new and needs a human
to adjudicate it -- deciding whether a record really carries the claim is a judgement
about biology, so this script reports candidates and never edits the tables.

    python tools/curation_pins.py report    # human-readable drift report
    python tools/curation_pins.py check     # exit nonzero on any drift
"""

import argparse
import sqlite3
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "plannotate" / "data"
CURATED_TABLES = ("selection_markers.csv", "ori_copy_number.csv")
PIN_ONLY_TABLES = ("feature_suppressions.csv", "fragment_suppression_regions.csv")
NESTED_OVERRIDE_TABLE = "nested_feature_overrides.csv"
COMPOSITE_REGION_TABLE = "composite_reference_regions.csv"

# Description databases keyed by the source name used in the tables' ``db`` column.
SOURCE_DATABASES = {
    "snapgene": (DATA / "BLAST_dbs" / "snapgene.db", "snapgene"),
    "swissprot": (DATA / "diamond_dbs" / "swissprot.db", "swissprot"),
    "fpbase": (DATA / "diamond_dbs" / "fpbase.db", "fpbase"),
}

# Records that share a curated row's name but were deliberately NOT pinned, and why.
# This is the durable home for those judgements: without it a later maintainer cannot
# tell a rejected record from one that simply has not been reviewed yet.
EXCLUSIONS: dict[tuple[str, str], str] = {
    ("swissprot", "Q9PT92"): "cat: CATA_DANRE, a catalase rather than a CAT",
    ("swissprot", "Q9PWF7"): "cat: CATA_GLARU, a catalase rather than a CAT",
    ("swissprot", "Q27710"): "cat: CATA_ONCVE, a catalase rather than a CAT",
    (
        "swissprot",
        "Q04938",
    ): "sacB: PTSB_LACLL, a PTS system protein, not a levansucrase",
    (
        "swissprot",
        "O68215",
    ): "sacB: SACB2_NEIMD, Neisseria capsule transport, not a levansucrase",
    (
        "swissprot",
        "Q83U59",
    ): "sacB: SACB3_NEIMD, Neisseria capsule transport, not a levansucrase",
    (
        "swissprot",
        "Q84CZ9",
    ): "sacB: SACB4_NEIMD, Neisseria capsule transport, not a levansucrase",
    (
        "swissprot",
        "Q84D00",
    ): "sacB: SACB5_NEIMD, Neisseria capsule transport, not a levansucrase",
    (
        "swissprot",
        "Q9JWW8",
    ): "sacB: SACB_NEIMA, Neisseria capsule transport, not a levansucrase",
    (
        "swissprot",
        "P32747",
    ): "ura3: PYRD (dihydroorotate dehydrogenase) is not the 5-FOA target",
    (
        "swissprot",
        "Q4VR96",
    ): "aadA: blurb states it does not confer streptomycin/spectinomycin resistance",
    ("swissprot", "Q2GPT7"): "ccdB: an unrelated fungal biosynthetic enzyme",
    ("swissprot", "P45709"): "ccdB: blurb does not support the gyrase-poison claim",
    ("swissprot", "Q7X2H8"): "codA: choline oxidase, not cytosine deaminase",
    # `leu2` carries no row at all: its only Swiss-Prot record is LEUC_SCHPO,
    # isopropylmalate dehydratase (LEU1 activity), not LEU2.
    ("swissprot", "O14289"): "leu2: LEUC_SCHPO is LEU1 activity; the row was dropped",
}


def _load_source(database: str) -> pd.DataFrame:
    path, table = SOURCE_DATABASES[database]
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} is missing; run 'plannotate setupdb' before auditing pins"
        )
    with sqlite3.connect(path) as connection:
        return pd.read_sql_query(
            f"select sseqid, name, blurb from {table}", connection
        ).fillna("")


def _curated_rows(filename: str) -> pd.DataFrame:
    return pd.read_csv(DATA / "data" / filename, dtype=str).fillna("")


@lru_cache(maxsize=None)
def _snapgene_sequence(accession: str) -> str:
    """Extract one installed SnapGene record for composite-region validation."""
    result = subprocess.run(
        [
            "blastdbcmd",
            "-db",
            str(DATA / "BLAST_dbs" / "snapgene"),
            "-entry",
            accession,
            "-outfmt",
            "%s",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return "".join(result.stdout.split()).upper()


def _pinned_accessions(row: pd.Series) -> set[str]:
    return {part.strip() for part in row["sseqid"].split(";") if part.strip()}


def audit() -> tuple[list[str], list[str]]:
    """Return (dead pins, unreviewed candidates) across both curated tables."""
    # one (sseqid, name, blurb) triple list per source, loaded on first use
    sources: dict[str, list[tuple[str, str, str]]] = {}
    dead: list[str] = []
    unreviewed: list[str] = []

    for filename in CURATED_TABLES:
        rows = _curated_rows(filename)
        # a display name may be split across rows when the accessions behind it are
        # not one feature (snapgene's oriV is both the incP and the F origin), so a
        # name is only unreviewed if no row claims the accession
        pins_by_name: dict[tuple[str, str], set[str]] = {}
        for _, row in rows.iterrows():
            pins_by_name.setdefault((row["db"].strip(), row["name"]), set()).update(
                _pinned_accessions(row)
            )

        reported: set[tuple[str, str]] = set()
        for _, row in rows.iterrows():
            database = row["db"].strip()
            if database not in sources:
                frame = _load_source(database)
                sources[database] = list(
                    zip(
                        frame["sseqid"].astype(str),
                        frame["name"].astype(str),
                        frame["blurb"].astype(str),
                        strict=True,
                    )
                )
            records = sources[database]
            known = {accession for accession, _, _ in records}

            for accession in sorted(_pinned_accessions(row) - known):
                dead.append(
                    f"{filename}: {row['name']} pins {database}:{accession}, "
                    "which is not in the installed bundle"
                )

            pinned = pins_by_name[(database, row["name"])]
            same_name = {
                accession: blurb
                for accession, name, blurb in records
                if name == row["name"] and accession not in pinned
            }
            for accession, blurb in sorted(same_name.items()):
                if (database, accession) in EXCLUSIONS:
                    continue
                if (database, accession) in reported:
                    continue
                reported.add((database, accession))
                unreviewed.append(
                    f"{filename}: {row['name']} does not pin {database}:{accession} "
                    f"-- {blurb.split(' - ', 1)[-1][:90]}"
                )

    # Global and fragment-region suppressions are exact records, not claims shared by
    # display name, so only dead pins matter; a same-name record must never inherit one.
    for filename in PIN_ONLY_TABLES:
        for _, row in _curated_rows(filename).iterrows():
            database = row["db"].strip()
            if database not in SOURCE_DATABASES:
                dead.append(f"{filename}: unknown source {database!r}")
                continue
            if database not in sources:
                frame = _load_source(database)
                sources[database] = list(
                    zip(
                        frame["sseqid"].astype(str),
                        frame["name"].astype(str),
                        frame["blurb"].astype(str),
                        strict=True,
                    )
                )
            known = {accession for accession, _, _ in sources[database]}
            for accession in sorted(_pinned_accessions(row) - known):
                dead.append(
                    f"{filename}: pins {database}:{accession}, which is not in "
                    "the installed bundle"
                )

    # Pair overrides pin both sides. Rfam has no SQLite descriptions database, so
    # its accessions are validated by the exhaustive nested audit instead.
    overrides = _curated_rows(NESTED_OVERRIDE_TABLE)
    for _, row in overrides.iterrows():
        for role, db_column, accession_column in (
            ("parent", "parent_db", "parent_sseqid"),
            ("child", "child_db", "child_sseqid"),
        ):
            database = row[db_column].strip()
            if database == "Rfam":
                continue
            if database not in SOURCE_DATABASES:
                dead.append(
                    f"{NESTED_OVERRIDE_TABLE}: {role} uses unknown source {database!r}"
                )
                continue
            if database not in sources:
                frame = _load_source(database)
                sources[database] = list(
                    zip(
                        frame["sseqid"].astype(str),
                        frame["name"].astype(str),
                        frame["blurb"].astype(str),
                        strict=True,
                    )
                )
            known = {accession for accession, _, _ in sources[database]}
            accessions = {
                part.strip()
                for part in row[accession_column].split(";")
                if part.strip()
            }
            for accession in sorted(accessions - known):
                dead.append(
                    f"{NESTED_OVERRIDE_TABLE}: {role} pins {database}:{accession}, "
                    "which is not in the installed bundle"
                )

    # Composite regions pin both the larger record and the component that explains
    # part of it. Validate source-bundle drift on both sides of the relationship.
    composite_regions = _curated_rows(COMPOSITE_REGION_TABLE)
    for _, row in composite_regions.iterrows():
        pins_valid = True
        for role, db_column, accession_column in (
            ("record", "db", "sseqid"),
            ("component", "component_db", "component_sseqid"),
        ):
            database = row[db_column].strip()
            accession = row[accession_column].strip()
            if database not in SOURCE_DATABASES:
                dead.append(
                    f"{COMPOSITE_REGION_TABLE}: {role} uses unknown source {database!r}"
                )
                pins_valid = False
                continue
            if database not in sources:
                frame = _load_source(database)
                sources[database] = list(
                    zip(
                        frame["sseqid"].astype(str),
                        frame["name"].astype(str),
                        frame["blurb"].astype(str),
                        strict=True,
                    )
                )
            known = {candidate for candidate, _, _ in sources[database]}
            if accession not in known:
                dead.append(
                    f"{COMPOSITE_REGION_TABLE}: {role} pins "
                    f"{database}:{accession}, which is not in the installed bundle"
                )
                pins_valid = False
        try:
            start = int(row["region_start"])
            end = int(row["region_end"])
        except ValueError:
            dead.append(
                f"{COMPOSITE_REGION_TABLE}: non-integer interval for "
                f"{row['db']}:{row['sseqid']}"
            )
            continue
        if start < 1 or end < start:
            dead.append(
                f"{COMPOSITE_REGION_TABLE}: invalid interval {start}-{end} for "
                f"{row['db']}:{row['sseqid']}"
            )
            continue
        # SnapGene records and components are nucleotide sequences, so validate the
        # strongest form of the claim: the curated interval is the exact component.
        if pins_valid and row["db"] == row["component_db"] == "snapgene":
            try:
                record_sequence = _snapgene_sequence(row["sseqid"])
                component_sequence = _snapgene_sequence(row["component_sseqid"])
            except (OSError, subprocess.CalledProcessError) as error:
                detail = getattr(error, "stderr", "") or str(error)
                dead.append(
                    f"{COMPOSITE_REGION_TABLE}: could not extract SnapGene "
                    f"sequences for validation: {detail.strip()}"
                )
                continue
            if end > len(record_sequence):
                dead.append(
                    f"{COMPOSITE_REGION_TABLE}: interval {start}-{end} exceeds "
                    f"{row['db']}:{row['sseqid']} length {len(record_sequence)}"
                )
            elif record_sequence[start - 1 : end] != component_sequence:
                dead.append(
                    f"{COMPOSITE_REGION_TABLE}: interval {start}-{end} of "
                    f"{row['db']}:{row['sseqid']} no longer equals component "
                    f"{row['component_db']}:{row['component_sseqid']}"
                )
        elif pins_valid:
            dead.append(
                f"{COMPOSITE_REGION_TABLE}: cannot fully validate {row['db']}:"
                f"{row['sseqid']} -> {row['component_db']}:"
                f"{row['component_sseqid']}; extend sequence validation before "
                "curating non-SnapGene composite regions"
            )
    return dead, unreviewed


def report(_: argparse.Namespace) -> int:
    dead, unreviewed = audit()
    print(f"dead pins: {len(dead)}")
    for line in dead:
        print(f"  {line}")
    print(f"unreviewed same-name records: {len(unreviewed)}")
    for line in unreviewed:
        print(f"  {line}")
    print(f"documented exclusions: {len(EXCLUSIONS)}")
    return 0


def check(_: argparse.Namespace) -> int:
    dead, unreviewed = audit()
    if not dead and not unreviewed:
        print("curated pins are consistent with the installed database bundle")
        return 0
    for line in (*dead, *unreviewed):
        print(line, file=sys.stderr)
    print(
        f"\n{len(dead)} dead pin(s), {len(unreviewed)} unreviewed record(s). "
        "Adjudicate each: pin it in the CSV, or add it to EXCLUSIONS with a reason.",
        file=sys.stderr,
    )
    return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    subparsers = parser.add_subparsers(dest="command", required=True)

    report_parser = subparsers.add_parser("report", help="print a drift report")
    report_parser.set_defaults(function=report)

    check_parser = subparsers.add_parser("check", help="exit nonzero on any drift")
    check_parser.set_defaults(function=check)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    raise SystemExit(args.function(args))


if __name__ == "__main__":
    main()
