#!/usr/bin/env python3
"""Find database features that pLannotate annotates inside other features.

The SnapGene database stores DNA, so its records can be annotated directly. FPbase
stores proteins; those records are deterministically back-translated so DIAMOND can
search their amino-acid content. Nucleotide and RNA hits in the synthetic FPbase DNA
are discarded because they have no biological meaning.

The parent feature's exact self-hit and other hits covering its complete sequence are
not nested features. Every strict subinterval retained by detailed mode is reported,
including fragments: short false positives are useful database-curation candidates.

    python tools/nested_feature_audit.py \
        --csv docs/nested-feature-audit.csv \
        --markdown docs/nested-feature-audit.md
"""

import argparse
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, cast

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plannotate import _package_data  # noqa: E402
from plannotate.annotate import annotate_batch  # noqa: E402

DATA = ROOT / "plannotate" / "data"
SNAPGENE_DB = DATA / "BLAST_dbs" / "snapgene"
FPBASE_DB = DATA / "diamond_dbs" / "fpbase.dmnd"
MANIFEST = DATA / "database-manifest.json"

REPORT_COLUMNS = [
    "parent_db",
    "parent_sseqid",
    "parent_name",
    "parent_type",
    "parent_length",
    "nested_db",
    "nested_sseqid",
    "nested_name",
    "nested_type",
    "nested_strand",
    "nested_start",
    "nested_end",
    "nested_length",
    "nested_subject_length",
    "nested_subject_start",
    "nested_subject_end",
    "percent_identity",
    "percent_match",
    "evalue",
    "fragment",
]


def _run(command: list[str]) -> str:
    """Run a sequence-database utility and return its standard output."""
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"Required executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip()
        raise RuntimeError(f"{' '.join(command)} failed: {detail}") from error
    return result.stdout


def _metadata(database: str, path: Path) -> dict[str, dict[str, str]]:
    """Load the stable id, display name, and type for one feature database."""
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing; run 'plannotate setupdb'")
    with sqlite3.connect(path) as connection:
        rows = connection.execute(
            f"SELECT sseqid, name, type FROM {database}"
        ).fetchall()
    return {
        str(sequence_id): {
            "name": str(name or sequence_id),
            "type": str(feature_type or "misc_feature"),
        }
        for sequence_id, name, feature_type in rows
    }


def load_snapgene_features() -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Extract every DNA record and its metadata from the installed bundle."""
    output = _run(
        [
            "blastdbcmd",
            "-db",
            str(SNAPGENE_DB),
            "-entry",
            "all",
            "-outfmt",
            "%a\t%s",
        ]
    )
    sequences: dict[str, str] = {}
    for line in output.splitlines():
        sequence_id, sequence = line.split("\t", 1)
        sequences[sequence_id] = sequence
    metadata = _align_metadata_ids(
        "snapgene",
        sequences,
        _metadata("snapgene", SNAPGENE_DB.with_suffix(".db")),
    )
    return sequences, metadata


def _back_translate(protein: str) -> str:
    """Create DNA whose six-frame translation contains the supplied protein."""
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    codons = table.back_table
    return "".join(
        "TAA" if residue == "*" else codons.get(residue.upper(), "NNN")
        for residue in protein
    )


def load_fpbase_features() -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Extract every protein and back-translate it for DIAMOND blastx."""
    if not FPBASE_DB.is_file():
        raise FileNotFoundError(f"{FPBASE_DB} is missing; run 'plannotate setupdb'")
    with tempfile.TemporaryDirectory(prefix="plannotate-nested-audit-") as directory:
        fasta = Path(directory) / "fpbase.faa"
        _run(["diamond", "getseq", "--db", str(FPBASE_DB), "--out", str(fasta)])
        sequences = {
            record.id: _back_translate(str(record.seq))
            for record in SeqIO.parse(fasta, "fasta")
        }
    metadata = _align_metadata_ids(
        "fpbase",
        sequences,
        _metadata("fpbase", FPBASE_DB.with_suffix(".db")),
    )
    return sequences, metadata


def _align_metadata_ids(
    database: str,
    sequences: Mapping[str, str],
    metadata: Mapping[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Key metadata by the sequence index ids, accepting unique case-only drift."""
    folded_metadata: dict[str, list[str]] = {}
    for metadata_id in metadata:
        folded_metadata.setdefault(metadata_id.casefold(), []).append(metadata_id)

    aligned: dict[str, dict[str, str]] = {}
    unmatched_sequences: list[str] = []
    consumed_metadata: set[str] = set()
    for sequence_id in sequences:
        lookup_id: str | None = sequence_id if sequence_id in metadata else None
        if lookup_id is None:
            candidates = folded_metadata.get(sequence_id.casefold(), [])
            if len(candidates) == 1:
                lookup_id = candidates[0]
        if lookup_id is None:
            unmatched_sequences.append(sequence_id)
            continue
        aligned[sequence_id] = metadata[lookup_id]
        consumed_metadata.add(lookup_id)

    unmatched_metadata = sorted(set(metadata) - consumed_metadata)
    if unmatched_sequences or unmatched_metadata:
        raise RuntimeError(
            f"{database} sequence/metadata ids differ: "
            f"{len(unmatched_metadata)} sequence(s) missing and "
            f"{len(unmatched_sequences)} metadata row(s) missing"
        )
    return aligned


def nested_rows(
    results: Mapping[str, pd.DataFrame],
    parents: Mapping[str, dict[str, Any]],
    source_methods: Mapping[str, str],
) -> tuple[pd.DataFrame, list[str]]:
    """Convert detailed annotations to strict parent/child containment rows."""
    rows: list[dict[str, Any]] = []
    missing_self_hits: list[str] = []
    for key, annotations in results.items():
        parent = parents[key]
        source_database = str(parent["db"])
        source_id = str(parent["sseqid"])
        parent_length = int(parent["length"])
        self_hit = (annotations["db"] == source_database) & (
            annotations["sseqid"] == source_id
        )
        if not self_hit.any():
            missing_self_hits.append(key)

        candidates = annotations.loc[~self_hit].copy()
        # A complete synonym/homolog is not physically nested in its parent.
        candidates = candidates.loc[
            (candidates["qstart"].astype(int) > 0)
            | (candidates["qend"].astype(int) < parent_length)
        ]
        if source_database == "fpbase":
            # The query DNA is synthetic. Only translated-protein searches reflect
            # information that exists in the FPbase source record.
            candidates = candidates.loc[
                candidates["db"].map(source_methods).eq("diamond")
            ]

        for _, child in candidates.iterrows():
            rows.append(
                {
                    "parent_db": source_database,
                    "parent_sseqid": source_id,
                    "parent_name": parent["name"],
                    "parent_type": parent["type"],
                    "parent_length": parent_length,
                    "nested_db": child["db"],
                    "nested_sseqid": child["sseqid"],
                    "nested_name": child["name"],
                    "nested_type": child["type"],
                    "nested_strand": int(child["sframe"]),
                    "nested_start": int(child["qstart"]),
                    "nested_end": int(child["qend"]),
                    "nested_length": int(child["length"]),
                    "nested_subject_length": int(child["slen"]),
                    "nested_subject_start": int(child["sstart"]),
                    "nested_subject_end": int(child["send"]),
                    "percent_identity": float(child["pident"]),
                    "percent_match": float(child["abs percmatch"]),
                    "evalue": float(child["evalue"]),
                    "fragment": bool(child["fragment"]),
                }
            )

    report = pd.DataFrame(rows, columns=REPORT_COLUMNS)
    if not report.empty:
        report = report.sort_values(
            [
                "parent_db",
                "parent_name",
                "parent_sseqid",
                "nested_start",
                "nested_end",
                "nested_name",
            ],
            kind="stable",
        ).reset_index(drop=True)
    return report, missing_self_hits


def _manifest_summary() -> str:
    if not MANIFEST.is_file():
        return "unknown database bundle"
    manifest = json.loads(MANIFEST.read_text())
    return f"{manifest.get('bundle', 'database bundle')} ({manifest.get('build_date', 'unknown date')})"


def write_markdown(
    path: Path,
    report: pd.DataFrame,
    source_counts: Mapping[str, int],
    missing_self_hits: Iterable[str],
) -> None:
    """Write a readable companion to the complete pairwise CSV report."""
    missing = list(missing_self_hits)
    problematic_counts = (
        report.groupby("parent_db")["parent_sseqid"].nunique().to_dict()
        if not report.empty
        else {}
    )
    lines = [
        "# Nested feature audit",
        "",
        f"Database bundle: `{_manifest_summary()}`.",
        "",
        (
            "Each installed SnapGene DNA feature and FPbase protein feature was "
            "annotated as a linear sequence in detailed mode. Exact self-hits and "
            "other full-length hits are omitted. FPbase records are protein-only, "
            "so only translated-protein results are meaningful and retained."
        ),
        "",
        (
            "These are review candidates, not automatic deletion decisions: some "
            "source records intentionally contain composite biological features."
        ),
        "",
        (
            "This Markdown file and its companion CSV are generated snapshots, not "
            "runtime policy inputs. Curated runtime decisions live in "
            "`plannotate/data/data/nested_feature_overrides.csv` and "
            "`feature_suppressions.csv`; see "
            "[`nested-feature-curation-policy.md`](nested-feature-curation-policy.md) "
            "for their semantics and the complete maintainer workflow."
        ),
        "",
        (
            "Regenerate with `python tools/nested_feature_audit.py --csv "
            "docs/nested-feature-audit.csv --markdown "
            "docs/nested-feature-audit.md`."
        ),
        "",
        f"Nested annotations: **{len(report)}**",
        "",
    ]
    for database in ("snapgene", "fpbase"):
        lines.append(
            f"- {database}: {problematic_counts.get(database, 0)} of "
            f"{source_counts.get(database, 0)} source features have nested hits"
        )
    lines.extend(
        [
            "",
            f"Missing exact self-hits: **{len(missing)}**",
            "",
            "The complete machine-readable result is in "
            "[`nested-feature-audit.csv`](nested-feature-audit.csv).",
            "",
            "## Problematic features",
            "",
        ]
    )

    if report.empty:
        lines.append("No nested annotations were found.")
    else:
        group_columns = [
            "parent_db",
            "parent_sseqid",
            "parent_name",
            "parent_type",
            "parent_length",
        ]
        for parent, children in report.groupby(group_columns, sort=False):
            parent_values = cast(tuple[object, object, object, object, object], parent)
            database, sequence_id, name, feature_type = map(str, parent_values[:4])
            length = int(str(parent_values[4]))
            lines.append(
                f"- `{database}:{sequence_id}` — {name} ({feature_type}, {length} bp)"
            )
            for _, child in children.iterrows():
                fragment = ", fragment" if child["fragment"] else ""
                lines.append(
                    "  - "
                    f"`{child['nested_db']}:{child['nested_sseqid']}` — "
                    f"{child['nested_name']} ({child['nested_type']}, "
                    f"{child['nested_start']}–{child['nested_end']}, "
                    f"{child['percent_identity']:.1f}% identity{fragment})"
                )

    if missing:
        lines.extend(["", "## Records missing an exact self-hit", ""])
        lines.extend(f"- `{key}`" for key in sorted(missing))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def audit(cores: int) -> tuple[pd.DataFrame, dict[str, int], list[str]]:
    """Load every source feature, annotate it, and return the nested-hit report."""
    snapgene_sequences, snapgene_metadata = load_snapgene_features()
    fpbase_sequences, fpbase_metadata = load_fpbase_features()

    sequences: dict[str, str] = {}
    parents: dict[str, dict[str, Any]] = {}
    for database, source_sequences, source_metadata in (
        ("snapgene", snapgene_sequences, snapgene_metadata),
        ("fpbase", fpbase_sequences, fpbase_metadata),
    ):
        for sequence_id, sequence in source_sequences.items():
            key = f"{database}:{sequence_id}"
            sequences[key] = sequence
            parents[key] = {
                "db": database,
                "sseqid": sequence_id,
                "name": source_metadata[sequence_id]["name"],
                "type": source_metadata[sequence_id]["type"],
                "length": len(sequence),
            }

    yaml_path = _package_data.get_yaml_path()
    source_methods = {
        name: str(config["method"])
        for name, config in _package_data.get_yaml(yaml_path).items()
    }
    results = annotate_batch(
        sequences,
        yaml_file=yaml_path,
        linear=True,
        is_detailed=True,
        cores=cores,
        apply_nested_policy=False,
    )
    report, missing_self_hits = nested_rows(results, parents, source_methods)
    counts = {
        "snapgene": len(snapgene_sequences),
        "fpbase": len(fpbase_sequences),
    }
    return report, counts, missing_self_hits


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "artifacts" / "nested-feature-audit.csv",
        help="complete parent/child report",
    )
    parser.add_argument(
        "--markdown",
        type=Path,
        help="optional human-readable grouped report",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="total search cores (default: all available)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.cores < 1:
        raise SystemExit("--cores must be at least 1")
    report, source_counts, missing_self_hits = audit(args.cores)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.csv, index=False)
    if args.markdown:
        write_markdown(
            args.markdown,
            report,
            source_counts,
            missing_self_hits,
        )
    problematic = report[["parent_db", "parent_sseqid"]].drop_duplicates()
    print(
        f"Annotated {sum(source_counts.values())} source features; found "
        f"{len(report)} nested hits in {len(problematic)} parent features."
    )
    print(f"CSV: {args.csv}")
    if args.markdown:
        print(f"Markdown: {args.markdown}")
    if missing_self_hits:
        print(
            f"Warning: {len(missing_self_hits)} source features had no exact self-hit.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
