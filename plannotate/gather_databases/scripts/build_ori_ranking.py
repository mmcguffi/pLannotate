#!/usr/bin/env python3
"""
Count origin-of-replication frequencies across a GenBank plasmid corpus.

This is a reference/curation aid, not part of the runtime path. The shipped
ranking (``plannotate/data/data/ori_ranking.csv``) is hand-curated by biological
prevalence rather than raw corpus counts, because public dumps are skewed by
individual lab collections. This script reports the corpus evidence behind that
curation and normalizes labels to the canonical snapgene names so the output
joins cleanly to the runtime ranking.

The corpus used originally was the Addgene plasmid dump that lives in this repo's
git history at ``25685f9:data/addgene_gbks/`` (4,989 plasmids). Recover it with:

    git archive 25685f9 data/addgene_gbks | tar -x -C <corpus_dir>

Usage:
    python3 build_ori_ranking.py --corpus <dir_of_gbks> \\
        --descriptions <path/to/descriptions.db> \\
        --output ori_corpus_counts.csv
"""

import argparse
import sqlite3
import unicodedata
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# ASCII spellings seen in the wild mapped to the canonical snapgene names.
_NAME_ALIASES = {
    "2u ori": "2μ ori",
    "r6k gamma ori": "R6K γ ori",
}


def _canonical_label(label: str) -> str:
    """Normalize a GenBank label toward its canonical snapgene name."""
    normalized = unicodedata.normalize("NFKC", label).strip()
    return _NAME_ALIASES.get(normalized.lower(), normalized)


def _snapgene_ori_names(descriptions_db: Path) -> set[str]:
    """Return the set of snapgene names classified as rep_origin."""
    connection = sqlite3.connect(descriptions_db)
    try:
        rows = connection.execute(
            "SELECT DISTINCT name FROM snapgene WHERE type = 'rep_origin'"
        ).fetchall()
    finally:
        connection.close()
    return {name for (name,) in rows if name}


def count_origins(corpus_dir: Path, valid_names: set[str]) -> Counter:
    """Count rep_origin labels (restricted to known snapgene names) in a corpus."""
    counts: Counter = Counter()
    for path in sorted(corpus_dir.rglob("*.gb*")):
        try:
            records = SeqIO.parse(path, "genbank")
            for record in records:
                for feature in record.features:
                    if feature.type != "rep_origin":
                        continue
                    labels = feature.qualifiers.get("label", [])
                    if not labels:
                        continue
                    name = _canonical_label(labels[0])
                    if name in valid_names:
                        counts[name] += 1
        except (ValueError, OSError) as exc:
            print(f"Skipping {path}: {exc}")
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Count rep_origin frequencies across a GenBank plasmid corpus"
    )
    parser.add_argument(
        "--corpus", required=True, help="Directory of GenBank plasmid files"
    )
    parser.add_argument(
        "--descriptions",
        required=True,
        help="Path to the snapgene descriptions.db (for canonical names)",
    )
    parser.add_argument(
        "--output",
        default="ori_corpus_counts.csv",
        help="Output CSV (name,count,rank). DEFAULT: ori_corpus_counts.csv",
    )
    args = parser.parse_args()

    corpus_dir = Path(args.corpus)
    if not corpus_dir.is_dir():
        parser.error(f"corpus directory does not exist: {corpus_dir}")

    valid_names = _snapgene_ori_names(Path(args.descriptions))
    counts = count_origins(corpus_dir, valid_names)

    ranked = counts.most_common()
    frame = pd.DataFrame(
        [
            {"name": name, "count": count, "rank": index + 1}
            for index, (name, count) in enumerate(ranked)
        ]
    )
    frame.to_csv(args.output, index=False)
    print(f"Wrote {len(frame)} origins to {args.output}")
    for name, count in ranked:
        print(f"  {count:>6}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
