#!/usr/bin/env python3
"""Classify nested-feature audit pairs using the packaged curation policy."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Re-export these names for audit scripts and tests that historically imported this
# command module. The executable policy itself belongs to the package so annotation
# and report generation cannot drift apart.
from plannotate._nested import (  # noqa: E402,F401
    Decision,
    as_bool,
    classify_report,
    classify_row,
    curated_decisions,
    pair_key,
)

DEFAULT_INPUT = ROOT / "docs" / "nested-feature-audit.csv"
DEFAULT_OUTPUT = ROOT / "docs" / "nested-feature-decisions.csv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    classified = classify_report(pd.read_csv(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    classified.to_csv(args.output, index=False)
    counts = classified["rule_status"].value_counts().sort_index().to_dict()
    print(f"Classified {len(classified)} nested calls: {counts}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
