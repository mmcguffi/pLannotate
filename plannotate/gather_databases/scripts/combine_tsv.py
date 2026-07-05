#!/usr/bin/env python3
"""Combine headerless TSV chunks without loading the dataset into memory."""

import argparse
import shutil
from pathlib import Path
from typing import Iterable


def combine_files(inputs: Iterable[Path], output: Path) -> None:
    """Concatenate input files in order while preserving them for Snakemake."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as destination:
        for input_path in inputs:
            with input_path.open("rb") as source:
                shutil.copyfileobj(source, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    combine_files(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
