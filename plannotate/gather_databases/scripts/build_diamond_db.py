#!/usr/bin/env python3
"""
Build DIAMOND database from TSV or FASTA file.

Usage:
    python3 build_diamond_db.py --input fpbase.tsv --output fpbase.dmnd --format tsv
    python3 build_diamond_db.py --input proteins.fasta --output proteins.dmnd --format fasta
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def tsv_to_fasta(tsv_file: Path, fasta_file: Path) -> bool:
    """Convert TSV file to FASTA format for DIAMOND database creation."""
    try:
        df = pd.read_csv(tsv_file, sep="\t")

        # Expected columns: name, slug, seq, blurb
        if "seq" not in df.columns or "name" not in df.columns:
            print("Error: TSV file must contain 'seq' and 'name' columns")
            return False

        with open(fasta_file, "w") as f:
            for _, row in df.iterrows():
                name = row["name"]
                seq = row["seq"]
                if pd.notna(seq) and pd.notna(name):
                    f.write(f">{name}\n{seq}\n")

        print(f"Converted {len(df)} proteins to FASTA format")
        return True
    except Exception as e:
        print(f"Error converting TSV to FASTA: {e}")
        return False


def build_diamond_database(fasta_file: Path, output_file: Path) -> bool:
    """Build DIAMOND database from FASTA file."""
    try:
        # Remove .dmnd extension from output for diamond makedb
        db_name = str(output_file).replace(".dmnd", "")

        cmd = ["diamond", "makedb", "--in", str(fasta_file), "--db", db_name]

        print(f"Building DIAMOND database: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print("DIAMOND database created successfully")
            return True
        else:
            print(f"Error creating DIAMOND database: {result.stderr}")
            return False
    except FileNotFoundError:
        print("Error: diamond command not found. Please ensure DIAMOND is installed.")
        return False
    except Exception as e:
        print(f"Error creating DIAMOND database: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Build DIAMOND database from TSV or FASTA file"
    )
    parser.add_argument("--input", required=True, help="Input file (TSV or FASTA)")
    parser.add_argument(
        "--output", required=True, help="Output DIAMOND database file (.dmnd)"
    )
    parser.add_argument(
        "--format", choices=["tsv", "fasta"], required=True, help="Input file format"
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    if not input_file.exists():
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if args.format == "tsv":
        # Convert TSV to FASTA first
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_fasta:
            temp_fasta_path = Path(tmp_fasta.name)

        if not tsv_to_fasta(input_file, temp_fasta_path):
            sys.exit(1)

        # Build database from temporary FASTA
        success = build_diamond_database(temp_fasta_path, output_file)

        # Clean up temporary file
        temp_fasta_path.unlink()

    elif args.format == "fasta":
        # Build database directly from FASTA
        success = build_diamond_database(input_file, output_file)

    if not success:
        sys.exit(1)

    print(f"DIAMOND database created: {output_file}")


if __name__ == "__main__":
    main()
