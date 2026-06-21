#!/usr/bin/env python3
"""
Convert CSV/TSV files to SQLite database for database descriptions.

Usage:
    python3 csv_to_sqlite.py --input file.csv --output descriptions.db --table table_name
    python3 csv_to_sqlite.py --input file.tsv --output descriptions.db --table table_name --delimiter tab
"""

import argparse
import sqlite3
from pathlib import Path
from typing import Mapping

import pandas as pd


def create_sqlite_from_csv(
    input_file: Path,
    output_file: Path,
    table_name: str,
    delimiter: str = ",",
    append: bool = False,
    no_header: bool = False,
    add_type: str | None = None,
    column_renames: Mapping[str, str] | None = None,
) -> bool:
    """Convert CSV/TSV file to SQLite database."""
    try:
        # Read the CSV/TSV file
        separator = "\t" if delimiter == "tab" else delimiter
        df = pd.read_csv(
            input_file, sep=separator, header=None if no_header else "infer"
        )

        # Enforce standard 4-column format: ["sseqid", "name", "type", "blurb"]
        if no_header:
            if len(df.columns) == 3:
                # 3-column format - will need type added via --add-type
                df.columns = ["sseqid", "name", "blurb"]
            elif len(df.columns) == 4:
                # Standard 4-column format
                df.columns = ["sseqid", "name", "type", "blurb"]
            else:
                raise ValueError(
                    f"Expected 3 or 4 columns, got {len(df.columns)}. Only standard format supported."
                )
        elif column_renames:
            df = df.rename(columns=column_renames)

        # Add type column if specified (for 3-column files)
        if add_type:
            if "type" in df.columns:
                print("Warning: Type column already exists, --add-type ignored")
            else:
                df.insert(2, "type", add_type)  # Insert type as 3rd column

        # Validate final format - must have exactly 4 columns in correct order
        expected_columns = ["sseqid", "name", "type", "blurb"]
        if list(df.columns) != expected_columns:
            raise ValueError(
                f"Final DataFrame must have columns {expected_columns}, got {list(df.columns)}"
            )

        print(f"Read {len(df)} rows from {input_file}")

        # Create SQLite database
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Remove existing database if not appending
        if not append and output_file.exists():
            output_file.unlink()

        conn = sqlite3.connect(output_file)

        # Write dataframe to SQLite
        if append:
            df.to_sql(table_name, conn, index=False, if_exists="append")
        else:
            df.to_sql(table_name, conn, index=False, if_exists="replace")

        # Create index on first column (usually the ID column)
        if len(df.columns) > 0:
            first_col = df.columns[0]
            conn.execute(
                f"CREATE INDEX IF NOT EXISTS idx_{table_name}_{first_col} ON {table_name} ({first_col})"
            )

        conn.close()

        action = "Appended to" if append else "Created"
        print(f"{action} SQLite database: {output_file}")
        print(f"  Table: {table_name}")
        print(f"  Columns: {list(df.columns)}")

        return True

    except Exception as e:
        print(f"Error converting {input_file} to SQLite: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Convert CSV/TSV to SQLite database")
    parser.add_argument("--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("--output", required=True, help="Output SQLite database file")
    parser.add_argument("--table", required=True, help="Table name in SQLite database")
    parser.add_argument(
        "--delimiter",
        default=",",
        help="CSV delimiter (default: comma, use 'tab' for TSV)",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing database instead of replacing",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="File has no header row, assign standard column names",
    )
    parser.add_argument(
        "--add-type",
        type=str,
        help="Add a type column with the specified value (e.g., 'CDS', 'ncRNA')",
    )
    parser.add_argument(
        "--rename",
        action="append",
        default=[],
        metavar="SOURCE=TARGET",
        help="Rename an input column; may be supplied more than once",
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    if not input_file.exists():
        print(f"Error: Input file does not exist: {input_file}")
        return 1

    try:
        column_renames = dict(rename.split("=", 1) for rename in args.rename)
    except ValueError:
        parser.error("--rename values must use SOURCE=TARGET syntax")

    success = create_sqlite_from_csv(
        input_file,
        output_file,
        args.table,
        args.delimiter,
        args.append,
        args.no_header,
        args.add_type,
        column_renames,
    )
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
