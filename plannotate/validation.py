#!/usr/bin/env python3
"""Sequence and file validation utilities.

Can be used both as a module for importing validation functions
or as a command-line script for validating sequence files (snakemake-able).
"""

import argparse
from io import StringIO
import os
import sys
from typing import Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Constants for validation
VALID_GENBANK_EXTS = [".gbk", ".gb", ".gbf", ".gbff"]
VALID_FASTA_EXTS = [".fa", ".fasta", ".fas", ".fna"]
MAX_PLAS_SIZE = 50000
IUPAC_NUCLEOTIDES = "GATCRYWSMKHBVDNgatcrywsmkhbvdn"


class InvalidSequenceError(ValueError):
    """Raised when sequence validation fails."""

    pass


def get_name_ext(file_loc: str) -> Tuple[str, str]:
    """Extract name and extension from file path."""
    base = os.path.basename(file_loc)
    name = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1]
    return name, ext


def validate_sequence(seq: str, max_length: int = MAX_PLAS_SIZE) -> None:
    """Validate DNA sequence content and length."""
    if not set(seq).issubset(IUPAC_NUCLEOTIDES):
        error = (
            "Sequence contains invalid characters -- must be ATCG "
            "and/or valid IUPAC nucleotide ambiguity code"
        )
        raise ValueError(error)

    if len(seq) > max_length:
        error = (
            f"Are you sure this is an engineered plasmid? Entry size is too large "
            f"-- must be {max_length} bases or less."
        )
        raise ValueError(error)


def validate_file(
    file: str,
    ext: str,
    max_length: int = MAX_PLAS_SIZE,
) -> SeqRecord:
    """
    Validate a sequence file and return the entry.

    Can raise InvalidSequenceError if not valid.
    """
    if ext in VALID_FASTA_EXTS:
        record = _validate_fasta_file(file)
    elif ext in VALID_GENBANK_EXTS:
        record = _validate_genbank_file(file)
    else:
        ERROR = "must be a FASTA or GenBank file"
        raise ValueError(ERROR)

    if len(record) != 1:
        ERROR = (
            "File contains multiple entries -- please submit a single sequence file."
        )
        raise InvalidSequenceError(ERROR)

    validate_sequence(str(record[0].seq), max_length)
    return record[0]


def _validate_fasta_file(file: str) -> list[SeqRecord]:
    """Validate FASTA file format and content."""
    try:
        with open(file) as handle:
            record = list(SeqIO.parse(handle, "fasta"))
        if not record:
            raise InvalidSequenceError(
                "Malformed fasta file --> please submit a fasta file in standard format"
            )

        # Ensure DNA molecule type annotation
        record[0].annotations["molecule_type"] = "DNA"

        # Round-trip through FASTA to validate format and normalize the records.
        with StringIO() as buffer:
            SeqIO.write(record, buffer, "fasta")
            buffer.seek(0)
            record = list(SeqIO.parse(buffer, "fasta"))

        return record
    except IndexError:
        raise InvalidSequenceError(
            "Malformed fasta file --> please submit a fasta file in standard format"
        )


def _validate_genbank_file(file: str) -> list[SeqRecord]:
    """Validate GenBank file format and extract sequence."""
    try:
        with open(file) as handle:
            record = list(SeqIO.parse(handle, "gb"))[0]

        # Convert to FASTA format for consistency
        with StringIO() as buffer:
            SeqIO.write(record, buffer, "fasta")
            buffer.seek(0)
            record = list(SeqIO.parse(buffer, "fasta"))

        return record
    except IndexError:
        raise InvalidSequenceError(
            "Malformed Genbank file --> please submit a Genbank file in standard format"
        )


def validate_and_write_fasta(
    input_file: str,
    output_file: str,
    max_length: int = MAX_PLAS_SIZE,
) -> SeqRecord:
    """
    Complete validation pipeline: validate input file and write cleaned FASTA output.

    This is the main function used by the validation script.
    """
    _, ext = get_name_ext(input_file)
    record = validate_file(input_file, ext, max_length)
    SeqIO.write(record, output_file, "fasta")
    return record


def main():
    """Command-line interface for validating input sequence files."""
    parser = argparse.ArgumentParser(description="Validate input sequence files")
    parser.add_argument("--input", required=True, help="Input sequence file path")
    parser.add_argument(
        "--output", required=True, help="Output validated FASTA file path"
    )
    parser.add_argument(
        "--max-length", type=int, default=50000, help="Maximum sequence length"
    )

    args = parser.parse_args()

    try:
        validate_and_write_fasta(args.input, args.output, args.max_length)
        name, ext = get_name_ext(args.input)
        print(f"Validated: {name}{ext}")
    except Exception as e:
        print(f"Error validating {args.input}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
