#!/usr/bin/env python3
"""Sequence and file validation utilities.

Can be used as a module or as a command-line sequence-file validator.
"""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Constants for validation
VALID_GENBANK_EXTS = [".gbk", ".gb", ".gbf", ".gbff"]
VALID_FASTA_EXTS = [".fa", ".fasta", ".fas", ".fna"]
MAX_PLAS_SIZE = 50000
IUPAC_NUCLEOTIDES = "GATCRYWSMKHBVDNgatcrywsmkhbvdn"


class InvalidSequenceError(ValueError):
    """Raised when sequence validation fails."""


def get_name_ext(file_loc: str | Path) -> tuple[str, str]:
    """Extract name and extension from file path."""
    path = Path(file_loc)
    return path.stem, path.suffix.lower()


def validate_sequence(seq: str, max_length: int | None = MAX_PLAS_SIZE) -> None:
    """Validate DNA sequence content and length."""
    if max_length is not None and max_length < 1:
        raise ValueError("max_length must be at least 1")
    if not seq:
        raise InvalidSequenceError("Sequence is empty")
    if not set(seq).issubset(IUPAC_NUCLEOTIDES):
        error = (
            "Sequence contains invalid characters -- must be ATCG "
            "and/or valid IUPAC nucleotide ambiguity code"
        )
        raise InvalidSequenceError(error)

    if max_length is not None and len(seq) > max_length:
        error = (
            f"Are you sure this is an engineered plasmid? Entry size is too large "
            f"-- must be {max_length} bases or less."
        )
        raise InvalidSequenceError(error)


def validate_file(
    file: str | Path,
    ext: str | None = None,
    max_length: int | None = MAX_PLAS_SIZE,
) -> SeqRecord:
    """
    Validate a sequence file and return the entry.

    Can raise InvalidSequenceError if not valid.
    """
    path = Path(file)
    ext = (ext or path.suffix).lower()
    if ext in VALID_FASTA_EXTS:
        record = _validate_fasta_file(path)
    elif ext in VALID_GENBANK_EXTS:
        record = _validate_genbank_file(path)
    else:
        raise ValueError("must be a FASTA or GenBank file")

    if len(record) != 1:
        error = (
            "File contains multiple entries -- please submit a single sequence file."
        )
        raise InvalidSequenceError(error)

    validate_sequence(str(record[0].seq), max_length)
    return record[0]


def _validate_fasta_file(file: Path) -> list[SeqRecord]:
    """Validate FASTA file format and content."""
    with file.open() as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise InvalidSequenceError(
            "Malformed FASTA file; submit a file in standard FASTA format"
        )
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


def _validate_genbank_file(file: Path) -> list[SeqRecord]:
    """Validate GenBank format while preserving record annotations."""
    with file.open() as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise InvalidSequenceError(
            "Malformed GenBank file; submit a file in standard GenBank format"
        )
    return records


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


def main() -> None:
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
    except Exception as exc:
        print(f"Error validating {args.input}: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
