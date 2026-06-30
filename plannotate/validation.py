"""Sequence and file validation utilities."""

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


def validate_records(
    file: str | Path,
    ext: str | None = None,
    max_length: int | None = MAX_PLAS_SIZE,
) -> list[SeqRecord]:
    """Validate every record in a sequence file and return them.

    Unlike :func:`validate_file`, this accepts multi-record FASTA/GenBank files so a
    whole batch can be annotated together. Each record's sequence is validated; an
    empty file raises ``InvalidSequenceError``.
    """
    path = Path(file)
    ext = (ext or path.suffix).lower()
    if ext in VALID_FASTA_EXTS:
        records = _validate_fasta_file(path)
    elif ext in VALID_GENBANK_EXTS:
        records = _validate_genbank_file(path)
    else:
        raise ValueError("must be a FASTA or GenBank file")

    for record in records:
        validate_sequence(str(record.seq), max_length)
    return records


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
