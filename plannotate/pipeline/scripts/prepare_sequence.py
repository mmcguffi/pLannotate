"""
Prepare and validate input sequence for pLannotate pipeline
"""
import json
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def validate_sequence(seq_str, max_length=50000):
    """Validate DNA sequence"""
    # Check sequence length
    if len(seq_str) > max_length:
        raise ValueError(f"Sequence length {len(seq_str)} exceeds maximum {max_length}")
    
    # Check for valid DNA characters
    valid_chars = set('ATCGNatcgn')
    seq_chars = set(seq_str)
    if not seq_chars.issubset(valid_chars):
        invalid = seq_chars - valid_chars
        raise ValueError(f"Invalid DNA characters found: {invalid}")
    
    return seq_str.upper()


def prepare_sequence(input_file, output_fasta, output_metadata, is_linear=False):
    """
    Prepare sequence for annotation pipeline
    - Validate input
    - Convert to FASTA if needed
    - Double sequence if circular (for origin-crossing features)
    - Save metadata
    """
    # Read input file
    try:
        # Try GenBank format first
        records = list(SeqIO.parse(input_file, "genbank"))
        if records:
            record = records[0]
            file_format = "genbank"
        else:
            # Try FASTA format
            records = list(SeqIO.parse(input_file, "fasta"))
            if records:
                record = records[0]
                file_format = "fasta"
            else:
                raise ValueError("Could not parse input file as GenBank or FASTA")
    except Exception as e:
        raise ValueError(f"Error reading input file: {e}")
    
    if len(records) > 1:
        raise ValueError("Input file contains multiple sequences. Please provide a single sequence.")
    
    # Get sequence and validate
    seq_str = str(record.seq)
    seq_str = validate_sequence(seq_str)
    original_length = len(seq_str)
    
    # Double sequence if circular
    if not is_linear:
        seq_str = seq_str + seq_str
    
    # Create new record for output
    out_record = SeqRecord(
        Seq(seq_str),
        id=record.id if record.id else "sequence",
        name=record.name if record.name else "sequence",
        description="Prepared for pLannotate annotation"
    )
    
    # Write FASTA output
    with open(output_fasta, 'w') as f:
        SeqIO.write(out_record, f, "fasta")
    
    # Save metadata
    metadata = {
        "original_file": input_file,
        "original_format": file_format,
        "original_id": record.id,
        "original_name": record.name,
        "original_length": original_length,
        "is_linear": is_linear,
        "prepared_length": len(seq_str),
        "doubled": not is_linear
    }
    
    with open(output_metadata, 'w') as f:
        json.dump(metadata, f, indent=2)


if __name__ == "__main__":
    # Snakemake inputs/outputs
    prepare_sequence(
        input_file=snakemake.input[0],
        output_fasta=snakemake.output.fasta,
        output_metadata=snakemake.output.metadata,
        is_linear=snakemake.params.linear
    )