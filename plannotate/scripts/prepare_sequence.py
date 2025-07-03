"""
Prepare and validate input sequence for pLannotate pipeline
"""

import json
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Get Snakemake variables
input_seq = snakemake.params.input_seq
is_linear = snakemake.params.linear
seq_file = snakemake.output.seq_file
seq_info_file = snakemake.output.seq_info

# Validate sequence using BioPython
with NamedTemporaryFile(mode='w', suffix='.fasta') as tmp:
    SeqIO.write(
        SeqRecord(Seq(input_seq), name="pLannotate", 
                 annotations={"molecule_type": "DNA"}),
        tmp.name, "fasta"
    )
    record = list(SeqIO.parse(tmp.name, "fasta"))[0]

# Get original sequence
original_seq = str(record.seq)

# Prepare sequence (double if circular)
if is_linear:
    prepared_seq = original_seq
else:
    prepared_seq = original_seq + original_seq

# Write prepared sequence
SeqIO.write(
    SeqRecord(Seq(prepared_seq), id="query", description=""),
    seq_file, "fasta"
)

# Save sequence info
seq_info = {
    "original_length": len(original_seq),
    "prepared_length": len(prepared_seq),
    "is_linear": is_linear,
    "original_seq": original_seq
}

with open(seq_info_file, 'w') as f:
    json.dump(seq_info, f, indent=2)