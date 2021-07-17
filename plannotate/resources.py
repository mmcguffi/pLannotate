
import os
from tempfile import NamedTemporaryFile

from Bio import SeqIO

valid_genbank_exts = ['.gbk', '.gb', '.gbf', '.gbff']
valid_fasta_exts = ['.fa', '.fasta']
maxPlasSize = 50000


def get_name_ext(file_loc):
    base = os.path.basename(file_loc)
    name = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1]
    return name,ext


def validate_file(file_loc, ext):
    if ext in valid_fasta_exts:
        #This catches errors on file uploads via Biopython
        fileloc = NamedTemporaryFile()
        record = list(SeqIO.parse(file_loc, "fasta"))
        try:
            record[0].annotations["molecule_type"] = "DNA"
        except IndexError:
            error = "Malformed fasta file --> please submit a fasta file in standard format"
            raise ValueError(error)
        SeqIO.write(record, fileloc.name, 'fasta')
        record = list(SeqIO.parse(fileloc.name, "fasta"))
        fileloc.close()

        if len(record)!=1:
            error = 'FASTA file contains many entries --> please submit a single FASTA file.'
            raise ValueError(error)

    elif ext in valid_genbank_exts:
        fileloc = NamedTemporaryFile()
        try:
            record = list(SeqIO.parse(file_loc, "gb"))[0]
        except IndexError:
            error = "Malformed Genbank file --> please submit a Genbank file in standard format"
            raise ValueError(error)
        # submitted_gbk = record # for combining -- not current imlementated
        SeqIO.write(record, fileloc.name, 'fasta')
        record = list(SeqIO.parse(fileloc.name, "fasta"))
        fileloc.close()
    
    else:
        error = 'must be a FASTA or GenBank file'
        raise ValueError(error)

    if len(record)!=1:
        error = 'FASTA file contains many entries --> please submit a single FASTA file.'
        raise ValueError(error)
        
    inSeq = str(record[0].seq)

    return inSeq


def validate_sequence(inSeq):
    IUPAC= 'GATCRYWSMKHBVDNgatcrywsmkhbvdn'
    if not set(inSeq).issubset(IUPAC):
        error = f'Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code'
        raise ValueError(error)

    if len(inSeq) > maxPlasSize:
        error = f'Are you sure this is an engineered plasmid? Entry size is too large -- must be {maxPlasSize} bases or less.'
        raise ValueError(error)
