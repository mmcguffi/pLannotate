"""
Snakemake workflow to gather and process all database files for pLannotate.

This workflow recreates the databases.yml file and downloads/processes all required
database files including BLAST databases, DIAMOND databases, and Infernal models.

Usage:
    snakemake -s gather.smk --cores 4
    snakemake -s gather.smk --cores 4 --dry-run  # Preview workflow

Output:
    - ../data/data/databases.yml (updated configuration)
    - ../data/data/BLAST_dbs/* (SnapGene BLAST databases)  
    - ../data/data/diamond_dbs/* (Swiss-Prot and FPbase DIAMOND databases)
    - ../data/data/infernal_dbs/* (Rfam covariance models)
"""

from datetime import datetime
from pathlib import Path
import os

# Configuration - save to gather directory
DATA_DIR = "gathered_data"
BLAST_DB_DIR = f"{DATA_DIR}/BLAST_dbs"
DIAMOND_DB_DIR = f"{DATA_DIR}/diamond_dbs"  
INFERNAL_DB_DIR = f"{DATA_DIR}/infernal_dbs"

# Ensure output directories exist
os.makedirs(BLAST_DB_DIR, exist_ok=True)
os.makedirs(DIAMOND_DB_DIR, exist_ok=True)
os.makedirs(INFERNAL_DB_DIR, exist_ok=True)

rule all:
    input:
        f"{DATA_DIR}/databases.yml",
        f"{BLAST_DB_DIR}/snapgene.nhr",
        f"{BLAST_DB_DIR}/snapgene.csv",
        f"{DIAMOND_DB_DIR}/fpbase.dmnd",
        f"{DIAMOND_DB_DIR}/swissprot.dmnd", 
        f"{INFERNAL_DB_DIR}/Rfam.cm"

rule gather_rfam:
    """Download and process Rfam database for Infernal searches."""
    output:
        cm_file=f"{INFERNAL_DB_DIR}/Rfam.cm",
        clan_file=f"{INFERNAL_DB_DIR}/Rfam.clanin",
        version_file=f"{INFERNAL_DB_DIR}/version.txt"
    log:
        "logs/gather_rfam.log"
    shell:
        """
        python3 rfam/gather_rfam.py \
            --output-dir {INFERNAL_DB_DIR} \
        >& {log}
        """

rule gather_snapgene:
    """Download SnapGene BLAST database and CSV from GitHub release."""
    output:
        csv_file=f"{BLAST_DB_DIR}/snapgene.csv",
        blast_nhr=f"{BLAST_DB_DIR}/snapgene.nhr",
        blast_nin=f"{BLAST_DB_DIR}/snapgene.nin", 
        blast_nsq=f"{BLAST_DB_DIR}/snapgene.nsq",
        version_file=f"{BLAST_DB_DIR}/version.txt"
    log:
        "logs/gather_snapgene.log"
    shell:
        """
        python3 snapgene/gather_snapgene.py \
            --output-dir {BLAST_DB_DIR} \
        >& {log}
        """

rule gather_fpbase:
    """Download and process FPbase database for DIAMOND searches."""
    output:
        tsv_file=f"{DIAMOND_DB_DIR}/fpbase.tsv",
        fasta_file=f"{DIAMOND_DB_DIR}/fpbase.fasta",
        version_file=f"{DIAMOND_DB_DIR}/fpbase_version.txt"
    log:
        "logs/gather_fpbase.log"
    shell:
        """
        cd fpbase && \
        bash gather_fpbase.sh \
        >& ../{log} && \
        mv fpbase-proteins_*.tsv ../{output.tsv_file} && \
        mv fpbase-proteins_*.fasta ../{output.fasta_file} && \
        echo "downloaded $(date +%Y-%m-%d)" > ../{output.version_file}
        """

rule build_fpbase_diamond_db:
    """Build DIAMOND database from FPbase FASTA file."""
    input:
        fasta_file=f"{DIAMOND_DB_DIR}/fpbase.fasta"
    output:
        dmnd_file=f"{DIAMOND_DB_DIR}/fpbase.dmnd"
    log:
        "logs/build_fpbase_diamond.log"
    shell:
        """
        python3 scripts/build_diamond_db.py \
            --input {input.fasta_file} \
            --output {output.dmnd_file} \
            --format fasta \
        >& {log}
        """

rule download_swissprot:
    """Download and split Swiss-Prot database."""
    output:
        fasta_file="swissprot/fresh_data/uniprot_sprot.fasta",
        dat_file="swissprot/fresh_data/uniprot_sprot.dat",
        split_dir=directory("swissprot/fresh_split")
    log:
        "logs/download_swissprot.log"
    shell:
        """
        cd swissprot && \
        python3 download_fresh_swissprot.py \
        >& ../{log}
        """

# Get list of chunk files dynamically
def get_chunk_files(wildcards):
    chunk_dir = Path("swissprot/fresh_split")
    if chunk_dir.exists():
        return [str(f) for f in chunk_dir.glob("chunk_*.dat")]
    else:
        # Return expected files for planning (will be created by download_swissprot)
        return [f"swissprot/fresh_split/chunk_{i:02d}.dat" for i in range(1, 18)]

rule process_swissprot_chunk:
    """Process individual Swiss-Prot chunk files."""
    input:
        chunk_file="swissprot/fresh_split/chunk_{chunk_num}.dat",
        parser_script="swissprot/parse_swissprot_file.py",
        split_dir="swissprot/fresh_split"  # Ensure download completes first
    output:
        normal_csv="swissprot/temp_results/swiss_description_chunk_{chunk_num}.csv",
        verbose_csv="swissprot/temp_results/swiss_description_verbose_chunk_{chunk_num}.csv"
    log:
        "logs/process_chunk_{chunk_num}.log"
    shell:
        """
        mkdir -p swissprot/temp_results && \
        cd swissprot && \
        
        # Process normal version
        python3 parse_swissprot_file.py chunk_{wildcards.chunk_num}.dat 2>/dev/null && \
        mv swiss_description.csv temp_results/swiss_description_chunk_{wildcards.chunk_num}.csv
        
        # Process verbose version  
        python3 parse_swissprot_file.py -v chunk_{wildcards.chunk_num}.dat 2>/dev/null && \
        mv swiss_description.csv temp_results/swiss_description_verbose_chunk_{wildcards.chunk_num}.csv
        
        >& ../{log}
        """

rule combine_swissprot_results:
    """Combine all processed chunk results."""
    input:
        normal_chunks=expand("swissprot/temp_results/swiss_description_chunk_{chunk_num}.csv", 
                           chunk_num=[f"{i:02d}" for i in range(1, 18)]),
        verbose_chunks=expand("swissprot/temp_results/swiss_description_verbose_chunk_{chunk_num}.csv",
                            chunk_num=[f"{i:02d}" for i in range(1, 18)])
    output:
        normal_csv="swissprot/fresh_results/swiss_description.csv.gz",
        verbose_csv="swissprot/fresh_results/swiss_description_verbose.csv.gz"
    log:
        "logs/combine_swissprot.log"
    shell:
        """
        mkdir -p swissprot/fresh_results && \
        
        # Combine normal results
        cat {input.normal_chunks} > swissprot/fresh_results/swiss_description.csv && \
        gzip swissprot/fresh_results/swiss_description.csv
        
        # Combine verbose results
        cat {input.verbose_chunks} > swissprot/fresh_results/swiss_description_verbose.csv && \
        gzip swissprot/fresh_results/swiss_description_verbose.csv
        
        # Clean up temp files
        rm -rf swissprot/temp_results
        
        >& {log}
        """

rule gather_swissprot:
    """Copy Swiss-Prot files to final location."""
    input:
        fasta_file="swissprot/fresh_data/uniprot_sprot.fasta",
        description_file="swissprot/fresh_results/swiss_description.csv.gz"
    output:
        fasta_file=f"{DIAMOND_DB_DIR}/swissprot.fasta",
        description_file=f"{DIAMOND_DB_DIR}/swiss_description.csv.gz",
        version_file=f"{DIAMOND_DB_DIR}/swissprot_version.txt"
    log:
        "logs/gather_swissprot.log"
    shell:
        """
        cp {input.fasta_file} {output.fasta_file} && \
        cp {input.description_file} {output.description_file} && \
        echo "Release $(date +%Y_%m)" > {output.version_file}
        """

rule build_swissprot_diamond_db:
    """Build DIAMOND database from Swiss-Prot FASTA file."""
    input:
        fasta_file=f"{DIAMOND_DB_DIR}/swissprot.fasta"
    output:
        dmnd_file=f"{DIAMOND_DB_DIR}/swissprot.dmnd"
    log:
        "logs/build_swissprot_diamond.log"
    shell:
        """
        diamond makedb \
            --in {input.fasta_file} \
            --db {output.dmnd_file} \
        >& {log}
        """

rule create_databases_yml:
    """Create the databases.yml configuration file."""
    input:
        rfam_version=f"{INFERNAL_DB_DIR}/version.txt",
        snapgene_version=f"{BLAST_DB_DIR}/version.txt", 
        fpbase_version=f"{DIAMOND_DB_DIR}/fpbase_version.txt",
        swissprot_version=f"{DIAMOND_DB_DIR}/swissprot_version.txt",
        template="templates/databases_template.yml"
    output:
        f"{DATA_DIR}/databases.yml"
    log:
        "logs/create_databases_yml.log"
    shell:
        """
        python3 scripts/create_databases_yml.py \
            --rfam-version {input.rfam_version} \
            --snapgene-version {input.snapgene_version} \
            --fpbase-version {input.fpbase_version} \
            --swissprot-version {input.swissprot_version} \
            --template {input.template} \
            --output {output} \
        >& {log}
        """