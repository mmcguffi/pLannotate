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

# Configuration
DATA_DIR = "../data/data"
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
    """Create SnapGene CSV and FASTA files from existing database."""
    output:
        csv_file=f"{BLAST_DB_DIR}/snapgene.csv",
        fasta_file=f"{BLAST_DB_DIR}/snapgene.fasta",
        version_file=f"{BLAST_DB_DIR}/snapgene_version.txt"
    log:
        "logs/gather_snapgene.log"
    shell:
        """
        python3 snapgene/gather_snapgene.py \
            --output-dir {BLAST_DB_DIR} \
        >& {log} && \
        echo "SnapGene processed $(date +%Y-%m-%d)" > {output.version_file}
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
        echo "FPbase downloaded $(date +%Y-%m-%d)" > ../{output.version_file}
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

rule gather_swissprot:
    """Download and process Swiss-Prot database."""
    output:
        fasta_file=f"{DIAMOND_DB_DIR}/swissprot.fasta",
        description_file=f"{DIAMOND_DB_DIR}/swiss_description.csv.gz",
        version_file=f"{DIAMOND_DB_DIR}/swissprot_version.txt"
    log:
        "logs/gather_swissprot.log"
    shell:
        """
        cd swissprot && \
        bash run_full_workflow.sh \
        >& ../{log} && \
        mv fresh_data/uniprot_sprot.fasta ../{output.fasta_file} && \
        mv fresh_results/swiss_description.csv.gz ../{output.description_file} && \
        echo "Swiss-Prot downloaded $(date +%Y-%m-%d)" > ../{output.version_file}
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
        snapgene_version=f"{BLAST_DB_DIR}/snapgene_version.txt", 
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