"""
Snakemake workflow to gather and process all database files for pLannotate.

This workflow recreates the databases.yml file and downloads/processes all required
database files including BLAST databases, DIAMOND databases, and Infernal models.

Usage:
    snakemake -s gather.smk --cores 4
    snakemake -s gather.smk --cores 4 --dry-run  # Preview workflow

Output:
    - plannotate-databases-v2.tar.gz
    - plannotate-databases-v2.tar.gz.sha256
"""

from datetime import datetime
from pathlib import Path
import os
import shlex
import sys

# Configuration - save to gather directory
DATA_DIR = "gathered_data"
BLAST_DB_DIR = f"{DATA_DIR}/BLAST_dbs"
DIAMOND_DB_DIR = f"{DATA_DIR}/diamond_dbs"  
INFERNAL_DB_DIR = f"{DATA_DIR}/infernal_dbs"
WORKFLOW_DIR = Path(workflow.basedir)
FPBASE_SCRIPT = shlex.quote(str(WORKFLOW_DIR / "fpbase/gather_fpbase.py"))
RFAM_SCRIPT = shlex.quote(str(WORKFLOW_DIR / "rfam/gather_rfam.py"))
SNAPGENE_SCRIPT = shlex.quote(str(WORKFLOW_DIR / "snapgene/gather_snapgene.py"))
SWISSPROT_DOWNLOAD = shlex.quote(
    str(WORKFLOW_DIR / "swissprot/download_fresh_swissprot.py")
)
SWISSPROT_PARSER = str(WORKFLOW_DIR / "swissprot/parse_swissprot_file.py")
BUILD_DIAMOND_SCRIPT = shlex.quote(
    str(WORKFLOW_DIR / "scripts/build_diamond_db.py")
)
COMBINE_TSV_SCRIPT = shlex.quote(str(WORKFLOW_DIR / "scripts/combine_tsv.py"))
CREATE_MANIFEST_SCRIPT = shlex.quote(
    str(WORKFLOW_DIR / "scripts/create_database_manifest.py")
)
CSV_TO_SQLITE_SCRIPT = shlex.quote(
    str(WORKFLOW_DIR / "scripts/csv_to_sqlite.py")
)
PACKAGE_BUNDLE_SCRIPT = shlex.quote(
    str(WORKFLOW_DIR / "scripts/package_database_bundle.py")
)

# Snakemake removes active Conda/virtualenv bin directories from rule PATHs.
# Resolve the project interpreter before jobs are launched so every Python rule
# uses the environment that invoked this workflow's dependencies.
ENV_PREFIX = os.environ.get("VIRTUAL_ENV") or os.environ.get("CONDA_PREFIX")
DEFAULT_PYTHON = Path(ENV_PREFIX) / "bin/python" if ENV_PREFIX else Path(sys.executable)
PYTHON = shlex.quote(config.get("python_executable", str(DEFAULT_PYTHON)))

# Ensure output directories exist
os.makedirs(BLAST_DB_DIR, exist_ok=True)
os.makedirs(DIAMOND_DB_DIR, exist_ok=True)
os.makedirs(INFERNAL_DB_DIR, exist_ok=True)

rule all:
    input:
        "plannotate-databases-v2.tar.gz",
        "plannotate-databases-v2.tar.gz.sha256"

rule package_database_bundle:
    """Package the runtime databases in the canonical 2.x archive layout."""
    input:
        snapgene_nhr=f"{BLAST_DB_DIR}/snapgene.nhr",
        snapgene_nin=f"{BLAST_DB_DIR}/snapgene.nin",
        snapgene_nsq=f"{BLAST_DB_DIR}/snapgene.nsq",
        snapgene_ndb=f"{BLAST_DB_DIR}/snapgene.ndb",
        snapgene_nog=f"{BLAST_DB_DIR}/snapgene.nog",
        snapgene_nos=f"{BLAST_DB_DIR}/snapgene.nos",
        snapgene_not=f"{BLAST_DB_DIR}/snapgene.not",
        snapgene_ntf=f"{BLAST_DB_DIR}/snapgene.ntf",
        snapgene_nto=f"{BLAST_DB_DIR}/snapgene.nto",
        snapgene_db=f"{BLAST_DB_DIR}/snapgene.db",
        fpbase=f"{DIAMOND_DB_DIR}/fpbase.dmnd",
        swissprot=f"{DIAMOND_DB_DIR}/swissprot.dmnd",
        fpbase_db=f"{DIAMOND_DB_DIR}/fpbase.db",
        swissprot_db=f"{DIAMOND_DB_DIR}/swissprot.db",
        rfam_cm=f"{INFERNAL_DB_DIR}/Rfam.cm",
        rfam_clanin=f"{INFERNAL_DB_DIR}/Rfam.clanin",
        rfam_i1f=f"{INFERNAL_DB_DIR}/Rfam.cm.i1f",
        rfam_i1i=f"{INFERNAL_DB_DIR}/Rfam.cm.i1i",
        rfam_i1m=f"{INFERNAL_DB_DIR}/Rfam.cm.i1m",
        rfam_i1p=f"{INFERNAL_DB_DIR}/Rfam.cm.i1p",
        manifest=f"{DATA_DIR}/database-manifest.json"
    output:
        archive="plannotate-databases-v2.tar.gz",
        checksum="plannotate-databases-v2.tar.gz.sha256"
    params:
        source=DATA_DIR
    log:
        "logs/package_database_bundle.log"
    shell:
        """
        {PYTHON} {PACKAGE_BUNDLE_SCRIPT} \
            --source {params.source} \
            --output {output.archive} \
        > {log} 2>&1
        """

rule gather_rfam:
    """Download and process Rfam database for Infernal searches."""
    output:
        cm_file=f"{INFERNAL_DB_DIR}/Rfam.cm",
        clan_file=f"{INFERNAL_DB_DIR}/Rfam.clanin",
        cm_i1f=f"{INFERNAL_DB_DIR}/Rfam.cm.i1f",
        cm_i1i=f"{INFERNAL_DB_DIR}/Rfam.cm.i1i",
        cm_i1m=f"{INFERNAL_DB_DIR}/Rfam.cm.i1m",
        cm_i1p=f"{INFERNAL_DB_DIR}/Rfam.cm.i1p",
        version_file=f"{INFERNAL_DB_DIR}/version.txt"
    log:
        "logs/gather_rfam.log"
    shell:
        """
        {PYTHON} {RFAM_SCRIPT} \
            --output-dir {INFERNAL_DB_DIR} \
        >& {log}
        """

rule gather_snapgene:
    """Download the canonical SnapGene BLAST database and metadata."""
    output:
        csv_file=f"{BLAST_DB_DIR}/snapgene.csv",
        blast_nhr=f"{BLAST_DB_DIR}/snapgene.nhr",
        blast_nin=f"{BLAST_DB_DIR}/snapgene.nin", 
        blast_nsq=f"{BLAST_DB_DIR}/snapgene.nsq",
        blast_ndb=f"{BLAST_DB_DIR}/snapgene.ndb",
        blast_nog=f"{BLAST_DB_DIR}/snapgene.nog",
        blast_nos=f"{BLAST_DB_DIR}/snapgene.nos",
        blast_not=f"{BLAST_DB_DIR}/snapgene.not",
        blast_ntf=f"{BLAST_DB_DIR}/snapgene.ntf",
        blast_nto=f"{BLAST_DB_DIR}/snapgene.nto",
        version_file=f"{BLAST_DB_DIR}/version.txt"
    log:
        "logs/gather_snapgene.log"
    shell:
        """
        {PYTHON} {SNAPGENE_SCRIPT} \
            --output-dir {BLAST_DB_DIR} \
        >& {log}
        """

rule create_snapgene_sqlite:
    """Convert SnapGene CSV to SQLite database."""
    input:
        csv_file=f"{BLAST_DB_DIR}/snapgene.csv"
    output:
        db_file=f"{BLAST_DB_DIR}/snapgene.db"
    log:
        "logs/create_snapgene_sqlite.log"
    shell:
        """
        {PYTHON} {CSV_TO_SQLITE_SCRIPT} \
            --input {input.csv_file} \
            --output {output.db_file} \
            --table snapgene \
            --rename Feature=name \
            --rename Type=type \
            --rename Description=blurb \
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
        {PYTHON} {FPBASE_SCRIPT} \
            --output {output.tsv_file} \
            --fasta {output.fasta_file} \
        >& {log}
        echo "downloaded $(date +%Y-%m-%d)" > {output.version_file}
        """

rule create_fpbase_sqlite:
    """Create the FPbase feature-description database (one source per file)."""
    input:
        fpbase_tsv=f"{DIAMOND_DB_DIR}/fpbase.tsv"
    output:
        db_file=f"{DIAMOND_DB_DIR}/fpbase.db"
    log:
        "logs/create_fpbase_sqlite.log"
    shell:
        """
        {PYTHON} {CSV_TO_SQLITE_SCRIPT} \
            --input {input.fpbase_tsv} \
            --output {output.db_file} \
            --table fpbase \
            --delimiter tab \
            --no-header \
            --add-type CDS \
        >& {log}
        """

rule create_swissprot_sqlite:
    """Create the Swiss-Prot feature-description database (one source per file)."""
    input:
        swissprot_tsv=f"{DIAMOND_DB_DIR}/swiss_description_verbose.tsv"
    output:
        db_file=f"{DIAMOND_DB_DIR}/swissprot.db"
    log:
        "logs/create_swissprot_sqlite.log"
    shell:
        """
        {PYTHON} {CSV_TO_SQLITE_SCRIPT} \
            --input {input.swissprot_tsv} \
            --output {output.db_file} \
            --table swissprot \
            --delimiter tab \
            --no-header \
        >& {log}
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
        {PYTHON} {BUILD_DIAMOND_SCRIPT} \
            --input {input.fasta_file} \
            --output {output.dmnd_file} \
            --format fasta \
        >& {log}
        """

checkpoint download_swissprot:
    """Download and split Swiss-Prot database."""
    output:
        fasta_file="gathered_data/swissprot/temp_data/uniprot_sprot.fasta",
        dat_file="gathered_data/swissprot/temp_data/uniprot_sprot.dat",
        version_file="gathered_data/swissprot/temp_data/version.txt",
        split_dir=directory("gathered_data/swissprot/temp_split")
    log:
        "logs/download_swissprot.log"
    shell:
        """
        {PYTHON} {SWISSPROT_DOWNLOAD} \
        >& {log}
        """

# Dynamic chunk processing using checkpoint pattern

rule process_swissprot_chunk:
    """Process individual Swiss-Prot chunk files."""
    input:
        chunk_file="gathered_data/swissprot/temp_split/chunk_{chunk_num}.dat",
        parser_script=SWISSPROT_PARSER,
        split_dir="gathered_data/swissprot/temp_split"  # Ensure download completes first
    output:
        verbose_tsv="gathered_data/swissprot/temp_results/swiss_description_verbose_chunk_{chunk_num}.tsv"
    log:
        "logs/process_chunk_{chunk_num}.log"
    shell:
        """
        mkdir -p gathered_data/swissprot/temp_results && \
        {PYTHON} {input.parser_script} \
            --output gathered_data/swissprot/temp_results/swiss_description_verbose_chunk_{wildcards.chunk_num}.tsv \
            gathered_data/swissprot/temp_split/chunk_{wildcards.chunk_num}.dat \
        >& {log}
        """

def get_chunk_numbers_from_split():
    """Get chunk numbers from split files."""
    chunk_dir = Path("gathered_data/swissprot/temp_split")
    
    if chunk_dir.exists():
        chunk_files = list(chunk_dir.glob("chunk_*.dat"))
        chunk_nums = []
        for f in chunk_files:
            try:
                num = f.stem.split('_')[1]  # Gets XX from chunk_XX
                chunk_nums.append(num)
            except (IndexError, ValueError):
                continue
        return sorted(chunk_nums)
    else:
        # Fallback for planning
        return [f"{i:02d}" for i in range(1, 18)]


def get_verbose_chunk_results(wildcards):
    """Get verbose chunk result files dynamically after checkpoint."""
    checkpoint_output = checkpoints.download_swissprot.get(**wildcards).output[0]
    chunk_nums = get_chunk_numbers_from_split()
    return expand("gathered_data/swissprot/temp_results/swiss_description_verbose_chunk_{chunk_num}.tsv",
                  chunk_num=chunk_nums)

rule combine_swissprot_results:
    """Combine all processed chunk results."""
    input:
        verbose_chunks=get_verbose_chunk_results
    output:
        verbose_tsv="gathered_data/swissprot/temp_results/swiss_description_verbose.tsv"
    log:
        "logs/combine_swissprot.log"
    shell:
        """
        {PYTHON} {COMBINE_TSV_SCRIPT} \
            --input {input.verbose_chunks} \
            --output {output.verbose_tsv} \
        >& {log}
        """

rule gather_swissprot:
    """Copy Swiss-Prot files to final location."""
    input:
        fasta_file="gathered_data/swissprot/temp_data/uniprot_sprot.fasta",
        verbose_description_file="gathered_data/swissprot/temp_results/swiss_description_verbose.tsv",
        version_file="gathered_data/swissprot/temp_data/version.txt"
    output:
        fasta_file=f"{DIAMOND_DB_DIR}/swissprot.fasta",
        verbose_description_file=f"{DIAMOND_DB_DIR}/swiss_description_verbose.tsv",
        version_file=f"{DIAMOND_DB_DIR}/swissprot_version.txt"
    log:
        "logs/gather_swissprot.log"
    shell:
        """
        cp {input.fasta_file} {output.fasta_file} && \
        cp {input.verbose_description_file} {output.verbose_description_file} && \
        cp {input.version_file} {output.version_file}
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

rule create_database_manifest:
    """Record database versions and checksums separately from search config."""
    input:
        rfam_version=f"{INFERNAL_DB_DIR}/version.txt",
        snapgene_version=f"{BLAST_DB_DIR}/version.txt",
        fpbase_version=f"{DIAMOND_DB_DIR}/fpbase_version.txt",
        swissprot_version=f"{DIAMOND_DB_DIR}/swissprot_version.txt",
        snapgene_nhr=f"{BLAST_DB_DIR}/snapgene.nhr",
        snapgene_nin=f"{BLAST_DB_DIR}/snapgene.nin",
        snapgene_nsq=f"{BLAST_DB_DIR}/snapgene.nsq",
        snapgene_ndb=f"{BLAST_DB_DIR}/snapgene.ndb",
        snapgene_nog=f"{BLAST_DB_DIR}/snapgene.nog",
        snapgene_nos=f"{BLAST_DB_DIR}/snapgene.nos",
        snapgene_not=f"{BLAST_DB_DIR}/snapgene.not",
        snapgene_ntf=f"{BLAST_DB_DIR}/snapgene.ntf",
        snapgene_nto=f"{BLAST_DB_DIR}/snapgene.nto",
        snapgene_db=f"{BLAST_DB_DIR}/snapgene.db",
        fpbase=f"{DIAMOND_DB_DIR}/fpbase.dmnd",
        swissprot=f"{DIAMOND_DB_DIR}/swissprot.dmnd",
        fpbase_db=f"{DIAMOND_DB_DIR}/fpbase.db",
        swissprot_db=f"{DIAMOND_DB_DIR}/swissprot.db",
        rfam_cm=f"{INFERNAL_DB_DIR}/Rfam.cm",
        rfam_clanin=f"{INFERNAL_DB_DIR}/Rfam.clanin",
        rfam_i1f=f"{INFERNAL_DB_DIR}/Rfam.cm.i1f",
        rfam_i1i=f"{INFERNAL_DB_DIR}/Rfam.cm.i1i",
        rfam_i1m=f"{INFERNAL_DB_DIR}/Rfam.cm.i1m",
        rfam_i1p=f"{INFERNAL_DB_DIR}/Rfam.cm.i1p"
    output:
        manifest=f"{DATA_DIR}/database-manifest.json"
    log:
        "logs/create_database_manifest.log"
    params:
        source=DATA_DIR
    shell:
        """
        {PYTHON} {CREATE_MANIFEST_SCRIPT} \
            --source {params.source} \
            --output {output.manifest} \
            --rfam-version-file {input.rfam_version} \
            --snapgene-version-file {input.snapgene_version} \
            --fpbase-version-file {input.fpbase_version} \
            --swissprot-version-file {input.swissprot_version} \
        >& {log}
        """
