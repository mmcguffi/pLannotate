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
        "databases_deployed.flag"

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
    """Copy SnapGene BLAST database and CSV from backup."""
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
        mkdir -p {BLAST_DB_DIR} && \
        cp ../data/temp_data_backup_csv/BLAST_dbs/snapgene.* {BLAST_DB_DIR}/ && \
        echo "Downloaded 2021-07-23" > {output.version_file} \
        >& {log}
        """

rule create_snapgene_sqlite:
    """Convert SnapGene CSV to SQLite database."""
    input:
        csv_file=f"{BLAST_DB_DIR}/snapgene.csv"
    output:
        db_file=f"{BLAST_DB_DIR}/descriptions.db"
    log:
        "logs/create_snapgene_sqlite.log"
    shell:
        """
        # Convert legacy format to standard format
        python3 -c "
import pandas as pd
df = pd.read_csv('{input.csv_file}')
df = df.rename(columns={{'Feature': 'name', 'Type': 'type', 'Description': 'blurb'}})
df.to_csv('{BLAST_DB_DIR}/snapgene_standard.csv', index=False, header=False)
" && \
        python3 scripts/csv_to_sqlite.py \
            --input {BLAST_DB_DIR}/snapgene_standard.csv \
            --output {output.db_file} \
            --table snapgene \
            --no-header \
        >& {log} && \
        rm {BLAST_DB_DIR}/snapgene_standard.csv
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

rule create_diamond_descriptions_sqlite:
    """Create SQLite database from FPbase TSV and Swiss-Prot verbose CSV."""
    input:
        fpbase_tsv=f"{DIAMOND_DB_DIR}/fpbase.tsv",
        swissprot_tsv=f"{DIAMOND_DB_DIR}/swiss_description_verbose.tsv"
    output:
        db_file=f"{DIAMOND_DB_DIR}/descriptions.db"
    log:
        "logs/create_diamond_sqlite.log"
    shell:
        """
        # Create FPbase table (no header, add CDS type)
        python3 scripts/csv_to_sqlite.py \
            --input {input.fpbase_tsv} \
            --output {output.db_file} \
            --table fpbase \
            --delimiter tab \
            --no-header \
            --add-type CDS \
        >& {log}
        
        # Add Swiss-Prot table (verbose version with CDS type in column 3)
        python3 scripts/csv_to_sqlite.py \
            --input {input.swissprot_tsv} \
            --output {output.db_file} \
            --table swissprot \
            --delimiter tab \
            --no-header \
            --append \
        >> {log} 2>&1
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

checkpoint download_swissprot:
    """Download and split Swiss-Prot database."""
    output:
        fasta_file="gathered_data/swissprot/temp_data/uniprot_sprot.fasta",
        dat_file="gathered_data/swissprot/temp_data/uniprot_sprot.dat",
        split_dir=directory("gathered_data/swissprot/temp_split")
    log:
        "logs/download_swissprot.log"
    shell:
        """
        python3 swissprot/download_fresh_swissprot.py \
        >& {log}
        """

# Dynamic chunk processing using checkpoint pattern

rule process_swissprot_chunk:
    """Process individual Swiss-Prot chunk files."""
    input:
        chunk_file="gathered_data/swissprot/temp_split/chunk_{chunk_num}.dat",
        parser_script="swissprot/parse_swissprot_file.py",
        split_dir="gathered_data/swissprot/temp_split"  # Ensure download completes first
    output:
        verbose_tsv="gathered_data/swissprot/temp_results/swiss_description_verbose_chunk_{chunk_num}.tsv"
    log:
        "logs/process_chunk_{chunk_num}.log"
    shell:
        """
        mkdir -p gathered_data/swissprot/temp_results && \
        python3 swissprot/parse_swissprot_file.py \
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
    run:
        import subprocess
        import os
        import pandas as pd
        
        # Get the input files 
        verbose_chunks = input.verbose_chunks
        
        # Create output directory
        os.makedirs("gathered_data/swissprot/temp_results", exist_ok=True)
        
        # Combine verbose results using pandas for proper TSV handling  
        verbose_dfs = []
        for chunk_file in verbose_chunks:
            df = pd.read_csv(chunk_file, header=None, sep='\t')
            verbose_dfs.append(df)
        
        if verbose_dfs:
            combined_verbose = pd.concat(verbose_dfs, ignore_index=True)
            combined_verbose.to_csv("gathered_data/swissprot/temp_results/swiss_description_verbose.tsv", 
                                  index=False, header=False, sep='\t')
        
        # No need to gzip - keep as TSV for direct use
        
        # Clean up chunk files
        for chunk_file in verbose_chunks:
            os.remove(chunk_file)

rule gather_swissprot:
    """Copy Swiss-Prot files to final location."""
    input:
        fasta_file="gathered_data/swissprot/temp_data/uniprot_sprot.fasta",
        verbose_description_file="gathered_data/swissprot/temp_results/swiss_description_verbose.tsv"
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

rule deploy_databases:
    """Copy essential database files to plannotate/data/data/databases/ for distribution."""
    input:
        databases_yml=f"{DATA_DIR}/databases.yml",
        # BLAST databases
        snapgene_nhr=f"{BLAST_DB_DIR}/snapgene.nhr",
        snapgene_nin=f"{BLAST_DB_DIR}/snapgene.nin", 
        snapgene_nsq=f"{BLAST_DB_DIR}/snapgene.nsq",
        snapgene_csv=f"{BLAST_DB_DIR}/snapgene.csv",
        blast_descriptions_db=f"{BLAST_DB_DIR}/descriptions.db",
        # DIAMOND databases  
        fpbase_dmnd=f"{DIAMOND_DB_DIR}/fpbase.dmnd",
        fpbase_tsv=f"{DIAMOND_DB_DIR}/fpbase.tsv",
        swissprot_dmnd=f"{DIAMOND_DB_DIR}/swissprot.dmnd",
        swiss_descriptions_tsv=f"{DIAMOND_DB_DIR}/swiss_description_verbose.tsv",
        diamond_descriptions_db=f"{DIAMOND_DB_DIR}/descriptions.db",
        # Infernal databases
        rfam_cm=f"{INFERNAL_DB_DIR}/Rfam.cm",
        rfam_clanin=f"{INFERNAL_DB_DIR}/Rfam.clanin"
    output:
        "databases_deployed.flag"
    log:
        "logs/deploy_databases.log"
    params:
        save_dir="../data/data/databases"
    shell:
        """
        # Create target directory
        mkdir -p {params.save_dir} && \

        # Copy databases.yml
        cp {input.databases_yml} {params.save_dir}/databases.yml && \

        # Copy BLAST databases
        cp {input.snapgene_nhr} {input.snapgene_nin} {input.snapgene_nsq} {params.save_dir}/ && \
        cp {input.snapgene_csv} {params.save_dir}/ && \
        cp {input.blast_descriptions_db} {params.save_dir}/snapgene_descriptions.db && \

        # Copy DIAMOND databases
        cp {input.fpbase_dmnd} {input.fpbase_tsv} {params.save_dir}/ && \
        cp {input.swissprot_dmnd} {params.save_dir}/ && \
        cp {input.swiss_descriptions_tsv} {params.save_dir}/ && \
        cp {input.diamond_descriptions_db} {params.save_dir}/swissprot_descriptions.db && \

        # Copy Infernal databases
        cp {input.rfam_cm} {input.rfam_clanin} {params.save_dir}/ && \

        # Create completion flag
        echo "Database deployment completed on $(date)" > {output}
        """