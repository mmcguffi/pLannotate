# Database Gathering Workflow

This directory contains a Snakemake workflow (`gather.smk`) to recreate the pLannotate databases and configuration file.

## Overview

The workflow downloads and processes all required databases:

- **Rfam**: RNA families for Infernal covariance model searches
- **SnapGene**: Plasmid features for BLAST nucleotide searches  
- **FPbase**: Fluorescent proteins for DIAMOND protein searches
- **Swiss-Prot**: Protein sequences for DIAMOND protein searches

## Files

### Main Workflow
- `gather.smk` - Snakemake workflow orchestrating all database gathering
- `validate_workflow.py` - Validation script to check workflow structure

### Helper Scripts
- `scripts/build_diamond_db.py` - Builds DIAMOND databases from TSV/FASTA files
- `scripts/create_databases_yml.py` - Creates the databases.yml configuration file from template

### Templates
- `templates/databases_template.yml` - YAML template with placeholders for version information

### Database-Specific Scripts
- `rfam/gather_rfam.py` - Downloads Rfam covariance models
- `snapgene/gather_snapgene.py` - Creates CSV and FASTA from existing SnapGene database
- `fpbase/gather_fpbase.py` - Downloads FPbase protein data via GraphQL API
- `fpbase/gather_fpbase.sh` - Shell wrapper for FPbase download
- `swissprot/download_fresh_swissprot.py` - Downloads Swiss-Prot database
- `swissprot/run_full_workflow.sh` - Complete Swiss-Prot processing workflow

## Usage

### Prerequisites
```bash
# Required external tools
conda install -c bioconda blast diamond-aligner infernal
# Or ensure these are available: makeblastdb, diamond, cmpress

# Required Python packages
pip install requests pandas biopython
```

### Running the Workflow
```bash
# Preview what will be executed
snakemake -s gather.smk --dry-run

# Run the full workflow (adjust cores as needed)
snakemake -s gather.smk --cores 4

# Run individual database gathering
snakemake -s gather.smk --cores 1 gather_rfam
snakemake -s gather.smk --cores 1 gather_snapgene
snakemake -s gather.smk --cores 1 gather_fpbase  
snakemake -s gather.smk --cores 1 gather_swissprot
```

### Validation
```bash
# Check workflow structure and scripts
python validate_workflow.py
```

## Outputs

The workflow creates:

```
../data/data/
├── databases.yml              # Updated configuration file
├── BLAST_dbs/                 # SnapGene BLAST databases
│   ├── snapgene.nhr
│   ├── snapgene.nin
│   ├── snapgene.nsq
│   ├── snapgene.csv           # Feature metadata (display_name, unique_name, details)
│   ├── snapgene.fasta         # Sequences with unique names as headers
│   └── version.txt
├── diamond_dbs/               # DIAMOND protein databases  
│   ├── fpbase.dmnd
│   ├── fpbase.tsv             # Protein metadata (name, slug, blurb)
│   ├── fpbase.fasta           # Protein sequences with slug headers
│   ├── swissprot.dmnd
│   ├── swissprot.fasta
│   ├── swiss_description.csv.gz
│   ├── fpbase_version.txt
│   └── swissprot_version.txt
└── infernal_dbs/              # Rfam covariance models
    ├── Rfam.cm
    ├── Rfam.clanin
    └── version.txt
```

## Notes

- The workflow follows the project's Bash command formatting guidelines
- All Python code is in separate scripts (no inline `python -c` commands)
- Each database gathering step is independent and can be run separately
- Version information is captured for reproducibility
- Log files are created in `logs/` directory for troubleshooting