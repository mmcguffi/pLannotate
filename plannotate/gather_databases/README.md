# Database Gathering Workflow

This directory contains a Snakemake workflow (`gather.smk`) to recreate the
pLannotate database bundle.

## Overview

The workflow downloads and processes all required databases:

- **Rfam**: RNA families for Infernal covariance model searches
- **SnapGene**: Plasmid features for BLAST nucleotide searches  
- **FPbase**: Fluorescent proteins for DIAMOND protein searches
- **Swiss-Prot**: Protein sequences for DIAMOND protein searches

## Files

### Main Workflow
- `gather.smk` - Snakemake workflow orchestrating all database gathering

### Helper Scripts
- `scripts/build_diamond_db.py` - Builds DIAMOND databases from TSV/FASTA files
- `scripts/combine_tsv.py` - Combines Swiss-Prot metadata chunks
- `scripts/create_database_manifest.py` - Records source versions and checksums
- `scripts/csv_to_sqlite.py` - Creates indexed metadata databases
- `scripts/package_database_bundle.py` - Creates the deterministic runtime archive

### Database-Specific Scripts
- `rfam/gather_rfam.py` - Downloads Rfam covariance models
- `snapgene/gather_snapgene.py` - Creates CSV and FASTA from existing SnapGene database
- `fpbase/gather_fpbase.py` - Downloads FPbase protein data via GraphQL API
- `swissprot/download_fresh_swissprot.py` - Downloads Swiss-Prot database
- `swissprot/parse_swissprot_file.py` - Creates Swiss-Prot metadata rows

## Usage

### Prerequisites
```bash
conda install -c conda-forge -c bioconda blast diamond infernal
pip install -e '.[databases]'
```

### Running the Workflow

From Python:

```python
from plannotate import build_databases

archive = build_databases("database-build", cores=4)
```

Or invoke Snakemake directly from this directory:

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

### Install and test the generated bundle

```bash
archive="$PWD/plannotate-databases-v2.tar.gz"
checksum=$(awk '{print $1}' "$archive.sha256")

cd ../../
PLANNOTATE_DATABASE_URL="file://$archive" \
PLANNOTATE_DATABASE_SHA256="$checksum" \
plannotate setupdb --force

python -m pytest --run-integration
```

## Outputs

The workflow creates a deterministic `plannotate-databases-v2.tar.gz` archive
and matching `.sha256` file. Search configuration remains in the package YAML;
database provenance is recorded in `database-manifest.json` inside the archive.
The archive has the runtime layout below:

```
plannotate-databases-v2.tar.gz
├── database-manifest.json     # Versions and per-file checksums
├── BLAST_dbs/                 # SnapGene BLAST databases
│   ├── snapgene.nhr
│   ├── snapgene.nin
│   ├── snapgene.nsq
│   └── descriptions.db        # Feature metadata
├── diamond_dbs/               # DIAMOND protein databases  
│   ├── fpbase.dmnd
│   ├── swissprot.dmnd
│   └── descriptions.db        # FPbase and Swiss-Prot metadata
└── infernal_dbs/              # Rfam covariance models
    ├── Rfam.cm
    ├── Rfam.clanin
    └── Rfam.cm.i1*
```

## Notes

- The workflow follows the project's Bash command formatting guidelines
- All Python code is in separate scripts (no inline `python -c` commands)
- Each database gathering step is independent and can be run separately
- Version information is captured for reproducibility
- Log files are created in `logs/` directory for troubleshooting
- Configure the repository variables `PLANNOTATE_DATABASE_URL` and
  `PLANNOTATE_DATABASE_SHA256` before publishing a release. The release workflow
  verifies and uploads this canonical archive.
