# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pLannotate is a command-line tool and web server for annotating engineered plasmid sequences. The tool performs multi-method searches using BLAST, DIAMOND, and Infernal to identify features in DNA sequences and generate interactive visualizations.

### Core Architecture

The main application is organized into focused modules:

- **`plannotate/main.py`** - CLI entry point using Typer with three main commands:
  - `plannotate batch` - Main annotation command
  - `plannotate setupdb` - Database setup
  - `plannotate yaml` - Configuration export
- **`plannotate/annotate.py`** - Core annotation engine with multi-method search capabilities
- **`plannotate/models.py`** - Data models including `Construct` class and `Feature` dataclass
- **`plannotate/resources.py`** - Resource management, validation, and database configuration
- **`plannotate/bokeh_plot.py`** - Interactive visualization using Bokeh
- **`plannotate/infernal.py`** - RNA structure search functionality
- **`plannotate/logging_config.py`** - Centralized logging configuration

### Database Architecture

The tool uses multiple annotation databases configured via YAML:
- **SnapGene** - Curated plasmid features (BLAST nucleotide search)
- **Swiss-Prot** - Protein sequences (DIAMOND protein search)  
- **FPbase** - Fluorescent proteins (DIAMOND protein search)
- **Rfam** - RNA families (Infernal covariance model search)

Database locations and search parameters are defined in `plannotate/data/data/databases.yml`.

### Key Data Flow

1. Input validation via `resources.validate_file()`
2. Sequence processing through `Construct` class initialization
3. Multi-method annotation via `annotate.annotate()`
4. Feature scoring and filtering via `annotate.calculate()` and `annotate.clean()`
5. Output generation via `Construct.to_genbank()`, `Construct.to_html()`, `Construct.to_csv()`

## Development Commands

### Environment Setup
```bash
# Create conda environment from file
conda env create -f environment.yml
conda activate plannotate

# Install from source
python setup.py install

# Download required databases
plannotate setupdb
```

### Testing
```bash
# Run unit tests
python -m pytest tests/test_units.py -v

# Run specific test
python -m pytest tests/test_units.py::test_BLAST -v
```

### Code Quality
```bash
# Format code (use ruff - it's in the environment)
ruff format .

# Lint code  
ruff check .
```

### Example Usage
```bash
# Basic annotation
plannotate batch -i input/plasmid.fa -o output/ --html

# Detailed annotation with custom database
plannotate batch -i input/plasmid.fa -o output/ --detailed --yaml_file custom_db.yaml

# Linear DNA annotation
plannotate batch -i input/linear.fa --linear --csv
```

## Pipeline Integration

The project includes a Snakemake pipeline (`plannotate.smk`) for batch processing:
- Processes multiple input files automatically
- Supports various input formats (FASTA, GenBank)
- Generates consistent outputs across samples
- Requires `config.yaml` for configuration

## Python API

The tool can be imported and used programmatically:

```python
from plannotate.annotate import annotate
from plannotate.models import Construct

# Direct annotation
hits_df = annotate(sequence_string, is_detailed=True, linear=False)

# Full pipeline with outputs
construct = Construct(seq=sequence, linear=False, detailed=True)
gbk_content = construct.to_genbank()
html_content = construct.to_html()
csv_df = construct.to_csv()
```

## Important Implementation Notes

### Database Dependencies
- Databases must be downloaded via `plannotate setupdb` before first use
- Custom databases can be configured by modifying the YAML configuration
- External tools required: BLAST+, DIAMOND, Infernal, tRNAscan-SE

### File Format Support
- Input: FASTA (.fa, .fasta, .fas, .fna), GenBank (.gbk, .gb, .gbf, .gbff)
- Output: GenBank, HTML (interactive Bokeh plots), CSV

### Performance Considerations
- `--detailed` mode increases sensitivity but also false positives
- Large sequences may require significant processing time
- DIAMOND searches are faster than BLAST for protein sequences

### Logging
- Centralized logging configuration in `logging_config.py`
- Use `--verbose` flag for debug-level logging
- All major operations are logged for troubleshooting

## Bash Command Formatting Style

When writing shell commands in Snakemake rules, follow these strict guidelines:

### Required Format
```bash
# Correct format - call Python scripts, never inline code
python3 scripts/search_database.py \
    --input {input.seq} \
    --output {output.hits} \
    --database {params.db_name} \
>& {log}
```

### Formatting Rules
- **NEVER use `python3 -c` with inline code in shell directives**
- All Python code must be in separate scripts in `scripts/` directory
- Use backslashes (`\`) for line continuation
- Align parameters vertically with proper indentation
- Place each parameter on its own line for long commands
- Use `>& {log}` for log redirection (combines stdout and stderr)
- Keep commands readable and well-structured

### Script Requirements
- All scripts in `scripts/` should use proper `argparse` for command-line arguments
- Scripts should have clear error handling and logging
- Each script should have a single, focused purpose
- Use docstrings and follow the project's Python style guidelines

This architecture emphasizes modularity, with each component having a single focused responsibility and clear interfaces between modules.