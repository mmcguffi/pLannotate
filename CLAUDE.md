# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pLannotate is a command-line tool and web server for annotating engineered plasmid sequences. The tool performs multi-method searches using BLAST, DIAMOND, and Infernal to identify features in DNA sequences and generate interactive visualizations.

### Core Architecture

The main application is organized into focused modules:

- **`plannotate/main.py`** - CLI entry point using Typer with commands:
  - `plannotate batch` - Main annotation command
  - `plannotate setupdb` - Database setup
  - `plannotate yaml` - Configuration export
  - `plannotate databases` - Print the installed database manifest
  - `plannotate streamlit` - Launch the optional web app (requires the `server` extra)
- **`plannotate/annotate.py`** - Candidate collection and final annotation pipeline
- **`plannotate/models.py`** - `Construct`, `Feature`, conversions, and output methods
- **`plannotate/_tools/`** - BLAST, DIAMOND, and Infernal integrations
- **`plannotate/_concurrency.py`** - Core allocation and ordered thread-pool execution
- **`plannotate/_package_data.py`** - Packaged assets and database configuration
- **`plannotate/bokeh_plot.py`** - Plot preparation, geometry, and Bokeh rendering
- **`plannotate/streamlit_app.py`** - Optional Streamlit web front end, built on `Construct`

### Database Architecture

The tool uses multiple annotation databases configured via YAML:
- **SnapGene** - Curated plasmid features (BLAST nucleotide search)
- **Swiss-Prot** - Protein sequences (DIAMOND protein search)  
- **FPbase** - Fluorescent proteins (DIAMOND protein search)
- **Rfam** - RNA families (Infernal covariance model search)

Database locations and search parameters are defined in `plannotate/data/data/databases.yml`.

### Key Data Flow

1. `validation.validate_file()` reads one FASTA or GenBank record.
2. `Construct` calls `annotate.annotate()`.
3. `annotate.annotate()` runs configured sources and finalizes their candidates.
4. `_filter.filter_and_clean_hits()` scores hits and resolves overlaps.
5. `Construct` exports GenBank, CSV, or optional Bokeh HTML.

## Development Commands

### Environment Setup
```bash
# Create conda environment from file
conda env create -f environment.yml
conda activate plannotate

# Install the package and development dependencies
pip install -e '.[test,lint]'

# Download required databases
plannotate setupdb
```

### Testing
```bash
# Fast suite
pytest

# Include external tools and downloaded databases
pytest --run-integration

# Static checks
mypy plannotate tests
ruff check .
ruff format --check .
```

### Code Quality
```bash
# Format code
ruff format .

# Lint code  
ruff check .
```

### Python comments

Python comments should be clear and concise, following the project's style guidelines. Use docstrings for module, class, and function documentation. Inline comments should explain complex logic or decisions; explain they "why" rather than the "what" of the code. Start inline comments with a lowercase letter and keep them brief. Use `# TODO` for tasks that need to be addressed later, and `# NOTE` for particularly thorny or important points that may not be immediately obvious.

### Example Usage
```bash
# Basic annotation
plannotate batch -i input/plasmid.fa -o output/ --html

# Detailed annotation with custom database
plannotate batch -i input/plasmid.fa -o output/ --detailed --yaml-file custom_db.yaml

# Linear DNA annotation
plannotate batch -i input/linear.fa --linear --csv
```

## Database Builds

Runtime annotation uses only the Python standard library for scheduling. Snakemake is
an optional database-build dependency: `pip install -e '.[databases]'`. The build
workflow lives under `plannotate/gather_databases/`; it is not part of the runtime
annotation path.

The supported Python entry point for rebuilding the bundle is
`plannotate.build_databases(output_directory, cores=...)`.

## Python API

The tool can be imported and used programmatically:

```python
from plannotate.annotate import annotate
from plannotate import Construct

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
- External tools required: BLAST+, DIAMOND, and Infernal

### File Format Support
- Input: FASTA (.fa, .fasta, .fas, .fna), GenBank (.gbk, .gb, .gbf, .gbff)
- Output: GenBank, HTML (interactive Bokeh plots), CSV

### Performance Considerations
- `--detailed` mode increases sensitivity but also false positives
- Large sequences may require significant processing time
- DIAMOND searches are faster than BLAST for protein sequences

### Logging
- Library modules use standard `logging.getLogger(__name__)` loggers.
- The CLI configures the `plannotate` logger; `--verbose` enables debug output.

## Bash Command Formatting Style

When writing shell commands in Snakemake rules, follow these strict guidelines:

### Required Format
```bash
# Correct format - call Python scripts, never inline code
python3 example.py \
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
