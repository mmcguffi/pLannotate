# pLannotate Snakemake Pipeline

This is a refactored version of pLannotate that uses Snakemake as the workflow engine, breaking the annotation process into modular, parallelizable steps.

## Overview

The pipeline performs the following steps:

1. **Sequence Preparation**: Validates input sequence and doubles it if circular
2. **Database Searches**: Runs BLAST/Diamond/Infernal searches in parallel across all databases
3. **Feature Details**: Adds descriptions and types to hits
4. **Score Calculation**: Computes scores based on identity, coverage, and database priority
5. **Database Merging**: Combines results from all databases
6. **Overlap Cleaning**: Removes overlapping features based on score priority
7. **Fragment Detection**: Identifies partial/fragment features
8. **Final Formatting**: Prepares final annotations with sequence extraction
9. **Output Generation**: Creates GenBank and CSV files, optionally with Bokeh visualization

## Installation

### Prerequisites

- Snakemake (install with `conda install -c bioconda snakemake` or `pip install snakemake`)
- BLAST+ tools
- Diamond (for protein searches)
- Infernal (for RNA searches)
- ripgrep (`rg` command for compressed file searching)
- Python packages: pandas, numpy, biopython, bokeh

### Setup

1. Ensure all pLannotate databases are downloaded (use existing `plannotate setupdb` command)
2. Configure the pipeline by editing `config.yaml`

## Usage

### Basic Usage

```bash
# Run the pipeline
snakemake -c 4  # Use 4 cores

# Dry run to see what will be executed
snakemake -n

# Generate a visualization of the workflow
snakemake --dag | dot -Tpng > workflow.png
```

### Configuration

Edit `config.yaml` to set:

- `input_sequence`: Path to your input FASTA or GenBank file
- `yaml_file`: Path to database configuration (usually default is fine)
- `linear`: Set to `true` for linear DNA, `false` for circular (default)
- `detailed`: Set to `true` for more detailed search with more false positives
- `work_dir`: Directory for intermediate files
- `results_dir`: Directory for final output files

### Example Config

```yaml
input_sequence: "my_plasmid.fasta"
yaml_file: "../data/databases.yml"
work_dir: "work"
results_dir: "results"
linear: false
detailed: false
threads: 4
```

## Output Files

The pipeline generates:

- `results/{sample}_annotations.csv`: User-friendly annotation table
- `results/{sample}_annotations.gbk`: GenBank file with all annotations
- `results/{sample}_plasmid_map.html`: Interactive Bokeh visualization (if rule is run)

## Intermediate Files

The pipeline creates intermediate files in the work directory:

- `work/prepared/`: Validated and prepared sequences
- `work/blast/`: Raw BLAST/Diamond/Infernal results
- `work/details/`: Hits with added feature descriptions
- `work/scored/`: Hits with calculated scores
- `work/merged/`: Combined results from all databases
- `work/cleaned/`: Results after overlap removal
- `work/fragments/`: Results with fragment detection
- `work/final/`: Final formatted annotations

## Advantages of Snakemake Refactoring

1. **Parallelization**: Database searches run in parallel, significantly speeding up annotation
2. **Modularity**: Each step is isolated, making debugging and modifications easier
3. **Reproducibility**: Snakemake tracks dependencies and only reruns necessary steps
4. **Scalability**: Easy to run on clusters or cloud environments
5. **Transparency**: Clear visualization of the workflow and dependencies

## Advanced Usage

### Running Specific Steps

```bash
# Only run BLAST searches
snakemake work/blast/my_plasmid_snapgene_raw.tsv

# Generate only the CSV output
snakemake results/my_plasmid_annotations.csv
```

### Cluster Execution

```bash
# Run on a SLURM cluster
snakemake --cluster "sbatch -p short -t 30:00" -j 10
```

### Custom Databases

To use custom databases, modify the YAML file specified in the config and ensure the database files are accessible.

## Troubleshooting

1. **Missing dependencies**: Ensure all required tools (BLAST, Diamond, Infernal) are in your PATH
2. **Database errors**: Run `plannotate setupdb` to ensure databases are downloaded
3. **Memory issues**: Reduce parallelization by using fewer cores
4. **Debugging**: Check log files in `.snakemake/log/`

## Integration with Original pLannotate

This pipeline can be integrated with the original pLannotate by:

1. Creating a new command-line option to use Snakemake backend
2. Wrapping Snakemake execution in the existing Python interface
3. Converting between the existing data structures and file-based approach

Example integration point in `pLannotate.py`:

```python
@main.command("pipeline")
@click.option("--input", "-i", required=True)
@click.option("--output", "-o", default="results")
@click.option("--cores", "-c", default=4)
def main_pipeline(input, output, cores):
    """Run pLannotate using Snakemake pipeline"""
    # Create config file
    # Run snakemake
    # Return results
```