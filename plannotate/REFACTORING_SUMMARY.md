# pLannotate Snakemake Refactoring Summary

## Overview

The pLannotate annotation pipeline has been refactored from a monolithic script (`annotate.py`) into a modular Snakemake pipeline. This refactoring improves maintainability, extensibility, and performance through parallel execution.

## Key Changes

### 1. Modular Architecture

The original `annotate.py` has been split into focused modules:

- **`pipeline/search.py`**: Handles BLAST, DIAMOND, and Infernal searches
- **`pipeline/details.py`**: Retrieves feature descriptions from databases
- **`pipeline/process.py`**: Score calculation, cleaning, and fragment detection
- **`pipeline/combine.py`**: Combines results from multiple databases

### 2. Snakemake Workflow

A new `Snakefile` orchestrates the pipeline with rules for:
- Parallel database searches
- Score calculation
- Result combination
- Overlap removal
- Final annotation generation

### 3. Backward Compatibility

The new `annotate.py` provides a wrapper that maintains the original API:

```python
# Original API still works
from plannotate import annotate
results = annotate.annotate(sequence)
```

### 4. New Features

- **Parallel Execution**: Database searches run in parallel
- **Intermediate Files**: Results saved at each step for debugging
- **Configurable**: YAML configuration for pipeline parameters
- **Extensible**: Easy to add new databases or search methods

## File Structure

```
plannotate/
├── Snakefile              # Main pipeline orchestration
├── config.yaml            # Configuration template
├── annotate.py            # Backward-compatible wrapper (new)
├── annotate_original.py   # Original monolithic script (backup)
├── pipeline/              # Modular pipeline components
│   ├── __init__.py
│   ├── search.py          # Database search module
│   ├── details.py         # Feature detail retrieval
│   ├── process.py         # Score calculation and cleaning
│   ├── combine.py         # Result combination
│   └── README.md          # Pipeline documentation
└── examples/
    └── use_pipeline.py    # Usage examples
```

## Benefits

1. **Maintainability**: Each module has a single, clear responsibility
2. **Performance**: Parallel database searches reduce runtime
3. **Debugging**: Intermediate results saved for inspection
4. **Extensibility**: Easy to add new databases or modify pipeline steps
5. **Testing**: Modules can be unit tested independently
6. **Reproducibility**: Snakemake ensures consistent execution

## Migration Guide

For most users, no changes are needed. The original API is preserved:

```python
# Old way (still works)
from plannotate.annotate import annotate
results = annotate(sequence)

# New way with more control
from plannotate import annotate
results = annotate.annotate(
    sequence,
    threads=8,              # Use 8 parallel threads
    output_dir="results",   # Keep intermediate files
    keep_temp=True          # Don't delete temp directory
)
```

## Running Snakemake Directly

Advanced users can run the pipeline directly:

```bash
snakemake -s plannotate/Snakefile \
    --config input_sequence="ATCG..." \
    --cores 4
```

## Adding New Databases

1. Add configuration to `databases.yml`:
```yaml
my_database:
  location: "path/to/db"
  method: "blastn"
  priority: 5
  parameters: ["-evalue", "1e-5"]
  details:
    location: "path/to/descriptions.csv"
    default_type: "misc_feature"
```

2. The pipeline automatically includes the new database

## Testing

The refactored pipeline produces identical results to the original implementation while providing better performance and maintainability.