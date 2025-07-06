# pLannotate Project Guidelines for Claude

## Project Overview

This is a command line tool to annotate plasmid sequences with the following design principles:

- **Fast**: Efficient annotation of plasmid sequences
- **Modular**: Clean separation of concerns across modules
- **Multi-format**: Handle multiple file formats (FASTA, GenBank)
- **User-friendly**: Easy to use with good UX
- **Accurate**: Deep and accurate plasmid sequence annotation

### Core Architecture

Main modules in `./plannotate/`:
- `annotate.py` - Core annotation logic (BLAST/DIAMOND/Infernal searches)
- `bokeh_plot.py` - Interactive visualization using Bokeh
- `infernal.py` - Infernal RNA search functionality
- `logging_config.py` - Centralized logging configuration
- `main.py` - CLI entry point with Typer
- `models.py` - Data models and Construct class
- `resources.py` - Resource management and validation

### Pipeline Scripts

Additional pipeline scripts in `./scripts/`:
- `validate_input.py` - Input sequence validation
- `search_database.py` - Database search wrapper
- `combine_and_filter.py` - Result combination and filtering
- `generate_outputs.py` - Output file generation
- `generate_plot.py` - Plot generation

## Python Style Guidelines

### General Principles
- Write concise, technical responses with accurate Python examples
- Use functional, declarative programming; avoid classes where possible
- Prefer iteration and modularization over code duplication
- Use descriptive variable names with auxiliary verbs (e.g., `is_active`, `has_permission`)

### Code Formatting
- Use Ruff formatting and generally abide by PEP8
    - acutally use `ruff` -- it is in the environment
- snake_case for variables and functions
- ALL_CAPS for constants
- TitleCase/CamelCase for objects
- No bare exceptions
- No long lines
- Use Python type hints for all functions

### Code Organization
- Minimize indents (aim for 1 indent max, maybe 2 depending on context)
- Each `.py` file should have a singular, focused purpose
- Files should not exceed 400 lines (ideally much smaller)
- No bare code outside functions
- No global variables or code in global scope

### Documentation
- All modules (`.py`) need a docstring
- Functions that need it should have a docstring
- First line of docstrings must be ONE line only (PEP8)
- If you can't describe it in one line, it's not focused enough
- Code should be "self documenting" with clear variable names
- Comments describe purpose (why), not syntax (how)

### Design Patterns
- Follow single responsibility principle
- Write functions that can be tested
- Prefer "pure" functions that only modify what they're given
- Avoid side effects (no printing, no in-place modification in pure functions)
- Avoid custom classes if simpler/prebuilt solutions exist
- Use classes for organizing state/data and dataclasses for data organization
- No magic numbers or magic strings

### Code Quality
- Do not leave large amounts of commented code
- Critical code should be tested with formal pytests
- Guidelines are good, but be smart enough to know when to break them if it simplifies things

## Snakemake Guidelines

When working with Snakemake files (`.smk`):

### Structure
- Only use `shell` directive (no `run`, etc)
- `shell` should call CLI programs or Python scripts
- Structure bash in shell for legibility with proper line breaks

### Command Formatting
```bash
my_tool \
    -p1 100 \
    -p2 200 \
    --flag \
>& {log}
```

### Rule Requirements
Every rule should generally have:
- `input`
- `output` 
- `shell`
- `log`
- `conda` (depending on needs)

### Best Practices
- Each rule should be well-contained and do a single task
- Use external Python scripts rather than inline Python code
- Call existing project functions rather than duplicating logic
- Use proper log redirection with `>& {log}`

## Testing Requirements

### Test Execution
- When any changes are implemented, always run tests after changes
- Confirm that changes do not break existing functionality
- Run pytest: `python -m pytest tests/test_units.py -v`

### Test Development  
- When new functionality is added, write appropriate new tests
- If you don't have access to appropriate data for testing, ask the human
- Critical code should be tested with formal pytests
- Tests should cover edge cases and error conditions

### Test Data
- Use existing test data in `tests/test_data/` when possible
- Test files include various formats (`.fa`, `.fasta`, `.gbk`, etc.)
- Sample sequences available for testing annotation functionality

## Development Workflow

1. **Understand the existing codebase** - Read relevant modules before making changes
2. **Follow the existing patterns** - Use established functions and conventions
3. **Make focused changes** - One responsibility per function/module
4. **Test thoroughly** - Run existing tests and add new ones as needed
5. **Document appropriately** - Clear docstrings and self-documenting code
6. **Use existing utilities** - Leverage functions in `resources.py`, `annotate.py`, etc.

## Key Functions to Understand

### Core Annotation Pipeline
- `annotate.BLAST()` - Multi-method search (BLAST/DIAMOND/Infernal)
- `annotate.get_details()` - Feature detail retrieval
- `annotate.calculate()` - Score calculation
- `annotate.clean()` - Result filtering

### Resource Management
- `resources.validate_file()` - Input validation
- `resources.get_yaml()` - Database configuration
- `resources.get_yaml_path()` - Configuration file path

### Visualization
- `bokeh_plot.get_bokeh()` - Interactive plot generation
- `models.Construct` - Main data model with output methods

This project emphasizes clean, modular code with comprehensive testing and clear separation of concerns.