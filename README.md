[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkab374.svg)](https://doi.org/10.1093/nar/gkab374)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/plannotate/README.html)


<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

Online Annotation
=================

pLannotate is web server for automatically annotating engineered plasmids.

Please visit http://plannotate.barricklab.org/


Local Installation
==================

To use pLannotate as a local server or a command line tool, please follow the installation instructions below.
### Quick install

The easiest way to install is via [conda](https://docs.conda.io/en/latest/):

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
```

Then activate the `plannotate` conda environment (`conda activate plannotate`) and proceed with using pLannotate (see **Using pLannotate locally** below).


### Installing from source

Conda is also recommended when installing from source.
Download the compressed source code from the [releases](https://github.com/barricklab/pLannotate/releases) page.
Uncompress the source code and move the directory to a location of your choice.

On the command line, navigate into the `pLannotate` folder.

Enter the following commands:
```
conda env create --name plannotate -f environment.yml
conda activate plannotate
pip install .[test]
```

After installation, run the following command to download the database files:
```
plannotate setupdb
```

Using pLannotate locally
=====
### Local server (GUI)

After installation, launch pLannotate as a local web server with:
```
plannotate streamlit
```

pLannotate should launch in your default web browser, or you may simply navigate to http://localhost:8501 in your web browser.

### Command Line Interface (batch mode)

To annotate FASTA or GenBank files and generate the interactive plasmid maps on the command line,
follow the above instructions to install pLannotate.

We can check the options using the following command:

`plannotate batch --help`

```
Usage: plannotate batch [OPTIONS]

  Annotates engineered DNA sequences, primarily plasmids. Accepts a FASTA file
  and outputs a gbk file with annotations, as well as an optional interactive
  plasmid map as an HTLM file.

Options:
  -i, --input TEXT      location of a FASTA or GBK file
  -o, --output TEXT     location of output folder. DEFAULT: current dir
  -f, --file_name TEXT  name of output file (do not add extension). DEFAULT:
                        input file name

  -s, --suffix TEXT     suffix appended to output files. Use '' for no suffix.
                        DEFAULT: '_pLann'

  -y, --yaml_file TEXT  path to YAML file for custom databases. DEFAULT:
                        builtin

  -l, --linear          enables linear DNA annotation
  -h, --html            creates an html plasmid map in specified path
  -c, --csv             creates a cvs file in specified path
  -d, --detailed        uses modified algorithm for a more-detailed search
                        with more false positives

  -x, --no_gbk          supresses GenBank output file
  --help                Show this message and exit.
  ```

Example usage:
```
plannotate batch -i ./plannotate/data/fastas/pUC19.fa --html --output ~/Desktop/ --file_name pLasmid
```

Custom databases can be added by supplying pLannotate a custom YAML file. To create the default YAML file, enter the following command:
```
plannotate yaml > plannotate_default.yaml
```

This configuration file can be edited to point to other external databases that you wish to use. When launching pLannotate, you can specify the path to your custom YAML file using the `--yaml_file` option.

### Using within Python

You can also directly import pLannotate as a Python module:

```python
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.resources import get_seq_record
from bokeh.io import show

# for inline plotting in jupyter
from bokeh.resources import INLINE
import bokeh.io
bokeh.io.output_notebook(INLINE)

seq = "tgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataggtctcaatccacgggtacgggtatggagaaacagtagagagttgcgataaaaagcgtcaggtagtatccgctaatcttatggataaaaatgctatggcatagcaaagtgtgacgccgtgcaaataatcaatgtggacttttctgccgtgattatagacacttttgttacgcgtttttgtcatggctttggtcccgctttgttacagaatgcttttaataagcggggttaccggtttggttagcgagaagagccagtaaaagacgcagtgacggcaatgtctgatgcaatatggacaattggtttcttgtaatcgttaatccgcaaataacgtaaaaacccgcttcggcgggtttttttatggggggagtttagggaaagagcatttgtcatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcgg"

# get pandas df of annotations
hits = annotate(seq, is_detailed = True, linear= True)

# get biopython SeqRecord object
seq_record = get_seq_record(hits, seq)

# show plot
show(get_bokeh(hits, linear=True))
```

This syntax will likely change in the future to be more user-friendly.

Testing
=====

Run the fast unit test suite with:

```bash
python -m pytest
```

Tests that require BLAST, DIAMOND, Infernal, and downloaded pLannotate databases
are marked as integration tests. They remain visible during pytest collection and
IDE discovery. They are deselected from normal command-line test runs by default,
but remain selectable in VSCode's Test UI. To run them locally from the command
line, install the test dependencies and databases, then pass `--run-integration`:

```bash
pip install .[test]
plannotate setupdb
python -m pytest --run-integration
```

The `test` extra includes `pytest-xdist`, so both default and integration test
runs can be parallelized:

```bash
python -m pytest -n auto
python -m pytest -n auto --run-integration
```

Each test has a timeout guard so a stuck external tool fails with a clear test
error instead of hanging indefinitely. Override the defaults with
`--test-timeout` and `--integration-timeout`, or use `0` to disable a timeout.

GitHub Actions runs the default unit suite first, then downloads databases and
runs the integration suite with `pytest -n auto --run-integration`.

The tests include a serialized annotation snapshot for the bundled example FASTA files in `plannotate/data/fastas`.
The snapshot is stored at `tests/test_data/example_fasta_annotations.json` and is checked by `tests/test_example_fasta_annotations.py`.
It records the cleaned annotation fields for each example sequence so future changes to annotation behavior are reviewed deliberately.
In addition to the default circular annotation mode for every bundled FASTA, the snapshot includes a few representative alternate-mode cases covering `linear=True`, `is_detailed=True`, and both together.
The snapshot check is split into one pytest case per FASTA file, allowing `pytest-xdist` to distribute the annotation work across workers.

If an annotation change is intentional, refresh the snapshot with:

```bash
PLANNOTATE_UPDATE_FASTA_ANNOTATION_SNAPSHOTS=1 python -m pytest tests/test_example_fasta_annotations.py --run-integration -q
```

About
=====
pLannotate was developed and is maintained by [Matt McGuffie](https://twitter.com/matt_mcguffie) at the [Barrick lab](https://barricklab.org/twiki/bin/view/Lab), University of Texas at Austin, Austin, Texas.
