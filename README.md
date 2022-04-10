[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkab374.svg)](https://doi.org/10.1093/nar/gkab374)


<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

Online Annotation
=================

pLannotate is web server for automatically annotating engineered plasmids.

Please visit http://plannotate.barricklab.org/


Local Installation
==================
To use pLannotate as a local server or a command line tool, please follow the installation instructions below.
### Quick install

The easiest way to install is via [conda](https://docs.conda.io/en/latest/), or for faster installation, [mamba](https://github.com/mamba-org/mamba):

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
```
or
```bash
mamba create -n plannotate -c conda-forge -c bioconda plannotate
```

Then activate the `plannotate` conda environment (`conda activate plannotate`) and proceed with using pLannotate (see **Using pLannotate locally** below).


### Installing from source
Installing from source also requires conda (or mamba), therefore the above method is recommended. If you still wish to install from source, download the compressed source code from the [releases](https://github.com/barricklab/pLannotate/releases) page. Uncompress the source code and move the directory to a location of your choice.

On the command line, navigate into the `pLannotate` folder.

Enter the following commands:
```
conda env create -f environment.yml
conda activate plannotate
python setup.py install
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
  -i, --input TEXT      location of a FASTA or GBK file; < 50,000 bases
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

About
=====
pLannotate is currently developed by [Matt McGuffie](https://twitter.com/matt_mcguffie) at the [Barrick lab](https://barricklab.org/twiki/bin/view/Lab), University of Texas at Austin, Austin, Texas.
