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

If you wish you to use pLannotate as a local server or to use as a command line tool, please follow the instructions below (requires [Conda](https://docs.conda.io/en/latest/)).

Download the source code as well as the compressed BLAST (and associated) databases from the [releases](https://github.com/barricklab/pLannotate/releases/tag/v1.1.0) page.

On the command line, navigate to the directory where you have placed `pLannotate` folder and compressed database folder.

Enter the following commands:
```
tar -zxf BLAST_dbs.tar.gz -C ./pLannotate/plannotate/data && rm BLAST_dbs.tar.gz
cd pLannotate
conda env create -f environment.yml
conda activate plannotate
python setup.py install
```

To launch pLannotate as a local web server:
```
plannotate streamlit
```

After execution of the final command, pLannotate should launch in your default web browser, or you may simply navigate to http://localhost:8501 in your web browser.

Command Line Interface (batch mode)
===================================

To annotate FASTA files and generate the interactive plasmid maps on the command line,
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
`plannotate batch -i ./plannotate/data/fastas/pUC19.fa --html --output ~/Desktop/ --file_name pLasmid`

About
=====
pLannotate is currently developed by [Matt McGuffie](https://twitter.com/matt_mcguffie) at the [Barrick lab](https://barricklab.org/twiki/bin/view/Lab), University of Texas at Austin, Austin, Texas.
