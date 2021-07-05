[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)

<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

Online Annotation
=================

pLannotate is web server for automatically annotating engineered plasmids.

Please visit http://plannotate.barricklab.org/


Local Installation
==================

If you wish you to use pLannotate as a local server, please follow the instructions below (requires [Conda](https://docs.conda.io/en/latest/)).

Download the source code as well as the compressed BLAST (and associated) databases from the [releases](https://github.com/barricklab/pLannotate/releases/tag/v1.0.0) page. Unzip the `BLAST_dbs` folder and place it within the `pLannotate` folder.

On the command line, navigate to the `pLannotate` folder and enter the following commands:

Create the Conda environment:
```
conda env create -f environment.yml
```
Activate the Conda environment:
```
conda activate pLannotate
```
Launch pLannotate:
```
plannotate streamlit [ --blast_db my_blast_db ]
```

After execution of the final command, pLannotate should launch in your default web browser, or you may simply navigate to http://localhost:8501 in your web browser.


Development Installation
========================

First setup a conda environment as above to obtain all the binaries (BLAST etc). Then install the Python code with:

```
python setup.py develop
```

You can now edit the files and run the CLI as above.


About
=====
pLannotate is currently developed by [Matt McGuffie](https://twitter.com/matt_mcguffie) at the [Barrick lab](https://barricklab.org/twiki/bin/view/Lab), University of Texas at Austin, Austin, Texas.
