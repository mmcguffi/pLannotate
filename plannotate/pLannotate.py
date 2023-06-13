import argparse
from dataclasses import dataclass
import sys

import streamlit.cli
from bokeh.embed import file_html
from bokeh.resources import CDN
import yaml

import plannotate.resources as rsc
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.streamlit_app import run_streamlit

import defopt
from pathlib import Path
from typing import Callable
from typing import List
from typing import Optional

#possible file structure for better containment
# plasmid = {
#     'fileloc': '',
#     'name': '',
#     'ext': '',
#     'blast_db': './BLAST_dbs/',
#     'linear': False,
#     'seq': '',
#     'raw_hits': pd.DataFrame(),
#     'hits': pd.DataFrame(),
#     'hits_detailed' : pd.DataFrame()
# }

# NOTE: streamline really wants us to use their entry point and give
#       them a script to run. Here we follow the hello world example
#       to bootstrap running of this file as a script (the streamlit_run
#       function). Unfortunately we have to buy in to using click as
#       the command-line parse in front of streamlit, but then also
#       use standard argparse to parse the final options in our script.


def streamlit(*, yaml_file: Path = Path(rsc.get_yaml_path())) -> None:
    """Launches pLannotate as an interactive web app.

    Args:
        yaml_file: the Path to the YAML file.
    """
    # taken from streamlit.cli.main_hello, @0.78.0
    #streamlit.cli._apply_config_options_from_cli(kwargs)
    # TODO: do this better?
    args = ['--yaml_file', yaml_file]
    
    if rsc.databases_exist():
        streamlit.cli._main_run(__file__, args)
    else:
        print("Databases not downloaded. Run 'plannotate setupdb' to download databases.")


def yaml() -> None:
    """Prints YAML file to stdout for custom database modification."""
    with open (rsc.get_yaml_path(), 'r') as stream:
        print(yaml.dump(yaml.load(stream, Loader = yaml.SafeLoader), default_flow_style=False))

        
def setupdb() -> None:
    """Downloads databases; required for use of pLannotate."""
        
    if rsc.databases_exist():
        print("Databases already downloaded.")
        print()
        
    else:
       rsc.download_databases()

    print("Run 'plannotate streamlit' or 'plannotate batch {arguments}' to launch pLannotate.")
    print("To get a list of available arguments for command line use, run 'plannotate batch --help'.")
    print("Please also consider citing: https://doi.org/10.1093/nar/gkab374 :)")


def batch(*,
          input: Path,
          output: Path = Path("./"),
          file_name: Optional[str] = None,
          suffix: str = "_pLann",
          yaml_file: Path = Path(rsc.get_yaml_path()),
          linear: bool = True,
          html: bool = True,
          csv: bool = True,
          detailed: bool = True,
          gbk: bool = False) -> None:
    """
    Annotates engineered DNA sequences, primarily plasmids. Accepts a FASTA or GenBank file and outputs
    a GenBank file with annotations, as well as an optional interactive plasmid map as an HTLM file.

    Args:
        input: path to the input FASTA or GBK file
        output: path to the output folder
        file_name: name of output file (do not add extension). If not provided, uses the input file name.
        suffix: suffix appended to output files. Use '' for no suffix.
        yaml: path to YAML file for custom databases. Defaults to the built-in databases. 
        linear: enables linear DNA annotation
        html: creates an html plasmid map in specified path
        csv: creates a cvs file in specified path
        detailed: uses modified algorithm for a more-detailed search with more false positives
        gkb: output a GenBank file
    """
    if not rsc.databases_exist():
        print("Databases not downloaded. Run 'plannotate setupdb' to download databases.")
        sys.exit()

    name, ext = rsc.get_name_ext(input)

    if file_name is None or file_name == "":
        file_name = name

    inSeq = rsc.validate_file(input, ext, max_length = float("inf"))

    recordDf = annotate(inSeq, yaml_file, linear, detailed)

    if gbk:
        gbk_rec = rsc.get_gbk(recordDf, inSeq, linear)
        with open(f"{output}/{file_name}{suffix}.gbk", "w") as handle:
            handle.write(gbk_rec)

    if html:
        bokeh_chart = get_bokeh(recordDf, linear)
        bokeh_chart.sizing_mode = "fixed"
        html = file_html(bokeh_chart, resources = CDN, title = f"{output}.html")
        with open(f"{output}/{file_name}{suffix}.html", "w") as handle:
            handle.write(html)

    if csv:
        csv_df = rsc.get_clean_csv_df(recordDf)
        csv_df.to_csv(f"{output}/{file_name}{suffix}.csv", index = None)


def main(argv: List[str] = sys.argv[1:]) -> None:
    tools: List[Callable] = [
        batch,
        setupdb,
        streamlit,
        yaml,
    ]
    if len(argv) != 0 and all(arg not in argv for arg in ["-h", "--help"]):
        print("Running command: client-tools " + " ".join(argv))
    try:
        defopt.run(funcs=tools, argv=argv)
        print("Completed successfully.")
    except Exception as e:
        sys.stderr.write("Failed on command: " + " ".join(argv) + "\n")
        raise e


if __name__ == '__main__':
    main()
