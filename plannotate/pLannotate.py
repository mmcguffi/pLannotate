import argparse
import sys

import click
import streamlit.cli
import yaml
from bokeh.embed import file_html
from bokeh.resources import CDN, INLINE

from plannotate import resources as rsc
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.streamlit_app import run_streamlit

# possible file structure for better containment
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


@click.group()
@click.version_option(prog_name=__package__)
def main():
    pass


@main.command("streamlit")
@streamlit.cli.configurator_options
@click.option(
    "--yaml_file",
    default=rsc.get_yaml_path(),
    help="path to YAML file.",
    type=click.Path(exists=True),
)
def main_streamlit(yaml_file, **kwargs):
    """Launches pLannotate as an interactive web app."""
    # taken from streamlit.cli.main_hello, @0.78.0
    # streamlit.cli._apply_config_options_from_cli(kwargs)
    # TODO: do this better?
    args = ["--yaml_file", yaml_file]

    if rsc.databases_exist():
        streamlit.cli._main_run(__file__, args)
    else:
        print(
            "Databases not downloaded. Run 'plannotate setupdb' to download databases."
        )


@main.command("yaml")
def main_yaml():
    """Prints YAML file to stdout for custom database modification."""
    with open(rsc.get_yaml_path(), "r") as stream:
        print(
            yaml.dump(
                yaml.load(stream, Loader=yaml.SafeLoader), default_flow_style=False
            )
        )


@main.command("setupdb")
def main_setupdb():
    """Downloads databases; required for use of pLannotate."""

    if rsc.databases_exist():
        print("Databases already downloaded.")
        print()

    else:
        rsc.download_databases()

    print(
        "Run 'plannotate streamlit' or 'plannotate batch {arguments}' to launch pLannotate."
    )
    print(
        "To get a list of available arguments for command line use, run 'plannotate batch --help'."
    )
    print("Please also consider citing: https://doi.org/10.1093/nar/gkab374 :)")


@main.command("batch")
@click.option(
    "--input",
    "-i",
    required=True,
    help="location of a FASTA or GBK file",
)
@click.option(
    "--output",
    "-o",
    default="./",
    help="location of output folder. DEFAULT: current dir",
)
@click.option(
    "--file_name",
    "-f",
    default="",
    help="name of output file (do not add extension). DEFAULT: input file name",
)
@click.option(
    "--suffix",
    "-s",
    default="_pLann",
    help="suffix appended to output files. Use '' for no suffix. DEFAULT: '_pLann'",
)
@click.option(
    "--yaml_file",
    "-y",
    default=rsc.get_yaml_path(),
    help="path to YAML file for custom databases. DEFAULT: builtin",
)
@click.option(
    "--linear",
    "-l",
    is_flag=True,
    help="enables linear DNA annotation",
)
@click.option(
    "--html",
    "-h",
    is_flag=True,
    help="creates an html plasmid map in specified path",
)
@click.option(
    "--htmlfull",
    "-hf",
    is_flag=True,
    help="creates an html plasmid map in specified path, with bokeh baked in",
)
@click.option(
    "--csv",
    "-c",
    is_flag=True,
    help="creates a cvs file in specified path",
)
@click.option(
    "--detailed",
    "-d",
    is_flag=True,
    help="uses modified algorithm for a more-detailed search with more false positives",
)
@click.option(
    "--no_gbk",
    "-x",
    is_flag=True,
    help="supresses GenBank output file",
)
def main_batch(**kwargs):
    """
    Annotates engineered DNA sequences, primarily plasmids. Accepts a FASTA or GenBank file and outputs
    a GenBank file with annotations, as well as an optional interactive plasmid map as an HTLM file.
    """
    if not rsc.databases_exist():
        print(
            "Databases not downloaded. Run 'plannotate setupdb' to download databases."
        )
        sys.exit()

    name, ext = rsc.get_name_ext(kwargs["input"])

    if kwargs["file_name"] == "":
        kwargs["file_name"] = name

    inSeq = rsc.validate_file(kwargs["input"], ext, max_length=float("inf"))

    recordDf = annotate(
        inSeq, kwargs["yaml_file"], kwargs["linear"], kwargs["detailed"]
    )

    if kwargs["no_gbk"] is False:
        gbk = rsc.get_gbk(recordDf, inSeq, kwargs["linear"])
        with open(
            f"{kwargs['output']}/{kwargs['file_name']}{kwargs['suffix']}.gbk", "w"
        ) as handle:
            handle.write(gbk)

    if kwargs["html"] or kwargs["htmlfull"]:
        bokeh_chart = get_bokeh(recordDf, kwargs["linear"])
        bokeh_chart.sizing_mode = "fixed"
        if kwargs["htmlfull"]:
            resource_type = INLINE
        else:
            resource_type = CDN
        html = file_html(
            bokeh_chart, resources=resource_type, title=f"{kwargs['output']}.html"
        )
        with open(
            f"{kwargs['output']}/{kwargs['file_name']}{kwargs['suffix']}.html", "w"
        ) as handle:
            handle.write(html)

    if kwargs["csv"]:
        csv_df = rsc.get_clean_csv_df(recordDf)
        csv_df.to_csv(
            f"{kwargs['output']}/{kwargs['file_name']}{kwargs['suffix']}.csv",
            index=None,
        )


def streamlit_run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml_file")
    args = parser.parse_args()

    run_streamlit(args)


if __name__ == "__main__":
    streamlit_run()
