import logging
from pathlib import Path

import typer
import yaml
from bokeh.embed import file_html
from bokeh.resources import CDN, INLINE

from plannotate import resources as rsc
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.logging_config import get_logger, setup_logging

# Set up logging once at module level
setup_logging(level=logging.INFO)
logger = get_logger(__name__)

app = typer.Typer()


@app.command("yaml")
def main_yaml():
    """Prints YAML file to stdout for custom database modification."""
    with open(rsc.get_yaml_path(), "r") as stream:
        print(
            yaml.dump(
                yaml.load(stream, Loader=yaml.SafeLoader), default_flow_style=False
            )
        )


@app.command("setupdb")
def main_setupdb():
    """Downloads databases; required for use of pLannotate."""

    if rsc.databases_exist():
        logger.info("Databases already downloaded.")

    else:
        rsc.download_databases()

    logger.info("Run 'plannotate batch {arguments}' to launch pLannotate.")
    logger.info(
        "To get a list of available arguments for command line use, run 'plannotate batch --help'."
    )
    logger.info("Please also consider citing: https://doi.org/10.1093/nar/gkab374 :)")


@app.command("batch")
def main_batch(
    input_file: Path = typer.Option(
        ...,
        "--input",
        "-i",
        help="location of a FASTA or GBK file",
        exists=True,
    ),
    output: Path = typer.Option(
        "./",
        "--output",
        "-o",
        help="location of output folder. DEFAULT: current dir",
    ),
    file_name: str = typer.Option(
        "",
        "--file_name",
        "-f",
        help="name of output file (do not add extension). DEFAULT: input file name",
    ),
    suffix: str = typer.Option(
        "_pLann",
        "--suffix",
        "-s",
        help="suffix appended to output files. Use '' for no suffix. DEFAULT: '_pLann'",
    ),
    yaml_file: Path = typer.Option(
        rsc.get_yaml_path(),
        "--yaml_file",
        "-y",
        help="path to YAML file for custom databases. DEFAULT: builtin",
        exists=True,
    ),
    linear: bool = typer.Option(
        False,
        "--linear",
        "-l",
        help="enables linear DNA annotation",
    ),
    html: bool = typer.Option(
        False,
        "--html",
        "-h",
        help="creates an html plasmid map in specified path",
    ),
    htmlfull: bool = typer.Option(
        False,
        "--htmlfull",
        "-hf",
        help="creates an html plasmid map in specified path, with bokeh baked in",
    ),
    csv: bool = typer.Option(
        False,
        "--csv",
        "-c",
        help="creates a cvs file in specified path",
    ),
    detailed: bool = typer.Option(
        False,
        "--detailed",
        "-d",
        help="uses modified algorithm for a more-detailed search with more false positives",
    ),
    no_gbk: bool = typer.Option(
        False,
        "--no_gbk",
        "-x",
        help="supresses GenBank output file",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="enable verbose logging",
    ),
):
    """
    Annotates engineered DNA sequences, primarily plasmids. Accepts a FASTA or GenBank file and outputs
    a GenBank file with annotations, as well as an optional interactive plasmid map as an HTML file.
    """
    # Update logging level if verbose is requested
    if verbose:
        setup_logging(level=logging.DEBUG)

    if not rsc.databases_exist():
        logger.error(
            "Databases not downloaded. Run 'plannotate setupdb' to download databases."
        )
        raise typer.Exit(1)

    name, ext = rsc.get_name_ext(str(input_file))

    if file_name == "":
        file_name = name

    INFINITY = 999_999_999_999
    seq = rsc.validate_file(str(input_file), ext, max_length=INFINITY)

    recordDf = annotate(seq, str(yaml_file), linear, detailed)

    if not no_gbk:
        gbk = rsc.get_gbk(recordDf, seq, linear)
        output_path = output / f"{file_name}{suffix}.gbk"
        with open(output_path, "w") as handle:
            handle.write(gbk)
        logger.info(f"Generated GenBank file: {output_path}")

    if html or htmlfull:
        bokeh_chart = get_bokeh(recordDf, linear)
        bokeh_chart.sizing_mode = "fixed"
        if htmlfull:
            resource_type = INLINE
        else:
            resource_type = CDN
        html_content = file_html(
            bokeh_chart, resources=resource_type, title=f"{output}.html"
        )
        html_path = output / f"{file_name}{suffix}.html"
        with open(html_path, "w") as handle:
            handle.write(html_content)
        logger.info(f"Generated HTML file: {html_path}")

    if csv:
        csv_df = rsc.get_clean_csv_df(recordDf)
        csv_path = output / f"{file_name}{suffix}.csv"
        csv_df.to_csv(csv_path, index=None)
        logger.info(f"Generated CSV file: {csv_path}")


def main():
    app()


if __name__ == "__main__":
    main()
