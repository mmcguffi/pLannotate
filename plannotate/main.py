"""
Main entry point for pLannotate, a plasmid annotation tool.

This module provides command-line interfaces for:
- Printing a YAML file for custom database modification.
- Setting up the database by downloading required files.
- Running batch annotations on plasmid sequences from FASTA or GenBank files.

Author: Matt McGuffie
"""

import json
import logging
import sys
from pathlib import Path

import typer
import yaml

from . import _package_data, validation
from .models import Construct

logger = logging.getLogger(__name__)
app = typer.Typer()
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def _configure_logging(level: int = logging.INFO) -> None:
    package_logger = logging.getLogger("plannotate")
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    package_logger.handlers.clear()
    package_logger.addHandler(handler)
    package_logger.setLevel(level)
    package_logger.propagate = False


@app.command("yaml")
def main_yaml():
    """Print the search configuration for custom database modification."""
    with _package_data.get_yaml_path().open() as stream:
        typer.echo(yaml.safe_dump(yaml.safe_load(stream), sort_keys=False))


@app.command("databases")
def main_databases():
    """Print versions and checksums for the installed database bundle."""
    try:
        manifest = _package_data.get_database_manifest()
    except FileNotFoundError as exc:
        logger.error(str(exc))
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(manifest, indent=2, sort_keys=True))


@app.command("setupdb")
def main_setupdb(
    force: bool = typer.Option(
        False,
        "--force",
        help="replace an existing database installation",
    ),
):
    """Downloads databases; required for use of pLannotate."""

    if _package_data.databases_exist() and not force:
        logger.info("Databases already downloaded.")

    else:
        _package_data.download_databases()

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
        "--file-name",
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
        _package_data.get_yaml_path(),
        "--yaml-file",
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
        help="creates a CSV file in specified path",
    ),
    detailed: bool = typer.Option(
        False,
        "--detailed",
        "-d",
        help="uses modified algorithm for a more-detailed search with more false positives",
    ),
    cores: int = typer.Option(
        1,
        "--cores",
        "-j",
        min=1,
        help="maximum annotation sources to run in parallel",
    ),
    no_gbk: bool = typer.Option(
        False,
        "--no-gbk",
        "--no_gbk",
        "-x",
        help="suppresses GenBank output file",
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
    if verbose:
        _configure_logging(logging.DEBUG)

    if not _package_data.databases_exist():
        logger.error(
            "Databases not downloaded. Run 'plannotate setupdb' to download databases."
        )
        raise typer.Exit(1)

    name, ext = validation.get_name_ext(str(input_file))

    if not file_name:
        file_name = name

    seqrecord = validation.validate_file(
        input_file,
        ext,
        max_length=None,
    )

    construct = Construct(
        seq=str(seqrecord.seq),
        linear=linear,
        detailed=detailed,
        db_options=yaml_file,
        prior_annotations=seqrecord if ext in validation.VALID_GENBANK_EXTS else None,
        cores=cores,
    )

    output.mkdir(parents=True, exist_ok=True)

    if not no_gbk:
        gbk = construct.to_genbank()
        output_path = output / f"{file_name}{suffix}.gbk"
        output_path.write_text(gbk)
        logger.info("Generated GenBank file: %s", output_path)

    if html or htmlfull:
        html_content = construct.to_html(htmlfull=htmlfull)
        html_path = output / f"{file_name}{suffix}.html"
        html_path.write_text(html_content)
        logger.info("Generated HTML file: %s", html_path)

    if csv:
        csv_df = construct.to_csv()
        csv_path = output / f"{file_name}{suffix}.csv"
        csv_df.to_csv(csv_path, index=False)
        logger.info("Generated CSV file: %s", csv_path)


def main() -> None:
    _configure_logging()
    app()


if __name__ == "__main__":
    main()
