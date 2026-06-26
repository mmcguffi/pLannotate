"""DIAMOND translated-nucleotide-search integration."""

import logging
import shlex
from typing import Any

import pandas as pd

from .._concurrency import parameters_with_threads
from .common import read_table, run_command, temporary_files

COLUMNS = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
logger = logging.getLogger(__name__)


def search(
    sequence: str,
    config: dict[str, Any],
    threads: int = 1,
) -> pd.DataFrame:
    """Search a protein database with DIAMOND blastx."""
    logger.info("Starting DIAMOND search")
    logger.debug("DIAMOND database=%s threads=%d", config["db_loc"], threads)
    parameters = parameters_with_threads(
        str(config["parameters"]), ("--threads", "-p"), "--threads", threads
    )
    with temporary_files(sequence) as (query_path, output_path):
        command = [
            "diamond",
            "blastx",
            "-d",
            str(config["db_loc"]),
            "-q",
            query_path,
            "-o",
            output_path,
            *shlex.split(parameters),
            "--outfmt",
            "6",
            *COLUMNS.split(),
        ]
        run_command(command, "diamond")
        dataframe = read_table(output_path, COLUMNS)

    if not dataframe.empty:
        sequence_ids = dataframe["sseqid"].astype(str)
        # only unwrap the accession (e.g. sp|P12345|NAME -> P12345) for ids that
        # actually carry it; ids without a second field keep their original value
        # instead of becoming NaN.
        has_accession = sequence_ids.str.contains(r"\|")
        if has_accession.any():
            accessions = sequence_ids.str.split("|", n=2).str.get(1)
            dataframe["sseqid"] = accessions.where(has_accession, sequence_ids)
    dataframe["sframe"] = (
        (dataframe["qstart"] < dataframe["qend"]).astype(int).replace(0, -1)
    )
    dataframe["slen"] *= 3
    dataframe["length"] = abs(dataframe["qend"] - dataframe["qstart"]) + 1
    logger.info("DIAMOND found %d candidate hits", len(dataframe))
    return dataframe
