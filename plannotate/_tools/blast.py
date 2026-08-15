"""BLAST nucleotide-search integration."""

import logging
import shlex
from collections.abc import Mapping
from typing import Any

import pandas as pd

from .._concurrency import parameters_with_threads
from .common import read_table, run_command, temporary_files

COLUMNS = (
    "qseqid qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue "
    "btop"
)
# btop is BLAST's compact alignment trace; it is alignment text, never a number
TEXT_COLUMNS = ("qseq", "btop")
logger = logging.getLogger(__name__)


def search(
    sequence: str | Mapping[str, str],
    config: dict[str, Any],
    threads: int = 1,
    executable: str = "blastn",
) -> pd.DataFrame:
    """Search a nucleotide database with BLAST.

    ``sequence`` may be one sequence or a ``{query_id: sequence}`` mapping; the
    returned frame carries a ``qseqid`` column identifying each hit's query.
    """
    logger.info("Starting BLAST search")
    logger.debug("BLAST database=%s threads=%d", config["db_loc"], threads)
    parameters = parameters_with_threads(
        str(config["parameters"]), ("-num_threads",), "-num_threads", threads
    )
    with temporary_files(sequence) as (query_path, output_path):
        command = [
            executable,
            "-task",
            "blastn-short",
            "-query",
            query_path,
            "-out",
            output_path,
            "-db",
            str(config["db_loc"]),
            *shlex.split(parameters),
            "-outfmt",
            f"6 {COLUMNS}",
        ]
        run_command(command, executable)
        dataframe = read_table(output_path, COLUMNS, TEXT_COLUMNS)
    # only a covariance-model search reports a consensus structure; the column exists
    # so every source shares one schema
    dataframe["structure"] = ""
    # BLAST's pident is a straight count of matching bases
    dataframe["sequence_identity"] = True
    logger.info("BLAST found %d candidate hits", len(dataframe))
    return dataframe
