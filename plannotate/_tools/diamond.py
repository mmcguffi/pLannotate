"""DIAMOND translated-nucleotide-search integration."""

import logging
import shlex
from collections.abc import Mapping
from typing import Any

import pandas as pd

from .._concurrency import parameters_with_threads
from .common import read_table, run_command, temporary_files

COLUMNS = (
    "qseqid qstart qend sseqid pident slen qseq length sstart send qlen evalue btop"
)
# btop is the compact alignment trace; it is alignment text, never a number
TEXT_COLUMNS = ("qseq", "btop")
logger = logging.getLogger(__name__)


def normalize_subject_ids(sequence_ids: pd.Series) -> pd.Series:
    """Unwrap a pipe-delimited protein id to its accession (sp|P12345|NAME -> P12345).

    Only ids that actually carry an accession are unwrapped; ids without a second
    ``|`` field keep their original value instead of becoming NaN. This is the
    single source of truth for how DIAMOND hit ids are keyed -- ``_database_builder``
    reuses it so synthesized descriptions match the reported ``sseqid``.
    """
    ids = sequence_ids.astype(str)
    has_accession = ids.str.contains(r"\|")
    if not has_accession.any():
        return ids
    accessions = ids.str.split("|", n=2).str.get(1)
    return accessions.where(has_accession, ids)


def search(
    sequence: str | Mapping[str, str],
    config: dict[str, Any],
    threads: int = 1,
) -> pd.DataFrame:
    """Search a protein database with DIAMOND blastx.

    ``sequence`` may be one sequence or a ``{query_id: sequence}`` mapping; the
    returned frame carries a ``qseqid`` column identifying each hit's query.
    """
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
        dataframe = read_table(output_path, COLUMNS, TEXT_COLUMNS)

    if not dataframe.empty:
        dataframe["sseqid"] = normalize_subject_ids(dataframe["sseqid"])
    dataframe["sframe"] = (
        (dataframe["qstart"] < dataframe["qend"]).astype(int).replace(0, -1)
    )
    dataframe["slen"] *= 3
    dataframe["length"] = abs(dataframe["qend"] - dataframe["qstart"]) + 1
    # only a covariance-model search reports a consensus structure; the column exists
    # so every source shares one schema
    dataframe["structure"] = ""
    # DIAMOND's pident is a straight count of identical residues
    dataframe["sequence_identity"] = True
    logger.info("DIAMOND found %d candidate hits", len(dataframe))
    return dataframe
