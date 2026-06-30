"""Infernal covariance-model-search integration."""

import logging
import shlex
from collections.abc import Mapping
from itertools import accumulate
from pathlib import Path
from typing import Any

import pandas as pd

from .._concurrency import parameters_with_threads
from .common import normalize_queries, run_command, temporary_files

REQUIRED_COLUMNS = [
    "#idx",
    "target name",
    "accession",
    "query name",
    "clan name",
    "seq from",
    "seq to",
    "mdl from",
    "mdl to",
    "strand",
    "score",
    "E-value",
    "description of target",
]
logger = logging.getLogger(__name__)


def search(
    sequence: str | Mapping[str, str],
    config: dict[str, Any],
    threads: int = 1,
) -> pd.DataFrame:
    """Search a covariance-model database with Infernal cmscan.

    ``sequence`` may be one sequence or a ``{query_id: sequence}`` mapping; the
    returned frame carries a ``qseqid`` column identifying each hit's query.
    """
    logger.info("Starting Infernal search")
    logger.debug(
        "Infernal database=%s clan=%s threads=%d",
        config.get("db_loc"),
        config.get("clanin_loc"),
        threads,
    )
    parameters = parameters_with_threads(
        str(config["parameters"]), ("--cpu",), "--cpu", threads
    )
    cm_path, clan_path = _configured_database_paths(config)
    with temporary_files(sequence) as (query_path, output_path):
        command = [
            "cmscan",
            "--cut_ga",
            "--rfam",
            "--noali",
            "--fmt",
            "2",
            *shlex.split(parameters),
            "--tblout",
            output_path,
            "--clanin",
            clan_path,
            cm_path,
            query_path,
        ]
        run_command(command, "cmscan")
        dataframe = parse_output(output_path)

    # qlen/qseq are derived from the matching query, looked up by qseqid, so a
    # batched multi-FASTA run resolves each hit against its own sequence. The qlen
    # set here is the searched (possibly doubled) length; the collector overrides it
    # with the true plasmid length, matching the blast/diamond path.
    queries = normalize_queries(sequence)
    dataframe["qlen"] = [len(queries[query_id]) for query_id in dataframe["qseqid"]]
    if not dataframe.empty:
        dataframe["qseq"] = dataframe.apply(
            lambda row: queries[row["qseqid"]][row["qstart"] - 1 : row["qend"]].upper(),
            axis=1,
        )
    logger.info("Infernal found %d candidate hits", len(dataframe))
    return dataframe


def database_paths(database_name: str, directory: Path) -> dict[str, str]:
    """Build the CM and clan paths for an Infernal database."""
    return {
        "db_loc": str(directory / f"{database_name}.cm"),
        "clanin_loc": str(directory / f"{database_name}.clanin"),
    }


def _configured_database_paths(config: dict[str, Any]) -> tuple[str, str]:
    """Read the configured covariance-model and clan paths."""
    if "db_loc" not in config or "clanin_loc" not in config:
        raise ValueError("Infernal configuration requires CM and clan database paths")
    return str(config["db_loc"]), str(config["clanin_loc"])


def parse_output(path: str | Path) -> pd.DataFrame:
    """Parse Infernal ``--tblout --fmt 2`` output into candidate hits."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("Infernal output is missing its column headers")

    widths = [len(field) + 1 for field in lines[1].split()]
    ends = list(accumulate(widths))
    ends[-1] += 100
    starts = [0, *ends[:-1]]
    positions = list(zip(starts, ends, strict=True))
    names = [lines[0][start:end].strip() for start, end in positions]

    try:
        dataframe = pd.read_fwf(path, comment="#", colspecs=positions, header=None)
        dataframe.columns = names
    except pd.errors.EmptyDataError:
        dataframe = pd.DataFrame(columns=names)

    dataframe = dataframe[REQUIRED_COLUMNS]
    dataframe = dataframe.loc[:, ~dataframe.columns.duplicated()]
    dataframe = dataframe.rename(
        columns={
            "target name": "name",
            "query name": "qseqid",
            "seq from": "qstart",
            "seq to": "qend",
            "mdl from": "sstart",
            "mdl to": "send",
            "E-value": "evalue",
            "strand": "sframe",
            "description of target": "blurb",
        }
    )

    dataframe["accession"] = dataframe["accession"].str.replace("-", " ")
    dataframe["clan name"] = dataframe["clan name"].str.replace("-", " ")
    dataframe["name"] = dataframe["name"].str.replace("_", " ")
    # Identify each hit by its stable Rfam accession (e.g. RF00162) instead of
    # cmscan's per-run hit ordinal, so the id is meaningful and reproducible across
    # runs. Fall back to the model name on the rare model without an accession.
    clean_accession = dataframe["accession"].str.strip()
    dataframe["sseqid"] = clean_accession.where(
        clean_accession != "", dataframe["name"]
    )
    dataframe["blurb"] = (
        "Accession: " + dataframe["accession"] + " - " + dataframe["blurb"]
    )
    dataframe["type"] = "ncRNA"
    dataframe["qseq"] = ""

    coordinates = dataframe[["qstart", "qend"]].apply(pd.to_numeric)
    dataframe["qstart"] = coordinates.min(axis=1).astype("int64")
    dataframe["qend"] = coordinates.max(axis=1).astype("int64")
    dataframe["sframe"] = dataframe["sframe"].map({"-": -1, "+": 1})
    dataframe["length"] = abs(dataframe["qend"] - dataframe["qstart"]) + 1
    dataframe["slen"] = abs(dataframe["send"] - dataframe["sstart"]) + 1
    dataframe["pident"] = 100
    return dataframe.drop(columns=["#idx", "accession", "clan name"])
