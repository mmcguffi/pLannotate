"""Database search functions and orchestration.

This module contains functions for running searches against different types of
databases (BLAST, DIAMOND, Infernal) and orchestrating searches across multiple databases.

Author: Matt McGuffie
"""

import shlex
import subprocess
from functools import lru_cache
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, List, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .infernal import parse_infernal
from .logging_config import get_logger

logger = get_logger(__name__)

# Type definitions for database configuration
DatabaseConfig = Dict[
    str,
    Union[str, int, List[str], Dict[str, Union[str, bool, List[str]]]],
]


def blast(seq: str, db: DatabaseConfig, mode: str = "blastn") -> pd.DataFrame:
    """Run BLASTn search against a nucleotide database."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"{mode} -task blastn-short -query {query.name} -out {tmp.name} "
        f'-db {db_loc} {parameters} -outfmt "6 {FLAGS}"'
    )

    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )

    with open(tmp.name, "r") as file_handle:
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=FLAGS.split())
    # Convert numeric columns, handling non-numeric values gracefully
    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col])
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass

    return inDf


def diamond(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run DIAMOND blastx search against a protein database."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"diamond blastx -d {db_loc} -q {query.name} -o {tmp.name} "
        f"{parameters} --outfmt 6 {FLAGS}"
    )
    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )

    with open(tmp.name, "r") as file_handle:
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=FLAGS.split())
    # Convert numeric columns, handling non-numeric values gracefully
    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col])
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass

    # DIAMOND-specific post-processing
    # Only apply pipe-splitting for databases that use pipe-separated identifiers (like SwissProt)
    if "sseqid" in inDf.columns and not inDf["sseqid"].empty:
        # Check if any sequence IDs contain pipes before attempting to split
        sseqid_str = inDf["sseqid"].astype(str)
        if sseqid_str.str.contains("\\|", regex=True).any():
            inDf["sseqid"] = sseqid_str.str.split("|", n=2).str.get(1)
    inDf["sframe"] = (inDf["qstart"] < inDf["qend"]).astype(int).replace(0, -1)
    inDf["slen"] = inDf["slen"] * 3
    inDf["length"] = abs(inDf["qend"] - inDf["qstart"]) + 1

    return inDf


def infernal(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run Infernal cmscan search against RNA covariance models."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
    cmd = f"cmscan {FLAGS} {parameters} --tblout {tmp.name} --clanin {db_loc} {query.name}"
    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )
    inDf = parse_infernal(tmp.name)

    inDf["qlen"] = len(seq)

    # manually gets DNA sequence from seq(x2)
    if not inDf.empty:
        inDf["qseq"] = inDf.apply(
            lambda x: (seq)[x["qstart"] : x["qend"] + 1].upper(), axis=1
        )

    tmp.close()
    query.close()

    return inDf


def _run_database_search(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run search against a database using the appropriate method."""
    task = db["method"]

    if task == "blastn":
        return blast(seq, db)
    elif task == "diamond":
        return diamond(seq, db)
    elif task == "infernal":
        return infernal(seq, db)
    else:
        raise ValueError(f"Unknown search method: {task}")


@lru_cache(maxsize=5)
def _search_all_databases_cached(
    query_sequence: str,
    is_linear: bool,
    yaml_file_str: str,
) -> pd.DataFrame:
    """Cached version of database search to avoid re-running expensive searches."""
    yaml_file = Path(yaml_file_str)
    logger.info("Starting annotation...")

    databases = rsc.get_yaml(yaml_file)
    all_database_hits = _search_individual_databases(
        query_sequence, databases, yaml_file, is_linear
    )

    if not all_database_hits:
        return pd.DataFrame()

    # Import here to avoid circular import
    from .filter_annotations import _combine_and_filter_results

    # Combine results using filtering module
    final_results = _combine_and_filter_results(all_database_hits)

    logger.info("Database search complete!")
    return final_results


def search_all_databases(
    query_sequence: str,
    is_linear: bool,
    yaml_file: Path,
) -> pd.DataFrame:
    """Search query sequence against all configured databases.
    
    This function includes automatic caching for identical inputs to improve performance
    when processing the same DNA sequence multiple times. The cache is applied at the
    database search level so that both detailed and non-detailed annotation modes can
    reuse the same raw search results.
    """
    # Convert inputs to hashable types for caching
    yaml_file_str = str(yaml_file)
    
    # Use cached version for computational efficiency
    return _search_all_databases_cached(query_sequence, is_linear, yaml_file_str)


def _search_individual_databases(
    query_sequence: str,
    databases: dict,
    yaml_file: Path,
    is_linear: bool,
) -> List[pd.DataFrame]:
    """Search each database individually and return processed hits."""

    database_results = []

    for i, (database_name, database_config) in enumerate(databases.items(), 1):
        logger.info(f"Processing database {i}/{len(databases)}: {database_name}")

        # Search this database
        hits = _run_database_search(seq=query_sequence, db=database_config)

        if hits.empty:
            continue

        # Process hits for this database
        processed_hits = _process_database_hits(
            hits, database_name, database_config, yaml_file, is_linear
        )

        if not processed_hits.empty:
            database_results.append(processed_hits)

    return database_results


def _process_database_hits(
    hits: pd.DataFrame,
    database_name: str,
    database_config: dict,
    yaml_file: Path,
    is_linear: bool,
) -> pd.DataFrame:
    """Process hits from a single database."""
    # Import here to avoid circular import
    from .filter_annotations import (
        load_feature_details,
        calculate_hit_scores,
        _merge_feature_details,
        _apply_database_priority,
    )

    # Add database metadata
    hits["db"] = database_name
    hits["sseqid"] = hits["sseqid"].astype(str)

    # Load feature details
    feature_details = load_feature_details(hits, yaml_file)

    # Merge with feature descriptions
    hits = _merge_feature_details(hits, feature_details)

    # Filter out primer binding sites
    hits = hits.loc[hits["Type"] != "primer_bind"]

    # Apply database priority
    hits = _apply_database_priority(hits, database_config)

    # Calculate scores
    hits = calculate_hit_scores(hits, is_linear)

    return hits
