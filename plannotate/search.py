"""Database search functions and orchestration.

This module contains functions for running searches against different types of
databases (BLAST, DIAMOND, Infernal) and orchestrating searches across multiple databases.

Author: Matt McGuffie
"""

import json
import os
import shlex
import shutil
import subprocess
import sys
from contextlib import contextmanager
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Dict, Iterator, List, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from . import sqlite_utils
from .infernal import parse_infernal
from .logging_config import get_logger

logger = get_logger(__name__)

# Type definitions for database configuration
DatabaseConfig = Dict[
    str,
    Union[str, int, List[str], Dict[str, Union[str, bool, List[str]]]],
]

# Constants for search

# Column definitions
BLAST_COLS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "score",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
    "db",
]

EXTRA_COLS = [
    "name",
    "blurb",
    "type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "fragment",
]

DEV_COLS = ["wiggle", "wstart", "wend", "kind", "qstart_dup", "qend_dup"]

DF_COLS = BLAST_COLS + EXTRA_COLS + DEV_COLS


def _run_external_command(command: str, tool: str) -> None:
    """Run a search tool and report failures with its diagnostic output."""
    try:
        result = subprocess.run(
            shlex.split(command),
            shell=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{tool} is not installed or is not available on PATH"
        ) from exc
    if result.returncode != 0:
        diagnostic = result.stderr.strip() or result.stdout.strip() or "no output"
        raise RuntimeError(
            f"{tool} failed with exit code {result.returncode}: {diagnostic}"
        )


@contextmanager
def _search_temp_files(seq: str) -> Iterator[tuple[str, str]]:
    """Provide a FASTA query and output path that close on every exit path."""
    with NamedTemporaryFile() as query, NamedTemporaryFile() as output:
        SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")
        yield query.name, output.name


def _parameters_with_threads(
    parameters: str,
    thread_options: tuple[str, ...],
    thread_option: str,
    threads: int,
) -> str:
    """Replace user-supplied thread flags with the scheduler allocation."""
    if threads < 1:
        raise ValueError("threads must be at least 1")

    tokens = shlex.split(parameters)
    filtered = []
    skip_value = False
    for token in tokens:
        if skip_value:
            skip_value = False
            continue
        if token in thread_options:
            skip_value = True
            continue
        if any(token.startswith(f"{option}=") for option in thread_options):
            continue
        filtered.append(token)
    filtered.extend((thread_option, str(threads)))
    return shlex.join(filtered)


def blast(
    seq: str,
    db: DatabaseConfig,
    mode: str = "blastn",
    threads: int = 1,
) -> pd.DataFrame:
    """Run BLASTn search against a nucleotide database."""
    parameters = _parameters_with_threads(
        str(db["parameters"]), ("-num_threads",), "-num_threads", threads
    )
    db_loc = db["db_loc"]
    flags = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
    with _search_temp_files(seq) as (query_path, output_path):
        cmd = (
            f"{mode} -task blastn-short -query {query_path} -out {output_path} "
            f'-db {db_loc} {parameters} -outfmt "6 {flags}"'
        )
        _run_external_command(cmd, mode)
        with open(output_path) as file_handle:
            align = file_handle.readlines()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=flags.split())
    # Convert numeric columns, handling non-numeric values gracefully
    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col])
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass

    return inDf


def diamond(seq: str, db: DatabaseConfig, threads: int = 1) -> pd.DataFrame:
    """Run DIAMOND blastx search against a protein database."""
    parameters = _parameters_with_threads(
        str(db["parameters"]), ("--threads", "-p"), "--threads", threads
    )
    db_loc = db["db_loc"]
    flags = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
    with _search_temp_files(seq) as (query_path, output_path):
        cmd = (
            f"diamond blastx -d {db_loc} -q {query_path} -o {output_path} "
            f"{parameters} --outfmt 6 {flags}"
        )
        _run_external_command(cmd, "diamond")
        with open(output_path) as file_handle:
            align = file_handle.readlines()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=flags.split())
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


def infernal(seq: str, db: DatabaseConfig, threads: int = 1) -> pd.DataFrame:
    """Run Infernal cmscan search against RNA covariance models."""
    parameters = _parameters_with_threads(
        str(db["parameters"]), ("--cpu",), "--cpu", threads
    )
    db_loc = db["db_loc"]
    flags = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
    with _search_temp_files(seq) as (query_path, output_path):
        cmd = (
            f"cmscan {flags} {parameters} --tblout {output_path} "
            f"--clanin {db_loc} {query_path}"
        )
        _run_external_command(cmd, "cmscan")
        inDf = parse_infernal(output_path)

    inDf["qlen"] = len(seq)

    # manually gets DNA sequence from seq(x2)
    if not inDf.empty:
        inDf["qseq"] = inDf.apply(
            lambda x: seq[x["qstart"] - 1 : x["qend"]].upper(), axis=1
        )

    return inDf


def _run_database_search(
    seq: str,
    db: DatabaseConfig,
    threads: int = 1,
) -> pd.DataFrame:
    """Run search against a database using the appropriate method."""
    task = db["method"]

    if task == "blastn":
        return blast(seq, db, threads=threads)
    elif task == "diamond":
        return diamond(seq, db, threads=threads)
    elif task == "infernal":
        return infernal(seq, db, threads=threads)
    else:
        raise ValueError(f"Unknown search method: {task}")


def _search_single_database(
    query: str,
    database_name: str,
    database_config: DatabaseConfig,
    yaml_file: Path,
    is_linear: bool,
    threads: int = 1,
) -> pd.DataFrame:
    """Search a single database and return processed hits."""
    if not is_linear:
        logger.debug("doubling query for circular sequence")
        query = query + query

    hits = _run_database_search(seq=query, db=database_config, threads=threads)

    if hits.empty:
        return pd.DataFrame()

    # Process hits for this database
    processed_hits = _process_database_hits(
        hits,
        database_name,
        database_config,
        yaml_file,
    )

    return processed_hits


def _process_database_hits(
    hits: pd.DataFrame,
    database_name: str,
    database_config: dict,
    yaml_file: Path,
) -> pd.DataFrame:
    """Process hits from a single database."""

    # Add database metadata
    hits["db"] = database_name
    hits["sseqid"] = hits["sseqid"].astype(str)

    # Load feature details
    feature_details = _load_feature_details(hits, yaml_file)

    # Merge with feature descriptions
    hits = _merge_feature_details(hits, feature_details)

    # Filter out primer binding sites (only if type column exists)
    if "type" in hits.columns:
        hits = hits.loc[hits["type"] != "primer_bind"]

    # Apply database priority
    hits = _apply_database_priority(hits, database_config)

    # Return raw hits with metadata - scoring will be done during filtering
    return hits


def _load_feature_details(hits: pd.DataFrame, yaml_file: Path) -> pd.DataFrame:
    """Load detailed feature information from database files."""
    database_name = _validate_single_database(hits)
    databases = rsc.get_yaml(yaml_file)
    database_config = databases[database_name]

    # Clean sequence IDs
    hits = _clean_sequence_ids(hits)

    # Get sequence IDs for lookup
    sequence_ids = _extract_sequence_ids(hits, database_name)

    # Load feature descriptions
    feature_details = _load_feature_descriptions(
        database_config, database_name, sequence_ids, hits
    )

    # Apply database-specific processing
    if database_name == "swissprot":
        feature_details = _process_swissprot_details(feature_details)

    # Apply default type if specified
    feature_details = _apply_default_type(feature_details, database_config)

    return feature_details


def _validate_single_database(hits: pd.DataFrame) -> str:
    """Ensure all hits come from the same database."""
    unique_databases = hits["db"].unique()
    if len(unique_databases) != 1:
        raise ValueError("All hits must be from the same database")
    return unique_databases[0]


def _clean_sequence_ids(hits: pd.DataFrame) -> pd.DataFrame:
    """Clean problematic sequence ID formats."""
    # Fix "pdb|3xHA|" format to extract middle part
    problematic_pattern = r"pdb\|(.*)\|"
    hits["sseqid"] = hits["sseqid"].str.replace(problematic_pattern, r"\1", regex=True)
    return hits


def _extract_sequence_ids(hits: pd.DataFrame, database_name: str) -> List[str]:
    """Extract and filter sequence IDs for database lookup."""
    sequence_ids = hits.loc[hits["db"] == database_name, "sseqid"].tolist()
    return [sid for sid in sequence_ids if sid]  # Remove empty values


def _load_feature_descriptions(
    database_config: dict,
    database_name: str,
    sequence_ids: List[str],
    hits: pd.DataFrame,
) -> pd.DataFrame:
    """Load feature descriptions from SQLite database."""

    db_details = database_config["details"]

    if db_details["location"] == "None":
        # Data already in dataframe (e.g., Rfam)
        return hits[["sseqid", "name", "type", "blurb"]].copy()

    # Use SQLite database for all other cases
    data_dir = rsc.get_data_directory()
    sequence_ids_set = set(sequence_ids)
    return sqlite_utils.load_descriptions_from_sqlite(
        database_name, data_dir, sequence_ids_set, database_config
    )


def _process_swissprot_details(feature_details: pd.DataFrame) -> pd.DataFrame:
    """Extract SwissProt protein existence levels for priority modification."""
    # Calculate priority modifications based on existence levels
    priority_modifications = [
        _extract_existence_level_priority(description)
        for description in feature_details["blurb"]
    ]

    feature_details["priority_mod"] = priority_modifications
    return feature_details


def _extract_existence_level_priority(description: str) -> int:
    """Extract priority modification from SwissProt existence level."""
    existence_marker = "existence level"
    marker_position = description.find(existence_marker)

    if marker_position == -1:
        return 0  # No existence level found

    # Extract the digit immediately after "existence level "
    level_start = marker_position + len(existence_marker) + 1
    level_end = level_start + 1

    try:
        existence_level = int(description[level_start:level_end])
        return existence_level - 1  # Convert to 0-based priority modification
    except (ValueError, IndexError):
        return 0


def _apply_default_type(
    feature_details: pd.DataFrame,
    database_config: dict,
) -> pd.DataFrame:
    """Apply default feature type if specified in database configuration."""
    default_type = database_config["details"].get("default_type")

    if default_type and default_type != "None":
        feature_details["type"] = default_type

    return feature_details


def _combine_results(database_results: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine results from all databases."""
    # Combine all results
    combined_results = pd.concat(database_results, ignore_index=True)

    # Note: Sorting by scores will happen in filter_annotations after scoring
    return combined_results


def _merge_feature_details(
    hits: pd.DataFrame,
    feature_details: pd.DataFrame,
) -> pd.DataFrame:
    """Merge hits with feature details, handling duplicate columns."""
    merged = hits.merge(
        feature_details, on="sseqid", how="left", suffixes=("_original", "")
    )

    # Remove duplicate columns from original data
    duplicate_columns = [col for col in merged.columns if col.endswith("_original")]
    merged = merged.drop(columns=duplicate_columns)

    return merged


def _apply_database_priority(hits: pd.DataFrame, database_config: dict) -> pd.DataFrame:
    """Apply database priority and any priority modifications."""
    hits["priority"] = database_config["priority"]

    # Apply priority modifications if available (e.g., from SwissProt existence levels)
    if "priority_mod" in hits.columns:
        hits["priority"] += hits["priority_mod"]
        hits = hits.drop("priority_mod", axis=1)

    return hits


@lru_cache(maxsize=5)
def _search_all_databases_cached(
    query_sequence: str,
    is_linear: bool,
    yaml_file_str: str,
    cores: int,
) -> pd.DataFrame:
    """Run the database-search workflow, caching identical invocations."""
    yaml_file = Path(yaml_file_str)
    logger.info("Starting annotation...")

    databases = rsc.get_yaml(yaml_file)
    if not databases:
        return pd.DataFrame()

    database_results = _run_annotation_workflow(
        query_sequence=query_sequence,
        is_linear=is_linear,
        yaml_file=yaml_file,
        cores=cores,
    )

    if not database_results:
        return pd.DataFrame()

    # Combine results from all databases
    final_results = _combine_results(database_results)

    logger.info("Database search complete!")
    return final_results


def _run_annotation_workflow(
    query_sequence: str,
    is_linear: bool,
    yaml_file: Path,
    cores: int,
) -> List[pd.DataFrame]:
    """Run one Snakemake job per configured annotation database."""
    if cores < 1:
        raise ValueError("cores must be at least 1")

    snakemake = shutil.which("snakemake")
    if snakemake is None:
        raise RuntimeError(
            "Snakemake is required for annotation; reinstall pLannotate or install "
            "snakemake>=7.18."
        )

    databases = rsc.get_yaml(yaml_file)
    database_count = len(databases)
    thread_allocations = _allocate_search_threads(databases, cores)

    workflow = Path(str(files(rsc.PACKAGE) / "workflows/annotation.smk"))
    with TemporaryDirectory(prefix="plannotate-annotation-") as temp_dir:
        work_dir = Path(temp_dir)
        query_file = work_dir / "query.txt"
        config_file = work_dir / "config.json"
        query_file.write_text(query_sequence)
        config_file.write_text(
            json.dumps(
                {
                    "database_count": database_count,
                    "linear": is_linear,
                    "python_executable": sys.executable,
                    "query_file": str(query_file),
                    "thread_allocations": thread_allocations,
                    "yaml_file": str(yaml_file.resolve()),
                }
            )
        )

        command = [
            snakemake,
            "--snakefile",
            str(workflow),
            "--directory",
            str(work_dir),
            "--configfile",
            str(config_file),
            "--cores",
            str(cores),
            "--rerun-triggers",
            "mtime",
            "--quiet",
            "all",
        ]
        environment = os.environ.copy()
        package_parent = str(Path(__file__).resolve().parent.parent)
        current_pythonpath = environment.get("PYTHONPATH")
        environment["PYTHONPATH"] = os.pathsep.join(
            part for part in (package_parent, current_pythonpath) if part
        )
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode != 0:
            diagnostic = result.stderr.strip() or result.stdout.strip() or "no output"
            raise RuntimeError(
                "Annotation workflow failed with exit code "
                f"{result.returncode}: {diagnostic}"
            )

        database_results = []
        for database_index in range(database_count):
            hits = pd.read_pickle(work_dir / "results" / f"{database_index}.pkl")
            if not hits.empty:
                database_results.append(hits)
        return database_results


def _allocate_search_threads(databases: dict, cores: int) -> List[int]:
    """Allocate cores across searches without exceeding the global core limit.

    Every search receives one thread first. Spare threads cycle through Infernal,
    DIAMOND, then BLAST jobs. Rfam is prioritized because it is consistently the
    slowest search across the bundled example plasmids.
    """
    if cores < 1:
        raise ValueError("cores must be at least 1")
    if not databases:
        return []

    allocations = [1] * len(databases)
    spare_cores = max(0, cores - len(databases))
    method_priority = {"infernal": 0, "diamond": 1, "blastn": 2, "blast": 2}
    database_configs = list(databases.values())
    prioritized_indices = sorted(
        range(len(databases)),
        key=lambda index: method_priority.get(
            str(database_configs[index]["method"]).lower(), 3
        ),
    )
    for offset in range(spare_cores):
        allocations[prioritized_indices[offset % len(prioritized_indices)]] += 1
    return allocations


def search_all_databases(
    query_sequence: str,
    is_linear: bool,
    yaml_file: Path,
    cores: int = 1,
) -> pd.DataFrame:
    """Search query sequence against all configured databases.

    This function includes automatic caching for identical inputs to improve performance
    when processing the same DNA sequence multiple times. The cache is applied at the
    database search level so that both detailed and non-detailed annotation modes can
    reuse the same raw search results.
    """
    # Convert inputs to hashable types for caching
    yaml_file_str = str(yaml_file)

    # Use cached implementation for computational efficiency
    return _search_all_databases_cached(query_sequence, is_linear, yaml_file_str, cores)
