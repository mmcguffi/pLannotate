"""Annotation filtering and processing functions.

This module contains functions for calculating scores, cleaning hits,
retrieving feature details, and processing raw annotation results.

Author: Matt McGuffie
"""

import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List

import numpy as np
import pandas as pd

from . import resources as rsc
from .logging_config import get_logger

logger = get_logger(__name__)

# Constants
WIGGLE_RATIO = 0.15  # Percent "trimmed" on either end for overlap detection
PERFECT_MATCH_BONUS = 10  # Bonus multiplier for 100% matches
MIN_QUALITY_THRESHOLD = 3  # Minimum pi_permatch to keep hits
EVALUE_THRESHOLD = 1.0  # Maximum e-value to keep hits

# Known problematic sequence IDs that cause overlap issues
BLACKLISTED_SSEQIDS = ["P03851", "P03845", "ISS", "P03846"]

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
    "Feature",
    "Description",
    "Type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "fragment",
]

DEV_COLS = ["wiggle", "wstart", "wend", "kind", "qstart_dup", "qend_dup"]

DF_COLS = BLAST_COLS + EXTRA_COLS + DEV_COLS


def calculate_hit_scores(hits: pd.DataFrame, is_linear: bool) -> pd.DataFrame:
    """Calculate quality metrics and scores for BLAST hits."""
    hits = hits.copy()

    # Convert to 0-based coordinates
    hits = _normalize_coordinates(hits)

    # Calculate basic quality metrics
    hits = _calculate_quality_metrics(hits)

    # Apply priority-based score adjustments
    hits = _apply_priority_scoring(hits)

    # Adjust for linear vs circular sequences
    if not is_linear:
        hits["qlen"] = (hits["qlen"] / 2).astype("int")

    # Apply perfect match bonus
    hits = _apply_perfect_match_bonus(hits)

    # Calculate wiggle room for overlap detection
    hits = _calculate_wiggle_boundaries(hits)

    return hits


def _normalize_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Convert 1-based coordinates to 0-based and ensure start < end."""
    hits["qstart"] = hits["qstart"] - 1
    hits["qend"] = hits["qend"] - 1

    # Ensure qstart is always the smaller coordinate
    hits["qstart"], hits["qend"] = (
        hits[["qstart", "qend"]].min(axis=1),
        hits[["qstart", "qend"]].max(axis=1),
    )
    return hits


def _calculate_quality_metrics(hits: pd.DataFrame) -> pd.DataFrame:
    """Calculate percentage match and weighted quality scores."""
    hits["percmatch"] = hits["length"] / hits["slen"] * 100
    hits["abs percmatch"] = 100 - abs(100 - hits["percmatch"])  # Cap at 100%
    hits["pi_permatch"] = (hits["pident"] * hits["abs percmatch"]) / 100
    hits["score"] = (hits["pi_permatch"] / 100) * hits["length"]
    return hits


def _apply_priority_scoring(hits: pd.DataFrame) -> pd.DataFrame:
    """Apply priority-based score adjustments (higher priority = less penalty)."""
    # Each priority increase halves the score penalty
    # Priority 1 = 1x, Priority 2 = 0.5x, Priority 3 = 0.25x, etc.
    priority_multiplier = 2 ** (-1 * hits["priority"].astype(float)) * 2
    hits["score"] = hits["score"] * priority_multiplier
    return hits


def _apply_perfect_match_bonus(hits: pd.DataFrame) -> pd.DataFrame:
    """Apply bonus scoring for perfect matches."""
    perfect_matches = hits["pi_permatch"] == 100
    if perfect_matches.any():
        bonus = (1 / hits.loc[perfect_matches, "priority"]) * PERFECT_MATCH_BONUS
        hits.loc[perfect_matches, "score"] *= bonus
    return hits


def _calculate_wiggle_boundaries(hits: pd.DataFrame) -> pd.DataFrame:
    """Calculate trimmed boundaries for overlap detection."""
    hits["wiggle"] = (hits["length"] * WIGGLE_RATIO).astype(int)
    hits["wstart"] = hits["qstart"] + hits["wiggle"]
    hits["wend"] = hits["qend"] - hits["wiggle"]
    return hits


def filter_and_clean_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """Clean BLAST hits by removing poor matches and resolving overlaps."""
    if hits.empty:
        return pd.DataFrame(columns=DF_COLS)

    # Store original coordinates before circular adjustments
    hits = _store_original_coordinates(hits)

    # Adjust coordinates for circular sequences
    hits = _adjust_circular_coordinates(hits)

    # Apply basic quality filters
    hits = _apply_quality_filters(hits)

    if hits.empty:
        return pd.DataFrame(columns=DF_COLS)

    # Remove overlapping hits using sequence space analysis
    hits = _remove_overlapping_hits(hits)

    return hits


def _store_original_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Store original coordinates before circular adjustments."""
    hits["qstart_dup"] = hits["qstart"]
    hits["qend_dup"] = hits["qend"]
    return hits


def _adjust_circular_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Adjust coordinates that exceed sequence length for circular sequences."""
    seq_length = hits["qlen"]

    # Adjust main coordinates
    hits["qstart"] = np.where(
        hits["qstart"] >= seq_length, hits["qstart"] - seq_length, hits["qstart"]
    )
    hits["qend"] = np.where(
        hits["qend"] >= seq_length, hits["qend"] - seq_length, hits["qend"]
    )

    # Adjust wiggle boundaries
    hits["wstart"] = np.where(
        hits["wstart"] >= seq_length, hits["wstart"] - seq_length, hits["wstart"]
    )
    hits["wend"] = np.where(
        hits["wend"] >= seq_length, hits["wend"] - seq_length, hits["wend"]
    )

    return hits


def _apply_quality_filters(hits: pd.DataFrame) -> pd.DataFrame:
    """Apply basic quality and blacklist filters."""
    # Remove blacklisted sequences
    hits = hits.loc[~hits["sseqid"].isin(BLACKLISTED_SSEQIDS)]

    # Filter by e-value
    hits = hits.loc[hits["evalue"] < EVALUE_THRESHOLD]

    # Remove very poor matches (usually artifacts)
    hits = hits.loc[hits["pi_permatch"] > MIN_QUALITY_THRESHOLD]

    # Remove duplicates
    hits = hits.drop_duplicates().reset_index(drop=True)

    return hits


def _remove_overlapping_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """Remove overlapping hits using sequence space analysis."""
    if hits.empty:
        return hits

    # Convert numeric columns to proper dtypes
    hits = _normalize_numeric_columns(hits)

    # Build sequence space representation
    sequence_space = _build_sequence_space_matrix(hits)

    # Find and remove overlapping hits
    overlapping_indices = _find_overlapping_indices_original(hits, sequence_space)
    sequence_space = sequence_space.drop(list(overlapping_indices))

    # Return filtered hits
    hits = hits.loc[sequence_space.index.get_level_values(0)]
    hits = hits.reset_index(drop=True)

    return hits


def _normalize_numeric_columns(hits: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric columns to proper integer dtypes."""
    for col in hits.columns:
        try:
            hits[col] = pd.to_numeric(hits[col], downcast="integer")
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass
    return hits


def _build_sequence_space_matrix(hits: pd.DataFrame) -> pd.DataFrame:
    """Build a matrix representing sequence space occupancy for each hit."""
    sequence_length = int(hits["qlen"].iloc[0])
    space_rows = []

    for i in hits.index:
        hit_row = hits.loc[i]
        wstart = hit_row["wstart"]
        wend = hit_row["wend"]
        sseqid = [hit_row["sseqid"]]
        kind = hit_row["kind"]

        # Create sequence space representation for this hit
        if wend < wstart:  # Hit crosses origin (circular sequence)
            left_segment = (wend + 1) * [kind]
            middle_segment = (wstart - wend - 1) * [None]
            right_segment = (sequence_length - wstart) * [kind]
        else:  # Normal linear hit
            left_segment = wstart * [None]
            middle_segment = (wend - wstart + 1) * [kind]
            right_segment = (sequence_length - wend - 1) * [None]

        space_rows.append(sseqid + left_segment + middle_segment + right_segment)

    # Create DataFrame with sequence positions as columns
    columns = ["sseqid"] + list(range(0, sequence_length))
    sequence_space = pd.DataFrame(space_rows, columns=columns)
    sequence_space = sequence_space.set_index([sequence_space.index, "sseqid"])

    return sequence_space


def _find_overlapping_indices_original(
    hits: pd.DataFrame, sequence_space: pd.DataFrame
) -> set:
    """Find indices of hits that overlap with higher-scoring hits (original logic)."""
    sequence_length = hits["qlen"].iloc[0]
    indices_to_drop = set()

    for i in range(len(sequence_space)):
        if sequence_space.iloc[i].name in indices_to_drop:
            continue

        # Get hit information
        hit_index = sequence_space.iloc[i].name[0]
        qstart = hits.loc[hit_index]["qstart"]
        qend = hits.loc[hit_index]["qend"]
        kind = hits.loc[hit_index]["kind"]

        # Get sequence positions occupied by this hit
        if qstart < qend:
            occupied_positions = list(range(qstart + 1, qend + 1))
        else:
            # Hit crosses origin
            occupied_positions = list(range(0, qend + 1)) + list(
                range(qstart, sequence_length)
            )

        # Find hits that overlap with this one in the same positions
        overlapping_mask = (sequence_space[occupied_positions] == kind).any(axis=1)
        overlapping_hits = sequence_space[overlapping_mask]

        # Mark lower-priority overlapping hits for removal
        indices_to_drop.update(overlapping_hits.loc[i + 1 :].index)

    return indices_to_drop


def load_feature_details(hits: pd.DataFrame, yaml_file: Path) -> pd.DataFrame:
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
    """Load feature descriptions from database files."""
    db_details = database_config["details"]

    if db_details["location"] == "None":
        # Data already in dataframe (e.g., Rfam)
        return hits[["sseqid", "Feature", "Description"]].copy()

    # Determine file location
    if db_details["location"] == "Default":
        details_file = rsc.get_details(database_name).with_suffix(".csv")
    else:
        details_file = Path(db_details["location"])

    # Load from compressed or uncompressed file
    if db_details.get("compressed", False):
        details_file = details_file.with_suffix(".csv.gz")
        return _load_compressed_details(sequence_ids, str(details_file))
    else:
        return pd.read_csv(details_file)


def _load_compressed_details(sequence_ids: List[str], file_path: str) -> pd.DataFrame:
    """Load feature details from compressed files using ripgrep."""
    search_pattern = "|".join(sequence_ids)

    with NamedTemporaryFile(suffix=".csv") as temp_file:
        # Use ripgrep to search compressed file
        subprocess.call(
            f'rg -z "{search_pattern}" {file_path} > {temp_file.name}', shell=True
        )

        return pd.read_csv(
            temp_file.name, header=None, names=["sseqid", "Feature", "Description"]
        )


def _process_swissprot_details(feature_details: pd.DataFrame) -> pd.DataFrame:
    """Extract SwissProt protein existence levels for priority modification."""
    # Calculate priority modifications based on existence levels
    priority_modifications = [
        _extract_existence_level_priority(description)
        for description in feature_details["Description"]
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
        feature_details["Type"] = default_type

    return feature_details


def _combine_and_filter_results(database_results: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine results from all databases (matching original behavior exactly)."""
    # Combine all results
    combined_results = pd.concat(database_results, ignore_index=True)

    # Sort by quality metrics (matching original get_raw_hits behavior)
    combined_results = combined_results.sort_values(
        by=["score", "length", "percmatch"], ascending=[False, False, False]
    )

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
