"""Annotation filtering and processing functions.

This module contains functions for calculating scores, cleaning hits,
retrieving feature details, and processing raw annotation results.

Author: Matt McGuffie
"""

import numpy as np
import pandas as pd

from .logging_config import get_logger
from .search import DF_COLS

logger = get_logger(__name__)

# Constants for filtering and scoring
WIGGLE_RATIO = 0.15  # Percent "trimmed" on either end for overlap detection
PERFECT_MATCH_BONUS = 10  # Bonus multiplier for 100% matches
MIN_QUALITY_THRESHOLD = 3  # Minimum pi_permatch to keep hits
EVALUE_THRESHOLD = 1.0  # Maximum e-value to keep hits

# Known problematic sequence IDs that cause overlap issues
BLACKLISTED_SSEQIDS = ["P03851", "P03845", "ISS", "P03846"]


def _store_original_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Store original coordinates before circular adjustments."""
    hits["qstart_dup"] = hits["qstart"]
    hits["qend_dup"] = hits["qend"]
    return hits


def _adjust_circular_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Adjust coordinates that exceed sequence length for circular sequences."""
    # By this point, qlen should already be adjusted to the original sequence length
    # We adjust coordinates that exceed this length (from the second half of doubled sequence)
    seq_length = hits["qlen"].iloc[0]

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
    logger.debug("Normalized numeric columns for overlap analysis")

    # Build sequence space representation
    sequence_space = _build_sequence_space_matrix(hits)
    logger.debug(f"Built sequence space matrix for {len(hits)} hits")

    # Find and remove overlapping hits
    overlapping_indices = _find_overlapping_indices_original(hits, sequence_space)
    logger.debug(f"Found {len(overlapping_indices)} overlapping hits to remove")
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


def calculate_hit_scores(hits: pd.DataFrame) -> pd.DataFrame:
    """Calculate quality metrics and scores for BLAST hits."""
    hits = hits.copy()

    # Convert to 0-based coordinates
    hits = _normalize_coordinates(hits)

    # Calculate basic quality metrics
    hits = _calculate_quality_metrics(hits)

    # Apply priority-based score adjustments
    hits = _apply_priority_scoring(hits)

    # Apply perfect match bonus
    hits = _apply_perfect_match_bonus(hits)

    # Calculate wiggle room for overlap detection
    hits = _calculate_wiggle_boundaries(hits)

    return hits


def filter_and_clean_hits(hits: pd.DataFrame, is_linear: bool = False) -> pd.DataFrame:
    """Clean BLAST hits by removing poor matches and resolving overlaps."""
    if hits.empty:
        return pd.DataFrame(columns=DF_COLS)

    initial_count = len(hits)
    logger.debug(f"Starting hit filtering with {initial_count} raw hits")

    # Calculate scores and quality metrics first
    hits = calculate_hit_scores(hits)
    logger.debug("Calculated hit scores and quality metrics")

    # Sort by quality metrics for consistent processing order
    hits = hits.sort_values(
        by=["score", "length", "percmatch"], ascending=[False, False, False]
    ).reset_index(drop=True)
    logger.debug("Sorted hits by score and quality metrics")

    # Store original coordinates before circular adjustments
    hits = _store_original_coordinates(hits)
    logger.debug("Stored original coordinates")

    # Adjust qlen for circular sequences before coordinate adjustment
    if not is_linear:
        hits["qlen"] = (hits["qlen"] / 2).astype("int")
        logger.debug("Adjusted qlen for circular sequences")
    
    # Adjust coordinates for circular sequences
    hits = _adjust_circular_coordinates(hits)
    logger.debug("Adjusted coordinates for circular sequences")

    # Apply basic quality filters
    hits = _apply_quality_filters(hits)
    after_quality_filter = len(hits)
    logger.debug(
        f"After quality filters: {after_quality_filter} hits ({initial_count - after_quality_filter} removed)"
    )

    if hits.empty:
        return pd.DataFrame(columns=DF_COLS)

    # Remove overlapping hits using sequence space analysis
    logger.info(f"Resolving overlaps among {after_quality_filter} hits...")
    hits = _remove_overlapping_hits(hits)
    final_count = len(hits)
    logger.debug(
        f"After overlap removal: {final_count} hits ({after_quality_filter - final_count} overlaps removed)"
    )

    return hits
