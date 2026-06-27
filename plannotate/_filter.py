"""Internal annotation filtering and scoring.

This module contains functions for calculating scores, cleaning hits,
retrieving feature details, and processing raw annotation results.

Author: Matt McGuffie
"""

import logging

import numpy as np
import pandas as pd

from ._schema import ANNOTATION_COLUMNS

logger = logging.getLogger(__name__)

# Constants for filtering and scoring
WIGGLE_RATIO = 0.15  # Percent "trimmed" on either end for overlap detection
PERFECT_MATCH_BONUS = 10  # Bonus multiplier for 100% matches
MIN_QUALITY_THRESHOLD = 3  # Minimum pi_permatch to keep hits
EVALUE_THRESHOLD = 1.0  # Maximum e-value to keep hits

# Known problematic sequence IDs that cause overlap issues
BLACKLISTED_SSEQIDS = ["P03851", "P03845", "ISS", "P03846"]


def _store_original_coordinates(hits: pd.DataFrame) -> pd.DataFrame:
    """Store original coordinates before circular adjustments."""
    # NOTE: these are the pre-wrap (doubled-sequence) coordinates, so the two copies
    # of a feature in a circular query keep distinct dup coordinates and survive the
    # drop_duplicates in _apply_quality_filters. They are collapsed afterwards by
    # same-kind overlap removal, which both copies trigger once wrapped onto the same
    # interval. Any change that lets the copies differ in "kind" would emit duplicates.
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
    hits = hits.loc[~hits["sseqid"].isin(BLACKLISTED_SSEQIDS)]
    hits = hits.loc[hits["evalue"] < EVALUE_THRESHOLD]
    hits = hits.loc[hits["pi_permatch"] > MIN_QUALITY_THRESHOLD]
    return hits.drop_duplicates().reset_index(drop=True)


def _remove_overlapping_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """Remove lower-scoring hits whose trimmed interval overlaps a better hit."""
    if hits.empty:
        return hits

    sequence_length = int(hits["qlen"].iloc[0])
    # hoist the hot columns into plain lists: pulling rows out of the DataFrame
    # inside this O(n^2) loop is dominated by pandas scalar-access overhead, which
    # turns complex plasmids (thousands of raw hits) into hundreds of milliseconds.
    qstart = hits["qstart"].astype(int).tolist()
    qend = hits["qend"].astype(int).tolist()
    wstart = hits["wstart"].astype(int).tolist()
    wend = hits["wend"].astype(int).tolist()
    kind = hits["kind"].tolist()

    count = len(hits)
    dropped: set[int] = set()
    for better_index in range(count):
        if better_index in dropped:
            continue
        better_kind = kind[better_index]
        occupied = _circular_segments(
            qstart[better_index],
            qend[better_index],
            sequence_length,
            exclude_start=True,
        )
        for candidate_index in range(better_index + 1, count):
            if candidate_index in dropped:
                continue
            if kind[candidate_index] != better_kind:
                continue
            trimmed = _circular_segments(
                wstart[candidate_index],
                wend[candidate_index],
                sequence_length,
            )
            if _segments_overlap(occupied, trimmed):
                dropped.add(candidate_index)

    logger.debug("Found %d overlapping hits to remove", len(dropped))
    return hits.drop(index=list(dropped)).reset_index(drop=True)


def _circular_segments(
    start: int,
    end: int,
    sequence_length: int,
    *,
    exclude_start: bool = False,
) -> list[tuple[int, int]]:
    """Represent an inclusive circular interval as one or two linear segments."""
    if exclude_start and start < end:
        start += 1
    if start <= end:
        return [(start, end)]
    segments = [(0, end), (start, sequence_length - 1)]
    return [
        (segment_start, segment_end)
        for segment_start, segment_end in segments
        if segment_start <= segment_end
    ]


def _segments_overlap(
    first: list[tuple[int, int]], second: list[tuple[int, int]]
) -> bool:
    return any(
        max(first_start, second_start) <= min(first_end, second_end)
        for first_start, first_end in first
        for second_start, second_end in second
    )


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
    hits = _normalize_coordinates(hits)
    hits = _calculate_quality_metrics(hits)
    hits = _apply_priority_scoring(hits)
    hits = _apply_perfect_match_bonus(hits)
    return _calculate_wiggle_boundaries(hits)


def filter_and_clean_hits(hits: pd.DataFrame, is_linear: bool = False) -> pd.DataFrame:
    """Clean BLAST hits by removing poor matches and resolving overlaps."""
    if hits.empty:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    initial_count = len(hits)
    logger.debug("Starting hit filtering with %d raw hits", initial_count)

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

        hits = _adjust_circular_coordinates(hits)
        logger.debug("Adjusted coordinates for circular sequences")

    # Apply basic quality filters
    hits = _apply_quality_filters(hits)
    after_quality_filter = len(hits)
    logger.debug(
        "After quality filters: %d hits (%d removed)",
        after_quality_filter,
        initial_count - after_quality_filter,
    )

    if hits.empty:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    # Remove overlapping hits using sequence space analysis
    logger.info("Resolving overlaps among %d hits...", after_quality_filter)
    hits = _remove_overlapping_hits(hits)
    final_count = len(hits)
    logger.debug(
        "After overlap removal: %d hits (%d overlaps removed)",
        final_count,
        after_quality_filter - final_count,
    )

    return hits
