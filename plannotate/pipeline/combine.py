"""
Combine module for pLannotate pipeline.
Handles combining results from multiple database searches.
"""

import pandas as pd

from .. import resources as rsc
from ..logging_config import get_logger

logger = get_logger(__name__)


def combine_results(result_files, is_detailed=False):
    """
    Combine results from multiple database searches.
    
    Args:
        result_files: List of paths to result CSV files
        is_detailed: Whether to use detailed annotation types
        
    Returns:
        Combined DataFrame sorted by score
    """
    if not result_files:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    all_results = []
    for file_path in result_files:
        try:
            df = pd.read_csv(file_path)
            if not df.empty:
                all_results.append(df)
        except Exception as e:
            logger.warning(f"Failed to read {file_path}: {e}")
    
    if not all_results:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Set kind column based on detail level
    if is_detailed:
        combined_df["kind"] = combined_df["Type"]
    else:
        combined_df["kind"] = 1
    
    # Sort by score and other metrics
    combined_df = combined_df.sort_values(
        by=["score", "length", "percmatch"], 
        ascending=[False, False, False]
    )
    
    return combined_df


def prepare_sequence(seq, is_linear):
    """
    Prepare sequence for searching.
    
    Args:
        seq: Input DNA sequence
        is_linear: Whether sequence is linear
        
    Returns:
        Prepared sequence (doubled if circular)
    """
    if is_linear:
        return seq
    else:
        # Double sequence for circular plasmids to catch origin-spanning features
        return seq + seq