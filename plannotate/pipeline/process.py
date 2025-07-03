"""
Process module for pLannotate pipeline.
Handles score calculation, hit cleaning, and fragment detection.
"""

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from .. import resources as rsc
from ..logging_config import get_logger

logger = get_logger(__name__)


def calculate_scores(df, is_linear):
    """
    Calculate scores and related metrics for hits.
    
    Args:
        df: DataFrame with search results
        is_linear: Whether the sequence is linear (vs circular)
        
    Returns:
        DataFrame with calculated scores
    """
    df = df.copy()
    
    # Adjust positions to 0-based
    df["qstart"] = df["qstart"] - 1
    df["qend"] = df["qend"] - 1
    
    # Ensure qstart < qend
    df["qstart"], df["qend"] = (
        df[["qstart", "qend"]].min(axis=1),
        df[["qstart", "qend"]].max(axis=1),
    )
    
    # Calculate match percentages
    df["percmatch"] = df["length"] / df["slen"] * 100
    df["abs percmatch"] = 100 - abs(100 - df["percmatch"])
    df["pi_permatch"] = (df["pident"] * df["abs percmatch"]) / 100
    df["score"] = (df["pi_permatch"] / 100) * df["length"]
    
    # Apply priority-based score adjustment
    df["score"] = df["score"] * (2 ** (-1 * df["priority"].astype(float)) * 2)
    
    # Adjust qlen for circular sequences
    if not is_linear:
        df["qlen"] = (df["qlen"] / 2).astype("int")
    
    # Apply bonus for perfect matches
    bonus = (1 / df["priority"]) * 10
    df.loc[df["pi_permatch"] == 100, "score"] = (
        df.loc[df["pi_permatch"] == 100, "score"] * bonus
    )
    
    # Calculate wiggle room for overlap detection
    wiggle_size = 0.15
    df["wiggle"] = (df["length"] * wiggle_size).astype(int)
    df["wstart"] = df["qstart"] + df["wiggle"]
    df["wend"] = df["qend"] - df["wiggle"]
    
    return df


def clean_hits(df):
    """
    Clean and filter hits, removing overlaps and low-quality matches.
    
    Args:
        df: DataFrame with scored hits
        
    Returns:
        Cleaned DataFrame
    """
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    df = df.copy()
    
    # Store original positions
    df["qstart_dup"] = df["qstart"]
    df["qend_dup"] = df["qend"]
    
    # Adjust positions that wrap around sequence
    df["qstart"] = np.where(
        df["qstart"] >= df["qlen"], df["qstart"] - df["qlen"], df["qstart"]
    )
    df["qend"] = np.where(
        df["qend"] >= df["qlen"], df["qend"] - df["qlen"], df["qend"]
    )
    df["wstart"] = np.where(
        df["wstart"] >= df["qlen"], df["wstart"] - df["qlen"], df["wstart"]
    )
    df["wend"] = np.where(
        df["wend"] >= df["qlen"], df["wend"] - df["qlen"], df["wend"]
    )
    
    # Remove known problematic hits
    problem_hits = ["P03851", "P03845", "ISS", "P03846"]
    df = df.loc[~df["sseqid"].isin(problem_hits)]
    
    # Filter by e-value and quality
    df = df.loc[df["evalue"] < 1]
    df = df.loc[df["pi_permatch"] > 3]
    
    df = df.drop_duplicates().reset_index(drop=True)
    
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Convert float columns to int where appropriate
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col], downcast="integer")
        except (ValueError, TypeError):
            pass
    
    # Remove overlapping hits
    df = _remove_overlaps(df)
    
    return df


def _remove_overlaps(df):
    """Remove overlapping hits using sequence space concept."""
    if df.empty:
        return df
        
    # Create conceptual sequence space
    seq_space = []
    end = int(df["qlen"].iloc[0])
    
    for i in df.index:
        wstart = df.loc[i]["wstart"]
        wend = df.loc[i]["wend"]
        sseqid = [df.loc[i]["sseqid"]]
        
        if wend < wstart:  # Hit crosses origin
            left = (wend + 1) * [df.loc[i]["kind"]]
            center = (wstart - wend - 1) * [None]
            right = (end - wstart) * [df.loc[i]["kind"]]
        else:  # Normal hit
            left = wstart * [None]
            center = (wend - wstart + 1) * [df.loc[i]["kind"]]
            right = (end - wend - 1) * [None]
        
        seq_space.append(sseqid + left + center + right)
    
    seq_space = pd.DataFrame(seq_space, columns=["sseqid"] + list(range(0, end)))
    seq_space = seq_space.set_index([seq_space.index, "sseqid"])
    
    # Find overlaps and determine which to drop
    to_drop = set()
    for i in range(len(seq_space)):
        if seq_space.iloc[i].name in to_drop:
            continue
            
        qstart = df.loc[seq_space.iloc[i].name[0]]["qstart"]
        qend = df.loc[seq_space.iloc[i].name[0]]["qend"]
        kind = df.loc[seq_space.iloc[i].name[0]]["kind"]
        
        if qstart < qend:
            column_slice = list(range(qstart + 1, qend + 1))
        else:
            column_slice = list(range(0, qend + 1)) + list(range(qstart, end))
        
        row_slice = (seq_space[column_slice] == kind).any(axis=1)
        to_drop = to_drop | set(seq_space[row_slice].loc[i + 1:].index)
    
    seq_space = seq_space.drop(to_drop)
    df = df.loc[seq_space.index.get_level_values(0)]
    df = df.reset_index(drop=True)
    
    return df


def detect_fragments(df):
    """
    Determine which features are fragments based on type and match quality.
    
    Args:
        df: DataFrame with cleaned hits
        
    Returns:
        DataFrame with fragment column added
    """
    def is_fragment(feature):
        if feature["Type"] == "CDS":
            if feature["pi_permatch"] == 100:
                return False
            elif ((feature["length"] % 3) == 0) & (feature["percmatch"] > 95):
                return False
            else:
                return True
        elif feature["Type"] != "CDS":
            if feature["percmatch"] < 95:
                return True
            else:
                return False
        else:
            logger.error("Fragment error.")
            return True
    
    df["fragment"] = df.apply(is_fragment, axis=1)
    return df


def finalize_annotations(df, seq):
    """
    Final processing of annotations including sequence extraction.
    
    Args:
        df: DataFrame with annotations
        seq: Original sequence (not doubled)
        
    Returns:
        Finalized DataFrame
    """
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Adjust end position for GBK format
    df["qend"] = df["qend"] + 1
    
    # Handle reverse complement for negative strand features
    df["qseq"] = df.apply(
        lambda x: (
            str(Seq(x["qseq"]).reverse_complement()) 
            if x["sframe"] == -1 
            else x["qseq"]
        ),
        axis=1,
    )
    
    # Fill in edge cases
    df["Feature"] = df["Feature"].fillna(df["sseqid"])
    df["Description"] = df["Description"].fillna("")
    df["Type"] = df["Type"].fillna("misc_feature")
    
    return df