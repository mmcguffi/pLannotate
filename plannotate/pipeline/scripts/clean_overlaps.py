"""
Clean overlapping features based on score priority
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
import resources as rsc


def clean_overlaps(df, is_detailed):
    """Remove overlapping features using sequence space algorithm"""
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Store original positions for circular sequences
    df["qstart_dup"] = df["qstart"]
    df["qend_dup"] = df["qend"]
    
    # Adjust positions if they exceed sequence length
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
    
    # Filter by e-value
    df = df.loc[df["evalue"] < 1]
    
    # Drop poor matches that are very small fragments
    df = df.loc[df["pi_permatch"] > 3]
    
    # Remove duplicates
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)
    
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Set kind column for overlap detection
    if is_detailed:
        df["kind"] = df["Type"]
    else:
        df["kind"] = 1
    
    # Create conceptual sequence space
    seq_length = int(df["qlen"].iloc[0])
    seq_space = []
    
    # Ensure integer types
    df = df.apply(pd.to_numeric, errors="ignore", downcast="integer")
    
    # Build sequence space representation
    for idx in df.index:
        wstart = df.loc[idx]["wstart"]
        wend = df.loc[idx]["wend"]
        sseqid = [df.loc[idx]["sseqid"]]
        
        if wend < wstart:  # Feature crosses origin
            left = (wend + 1) * [df.loc[idx]["kind"]]
            center = (wstart - wend - 1) * [None]
            right = (seq_length - wstart) * [df.loc[idx]["kind"]]
        else:  # Normal feature
            left = wstart * [None]
            center = (wend - wstart + 1) * [df.loc[idx]["kind"]]
            right = (seq_length - wend - 1) * [None]
        
        seq_space.append(sseqid + left + center + right)
    
    # Create multi-indexed dataframe
    seq_space_df = pd.DataFrame(seq_space, columns=["sseqid"] + list(range(seq_length)))
    seq_space_df = seq_space_df.set_index([seq_space_df.index, "sseqid"])
    
    # Filter overlapping features
    to_drop = set()
    
    for i in range(len(seq_space_df)):
        if seq_space_df.iloc[i].name in to_drop:
            continue
        
        qstart = df.loc[seq_space_df.iloc[i].name[0]]["qstart"]
        qend = df.loc[seq_space_df.iloc[i].name[0]]["qend"]
        kind = df.loc[seq_space_df.iloc[i].name[0]]["kind"]
        
        # Get column slice for this feature
        if qstart < qend:
            column_slice = list(range(qstart + 1, qend + 1))
        else:  # Crosses origin
            column_slice = list(range(0, qend + 1)) + list(range(qstart, seq_length))
        
        # Find overlapping features
        row_slice = (seq_space_df[column_slice] == kind).any(axis=1)
        overlapping = seq_space_df[row_slice].loc[i + 1:]
        to_drop.update(overlapping.index)
    
    # Drop overlapping features
    seq_space_df = seq_space_df.drop(to_drop)
    
    # Filter original dataframe
    keep_indices = seq_space_df.index.get_level_values(0)
    df = df.loc[keep_indices]
    df = df.reset_index(drop=True)
    
    return df


def main(input_file, metadata_file, output_file, is_detailed):
    """Main function to clean overlaps"""
    # Read hits
    df = pd.read_csv(input_file, sep="\t")
    
    # Read metadata
    with open(metadata_file) as f:
        metadata = json.load(f)
    
    # Clean overlaps
    cleaned_df = clean_overlaps(df, is_detailed)
    
    # Save results
    cleaned_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_file=snakemake.input.all_hits,
        metadata_file=snakemake.input.metadata,
        output_file=snakemake.output[0],
        is_detailed=snakemake.params.detailed
    )