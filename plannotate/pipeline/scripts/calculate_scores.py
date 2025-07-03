"""
Calculate scores and adjusted positions for pLannotate hits
"""
import pandas as pd
import json


def calculate_scores(df, is_linear):
    """Calculate scores and positions for hits"""
    if df.empty:
        return df
    
    # Adjust to 0-based coordinates
    df["qstart"] = df["qstart"] - 1
    df["qend"] = df["qend"] - 1
    
    # Ensure start < end
    df["qstart"], df["qend"] = (
        df[["qstart", "qend"]].min(axis=1),
        df[["qstart", "qend"]].max(axis=1),
    )
    
    # Calculate match percentages
    df["percmatch"] = df["length"] / df["slen"] * 100
    df["abs percmatch"] = 100 - abs(100 - df["percmatch"])
    df["pi_permatch"] = (df["pident"] * df["abs percmatch"]) / 100
    df["score"] = (df["pi_permatch"] / 100) * df["length"]
    
    # Score adjustment based on priority
    # Higher priority = less score deduction
    df["score"] = df["score"] * (2 ** (-1 * df["priority"].astype(float)) * 2)
    
    # Adjust qlen for circular sequences
    if not is_linear:
        df["qlen"] = (df["qlen"] / 2).astype("int")
    
    # Apply bonus for perfect matches
    bonus = (1 / df["priority"]) * 10
    perfect_mask = df["pi_permatch"] == 100
    df.loc[perfect_mask, "score"] = df.loc[perfect_mask, "score"] * bonus
    
    # Calculate wiggle room for overlap detection
    wiggle_size = 0.15  # 15% trimmed on each end
    df["wiggle"] = (df["length"] * wiggle_size).astype(int)
    df["wstart"] = df["qstart"] + df["wiggle"]
    df["wend"] = df["qend"] - df["wiggle"]
    
    return df


def main(input_file, metadata_file, output_file, db_config, is_linear):
    """Main function to calculate scores"""
    # Read hits
    df = pd.read_csv(input_file, sep="\t")
    
    if df.empty:
        df.to_csv(output_file, sep="\t", index=False)
        return
    
    # Read metadata for sequence length info
    with open(metadata_file) as f:
        metadata = json.load(f)
    
    # Calculate scores
    scored_df = calculate_scores(df, is_linear)
    
    # Save results
    scored_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_file=snakemake.input.detailed_hits,
        metadata_file=snakemake.input.metadata,
        output_file=snakemake.output[0],
        db_config=snakemake.params.db_config,
        is_linear=snakemake.params.linear
    )