"""
Merge results from all databases and sort by score
"""
import pandas as pd
from pathlib import Path


def merge_database_results(input_files):
    """Merge all database results into a single dataframe"""
    all_dfs = []
    
    for file_path in input_files:
        if Path(file_path).exists():
            df = pd.read_csv(file_path, sep="\t")
            if not df.empty:
                all_dfs.append(df)
    
    if not all_dfs:
        # Return empty dataframe with expected columns
        return pd.DataFrame()
    
    # Concatenate all results
    merged_df = pd.concat(all_dfs, ignore_index=True)
    
    # Sort by score (descending), then by length and percmatch
    merged_df = merged_df.sort_values(
        by=["score", "length", "percmatch"],
        ascending=[False, False, False]
    )
    
    return merged_df


def main(input_files, output_file):
    """Main function to merge database results"""
    # Merge all input files
    merged_df = merge_database_results(input_files)
    
    # Save results
    merged_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_files=snakemake.input,
        output_file=snakemake.output[0]
    )