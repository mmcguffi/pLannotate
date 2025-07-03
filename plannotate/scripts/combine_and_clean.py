"""
Combine all database results and clean overlapping hits
"""

import json
import pandas as pd

# Import pipeline modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import combine, process
import resources as rsc

# Get Snakemake variables
result_files = snakemake.input.results
seq_info_file = snakemake.input.seq_info
output_file = snakemake.output.combined
is_detailed = snakemake.params.is_detailed

# Read sequence info
with open(seq_info_file, 'r') as f:
    seq_info = json.load(f)

# Combine all results
print(f"Combining results from {len(result_files)} databases...")
combined_df = combine.combine_results(result_files, is_detailed)

if combined_df.empty:
    print("No annotations found in any database")
    combined_df.to_csv(output_file, index=False)
else:
    print(f"Combined {len(combined_df)} total annotations")
    
    # Clean overlapping hits
    print("Cleaning overlapping hits...")
    cleaned_df = process.clean_hits(combined_df)
    
    print(f"Retained {len(cleaned_df)} annotations after cleaning")
    
    # Save cleaned results
    cleaned_df.to_csv(output_file, index=False)