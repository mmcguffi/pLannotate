"""
Finalize annotations with fragment detection and sequence extraction
"""

import json
import pandas as pd

# Import pipeline modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import process
import resources as rsc

# Get Snakemake variables
annotations_file = snakemake.input.annotations
seq_info_file = snakemake.input.seq_info
output_file = snakemake.output.final

# Read data
annotations_df = pd.read_csv(annotations_file)

with open(seq_info_file, 'r') as f:
    seq_info = json.load(f)

if annotations_df.empty:
    print("No annotations to finalize")
    final_df = pd.DataFrame(columns=rsc.DF_COLS)
else:
    print(f"Finalizing {len(annotations_df)} annotations...")
    
    # Detect fragments
    print("Detecting fragments...")
    annotations_df = process.detect_fragments(annotations_df)
    
    num_fragments = annotations_df['fragment'].sum()
    print(f"Found {num_fragments} fragment annotations")
    
    # Finalize annotations
    print("Finalizing annotations...")
    final_df = process.finalize_annotations(annotations_df, seq_info["original_seq"])
    
    print(f"Final annotation count: {len(final_df)}")

# Save final annotations
final_df.to_csv(output_file, index=False)