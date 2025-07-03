"""
Generate GenBank file and optional visualization
"""

import json
import pandas as pd

# Import required modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import resources as rsc
from bokeh_plot import get_bokeh
from bokeh.plotting import output_file, save

# Get Snakemake variables
annotations_file = snakemake.input.annotations
seq_info_file = snakemake.input.seq_info
gbk_file = snakemake.output.gbk
html_file = snakemake.output.html
create_plot = snakemake.params.create_plot

# Read data
annotations_df = pd.read_csv(annotations_file)

with open(seq_info_file, 'r') as f:
    seq_info = json.load(f)

# Create GenBank file
print("Creating GenBank file...")
gbk_content = rsc.get_gbk(
    annotations_df, 
    seq_info["original_seq"], 
    is_linear=seq_info["is_linear"]
)

with open(gbk_file, 'w') as f:
    f.write(gbk_content)
print(f"GenBank file saved to: {gbk_file}")

# Create visualization if requested
if create_plot:
    print("Creating visualization...")
    plot = get_bokeh(annotations_df, linear=seq_info["is_linear"])
    
    # Save to HTML
    output_file(html_file)
    save(plot)
    print(f"Visualization saved to: {html_file}")
else:
    # Create empty file to satisfy Snakemake
    with open(html_file, 'w') as f:
        f.write("<html><body>Visualization disabled</body></html>")