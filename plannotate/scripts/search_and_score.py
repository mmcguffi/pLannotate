"""
Search database and calculate scores with feature details
"""

import os
import json
import pandas as pd
from Bio import SeqIO

# Import pipeline modules
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import search, process, details
import resources as rsc

# Get Snakemake variables
seq_file = snakemake.input.seq
seq_info_file = snakemake.input.seq_info
output_file = snakemake.output.result
db_name = snakemake.params.db_name
yaml_file = snakemake.params.yaml_file

# Create output directory
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Read sequence
record = list(SeqIO.parse(seq_file, "fasta"))[0]
seq = str(record.seq)

# Read sequence info
with open(seq_info_file, 'r') as f:
    seq_info = json.load(f)

# Get database configuration
databases = rsc.get_yaml(yaml_file)
db_config = databases[db_name]

# Run search
print(f"Searching database: {db_name}")
results = search.search_database(seq, db_config, db_name)

if results.empty:
    # Save empty results
    results.to_csv(output_file, index=False)
else:
    # Get feature details
    detailed_results = details.get_feature_details(results, yaml_file)
    
    # Calculate scores
    scored_results = process.calculate_scores(detailed_results, seq_info["is_linear"])
    
    # Save results
    scored_results.to_csv(output_file, index=False)
    print(f"Found {len(scored_results)} annotations in {db_name}")