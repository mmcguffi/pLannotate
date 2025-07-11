#!/bin/bash
# Process Swiss-Prot split files to extract descriptions using parse_swissprot_file.py
# Usage: ./process_splits.sh <split_dir> <output_dir>
# 
# This script:
# 1. Processes each split file with parse_swissprot_file.py (both verbose and non-verbose)
# 2. Combines results into swiss_description.csv and swiss_description_verbose.csv
# 3. Compresses final CSV files

SPLIT_DIR="${1:-./fresh_split}"
OUTPUT_DIR="${2:-./fresh_results}"

echo "Processing split files from $SPLIT_DIR to $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Clean up any existing files
rm -f "$OUTPUT_DIR"/swiss_description*.csv

# Process each split file
for file in "$SPLIT_DIR"/chunk_*.dat; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        echo "Processing $filename..."
        
        # Non-verbose
        python parse_swissprot_file.py "$filename" 2>/dev/null || echo "Warning: Failed $filename"
        if [ -f "./swiss_description.csv" ]; then
            cat "./swiss_description.csv" >> "$OUTPUT_DIR/swiss_description.csv"
            rm "./swiss_description.csv"
        fi
        
        # Verbose
        python parse_swissprot_file.py -v "$filename" 2>/dev/null || echo "Warning: Failed $filename verbose"
        if [ -f "./swiss_description.csv" ]; then
            cat "./swiss_description.csv" >> "$OUTPUT_DIR/swiss_description_verbose.csv"
            rm "./swiss_description.csv"
        fi
    fi
done

# Compress results
echo "Compressing results..."
[ -f "$OUTPUT_DIR/swiss_description.csv" ] && gzip "$OUTPUT_DIR/swiss_description.csv"
[ -f "$OUTPUT_DIR/swiss_description_verbose.csv" ] && gzip "$OUTPUT_DIR/swiss_description_verbose.csv"

echo "Done! Results in $OUTPUT_DIR/"