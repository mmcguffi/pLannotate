#!/bin/bash
# Master script to run the complete Swiss-Prot processing workflow
# Usage: ./run_full_workflow.sh [--download-only | --skip-download]
#
# Required scripts: download_fresh_swissprot.py, process_splits.sh, parse_swissprot_file.py
# Required packages: pip install biopython pandas

set -e

# Handle arguments
DOWNLOAD_ONLY=false
SKIP_DOWNLOAD=false

case "${1:-}" in
    --download-only) DOWNLOAD_ONLY=true ;;
    --skip-download) SKIP_DOWNLOAD=true ;;
    --help|-h) 
        echo "Usage: $0 [--download-only | --skip-download]"
        echo "  --download-only  Only download and split"
        echo "  --skip-download  Only process existing files"
        exit 0
        ;;
    "") ;; # No args, run full workflow
    *) echo "Unknown option: $1"; exit 1 ;;
esac

# Check requirements
python -c "import Bio; import pandas" 2>/dev/null || { echo "Missing packages. Run: pip install biopython pandas"; exit 1; }

# Step 1: Download and split
if [ "$SKIP_DOWNLOAD" = false ]; then
    echo "=== Downloading and splitting Swiss-Prot database ==="
    if [ "$DOWNLOAD_ONLY" = true ]; then
        python download_fresh_swissprot.py --download-only
    else
        python download_fresh_swissprot.py
    fi
fi

# Step 2: Process descriptions
if [ "$DOWNLOAD_ONLY" = false ]; then
    echo "=== Processing split files ==="
    chmod +x process_splits.sh
    ./process_splits.sh ./fresh_split ./fresh_results
fi

echo "=== Complete ==="
echo "FASTA sequences: ./fresh_data/uniprot_sprot.fasta"
[ "$DOWNLOAD_ONLY" = false ] && echo "Descriptions: ./fresh_results/swiss_description*.csv.gz"