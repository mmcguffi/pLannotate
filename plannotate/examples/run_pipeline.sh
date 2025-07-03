#!/bin/bash
# Example script for running pLannotate Snakemake pipeline

# Example 1: Run with a FASTA file
echo "Example 1: Running pipeline with FASTA file"
SEQUENCE=$(grep -v ">" ../data/fastas/pUC19.fa | tr -d '\n')
snakemake -s ../Snakefile \
    --config input_sequence="$SEQUENCE" \
             output_dir="pUC19_results" \
             linear=False \
    --cores 4

# Example 2: Run with command line sequence
echo -e "\nExample 2: Running pipeline with command line sequence"
snakemake -s ../Snakefile \
    --config input_sequence="ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA" \
             output_dir="test_results" \
             linear=True \
             create_plot=False \
    --cores 2

# Example 3: Run with config file
echo -e "\nExample 3: Running pipeline with config file"
cat > example_config.yaml << EOF
input_sequence: "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTG"
linear: False
is_detailed: True
output_dir: "detailed_results"
create_plot: True
threads: 8
EOF

snakemake -s ../Snakefile \
    --configfile example_config.yaml \
    --cores 8

echo -e "\nPipeline runs complete! Check the output directories for results."