"""
Finalize annotations with sequence extraction and formatting
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
import resources as rsc


def finalize_annotations(df, sequence):
    """Final formatting and sequence extraction"""
    if df.empty:
        return pd.DataFrame(columns=rsc.DF_COLS)
    
    # Adjust end position for GenBank format
    df["qend"] = df["qend"] + 1
    
    # Extract DNA sequences
    df["qseq"] = df.apply(
        lambda row: extract_sequence(row, sequence), axis=1
    )
    
    # Fill in missing values
    df["Feature"] = df["Feature"].fillna(df["sseqid"])
    df["Description"] = df["Description"].fillna("")
    df["Type"] = df["Type"].fillna("misc_feature")
    
    return df


def extract_sequence(row, full_sequence):
    """Extract sequence for a feature, handling reverse complement"""
    seq = row["qseq"]
    
    # If on negative strand, reverse complement
    if row["sframe"] == -1:
        seq = str(Seq(seq).reverse_complement())
    
    return seq


def create_clean_csv(df):
    """Create clean CSV output for users"""
    if df.empty:
        return pd.DataFrame()
    
    # Select and rename columns for user-friendly output
    csv_columns = {
        "Feature": "Feature",
        "Type": "Type",
        "qstart": "Start",
        "qend": "End",
        "sframe": "Strand",
        "Description": "Description",
        "db": "Database",
        "pident": "Identity%",
        "percmatch": "Coverage%",
        "fragment": "Fragment",
        "evalue": "E-value"
    }
    
    clean_df = df[list(csv_columns.keys())].copy()
    clean_df = clean_df.rename(columns=csv_columns)
    
    # Convert strand to +/-
    clean_df["Strand"] = clean_df["Strand"].replace({1: "+", -1: "-"})
    
    # Adjust positions to 1-based for user display
    clean_df["Start"] = clean_df["Start"] + 1
    
    return clean_df


def main(input_fragments, input_fasta, output_csv, output_internal):
    """Main function to finalize annotations"""
    # Read fragments dataframe
    df = pd.read_csv(input_fragments, sep="\t")
    
    # Read original sequence (doubled if circular)
    record = SeqIO.read(input_fasta, "fasta")
    sequence = str(record.seq)
    
    # Finalize annotations
    final_df = finalize_annotations(df, sequence)
    
    # Save internal format
    final_df.to_csv(output_internal, sep="\t", index=False)
    
    # Create and save user-friendly CSV
    clean_df = create_clean_csv(final_df)
    clean_df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    main(
        input_fragments=snakemake.input.fragments,
        input_fasta=snakemake.input.fasta,
        output_csv=snakemake.output.csv,
        output_internal=snakemake.output.internal
    )