"""
Create GenBank file with pLannotate annotations
"""
import pandas as pd
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import date
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
from plannotate import __version__ as plannotate_version


def create_feature_location(row):
    """Create a FeatureLocation object, handling origin-spanning features"""
    if row["qend"] > row["qstart"]:
        # Normal feature
        return FeatureLocation(row["qstart"], row["qend"], row["sframe"])
    else:
        # Origin-spanning feature
        first = FeatureLocation(row["qstart"], row["qlen"], row["sframe"])
        second = FeatureLocation(0, row["qend"], row["sframe"])
        if row["sframe"] == 1 or row["sframe"] == 0:
            return first + second
        else:  # sframe == -1
            return second + first


def add_fragment_label(row):
    """Add (fragment) suffix to fragment features"""
    if row["fragment"]:
        return f"{row['Feature']} (fragment)"
    else:
        return row["Feature"]


def create_genbank(annotations_df, original_seq, metadata, is_linear):
    """Create GenBank record with annotations"""
    # Create base record
    record = SeqRecord(
        seq=Seq(original_seq),
        id=metadata.get("original_id", "sequence"),
        name=metadata.get("original_name", "sequence"),
        description="Annotated by pLannotate"
    )
    
    # Add annotations
    record.annotations["molecule_type"] = "DNA"
    record.annotations["data_file_division"] = "SYN"
    record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()
    record.annotations["accession"] = "."
    record.annotations["version"] = "."
    record.annotations["comment"] = f"Annotated with pLannotate v{plannotate_version}"
    record.annotations["topology"] = "linear" if is_linear else "circular"
    
    if annotations_df.empty:
        return record
    
    # Add feature locations
    annotations_df["feat_loc"] = annotations_df.apply(create_feature_location, axis=1)
    
    # Add fragment labels
    annotations_df["Feature"] = annotations_df.apply(add_fragment_label, axis=1)
    
    # Fix type names
    annotations_df["Type"] = annotations_df["Type"].str.replace("origin of replication", "rep_origin")
    
    # Add features to record
    for _, row in annotations_df.iterrows():
        feature = SeqFeature(
            location=row["feat_loc"],
            type=row["Type"],
            qualifiers={
                "note": "pLannotate",
                "label": row["Feature"],
                "database": row["db"],
                "identity": round(row["pident"], 1),
                "match_length": round(row["percmatch"], 1),
                "fragment": str(row["fragment"]),
                "other": row["Type"]
            }
        )
        record.features.append(feature)
    
    return record


def main(annotations_file, fasta_file, metadata_file, output_file, is_linear):
    """Main function to create GenBank file"""
    # Read annotations
    annotations_df = pd.read_csv(annotations_file, sep="\t")
    
    # Read metadata
    with open(metadata_file) as f:
        metadata = json.load(f)
    
    # Read original sequence (not doubled)
    original_length = metadata["original_length"]
    record = SeqIO.read(fasta_file, "fasta")
    
    # Get original sequence (undouble if needed)
    if metadata["doubled"]:
        original_seq = str(record.seq)[:original_length]
    else:
        original_seq = str(record.seq)
    
    # Create GenBank record
    gbk_record = create_genbank(annotations_df, original_seq, metadata, is_linear)
    
    # Write GenBank file
    with open(output_file, "w") as f:
        SeqIO.write(gbk_record, f, "genbank")


if __name__ == "__main__":
    main(
        annotations_file=snakemake.input.annotations,
        fasta_file=snakemake.input.fasta,
        metadata_file=snakemake.input.metadata,
        output_file=snakemake.output[0],
        is_linear=snakemake.params.linear
    )