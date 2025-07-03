"""
Add feature details and descriptions to BLAST hits
"""
import pandas as pd
import subprocess
import yaml
from tempfile import NamedTemporaryFile
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
import resources as rsc


def parse_gz(sseqids, gz_loc):
    """Parse compressed detail files"""
    hits = "|".join(sseqids)
    output = NamedTemporaryFile(suffix="csv", delete=False)
    try:
        subprocess.call(f'rg -z "{hits}" {gz_loc} > {output.name}', shell=True)
        gz_details = pd.read_csv(
            output.name, header=None, names=["sseqid", "Feature", "Description"]
        )
    finally:
        output.close()
        Path(output.name).unlink(missing_ok=True)
    return gz_details


def get_feature_details(hits_df, db_config, db_name, yaml_file):
    """Get feature details for hits from database configuration"""
    if hits_df.empty:
        return hits_df
    
    # Get unique sequence IDs
    sseqids = hits_df["sseqid"].unique().tolist()
    sseqids = [s for s in sseqids if s]  # Remove empty values
    
    # Fix problematic IDs (e.g., "pdb|3xHA|" -> "3xHA")
    problem_pattern = r"pdb\|(.*)\|"
    hits_df["sseqid"] = hits_df["sseqid"].str.replace(problem_pattern, r"\1", regex=True)
    
    db_details = db_config["details"]
    
    if db_details["location"] == "None":
        # Details already in dataframe (e.g., from Infernal)
        feat_desc = hits_df[["sseqid", "Feature", "Description"]].drop_duplicates()
    else:
        # Load details from file
        if db_details["location"] == "Default":
            details_file = rsc.get_details(db_name) + ".csv"
        else:
            details_file = db_details["location"]
        
        # Handle compressed files
        if db_details.get("compressed", False):
            details_file += ".gz"
            feat_desc = parse_gz(sseqids, details_file)
        else:
            feat_desc = pd.read_csv(details_file)
            # Filter to only needed IDs
            feat_desc = feat_desc[feat_desc["sseqid"].isin(sseqids)]
        
        # Special handling for SwissProt protein existence levels
        if db_name == "swissprot":
            level_pos = feat_desc["Description"].str.find("existence level") + 16
            feat_desc["priority_mod"] = feat_desc.apply(
                lambda row: 0 if level_pos[row.name] == 15 else int(
                    row["Description"][level_pos[row.name]:level_pos[row.name] + 1]
                ) - 1,
                axis=1
            )
    
    # Add default type if specified
    if db_details.get("default_type", "None") != "None":
        feat_desc["Type"] = db_details["default_type"]
    
    # Merge details with hits
    result = hits_df.merge(feat_desc, on="sseqid", how="left", suffixes=("_x", None))
    
    # Drop duplicate columns
    result = result[result.columns.drop(list(result.filter(regex="_x")))]
    
    # Remove primer binding sites
    result = result[result["Type"] != "primer_bind"]
    
    # Add priority
    result["priority"] = db_config["priority"]
    if "priority_mod" in result.columns:
        result["priority"] = result["priority"] + result["priority_mod"]
        result = result.drop("priority_mod", axis=1)
    
    return result


def main(input_file, output_file, db_config, db_name, yaml_file):
    """Main function to add feature details"""
    # Read BLAST hits
    hits_df = pd.read_csv(input_file, sep="\t")
    
    if hits_df.empty:
        hits_df.to_csv(output_file, sep="\t", index=False)
        return
    
    # Add feature details
    detailed_df = get_feature_details(hits_df, db_config, db_name, yaml_file)
    
    # Save results
    detailed_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_file=snakemake.input.blast_hits,
        output_file=snakemake.output[0],
        db_config=snakemake.params.db_config,
        db_name=snakemake.params.db_name,
        yaml_file=snakemake.params.yaml_file
    )