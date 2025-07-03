"""
Details module for pLannotate pipeline.
Handles retrieving feature descriptions from various databases.
"""

import subprocess
from tempfile import NamedTemporaryFile

import pandas as pd

from .. import resources as rsc
from ..logging_config import get_logger

logger = get_logger(__name__)


def get_feature_details(hits_df, yaml_file_loc):
    """
    Get feature details for a set of hits from a single database.
    
    Args:
        hits_df: DataFrame with hits from a single database
        yaml_file_loc: Path to databases YAML file
        
    Returns:
        DataFrame with feature details merged
    """
    if hits_df.empty:
        return hits_df
    
    # Ensure all hits are from the same database
    db_names = hits_df["db"].unique()
    if len(db_names) != 1:
        raise ValueError("All hits must be from the same database")
    
    db_name = db_names[0]
    databases = rsc.get_yaml(yaml_file_loc)
    database = databases[db_name]
    
    # Get unique sequence IDs
    sseqids = hits_df["sseqid"].tolist()
    sseqids = [s for s in sseqids if s]  # Remove empty values
    
    # Fix problematic sequence IDs (e.g., "pdb|3xHA|" -> "3xHA")
    problem_pattern = r"pdb\|(.*)\|"
    hits_df["sseqid"] = hits_df["sseqid"].str.replace(
        problem_pattern, r"\1", regex=True
    )
    
    # Get feature descriptions
    feat_desc = _retrieve_descriptions(
        database, db_name, sseqids, hits_df
    )
    
    # Merge descriptions with hits
    hits_df = hits_df.merge(
        feat_desc, on="sseqid", how="left", suffixes=("_x", None)
    )
    
    # Drop duplicate columns ending with _x
    hits_df = hits_df[hits_df.columns.drop(list(hits_df.filter(regex="_x")))]
    
    # Add default type if specified
    db_details = database["details"]
    if db_details["default_type"] != "None":
        hits_df["Type"] = db_details["default_type"]
    
    # Remove primer binding sites
    hits_df = hits_df.loc[hits_df["Type"] != "primer_bind"]
    
    # Add priority and adjust for SwissProt if needed
    hits_df["priority"] = database["priority"]
    if "priority_mod" in hits_df.columns:
        hits_df["priority"] = hits_df["priority"] + hits_df["priority_mod"]
        hits_df = hits_df.drop("priority_mod", axis=1)
    
    return hits_df


def _retrieve_descriptions(database, db_name, sseqids, hits_df):
    """Retrieve descriptions from database files or embedded data."""
    db_details = database["details"]
    
    if db_details["location"] == "None":
        # Data is already in the dataframe
        return hits_df[["sseqid", "Feature", "Description"]]
    
    # Determine file location
    if db_details["location"] == "Default":
        details_file_loc = rsc.get_details(db_name) + ".csv"
    else:
        details_file_loc = db_details["location"]
    
    # Handle compressed files
    if db_details["compressed"] is True:
        details_file_loc += ".gz"
        feat_desc = _parse_compressed_details(sseqids, details_file_loc)
    else:
        feat_desc = pd.read_csv(details_file_loc)
    
    # Special handling for SwissProt protein existence levels
    if db_name == "swissprot":
        feat_desc = _add_swissprot_priority(feat_desc)
    
    return feat_desc


def _parse_compressed_details(sseqids, gz_loc):
    """Parse compressed detail files using ripgrep."""
    hits_pattern = "|".join(sseqids)
    output = NamedTemporaryFile(suffix="csv")
    
    subprocess.call(
        f'rg -z "{hits_pattern}" {gz_loc} > {output.name}', 
        shell=True
    )
    
    gz_details = pd.read_csv(
        output.name, 
        header=None, 
        names=["sseqid", "Feature", "Description"]
    )
    
    output.close()
    return gz_details


def _add_swissprot_priority(feat_desc):
    """Add priority modifier based on SwissProt protein existence level."""
    # Find existence level in description
    level_pos = feat_desc["Description"].str.find("existence level") + 16
    feat_desc["s"] = level_pos
    feat_desc["e"] = level_pos + 1
    
    def calc_priority_mod(desc, start, end):
        # Default priority if no existence level found
        if start == 15 and end == 16:
            return 0
        else:
            return int(desc[start:end]) - 1
    
    feat_desc["priority_mod"] = [
        calc_priority_mod(d, s, e)
        for d, s, e in zip(
            feat_desc["Description"], 
            feat_desc["s"], 
            feat_desc["e"]
        )
    ]
    
    return feat_desc.drop(columns=["s", "e"])