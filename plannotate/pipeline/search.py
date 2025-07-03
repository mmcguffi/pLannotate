"""
Search module for pLannotate pipeline.
Handles BLAST, DIAMOND, and Infernal searches against individual databases.
"""

import shlex
import subprocess
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..infernal import parse_infernal
from ..logging_config import get_logger

logger = get_logger(__name__)


def search_database(seq, db, db_name):
    """
    Run search against a single database using appropriate method.
    
    Args:
        seq: DNA sequence string
        db: Database configuration dictionary
        db_name: Name of the database
        
    Returns:
        DataFrame with search results
    """
    task = db["method"]
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    
    logger.info(f"Searching {db_name} using {task}")
    
    query = NamedTemporaryFile(delete=False)
    tmp = NamedTemporaryFile(delete=False)
    
    try:
        # Write query sequence
        SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")
        
        if task == "blastn":
            results = _run_blastn(query.name, tmp.name, db_loc, parameters)
        elif task == "diamond":
            results = _run_diamond(query.name, tmp.name, db_loc, parameters, seq)
        elif task == "infernal":
            results = _run_infernal(query.name, tmp.name, db_loc, parameters, seq)
        else:
            raise ValueError(f"Unknown search method: {task}")
            
        if not results.empty:
            results["db"] = db_name
            results["sseqid"] = results["sseqid"].astype(str)
            logger.info(f"Found {len(results)} hits in {db_name}")
        else:
            logger.info(f"No hits found in {db_name}")
            
        return results
        
    except Exception as e:
        logger.error(f"Error searching {db_name}: {e}")
        return pd.DataFrame()
        
    finally:
        query.close()
        tmp.close()


def _run_blastn(query_file, output_file, db_loc, parameters):
    """Run BLASTN search."""
    flags = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"blastn -task blastn-short -query {query_file} -out {output_file} "
        f'-db {db_loc} {parameters} -outfmt "6 {flags}"'
    )
    
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        logger.warning(f"BLASTN error: {result.stderr}")
    
    with open(output_file, "r") as f:
        align = f.readlines()
    
    if not align:
        return pd.DataFrame()
        
    df = pd.DataFrame([line.split() for line in align], columns=flags.split())
    
    # Convert numeric columns
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass
            
    return df


def _run_diamond(query_file, output_file, db_loc, parameters, seq):
    """Run DIAMOND BLASTX search."""
    flags = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"diamond blastx -d {db_loc} -q {query_file} -o {output_file} "
        f"{parameters} --outfmt 6 {flags}"
    )
    
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        logger.warning(f"DIAMOND error: {result.stderr}")
    
    with open(output_file, "r") as f:
        align = f.readlines()
    
    if not align:
        return pd.DataFrame()
        
    df = pd.DataFrame([line.split() for line in align], columns=flags.split())
    
    # Convert numeric columns
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass
    
    # Process DIAMOND-specific fields
    if "sseqid" in df.columns and not df["sseqid"].empty:
        df["sseqid"] = df["sseqid"].astype(str).str.split("|", n=2).str.get(1)
    df["sframe"] = (df["qstart"] < df["qend"]).astype(int).replace(0, -1)
    df["slen"] = df["slen"] * 3
    df["length"] = abs(df["qend"] - df["qstart"]) + 1
    
    return df


def _run_infernal(query_file, output_file, db_loc, parameters, seq):
    """Run Infernal cmscan search."""
    flags = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
    cmd = f"cmscan {flags} {parameters} --tblout {output_file} --clanin {db_loc} {query_file}"
    
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        logger.warning(f"Infernal error: {result.stderr}")
    
    df = parse_infernal(output_file)
    
    if not df.empty:
        df["qlen"] = len(seq)
        df["qseq"] = df.apply(
            lambda x: seq[x["qstart"]:x["qend"] + 1].upper(), axis=1
        )
    
    return df