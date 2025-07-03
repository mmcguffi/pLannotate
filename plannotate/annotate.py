"""
Wrapper for pLannotate Snakemake pipeline.
Provides the same interface as the original annotate() function.
"""

import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd
import snakemake
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .logging_config import get_logger

logger = get_logger(__name__)


def annotate(in_seq, yaml_file=None, linear=False, is_detailed=False, 
             output_dir=None, threads=4, keep_temp=False):
    """
    Annotate a plasmid sequence using the Snakemake pipeline.
    
    Args:
        in_seq: Input DNA sequence string
        yaml_file: Path to custom databases YAML (default: use bundled)
        linear: Whether sequence is linear (default: False)
        is_detailed: Use detailed annotation types (default: False)
        output_dir: Directory for output files (default: temporary)
        threads: Number of threads for parallel processing
        keep_temp: Keep temporary files (default: False)
        
    Returns:
        DataFrame with annotations
    """
    # Use default yaml file if not provided
    if yaml_file is None:
        yaml_file = rsc.get_yaml_path()
    
    # Validate sequence
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    try:
        SeqIO.write(
            SeqRecord(Seq(in_seq), name="pLannotate", 
                     annotations={"molecule_type": "DNA"}),
            temp_fasta.name, "fasta"
        )
        record = list(SeqIO.parse(temp_fasta.name, "fasta"))[0]
    finally:
        temp_fasta.close()
        os.unlink(temp_fasta.name)
    
    # Create temporary directory if output_dir not specified
    temp_dir = None
    if output_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="plannotate_")
        output_dir = temp_dir
    
    try:
        # Get path to Snakefile
        snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
        
        # Create config for this run
        config = {
            "input_sequence": str(record.seq),
            "linear": linear,
            "is_detailed": is_detailed,
            "yaml_file": yaml_file,
            "output_dir": output_dir
        }
        
        # Run Snakemake pipeline
        logger.info("Starting pLannotate pipeline...")
        success = snakemake.snakemake(
            snakefile=snakefile_path,
            config=config,
            cores=threads,
            quiet=True,
            log_handler=logger.handlers
        )
        
        if not success:
            raise RuntimeError("Snakemake pipeline failed")
        
        # Read final results
        final_annotations = os.path.join(output_dir, "final_annotations.csv")
        if os.path.exists(final_annotations):
            df = pd.read_csv(final_annotations)
        else:
            df = pd.DataFrame(columns=rsc.DF_COLS)
        
        return df
        
    finally:
        # Clean up temporary directory if created
        if temp_dir and not keep_temp:
            shutil.rmtree(temp_dir)


def annotate_to_gbk(in_seq, output_gbk, yaml_file=None, linear=False, 
                    is_detailed=False, threads=4):
    """
    Annotate a plasmid and save as GenBank file.
    
    Args:
        in_seq: Input DNA sequence string
        output_gbk: Path for output GenBank file
        yaml_file: Path to custom databases YAML (default: use bundled)
        linear: Whether sequence is linear (default: False)
        is_detailed: Use detailed annotation types (default: False)
        threads: Number of threads for parallel processing
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        # Run annotation
        df = annotate(
            in_seq, 
            yaml_file=yaml_file,
            linear=linear,
            is_detailed=is_detailed,
            output_dir=temp_dir,
            threads=threads,
            keep_temp=True
        )
        
        # Copy GenBank file to desired location
        temp_gbk = os.path.join(temp_dir, "final_annotations.gbk")
        if os.path.exists(temp_gbk):
            shutil.copy(temp_gbk, output_gbk)
        else:
            # Create GenBank file if pipeline didn't
            gbk_content = rsc.get_gbk(df, in_seq, is_linear=linear)
            with open(output_gbk, 'w') as f:
                f.write(gbk_content)


def run_pipeline(config_file=None, threads=4):
    """
    Run the Snakemake pipeline with a config file.
    
    Args:
        config_file: Path to config YAML file
        threads: Number of threads for parallel processing
        
    Returns:
        Boolean indicating success
    """
    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    
    if config_file is None:
        config_file = os.path.join(os.path.dirname(__file__), "config.yaml")
    
    return snakemake.snakemake(
        snakefile=snakefile_path,
        configfiles=[config_file],
        cores=threads
    )