"""
Example integration of Snakemake pipeline with existing pLannotate interface

This shows how the refactored pipeline could be integrated as an alternative
backend while maintaining compatibility with the existing API.
"""
import pandas as pd
from pathlib import Path
from .run_pipeline import run_pipeline


def annotate_with_pipeline(inSeq, yaml_file=None, linear=False, is_detailed=False, cores=4):
    """
    Drop-in replacement for the original annotate() function using Snakemake pipeline
    
    This function provides the same interface as the original annotate() but uses
    the Snakemake pipeline for parallelized execution.
    
    Args:
        inSeq: DNA sequence string
        yaml_file: Path to database configuration YAML
        linear: Whether the DNA is linear (False = circular)
        is_detailed: Whether to use detailed search mode
        cores: Number of cores to use for parallel execution
    
    Returns:
        pandas.DataFrame with same columns as original annotate()
    """
    import tempfile
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_input:
        record = SeqRecord(Seq(inSeq), id="temp_sequence", description="")
        SeqIO.write(record, tmp_input, "fasta")
        input_file = tmp_input.name
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        try:
            # Run pipeline
            results = run_pipeline(
                input_sequence=input_file,
                output_dir=output_dir,
                yaml_file=yaml_file,
                linear=linear,
                detailed=is_detailed,
                cores=cores
            )
            
            # Get annotations DataFrame
            if "annotations_csv" in results:
                # Read the internal TSV format for compatibility
                internal_file = Path(output_dir).parent / "work" / "final" / "temp_sequence_annotations.tsv"
                if internal_file.exists():
                    df = pd.read_csv(internal_file, sep="\t")
                else:
                    # Fall back to CSV and convert
                    df = results["annotations_csv"]
                    # Would need column mapping here for full compatibility
                
                return df
            else:
                # Return empty DataFrame with expected columns
                from .. import resources as rsc
                return pd.DataFrame(columns=rsc.DF_COLS)
                
        finally:
            # Clean up temp input file
            Path(input_file).unlink(missing_ok=True)


def integrate_with_cli():
    """
    Example of how to add pipeline option to existing CLI
    
    This could be added to pLannotate.py
    """
    # This would be added to the existing click command group
    # @main.command("pipeline")
    # @click.option("--input", "-i", required=True, help="Input FASTA/GenBank file")
    # @click.option("--output", "-o", default="./", help="Output directory")
    # @click.option("--cores", "-c", default=4, help="Number of cores for parallel execution")
    # @click.option("--linear", "-l", is_flag=True, help="Linear DNA mode")
    # @click.option("--detailed", "-d", is_flag=True, help="Detailed search mode")
    # def main_pipeline(input, output, cores, linear, detailed):
    #     """Run pLannotate using parallelized Snakemake pipeline"""
    #     results = run_pipeline(
    #         input_sequence=input,
    #         output_dir=output,
    #         cores=cores,
    #         linear=linear,
    #         detailed=detailed
    #     )
    #     click.echo(f"Annotation complete. Results saved to {output}/")
    pass


def integrate_with_streamlit():
    """
    Example of how to use pipeline in Streamlit app
    
    This could replace the existing annotate() call in the Streamlit interface
    """
    # In the Streamlit app, replace:
    #   recordDf = annotate(inSeq, st.session_state.local_db_loc, linear)
    # 
    # With:
    #   if st.checkbox("Use parallel pipeline (faster for large sequences)"):
    #       recordDf = annotate_with_pipeline(
    #           inSeq, 
    #           st.session_state.local_db_loc, 
    #           linear,
    #           cores=st.slider("CPU cores", 1, 8, 4)
    #       )
    #   else:
    #       recordDf = annotate(inSeq, st.session_state.local_db_loc, linear)
    pass


def benchmark_comparison():
    """
    Example code to benchmark original vs pipeline implementation
    """
    import time
    from ..annotate import annotate as original_annotate
    
    # Example sequence (would use real plasmid sequence)
    test_sequence = "ATGC" * 1000
    
    # Time original implementation
    start = time.time()
    original_result = original_annotate(test_sequence)
    original_time = time.time() - start
    
    # Time pipeline implementation
    start = time.time()
    pipeline_result = annotate_with_pipeline(test_sequence, cores=4)
    pipeline_time = time.time() - start
    
    print(f"Original implementation: {original_time:.2f} seconds")
    print(f"Pipeline implementation: {pipeline_time:.2f} seconds")
    print(f"Speedup: {original_time / pipeline_time:.2f}x")
    
    # Verify results are equivalent
    # (would need more sophisticated comparison in practice)
    assert len(original_result) == len(pipeline_result)


if __name__ == "__main__":
    # Example usage
    sequence = "ATGCGATCGTAGCTAGCTAGCTAGCTAGCTAGC"
    
    # Use pipeline version
    result = annotate_with_pipeline(sequence, linear=False, cores=4)
    print(f"Found {len(result)} annotations")
    print(result[["Feature", "Type", "qstart", "qend"]].head())