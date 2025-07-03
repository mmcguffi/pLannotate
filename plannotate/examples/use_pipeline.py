#!/usr/bin/env python3
"""
Example usage of the pLannotate Snakemake pipeline
"""

from plannotate import annotate

# Example 1: Basic usage (same as before)
def example_basic():
    # Example sequence (portion of pUC19)
    sequence = """
    TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTG
    TAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCT
    """
    sequence = sequence.replace('\n', '').replace(' ', '')
    
    print("Running basic annotation...")
    results = annotate.annotate(sequence)
    print(f"Found {len(results)} annotations")
    
    return results


# Example 2: With custom options
def example_with_options():
    # Example linear sequence
    sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    
    print("\nRunning annotation with options...")
    results = annotate.annotate(
        sequence,
        linear=True,           # Linear sequence
        is_detailed=True,      # Detailed annotation types
        threads=8,             # Use 8 threads
        keep_temp=False        # Don't keep temporary files
    )
    
    return results


# Example 3: Save to GenBank file
def example_save_genbank():
    sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    
    print("\nSaving annotation to GenBank file...")
    annotate.annotate_to_gbk(
        sequence,
        "output.gbk",
        linear=True,
        threads=4
    )
    print("Saved to output.gbk")


# Example 4: Keep intermediate files for debugging
def example_with_output_dir():
    sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    
    print("\nRunning with specified output directory...")
    results = annotate.annotate(
        sequence,
        output_dir="my_results",  # Keep all intermediate files here
        keep_temp=True,            # Don't delete the directory
        threads=4
    )
    
    print("Results saved in my_results/")
    print("Check my_results/searches/ for per-database results")
    
    return results


# Example 5: Using custom database configuration
def example_custom_databases():
    sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    
    print("\nRunning with custom database configuration...")
    results = annotate.annotate(
        sequence,
        yaml_file="custom_databases.yml",  # Path to custom YAML
        threads=4
    )
    
    return results


# Example 6: Direct Snakemake usage
def example_snakemake_direct():
    """
    For advanced users who want to run Snakemake directly
    """
    import yaml
    
    # Create custom config
    config = {
        "input_sequence": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA",
        "linear": True,
        "is_detailed": False,
        "output_dir": "snakemake_results"
    }
    
    # Save config
    with open("my_config.yaml", "w") as f:
        yaml.dump(config, f)
    
    print("\nRunning Snakemake directly...")
    success = annotate.run_pipeline("my_config.yaml", threads=4)
    
    if success:
        print("Pipeline completed successfully")
    else:
        print("Pipeline failed")
    
    return success


if __name__ == "__main__":
    # Run examples
    print("pLannotate Pipeline Examples")
    print("=" * 50)
    
    # Example 1
    results1 = example_basic()
    if not results1.empty:
        print(results1[['Feature', 'Type', 'qstart', 'qend']].head())
    
    # Example 2
    results2 = example_with_options()
    
    # Example 3
    example_save_genbank()
    
    # Example 4
    results4 = example_with_output_dir()
    
    print("\nDone!")