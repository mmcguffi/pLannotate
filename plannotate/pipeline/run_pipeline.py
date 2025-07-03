"""
Wrapper script to run pLannotate Snakemake pipeline programmatically
"""
import os
import yaml
import subprocess
import tempfile
import shutil
from pathlib import Path
import pandas as pd


class PlannotatePipeline:
    """Wrapper class for running pLannotate Snakemake pipeline"""
    
    def __init__(self, work_dir=None):
        """Initialize pipeline wrapper"""
        self.pipeline_dir = Path(__file__).parent
        self.work_dir = work_dir or tempfile.mkdtemp(prefix="plannotate_")
        self.config_file = None
    
    def create_config(self, input_sequence, output_dir="results", 
                     yaml_file=None, linear=False, detailed=False, threads=4):
        """Create configuration file for pipeline"""
        if yaml_file is None:
            yaml_file = str(self.pipeline_dir.parent / "data" / "databases.yml")
        
        config = {
            "input_sequence": str(input_sequence),
            "yaml_file": str(yaml_file),
            "work_dir": os.path.join(self.work_dir, "work"),
            "results_dir": str(output_dir),
            "linear": linear,
            "detailed": detailed,
            "threads": threads
        }
        
        # Write config file
        self.config_file = os.path.join(self.work_dir, "config.yaml")
        with open(self.config_file, 'w') as f:
            yaml.dump(config, f)
        
        return self.config_file
    
    def run(self, targets=None, cores=4, dryrun=False, quiet=False):
        """Run the Snakemake pipeline"""
        cmd = [
            "snakemake",
            "-s", str(self.pipeline_dir / "Snakefile"),
            "--configfile", self.config_file,
            "-c", str(cores)
        ]
        
        if targets:
            cmd.extend(targets)
        
        if dryrun:
            cmd.append("-n")
        
        if quiet:
            cmd.append("-q")
        
        # Change to pipeline directory to find scripts
        original_dir = os.getcwd()
        try:
            os.chdir(self.pipeline_dir)
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"Snakemake failed:\n{result.stderr}")
            
            return result
        finally:
            os.chdir(original_dir)
    
    def get_results(self, output_dir):
        """Get results from pipeline run"""
        results = {}
        
        # Find result files
        result_files = Path(output_dir).glob("*")
        for file in result_files:
            if file.suffix == ".csv":
                results["annotations_csv"] = pd.read_csv(file)
            elif file.suffix == ".gbk":
                with open(file) as f:
                    results["genbank"] = f.read()
            elif file.suffix == ".html":
                with open(file) as f:
                    results["bokeh_html"] = f.read()
        
        return results
    
    def cleanup(self):
        """Clean up temporary work directory"""
        if self.work_dir and self.work_dir.startswith(tempfile.gettempdir()):
            shutil.rmtree(self.work_dir, ignore_errors=True)


def run_pipeline(input_sequence, output_dir="results", **kwargs):
    """
    Convenience function to run pLannotate pipeline
    
    Args:
        input_sequence: Path to input FASTA/GenBank file
        output_dir: Directory for output files
        **kwargs: Additional arguments (linear, detailed, cores, etc.)
    
    Returns:
        Dictionary with results
    """
    pipeline = PlannotatePipeline()
    
    try:
        # Create config
        pipeline.create_config(
            input_sequence=input_sequence,
            output_dir=output_dir,
            linear=kwargs.get("linear", False),
            detailed=kwargs.get("detailed", False),
            threads=kwargs.get("cores", 4),
            yaml_file=kwargs.get("yaml_file")
        )
        
        # Run pipeline
        pipeline.run(cores=kwargs.get("cores", 4))
        
        # Get results
        results = pipeline.get_results(output_dir)
        
        return results
        
    finally:
        pipeline.cleanup()


if __name__ == "__main__":
    # Example usage
    import argparse
    
    parser = argparse.ArgumentParser(description="Run pLannotate Snakemake pipeline")
    parser.add_argument("input", help="Input FASTA or GenBank file")
    parser.add_argument("-o", "--output", default="results", help="Output directory")
    parser.add_argument("-c", "--cores", type=int, default=4, help="Number of cores")
    parser.add_argument("--linear", action="store_true", help="Linear DNA mode")
    parser.add_argument("--detailed", action="store_true", help="Detailed search mode")
    
    args = parser.parse_args()
    
    results = run_pipeline(
        args.input,
        args.output,
        cores=args.cores,
        linear=args.linear,
        detailed=args.detailed
    )
    
    print(f"Pipeline completed. Results saved to {args.output}/")