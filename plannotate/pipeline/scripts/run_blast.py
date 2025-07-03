"""
Run BLAST/Diamond/Infernal searches for pLannotate pipeline
"""
import json
import shlex
import subprocess
import pandas as pd
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from pathlib import Path
import sys

# Add parent directory to path to import infernal parser
sys.path.append(str(Path(__file__).parent.parent.parent))
from infernal import parse_infernal


def run_blast(seq, db_config):
    """Run BLAST/Diamond/Infernal search based on database configuration"""
    task = db_config["method"]
    parameters = db_config["parameters"]
    db_loc = db_config["db_loc"]
    
    query = NamedTemporaryFile(delete=False, suffix=".fasta")
    tmp = NamedTemporaryFile(delete=False)
    
    # Write query sequence
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")
    
    try:
        if task == "blastn":
            flags = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
            cmd = (
                f"blastn -task blastn-short -query {query.name} -out {tmp.name} "
                f'-db {db_loc} {parameters} -outfmt "6 {flags}"'
            )
            result = subprocess.run(
                shlex.split(cmd),
                shell=False,
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"BLAST failed: {result.stderr}")
            
            # Parse BLAST output
            with open(tmp.name, "r") as f:
                align = f.readlines()
            
            if not align:
                return pd.DataFrame()
            
            df = pd.DataFrame([line.split() for line in align], columns=flags.split())
            df = df.apply(pd.to_numeric, errors="ignore")
            
        elif task == "diamond":
            flags = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
            cmd = (
                f"diamond blastx -d {db_loc} -q {query.name} -o {tmp.name} "
                f"{parameters} --outfmt 6 {flags}"
            )
            result = subprocess.run(
                shlex.split(cmd),
                shell=False,
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"Diamond failed: {result.stderr}")
            
            # Parse Diamond output
            with open(tmp.name, "r") as f:
                align = f.readlines()
            
            if not align:
                return pd.DataFrame()
            
            df = pd.DataFrame([line.split() for line in align], columns=flags.split())
            df = df.apply(pd.to_numeric, errors="ignore")
            
            # Diamond-specific processing
            try:
                df["sseqid"] = df["sseqid"].str.split("|", n=2, expand=True)[1]
            except (ValueError, KeyError):
                pass
            df["sframe"] = (df["qstart"] < df["qend"]).astype(int).replace(0, -1)
            df["slen"] = df["slen"] * 3
            df["length"] = abs(df["qend"] - df["qstart"]) + 1
            
        elif task == "infernal":
            flags = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
            cmd = f"cmscan {flags} {parameters} --tblout {tmp.name} --clanin {db_loc} {query.name}"
            result = subprocess.run(
                shlex.split(cmd),
                shell=False,
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"Infernal failed: {result.stderr}")
            
            # Parse Infernal output
            df = parse_infernal(tmp.name)
            df["qlen"] = len(seq)
            
            # Extract DNA sequence
            if not df.empty:
                df["qseq"] = df.apply(
                    lambda x: seq[x["qstart"]:x["qend"] + 1].upper(), axis=1
                )
        else:
            raise ValueError(f"Unknown search method: {task}")
            
    finally:
        # Clean up temp files
        query.close()
        tmp.close()
        Path(query.name).unlink(missing_ok=True)
        Path(tmp.name).unlink(missing_ok=True)
    
    return df


def main(input_fasta, input_metadata, output_file, db_config, db_name):
    """Main function to run BLAST search"""
    # Read sequence
    record = SeqIO.read(input_fasta, "fasta")
    seq = str(record.seq)
    
    # Read metadata
    with open(input_metadata) as f:
        metadata = json.load(f)
    
    # Run search
    df = run_blast(seq, db_config)
    
    # Add database name
    if not df.empty:
        df["db"] = db_name
        df["sseqid"] = df["sseqid"].astype(str)
    
    # Save results
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_fasta=snakemake.input.fasta,
        input_metadata=snakemake.input.metadata,
        output_file=snakemake.output[0],
        db_config=snakemake.params.db_config,
        db_name=snakemake.params.db_name
    )