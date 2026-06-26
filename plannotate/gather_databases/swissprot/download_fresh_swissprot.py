#!/usr/bin/env python3
"""
Download fresh Swiss-Prot database and split into manageable chunks

Usage examples:
    # Download and split everything (recommended)
    python download_temp_swissprot.py

    # Only download files, don't split
    python download_temp_swissprot.py --download-only

    # Custom split size (default is 1000 lines per file)
    python download_temp_swissprot.py --split-lines 2000

    # Custom output directory
    python download_fresh_swissprot.py --output-dir ../gathered_data/swissprot_temp

What it does:
    1. Downloads latest uniprot_sprot.dat.gz and uniprot_sprot.fasta.gz
    2. Uncompresses them to ../gathered_data/swissprot/temp_data/
    3. Splits the .dat file into chunks in ../gathered_data/swissprot/temp_split/

After running, description processing is handled by Snakemake workflow.

Output files:
    - ../gathered_data/swissprot/temp_data/uniprot_sprot.fasta (protein sequences)
    - ../gathered_data/swissprot/temp_split/chunk_*.dat (split files for processing)
    - Processed results go to ../gathered_data/swissprot/temp_results/

Requirements:
    - Internet connection
    - Sufficient disk space (~2GB+)
    - process_splits.sh and parse_swissprot_file.py for description extraction
"""

import argparse
import gzip
import os
import shutil
import urllib.request
from pathlib import Path

SWISSPROT_RELEASE_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/reldate.txt"
)


def download_release_version(output_dir):
    """Record the exact Swiss-Prot release represented by the downloaded data."""
    release_path = Path(output_dir) / "reldate.txt"
    urllib.request.urlretrieve(SWISSPROT_RELEASE_URL, release_path)
    release_line = next(
        (
            line.strip()
            for line in release_path.read_text().splitlines()
            if line.startswith("UniProtKB/Swiss-Prot Release")
        ),
        "Unknown Swiss-Prot release",
    )
    version_path = Path(output_dir) / "version.txt"
    version_path.write_text(release_line + "\n")
    return version_path


def download_swissprot(output_dir="gathered_data/swissprot/temp_data"):
    """Download latest Swiss-Prot database from UniProt"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # URLs for Swiss-Prot database
    urls = {
        "uniprot_sprot.dat.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
        "uniprot_sprot.fasta.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    }

    downloaded_files = {}

    for filename, url in urls.items():
        output_path = os.path.join(output_dir, filename)
        uncompressed_path = output_path.replace(".gz", "")

        print(f"Downloading {filename}...")
        try:
            urllib.request.urlretrieve(url, output_path)
            print(f"Downloaded: {output_path}")

            # Uncompress the file
            print(f"Uncompressing {filename}...")
            with gzip.open(output_path, "rb") as f_in:
                with open(uncompressed_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove compressed file to save space
            os.remove(output_path)
            downloaded_files[filename.replace(".gz", "")] = uncompressed_path
            print(f"Uncompressed: {uncompressed_path}")

        except Exception as e:
            print(f"Error downloading {filename}: {e}")

    return downloaded_files


def split_swissprot_file(
    input_file, output_dir="gathered_data/swissprot/temp_split", records_per_file=35000
):
    """Split Swiss-Prot .dat file into smaller chunks on record boundaries (//)"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Splitting {input_file} into chunks of ~{records_per_file} records...")
    return split_by_records_python(input_file, output_dir, records_per_file)


def split_by_records_python(input_file, output_dir, records_per_file=35000):
    """Fallback Python-based splitting by Swiss-Prot records"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open(input_file, "r") as f:
        record_count = 0
        file_count = 0
        current_record = []
        output_file = None

        for line in f:
            current_record.append(line)

            if line.strip() == "//":
                # End of record
                if record_count % records_per_file == 0:
                    # Start new file
                    if output_file:
                        output_file.close()
                    file_count += 1
                    output_file = open(f"{output_dir}/chunk_{file_count:02d}.dat", "w")

                # Write complete record
                if output_file:
                    output_file.writelines(current_record)
                current_record = []
                record_count += 1

        if output_file:
            output_file.close()

    print(f"Created {file_count} split files in {output_dir}")
    return output_dir


def main():
    parser = argparse.ArgumentParser(
        description="Download and process fresh Swiss-Prot database"
    )
    parser.add_argument(
        "--download-only", action="store_true", help="Only download, don't split"
    )
    parser.add_argument(
        "--split-records",
        type=int,
        default=35000,
        help="Records per split file (default: 35K for ~16 chunks)",
    )
    parser.add_argument(
        "--output-dir",
        default="gathered_data/swissprot/temp_data",
        help="Output directory for downloads",
    )

    args = parser.parse_args()

    # Download fresh data
    print("=== Downloading fresh Swiss-Prot database ===")
    downloaded_files = download_swissprot(args.output_dir)
    download_release_version(args.output_dir)

    if not args.download_only and "uniprot_sprot.dat" in downloaded_files:
        # Split the .dat file
        print("\n=== Splitting Swiss-Prot database ===")
        split_dir = split_swissprot_file(
            downloaded_files["uniprot_sprot.dat"], records_per_file=args.split_records
        )

        print("\nNext steps:")
        print("1. Processing handled by Snakemake workflow")
        print(f"2. Split files created in: {split_dir}")
        print(
            f"2. Fresh FASTA sequences are available at: {downloaded_files.get('uniprot_sprot.fasta', 'Not downloaded')}"
        )


if __name__ == "__main__":
    main()
