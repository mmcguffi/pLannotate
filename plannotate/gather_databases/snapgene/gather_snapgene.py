#!/usr/bin/env python3
"""
Script to download SnapGene database from plannotate GitHub release 1.2.0.

Usage:
    python3 gather_snapgene.py --output-dir /path/to/output
    python3 gather_snapgene.py --output-dir /path/to/output --keep-archive
"""

import argparse
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Optional

import requests


def get_release_info() -> Optional[dict]:
    """Get release information from GitHub API."""
    try:
        response = requests.get(
            "https://api.github.com/repos/mmcguffi/pLannotate/releases"
        )
        if response.status_code == 200:
            releases = response.json()
            for release in releases:
                if release["tag_name"] == "v1.2.0":
                    return release
        return None
    except Exception as e:
        print(f"Error getting release info: {e}")
        return None


def download_file(url: str, output_path: Path) -> bool:
    """Download a file from URL to output path."""
    try:
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"Downloaded to {output_path}")
        return True
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False


def extract_snapgene_files(tar_path: Path, output_dir: Path) -> bool:
    """Extract SnapGene database files from tar archive."""
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            # Extract snapgene files (they should be in the BLAST_dbs directory)
            for member in tar.getmembers():
                if "snapgene" in member.name.lower():
                    # Flatten the path to just the filename
                    filename = member.name.split("/")[-1]
                    if filename:  # Skip directories
                        member.name = filename
                        tar.extract(member, output_dir)
                        print(f"Extracted {filename}")
        return True
    except Exception as e:
        print(f"Error extracting {tar_path}: {e}")
        return False


def check_blast_db(output_dir: Path) -> bool:
    """Check if BLAST database files are present."""
    # Check for typical BLAST database files
    blast_extensions = [".nhr", ".nin", ".nsq"]
    snapgene_base = None

    for file in output_dir.glob("snapgene.*"):
        if any(file.name.endswith(ext) for ext in blast_extensions):
            snapgene_base = file.name.split(".")[0]
            break

    if snapgene_base:
        print(f"Found BLAST database files for {snapgene_base}")
        return True

    # If no pre-built database, try to make one from FASTA
    fasta_file = None
    for file in output_dir.glob("*.fasta"):
        if "snapgene" in file.name.lower():
            fasta_file = file
            break

    if not fasta_file:
        print("No SnapGene FASTA file found")
        return False

    try:
        print("Creating BLAST database from FASTA...")
        result = subprocess.run(
            ["makeblastdb", "-in", str(fasta_file), "-dbtype", "nucl"],
            cwd=output_dir,
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            print("BLAST database created successfully")
            return True
        else:
            print(f"Error creating BLAST database: {result.stderr}")
            return False
    except FileNotFoundError:
        print(
            "Error: makeblastdb command not found. Please ensure BLAST+ is installed."
        )
        return False
    except Exception as e:
        print(f"Error creating BLAST database: {e}")
        return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download SnapGene database from plannotate release 1.2.0"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for downloaded files",
    )
    parser.add_argument(
        "--keep-archive",
        action="store_true",
        help="Keep the downloaded tar.gz file",
    )
    parser.add_argument(
        "--no-makeblastdb",
        action="store_true",
        help="Skip creating BLAST database",
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get release information
    release_info = get_release_info()
    if not release_info:
        print("Could not get release information")
        sys.exit(1)

    # Find the tar.gz asset
    asset_url = None
    for asset in release_info.get("assets", []):
        if asset["name"].endswith(".tar.gz"):
            asset_url = asset["browser_download_url"]
            break

    if not asset_url:
        print("Could not find tar.gz asset in release 1.2.0")
        sys.exit(1)

    # Download the release
    tar_path = output_dir / "plannotate-1.2.0.tar.gz"
    if not download_file(asset_url, tar_path):
        sys.exit(1)

    # Extract SnapGene files
    if not extract_snapgene_files(tar_path, output_dir):
        sys.exit(1)

    # Check/create BLAST database unless skipped
    if not args.no_makeblastdb:
        if not check_blast_db(output_dir):
            print("Warning: BLAST database not found or could not be created")

    # Write version info
    version_file = output_dir / "version.txt"
    with open(version_file, "w") as f:
        f.write("SnapGene from plannotate release 1.2.0\n")

    print(f"Version information written to {version_file}")

    # Clean up archive if requested
    if not args.keep_archive:
        tar_path.unlink()
        print("Cleaned up archive file")

    print("SnapGene database download complete")


if __name__ == "__main__":
    main()
