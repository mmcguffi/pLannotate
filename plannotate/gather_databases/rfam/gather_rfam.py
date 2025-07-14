#!/usr/bin/env python3

"""
Script to download the latest Rfam release and capture version information.

"no-press" is an option to skip pressing (indexing) covariance models.
Indexing speeds searches and is not required for use of infernal.

Usage:
    python3 gather_rfam.py --output-dir /path/to/output
    python3 gather_rfam.py --output-dir /path/to/output --keep-archive
    python3 gather_rfam.py --output-dir /path/to/output --no-press
"""

import argparse
import gzip
import subprocess
import sys
from pathlib import Path
from typing import Optional

import requests


def get_rfam_version() -> Optional[str]:
    """Get the latest Rfam release number from the FTP site."""
    try:
        response = requests.get(
            "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/README"
        )
        if response.status_code == 200:
            for line in response.text.split("\n"):
                if line.startswith("Release "):
                    # Extract version from line like "Release 15.0"
                    release_num = line.replace("Release ", "").strip()
                    print(f"Found release: {release_num}")
                    return release_num
        return None
    except Exception as e:
        print(f"Error getting version: {e}")
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


def extract_cm_file(cm_gz_path: Path, output_dir: Path) -> bool:
    """Extract covariance model file from gzipped archive."""
    try:
        cm_path = output_dir / "Rfam.cm"
        print(f"Extracting {cm_gz_path} to {cm_path}")

        with gzip.open(cm_gz_path, "rb") as f_in:
            with open(cm_path, "wb") as f_out:
                f_out.write(f_in.read())

        print("Extracted Rfam.cm")
        return True
    except Exception as e:
        print(f"Error extracting {cm_gz_path}: {e}")
        return False


def press_covariance_models(output_dir: Path) -> bool:
    """Press covariance models using cmpress for Infernal."""
    cm_file = output_dir / "Rfam.cm"
    if not cm_file.exists():
        print(f"Error: {cm_file} not found")
        return False

    try:
        print("Pressing covariance models...")
        result = subprocess.run(
            ["cmpress", "Rfam.cm"], cwd=output_dir, capture_output=True, text=True
        )

        if result.returncode == 0:
            print("Covariance models pressed successfully")
            print(result.stdout)
            return True
        else:
            print(f"Error pressing models (return code {result.returncode}): {result.stderr}")
            if result.stdout:
                print(f"stdout: {result.stdout}")
            return False
    except FileNotFoundError:
        print("Error: cmpress command not found. Please ensure Infernal is installed.")
        return False
    except Exception as e:
        print(f"Error pressing models: {e}")
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Download latest Rfam release")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for downloaded files",
    )
    parser.add_argument(
        "--keep-archive",
        action="store_true",
        help="Keep the downloaded .gz file",
    )
    parser.add_argument(
        "--no-press",
        action="store_true",
        help="Skip pressing (indexing) covariance models",
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get version information
    version = get_rfam_version()
    if not version:
        print("Could not determine Rfam version")
        sys.exit(1)

    # Download the consolidated covariance models file
    rfam_cm_url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    cm_gz_path = output_dir / "Rfam.cm.gz"

    if not download_file(rfam_cm_url, cm_gz_path):
        sys.exit(1)

    # Download the clan information file
    rfam_clan_url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
    clan_path = output_dir / "Rfam.clanin"

    if not download_file(rfam_clan_url, clan_path):
        sys.exit(1)

    # Extract the CM file
    if not extract_cm_file(cm_gz_path, output_dir):
        sys.exit(1)

    # Press covariance models unless skipped
    if not args.no_press:
        if not press_covariance_models(output_dir):
            print("Warning: Failed to press covariance models")

    # Write version info in the correct format (just "release X.Y")
    version_file = output_dir / "version.txt"
    with open(version_file, "w") as f:
        f.write(f"release {version}\n")

    print(f"Version information written to {version_file}")

    # Clean up archive if requested
    if not args.keep_archive:
        cm_gz_path.unlink()
        print("Cleaned up archive file")

    print(f"Rfam database download complete: {version}")


if __name__ == "__main__":
    main()
