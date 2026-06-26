"""Programmatic API for rebuilding pLannotate's database bundle."""

import shutil
import subprocess
import sys
from pathlib import Path

ARCHIVE_NAME = "plannotate-databases-v2.tar.gz"


def build_databases(
    output_directory: str | Path = ".",
    *,
    cores: int = 1,
    dry_run: bool = False,
) -> Path:
    """Run the packaged database-build workflow and return its archive path.

    This requires the ``databases`` optional dependencies and the BLAST, DIAMOND,
    and Infernal executables used by the selected targets.
    """
    if cores < 1:
        raise ValueError("cores must be at least 1")

    snakemake = shutil.which("snakemake")
    if snakemake is None:
        raise RuntimeError(
            "Database building requires Snakemake; install pLannotate[databases]"
        )

    output_directory = Path(output_directory).expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    workflow = Path(__file__).with_name("gather.smk")
    command = [
        snakemake,
        "--snakefile",
        str(workflow),
        "--directory",
        str(output_directory),
        "--cores",
        str(cores),
        "--config",
        f"python_executable={sys.executable}",
    ]
    if dry_run:
        command.append("--dry-run")

    result = subprocess.run(command)
    if result.returncode != 0:
        raise RuntimeError(f"Database build failed with exit code {result.returncode}")
    return output_directory / ARCHIVE_NAME
