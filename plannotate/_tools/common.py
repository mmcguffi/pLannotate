"""Execution and file helpers shared by external annotation tools."""

import logging
import shlex
import subprocess
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def run_command(command: str | Sequence[str], tool: str) -> None:
    """Run an external tool and include its diagnostic output in failures."""
    arguments = shlex.split(command) if isinstance(command, str) else list(command)
    logger.debug("Executing %s command: %s", tool, shlex.join(arguments))
    try:
        result = subprocess.run(arguments, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"{tool} is not installed or available on PATH") from exc
    if result.returncode != 0:
        diagnostic = result.stderr.strip() or result.stdout.strip() or "no output"
        raise RuntimeError(
            f"{tool} failed with exit code {result.returncode}: {diagnostic}"
        )
    logger.debug("%s command completed successfully", tool)


@contextmanager
def temporary_files(sequence: str) -> Iterator[tuple[str, str]]:
    """Provide temporary FASTA input and tabular output paths."""
    with TemporaryDirectory(prefix="plannotate-tool-") as temp_dir:
        query_path = Path(temp_dir) / "query.fasta"
        output_path = Path(temp_dir) / "results.tsv"
        SeqIO.write(SeqRecord(Seq(sequence), id="query"), query_path, "fasta")
        output_path.touch()
        yield str(query_path), str(output_path)


def read_table(path: str, columns: str) -> pd.DataFrame:
    """Read whitespace-separated output and infer numeric columns."""
    rows = [line.split() for line in Path(path).read_text().splitlines()]
    dataframe = pd.DataFrame(rows, columns=columns.split())
    for column in dataframe.columns:
        try:
            dataframe[column] = pd.to_numeric(dataframe[column])
        except (TypeError, ValueError):
            pass
    return dataframe
