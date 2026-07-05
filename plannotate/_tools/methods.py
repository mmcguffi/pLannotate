"""Annotation-method registry and external-tool dispatch.

To add a new annotation tool (e.g. an ORF finder or IS-element detector):

1. Implement ``search(sequence, config, threads) -> pandas.DataFrame`` in a new
   ``_tools`` module that returns the columns in ``_schema.ADAPTER_COLUMNS``.
2. Add one entry to ``METHODS`` below. That is the only registration step:
   ``Method(search_fn, database_directory)`` for a tool with a packaged database,
   or ``database_directory=None`` for a database-free detector. Pass a custom
   ``build_paths`` only when the tool needs more than a single database file.
"""

import logging
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd

from . import blast, diamond, infernal

logger = logging.getLogger(__name__)
Search = Callable[[str | Mapping[str, str], dict[str, Any], int], pd.DataFrame]
PathBuilder = Callable[[str, Path], dict[str, str]]


def _single_database_path(database_name: str, directory: Path) -> dict[str, str]:
    return {"db_loc": str(directory / database_name)}


class Method(NamedTuple):
    """How to run one annotation method and locate its database, if any."""

    search: Search
    # packaged database subdirectory, or None for database-free detectors
    database_directory: str | None
    # builds tool-specific database paths from a base location
    build_paths: PathBuilder = _single_database_path


METHODS: dict[str, Method] = {
    "blast": Method(blast.search, "BLAST_dbs"),
    "blastn": Method(blast.search, "BLAST_dbs"),
    "diamond": Method(diamond.search, "diamond_dbs"),
    "infernal": Method(infernal.search, "infernal_dbs", infernal.database_paths),
}


def _method(name: str) -> Method:
    normalized_name = name.lower()
    try:
        return METHODS[normalized_name]
    except KeyError:
        raise ValueError(f"Unsupported database method: {normalized_name}") from None


def database_directory(name: str) -> str | None:
    """Return the packaged database directory for a method, if it needs one."""
    return _method(name).database_directory


def database_paths(name: str, database_name: str, directory: Path) -> dict[str, str]:
    """Build tool-specific database paths for a configured source."""
    return _method(name).build_paths(database_name, directory)


def run(
    sequence: str | Mapping[str, str], config: dict[str, Any], threads: int = 1
) -> pd.DataFrame:
    """Run the tool selected by an annotation-source configuration."""
    name = str(config["method"])
    logger.debug("Dispatching annotation method=%s threads=%d", name.lower(), threads)
    return _method(name).search(sequence, config, threads)
