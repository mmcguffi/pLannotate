"""Configured annotation methods and external-tool dispatch."""

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd

from . import blast, diamond, infernal

logger = logging.getLogger(__name__)
Search = Callable[[str, dict[str, Any], int], pd.DataFrame]
PathBuilder = Callable[[str, Path], dict[str, str]]

SEARCHERS: dict[str, Search] = {
    "blast": blast.search,
    "blastn": blast.search,
    "diamond": diamond.search,
    "infernal": infernal.search,
}

DATABASE_DIRECTORIES: dict[str, str | None] = {
    "blast": "BLAST_dbs",
    "blastn": "BLAST_dbs",
    "diamond": "diamond_dbs",
    "infernal": "infernal_dbs",
}

PATH_BUILDERS: dict[str, PathBuilder] = {
    "infernal": infernal.database_paths,
}


def _method_name(name: str) -> str:
    normalized_name = name.lower()
    if normalized_name not in SEARCHERS:
        raise ValueError(f"Unsupported database method: {normalized_name}")
    return normalized_name


def database_directory(name: str) -> str | None:
    """Return the packaged database directory for a method, if it needs one."""
    return DATABASE_DIRECTORIES[_method_name(name)]


def database_paths(name: str, database_name: str, directory: Path) -> dict[str, str]:
    """Build tool-specific database paths for a configured source."""
    method = _method_name(name)
    path_builder = PATH_BUILDERS.get(method, _single_database_path)
    return path_builder(database_name, directory)


def _single_database_path(database_name: str, directory: Path) -> dict[str, str]:
    return {"db_loc": str(directory / database_name)}


def run(sequence: str, config: dict[str, Any], threads: int = 1) -> pd.DataFrame:
    """Run the tool selected by an annotation-source configuration."""
    method = _method_name(str(config["method"]))
    logger.debug("Dispatching annotation method=%s threads=%d", method, threads)
    search = SEARCHERS[method]
    return search(sequence, config, threads)
