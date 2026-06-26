"""Feature-description lookup in SQLite databases."""

import logging
import re
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)
DESCRIPTION_COLUMNS = ["sseqid", "name", "type", "blurb"]


def _validated_table_name(database_name: str) -> str:
    """Return a safe SQLite identifier for a configured database name."""
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", database_name):
        raise ValueError(f"Invalid database name: {database_name!r}")
    return database_name


def _description_query(connection: sqlite3.Connection, table_name: str) -> str:
    """Build a normalized description query for current and legacy schemas."""
    columns = {
        row[1] for row in connection.execute(f'PRAGMA table_info("{table_name}")')
    }
    if not columns:
        raise ValueError(f"Description table {table_name!r} does not exist")
    if {"sseqid", "name", "type", "blurb"} <= columns:
        return f'SELECT sseqid, name, type, blurb FROM "{table_name}"'
    if {"sseqid", "Feature", "Type", "Description"} <= columns:
        return (
            f"SELECT sseqid, Feature AS name, Type AS type, "
            f'Description AS blurb FROM "{table_name}"'
        )
    if {"sseqid", "name", "blurb"} <= columns:
        return (
            f"SELECT sseqid, name, 'misc_feature' AS type, blurb FROM \"{table_name}\""
        )
    if {"sseqid", "Feature", "Description"} <= columns:
        return (
            f"SELECT sseqid, Feature AS name, 'misc_feature' AS type, "
            f'Description AS blurb FROM "{table_name}"'
        )
    raise ValueError(
        f"Description table {table_name!r} has unsupported columns: {sorted(columns)}"
    )


def _resolve_database_config(
    database_name: str, db_config: dict[str, Any] | None
) -> dict[str, Any]:
    if db_config is not None:
        return db_config

    from . import _package_data

    return _package_data.get_yaml(_package_data.get_yaml_path()).get(database_name, {})


def get_descriptions_db_path(
    database_name: str,
    db_config: dict[str, Any] | None = None,
) -> Path:
    """Get the path to the SQLite descriptions database for a given database."""
    config = _resolve_database_config(database_name, db_config)
    if not config:
        raise ValueError(f"Database '{database_name}' not found in configuration")
    details_location = config.get("details", {}).get("location")
    if details_location not in (None, "None", "Default"):
        details_path = Path(str(details_location))
        return (
            details_path
            if details_path.suffix == ".db"
            else details_path / "descriptions.db"
        )

    database_path = config.get("db_loc")
    if not isinstance(database_path, str) or not database_path:
        raise ValueError(f"Database {database_name!r} has no resolved path")
    return Path(database_path).parent / "descriptions.db"


def load_descriptions_from_sqlite(
    database_name: str,
    sseqids: set[str] | None = None,
    db_config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """
    Load feature descriptions from SQLite database.

    Args:
        database_name: Name of the database (fpbase, swissprot, snapgene)
        sseqids: Optional set of sequence IDs to filter results
        db_config: Optional database configuration dict

    Returns:
        DataFrame with columns: sseqid, name, type, blurb
    """
    db_path = get_descriptions_db_path(database_name, db_config)

    if not db_path.is_file():
        raise FileNotFoundError(f"SQLite description database not found: {db_path}")
    if sseqids is not None and not sseqids:
        return pd.DataFrame(columns=DESCRIPTION_COLUMNS)

    table_name = _validated_table_name(database_name)
    with closing(sqlite3.connect(db_path)) as connection:
        query = _description_query(connection, table_name)
        params: tuple[str, ...] = ()
        if sseqids:
            params = tuple(sorted(sseqids))
            placeholders = ", ".join("?" for _ in params)
            query = f"{query} WHERE sseqid IN ({placeholders})"
        descriptions = pd.read_sql_query(query, connection, params=params)

    logger.debug(
        "Loaded %d descriptions from %s SQLite database",
        len(descriptions),
        database_name,
    )
    return descriptions


def query_description_by_sseqid(
    database_name: str,
    sseqid: str,
    db_config: dict[str, Any] | None = None,
) -> dict[str, str] | None:
    """
    Query a single description by sseqid.

    Args:
        database_name: Name of the database
        sseqid: Sequence ID to query
        db_config: Optional database configuration dict

    Returns:
        Dictionary with name, type, and blurb, or None if not found
    """
    db_path = get_descriptions_db_path(database_name, db_config)

    if not db_path.exists():
        return None

    descriptions = load_descriptions_from_sqlite(
        database_name,
        {sseqid},
        db_config,
    )
    if descriptions.empty:
        return None
    row = descriptions.iloc[0]
    return {column: str(row[column]) for column in ("name", "type", "blurb")}


def check_sqlite_database(
    database_name: str,
    db_config: dict[str, Any] | None = None,
) -> bool:
    """Check if SQLite database exists and is accessible."""
    try:
        db_path = get_descriptions_db_path(database_name, db_config)
        if not db_path.exists():
            return False

        table_name = _validated_table_name(database_name)
        with closing(sqlite3.connect(db_path)) as connection:
            result = connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (table_name,),
            ).fetchone()

        return result is not None

    except (OSError, sqlite3.Error, TypeError, ValueError):
        return False
