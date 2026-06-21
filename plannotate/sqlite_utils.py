"""
Utilities for working with SQLite description databases.

Replaces the CSV-based description system with SQLite for better performance
and to eliminate dependencies on ripgrep for compressed file searching.
"""

import sqlite3
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd


def is_sqlite_enabled(database_name: str, db_config: Optional[Dict] = None) -> bool:
    """Check if SQLite is enabled for a given database."""
    if db_config is None:
        # Import here to avoid circular imports
        from . import resources as rsc

        yaml_config = rsc.get_yaml(rsc.get_yaml_path())
        db_config = yaml_config.get(database_name, {})

    if db_config is None:
        return False

    return db_config.get("details", {}).get("sqlite", False)


def get_descriptions_db_path(
    database_name: str,
    data_dir: Path,
    db_config: Optional[Dict] = None,
) -> Path:
    """Get the path to the SQLite descriptions database for a given database."""
    if db_config is None:
        # Import here to avoid circular imports
        from . import resources as rsc

        yaml_config = rsc.get_yaml(rsc.get_yaml_path())
        db_config = yaml_config.get(database_name, {})

    if db_config is None:
        raise ValueError(f"Database '{database_name}' not found in configuration")

    # Get method from database configuration to determine subdirectory
    method = db_config.get("method", "").lower()

    if method == "diamond":
        return data_dir / "diamond_dbs" / "descriptions.db"
    elif method in ["blastn", "blast"]:
        return data_dir / "BLAST_dbs" / "descriptions.db"
    elif method == "infernal":
        return data_dir / "infernal_dbs" / "descriptions.db"
    else:
        raise ValueError(
            f"Unknown method '{method}' for database '{database_name}'. Cannot determine SQLite path."
        )


def load_descriptions_from_sqlite(
    database_name: str,
    data_dir: Path,
    sseqids: Optional[Set[str]] = None,
    db_config: Optional[Dict] = None,
) -> pd.DataFrame:
    """
    Load feature descriptions from SQLite database.

    Args:
        database_name: Name of the database (fpbase, swissprot, snapgene)
        data_dir: Path to the data directory
        sseqids: Optional set of sequence IDs to filter results
        db_config: Optional database configuration dict

    Returns:
        DataFrame with columns: sseqid, name, type, blurb
    """
    db_path = get_descriptions_db_path(database_name, data_dir, db_config)

    if not db_path.exists():
        print(f"Warning: SQLite database not found: {db_path}")
        # Return empty dataframe with expected columns
        return pd.DataFrame(columns=["sseqid", "name", "type", "blurb"])

    try:
        conn = sqlite3.connect(db_path)

        # Build query based on database type and available columns
        table_name = database_name

        # Check what columns exist in the table
        cursor = conn.execute(f"PRAGMA table_info({table_name})")
        column_info = cursor.fetchall()
        available_columns = [col[1] for col in column_info]  # col[1] is the column name

        # Build query based on available columns, using standardized names
        # All databases should now have type column
        if "type" in available_columns:
            # Standard format with type column
            base_query = f"SELECT sseqid, name, type, blurb FROM {table_name}"
            columns = ["sseqid", "name", "type", "blurb"]
        elif "Type" in available_columns:
            # Legacy format with Type column
            base_query = f"SELECT sseqid, Feature as name, Type as type, Description as blurb FROM {table_name}"
            columns = ["sseqid", "name", "type", "blurb"]
        elif "name" in available_columns and "blurb" in available_columns:
            # Legacy format without type - add default type
            base_query = (
                f"SELECT sseqid, name, 'misc_feature' as type, blurb FROM {table_name}"
            )
            columns = ["sseqid", "name", "type", "blurb"]
        else:
            # Legacy format - map old column names to new ones and add default type
            base_query = f"SELECT sseqid, Feature as name, 'misc_feature' as type, Description as blurb FROM {table_name}"
            columns = ["sseqid", "name", "type", "blurb"]

        # Add WHERE clause if filtering by sseqids
        if sseqids:
            # Convert set to SQL IN clause
            sseqid_list = "','".join(sseqids)
            query = f"{base_query} WHERE sseqid IN ('{sseqid_list}')"
        else:
            query = base_query

        df = pd.read_sql_query(query, conn)
        conn.close()

        print(f"Loaded {len(df)} descriptions from {database_name} SQLite database")
        return df

    except Exception as e:
        print(f"Error loading descriptions from SQLite: {e}")
        return pd.DataFrame(columns=columns)


def query_description_by_sseqid(
    database_name: str,
    data_dir: Path,
    sseqid: str,
    db_config: Optional[Dict] = None,
) -> Optional[Dict[str, str]]:
    """
    Query a single description by sseqid.

    Args:
        database_name: Name of the database
        data_dir: Path to the data directory
        sseqid: Sequence ID to query
        db_config: Optional database configuration dict

    Returns:
        Dictionary with name, type, and blurb, or None if not found
    """
    db_path = get_descriptions_db_path(database_name, data_dir, db_config)

    if not db_path.exists():
        return None

    try:
        conn = sqlite3.connect(db_path)

        table_name = database_name

        # Check what columns exist in the table
        cursor = conn.execute(f"PRAGMA table_info({table_name})")
        column_info = cursor.fetchall()
        available_columns = [col[1] for col in column_info]

        if "type" in available_columns:
            query = f"SELECT name, type, blurb FROM {table_name} WHERE sseqid = ?"
            cursor = conn.execute(query, (sseqid,))
            result = cursor.fetchone()
            conn.close()

            if result:
                return {
                    "name": result[0],
                    "type": result[1],
                    "blurb": result[2],
                }
        elif "Type" in available_columns:
            # Legacy format
            query = (
                f"SELECT Feature, Type, Description FROM {table_name} WHERE sseqid = ?"
            )
            cursor = conn.execute(query, (sseqid,))
            result = cursor.fetchone()
            conn.close()

            if result:
                return {
                    "name": result[0],
                    "type": result[1],
                    "blurb": result[2],
                }
        elif "name" in available_columns and "blurb" in available_columns:
            query = f"SELECT name, 'misc_feature' as type, blurb FROM {table_name} WHERE sseqid = ?"
            cursor = conn.execute(query, (sseqid,))
            result = cursor.fetchone()
            conn.close()

            if result:
                return {"name": result[0], "type": result[1], "blurb": result[2]}
        else:
            # Legacy format
            query = f"SELECT Feature, 'misc_feature' as type, Description FROM {table_name} WHERE sseqid = ?"
            cursor = conn.execute(query, (sseqid,))
            result = cursor.fetchone()
            conn.close()

            if result:
                return {"name": result[0], "type": result[1], "blurb": result[2]}

        return None

    except Exception as e:
        print(f"Error querying description for {sseqid}: {e}")
        return None


def check_sqlite_database(
    database_name: str,
    data_dir: Path,
    db_config: Optional[Dict] = None,
) -> bool:
    """Check if SQLite database exists and is accessible."""
    try:
        db_path = get_descriptions_db_path(database_name, data_dir, db_config)
        if not db_path.exists():
            return False

        conn = sqlite3.connect(db_path)

        # Check if the table exists
        cursor = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (database_name,),
        )
        result = cursor.fetchone()
        conn.close()

        return result is not None

    except Exception:
        return False
