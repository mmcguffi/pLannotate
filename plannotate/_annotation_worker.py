"""Execute one database-search job for the annotation Snakemake workflow."""

import argparse
from pathlib import Path

from . import resources as rsc
from .search import _search_single_database


def run_database_search(
    query_file: Path,
    yaml_file: Path,
    database_index: int,
    output_file: Path,
    linear: bool = False,
    threads: int = 1,
) -> None:
    """Search one configured database and serialize its processed hits."""
    databases = rsc.get_yaml(yaml_file)
    try:
        database_name, database_config = list(databases.items())[database_index]
    except IndexError as exc:
        raise ValueError(f"Unknown database index: {database_index}") from exc

    hits = _search_single_database(
        query_file.read_text().strip(),
        database_name,
        database_config,
        yaml_file,
        linear,
        threads,
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    hits.to_pickle(output_file)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--yaml", type=Path, required=True)
    parser.add_argument("--database-index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--linear", action="store_true")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    run_database_search(
        query_file=args.query,
        yaml_file=args.yaml,
        database_index=args.database_index,
        output_file=args.output,
        linear=args.linear,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
