"""Build ad-hoc annotation databases from a FASTA and optional descriptions CSV.

This is the machinery behind ``plannotate makedb``. It turns a user-supplied FASTA
into a BLAST or DIAMOND search index, an optional SQLite descriptions database, and
a ready-to-run YAML configuration that layers the new source on top of the builtin
databases -- so an end user can add a custom database without touching the internal
build (Snakemake) workflow.

The generated layout mirrors what a configured source expects (see
``_package_data.get_yaml``): the index and its ``<name>.db`` descriptions file share
one output directory, the index basename equals the source name, and the emitted
YAML points ``location`` at that directory with ``details.location: Default`` so the
descriptions file is found beside the index.
"""

import logging
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from . import _package_data, _sqlite
from ._tools.common import run_command

logger = logging.getLogger(__name__)

# Map the user-facing method spellings onto the canonical YAML ``method`` values.
# ``blastn`` searches a nucleotide index; ``diamond`` searches a translated protein
# index. Aliases keep the CLI forgiving without widening the stored vocabulary.
METHOD_ALIASES: dict[str, str] = {
    "blast": "blastn",
    "blastn": "blastn",
    "nucl": "blastn",
    "nucleotide": "blastn",
    "diamond": "diamond",
    "prot": "diamond",
    "protein": "diamond",
}

# Tool binaries each method shells out to; checked up front for a clearer error than
# a mid-build failure.
_METHOD_TOOL: dict[str, str] = {"blastn": "makeblastdb", "diamond": "diamond"}

# Sensible starting parameters per method, mirroring the builtin sources so a fresh
# custom database behaves reasonably out of the box. Users can tune them by editing
# the generated YAML.
_DEFAULT_PARAMETERS: dict[str, list[str]] = {
    "blastn": [
        "-perc_identity 95",
        "-max_target_seqs 20000",
        "-culling_limit 25",
        "-word_size 12",
    ],
    "diamond": [
        "-k 0",
        "--min-orf 1",
        "--matrix BLOSUM90",
        "--gapopen 10",
        "--gapextend 1",
        "--algo ctg",
        "--id 50",
        "--max-hsps 10",
        "--culling-overlap 200",
        "--seed-cut .001",
        "--comp-based-stats 0",
    ],
}


def normalize_method(method: str) -> str:
    """Resolve a user-supplied method spelling to a canonical YAML method."""
    canonical = METHOD_ALIASES.get(method.strip().lower())
    if canonical is None:
        supported = ", ".join(sorted(set(METHOD_ALIASES)))
        raise ValueError(f"Unsupported method {method!r}; choose one of: {supported}")
    return canonical


def _index_command(method: str, fasta: Path, out_base: Path) -> list[str]:
    """Build the makeblastdb/diamond command that writes an index named ``out_base``."""
    if method == "blastn":
        return [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(out_base),
        ]
    # diamond appends the .dmnd suffix itself, so pass the bare base path
    return ["diamond", "makedb", "--in", str(fasta), "--db", str(out_base)]


def build_search_index(fasta: Path, name: str, method: str, out_dir: Path) -> Path:
    """Create a BLAST or DIAMOND index for ``fasta`` named ``name`` under ``out_dir``.

    Returns the index base path (no tool-specific suffix), which is what the runtime
    stores as ``db_loc``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_base = out_dir / name
    logger.info("Building %s index %r from %s", method, name, fasta)
    run_command(_index_command(method, fasta, out_base), _METHOD_TOOL[method])
    return out_base


def build_descriptions(name: str, csv: Path, out_dir: Path) -> Path:
    """Build the ``<name>.db`` descriptions database beside the index from a CSV."""
    frame = pd.read_csv(csv, sep=None, engine="python")
    db_path = out_dir / f"{name}.db"
    count = _sqlite.write_descriptions_to_sqlite(name, frame, db_path)
    logger.info("Wrote %d feature descriptions to %s", count, db_path)
    return db_path


def build_source_config(
    name: str, method: str, out_dir: Path, priority: int, has_descriptions: bool
) -> dict[str, Any]:
    """Build the YAML source entry for a freshly created custom database."""
    if has_descriptions:
        # descriptions.db sits beside the index, so Default resolves to it; keeping
        # default_type None lets the CSV's own types through unchanged.
        details = {"default_type": None, "location": "Default"}
    else:
        # no descriptions file: synthesize feature names from hit ids, and give
        # protein hits a CDS default type since they are almost always coding.
        details = {
            "default_type": "CDS" if method == "diamond" else None,
            "location": None,
        }
    return {
        "method": method,
        "location": str(out_dir.resolve()),
        "priority": priority,
        "cost": 1.0,
        "parameters": list(_DEFAULT_PARAMETERS[method]),
        "details": details,
    }


def build_full_config(
    name: str,
    method: str,
    out_dir: Path,
    priority: int,
    has_descriptions: bool,
    include_builtins: bool,
) -> dict[str, Any]:
    """Assemble the complete source configuration for the new database.

    When ``include_builtins`` is set, the packaged sources are carried through
    verbatim (still pointing at their Default locations) and the custom source is
    appended, so a single ``--yaml-file`` annotates with builtins plus the new
    database. Otherwise the returned config contains only the custom source.
    """
    config: dict[str, Any] = {}
    if include_builtins:
        with _package_data.get_yaml_path().open() as handle:
            builtins = yaml.safe_load(handle)
        if not isinstance(builtins, dict):
            raise ValueError("Builtin database configuration is malformed")
        if name in builtins:
            raise ValueError(
                f"Database name {name!r} collides with a builtin source; "
                "choose a different --name"
            )
        config.update(builtins)
    config[name] = build_source_config(
        name, method, out_dir, priority, has_descriptions
    )
    return config


def make_database(
    fasta: Path,
    name: str,
    method: str,
    out_dir: Path,
    csv: Path | None = None,
    priority: int = 1,
    include_builtins: bool = True,
) -> dict[str, Any]:
    """Build index, optional descriptions, and the YAML config for a custom database.

    ``method`` is normalized here, so callers may pass any accepted alias. Returns
    the full source configuration (see :func:`build_full_config`); writing it to disk
    is left to the caller.
    """
    canonical = normalize_method(method)
    name = _sqlite._validated_table_name(name)
    tool = _METHOD_TOOL[canonical]
    if shutil.which(tool) is None:
        raise RuntimeError(
            f"{tool!r} was not found on PATH; install it to build a {canonical} "
            "database"
        )

    build_search_index(fasta, name, canonical, out_dir)
    has_descriptions = csv is not None
    if csv is not None:
        build_descriptions(name, csv, out_dir)
    return build_full_config(
        name, canonical, out_dir, priority, has_descriptions, include_builtins
    )
