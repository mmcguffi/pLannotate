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
from Bio import SeqIO

from . import _package_data, _sqlite, annotate
from ._tools import diamond
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


def _normalize_subject_ids(sequence_ids: pd.Series, method: str) -> pd.Series:
    """Key descriptions by the id the search pipeline reports for a hit.

    A synthesized descriptions table must be keyed by the *normalized* subject id or
    the per-hit merge in ``annotate._enrich_hits`` misses. Rather than reimplement
    the normalization (and risk silent drift), reuse the runtime's own steps:
    ``diamond.normalize_subject_ids`` unwraps a pipe-delimited protein id to its
    accession (``sp|P12345|NAME`` -> ``P12345``), and ``annotate.strip_pdb_wrapper``
    strips a PDB wrapper (``pdb|1ABC|`` -> ``1ABC``) for every method.
    """
    normalized = sequence_ids.astype(str)
    if method == "diamond":
        normalized = diamond.normalize_subject_ids(normalized)
    return annotate.strip_pdb_wrapper(normalized)


def descriptions_from_fasta(fasta: Path, method: str) -> pd.DataFrame:
    """Synthesize a descriptions table from FASTA headers when no CSV is supplied.

    The record id (normalized to match how the search tool reports the hit's
    ``sseqid``, see :func:`_normalize_subject_ids`) is the key; the original id
    becomes the feature ``name`` and any text after it on the header line becomes
    the ``blurb``. Protein features default to ``CDS`` (they are almost always
    coding) and nucleotide features to ``misc_feature``. Always producing a
    descriptions table keeps the ``--csv`` argument genuinely optional: the
    annotation path needs the table to exist.
    """
    default_type = "CDS" if method == "diamond" else "misc_feature"
    rows = []
    for record in SeqIO.parse(str(fasta), "fasta"):
        # everything after the first whitespace-delimited token is the blurb; split
        # rather than slice off the id so a header like "> feature" (leading space,
        # which Biopython keeps in the description) does not yield a garbage blurb.
        header_parts = record.description.split(None, 1)
        blurb = header_parts[1] if len(header_parts) > 1 else ""
        rows.append(
            {
                "sseqid": record.id,
                "name": record.id,  # original id kept as the label
                "type": default_type,
                "blurb": blurb,
            }
        )
    if not rows:
        raise ValueError(f"No sequences found in {fasta}")
    frame = pd.DataFrame(rows)
    return frame.assign(sseqid=_normalize_subject_ids(frame["sseqid"], method))


def _descriptions_frame(fasta: Path, csv: Path | None, method: str) -> pd.DataFrame:
    """Build the descriptions table for a source from a CSV or the FASTA headers.

    Either way the returned ``sseqid`` is keyed by the id the search tool actually
    reports for a hit (see :func:`_normalize_subject_ids`): DIAMOND unwraps a
    pipe-delimited id to its accession, so a raw ``sp|P12345|NAME`` id -- whether it
    came from a UniProt-style FASTA header or a user CSV whose ids match it -- must
    be normalized here or the per-hit description merge silently misses.
    """
    if csv is None:
        # descriptions_from_fasta already keys by the normalized id
        return descriptions_from_fasta(fasta, method)
    frame = _sqlite._normalize_description_frame(
        pd.read_csv(csv, sep=None, engine="python")
    )
    return frame.assign(sseqid=_normalize_subject_ids(frame["sseqid"], method))


def build_descriptions(name: str, frame: pd.DataFrame, out_dir: Path) -> Path:
    """Write the ``<name>.db`` descriptions database beside the index."""
    db_path = out_dir / f"{name}.db"
    count = _sqlite.write_descriptions_to_sqlite(name, frame, db_path)
    logger.info("Wrote %d feature descriptions to %s", count, db_path)
    return db_path


def build_source_config(
    name: str, method: str, out_dir: Path, priority: int
) -> dict[str, Any]:
    """Build the YAML source entry for a freshly created custom database.

    A descriptions database is always written beside the index (from the CSV or
    synthesized from the FASTA), so ``details.location`` is ``Default`` -- it
    resolves to that file -- and ``default_type`` stays ``None`` so the table's own
    types are respected.
    """
    return {
        "method": method,
        "location": str(out_dir.resolve()),
        "priority": priority,
        "cost": 1.0,
        "parameters": list(_DEFAULT_PARAMETERS[method]),
        "details": {"default_type": None, "location": "Default"},
    }


def build_full_config(
    name: str,
    method: str,
    out_dir: Path,
    priority: int,
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
    config[name] = build_source_config(name, method, out_dir, priority)
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
    """Build the index, descriptions, and YAML config for a custom database.

    ``method`` is normalized here, so callers may pass any accepted alias. When
    ``csv`` is omitted the descriptions are synthesized from the FASTA headers (see
    :func:`descriptions_from_fasta`). Returns the full source configuration (see
    :func:`build_full_config`); writing it to disk is left to the caller.
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
    build_descriptions(name, _descriptions_frame(fasta, csv, canonical), out_dir)
    return build_full_config(name, canonical, out_dir, priority, include_builtins)
