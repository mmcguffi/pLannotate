"""Collect, rank, and finalize annotations for a DNA sequence."""

import logging
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Seq import Seq

from . import _concurrency, _package_data, _sqlite
from ._filter import filter_and_clean_hits
from ._schema import ADAPTER_COLUMNS, ANNOTATION_COLUMNS
from ._tools.methods import run as run_tool

logger = logging.getLogger(__name__)


def _collect_source_hits(
    query: str,
    source_name: str,
    source_config: dict[str, Any],
    is_linear: bool,
    threads: int = 1,
) -> pd.DataFrame:
    """Collect and enrich hits from one configured annotation source."""
    search_query = query if is_linear else query + query
    hits = run_tool(search_query, source_config, threads=threads)
    if hits.empty:
        return pd.DataFrame()

    missing_columns = set(ADAPTER_COLUMNS) - set(hits.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Source {source_name!r} did not return columns: {missing}")
    return _enrich_hits(hits, source_name, source_config)


def _enrich_hits(
    hits: pd.DataFrame,
    source_name: str,
    source_config: dict[str, Any],
) -> pd.DataFrame:
    """Attach descriptions, feature types, and priority to raw hits."""
    enriched = hits.copy()
    enriched["db"] = source_name
    enriched["sseqid"] = (
        enriched["sseqid"].astype(str).str.replace(r"pdb\|(.*)\|", r"\1", regex=True)
    )
    details = _load_feature_details(enriched, source_name, source_config)
    enriched = enriched.merge(
        details,
        on="sseqid",
        how="left",
        suffixes=("_original", ""),
    )
    duplicate_columns = [
        column for column in enriched if str(column).endswith("_original")
    ]
    enriched = enriched.drop(columns=duplicate_columns)
    if "type" in enriched.columns:
        enriched = enriched.loc[enriched["type"] != "primer_bind"]

    enriched["priority"] = source_config["priority"]
    if "priority_mod" in enriched.columns:
        # a left-merge miss (hit absent from the descriptions DB) leaves priority_mod
        # NaN; treat it as no penalty so it cannot poison the score downstream.
        priority_mod = enriched.pop("priority_mod").fillna(0)
        enriched["priority"] = (enriched["priority"] + priority_mod).astype(int)
    return enriched


def _load_feature_details(
    hits: pd.DataFrame,
    source_name: str,
    source_config: dict[str, Any],
) -> pd.DataFrame:
    detail_config = source_config["details"]
    detail_location = detail_config["location"]
    if detail_location is None or detail_location == "None":
        details = hits[["sseqid", "name", "type", "blurb"]].copy()
    else:
        sequence_ids = {
            identifier for identifier in hits["sseqid"].tolist() if identifier
        }
        details = _sqlite.load_descriptions_from_sqlite(
            source_name,
            sequence_ids,
            source_config,
        )

    if detail_config.get("priority_from_protein_existence", False):
        details["priority_mod"] = details["blurb"].map(_existence_level_priority)
    default_type = detail_config.get("default_type")
    if default_type and default_type != "None":
        details["type"] = default_type
    return details


def _existence_level_priority(description: object) -> int:
    """Translate a Swiss-Prot existence level into a priority penalty."""
    if not isinstance(description, str):
        return 0
    marker = "existence level"
    marker_position = description.find(marker)
    if marker_position == -1:
        return 0
    try:
        return int(description[marker_position + len(marker) + 1]) - 1
    except (ValueError, IndexError):
        return 0


# Cores only affects search parallelism, never the resulting hits, so it is kept out
# of the cache key below. The value in effect for a cache miss is read from here.
_active_cores = 1


@lru_cache(maxsize=5)
def _collect_hits_cached(
    query_sequence: str,
    is_linear: bool,
    yaml_file_str: str,
    yaml_modified_ns: int,
) -> pd.DataFrame:
    """Cache an immutable snapshot of candidate hits for identical inputs."""
    # yaml_modified_ns is part of the cache key so edits to the YAML invalidate it.
    del yaml_modified_ns
    sources = _package_data.get_yaml(Path(yaml_file_str))
    if not sources:
        return pd.DataFrame()

    results = _concurrency.run_sources(
        _collect_source_hits,
        query_sequence,
        is_linear,
        sources,
        _active_cores,
    )
    results = [result for result in results if not result.empty]
    if not results:
        return pd.DataFrame()
    return pd.concat(results, ignore_index=True)


def _collect_hits(
    query_sequence: str,
    is_linear: bool,
    yaml_file: Path,
    cores: int,
) -> pd.DataFrame:
    """Collect an independent copy of all candidate hits."""
    global _active_cores
    yaml_file = yaml_file.resolve()
    _active_cores = cores
    cached = _collect_hits_cached(
        query_sequence,
        is_linear,
        str(yaml_file),
        yaml_file.stat().st_mtime_ns,
    )
    return cached.copy(deep=True)


def _is_fragment(feature: pd.Series) -> bool:
    """Determine if a feature is a fragment based on type and match quality."""
    if "type" not in feature.index:
        return False
    if feature["type"] != "CDS":
        return bool(feature["percmatch"] < 95)
    is_complete_cds = feature["pi_permatch"] == 100 or (
        feature["length"] % 3 == 0 and feature["percmatch"] > 95
    )
    return not is_complete_cds


def _empty_annotations() -> pd.DataFrame:
    return pd.DataFrame(columns=ANNOTATION_COLUMNS)


def _orient_query_sequence(feature: pd.Series) -> str:
    query_sequence = str(feature["qseq"])
    if feature["sframe"] == -1:
        return str(Seq(query_sequence).reverse_complement())
    return query_sequence


def annotate(
    seq: str | Seq,
    yaml_file: Path | None = None,
    linear: bool = False,
    is_detailed: bool = False,
    cores: int = 1,
) -> pd.DataFrame:
    """Annotate a DNA sequence and return results as a DataFrame."""
    yaml_file = (
        Path(yaml_file) if yaml_file is not None else _package_data.get_yaml_path()
    )
    sequence = Seq(seq)

    logger.info("Collecting candidate annotations")
    hits = _collect_hits(str(sequence), linear, yaml_file, cores)
    if hits.empty:
        return _empty_annotations()

    hits["kind"] = hits["type"] if is_detailed else 1
    hits = filter_and_clean_hits(hits, linear)
    if hits.empty:
        return _empty_annotations()

    hits["fragment"] = hits.apply(_is_fragment, axis=1)
    hits["qend"] += 1
    hits["qseq"] = hits.apply(_orient_query_sequence, axis=1)
    hits["name"] = hits["name"].fillna(hits["sseqid"])
    hits["blurb"] = hits["blurb"].fillna("")
    if "type" in hits.columns:
        hits["type"] = hits["type"].fillna("misc_feature")
    else:
        hits["type"] = "misc_feature"

    logger.info("Annotation complete: %d features identified", len(hits))
    return hits
