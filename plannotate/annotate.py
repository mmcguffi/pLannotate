"""Collect, rank, and finalize annotations for a DNA sequence."""

import logging
from functools import lru_cache, partial
from pathlib import Path
from typing import Any, cast

import pandas as pd
from Bio.Seq import Seq

from . import _concurrency, _package_data, _sqlite
from ._filter import filter_and_clean_hits
from ._schema import ADAPTER_COLUMNS, ANNOTATION_COLUMNS
from ._tools.methods import run as run_tool

logger = logging.getLogger(__name__)

# Sources retained in fast mode: the nucleotide snapgene search and the tiny
# fpbase diamond db are the two cheapest sources (cost 0.45 and 0.03), yet cover
# the bulk of common engineered features. Fast mode skips the costly swissprot
# and Rfam scans, which together account for nearly all annotation runtime.
FAST_SOURCES = frozenset({"snapgene", "fpbase"})


def _circular_search_query(query: str) -> str:
    """Double a circular query so features spanning the linear seam are recovered.

    The whole sequence is searched twice end-to-end, which is lossless regardless
    of where the arbitrary circular->linear seam falls. Trimming the second copy
    to a window was tried (both per-source and at a rotated-ori seam) but perturbs
    the genome-wide overlap/culling resolution non-locally -- it can drop or swap
    low-quality fragment calls far from the seam -- so the full second copy is kept.
    """
    return query + query


def _collect_source_hits(
    query: str,
    source_name: str,
    source_config: dict[str, Any],
    is_linear: bool,
    threads: int = 1,
    fast: bool = False,
) -> pd.DataFrame:
    """Collect and enrich hits from one configured annotation source.

    A circular search normally doubles the query so a feature crossing the linear
    seam is recovered as one alignment. In ``fast`` mode the query is searched
    once and seam-spanning features are reconstructed afterwards by
    :func:`_stitch_seam_hits`, halving the search at the cost of a little seam
    fidelity (see :data:`FAST_SOURCES`).
    """
    if is_linear or fast:
        search_query = query
    else:
        search_query = _circular_search_query(query)
    hits = run_tool(search_query, source_config, threads=threads)
    if hits.empty:
        return pd.DataFrame()

    missing_columns = set(ADAPTER_COLUMNS) - set(hits.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Source {source_name!r} did not return columns: {missing}")
    # Record the true sequence length rather than the doubled search length so the
    # circular wrap in _filter and the construct geometry use the real size.
    hits["qlen"] = len(query)
    if fast and not is_linear:
        hits = _stitch_seam_hits(hits)
    return _enrich_hits(hits, source_name, source_config)


# Tolerance, in subject units, for treating two terminal fragments as one
# seam-spanning feature. Small so distinct features are never fused, but >0 so a
# ragged blast end or a diamond codon split exactly at the seam does not block a
# real merge.
_SEAM_SUBJECT_TOLERANCE = 3


def _stitch_seam_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """Fuse terminal fragment pairs into origin-spanning hits (fast mode).

    Without the doubled query, a feature crossing the circular seam is reported as
    two partial alignments -- one ending at the last base, one starting at the
    first. When such a pair shares a subject, a strand, and contiguous subject
    coordinates, it is one feature split by the seam, so it is merged into a single
    hit whose query end is pushed past ``qlen``. The existing circular wrap in
    :mod:`._filter` then folds that back into an origin-spanning interval, exactly
    as the doubled-query path does.
    """
    if hits.empty:
        return hits

    qlen = int(hits["qlen"].iloc[0])
    q_lo = hits[["qstart", "qend"]].min(axis=1)
    q_hi = hits[["qstart", "qend"]].max(axis=1)
    s_lo = hits[["sstart", "send"]].min(axis=1)
    s_hi = hits[["sstart", "send"]].max(axis=1)
    right_fragments = hits.index[q_hi == qlen].tolist()  # touch the last base
    left_fragments = hits.index[q_lo == 1].tolist()  # touch the first base
    if not right_fragments or not left_fragments:
        return hits

    merged_rows: list[dict[str, Any]] = []
    consumed: set[int] = set()
    for right in right_fragments:
        if right in consumed:
            continue
        best_left: int | None = None
        best_gap = _SEAM_SUBJECT_TOLERANCE + 1
        for left in left_fragments:
            if left == right or left in consumed:
                continue
            if hits.at[right, "sseqid"] != hits.at[left, "sseqid"]:
                continue
            if hits.at[right, "sframe"] != hits.at[left, "sframe"]:
                continue
            # the subject must continue where the partner stops; check both
            # orderings so forward and reverse strands are handled alike
            gap = min(
                abs(int(s_lo[left]) - (int(s_hi[right]) + 1)),
                abs(int(s_lo[right]) - (int(s_hi[left]) + 1)),
            )
            if gap < best_gap:
                best_gap = gap
                best_left = left
        if best_left is None:
            continue
        consumed.add(right)
        consumed.add(best_left)
        merged_rows.append(_merge_seam_pair(hits.loc[right], hits.loc[best_left], qlen))

    if not merged_rows:
        return hits
    kept = hits.drop(index=list(consumed))
    return pd.concat([kept, pd.DataFrame(merged_rows)], ignore_index=True)


def _merge_seam_pair(right: pd.Series, left: pd.Series, qlen: int) -> dict[str, Any]:
    """Combine a query-3' fragment and a query-5' fragment into one wrapped hit."""
    right_start = min(int(right["qstart"]), int(right["qend"]))
    left_end = max(int(left["qstart"]), int(left["qend"]))
    total_length = int(right["length"]) + int(left["length"])
    # length-weighted identity so a short, weaker fragment cannot dominate
    pident = (
        right["pident"] * int(right["length"]) + left["pident"] * int(left["length"])
    ) / total_length

    merged = right.copy()
    # span runs from the 3'-end fragment through the seam (qlen) into the 5'-start
    # fragment; the >qlen end is wrapped later by _adjust_circular_coordinates
    merged["qstart"] = right_start
    merged["qend"] = qlen + left_end
    merged["length"] = total_length
    merged["pident"] = pident
    merged["evalue"] = min(float(right["evalue"]), float(left["evalue"]))
    # subject coords are informational after filtering; keep a sane combined view
    merged["sstart"] = min(int(right["sstart"]), int(left["sstart"]))
    merged["send"] = max(int(right["send"]), int(left["send"]))
    # plus-strand query order across the seam is the 3'-end fragment then the 5'
    merged["qseq"] = str(right["qseq"]) + str(left["qseq"])
    return cast(dict[str, Any], merged.to_dict())


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
        # details are synthesized from the hits themselves, so collapse repeats to one
        # row per sseqid; otherwise the merge below fans out when an id occurs twice.
        details = (
            hits[["sseqid", "name", "type", "blurb"]]
            .drop_duplicates(subset="sseqid")
            .copy()
        )
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
    fast: bool,
) -> pd.DataFrame:
    """Cache an immutable snapshot of candidate hits for identical inputs."""
    # yaml_modified_ns is part of the cache key so edits to the YAML invalidate it.
    del yaml_modified_ns
    sources = _package_data.get_yaml(Path(yaml_file_str))
    if fast:
        sources = {
            name: config for name, config in sources.items() if name in FAST_SOURCES
        }
    if not sources:
        return pd.DataFrame()

    results = _concurrency.run_sources(
        partial(_collect_source_hits, fast=fast),
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
    fast: bool = False,
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
        fast,
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
    fast: bool = False,
) -> pd.DataFrame:
    """Annotate a DNA sequence and return results as a DataFrame.

    Circular sequences are fully doubled so origin-spanning features are never
    missed. ``fast`` restricts the search to the cheapest sources (see
    :data:`FAST_SOURCES`) for a quicker, lower-coverage annotation.
    """
    yaml_file = (
        Path(yaml_file) if yaml_file is not None else _package_data.get_yaml_path()
    )
    sequence = Seq(seq)

    logger.info("Collecting candidate annotations")
    hits = _collect_hits(str(sequence), linear, yaml_file, cores, fast)
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
