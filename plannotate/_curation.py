"""Curated knowledge that enriches annotations beyond what the search databases hold.

The search databases say *what* a feature is; they do not say what a cloner needs to
know about it. Five lookups fill that gap:

* :func:`selection_marker` -- how a marker gene is selected for (agent, host range).
* :func:`origin_copy_number` -- the plasmid copy number an origin of replication sets.
* :func:`suppressed_feature_accessions` -- source-pinned records excluded from search
  results because they are known global false positives.
* :func:`composite_reference_regions` -- intervals of a composite source record that
  are fully explained by a smaller embedded component.
* :func:`fragment_suppression_regions` -- source intervals known to produce recurring
  low-specificity fragment labels.

A feature name is NOT a safe key: it is a display label, not an identifier, and the
same string can mean different things even within one source. ``cat`` is used by
Swiss-Prot for both chloramphenicol acetyltransferases and catalases, for example.
Every curated row therefore declares the explicit source accession(s) it covers, and
lookups are keyed only on ``(database, sseqid)``. The human-readable ``name`` column is
documentation, not a fallback join key.

NOTE: matching is exact and case-sensitive on purpose because source accessions are
opaque identifiers and must not be normalized beyond surrounding whitespace.
"""

import logging
from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from . import _package_data

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SelectionMarker:
    """How a selection or resistance marker is selected for.

    ``domain`` is the controlled coarse split (``bacterial``/``eukaryotic``/``both``)
    and ``host_range`` the finer breadth, written as ``broad``/``narrow`` plus the
    taxa. ``reference`` carries the PMID a marker's assignment came from, and is
    empty where no primary source was confirmed.
    """

    marker_class: str
    selection_agent: str
    domain: str
    host_range: str
    reference: str


@dataclass(frozen=True)
class OriginCopyNumber:
    """The copy number an origin of replication confers on its plasmid.

    Copy number is strain-, medium-, and growth-rate-dependent, so ``copy_number``
    is a published figure rather than a guarantee and is left empty for origins with
    no measurement behind them -- ``copy_class`` is then the only claim made, and is
    itself ``unreported`` where the literature does not support even that.
    ``reference`` carries the PMIDs the figure came from.
    """

    copy_number: str
    copy_class: str
    domain: str
    host_range: str
    note: str
    reference: str


@dataclass(frozen=True)
class CompositeReferenceRegion:
    """A component interval embedded in a larger source record.

    Coordinates are one-based and inclusive in the nucleotide-equivalent subject
    coordinate system emitted by the annotation adapters. A fragment confined to
    this interval supports the component, not the larger record's label.
    """

    start: int
    end: int
    component_db: str
    component_sseqid: str
    rationale: str
    source: str


@dataclass(frozen=True)
class FragmentSuppressionRegion:
    """A source interval that produces a curated low-specificity fragment artifact."""

    subject_length: int
    start: int
    end: int
    max_identity: float
    name: str
    rationale: str
    source: str


_Key = tuple[str, str]
_Table = dict[_Key, tuple[str, ...]]


def _load_table(filename: str, columns: tuple[str, ...]) -> _Table:
    """Load a curated table into an accession-keyed lookup.

    ``sseqid`` is a ``;``-separated set, so one row covers every accession that carries
    the same claim without duplicating its metadata. ``db`` names exactly one source:
    an accession belongs to a single source vocabulary, so a row spanning two would
    have to pair each accession with a database it does not live in. A claim that holds
    in two sources is therefore two rows. Blank metadata cells read back as NaN, so
    every value is normalized to a string here; callers treat an empty string as
    "not curated" and omit it from output.
    """
    frame = pd.read_csv(_package_data.get_resource("data", filename), dtype=str)
    frame = frame.fillna("")
    table: _Table = {}
    for _, row in frame.iterrows():
        values = tuple(str(row[column]).strip() for column in columns)
        database = str(row["db"]).strip()
        for accession in str(row["sseqid"]).split(";"):
            accession = accession.strip()
            if database and accession:
                table[(database, accession)] = values
    return table


@lru_cache(maxsize=1)
def _selection_markers() -> _Table:
    return _load_table(
        "selection_markers.csv",
        ("marker_class", "selection_agent", "domain", "host_range", "reference"),
    )


@lru_cache(maxsize=1)
def _origin_copy_numbers() -> _Table:
    return _load_table(
        "ori_copy_number.csv",
        ("copy_number", "copy_class", "domain", "host_range", "note", "reference"),
    )


@lru_cache(maxsize=1)
def suppressed_feature_accessions() -> frozenset[_Key]:
    """Return source-pinned accessions excluded from every annotation result.

    This replaces the historical unscoped list in :mod:`._filter`. Keeping the source
    database in the key prevents an identifier used by a custom or future database
    from inheriting an unrelated suppression.
    """
    table = _load_table(
        "feature_suppressions.csv",
        ("name", "rationale", "reference"),
    )
    return frozenset(table)


@lru_cache(maxsize=1)
def composite_reference_regions() -> dict[_Key, tuple[CompositeReferenceRegion, ...]]:
    """Return curated embedded-component intervals keyed by source record."""
    frame = pd.read_csv(
        _package_data.get_resource("data", "composite_reference_regions.csv"),
        dtype=str,
    ).fillna("")
    regions: dict[_Key, list[CompositeReferenceRegion]] = {}
    seen: set[tuple[str, str, int, int, str, str]] = set()
    for _, row in frame.iterrows():
        database = str(row["db"]).strip()
        accession = str(row["sseqid"]).strip()
        component_database = str(row["component_db"]).strip()
        component_accession = str(row["component_sseqid"]).strip()
        rationale = str(row["rationale"]).strip()
        source = str(row["source"]).strip()
        try:
            start = int(str(row["region_start"]).strip())
            end = int(str(row["region_end"]).strip())
        except ValueError as error:
            raise ValueError(
                "Composite reference coordinates must be integers"
            ) from error
        if min(start, end) < 1 or start > end:
            raise ValueError(f"Invalid composite reference interval: {start}-{end}")
        if not all(
            (
                database,
                accession,
                component_database,
                component_accession,
                rationale,
                source,
            )
        ):
            raise ValueError(
                "Composite reference regions require pinned ids and provenance"
            )
        unique_key = (
            database,
            accession,
            start,
            end,
            component_database,
            component_accession,
        )
        if unique_key in seen:
            raise ValueError(f"Duplicate composite reference region: {unique_key!r}")
        seen.add(unique_key)
        regions.setdefault((database, accession), []).append(
            CompositeReferenceRegion(
                start,
                end,
                component_database,
                component_accession,
                rationale,
                source,
            )
        )
    return {
        key: tuple(sorted(values, key=lambda region: (region.start, region.end)))
        for key, values in regions.items()
    }


@lru_cache(maxsize=1)
def fragment_suppression_regions() -> dict[_Key, tuple[FragmentSuppressionRegion, ...]]:
    """Return manually curated low-specificity fragment intervals by source id."""
    frame = pd.read_csv(
        _package_data.get_resource("data", "fragment_suppression_regions.csv"),
        dtype=str,
    ).fillna("")
    regions: dict[_Key, list[FragmentSuppressionRegion]] = {}
    seen: set[tuple[str, str, int, int, int, float]] = set()
    for _, row in frame.iterrows():
        database = str(row["db"]).strip()
        accession = str(row["sseqid"]).strip()
        name = str(row["name"]).strip()
        rationale = str(row["rationale"]).strip()
        source = str(row["source"]).strip()
        try:
            subject_length = int(str(row["subject_length"]).strip())
            start = int(str(row["region_start"]).strip())
            end = int(str(row["region_end"]).strip())
            max_identity = float(str(row["max_identity"]).strip())
        except ValueError as error:
            raise ValueError(
                "Fragment suppression geometry and identity must be numeric"
            ) from error
        if min(subject_length, start, end) < 1 or start > end or end > subject_length:
            raise ValueError(
                "Invalid fragment suppression geometry: "
                f"length {subject_length}, interval {start}-{end}"
            )
        if not 0 <= max_identity <= 100:
            raise ValueError(f"Invalid fragment suppression identity: {max_identity}")
        if not all((database, accession, name, rationale, source)):
            raise ValueError(
                "Fragment suppression regions require pinned ids and provenance"
            )
        unique_key = (
            database,
            accession,
            subject_length,
            start,
            end,
            max_identity,
        )
        if unique_key in seen:
            raise ValueError(f"Duplicate fragment suppression region: {unique_key!r}")
        seen.add(unique_key)
        regions.setdefault((database, accession), []).append(
            FragmentSuppressionRegion(
                subject_length,
                start,
                end,
                max_identity,
                name,
                rationale,
                source,
            )
        )
    return {
        key: tuple(sorted(values, key=lambda region: (region.start, region.end)))
        for key, values in regions.items()
    }


def _lookup(table: _Table, database: str, sseqid: str) -> tuple[str, ...] | None:
    """Resolve a hit against a curated table by source accession."""
    return table.get((database.strip(), sseqid.strip()))


def selection_marker(database: str, sseqid: str) -> SelectionMarker | None:
    """Return marker details for a source record, or None if it is not a marker.

    ``database`` and ``sseqid`` identify the packaged record. The feature name is
    deliberately not a parameter: accepting one that is never read invites callers to
    believe it selects something.
    """
    values = _lookup(_selection_markers(), database, sseqid)
    return None if values is None else SelectionMarker(*values)


def origin_copy_number(database: str, sseqid: str) -> OriginCopyNumber | None:
    """Return copy-number details for a source record, or None if it is not an origin."""
    values = _lookup(_origin_copy_numbers(), database, sseqid)
    return None if values is None else OriginCopyNumber(*values)
