"""Curated knowledge that enriches annotations beyond what the search databases hold.

The search databases say *what* a feature is; they do not say what a cloner needs to
know about it. Two lookups fill that gap:

* :func:`selection_marker` -- how a marker gene is selected for (agent, host range).
* :func:`origin_copy_number` -- the plasmid copy number an origin of replication sets.

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
