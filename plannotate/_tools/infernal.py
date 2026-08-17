"""Infernal covariance-model-search integration."""

import logging
import re
import shlex
from collections.abc import Iterator, Mapping
from functools import lru_cache
from itertools import accumulate
from pathlib import Path
from typing import Any

import pandas as pd

from .._concurrency import parameters_with_threads
from .common import normalize_queries, run_command, temporary_files

REQUIRED_COLUMNS = [
    "#idx",
    "target name",
    "accession",
    "query name",
    "clan name",
    "seq from",
    "seq to",
    "mdl from",
    "mdl to",
    "strand",
    "score",
    "E-value",
    "description of target",
]
logger = logging.getLogger(__name__)
# the per-hit row inside an alignment block, e.g. "  (1) !   5.5e-13   73.7 ..."
_HIT_HEADER = re.compile(r"\s+\(\d+\) [!?]\s")


def search(
    sequence: str | Mapping[str, str],
    config: dict[str, Any],
    threads: int = 1,
) -> pd.DataFrame:
    """Search a covariance-model database with Infernal cmscan.

    ``sequence`` may be one sequence or a ``{query_id: sequence}`` mapping; the
    returned frame carries a ``qseqid`` column identifying each hit's query.
    """
    logger.info("Starting Infernal search")
    logger.debug(
        "Infernal database=%s clan=%s threads=%d",
        config.get("db_loc"),
        config.get("clanin_loc"),
        threads,
    )
    parameters = parameters_with_threads(
        str(config["parameters"]), ("--cpu",), "--cpu", threads
    )
    cm_path, clan_path = _configured_database_paths(config)
    with temporary_files(sequence) as (query_path, output_path):
        # the alignment report goes beside the table in the same temporary directory;
        # cmscan computes the alignment either way, so writing it out rather than
        # discarding it with --noali costs nothing but buys the consensus structure
        alignment_path = str(Path(output_path).with_name("alignments.txt"))
        command = [
            "cmscan",
            "--cut_ga",
            "--rfam",
            "--fmt",
            "2",
            *shlex.split(parameters),
            "-o",
            alignment_path,
            "--tblout",
            output_path,
            "--clanin",
            clan_path,
            cm_path,
            query_path,
        ]
        run_command(command, "cmscan")
        dataframe = parse_output(output_path, alignment_path, cm_path)

    # qlen/qseq are derived from the matching query, looked up by qseqid, so a
    # batched multi-FASTA run resolves each hit against its own sequence. The qlen
    # set here is the searched (possibly doubled) length; the collector overrides it
    # with the true plasmid length, matching the blast/diamond path.
    queries = normalize_queries(sequence)
    dataframe["qlen"] = [len(queries[query_id]) for query_id in dataframe["qseqid"]]
    if not dataframe.empty:
        dataframe["qseq"] = dataframe.apply(
            lambda row: queries[row["qseqid"]][row["qstart"] - 1 : row["qend"]].upper(),
            axis=1,
        )
    logger.info("Infernal found %d candidate hits", len(dataframe))
    return dataframe


def database_paths(database_name: str, directory: Path) -> dict[str, str]:
    """Build the CM and clan paths for an Infernal database."""
    return {
        "db_loc": str(directory / f"{database_name}.cm"),
        "clanin_loc": str(directory / f"{database_name}.clanin"),
    }


def _configured_database_paths(config: dict[str, Any]) -> tuple[str, str]:
    """Read the configured covariance-model and clan paths."""
    if "db_loc" not in config or "clanin_loc" not in config:
        raise ValueError("Infernal configuration requires CM and clan database paths")
    return str(config["db_loc"]), str(config["clanin_loc"])


def _coordinate_text(column: pd.Series) -> list[str]:
    """Render a coordinate column the way cmscan prints it in its alignment report."""
    return ["" if pd.isna(value) else str(int(float(value))) for value in column]


def _alignment_records(text: str) -> Iterator[tuple[str, str, list[str]]]:
    """Yield each ``>>`` block in cmscan's report as (query name, model name, lines).

    The query is tracked from the enclosing ``Query:`` header rather than the block
    itself, because a block names only its model -- and one report covers every
    sequence in a batched search.
    """
    query = ""
    block: list[str] = []
    model = ""
    for line in text.splitlines():
        if line.startswith("Query:"):
            if block:
                yield query, model, block
                block = []
            fields = line.split()
            query = fields[1] if len(fields) > 1 else ""
            continue
        if line.startswith(">> "):
            if block:
                yield query, model, block
            model = line[3:].split()[0]
            block = []
            continue
        if model:
            block.append(line)
    if block:
        yield query, model, block


def parse_alignments(
    path: str | Path,
) -> dict[tuple[str, str, str, str], tuple[str, float]]:
    """Map each hit to its consensus structure and posterior-probability accuracy.

    cmscan's alignment report is a human-readable format that Infernal makes no
    stability promise about, unlike ``--tblout``. Every failure therefore degrades to
    "no structure for this hit" rather than failing the search: the structure is an
    enrichment, and losing it must never cost an annotation. Each block is parsed
    independently so one unrecognised hit cannot discard the hits around it.

    Hits are keyed by ``(query, model, seq from, seq to)``: a batched search puts every
    sequence through one cmscan, and two plasmids sharing a backbone can hit the same
    model at the same offset, so the query is part of a hit's identity. The coordinates
    are the raw ones cmscan prints -- a minus-strand hit reports them descending,
    before the table parser sorts them.
    """
    try:
        text = Path(path).read_text()
    except OSError:
        logger.debug("No cmscan alignment report at %s", path)
        return {}

    alignments: dict[tuple[str, str, str, str], tuple[str, float]] = {}
    for query, model, block in _alignment_records(text):
        header = next((line for line in block if _HIT_HEADER.match(line)), None)
        if header is None:
            continue
        fields = header.split()
        try:
            # the row ends "... <seq from> <seq to> <strand> <trunc> <acc> <trunc> <gc>",
            # counted from the right because the leading columns differ between a
            # covariance-model hit and an HMM-only one
            accuracy = float(fields[-3])
            key = (query, model, fields[-7], fields[-6])
        except (IndexError, ValueError) as exc:
            logger.debug("Skipping an unparsable cmscan alignment header: %s", exc)
            continue
        alignments[key] = (
            "".join(line[:-3].strip() for line in block if line.endswith(" CS")),
            accuracy,
        )
    return alignments


@lru_cache(maxsize=4)
def _model_lengths(path: str, _fingerprint: tuple[int, int]) -> dict[str, int]:
    """Map every covariance model's accession and name to its consensus length.

    ``CLEN`` is the model's full width, which is what "full length of feature in db"
    means for a CM hit. It appears only in the model file, never in cmscan's output,
    so the whole file is scanned once per process and cached against its size and
    mtime.

    NOTE: a pressed database follows each covariance model with the HMM filter built
    from it, which repeats the same NAME and ACC but reports ``LENG`` rather than
    ``CLEN``. Restarting the key list at every NAME keeps those repeats from being
    carried forward and bound to the *next* model's length.
    """
    lengths: dict[str, int] = {}
    keys: list[str] = []
    try:
        with open(path, "rb") as handle:
            for raw in handle:
                if raw.startswith(b"NAME "):
                    keys = [raw.split()[1].decode()]
                elif raw.startswith(b"ACC "):
                    keys.append(raw.split()[1].decode())
                elif raw.startswith(b"CLEN "):
                    lengths.update(dict.fromkeys(keys, int(raw.split()[1])))
                    keys = []
    except (OSError, IndexError, ValueError) as exc:
        logger.debug("Could not read consensus lengths from %s: %s", path, exc)
        return {}
    return lengths


def model_lengths(path: str | Path) -> dict[str, int]:
    """Return the accession/name to consensus-length map for a CM database.

    The cached map is copied out so a caller mutating the result cannot corrupt every
    later lookup in the process.
    """
    try:
        stat = Path(path).stat()
    except OSError:
        return {}
    return dict(_model_lengths(str(path), (stat.st_size, stat.st_mtime_ns)))


def parse_output(
    path: str | Path,
    alignment_path: str | Path | None = None,
    cm_path: str | Path | None = None,
) -> pd.DataFrame:
    """Parse Infernal ``--tblout --fmt 2`` output into candidate hits."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("Infernal output is missing its column headers")

    widths = [len(field) + 1 for field in lines[1].split()]
    ends = list(accumulate(widths))
    ends[-1] += 100
    starts = [0, *ends[:-1]]
    positions = list(zip(starts, ends, strict=True))
    names = [lines[0][start:end].strip() for start, end in positions]

    try:
        dataframe = pd.read_fwf(path, comment="#", colspecs=positions, header=None)
        dataframe.columns = names
    except pd.errors.EmptyDataError:
        dataframe = pd.DataFrame(columns=names)

    dataframe = dataframe[REQUIRED_COLUMNS]
    dataframe = dataframe.loc[:, ~dataframe.columns.duplicated()]
    dataframe = dataframe.rename(
        columns={
            "target name": "name",
            "query name": "qseqid",
            "seq from": "qstart",
            "seq to": "qend",
            "mdl from": "sstart",
            "mdl to": "send",
            "E-value": "evalue",
            "strand": "sframe",
            "description of target": "blurb",
        }
    )

    # The alignment report shares only the model name and the raw sequence coordinates
    # with the table, so the join key is taken before the columns below are reformatted
    # for display: the name keeps its underscores and the coordinates stay in cmscan's
    # own order, which runs descending on the minus strand.
    raw_name = dataframe["name"].astype(str)
    raw_from = _coordinate_text(dataframe["qstart"])
    raw_to = _coordinate_text(dataframe["qend"])

    dataframe["accession"] = dataframe["accession"].str.replace("-", " ")
    dataframe["clan name"] = dataframe["clan name"].str.replace("-", " ")
    dataframe["name"] = dataframe["name"].str.replace("_", " ")
    # Identify each hit by its stable Rfam accession (e.g. RF00162) instead of
    # cmscan's per-run hit ordinal, so the id is meaningful and reproducible across
    # runs. Fall back to the model name on the rare model without an accession.
    clean_accession = dataframe["accession"].str.strip()
    dataframe["sseqid"] = clean_accession.where(
        clean_accession != "", dataframe["name"]
    )
    dataframe["blurb"] = (
        "Accession: " + dataframe["accession"] + " - " + dataframe["blurb"]
    )
    dataframe["type"] = "ncRNA"
    dataframe["qseq"] = ""
    # cmscan reports covariance-model hits, which have no base-by-base traceback to
    # summarize; the column exists so every source shares one schema.
    dataframe["btop"] = ""

    alignments = parse_alignments(alignment_path) if alignment_path else {}
    matched = [
        alignments.get(key)
        for key in zip(
            dataframe["qseqid"].astype(str), raw_name, raw_from, raw_to, strict=True
        )
    ]
    if alignments and any(hit is None for hit in matched):
        # every table row should have a block; a gap means the two outputs disagree,
        # which is worth surfacing because the identity below then falls back
        logger.warning(
            "cmscan reported %d hit(s) with no alignment block; "
            "their identity falls back to 100",
            sum(hit is None for hit in matched),
        )
    # The consensus secondary structure in WUSS notation: what a covariance model
    # actually matched on, and the closest thing a CM hit has to a traceback.
    # NOTE: this is indexed by alignment column, not by model position -- cmscan's CS
    # line includes the insertion columns, so it can be longer than the model span.
    dataframe["structure"] = ["" if hit is None else hit[0] for hit in matched]

    coordinates = dataframe[["qstart", "qend"]].apply(pd.to_numeric)
    dataframe["qstart"] = coordinates.min(axis=1).astype("int64")
    dataframe["qend"] = coordinates.max(axis=1).astype("int64")
    dataframe["sframe"] = dataframe["sframe"].map({"-": -1, "+": 1})
    dataframe["length"] = abs(dataframe["qend"] - dataframe["qstart"]) + 1

    # CLEN is the model's full width; the aligned span covers only the matched part, so
    # measuring against it would report every partial hit as a full-length match.
    lengths = model_lengths(cm_path) if cm_path else {}
    aligned_span = abs(dataframe["send"] - dataframe["sstart"]) + 1
    dataframe["slen"] = [
        lengths.get(accession) or lengths.get(name) or span
        for accession, name, span in zip(
            dataframe["sseqid"].astype(str), raw_name, aligned_span, strict=True
        )
    ]
    # Infernal scores a covariance model on how well a sequence fits the family's
    # structure, not on base identity: a canonical 5S rRNA is only ~64% identical to
    # its consensus. Reporting that as identity would rank every structured RNA far
    # below its real quality, so /identity carries the alignment's average posterior
    # probability instead -- the model's own confidence, already on a 0-100 scale.
    # A hit whose alignment could not be read keeps the historical 100, which is what
    # every Rfam hit reported before accuracy was available.
    dataframe["pident"] = [
        100.0 if hit is None else round(hit[1] * 100, 1) for hit in matched
    ]
    # That identity is a confidence, not a count of matching bases, so it must not earn
    # the exact-match bonus: accuracy is printed to two decimals, and rewarding 1.00
    # while passing over 0.99 would be a tenfold score step off a rounding boundary.
    dataframe["sequence_identity"] = False
    return dataframe.drop(columns=["#idx", "accession", "clan name"])
