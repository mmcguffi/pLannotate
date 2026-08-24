"""Conservative policy for nested and composite-reference fragment annotations.

The general nesting rules deliberately fail open: only a fragment with clear evidence
of being an incidental contained match receives ``suppress_child``. A separate unary
rule removes fragment labels whose subject interval is almost entirely explained by a
curated component of a composite reference. Whole features, structured RNAs,
compound-feature components, near-complete matches, edge-clipped matches, and
high-confidence CDS-derived sequence survive. Source-specific pair exceptions and
component intervals and recurring low-specificity fragment regions live in packaged
CSVs so they can be reviewed without editing executable code.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from math import ceil, isfinite
from typing import NamedTuple, cast

import pandas as pd

from . import _curation, _package_data

PairKey = tuple[str, str, str, str]

MIN_DERIVED_CDS_IDENTITY = 95.0
MIN_DERIVED_CDS_LENGTH_NT = 90
MIN_NEAR_COMPLETE_MATCH = 80.0
MIN_NEAR_COMPLETE_IDENTITY = 70.0
MIN_SAME_SOURCE_IDENTITY = 95.0
MIN_EXACT_ELEMENT_MATCH = 30.0
MIN_EXACT_ELEMENT_IDENTITY = 98.0
MAX_STRONG_EVALUE = 1e-10
BOUNDARY_SLOP_NT = 3
MAX_COMPOSITE_INFORMATIVE_NT = 3
MIN_BOUNDARY_OVERHANG_NT = 30
MIN_BOUNDARY_OVERHANG_FRACTION = 0.10
VALID_STATUSES = frozenset({"good", "bad", "review"})
VALID_ACTIONS = frozenset(
    {
        "keep",
        "suppress_child",
        "replace_parent",
        "trim_parent",
        "replace_child",
        "relabel_child",
        "review",
    }
)


class Decision(NamedTuple):
    status: str
    action: str
    rationale: str
    source: str


def as_bool(value: object) -> bool:
    """Interpret CSV-compatible boolean values."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def pair_key(row: Mapping[str, object]) -> PairKey:
    """Return the stable source/accession tuple used by curated overrides."""
    return (
        str(row["parent_db"]),
        str(row["parent_sseqid"]),
        str(row["nested_db"]),
        str(row["nested_sseqid"]),
    )


@lru_cache(maxsize=1)
def curated_decisions() -> dict[PairKey, Decision]:
    """Load and expand the packaged parent/child overrides."""
    frame = pd.read_csv(
        _package_data.get_resource("data", "nested_feature_overrides.csv"),
        dtype=str,
    ).fillna("")
    decisions: dict[PairKey, Decision] = {}
    for _, row in frame.iterrows():
        decision = Decision(
            str(row["status"]).strip(),
            str(row["action"]).strip(),
            str(row["rationale"]).strip(),
            str(row["source"]).strip(),
        )
        if decision.status not in VALID_STATUSES:
            raise ValueError(f"Unknown nested-feature status: {decision.status!r}")
        if decision.action not in VALID_ACTIONS:
            raise ValueError(f"Unknown nested-feature action: {decision.action!r}")
        if not decision.rationale or not decision.source:
            raise ValueError("Nested-feature overrides require rationale and source")
        for parent in str(row["parent_sseqid"]).split(";"):
            for child in str(row["child_sseqid"]).split(";"):
                key = (
                    str(row["parent_db"]).strip(),
                    parent.strip(),
                    str(row["child_db"]).strip(),
                    child.strip(),
                )
                if not all(key):
                    raise ValueError(
                        f"Nested-feature override has a blank key component: {key!r}"
                    )
                if key in decisions:
                    raise ValueError(f"Duplicate nested-feature override: {key!r}")
                decisions[key] = decision
    return decisions


def _number(row: Mapping[str, object], key: str, default: float = 0.0) -> float:
    value = row.get(key)
    if value is None:
        return default
    try:
        number = float(str(value))
    except (TypeError, ValueError):
        return default
    return default if number != number else number


def _strong_evalue(row: Mapping[str, object]) -> bool:
    """Return whether the row has strong statistical support, if recorded."""
    value = row.get("evalue")
    if value is None:
        # Callers predating the audit e-value column can still exercise deterministic
        # identity/coverage rules. Fresh annotation and reports always provide it.
        return True
    try:
        evalue = float(str(value))
    except (TypeError, ValueError):
        return True
    return evalue != evalue or evalue <= MAX_STRONG_EVALUE


def _coordinate(value: object) -> int:
    """Return an integral subject coordinate, or zero for missing/invalid input."""
    try:
        number = float(str(value))
    except (TypeError, ValueError):
        return 0
    return int(number) if isfinite(number) else 0


def _composite_fragment_decision(
    database: object,
    accession: object,
    subject_start: object,
    subject_end: object,
    fragment: object,
) -> Decision | None:
    """Suppress a fragment explained by curated components of its source record."""
    if not as_bool(fragment):
        return None
    regions = _curation.composite_reference_regions().get(
        (str(database).strip(), str(accession).strip()), ()
    )
    if not regions:
        return None
    start = _coordinate(subject_start)
    end = _coordinate(subject_end)
    if min(start, end) < 1:
        return None
    interval_start, interval_end = sorted((start, end))
    covered = 0
    cursor = interval_start
    used = []
    for region in regions:
        overlap_start = max(cursor, region.start)
        overlap_end = min(interval_end, region.end)
        if overlap_start > overlap_end:
            continue
        covered += overlap_end - overlap_start + 1
        cursor = overlap_end + 1
        used.append(region)
        if cursor > interval_end:
            break
    informative = interval_end - interval_start + 1 - covered
    if not used or informative > MAX_COMPOSITE_INFORMATIVE_NT:
        return None
    sources = ";".join(dict.fromkeys(region.source for region in used))
    components = ", ".join(
        dict.fromkeys(
            f"{region.component_db}:{region.component_sseqid}" for region in used
        )
    )
    return Decision(
        "bad",
        "suppress_child",
        f"The fragment has only {informative} aligned subject base(s) outside "
        f"curated embedded component region(s) ({components}), so it supports the "
        "component rather than the larger source-record label.",
        sources,
    )


def _curated_fragment_suppression_decision(
    database: object,
    accession: object,
    subject_length: object,
    subject_start: object,
    subject_end: object,
    percent_identity: object,
    fragment: object,
) -> Decision | None:
    """Suppress a weak fragment in a manually adjudicated source interval."""
    if not as_bool(fragment):
        return None
    regions = _curation.fragment_suppression_regions().get(
        (str(database).strip(), str(accession).strip()), ()
    )
    if not regions:
        return None
    length = _coordinate(subject_length)
    start = _coordinate(subject_start)
    end = _coordinate(subject_end)
    try:
        identity = float(str(percent_identity))
    except (TypeError, ValueError):
        return None
    if min(length, start, end) < 1 or not isfinite(identity):
        return None
    interval_start, interval_end = sorted((start, end))
    for region in regions:
        if length != region.subject_length or identity > region.max_identity:
            continue
        overlap = max(
            0,
            min(interval_end, region.end) - max(interval_start, region.start) + 1,
        )
        outside = interval_end - interval_start + 1 - overlap
        if overlap and outside <= BOUNDARY_SLOP_NT:
            return Decision(
                "bad",
                "suppress_child",
                region.rationale,
                region.source,
            )
    return None


def _is_near_complete(row: Mapping[str, object]) -> bool:
    """Keep a fragment that covers most of its reference with credible support."""
    if _number(row, "percent_match") < MIN_NEAR_COMPLETE_MATCH:
        return False
    identity = _number(row, "percent_identity")
    same_source = str(row.get("parent_db", "")) == str(row.get("nested_db", ""))
    return (same_source and identity >= MIN_SAME_SOURCE_IDENTITY) or (
        identity >= MIN_NEAR_COMPLETE_IDENTITY and _strong_evalue(row)
    )


def _is_boundary_extension(row: Mapping[str, object]) -> bool:
    """Identify a strong child whose unaligned reference tail crosses a parent edge."""
    parent_length = int(_number(row, "parent_length"))
    subject_length = int(_number(row, "nested_subject_length"))
    subject_start = int(_number(row, "nested_subject_start"))
    subject_end = int(_number(row, "nested_subject_end"))
    if min(parent_length, subject_length, subject_start, subject_end) <= 0:
        return False

    prefix = min(subject_start, subject_end) - 1
    suffix = subject_length - max(subject_start, subject_end)
    strand = int(_number(row, "nested_strand", 1))
    touches_left = int(_number(row, "nested_start", parent_length)) <= BOUNDARY_SLOP_NT
    touches_right = int(_number(row, "nested_end")) >= parent_length - BOUNDARY_SLOP_NT
    # Query direction reverses which reference tail would continue past an edge.
    left_tail = prefix if strand >= 0 else suffix
    right_tail = suffix if strand >= 0 else prefix
    required_tail = max(
        MIN_BOUNDARY_OVERHANG_NT,
        ceil(subject_length * MIN_BOUNDARY_OVERHANG_FRACTION),
    )
    extends_outside = (touches_left and left_tail >= required_tail) or (
        touches_right and right_tail >= required_tail
    )
    return (
        extends_outside
        and int(_number(row, "nested_length")) >= MIN_DERIVED_CDS_LENGTH_NT
        and _number(row, "percent_identity") >= MIN_DERIVED_CDS_IDENTITY
        and _strong_evalue(row)
    )


def classify_row(row: Mapping[str, object]) -> Decision:
    """Apply ordered nested-feature policy rules to one parent/child pair."""
    composite = _composite_fragment_decision(
        row.get("nested_db", ""),
        row.get("nested_sseqid", ""),
        row.get("nested_subject_start"),
        row.get("nested_subject_end"),
        row.get("fragment", False),
    )
    if composite is not None:
        return composite
    curated_fragment = _curated_fragment_suppression_decision(
        row.get("nested_db", ""),
        row.get("nested_sseqid", ""),
        row.get("nested_subject_length"),
        row.get("nested_subject_start"),
        row.get("nested_subject_end"),
        row.get("percent_identity"),
        row.get("fragment", False),
    )
    if curated_fragment is not None:
        return curated_fragment
    curated = curated_decisions().get(pair_key(row))
    if curated is not None:
        return curated
    if str(row.get("nested_type", "")) == "ncRNA":
        return Decision(
            "good",
            "keep",
            "A structured ncRNA annotation can legitimately be nested inside "
            "another sequence feature.",
            "rule:ncrna_keep",
        )
    if (
        not as_bool(row.get("fragment", False))
        and str(row.get("parent_type", "")) == "gene"
        and str(row.get("nested_type", ""))
        in {"CDS", "promoter", "terminator", "intron", "polyA_signal"}
    ):
        return Decision(
            "good",
            "keep",
            "A whole functional component inside a gene-level container is part "
            "of a compound cassette, as in the MX marker families.",
            "rule:compound_gene_component",
        )
    if as_bool(row.get("fragment", False)):
        if str(row.get("parent_type", "")) == str(row.get("nested_type", "")):
            return Decision(
                "review",
                "review",
                "This fragment has the same annotation kind as its parent, so it "
                "does not satisfy the automatic different-kind suppression rule.",
                "rule:same_kind_fragment_review",
            )
        if (
            str(row.get("nested_type", "")) != "CDS"
            and _number(row, "percent_identity") >= MIN_EXACT_ELEMENT_IDENTITY
            and _number(row, "percent_match") >= MIN_EXACT_ELEMENT_MATCH
        ):
            return Decision(
                "good",
                "keep",
                "An exact or near-exact match covers at least 30% of a non-CDS "
                "reference. Preserve short functional elements such as promoters, "
                "operators, recombination sites, repeats, and small RNAs.",
                "rule:exact_short_element_keep",
            )
        if _is_near_complete(row):
            return Decision(
                "good",
                "keep",
                "At least 80% of the child reference is aligned with credible "
                "identity and statistical support. Keep it explicitly labeled as "
                "a fragment rather than discarding genuine feature provenance.",
                "rule:near_complete_fragment_keep",
            )
        if _is_boundary_extension(row):
            return Decision(
                "good",
                "keep",
                "A strong alignment reaches a parent boundary and the unaligned "
                "tail of its reference continues past that edge. Treat it as an "
                "edge-clipped feature rather than contained fragment noise.",
                "rule:boundary_extension_keep",
            )
        if (
            str(row.get("nested_type", "")) == "CDS"
            and _number(row, "percent_identity") >= MIN_DERIVED_CDS_IDENTITY
            and int(_number(row, "nested_length")) >= MIN_DERIVED_CDS_LENGTH_NT
            and _strong_evalue(row)
        ):
            return Decision(
                "good",
                "keep",
                "A translated alignment of at least 30 amino acids at 95% or "
                "greater identity is strong evidence of genuine CDS-derived "
                "sequence. Keep it labeled as a fragment; do not infer a functional "
                "full-length protein.",
                "rule:high_confidence_cds_fragment",
            )
        return Decision(
            "bad",
            "suppress_child",
            "Low-coverage contained fragment of a larger reference feature with no "
            "independent structural, boundary, or high-confidence CDS evidence.",
            "rule:contained_fragment",
        )
    return Decision(
        "good",
        "keep",
        "This child is sufficiently complete to be treated as a feature call. Whole "
        "nested relationships survive by default; functional uncertainty or display "
        "redundancy must be handled without discarding the sequence match.",
        "rule:whole_child_keep",
    )


def classify_report(report: pd.DataFrame) -> pd.DataFrame:
    """Return the audit with policy decision columns appended."""
    rows = cast(list[dict[str, object]], report.to_dict("records"))
    decisions = [classify_row(row) for row in rows]
    classified = report.copy()
    classified["rule_status"] = [decision.status for decision in decisions]
    classified["suggested_action"] = [decision.action for decision in decisions]
    classified["decision_rationale"] = [decision.rationale for decision in decisions]
    classified["decision_source"] = [decision.source for decision in decisions]
    return classified


def _interval_length(start: int, end: int, sequence_length: int) -> int:
    """Return a half-open feature length on a circular coordinate system."""
    length = (end - start) % sequence_length
    return sequence_length if length == 0 and start != end else length


def _relative_containment(
    parent_start: int,
    parent_end: int,
    child_start: int,
    child_end: int,
    sequence_length: int,
) -> tuple[int, int, int] | None:
    """Return child start/end relative to a strict containing parent, if any."""
    parent_length = _interval_length(parent_start, parent_end, sequence_length)
    child_length = _interval_length(child_start, child_end, sequence_length)
    relative_start = (child_start - parent_start) % sequence_length
    relative_end = relative_start + child_length
    if child_length <= 0 or parent_length <= child_length:
        return None
    if relative_end > parent_length:
        return None
    return relative_start, relative_end, parent_length


def suppress_nested_fragments(hits: pd.DataFrame) -> pd.DataFrame:
    """Drop only detailed-mode children unanimously classified ``suppress_child``.

    Parents must be whole, score at least as well as the child, and strictly contain
    it. If multiple credible parents contain a child, any keep/review decision wins;
    this fail-open behavior prevents one broad annotation from erasing a legitimate
    nested feature. Curated ``replace_parent``/``replace_child`` actions are report
    guidance and do not silently mutate accessions during annotation.
    """
    if hits.empty or len(hits) < 2:
        return hits

    sequence_length = int(hits["qlen"].iloc[0])
    records = cast(list[dict[str, object]], hits.to_dict("records"))
    suppress: set[int] = set()
    for child_index, child in enumerate(records):
        if not as_bool(child.get("fragment", False)):
            continue
        decisions: list[Decision] = []
        for parent_index, parent in enumerate(records):
            if parent_index == child_index or as_bool(parent.get("fragment", False)):
                continue
            if _number(parent, "score") < _number(child, "score"):
                continue
            geometry = _relative_containment(
                int(str(parent["qstart"])),
                int(str(parent["qend"])),
                int(str(child["qstart"])),
                int(str(child["qend"])),
                sequence_length,
            )
            if geometry is None:
                continue
            relative_start, relative_end, parent_length = geometry
            decisions.append(
                classify_row(
                    {
                        "parent_db": parent.get("db", ""),
                        "parent_sseqid": parent.get("sseqid", ""),
                        "parent_type": parent.get("type", ""),
                        "parent_length": parent_length,
                        "nested_db": child.get("db", ""),
                        "nested_sseqid": child.get("sseqid", ""),
                        "nested_type": child.get("type", ""),
                        "nested_strand": child.get("sframe", 1),
                        "nested_start": relative_start,
                        "nested_end": relative_end,
                        "nested_length": child.get("length", 0),
                        "nested_subject_length": child.get("slen", 0),
                        "nested_subject_start": child.get("sstart", 0),
                        "nested_subject_end": child.get("send", 0),
                        "percent_identity": child.get("pident", 0),
                        "percent_match": child.get("abs percmatch", 0),
                        "evalue": child.get("evalue"),
                        "fragment": child.get("fragment", False),
                    }
                )
            )
        if decisions and all(
            decision.action == "suppress_child" for decision in decisions
        ):
            suppress.add(child_index)

    if not suppress:
        return hits
    return hits.drop(index=hits.index[list(sorted(suppress))]).reset_index(drop=True)


def suppress_uninformative_composite_fragments(hits: pd.DataFrame) -> pd.DataFrame:
    """Drop fragment labels explained by curated components of their references."""
    if hits.empty:
        return hits
    records = cast(list[dict[str, object]], hits.to_dict("records"))
    suppress = [
        index
        for index, row in enumerate(records)
        if _composite_fragment_decision(
            row.get("db", ""),
            row.get("sseqid", ""),
            row.get("sstart"),
            row.get("send"),
            row.get("fragment", False),
        )
        is not None
    ]
    if not suppress:
        return hits
    return hits.drop(index=hits.index[suppress]).reset_index(drop=True)


def suppress_curated_fragment_artifacts(hits: pd.DataFrame) -> pd.DataFrame:
    """Drop source-region fragment artifacts recorded by manual curation."""
    if hits.empty:
        return hits
    records = cast(list[dict[str, object]], hits.to_dict("records"))
    suppress = [
        index
        for index, row in enumerate(records)
        if _curated_fragment_suppression_decision(
            row.get("db", ""),
            row.get("sseqid", ""),
            row.get("slen"),
            row.get("sstart"),
            row.get("send"),
            row.get("pident"),
            row.get("fragment", False),
        )
        is not None
    ]
    if not suppress:
        return hits
    return hits.drop(index=hits.index[suppress]).reset_index(drop=True)
