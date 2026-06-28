"""Rotate a circular plasmid to a canonical origin-of-replication frame.

The public entry point is :func:`rotate_to_origin`. It detects the origin of
replication with a cheap snapgene-only search, then rotates (and, if needed,
reverse-complements) the raw sequence so the chosen ori begins at base 1 on the
forward strand. Annotation then runs once, downstream, on the rotated sequence,
so no feature-coordinate transform is needed here.

When no qualifying bacterial ori is found, the sequence is placed in a
deterministic, rotation- and strand-invariant frame via a min-hash k-mer rule so
that the same circular molecule always yields the same output.
"""

import logging
import zlib
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from . import _package_data

logger = logging.getLogger(__name__)

# k-mer width for the deterministic fallback. Clamped to the sequence length for
# inputs shorter than this; an odd width >10 makes hash ties rare in practice and
# avoids the self-overlap symmetries of an even k on repetitive sequence.
_FALLBACK_K = 17


@dataclass(frozen=True)
class RotationResult:
    """Outcome of a rotation: the new sequence and how it was derived."""

    rotated_seq: str
    offset: int  # rotation applied, expressed in the final (possibly flipped) frame
    flipped: bool  # whole plasmid reverse-complemented to put the ori on the top strand
    ori_name: str | None  # snapgene name of the chosen ori, else None
    rank: int | None  # ranking position of the chosen ori, else None
    fallback_used: bool  # True when the deterministic minimizer fallback was used

    def as_seq(self) -> Seq:
        """Return the rotated sequence as a Biopython :class:`~Bio.Seq.Seq`."""
        return Seq(self.rotated_seq)


def _load_ori_ranking() -> pd.DataFrame:
    """Load the curated, prevalence-ranked ori table shipped with the package.

    Columns: ``name`` (matches the snapgene name), ``rank`` (prevalence order,
    populated for bacterial origins only), and ``type`` (origin category, e.g.
    ``bacterial``, ``phage``, ``viral``, ``yeast``, ``fungal``).
    """
    path = _package_data.get_resource("data", "ori_ranking.csv")
    ranking = pd.read_csv(path)
    ranking["type"] = ranking["type"].astype(str).str.strip()
    return ranking


def _detect_origins(seq: str, yaml_file: Path | None, cores: int) -> pd.DataFrame:
    """Run a snapgene-only search and return non-fragment rep_origin candidates."""
    # Imported lazily so importing this module does not pull in the search engines.
    from . import _filter, annotate

    sources = _package_data.get_yaml(yaml_file or _package_data.get_yaml_path())
    config = sources.get("snapgene")
    if config is None:
        logger.warning("No 'snapgene' source configured; cannot detect origins")
        return pd.DataFrame()

    hits = annotate._collect_source_hits(
        seq, "snapgene", config, is_linear=False, threads=cores
    )
    if hits.empty or "type" not in hits.columns:
        return pd.DataFrame()

    hits = hits.loc[hits["type"] == "rep_origin"].copy()
    if hits.empty:
        return hits

    # filter_and_clean_hits scores hits (percmatch/pi_permatch), wraps circular
    # coordinates, and collapses the doubled-query duplicates by same-kind overlap.
    hits["kind"] = 1
    hits = _filter.filter_and_clean_hits(hits, is_linear=False)
    if hits.empty:
        return hits

    hits["fragment"] = hits.apply(annotate._is_fragment, axis=1)
    hits["qend"] = hits["qend"] + 1  # exclusive end, matching annotate.annotate
    return hits.loc[~hits["fragment"]].reset_index(drop=True)


def _select_origin(candidates: pd.DataFrame, ranking: pd.DataFrame) -> pd.Series | None:
    """Pick the best bacterial ori candidate, or None if none are ranked."""
    bacterial = ranking.loc[ranking["type"] == "bacterial"]
    rank_by_name = dict(zip(bacterial["name"], bacterial["rank"], strict=True))

    scored = candidates.copy()
    scored["rank"] = scored["name"].map(rank_by_name)
    scored = scored.dropna(subset=["rank"])
    if scored.empty:
        return None

    # deterministic: best (lowest) rank, then earliest position, then forward strand
    scored = scored.sort_values(
        ["rank", "qstart", "sframe"], ascending=[True, True, False]
    )
    return scored.iloc[0]


def _rotate_to_anchor(
    seq: str, qstart: int, qend: int, sframe: int, n: int
) -> tuple[str, int, bool]:
    """Rotate so the feature starts at base 1 on the forward strand.

    ``qend`` is exclusive (as produced by annotate). For a forward feature the
    rotation anchor is ``qstart``. For a reverse feature the biological 5' end is
    forward index ``qend - 1``, which maps to index ``n - qend`` after the whole
    sequence is reverse-complemented; rotating there leaves the ori reading 5'->3'
    from base 1. Both branches handle origin-spanning features (qstart > qend)
    because the anchor is a single coordinate, not the interval.
    """
    if sframe >= 0:
        offset = qstart % n
        return seq[offset:] + seq[:offset], offset, False

    reverse = str(Seq(seq).reverse_complement())
    offset = (n - qend) % n
    return reverse[offset:] + reverse[:offset], offset, True


def _canonical_fallback(seq: str) -> tuple[str, int, bool]:
    """Place the sequence in a deterministic, rotation/strand-invariant frame.

    The candidate frames are every rotation of ``seq`` and of its reverse
    complement. Each is named by ``(strand, start index)`` and scored by the
    crc32 of the k-mer beginning at that index. The minimum is chosen; ties are
    broken by the full rotation string (guaranteeing a unique result even under
    hash collisions) and finally by strand (resolving true palindromes).

    Invariance: rotating the input only renames the index of a given rotation,
    and the reverse complement of any rotated input is itself a rotation of
    ``RC(seq)`` -- so the candidate *set* is unchanged by how the user entered the
    molecule. The argmin therefore returns the identical string every time.
    """
    n = len(seq)
    if n == 0:
        return seq, 0, False

    reverse = str(Seq(seq).reverse_complement())
    k = min(_FALLBACK_K, n)

    # pass 1: find the minimal k-mer hash across both strands
    best_hash: int | None = None
    tied: list[tuple[bool, int]] = []
    for flipped, strand in ((False, seq), (True, reverse)):
        doubled = strand + strand
        for i in range(n):
            digest = zlib.crc32(doubled[i : i + k].encode("ascii"))
            if best_hash is None or digest < best_hash:
                best_hash = digest
                tied = [(flipped, i)]
            elif digest == best_hash:
                tied.append((flipped, i))

    # pass 2: resolve ties on the full rotation string, then strand
    def rotation(flipped: bool, i: int) -> str:
        strand = reverse if flipped else seq
        return (strand + strand)[i : i + n]

    flipped, index = min(tied, key=lambda fi: (rotation(*fi), fi[0]))
    return rotation(flipped, index), index, flipped


def rotate_to_origin(
    seq: str | Seq, yaml_file: Path | None = None, cores: int = 1
) -> RotationResult:
    """Rotate a circular sequence to a canonical origin-of-replication frame.

    Accepts a plain ``str`` or a Biopython :class:`~Bio.Seq.Seq`. The rotated
    sequence is available as ``result.rotated_seq`` (``str``) or via
    ``result.as_seq()`` (``Seq``).
    """
    seq = str(seq)
    n = len(seq)
    candidates = _detect_origins(seq, yaml_file, cores)
    if not candidates.empty:
        chosen = _select_origin(candidates, _load_ori_ranking())
        if chosen is not None:
            rotated, offset, flipped = _rotate_to_anchor(
                seq,
                int(chosen["qstart"]),
                int(chosen["qend"]),
                int(chosen["sframe"]),
                n,
            )
            logger.info(
                "Rotated to ori %r (rank %d)%s",
                chosen["name"],
                int(chosen["rank"]),
                " on reverse strand (flipped)" if flipped else "",
            )
            return RotationResult(
                rotated_seq=rotated,
                offset=offset,
                flipped=flipped,
                ori_name=str(chosen["name"]),
                rank=int(chosen["rank"]),
                fallback_used=False,
            )

    rotated, offset, flipped = _canonical_fallback(seq)
    logger.info("No qualifying ori found; applied deterministic minimizer rotation")
    return RotationResult(
        rotated_seq=rotated,
        offset=offset,
        flipped=flipped,
        ori_name=None,
        rank=None,
        fallback_used=True,
    )
