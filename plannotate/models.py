"""Construct dataclass for holding plasmid annotation metadata and providing convenient methods."""

from dataclasses import dataclass, field
from typing import List, Optional, Union

import pandas as pd
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .annotate import annotate
from .logging_config import get_logger

logger = get_logger(__name__)


@dataclass
class Feature:
    """Represents a single annotation feature on a DNA construct."""

    # Core identification
    sseqid: str  # Subject sequence ID from database
    feature_name: str  # Human-readable feature name
    description: str  # Feature description
    feature_type: str  # Type of feature (CDS, promoter, etc.)
    database: str  # Source database

    # Position information
    qstart: int  # Query start position (0-based)
    qend: int  # Query end position (0-based)
    qlen: int  # Query sequence length
    sstart: int  # Subject start position
    send: int  # Subject end position
    sframe: int  # Strand frame (1 for forward, -1 for reverse)

    # Sequence information
    qseq: str  # Query sequence
    length: int  # Feature length
    slen: int  # Subject sequence length

    # Quality metrics
    pident: float  # Percent identity
    percmatch: float  # Percent match
    abs_percmatch: float  # Absolute percent match
    pi_permatch: float  # Priority-weighted percent match
    evalue: float  # E-value
    score: float  # Calculated score

    # Processing information
    priority: int  # Database priority
    kind: Union[str, int]  # Kind of feature (for detailed mode)
    fragment: bool  # Whether this is a fragment

    # Wiggle room for overlap resolution
    wiggle: int  # Wiggle size
    wstart: int  # Wiggle start
    wend: int  # Wiggle end

    # Duplicate positions for origin crossing
    qstart_dup: Optional[int] = None
    qend_dup: Optional[int] = None

    def __post_init__(self) -> None:
        """Set default values and validate data."""
        if self.qstart_dup is None:
            self.qstart_dup = self.qstart
        if self.qend_dup is None:
            self.qend_dup = self.qend

    @property
    def is_forward_strand(self) -> bool:
        """Check if feature is on forward strand."""
        return self.sframe >= 0

    @property
    def is_reverse_strand(self) -> bool:
        """Check if feature is on reverse strand."""
        return self.sframe < 0

    @property
    def spans_origin(self) -> bool:
        """Check if feature spans the origin (qend < qstart)."""
        return self.qend < self.qstart

    @property
    def feature_location(self) -> FeatureLocation:
        """Get Biopython FeatureLocation object."""
        if self.qend > self.qstart:
            return FeatureLocation(self.qstart, self.qend, self.sframe)
        else:
            # Handle origin-spanning features
            first = FeatureLocation(self.qstart, self.qlen, self.sframe)
            second = FeatureLocation(0, self.qend, self.sframe)
            if self.sframe >= 0:
                return first + second
            else:
                return second + first

    @property
    def seqfeature(self) -> SeqFeature:
        """Get Biopython SeqFeature object."""
        return SeqFeature(
            self.feature_location,
            type=self.feature_type,
            qualifiers={
                "note": "pLannotate",
                "label": self.feature_name,
                "database": self.database,
                "identity": round(self.pident, 1),
                "match_length": round(self.percmatch, 1),
                "fragment": self.fragment,
                "other": self.feature_type,
            },
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for DataFrame conversion."""
        return {
            "sseqid": self.sseqid,
            "qstart": self.qstart,
            "qend": self.qend,
            "sstart": self.sstart,
            "send": self.send,
            "sframe": self.sframe,
            "score": self.score,
            "evalue": self.evalue,
            "qseq": self.qseq,
            "length": self.length,
            "slen": self.slen,
            "pident": self.pident,
            "qlen": self.qlen,
            "db": self.database,
            "Feature": self.feature_name,
            "Description": self.description,
            "Type": self.feature_type,
            "priority": self.priority,
            "percmatch": self.percmatch,
            "abs percmatch": self.abs_percmatch,
            "pi_permatch": self.pi_permatch,
            "wiggle": self.wiggle,
            "wstart": self.wstart,
            "wend": self.wend,
            "kind": self.kind,
            "qstart_dup": self.qstart_dup,
            "qend_dup": self.qend_dup,
            "fragment": self.fragment,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Feature":
        """Create Feature from dictionary."""
        return cls(
            sseqid=data["sseqid"],
            feature_name=data["Feature"],
            description=data["Description"],
            feature_type=data["Type"],
            database=data["db"],
            qstart=data["qstart"],
            qend=data["qend"],
            qlen=data["qlen"],
            sstart=data["sstart"],
            send=data["send"],
            sframe=data["sframe"],
            qseq=data["qseq"],
            length=data["length"],
            slen=data["slen"],
            pident=data["pident"],
            percmatch=data["percmatch"],
            abs_percmatch=data["abs percmatch"],
            pi_permatch=data["pi_permatch"],
            evalue=data["evalue"],
            score=data["score"],
            priority=data["priority"],
            kind=data["kind"],
            fragment=data["fragment"],
            wiggle=data["wiggle"],
            wstart=data["wstart"],
            wend=data["wend"],
            qstart_dup=data.get("qstart_dup"),
            qend_dup=data.get("qend_dup"),
        )

    def __str__(self) -> str:
        """String representation."""
        return f"{self.feature_name} ({self.feature_type}) at {self.qstart}-{self.qend}"

    def __repr__(self) -> str:
        """Detailed representation."""
        return (
            f"Feature('{self.feature_name}', type='{self.feature_type}', "
            f"pos={self.qstart}-{self.qend}, db='{self.database}')"
        )


def df_to_features(df: pd.DataFrame) -> List[Feature]:
    """Convert pandas DataFrame to list of Feature objects."""
    if df.empty:
        return []

    features = []
    for _, row in df.iterrows():
        try:
            feature = Feature.from_dict(row.to_dict())
            features.append(feature)
        except Exception as e:
            logger.warning(f"Failed to create Feature from row: {e}")
            continue

    return features


def features_to_df(features: List[Feature]) -> pd.DataFrame:
    """Convert list of Feature objects to pandas DataFrame."""
    if not features:
        return pd.DataFrame(columns=rsc.DF_COLS)

    return pd.DataFrame([feature.to_dict() for feature in features])


@dataclass
class Construct:
    """Holds plasmid annotation metadata and provides convenient methods for interaction."""

    seq: str
    linear: bool = False
    detailed: bool = False
    anno_options: str = field(default_factory=rsc.get_yaml_path)
    prior_features: Optional[List[Feature]] = None  # holds pre-existing annos
    name: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate and set default values after initialization."""
        if self.name is None:
            self.name = "construct"

        features: pd.DataFrame = annotate(
            self.seq,
            self.anno_options,
            self.linear,
            self.detailed,
        )
        self.features: List[Feature] = df_to_features(features)

    @property
    def length(self) -> int:
        """Get the length of the DNA sequence."""
        return len(self.seq)

    @property
    def is_circular(self) -> bool:
        """Check if the construct is circular (not linear)."""
        return not self.linear

    @property
    def has_annotations(self) -> bool:
        """Check if the construct has any annotations."""
        return len(self.features) > 0

    @property
    def annotation_count(self) -> int:
        """Get the number of annotations."""
        return len(self.features)

    @property
    def annotations_df(self) -> pd.DataFrame:
        """Get annotations as pandas DataFrame for compatibility."""
        if not self.features:
            return pd.DataFrame(columns=rsc.DF_COLS)
        return pd.DataFrame([feature.to_dict() for feature in self.features])

    def get_features_by_type(self, feature_type: str) -> List[Feature]:
        """Get all features of a specific type."""
        return [f for f in self.features if f.feature_type == feature_type]

    def get_features_by_database(self, database: str) -> List[Feature]:
        """Get all features from a specific database."""
        return [f for f in self.features if f.database == database]

    def get_fragments(self) -> List[Feature]:
        """Get all fragment features."""
        return [f for f in self.features if f.fragment]

    def get_complete_features(self) -> List[Feature]:
        """Get all complete (non-fragment) features."""
        return [f for f in self.features if not f.fragment]

    def get_high_quality_features(self, min_identity: float = 95.0) -> List[Feature]:
        """Get features with identity above threshold."""
        return [f for f in self.features if f.pident >= min_identity]

    def to_genbank(self): ...

    def to_csv(self): ...

    def plot(self, htmlfull=False): ...

    def to_seqrecord(self) -> SeqRecord:
        """Convert to Biopython SeqRecord object."""
        return rsc.get_seq_record(self.features, self.seq, self.linear)

    def summary(self) -> dict:
        """Get a summary of the construct."""
        return {
            "name": self.name,
            "length": self.length,
            "topology": "linear" if self.linear else "circular",
            "annotation_count": self.annotation_count,
            "databases_used": list(set(f.database for f in self.features)),
            "feature_types": list(set(f.feature_type for f in self.features)),
            "fragment_count": len(self.get_fragments()),
            "complete_feature_count": len(self.get_complete_features()),
        }

    def __str__(self) -> str:
        """String representation of the construct."""
        summary = self.summary()
        return (
            f"Construct '{summary['name']}' ({summary['length']} bp, "
            f"{summary['topology']}, {summary['annotation_count']} annotations)"
        )

    def __repr__(self) -> str:
        """Detailed representation of the construct."""
        return (
            f"Construct(sequence='{self.seq[:50]}{'...' if len(self.seq) > 50 else ''}', "
            f"features={self.annotation_count} features, "
            f"linear={self.linear}, detailed={self.detailed})"
        )
