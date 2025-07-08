"""Construct dataclass for holding plasmid annotation metadata and providing convenient methods."""

from dataclasses import dataclass, field, fields
from datetime import date
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List, Optional, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from rich.console import _is_jupyter

from . import __version__ as plannotate_version
from . import annotate, bokeh_plot
from . import resources as rsc
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
        return pd.DataFrame(columns=annotate.DF_COLS)

    return pd.DataFrame([feature.to_dict() for feature in features])


@dataclass
class Construct:
    """Holds plasmid annotation metadata and provides convenient methods for interaction."""

    seq: str | Seq
    linear: bool = False
    detailed: bool = False
    # TODO: formalize db_options
    db_options: Path = field(default_factory=rsc.get_yaml_path)
    prior_annotations: Optional[SeqRecord] = None  # holds pre-existing annos
    name: Optional[str] = None
    features: List[Feature] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate and set default values after initialization."""
        self.seq = Seq(self.seq)
        if self.name is None:
            self.name = "construct"

        features: pd.DataFrame = annotate.annotate(
            self.seq,
            self.db_options,
            self.linear,
            self.detailed,
        )
        self.features: List[Feature] = df_to_features(features)

    def __len__(self) -> int:
        """Get the length of the DNA sequence."""
        return len(self.seq)

    @property
    def annotations_df(self) -> pd.DataFrame:
        """Get annotations as pandas DataFrame for compatibility."""
        if not self.features:
            return pd.DataFrame(columns=annotate.DF_COLS)
        return features_to_df(self.features)

    def to_genbank(self) -> str:
        seq_record = self.to_seqrecord()

        outfileloc = NamedTemporaryFile()
        with open(outfileloc.name, "w") as handle:
            seq_record.annotations["molecule_type"] = "DNA"
            SeqIO.write(seq_record, handle, "genbank")
        with open(outfileloc.name) as handle:
            gbk_text = handle.read()
        outfileloc.close()

        return gbk_text

    def to_csv(self) -> pd.DataFrame:
        """Convert to clean CSV DataFrame."""
        CSV_COLS = [
            "sseqid",
            "qstart",
            "qend",
            "sframe",
            "pident",
            "slen",
            "length",
            "abs percmatch",
            "fragment",
            "db",
            "Feature",
            "Type",
            "Description",
            "qseq",
        ]
        REPLACEMENTS = {
            "qstart": "start location",
            "qend": "end location",
            "sframe": "strand",
            "pident": "percent identity",
            "slen": "full length of feature in db",
            "qseq": "sequence",
            "length": "length of found feature",
            "abs percmatch": "percent match length",
            "db": "database",
        }
        return self.annotations_df[CSV_COLS].rename(columns=REPLACEMENTS)

    def plot(self, linear: Optional[bool] = None):
        """Render plot in Jupyter notebook and returns Figure object."""
        try:
            import bokeh.io
            from bokeh.resources import INLINE
        except ImportError:
            logger.warning(
                "Bokeh is not installed. Please install it to use this feature."
            )
            return None

        if linear is None:
            linear = self.linear

        if _is_jupyter():  # for inline plotting in jupyter
            bokeh.io.output_notebook(INLINE)
            bokeh.io.show(bokeh_plot.get_bokeh(self.annotations_df, linear=linear))

        return bokeh_plot.get_bokeh(self.annotations_df, linear=linear)

    def to_html(self, htmlfull: bool = False) -> str:
        """Generate HTML file content from Bokeh plot."""
        try:
            from bokeh.embed import file_html
            from bokeh.resources import CDN, INLINE
        except ImportError:
            logger.warning(
                "Bokeh is not installed. Please install it to use this feature."
            )
            return ""

        bokeh_chart = self.plot()
        bokeh_chart.sizing_mode = "fixed"  # type: ignore # NOTE: version error?

        if htmlfull:
            resource_type = INLINE
        else:
            resource_type = CDN

        if bokeh_chart is not None:
            return file_html(bokeh_chart, resources=resource_type, title="pLannotate")  # type: ignore # noqa: E501
        else:
            return ""

    def to_seqrecord(self) -> SeqRecord:
        """Convert to Biopython SeqRecord object."""
        # this could be passed a more annotated df
        inDf = self.annotations_df.reset_index(drop=True)

        if inDf.empty:
            inDf = pd.DataFrame(columns=annotate.DF_COLS)

        def FeatureLocation_smart(r: pd.Series) -> FeatureLocation:
            # creates compound locations if needed
            if r.qend > r.qstart:
                return FeatureLocation(r.qstart, r.qend, r.sframe)
            elif r.qstart > r.qend:
                first = FeatureLocation(r.qstart, r.qlen, r.sframe)
                second = FeatureLocation(0, r.qend, r.sframe)
                if r.sframe == 1 or r.sframe == 0:
                    return first + second
                elif r.sframe == -1:
                    return second + first
            # fallback: return a zero-length feature at qstart
            return FeatureLocation(r.qstart, r.qstart, r.sframe)

        # adds a FeatureLocation object so it can be used in gbk construction
        inDf["feat loc"] = [FeatureLocation_smart(row) for _, row in inDf.iterrows()]

        # make a record
        record = SeqRecord(seq=Seq(self.seq), name="plasmid")

        record.annotations["data_file_division"] = "SYN"

        if "comment" not in record.annotations:
            record.annotations["comment"] = (
                f"Annotated with pLannotate v{plannotate_version}"
            )
        else:
            record.annotations["comment"] = (
                f"Annotated with pLannotate v{plannotate_version}. {record.annotations['comment']}"
            )

        if "date" not in record.annotations:
            record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()

        if "accession" not in record.annotations:
            record.annotations["accession"] = "."

        if "version" not in record.annotations:
            record.annotations["version"] = "."

        if self.linear:
            record.annotations["topology"] = "linear"
        else:
            record.annotations["topology"] = "circular"

        # this adds "(fragment)" to the end of a feature name
        # if it is a fragment. Maybe a better way show this data in the gbk
        # for downstream analysis, though this may suffice. change type to
        # non-canonical `fragment`?
        def append_frag(row: pd.Series) -> str:
            if row["fragment"] is True:
                return f"{row['Feature']} (fragment)"
            else:
                return f"{row['Feature']}"

        inDf["Feature"] = inDf.apply(lambda x: append_frag(x), axis=1)

        inDf["Type"] = inDf["Type"].str.replace("origin of replication", "rep_origin")
        for index in inDf.index:
            record.features.append(
                SeqFeature(
                    inDf.loc[index]["feat loc"],
                    type=inDf.loc[index]["Type"],  # maybe change 'Type'
                    qualifiers={
                        "note": "pLannotate",
                        "label": inDf.loc[index]["Feature"],
                        "database": inDf.loc[index]["db"],
                        "identity": round(inDf.loc[index]["pident"], 1),
                        "match_length": round(inDf.loc[index]["percmatch"], 1),
                        "fragment": inDf.loc[index]["fragment"],
                        "other": inDf.loc[index]["Type"],
                    },
                )
            )  # maybe change 'Type'

        return record

    def get_fragments(self) -> List[Feature]:
        """Get all fragment features."""
        return [f for f in self.features if f.fragment]

    def get_complete_features(self) -> List[Feature]:
        """Get all complete (non-fragment) features."""
        return [f for f in self.features if not f.fragment]

    @property
    def annotation_count(self) -> int:
        """Get the number of annotations."""
        return len(self.features)

    def summary(self) -> dict:
        """Get a summary of the construct."""
        return {
            "name": self.name,
            "length": len(self),
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
        class_name = self.__class__.__name__
        attrs_str = "\n ".join(
            f"{_.name}: {getattr(self, _.name)!r}" for _ in fields(self) if _.repr
        )
        return f"{class_name}(\n {attrs_str}\n)"
