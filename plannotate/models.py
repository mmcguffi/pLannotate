"""Public domain models for annotated DNA constructs."""

import logging
from copy import deepcopy
from dataclasses import dataclass, field, fields
from datetime import date
from io import StringIO
from pathlib import Path
from typing import Any, cast

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from rich import get_console

from . import __version__ as plannotate_version
from . import _package_data, validation
from ._schema import (
    ANNOTATION_COLUMNS,
    CSV_COLUMN_NAMES,
    CSV_COLUMNS,
)

logger = logging.getLogger(__name__)

# Annotation columns whose Feature field name differs from the column name. Every
# other column in ANNOTATION_COLUMNS maps to an identically-named field.
_COLUMN_TO_FIELD = {
    "db": "database",
    "name": "feature_name",
    "blurb": "description",
    "type": "feature_type",
    "abs percmatch": "abs_percmatch",
}
# Fields with dataclass defaults; absent columns fall back to those defaults.
_OPTIONAL_FIELDS = {"qstart_dup", "qend_dup"}


def _field(column: str) -> str:
    """Return the Feature field name for an annotation column."""
    return _COLUMN_TO_FIELD.get(column, column)


@dataclass
class Feature:
    """A single database annotation on a DNA construct."""

    sseqid: str
    feature_name: str
    description: str
    feature_type: str
    database: str
    qstart: int
    qend: int
    qlen: int
    sstart: int
    send: int
    sframe: int
    qseq: str
    length: int
    slen: int
    pident: float
    percmatch: float
    abs_percmatch: float
    pi_permatch: float
    evalue: float
    score: float
    priority: int
    kind: str | int
    fragment: bool
    wiggle: int
    wstart: int
    wend: int
    qstart_dup: int | None = None
    qend_dup: int | None = None

    def __post_init__(self) -> None:
        if self.qstart_dup is None:
            self.qstart_dup = self.qstart
        if self.qend_dup is None:
            self.qend_dup = self.qend

    @property
    def is_forward_strand(self) -> bool:
        return self.sframe >= 0

    @property
    def is_reverse_strand(self) -> bool:
        return self.sframe < 0

    @property
    def feature_location(self) -> FeatureLocation:
        """Return a simple or origin-spanning Biopython location."""
        if self.qend > self.qstart:
            return FeatureLocation(self.qstart, self.qend, self.sframe)
        if self.qend == self.qstart:
            return FeatureLocation(self.qstart, self.qstart, self.sframe)

        first = FeatureLocation(self.qstart, self.qlen, self.sframe)
        second = FeatureLocation(0, self.qend, self.sframe)
        return first + second if self.is_forward_strand else second + first

    @property
    def seqfeature(self) -> SeqFeature:
        """Return the GenBank representation of this annotation."""
        feature_type = self.feature_type.replace("origin of replication", "rep_origin")
        label = (
            f"{self.feature_name} (fragment)" if self.fragment else self.feature_name
        )
        return SeqFeature(
            self.feature_location,
            type=feature_type,
            qualifiers={
                "note": "pLannotate",
                "label": label,
                "database": self.database,
                "identity": round(self.pident, 1),
                "match_length": round(self.percmatch, 1),
                "fragment": self.fragment,
                "other": feature_type,
            },
        )

    def to_dict(self) -> dict[str, Any]:
        return {column: getattr(self, _field(column)) for column in ANNOTATION_COLUMNS}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Feature":
        values: dict[str, Any] = {}
        for column in ANNOTATION_COLUMNS:
            field = _field(column)
            values[field] = (
                data.get(column) if field in _OPTIONAL_FIELDS else data[column]
            )
        return cls(**values)

    def __str__(self) -> str:
        return f"{self.feature_name} ({self.feature_type}) at {self.qstart}-{self.qend}"

    def __repr__(self) -> str:
        return (
            f"Feature('{self.feature_name}', type='{self.feature_type}', "
            f"pos={self.qstart}-{self.qend}, db='{self.database}')"
        )


def df_to_features(dataframe: pd.DataFrame) -> list[Feature]:
    """Convert valid DataFrame rows into feature objects."""
    features = []
    for index, row in dataframe.iterrows():
        try:
            data = cast(dict[str, Any], row.to_dict())
            features.append(Feature.from_dict(data))
        except (KeyError, TypeError, ValueError) as exc:
            logger.warning("Skipping invalid feature row %s: %s", index, exc)
    return features


def features_to_df(features: list[Feature]) -> pd.DataFrame:
    """Convert feature objects into the canonical annotation DataFrame."""
    if not features:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    return pd.DataFrame(feature.to_dict() for feature in features)


@dataclass
class Construct:
    """A DNA sequence and its pLannotate features."""

    seq: str | Seq
    linear: bool = False
    detailed: bool = False
    db_options: str | Path = field(default_factory=_package_data.get_yaml_path)
    prior_annotations: SeqRecord | None = None
    name: str | None = None
    features: list[Feature] = field(default_factory=list)
    _skip_annotation: bool = False
    cores: int = 1

    def __post_init__(self) -> None:
        self.seq = Seq(self.seq)
        validation.validate_sequence(str(self.seq), max_length=None)
        if self.cores < 1:
            raise ValueError("cores must be at least 1")
        self.db_options = Path(self.db_options)
        self.name = self.name or self._prior_record_name() or "construct"
        if (
            self.prior_annotations is not None
            and self.prior_annotations.seq != self.seq
        ):
            raise ValueError(
                "prior_annotations sequence does not match construct sequence"
            )
        if self._skip_annotation:
            return

        # Imported lazily so importing the domain models does not load search engines.
        from .annotate import annotate

        annotations = annotate(
            self.seq,
            self.db_options,
            self.linear,
            self.detailed,
            self.cores,
        )
        self.features = df_to_features(annotations)

    def _prior_record_name(self) -> str | None:
        if self.prior_annotations is None:
            return None
        name = self.prior_annotations.name
        return None if name == "<unknown name>" else name

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def annotations_df(self) -> pd.DataFrame:
        return features_to_df(self.features)

    def to_genbank(self, prior_annotations: SeqRecord | None = None) -> str:
        output = StringIO()
        SeqIO.write(self.to_seqrecord(prior_annotations), output, "genbank")
        return output.getvalue()

    def to_csv(self) -> pd.DataFrame:
        return self.annotations_df[CSV_COLUMNS].rename(columns=CSV_COLUMN_NAMES)

    def _create_plot(self, linear: bool | None = None) -> Any:
        try:
            from . import bokeh_plot
        except ImportError as exc:
            raise RuntimeError(
                "Plotting requires Bokeh; install pLannotate with the 'plot' extra."
            ) from exc
        return bokeh_plot.get_bokeh(
            self.annotations_df,
            linear=self.linear if linear is None else linear,
        )

    def plot(self, linear: bool | None = None) -> Any:
        """Create a Bokeh figure and display it when running in Jupyter."""
        figure = self._create_plot(linear)

        if get_console().is_jupyter:
            import bokeh.io
            from bokeh.resources import INLINE

            bokeh.io.output_notebook(INLINE)
            bokeh.io.show(figure)
        return figure

    def to_html(self, htmlfull: bool = False) -> str:
        try:
            from bokeh.embed import file_html
            from bokeh.resources import CDN, INLINE
        except ImportError as exc:
            raise RuntimeError(
                "HTML output requires Bokeh; install pLannotate with the 'plot' extra."
            ) from exc

        chart = self._create_plot()
        chart.sizing_mode = "fixed"
        return file_html(
            chart,
            resources=INLINE if htmlfull else CDN,
            title="pLannotate",
        )

    def to_seqrecord(self, prior_annotations: SeqRecord | None = None) -> SeqRecord:
        # an explicit base record overrides the construct's own prior annotations,
        # letting callers layer features onto a different header without mutating self.
        base_record = (
            prior_annotations
            if prior_annotations is not None
            else self.prior_annotations
        )
        record = (
            deepcopy(base_record)
            if base_record is not None
            else SeqRecord(
                seq=Seq(self.seq),
                id=self.name or "construct",
                name=self.name or "construct",
            )
        )
        record.seq = Seq(self.seq)
        annotation_note = f"Annotated with pLannotate v{plannotate_version}"
        existing_comment = record.annotations.get("comment")
        record.annotations.update(
            {
                "accession": record.annotations.get("accession", "."),
                "comment": (
                    f"{annotation_note}. {existing_comment}"
                    if existing_comment
                    else annotation_note
                ),
                "data_file_division": record.annotations.get(
                    "data_file_division", "SYN"
                ),
                "date": record.annotations.get(
                    "date", date.today().strftime("%d-%b-%Y").upper()
                ),
                "molecule_type": "DNA",
                "topology": "linear" if self.linear else "circular",
                "version": record.annotations.get("version", "."),
            }
        )
        record.features.extend(feature.seqfeature for feature in self.features)
        return record

    def get_fragments(self) -> list[Feature]:
        return [feature for feature in self.features if feature.fragment]

    def get_complete_features(self) -> list[Feature]:
        return [feature for feature in self.features if not feature.fragment]

    @property
    def annotation_count(self) -> int:
        return len(self.features)

    def summary(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "length": len(self),
            "topology": "linear" if self.linear else "circular",
            "annotation_count": self.annotation_count,
            "databases_used": sorted({f.database for f in self.features}),
            "feature_types": sorted({f.feature_type for f in self.features}),
            "fragment_count": len(self.get_fragments()),
            "complete_feature_count": len(self.get_complete_features()),
        }

    def __str__(self) -> str:
        summary = self.summary()
        return (
            f"Construct '{summary['name']}' ({summary['length']} bp, "
            f"{summary['topology']}, {summary['annotation_count']} annotations)"
        )

    def __repr__(self) -> str:
        attrs = "\n ".join(
            f"{item.name}: {getattr(self, item.name)!r}"
            for item in fields(self)
            if item.repr
        )
        return f"{type(self).__name__}(\n {attrs}\n)"
