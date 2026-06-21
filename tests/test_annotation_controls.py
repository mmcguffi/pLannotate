import json
import re
import subprocess
import warnings
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from pandas.testing import assert_frame_equal

from plannotate import resources
from plannotate.models import Construct


pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[1]
FASTA_DIR = ROOT / "plannotate" / "data" / "fastas"
CONTROL_DIR = Path(__file__).parent / "test_data" / "annotation_controls"
FASTA_PATHS = sorted(FASTA_DIR.glob("*.fa"))
TOOLCHAIN = {
    "blastn": (["blastn", "-version"], r"blastn: ([^+\s]+)"),
    "diamond": (["diamond", "version"], r"diamond version ([^\s]+)"),
    "infernal": (["cmscan", "-h"], r"INFERNAL ([^\s]+)"),
}


class AnnotationControlChangedWarning(UserWarning):
    """An annotation differs from its non-blocking master control."""


@dataclass(frozen=True)
class ControlCase:
    mode: str
    fasta_path: Path
    linear: bool = False
    detailed: bool = False

    @property
    def id(self):
        return f"{self.mode}:{self.fasta_path.name}"

    @property
    def control_stem(self):
        return CONTROL_DIR / self.mode / self.fasta_path.stem


CONTROL_CASES = [
    ControlCase(mode, fasta_path, linear=linear, detailed=detailed)
    for fasta_path in FASTA_PATHS
    for mode, linear, detailed in (
        ("regular", False, False),
        ("detailed", False, True),
        ("linear", True, False),
    )
]
CONTROL_CASES.append(
    ControlCase(
        "detailed-linear",
        FASTA_DIR / "pXampl3.fa",
        linear=True,
        detailed=True,
    )
)


def _canonical_sseqid(value: object) -> str:
    """Ignore pandas' historical 5.0 versus 5 formatting for numeric Rfam IDs."""
    text = str(value)
    try:
        numeric = float(text)
    except ValueError:
        return text
    return str(int(numeric)) if numeric.is_integer() else text


def _normalized_csv(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.reset_index(drop=True).copy()
    frame["sseqid"] = frame["sseqid"].map(_canonical_sseqid)
    return frame


def _read_genbank(text: str):
    return SeqIO.read(StringIO(text), "genbank")


def _feature_signature(feature):
    qualifiers = tuple(
        (key, tuple(value)) for key, value in sorted(feature.qualifiers.items())
    )
    return str(feature.location), feature.type, qualifiers


def _feature_label(signature):
    location, feature_type, qualifiers = signature
    qualifier_map = dict(qualifiers)
    label = qualifier_map.get("label", ("unlabelled",))[0]
    return f"{label} ({feature_type}, {location})"


@lru_cache(maxsize=1)
def _installed_tool_versions():
    versions = {}
    for name, (command, pattern) in TOOLCHAIN.items():
        try:
            result = subprocess.run(
                command, check=True, capture_output=True, text=True
            )
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            versions[name] = f"unavailable: {exc}"
            continue
        match = re.search(pattern, result.stdout + result.stderr)
        versions[name] = match.group(1) if match else "unknown"
    return versions


@lru_cache(maxsize=1)
def _context_changes():
    changes = []
    expected_manifest = json.loads(
        (CONTROL_DIR / "current-database-manifest.json").read_text()
    )
    try:
        installed_manifest = resources.get_database_manifest()
    except FileNotFoundError as exc:
        changes.append(str(exc))
    else:
        if installed_manifest != expected_manifest:
            changes.append("installed database manifest differs from the control")

    context = json.loads((CONTROL_DIR / "regression-context.json").read_text())
    installed_tools = _installed_tool_versions()
    expected_tools = context["search_tools"]
    if installed_tools != expected_tools:
        changes.append(
            "search-tool versions differ: "
            f"expected {expected_tools}, found {installed_tools}"
        )
    return tuple(changes)


def _report_change(message, config):
    print(f"\nANNOTATION CONTROL DIFFERENCE\n{message}\n", flush=True)
    warnings.warn(message, AnnotationControlChangedWarning, stacklevel=2)
    if config.getoption("--strict-annotation-controls"):
        pytest.fail(message, pytrace=False)
    pytest.xfail(message)


def _compare_csv(actual, expected):
    actual = _normalized_csv(actual)
    expected = _normalized_csv(expected)
    try:
        assert_frame_equal(
            actual,
            expected,
            check_dtype=False,
            check_exact=False,
            rtol=1e-9,
            atol=1e-9,
        )
    except AssertionError as exc:
        headline = str(exc).splitlines()[0]
        display_columns = [column for column in actual.columns if column != "sequence"]

        def row_counts(frame):
            rows = []
            for row in frame[display_columns].itertuples(index=False, name=None):
                rows.append(
                    tuple(
                        None
                        if pd.isna(value)
                        else round(value, 9)
                        if isinstance(value, float)
                        else value
                        for value in row
                    )
                )
            return Counter(rows)

        actual_counts = row_counts(actual)
        expected_counts = row_counts(expected)

        def format_rows(rows):
            formatted = []
            for row, count in rows.items():
                values = {
                    column: value
                    for column, value in zip(display_columns, row)
                    if value is not None
                }
                prefix = f"{count} x " if count > 1 else ""
                formatted.append(prefix + json.dumps(values, default=str))
            return formatted

        added = format_rows(actual_counts - expected_counts)
        removed = format_rows(expected_counts - actual_counts)
        return (
            f"CSV changed ({len(expected)} control annotations, "
            f"{len(actual)} current annotations): {headline}; "
            f"added rows={added}; removed rows={removed}"
        )
    return None


def _compare_genbank(actual, expected):
    changes = []
    if actual.seq != expected.seq:
        changes.append("sequence changed")
    if actual.annotations.get("topology") != expected.annotations.get("topology"):
        changes.append(
            "topology changed from "
            f"{expected.annotations.get('topology')} to "
            f"{actual.annotations.get('topology')}"
        )

    actual_features = [_feature_signature(feature) for feature in actual.features]
    expected_features = [_feature_signature(feature) for feature in expected.features]
    if actual_features != expected_features:
        actual_counts = Counter(actual_features)
        expected_counts = Counter(expected_features)
        added = sorted(
            _feature_label(item) for item in (actual_counts - expected_counts).elements()
        )
        removed = sorted(
            _feature_label(item) for item in (expected_counts - actual_counts).elements()
        )
        changes.append(
            f"features changed ({len(expected_features)} control, "
            f"{len(actual_features)} current); added={added}; removed={removed}"
        )
    return "; ".join(changes) or None


def test_annotation_control_files_cover_requested_modes():
    fasta_stems = {path.stem for path in FASTA_PATHS}
    for mode in ("regular", "detailed", "linear"):
        mode_dir = CONTROL_DIR / mode
        assert {path.stem for path in mode_dir.glob("*.csv")} == fasta_stems
        assert {path.stem for path in mode_dir.glob("*.gbk")} == fasta_stems

    combined_dir = CONTROL_DIR / "detailed-linear"
    assert {path.stem for path in combined_dir.glob("*.csv")} == {"pXampl3"}
    assert {path.stem for path in combined_dir.glob("*.gbk")} == {"pXampl3"}


def test_annotation_control_context(request):
    if changes := _context_changes():
        _report_change(
            "Annotation control context changed: " + "; ".join(changes),
            request.config,
        )


@pytest.mark.parametrize("case", CONTROL_CASES, ids=lambda case: case.id)
def test_annotation_output_matches_control(case, request):
    if changes := _context_changes():
        warnings.warn(
            "Annotation controls cannot be compared reliably: " + "; ".join(changes),
            AnnotationControlChangedWarning,
            stacklevel=2,
        )
        pytest.skip("annotation control context differs")

    sequence = SeqIO.read(case.fasta_path, "fasta").seq
    construct = Construct(
        seq=sequence,
        linear=case.linear,
        detailed=case.detailed,
    )

    expected_csv = pd.read_csv(f"{case.control_stem}.csv")
    expected_gbk = SeqIO.read(f"{case.control_stem}.gbk", "genbank")
    actual_gbk = _read_genbank(construct.to_genbank())

    changes = [
        change
        for change in (
            _compare_csv(construct.to_csv(), expected_csv),
            _compare_genbank(actual_gbk, expected_gbk),
        )
        if change
    ]
    if changes:
        _report_change(
            f"Annotation control changed for {case.id}: " + " | ".join(changes),
            request.config,
        )
