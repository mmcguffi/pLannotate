import json
import re
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from functools import lru_cache
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from pandas.testing import assert_frame_equal

from plannotate import resources
from plannotate.models import Construct


ROOT = Path(__file__).resolve().parents[1]
FASTA_DIR = ROOT / "plannotate" / "data" / "fastas"
CONTROL_DIR = Path(__file__).parent / "test_data" / "annotation_controls"
FASTA_PATHS = sorted(FASTA_DIR.glob("*.fa"))
TOOLCHAIN = {
    "blastn": (["blastn", "-version"], r"blastn: ([^+\s]+)"),
    "diamond": (["diamond", "version"], r"diamond version ([^\s]+)"),
    "infernal": (["cmscan", "-h"], r"INFERNAL ([^\s]+)"),
}


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


@dataclass
class CaseResult:
    plasmid: str
    mode: str
    status: str
    control_annotations: int | None
    current_annotations: int | None
    reason: str


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


def canonical_sseqid(value: object) -> str:
    """Ignore pandas' historical 5.0 versus 5 formatting for numeric Rfam IDs."""
    text = str(value)
    try:
        numeric = float(text)
    except ValueError:
        return text
    return str(int(numeric)) if numeric.is_integer() else text


def normalized_csv(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.reset_index(drop=True).copy()
    frame["sseqid"] = frame["sseqid"].map(canonical_sseqid)
    return frame


def read_genbank(text: str):
    return SeqIO.read(StringIO(text), "genbank")


def feature_signature(feature):
    qualifiers = tuple(
        (key, tuple(value)) for key, value in sorted(feature.qualifiers.items())
    )
    return str(feature.location), feature.type, qualifiers


def feature_label(signature):
    location, feature_type, qualifiers = signature
    qualifier_map = dict(qualifiers)
    label = qualifier_map.get("label", ("unlabelled",))[0]
    return f"{label} ({feature_type}, {location})"


@lru_cache(maxsize=1)
def installed_tool_versions():
    versions = {}
    for name, (command, pattern) in TOOLCHAIN.items():
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            versions[name] = f"unavailable: {exc}"
            continue
        match = re.search(pattern, result.stdout + result.stderr)
        versions[name] = match.group(1) if match else "unknown"
    return versions


@lru_cache(maxsize=1)
def context_changes():
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
    installed_tools = installed_tool_versions()
    expected_tools = context["search_tools"]
    if installed_tools != expected_tools:
        changes.append(
            "search-tool versions differ: "
            f"expected {expected_tools}, found {installed_tools}"
        )
    return tuple(changes)


def compare_csv(actual, expected):
    actual = normalized_csv(actual)
    expected = normalized_csv(expected)
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


def compare_genbank(actual, expected):
    changes = []
    if actual.seq != expected.seq:
        changes.append("sequence changed")
    if actual.annotations.get("topology") != expected.annotations.get("topology"):
        changes.append(
            "topology changed from "
            f"{expected.annotations.get('topology')} to "
            f"{actual.annotations.get('topology')}"
        )

    actual_features = [feature_signature(feature) for feature in actual.features]
    expected_features = [feature_signature(feature) for feature in expected.features]
    if actual_features != expected_features:
        actual_counts = Counter(actual_features)
        expected_counts = Counter(expected_features)
        added = sorted(
            feature_label(item) for item in (actual_counts - expected_counts).elements()
        )
        removed = sorted(
            feature_label(item) for item in (expected_counts - actual_counts).elements()
        )
        changes.append(
            f"features changed ({len(expected_features)} control, "
            f"{len(actual_features)} current); added={added}; removed={removed}"
        )
    return "; ".join(changes) or None


def evaluate_case(case: ControlCase):
    expected_csv = pd.read_csv(f"{case.control_stem}.csv")
    expected_gbk = SeqIO.read(f"{case.control_stem}.gbk", "genbank")
    sequence = SeqIO.read(case.fasta_path, "fasta").seq
    construct = Construct(
        seq=sequence,
        linear=case.linear,
        detailed=case.detailed,
    )
    actual_csv = construct.to_csv()
    actual_gbk_text = construct.to_genbank()
    actual_gbk = read_genbank(actual_gbk_text)
    changes = [
        change
        for change in (
            compare_csv(actual_csv, expected_csv),
            compare_genbank(actual_gbk, expected_gbk),
        )
        if change
    ]
    result = CaseResult(
        plasmid=case.fasta_path.stem,
        mode=case.mode,
        status="changed" if changes else "passed",
        control_annotations=len(expected_csv),
        current_annotations=len(actual_csv),
        reason=" | ".join(changes) if changes else "CSV and GenBank match",
    )
    return result, actual_csv, actual_gbk_text


def error_result(case: ControlCase, exc: Exception):
    return CaseResult(
        plasmid=case.fasta_path.stem,
        mode=case.mode,
        status="error",
        control_annotations=None,
        current_annotations=None,
        reason=f"{type(exc).__name__}: {exc}",
    )


def _table_text(value: object, limit: int = 180):
    text = str(value).replace("|", "\\|").replace("\n", " ")
    return text if len(text) <= limit else text[: limit - 1] + "…"


def _result_summary(result: CaseResult):
    if result.status == "passed":
        return result.reason
    if result.status == "error":
        return _table_text(result.reason)
    changes = []
    if result.control_annotations != result.current_annotations:
        changes.append(
            f"annotation count {result.control_annotations} → {result.current_annotations}"
        )
    if "CSV changed" in result.reason:
        changes.append("CSV rows or metadata differ")
    if "features changed" in result.reason:
        changes.append("GenBank features or qualifiers differ")
    if "sequence changed" in result.reason:
        changes.append("sequence differs")
    if "topology changed" in result.reason:
        changes.append("topology differs")
    return "; ".join(changes) or _table_text(result.reason)


def render_markdown_report(results, context_warnings=()):
    counts = Counter(result.status for result in results)
    lines = [
        "# Annotation control comparison",
        "",
        f"- Total: {len(results)}",
        f"- Passed: {counts['passed']}",
        f"- Changed: {counts['changed']}",
        f"- Errors: {counts['error']}",
        "- Annotation changes are informational unless strict mode is requested.",
        "",
    ]
    if context_warnings:
        lines.extend(
            [
                "## Context warnings",
                "",
                *[f"- {warning}" for warning in context_warnings],
                "",
            ]
        )
    lines.extend(
        [
            "## Results",
            "",
            "| Plasmid | Mode | Result | Control | Current | Why |",
            "| --- | --- | --- | ---: | ---: | --- |",
        ]
    )
    for result in results:
        control = (
            result.control_annotations
            if result.control_annotations is not None
            else "—"
        )
        current = (
            result.current_annotations
            if result.current_annotations is not None
            else "—"
        )
        lines.append(
            f"| {result.plasmid} | {result.mode} | {result.status.upper()} | "
            f"{control} | {current} | {_result_summary(result)} |"
        )

    changed = [result for result in results if result.status != "passed"]
    if changed:
        lines.extend(["", "## Change details", ""])
        for result in changed:
            lines.extend(
                [
                    f"### {result.plasmid} — {result.mode}",
                    "",
                    result.reason,
                    "",
                ]
            )
    return "\n".join(lines).rstrip() + "\n"


def write_report_artifact(output_dir: Path, results, context_warnings=()):
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).isoformat()
    counts = Counter(result.status for result in results)
    payload = {
        "generated_at": generated_at,
        "context_warnings": list(context_warnings),
        "summary": {
            "total": len(results),
            "passed": counts["passed"],
            "changed": counts["changed"],
            "errors": counts["error"],
        },
        "results": [asdict(result) for result in results],
    }
    (output_dir / "report.json").write_text(json.dumps(payload, indent=2) + "\n")
    (output_dir / "report.md").write_text(
        render_markdown_report(results, context_warnings)
    )
