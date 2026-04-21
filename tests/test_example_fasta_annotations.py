import json
import os
from pathlib import Path

import pytest
from Bio import SeqIO

from plannotate import annotate, resources


ROOT = Path(__file__).resolve().parents[1]
FASTA_DIR = ROOT / "plannotate" / "data" / "fastas"
SNAPSHOT_PATH = Path(__file__).parent / "test_data" / "example_fasta_annotations.json"
UPDATE_SNAPSHOT_ENV = "PLANNOTATE_UPDATE_FASTA_ANNOTATION_SNAPSHOTS"
SNAPSHOT_SCHEMA_VERSION = 2
FASTA_EXTENSIONS = tuple(resources.valid_fasta_exts)
_UPDATED_SNAPSHOT = False
MODE_SNAPSHOT_CASES = [
    {
        "id": "pUC19_linear",
        "fasta": "pUC19.fa",
        "linear": True,
        "is_detailed": False,
    },
    {
        "id": "pUC19_detailed",
        "fasta": "pUC19.fa",
        "linear": False,
        "is_detailed": True,
    },
    {
        "id": "pXampl3_linear_detailed",
        "fasta": "pXampl3.fa",
        "linear": True,
        "is_detailed": True,
    },
]


def _get_example_fasta_paths():
    return sorted(
        path
        for path in FASTA_DIR.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def _normalize_record(record):
    normalized = {}
    for key, value in record.items():
        if isinstance(value, float):
            value = round(value, 6)
        normalized[key] = value
    return normalized


def _get_annotation_records(sequence, *, linear=False, is_detailed=False):
    hits = annotate.annotate(sequence, linear=linear, is_detailed=is_detailed)
    clean_hits = resources.get_clean_csv_df(hits)
    records = json.loads(clean_hits.to_json(orient="records"))
    return [_normalize_record(record) for record in records]


def _build_example_fasta_annotation_snapshot():
    examples = {}

    for fasta_path in _get_example_fasta_paths():
        examples[fasta_path.name] = _build_example_fasta_annotation(fasta_path)

    return {
        "schema_version": SNAPSHOT_SCHEMA_VERSION,
        "examples": examples,
        "mode_examples": _build_mode_annotation_snapshots(),
    }


def _build_example_fasta_annotation(fasta_path, *, linear=False, is_detailed=False):
    sequence = str(SeqIO.read(fasta_path, "fasta").seq)
    annotations = _get_annotation_records(
        sequence,
        linear=linear,
        is_detailed=is_detailed,
    )
    return {
        "sequence_length": len(sequence),
        "linear": linear,
        "is_detailed": is_detailed,
        "annotation_count": len(annotations),
        "annotations": annotations,
    }


def _build_mode_annotation_snapshots():
    return {
        case["id"]: _build_mode_annotation_snapshot(case)
        for case in MODE_SNAPSHOT_CASES
    }


def _build_mode_annotation_snapshot(case):
    fasta_path = FASTA_DIR / case["fasta"]
    snapshot = _build_example_fasta_annotation(
        fasta_path,
        linear=case["linear"],
        is_detailed=case["is_detailed"],
    )
    snapshot["fasta"] = case["fasta"]
    return snapshot


def _load_snapshot():
    if not SNAPSHOT_PATH.exists():
        raise AssertionError(
            f"Missing {SNAPSHOT_PATH}. Run with {UPDATE_SNAPSHOT_ENV}=1 to create it."
        )

    return json.loads(SNAPSHOT_PATH.read_text())


def _write_snapshot(snapshot):
    SNAPSHOT_PATH.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")


def _update_snapshot_if_requested():
    global _UPDATED_SNAPSHOT

    if os.environ.get(UPDATE_SNAPSHOT_ENV) and not _UPDATED_SNAPSHOT:
        _write_snapshot(_build_example_fasta_annotation_snapshot())
        _UPDATED_SNAPSHOT = True


def test_example_fasta_snapshot_lists_current_fastas():
    _update_snapshot_if_requested()

    current_fasta_names = {path.name for path in _get_example_fasta_paths()}
    expected = _load_snapshot()

    assert expected["schema_version"] == SNAPSHOT_SCHEMA_VERSION
    assert set(expected["examples"]) == current_fasta_names
    assert set(expected["mode_examples"]) == {
        case["id"] for case in MODE_SNAPSHOT_CASES
    }


@pytest.mark.skipif(
    bool(os.environ.get(UPDATE_SNAPSHOT_ENV)),
    reason="snapshot refresh is handled once by test_example_fasta_snapshot_lists_current_fastas",
)
@pytest.mark.parametrize(
    "fasta_path",
    _get_example_fasta_paths(),
    ids=lambda path: path.name,
)
def test_example_fasta_annotation_matches_snapshot(fasta_path):
    expected = _load_snapshot()
    expected_example = expected["examples"].get(fasta_path.name)

    assert expected_example is not None, (
        f"Missing snapshot entry for {fasta_path.name}. "
        f"Run with {UPDATE_SNAPSHOT_ENV}=1 to refresh it."
    )

    current_example = _build_example_fasta_annotation(fasta_path)

    assert current_example["sequence_length"] == expected_example["sequence_length"]
    assert current_example["annotation_count"] == expected_example["annotation_count"]
    assert current_example["annotations"] == expected_example["annotations"]


@pytest.mark.skipif(
    bool(os.environ.get(UPDATE_SNAPSHOT_ENV)),
    reason="snapshot refresh is handled once by test_example_fasta_snapshot_lists_current_fastas",
)
@pytest.mark.parametrize(
    "case",
    MODE_SNAPSHOT_CASES,
    ids=lambda case: case["id"],
)
def test_example_fasta_annotation_mode_matches_snapshot(case):
    expected = _load_snapshot()
    expected_example = expected["mode_examples"].get(case["id"])

    assert expected_example is not None, (
        f"Missing snapshot entry for {case['id']}. "
        f"Run with {UPDATE_SNAPSHOT_ENV}=1 to refresh it."
    )

    current_example = _build_mode_annotation_snapshot(case)

    assert current_example["fasta"] == expected_example["fasta"]
    assert current_example["linear"] == expected_example["linear"]
    assert current_example["is_detailed"] == expected_example["is_detailed"]
    assert current_example["sequence_length"] == expected_example["sequence_length"]
    assert current_example["annotation_count"] == expected_example["annotation_count"]
    assert current_example["annotations"] == expected_example["annotations"]
