"""Tests for the branch-to-branch annotation diff helper."""

import json
from argparse import Namespace

from Bio import SeqIO

from tools.annotation_diff import comparison_policy_note, generate


def test_comparison_policy_note_discloses_legacy_base_normalization():
    note = comparison_policy_note(
        {"legacy_detailed_applied": True},
        {"legacy_detailed_applied": False},
    )

    assert "base was explicitly run" in note
    assert "does not show the user-visible default change" in note


def test_comparison_policy_note_is_empty_without_normalization():
    assert comparison_policy_note({}, {}) == ""


def test_generate_records_metadata_when_every_case_errors(monkeypatch, tmp_path):
    fastas = tmp_path / "fastas"
    fastas.mkdir()
    (fastas / "invalid.fa").write_text(">record\nACGT\n")
    output = tmp_path / "output"
    monkeypatch.setattr(
        SeqIO,
        "read",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("invalid FASTA")),
    )

    result = generate(
        Namespace(
            fastas=fastas,
            out=output,
            legacy_detailed_if_supported=True,
        )
    )

    assert result == 0
    assert json.loads((output / "errors.json").read_text())
    metadata = json.loads((output / "run-metadata.json").read_text())
    assert metadata["legacy_detailed_requested"] is True
    assert metadata["legacy_detailed_applied"] is False
