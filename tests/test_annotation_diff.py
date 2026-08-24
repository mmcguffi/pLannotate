"""Tests for the branch-to-branch annotation diff helper."""

import json
from argparse import Namespace
from pathlib import Path

from Bio import SeqIO

from tools.annotation_diff import cases_for, comparison_policy_note, generate


def test_cases_cover_both_supported_topologies_for_every_fasta():
    paths = [Path("alpha.fa"), Path("beta.fa")]

    assert [case.id for case in cases_for(paths)] == [
        "default:alpha",
        "linear:alpha",
        "default:beta",
        "linear:beta",
    ]


def test_generate_automatically_uses_detailed_on_legacy_checkout(monkeypatch, tmp_path):
    fastas = tmp_path / "fastas"
    fastas.mkdir()
    (fastas / "record.fa").write_text(">record\nACGT\n")
    output = tmp_path / "output"
    calls = []

    class FakeTable:
        def to_csv(self, path, index=False):
            assert index is False
            Path(path).write_text("annotation\n")

    class LegacyConstruct:
        def __init__(self, seq, linear=False, detailed=False):
            calls.append((str(seq), linear, detailed))

        def to_csv(self):
            return FakeTable()

        def to_genbank(self):
            return "LOCUS legacy\n"

    monkeypatch.setattr("plannotate.models.Construct", LegacyConstruct)

    assert generate(Namespace(fastas=fastas, out=output)) == 0
    assert calls == [("ACGT", False, True), ("ACGT", True, True)]
    metadata = json.loads((output / "run-metadata.json").read_text())
    assert metadata == {"legacy_detailed_applied": True}
    assert (output / "default/record.csv").is_file()
    assert (output / "linear/record.csv").is_file()


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
        )
    )

    assert result == 0
    assert json.loads((output / "errors.json").read_text())
    metadata = json.loads((output / "run-metadata.json").read_text())
    assert metadata["legacy_detailed_applied"] is False
