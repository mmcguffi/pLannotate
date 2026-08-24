"""Tests for annotation-control report generation."""

import json
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from tests.annotation_control_utils import (
    CaseResult,
    compare_genbank,
    render_markdown_report,
    write_report_artifact,
)


def _record(features):
    record = SeqRecord(Seq("ATGC" * 25), id="test", name="test")
    record.annotations["topology"] = "circular"
    record.annotations["molecule_type"] = "DNA"
    record.features = features
    return record


def _feature(label="AmpR", **qualifiers):
    # Biopython represents every qualifier as a list of values, which is what the
    # comparison reads out of a parsed control file
    values = {key: [value] for key, value in qualifiers.items()}
    return SeqFeature(
        FeatureLocation(0, 30, 1), type="CDS", qualifiers={"label": [label], **values}
    )


def test_annotation_report_summarizes_results_and_reasons(tmp_path):
    results = [
        CaseResult("matching", "default", "passed", 2, 2, "CSV and GenBank match"),
        CaseResult(
            "changed",
            "linear",
            "changed",
            3,
            2,
            "CSV changed: removed duplicate | features changed: removed duplicate",
        ),
    ]

    report = render_markdown_report(results, ["tool version differs"])

    assert "- Passed: 1" in report
    assert "- Changed: 1" in report
    assert "annotation count 3 → 2" in report
    assert "CSV rows or metadata differ" in report
    assert "GenBank features or qualifiers differ" in report
    assert "CSV changed: removed duplicate" in report

    write_report_artifact(tmp_path, results, ["tool version differs"])
    payload = json.loads((tmp_path / "report.json").read_text())
    assert payload["summary"] == {
        "total": 2,
        "passed": 1,
        "changed": 1,
        "errors": 0,
    }
    assert payload["context_warnings"] == ["tool version differs"]
    assert (tmp_path / "report.md").read_text() == report


def test_a_qualifier_only_change_is_reported_as_one():
    # An annotation that gained a qualifier sits at the same place under the same
    # label, so reporting it as an addition and a removal prints the same string
    # twice and hides what actually moved -- exactly what a run that adds a qualifier
    # to every feature produces.
    expected = _record([_feature()])
    actual = _record([_feature(copy_number="~5", identity="99.0")])

    summary = compare_genbank(actual, expected) or ""

    assert "added=" not in summary
    assert "removed=" not in summary
    assert "+copy_number=~5" in summary
    assert "+identity=99.0" in summary
    assert "AmpR (CDS, [0:30](+))" in summary


def test_a_changed_qualifier_value_reports_both_sides():
    expected = _record([_feature(identity="98.0")])
    actual = _record([_feature(identity="99.0")])

    summary = compare_genbank(actual, expected) or ""

    assert "identity: 98.0 -> 99.0" in summary


def test_a_genuinely_new_annotation_is_still_reported_as_added():
    expected = _record([_feature()])
    actual = _record([_feature(), _feature(label="KanR")])

    summary = compare_genbank(actual, expected) or ""

    assert "added=['KanR (CDS, [0:30](+))']" in summary
    assert "qualifiers=" not in summary


def test_matching_records_report_no_change():
    assert compare_genbank(_record([_feature()]), _record([_feature()])) is None


def test_written_genbank_round_trips_through_the_comparison():
    # the comparison reads controls from disk, so it must agree with itself across a
    # write and a read rather than only over in-memory records
    record = _record([_feature(copy_number="~5")])
    handle = StringIO()
    SeqIO.write(record, handle, "genbank")

    reparsed = SeqIO.read(StringIO(handle.getvalue()), "genbank")

    assert compare_genbank(reparsed, record) is None
