"""Tests for annotation-control report generation."""

import json

from tests.annotation_control_utils import (
    CaseResult,
    render_markdown_report,
    write_report_artifact,
)


def test_annotation_report_summarizes_results_and_reasons(tmp_path):
    results = [
        CaseResult("matching", "regular", "passed", 2, 2, "CSV and GenBank match"),
        CaseResult(
            "changed",
            "detailed",
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
