"""Golden-output regression tests for annotation results."""

import warnings

import pytest

from tests.annotation_control_utils import (
    CONTROL_CASES,
    CONTROL_DIR,
    context_changes,
    evaluate_case,
)

pytestmark = pytest.mark.integration


class AnnotationControlChangedWarning(UserWarning):
    """An annotation differs from its non-blocking published control."""


def _report_change(message, config):
    print(f"\nANNOTATION CONTROL DIFFERENCE\n{message}\n", flush=True)
    warnings.warn(message, AnnotationControlChangedWarning, stacklevel=2)
    if config.getoption("--strict-annotation-controls"):
        pytest.fail(message, pytrace=False)
    pytest.xfail(message)


def test_annotation_control_files_cover_default_behavior():
    expected_by_control = {
        control_mode: {
            case.fasta_path.stem
            for case in CONTROL_CASES
            if case.control_mode == control_mode
        }
        for control_mode in {case.control_mode for case in CONTROL_CASES}
    }
    assert set(expected_by_control) == {"detailed", "detailed-linear"}
    for control_mode, expected_stems in expected_by_control.items():
        mode_dir = CONTROL_DIR / control_mode
        assert {path.stem for path in mode_dir.glob("*.csv")} == expected_stems
        assert {path.stem for path in mode_dir.glob("*.gbk")} == expected_stems


def test_annotation_control_context(request):
    if changes := context_changes():
        _report_change(
            "Annotation control context changed: " + "; ".join(changes),
            request.config,
        )


@pytest.mark.parametrize("case", CONTROL_CASES, ids=lambda case: case.id)
def test_annotation_output_matches_control(case, request):
    if changes := context_changes():
        warnings.warn(
            "Annotation controls cannot be compared reliably: " + "; ".join(changes),
            AnnotationControlChangedWarning,
            stacklevel=2,
        )
        pytest.skip("annotation control context differs")

    result, _, _ = evaluate_case(case)
    if result.status == "changed":
        _report_change(
            f"Annotation control changed for {case.id}: {result.reason}",
            request.config,
        )
