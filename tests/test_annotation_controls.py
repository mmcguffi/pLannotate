import warnings

import pytest

from tests.annotation_control_utils import (
    CONTROL_CASES,
    CONTROL_DIR,
    FASTA_PATHS,
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
