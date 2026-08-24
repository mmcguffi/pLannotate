"""Fast structural checks for the checked-in annotation controls."""

from tests.annotation_control_utils import CONTROL_CASES, CONTROL_DIR


def test_annotation_control_files_cover_default_behavior():
    assert {path.name for path in CONTROL_DIR.iterdir() if path.is_dir()} == {
        "default",
        "linear",
    }
    expected_by_mode = {
        mode: {case.fasta_path.stem for case in CONTROL_CASES if case.mode == mode}
        for mode in {case.mode for case in CONTROL_CASES}
    }
    assert set(expected_by_mode) == {"default", "linear"}
    for mode, expected_stems in expected_by_mode.items():
        mode_dir = CONTROL_DIR / mode
        assert {path.stem for path in mode_dir.glob("*.csv")} == expected_stems
        assert {path.stem for path in mode_dir.glob("*.gbk")} == expected_stems
