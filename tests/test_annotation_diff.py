"""Tests for the branch-to-branch annotation diff helper."""

from tools.annotation_diff import comparison_policy_note


def test_comparison_policy_note_discloses_legacy_base_normalization():
    note = comparison_policy_note(
        {"legacy_detailed_applied": True},
        {"legacy_detailed_applied": False},
    )

    assert "base was explicitly run" in note
    assert "does not show the user-visible default change" in note


def test_comparison_policy_note_is_empty_without_normalization():
    assert comparison_policy_note({}, {}) == ""
