from subprocess import CompletedProcess

import pytest

from plannotate import search


def test_external_search_failure_includes_tool_diagnostic(monkeypatch):
    monkeypatch.setattr(
        search.subprocess,
        "run",
        lambda *args, **kwargs: CompletedProcess(args, 1, stderr="database missing"),
    )

    with pytest.raises(RuntimeError, match="cmscan.*database missing"):
        search._run_external_command("cmscan query.fa", "cmscan")
