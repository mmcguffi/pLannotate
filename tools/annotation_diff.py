#!/usr/bin/env python3
"""Generate and diff annotation outputs between two pLannotate checkouts.

This backs the ``@plannotate-bot`` PR command (see
``.github/workflows/annotation-diff.yml``): it annotates a fixed set of FASTAs
with whatever ``plannotate`` is installed, then compares the ``main`` and PR
outputs into a Markdown table for manual inspection.

The tool has two subcommands with deliberately different dependency footprints:

- ``generate`` runs once under *each* checkout's environment. It imports only
  ``plannotate`` (plus Biopython/pandas), so it works on ``main`` and on the PR
  branch even though the branches expose different CLIs and helper modules. It
  is copied to a stable location before the workflow switches branches.
- ``report`` runs once, under the PR branch's environment. It reuses the
  branch's tested ``compare_csv``/``compare_genbank`` helpers, imported lazily
  so ``generate`` never needs them.
"""

import argparse
import inspect
import json
import sys
from dataclasses import dataclass
from pathlib import Path


# The case matrix mirrors tests/annotation_control_utils.CONTROL_CASES, but is
# duplicated here (rather than imported) so `generate` stays importable under a
# `main` checkout that predates the annotation-controls test module.
@dataclass(frozen=True)
class Case:
    mode: str
    stem: str
    linear: bool = False

    @property
    def id(self) -> str:
        return f"{self.mode}:{self.stem}"


PER_FASTA_MODES = (
    ("default", False),
    ("linear", True),
)


def cases_for(fasta_paths):
    cases = [
        Case(mode, path.stem, linear=linear)
        for path in fasta_paths
        for mode, linear in PER_FASTA_MODES
    ]
    return cases


def _fasta_paths(fasta_dir: Path):
    paths = sorted(fasta_dir.glob("*.fa"))
    if not paths:
        raise SystemExit(f"no *.fa inputs found in {fasta_dir}")
    return paths


def generate(args) -> int:
    """Annotate every case with the installed plannotate; write csv + gbk."""
    from Bio import SeqIO  # noqa: PLC0415 — deferred so `report` needn't import Bio

    from plannotate.models import Construct  # noqa: PLC0415

    fasta_dir = args.fastas.resolve()
    out_dir = args.out.resolve()
    fasta_paths = _fasta_paths(fasta_dir)
    by_stem = {path.stem: path for path in fasta_paths}
    cases = cases_for(fasta_paths)

    errors: dict[str, str] = {}
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] {case.id}", flush=True)
        case_dir = out_dir / case.mode
        case_dir.mkdir(parents=True, exist_ok=True)
        try:
            sequence = SeqIO.read(by_stem[case.stem], "fasta").seq
            kwargs = {"seq": sequence, "linear": case.linear}
            # This script is copied before CI checks out the base revision. Opt in
            # there so both old and new revisions exercise the same behavior.
            if "detailed" in inspect.signature(Construct).parameters:
                kwargs["detailed"] = True
            construct = Construct(**kwargs)
            construct.to_csv().to_csv(case_dir / f"{case.stem}.csv", index=False)
            (case_dir / f"{case.stem}.gbk").write_text(construct.to_genbank())
        except Exception as exc:  # noqa: BLE001 — record and continue past bad cases
            errors[case.id] = f"{type(exc).__name__}: {exc}"
            print(f"    ERROR: {errors[case.id]}", file=sys.stderr, flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "errors.json").write_text(json.dumps(errors, indent=2) + "\n")
    print(f"Generated {len(cases) - len(errors)}/{len(cases)} cases into {out_dir}")
    return 0


def _load_errors(directory: Path) -> dict[str, str]:
    path = directory / "errors.json"
    if path.is_file():
        return json.loads(path.read_text())
    return {}


def report(args) -> int:
    """Compare base (main) and head (PR) output dirs into a Markdown table."""
    import pandas as pd  # noqa: PLC0415
    from Bio import SeqIO  # noqa: PLC0415

    # Reuse the branch's tested diff/format helpers; only available on the PR
    # checkout, which is where `report` always runs.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from tests.annotation_control_utils import (  # noqa: PLC0415
        CaseResult,
        compare_csv,
        compare_genbank,
        read_genbank,
        render_markdown_report,
    )

    base_dir = args.base.resolve()
    head_dir = args.head.resolve()
    base_errors = _load_errors(base_dir)
    head_errors = _load_errors(head_dir)

    # Discover cases from the union of what either side produced.
    discovered: dict[str, tuple[str, str]] = {}
    for root in (base_dir, head_dir):
        for csv_path in root.glob("*/*.csv"):
            discovered[f"{csv_path.parent.name}:{csv_path.stem}"] = (
                csv_path.parent.name,
                csv_path.stem,
            )
    all_ids = sorted(set(discovered) | set(base_errors) | set(head_errors))

    results: list[CaseResult] = []
    for case_id in all_ids:
        mode, stem = discovered.get(case_id, tuple(case_id.split(":", 1)))
        if case_id in base_errors or case_id in head_errors:
            reason = "; ".join(
                filter(
                    None,
                    (
                        f"main errored: {base_errors[case_id]}"
                        if case_id in base_errors
                        else "",
                        f"this branch errored: {head_errors[case_id]}"
                        if case_id in head_errors
                        else "",
                    ),
                )
            )
            results.append(CaseResult(stem, mode, "error", None, None, reason))
            continue

        base_csv_path = base_dir / mode / f"{stem}.csv"
        head_csv_path = head_dir / mode / f"{stem}.csv"
        if not base_csv_path.is_file() or not head_csv_path.is_file():
            missing = "main" if not base_csv_path.is_file() else "this branch"
            results.append(
                CaseResult(
                    stem, mode, "error", None, None, f"missing output on {missing}"
                )
            )
            continue

        try:
            base_csv = pd.read_csv(base_csv_path)
            head_csv = pd.read_csv(head_csv_path)
            base_gbk = SeqIO.read(base_dir / mode / f"{stem}.gbk", "genbank")
            head_gbk = read_genbank((head_dir / mode / f"{stem}.gbk").read_text())
            # compare_*(actual, expected): actual = this branch, expected = main.
            changes = [
                change
                for change in (
                    compare_csv(head_csv, base_csv),
                    compare_genbank(head_gbk, base_gbk),
                )
                if change
            ]
        except Exception as exc:  # noqa: BLE001 — a bad case must not sink the report
            results.append(
                CaseResult(
                    stem, mode, "error", None, None, f"{type(exc).__name__}: {exc}"
                )
            )
            continue

        results.append(
            CaseResult(
                plasmid=stem,
                mode=mode,
                status="changed" if changes else "passed",
                control_annotations=len(base_csv),
                current_annotations=len(head_csv),
                reason=" | ".join(changes) if changes else "identical",
            )
        )

    markdown = render_markdown_report(results)
    # The reused renderer is written for control comparisons; relabel its
    # headings and columns so the table reads as a main-vs-branch diff.
    markdown = (
        markdown.replace(
            "# Annotation control comparison",
            "# Annotation diff: `main` vs this branch",
        )
        .replace(
            "- Annotation changes are informational unless strict mode is requested.\n",
            "",
        )
        .replace("| Result | Control | Current |", "| Result | main | branch |")
    )
    markdown = (
        "<!-- annotation-diff-bot -->\n"
        "> Full annotation pipeline run on `main` and on this PR branch over the "
        "packaged FASTAs. **Control = `main`, Current = this branch.** Rows marked "
        "CHANGED differ — inspect them manually.\n\n" + markdown
    )
    args.out.write_text(markdown)
    changed = sum(result.status == "changed" for result in results)
    errored = sum(result.status == "error" for result in results)
    print(
        f"Wrote {args.out}: {len(results) - changed - errored} identical, "
        f"{changed} changed, {errored} errors"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate_parser = subparsers.add_parser(
        "generate", help="annotate the fixed FASTA set with the installed plannotate"
    )
    generate_parser.add_argument("--fastas", type=Path, required=True)
    generate_parser.add_argument("--out", type=Path, required=True)
    generate_parser.set_defaults(function=generate)

    report_parser = subparsers.add_parser(
        "report", help="diff two generated output dirs into a Markdown table"
    )
    report_parser.add_argument("--base", type=Path, required=True, help="main outputs")
    report_parser.add_argument("--head", type=Path, required=True, help="PR outputs")
    report_parser.add_argument("--out", type=Path, required=True, help="report.md path")
    report_parser.set_defaults(function=report)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    raise SystemExit(args.function(args))


if __name__ == "__main__":
    main()
