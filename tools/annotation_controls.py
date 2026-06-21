#!/usr/bin/env python3
"""Regenerate and compare the published annotation control outputs."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests.annotation_control_utils import (  # noqa: E402
    CONTROL_CASES,
    CONTROL_DIR,
    CaseResult,
    context_changes,
    error_result,
    evaluate_case,
    write_report_artifact,
)


BASELINE_ENVIRONMENT = CONTROL_DIR / "baseline-environment.yml"
DEFAULT_PREFIX = ROOT / ".annotation-controls" / "plannotate-1.2.5"
DEFAULT_ARTIFACT_DIR = ROOT / "artifacts" / "annotation-controls"
EXPECTED_BASELINE_PACKAGES = {
    "plannotate": "1.2.5",
    "python": "3.12.13",
    "blast": "2.16.0",
    "diamond": "2.2.2",
    "infernal": "1.1.5",
    "biopython": "1.87",
    "numpy": "1.26.4",
    "pandas": "2.3.3",
    "bokeh": "2.4.3",
    "streamlit": "1.41.1",
}


def _run(command, **kwargs):
    print("+", " ".join(str(part) for part in command), flush=True)
    try:
        return subprocess.run(command, check=True, **kwargs)
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            print(exc.stdout, file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        raise


def _conda_packages(conda, prefix):
    result = _run(
        [conda, "list", "--prefix", str(prefix), "--json"],
        capture_output=True,
        text=True,
    )
    return {package["name"]: package for package in json.loads(result.stdout)}


def _ensure_baseline_environment(args):
    prefix = args.prefix.resolve()
    if args.recreate and prefix.exists():
        shutil.rmtree(prefix)
    if not (prefix / "conda-meta" / "history").exists():
        prefix.parent.mkdir(parents=True, exist_ok=True)
        _run(
            [
                args.conda,
                "env",
                "create",
                "--prefix",
                str(prefix),
                "--file",
                str(BASELINE_ENVIRONMENT),
            ]
        )

    packages = _conda_packages(args.conda, prefix)
    mismatches = {
        name: (version, packages.get(name, {}).get("version"))
        for name, version in EXPECTED_BASELINE_PACKAGES.items()
        if packages.get(name, {}).get("version") != version
    }
    package = packages.get("plannotate", {})
    if package.get("build_string") != "pyhdfd78af_0":
        mismatches["plannotate build"] = (
            "pyhdfd78af_0",
            package.get("build_string"),
        )
    if mismatches:
        raise RuntimeError(
            f"baseline environment differs from its pins: {mismatches}; use --recreate"
        )
    return prefix


def regenerate(args):
    prefix = _ensure_baseline_environment(args)
    executable = prefix / ("Scripts/plannotate.exe" if os.name == "nt" else "bin/plannotate")
    path_dir = executable.parent
    environment = os.environ.copy()
    environment["PATH"] = os.pathsep.join([str(path_dir), environment.get("PATH", "")])

    with tempfile.TemporaryDirectory(prefix="plannotate-controls-") as temp_name:
        staging = Path(temp_name)
        for case in CONTROL_CASES:
            output_dir = staging / case.mode
            output_dir.mkdir(parents=True, exist_ok=True)
            command = [
                str(executable),
                "batch",
                "--input",
                str(case.fasta_path),
                "--output",
                str(output_dir),
                "--file_name",
                case.fasta_path.stem,
                "--suffix",
                "",
                "--csv",
            ]
            if case.linear:
                command.append("--linear")
            if case.detailed:
                command.append("--detailed")
            _run(command, env=environment, capture_output=True, text=True)

        for case in CONTROL_CASES:
            destination = args.output_dir / case.mode
            destination.mkdir(parents=True, exist_ok=True)
            for extension in ("csv", "gbk"):
                source = staging / case.mode / f"{case.fasta_path.stem}.{extension}"
                if not source.is_file():
                    raise RuntimeError(f"baseline did not create {source}")
                shutil.copy2(source, destination / source.name)

    print(f"Regenerated {len(CONTROL_CASES)} control pairs in {args.output_dir}")
    return 0


def compare(args):
    output_dir = args.output_dir.resolve()
    current_dir = output_dir / "current"
    baseline_dir = output_dir / "baseline"
    for generated_dir in (current_dir, baseline_dir):
        if generated_dir.exists():
            shutil.rmtree(generated_dir)
    results: list[CaseResult] = []
    warnings = context_changes()
    for index, case in enumerate(CONTROL_CASES, start=1):
        print(f"[{index}/{len(CONTROL_CASES)}] {case.id}", flush=True)
        baseline_case_dir = baseline_dir / case.mode
        baseline_case_dir.mkdir(parents=True, exist_ok=True)
        for extension in ("csv", "gbk"):
            source = Path(f"{case.control_stem}.{extension}")
            shutil.copy2(source, baseline_case_dir / source.name)
        try:
            result, actual_csv, actual_gbk = evaluate_case(case)
        except Exception as exc:
            result = error_result(case, exc)
        else:
            case_dir = current_dir / case.mode
            case_dir.mkdir(parents=True, exist_ok=True)
            actual_csv.to_csv(case_dir / f"{case.fasta_path.stem}.csv", index=False)
            (case_dir / f"{case.fasta_path.stem}.gbk").write_text(actual_gbk)
        results.append(result)

    write_report_artifact(output_dir, results, warnings)
    changed = sum(result.status == "changed" for result in results)
    errors = sum(result.status == "error" for result in results)
    print(
        f"Wrote {output_dir / 'report.md'}: "
        f"{len(results) - changed - errors} passed, {changed} changed, {errors} errors"
    )
    return 1 if args.strict and (changed or errors or warnings) else 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    regenerate_parser = subparsers.add_parser(
        "regenerate", help="recreate controls with vanilla pLannotate 1.2.5"
    )
    regenerate_parser.add_argument(
        "--conda",
        default=shutil.which("mamba") or shutil.which("conda") or "conda",
        help="conda-compatible executable (defaults to mamba when available)",
    )
    regenerate_parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
    regenerate_parser.add_argument("--output-dir", type=Path, default=CONTROL_DIR)
    regenerate_parser.add_argument(
        "--recreate", action="store_true", help="recreate the pinned conda environment"
    )
    regenerate_parser.set_defaults(function=regenerate)

    compare_parser = subparsers.add_parser(
        "compare", help="compare the current checkout and write an artifact"
    )
    compare_parser.add_argument("--output-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
    compare_parser.add_argument(
        "--strict", action="store_true", help="exit nonzero for changes or errors"
    )
    compare_parser.set_defaults(function=compare)
    return parser


def main():
    args = build_parser().parse_args()
    raise SystemExit(args.function(args))


if __name__ == "__main__":
    main()
