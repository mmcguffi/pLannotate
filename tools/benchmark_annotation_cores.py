#!/usr/bin/env python3
"""Benchmark end-to-end annotation scaling and render a summary figure.

Each timing observation runs in a fresh Python subprocess. That keeps repeated
runs independent of pLannotate's in-process annotation cache without reaching
into private implementation details.
"""

import argparse
import csv
import json
import logging
import os
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO  # noqa: E402

from plannotate.models import Construct  # noqa: E402


def _run_annotation(fasta: Path, cores: int) -> dict[str, float | int]:
    """Annotate one construct and return a timing observation."""
    sequence = str(SeqIO.read(fasta, "fasta").seq)
    logging.getLogger("plannotate").setLevel(logging.WARNING)
    started = time.perf_counter()
    annotations = len(Construct(sequence, cores=cores).features)
    elapsed = time.perf_counter() - started
    return {"seconds": elapsed, "annotations": annotations}


def _run_annotation_subprocess(
    fasta: Path,
    cores: int,
    timeout_seconds: int,
) -> dict[str, float | int]:
    """Run one annotation timing in a fresh Python process."""
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--fasta",
        str(fasta),
        "--cores",
        str(cores),
    ]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
        text=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.communicate()
        raise TimeoutError(
            f"Annotation subprocess timed out after {timeout_seconds}s at cores={cores}"
        ) from exc

    if process.returncode != 0:
        raise RuntimeError(
            f"Annotation subprocess failed at cores={cores} "
            f"with exit code {process.returncode}:\n{stderr}"
        )
    return json.loads(stdout)


def benchmark(
    fasta: Path,
    max_cores: int,
    repeats: int,
    seed: int,
    timeout_seconds: int,
    retries: int,
) -> list[dict]:
    """Return randomized, cache-independent timing observations."""
    schedule = [core for core in range(1, max_cores + 1) for _ in range(repeats)]
    random.Random(seed).shuffle(schedule)
    observations = []
    repeat_by_core = {core: 0 for core in range(1, max_cores + 1)}
    expected_annotations = None
    for core in schedule:
        repeat_by_core[core] += 1
        result = None
        for attempt in range(retries + 1):
            try:
                result = _run_annotation_subprocess(
                    fasta.resolve(),
                    core,
                    timeout_seconds,
                )
                break
            except TimeoutError:
                if attempt == retries:
                    raise
                print(
                    f"cores={core:2d} repeat={repeat_by_core[core]} "
                    f"timed out; retrying",
                    flush=True,
                )
        if result is None:
            raise RuntimeError(f"No timing result collected for cores={core}")
        annotations = int(result["annotations"])
        elapsed = float(result["seconds"])
        if expected_annotations is None:
            expected_annotations = annotations
        elif annotations != expected_annotations:
            raise RuntimeError(
                f"Annotation count changed at {core} cores: "
                f"expected {expected_annotations}, found {annotations}"
            )
        observations.append(
            {
                "cores": core,
                "repeat": repeat_by_core[core],
                "seconds": elapsed,
                "annotations": annotations,
            }
        )
        print(
            f"cores={core:2d} repeat={repeat_by_core[core]} seconds={elapsed:.3f}",
            flush=True,
        )
    return observations


def summarize(observations: list[dict], max_cores: int) -> list[dict]:
    """Summarize timing observations by core count."""
    elapsed_by_core = {
        core: [row["seconds"] for row in observations if row["cores"] == core]
        for core in range(1, max_cores + 1)
    }
    baseline = statistics.mean(elapsed_by_core[1])
    summary = []
    for core, elapsed in elapsed_by_core.items():
        average = statistics.mean(elapsed)
        summary.append(
            {
                "cores": core,
                "average_seconds": average,
                "median_seconds": statistics.median(elapsed),
                "min_seconds": min(elapsed),
                "max_seconds": max(elapsed),
                "speedup": baseline / average,
                "efficiency": baseline / average / core,
            }
        )
    return summary


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def render_comparison_figure(
    path: Path,
    summaries: dict[str, list[dict]],
    repeats: int,
) -> None:
    """Render runtime and speedup curves for several plasmids on shared axes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, (runtime_axis, speedup_axis) = plt.subplots(1, 2, figsize=(13, 5.2))
    palette = ["#2864A5", "#C44E52", "#55A868", "#8172B3"]

    for index, (name, summary) in enumerate(summaries.items()):
        color = palette[index % len(palette)]
        cores = [row["cores"] for row in summary]
        runtime_axis.plot(
            cores,
            [row["average_seconds"] for row in summary],
            color=color,
            marker="o",
            linewidth=2.5,
            markersize=6,
            label=name,
        )
        speedup_axis.plot(
            cores,
            [row["speedup"] for row in summary],
            color=color,
            marker="o",
            linewidth=2.5,
            markersize=6,
            label=name,
        )

    any_summary = next(iter(summaries.values()))
    cores = [row["cores"] for row in any_summary]
    speedup_axis.plot(
        cores, cores, color="#888888", linestyle="--", linewidth=1.5, label="Ideal"
    )

    runtime_axis.set(
        title="End-to-end annotation runtime",
        xlabel="Allocated cores",
        ylabel="Seconds",
    )
    speedup_axis.set(
        title="Speedup over one core",
        xlabel="Allocated cores",
        ylabel="Speedup (x)",
    )
    for axis in (runtime_axis, speedup_axis):
        axis.set_xticks(cores)
        axis.legend(frameon=False)

    figure.suptitle(
        f"pLannotate core scaling - {', '.join(summaries)}\n"
        f"{repeats} runs per core count",
        fontsize=13,
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def copy_docs_image(source: Path, docs_image_dir: Path) -> Path:
    """Copy the generated benchmark plot to the documentation image directory."""
    docs_image_dir.mkdir(parents=True, exist_ok=True)
    destination = docs_image_dir / source.name
    shutil.copy2(source, destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--worker",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        nargs="+",
        default=[Path("plannotate/data/fastas/pUC19.fa")],
        help="One or more FASTA files; multiple files add a combined comparison plot.",
    )
    parser.add_argument("--cores", type=int, default=1, help=argparse.SUPPRESS)
    parser.add_argument("--max-cores", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260621)
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=120,
        help="Per-observation timeout for a single annotation subprocess.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=1,
        help="Retries for an annotation subprocess that exceeds the timeout.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/benchmarks"),
    )
    parser.add_argument(
        "--docs-image-dir",
        type=Path,
        default=Path("docs/images"),
        help="Directory that receives a copy of the rendered benchmark plot.",
    )
    parser.add_argument(
        "--skip-docs-copy",
        action="store_true",
        help="Only write the benchmark outputs under --output-dir.",
    )
    args = parser.parse_args()
    if args.worker:
        if args.cores < 1:
            parser.error("--cores must be positive")
        print(json.dumps(_run_annotation(args.fasta[0], args.cores)))
        return

    if args.max_cores < 1 or args.repeats < 1 or args.timeout_seconds < 1:
        parser.error("--max-cores, --repeats, and --timeout-seconds must be positive")
    if args.retries < 0:
        parser.error("--retries must be zero or positive")

    summaries: dict[str, list[dict]] = {}
    for fasta in args.fasta:
        print(f"Benchmarking {fasta.name}", flush=True)
        observations = benchmark(
            fasta,
            args.max_cores,
            args.repeats,
            args.seed,
            args.timeout_seconds,
            args.retries,
        )
        summary = summarize(observations, args.max_cores)
        summaries[fasta.stem] = summary
        stem = f"{fasta.stem}-core-scaling"
        write_csv(args.output_dir / f"{stem}-raw.csv", observations)
        write_csv(args.output_dir / f"{stem}.csv", summary)

    comparison_path = args.output_dir / "core-scaling-comparison.png"
    render_comparison_figure(comparison_path, summaries, args.repeats)
    if not args.skip_docs_copy:
        copy_docs_image(comparison_path, args.docs_image_dir)


if __name__ == "__main__":
    main()
