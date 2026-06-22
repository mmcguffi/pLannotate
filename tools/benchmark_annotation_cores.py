#!/usr/bin/env python3
"""Benchmark end-to-end annotation scaling and render a summary figure."""

import argparse
import csv
import logging
import random
import statistics
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import SeqIO  # noqa: E402

from plannotate.models import Construct  # noqa: E402
from plannotate.search import _search_all_databases_cached  # noqa: E402


def benchmark(fasta: Path, max_cores: int, repeats: int) -> list[dict]:
    """Return randomized, cache-independent timing observations."""
    sequence = str(SeqIO.read(fasta, "fasta").seq)
    logging.getLogger("plannotate").setLevel(logging.WARNING)

    # Warm database files and executables before collecting timings.
    _search_all_databases_cached.cache_clear()
    expected_annotations = len(Construct(sequence, cores=max_cores).features)

    schedule = [core for core in range(1, max_cores + 1) for _ in range(repeats)]
    random.Random(20260621).shuffle(schedule)
    observations = []
    repeat_by_core = {core: 0 for core in range(1, max_cores + 1)}
    for core in schedule:
        repeat_by_core[core] += 1
        _search_all_databases_cached.cache_clear()
        started = time.perf_counter()
        annotations = len(Construct(sequence, cores=core).features)
        elapsed = time.perf_counter() - started
        if annotations != expected_annotations:
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
        print(f"cores={core:2d} repeat={repeat_by_core[core]} seconds={elapsed:.3f}")
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


def render_figure(
    path: Path,
    observations: list[dict],
    summary: list[dict],
    fasta: Path,
    repeats: int,
) -> None:
    """Render individual runtimes and their average by core count."""
    cores = [row["cores"] for row in summary]
    averages = [row["average_seconds"] for row in summary]

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, runtime_axis = plt.subplots(figsize=(7.5, 5.2))
    color = "#2864A5"
    runtime_axis.scatter(
        [row["cores"] for row in observations],
        [row["seconds"] for row in observations],
        color=color,
        alpha=0.35,
        s=34,
        edgecolors="none",
        label="Individual runs",
        zorder=2,
    )
    runtime_axis.plot(
        cores,
        averages,
        color=color,
        marker="o",
        linewidth=3.5,
        markersize=7,
        label="Mean runtime",
        zorder=3,
    )
    runtime_axis.set(
        title="End-to-end annotation runtime",
        xlabel="Allocated cores",
        ylabel="Seconds",
    )
    runtime_axis.set_xticks(cores)
    runtime_axis.legend(frameon=False)

    best = max(summary, key=lambda row: row["speedup"])
    figure.suptitle(
        f"pLannotate core scaling — {fasta.name}\n"
        f"{repeats} runs per core count; {best['speedup']:.2f}× peak mean speedup",
        fontsize=13,
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("plannotate/data/fastas/pUC19.fa"),
    )
    parser.add_argument("--max-cores", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/benchmarks"),
    )
    args = parser.parse_args()
    if args.max_cores < 1 or args.repeats < 1:
        parser.error("--max-cores and --repeats must be positive")

    observations = benchmark(args.fasta, args.max_cores, args.repeats)
    summary = summarize(observations, args.max_cores)
    stem = f"{args.fasta.stem}-core-scaling"
    write_csv(args.output_dir / f"{stem}-raw.csv", observations)
    write_csv(args.output_dir / f"{stem}.csv", summary)
    render_figure(
        args.output_dir / f"{stem}.png",
        observations,
        summary,
        args.fasta,
        args.repeats,
    )


if __name__ == "__main__":
    main()
