#!/usr/bin/env python3
"""Visualize the installed pLannotate annotation databases.

Reads the database bundle that ``plannotate setupdb`` installs and renders a
summary figure describing each annotation source: how many entries it holds,
how much disk it occupies, and -- for the curated SnapGene set -- the mix of
feature types it can annotate.

Run it after downloading the databases to inspect a bundle:

    python scripts/plot_databases.py --output database_overview.png

By default it inspects the bundle reported by ``plannotate`` itself; pass
``--data-dir`` to point at a different installation.
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless by default; safe for CI and servers
import matplotlib.pyplot as plt

from plannotate import _package_data
from plannotate._sqlite import get_descriptions_db_path

logger = logging.getLogger("plannotate.plot_databases")

# stable color per search method so the same source keeps its hue across plots
METHOD_COLORS = {
    "blastn": "#4C72B0",
    "diamond": "#DD8452",
    "infernal": "#55A868",
}
METHOD_LABELS = {
    "blastn": "BLAST (nucleotide)",
    "diamond": "DIAMOND (protein)",
    "infernal": "Infernal (ncRNA)",
}


@dataclass
class SourceStats:
    """Descriptive statistics for a single annotation source."""

    name: str
    method: str
    entries: int
    size_bytes: int
    version: str
    type_counts: dict[str, int] = field(default_factory=dict)


def _source_size_bytes(db_loc: Path) -> int:
    """Sum every on-disk file that backs a source (db plus search indexes)."""
    base = db_loc.parent
    stem = db_loc.stem if db_loc.suffix else db_loc.name
    # NOTE: the per-source descriptions DB is named ``{source}.db`` beside the search
    # index, so the same glob folds it into the footprint automatically.
    return sum(p.stat().st_size for p in base.glob(f"{stem}.*") if p.is_file())


def _sqlite_counts(descriptions_db: Path, table: str) -> tuple[int, dict[str, int]]:
    """Return (total entries, per-feature-type counts) from a descriptions DB."""
    with sqlite3.connect(f"file:{descriptions_db}?mode=ro", uri=True) as con:
        total = con.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        rows = con.execute(
            f"SELECT type, count(*) FROM {table} GROUP BY type ORDER BY count(*) DESC"
        ).fetchall()
    return total, {str(feature_type): n for feature_type, n in rows}


def _covariance_model_count(cm_path: Path) -> int:
    """Count covariance models in an Infernal ``.cm`` file.

    Each family is stored as two stacked records: the covariance model and an
    embedded HMM filter Infernal uses to pre-screen candidates. Both carry their
    own ``NAME``/``ACC``/``//`` lines, so counting any of those double-counts
    families. Only the covariance model opens with an ``INFERNAL`` header, so
    counting those headers yields the true family total.
    """
    count = 0
    with cm_path.open("rt", errors="ignore") as handle:
        for line in handle:
            if line.startswith("INFERNAL"):
                count += 1
    return count


def collect_stats(data_dir: Path) -> list[SourceStats]:
    """Gather descriptive statistics for every configured annotation source."""
    config = _package_data.get_yaml(_package_data.get_yaml_path())
    try:
        manifest = _package_data.get_database_manifest()
    except (FileNotFoundError, ValueError):
        manifest = {}
    versions = {
        name: entry.get("version", "unknown")
        for name, entry in manifest.get("databases", {}).items()
        if isinstance(entry, dict)
    }

    stats: list[SourceStats] = []
    for name, source in config.items():
        method = source["method"]
        db_loc = Path(source["db_loc"])
        size = _source_size_bytes(db_loc)

        entries = 0
        type_counts: dict[str, int] = {}
        if method == "infernal":
            if db_loc.is_file():
                entries = _covariance_model_count(db_loc)
        else:
            # resolve the descriptions DB exactly like the runtime does
            descriptions_db = get_descriptions_db_path(name, source)
            if descriptions_db.is_file():
                entries, type_counts = _sqlite_counts(descriptions_db, name)

        # manifest keys are lowercased for some sources; fall back gracefully
        version = versions.get(name) or versions.get(name.lower(), "unknown")
        stats.append(
            SourceStats(
                name=name,
                method=method,
                entries=entries,
                size_bytes=size,
                version=version,
                type_counts=type_counts,
            )
        )

    # largest databases first reads most naturally on horizontal bars
    stats.sort(key=lambda s: s.entries, reverse=True)
    return stats


def _bar_color(method: str) -> str:
    return METHOD_COLORS.get(method, "#937860")


def _human_bytes(num: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if num < 1024 or unit == "GB":
            return f"{num:.0f} {unit}" if unit != "B" else f"{num:.0f} B"
        num /= 1024
    return f"{num:.0f} GB"


def make_figure(stats: list[SourceStats], build_date: str | None) -> plt.Figure:
    """Render the database overview figure."""
    names = [s.name for s in stats]
    colors = [_bar_color(s.method) for s in stats]

    fig, (ax_entries, ax_size, ax_types) = plt.subplots(1, 3, figsize=(15, 5))

    # panel 1: number of annotatable entries per source (log scale; spans 1k -> 575k)
    entries = [max(s.entries, 1) for s in stats]
    ax_entries.barh(names, entries, color=colors)
    ax_entries.set_xscale("log")
    ax_entries.set_xlabel("entries (log scale)")
    ax_entries.set_title("Annotatable entries per database")
    ax_entries.invert_yaxis()
    for y, value in enumerate(entries):
        ax_entries.text(value, y, f" {value:,}", va="center", fontsize=9)

    # panel 2: on-disk footprint per source
    sizes_mb = [s.size_bytes / 1e6 for s in stats]
    ax_size.barh(names, sizes_mb, color=colors)
    ax_size.set_xscale("log")
    ax_size.set_xlabel("on-disk size, MB (log scale)")
    ax_size.set_title("Bundle size per database")
    ax_size.invert_yaxis()
    for y, s in enumerate(stats):
        ax_size.text(
            max(s.size_bytes / 1e6, 0.01),
            y,
            f" {_human_bytes(s.size_bytes)}",
            va="center",
            fontsize=9,
        )

    # panel 3: feature-type composition of the most diverse (curated) source
    diverse = max(stats, key=lambda s: len(s.type_counts), default=None)
    if diverse and diverse.type_counts:
        top = list(diverse.type_counts.items())[:12]
        labels = [t for t, _ in top]
        counts = [n for _, n in top]
        ax_types.barh(labels, counts, color=_bar_color(diverse.method))
        ax_types.set_xlabel("feature count")
        ax_types.set_title(f"{diverse.name} feature types")
        ax_types.invert_yaxis()
        for y, value in enumerate(counts):
            ax_types.text(value, y, f" {value:,}", va="center", fontsize=9)
    else:
        ax_types.axis("off")

    # shared legend mapping color -> search method
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=color)
        for method, color in METHOD_COLORS.items()
        if method in {s.method for s in stats}
    ]
    labels = [
        METHOD_LABELS.get(method, method)
        for method in METHOD_COLORS
        if method in {s.method for s in stats}
    ]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, 0.0),
    )

    title = "pLannotate annotation databases"
    if build_date:
        title += f"  (bundle built {build_date})"
    fig.suptitle(title, fontsize=14, fontweight="bold")
    # reserve top room for the title and bottom room for the shared legend
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="database directory to inspect (default: the installed bundle)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("database_overview.png"),
        help="image file to write (extension sets the format; default: %(default)s)",
    )
    parser.add_argument(
        "--dpi", type=int, default=150, help="output resolution (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="print per-database statistics"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )

    data_dir = args.data_dir or _package_data.get_data_directory()
    logger.info("Inspecting database bundle in %s", data_dir)

    stats = collect_stats(Path(data_dir))
    if not stats:
        raise SystemExit("No databases found; run 'plannotate setupdb' first.")

    for s in stats:
        logger.debug(
            "%-10s method=%-9s entries=%-8d size=%-9s version=%s",
            s.name,
            s.method,
            s.entries,
            _human_bytes(s.size_bytes),
            s.version,
        )

    try:
        manifest = _package_data.get_database_manifest()
        build_date = manifest.get("build_date")
    except (FileNotFoundError, ValueError):
        build_date = None

    fig = make_figure(stats, build_date)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    logger.info("Wrote %s", args.output.resolve())


if __name__ == "__main__":
    main()
