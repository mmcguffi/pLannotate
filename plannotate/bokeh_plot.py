"""Interactive Bokeh rendering for annotated DNA constructs."""

from collections.abc import Mapping
from math import pi
from typing import Any

import numpy as np
import pandas as pd

from . import _package_data

try:
    from bokeh.models import ColumnDataSource, HoverTool, Range1d, WheelZoomTool
    from bokeh.models.annotations import Label
    from bokeh.plotting import figure
except ImportError as exc:
    raise ImportError(
        "Bokeh is not installed. Install pLannotate with the 'plot' extra."
    ) from exc

BASE_RADIUS = 0.18
FEATURE_THICKNESS = 0.017
PLOT_SIZE = 0.35
PLOT_DIMENSIONS = 800
GLYPH_COLUMNS = [
    "x",
    "y",
    "Lx1",
    "Ly1",
    "annoLineColor",
    "lineX",
    "lineY",
    "theta",
    "text_align",
]
TEXT_OFFSETS: Mapping[str, tuple[str, int, int]] = {
    "right": ("left", 3, 8),
    "left": ("right", -5, 8),
    "b_center": ("center", 0, 15),
    "t_center": ("center", 0, 0),
}


def _text_position(theta: float, position: str = "outer") -> str:
    if position == "inner":
        theta -= pi
    theta %= 2 * pi
    if 4 * pi / 3 <= theta <= 5 * pi / 3:
        return "b_center"
    if pi / 3 <= theta <= 2 * pi / 3:
        return "t_center"
    if theta <= pi / 3 or theta >= 5 * pi / 3:
        return "right"
    return "left"


def _feature_glyph(feature: pd.Series) -> pd.Series:
    """Calculate polygon and leader-line coordinates for one feature."""
    end_angle = feature["rend"]
    start_angle = feature["rstart"]
    radius = BASE_RADIUS + FEATURE_THICKNESS * 2.3 * feature["level"]
    segment_length = end_angle - start_angle
    if feature["sframe"] == 1:
        end_angle, start_angle = start_angle, end_angle

    shift = pi / 2
    theta = np.linspace(
        shift - end_angle,
        shift - start_angle,
        int(25 * segment_length) + 3,
    )
    outer_x = (radius + FEATURE_THICKNESS) * np.cos(theta)
    outer_y = (radius + FEATURE_THICKNESS) * np.sin(theta)
    leader_angle = np.mean([end_angle, start_angle])
    if feature["has_orientation"] is True:
        outer_x = outer_x[:-2]
        outer_y = outer_y[:-2]
        leader_angle = np.arctan2(np.mean(outer_x), np.mean(outer_y))
        outer_x = np.append(outer_x, radius * np.cos(shift - start_angle))
        outer_y = np.append(outer_y, radius * np.sin(shift - start_angle))

    inner_x = (radius - FEATURE_THICKNESS) * np.cos(theta[::-1])
    inner_y = (radius - FEATURE_THICKNESS) * np.sin(theta[::-1])
    if feature["has_orientation"] is True:
        inner_x = inner_x[2:]
        inner_y = inner_y[2:]

    text_angle = (pi / 2) - leader_angle
    line_start_x = np.cos(text_angle) * (radius + FEATURE_THICKNESS)
    line_start_y = np.sin(text_angle) * (radius + FEATURE_THICKNESS)
    line_end_x = np.cos(text_angle) * radius * 1.3
    line_end_y = np.sin(text_angle) * radius * 1.3
    line_color = feature["fill_color"]
    if line_color == "#ffffff":
        line_color = feature["line_color"]
    return pd.Series(
        [
            list(np.hstack((outer_x, inner_x))),
            list(np.hstack((outer_y, inner_y))),
            line_end_x,
            line_end_y,
            line_color,
            [line_start_x, line_end_x],
            [line_start_y, line_end_y],
            text_angle,
            _text_position(float(text_angle)),
        ]
    )


def _sequence_markers(sequence_length: int) -> pd.DataFrame:
    chunk_size = max(500, round((sequence_length // 5) / 500) * 500)
    chunks = pd.Series(range(0, sequence_length - chunk_size // 2, chunk_size))
    chunks = chunks.loc[chunks.lt(sequence_length)].replace(0, 1)
    theta = (pi / 2) - (chunks / sequence_length) * 2 * pi
    offset = 0.155
    start_x = np.cos(theta) * offset
    start_y = np.sin(theta) * offset
    end_x = np.cos(theta) * offset / 1.08
    end_y = np.sin(theta) * offset / 1.08
    markers = pd.DataFrame(
        {
            "lineX": list(zip(start_x, end_x, strict=True)),
            "lineY": list(zip(start_y, end_y, strict=True)),
            "theta": theta,
            "bp": chunks,
            "Lx1": end_x,
            "Ly1": end_y,
            "size": "12px",
        }
    )
    markers["text_align"] = markers["theta"].apply(_text_position, position="inner")
    return markers


def _assign_feature_levels(annotations: pd.DataFrame) -> pd.DataFrame:
    annotations = annotations.sort_values(by="score", ascending=False).copy()
    if annotations.empty:
        annotations["level"] = pd.Series(dtype=int)
        return annotations

    intervals = annotations[["qstart", "qend", "qlen"]].copy()
    intervals["qstart"] = np.where(
        intervals["qstart"] >= intervals["qend"],
        intervals["qstart"] - intervals["qlen"],
        intervals["qstart"],
    )
    placed: list[tuple[Any, float, float, int]] = []
    for index, row in intervals.iterrows():
        start, end = row["qstart"], row["qend"]
        occupied = {
            level
            for _, placed_start, placed_end, level in placed
            if placed_start < end and start < placed_end
        }
        level = 0
        while level in occupied:
            level += 1
        placed.append((index, start, end, level))

    levels = pd.Series(
        {index: level for index, _, _, level in placed},
        name="level",
        dtype=int,
    )
    return annotations.join(levels)


def _prepare_plot_data(annotations: pd.DataFrame) -> pd.DataFrame:
    """Return a render-ready copy of an annotation table."""
    dataframe = _assign_feature_levels(annotations)
    if dataframe.empty:
        return dataframe

    dataframe["pi_permatch_int"] = (
        dataframe["pi_permatch"].astype(int).astype(str) + "%"
    )
    dataframe.loc[dataframe["db"] == "Rfam", "pi_permatch_int"] = ""
    dataframe["rstart"] = (dataframe["qstart"] / dataframe["qlen"]) * 2 * pi
    dataframe["rend"] = (dataframe["qend"] / dataframe["qlen"]) * 2 * pi
    dataframe["rstart"] %= 2 * pi
    dataframe["rend"] %= 2 * pi
    dataframe["rend"] = np.where(
        dataframe["rend"] < dataframe["rstart"],
        dataframe["rend"] + 2 * pi,
        dataframe["rend"],
    )
    dataframe["type"] = dataframe["type"].str.replace(
        "rep_origin", "origin of replication"
    )

    colors = pd.read_csv(_package_data.get_resource("data", "colors.csv"), index_col=0)
    colors = colors.rename(columns={"Type": "type"})
    fragment_colors = colors.copy()
    fragment_colors[["fill_color", "line_color"]] = fragment_colors[
        ["line_color", "fill_color"]
    ]
    fragment_colors["fill_color"] = "#ffffff"

    complete = dataframe.loc[~dataframe["fragment"]].merge(
        colors, how="left", on="type"
    )
    complete = complete.fillna(
        {"color": "grey", "fill_color": "#808080", "line_color": "#000000"}
    )
    fragments = dataframe.loc[dataframe["fragment"]].merge(
        fragment_colors, how="left", on="type"
    )
    fragments = fragments.fillna(
        {"color": "grey", "fill_color": "#ffffff", "line_color": "#808080"}
    )
    dataframe = pd.concat([complete, fragments], ignore_index=True)

    orientation = pd.read_csv(
        _package_data.get_resource("data", "feature_orientation.csv"),
        header=None,
        names=["type", "has_orientation"],
    )
    orientation["has_orientation"] = orientation["has_orientation"].eq("T")
    dataframe = dataframe.merge(orientation, on="type", how="left")
    dataframe["has_orientation"] = dataframe["has_orientation"].fillna(False)
    dataframe["type"] = dataframe["type"].str.replace("_", " ")
    dataframe[GLYPH_COLUMNS] = dataframe.apply(_feature_glyph, axis=1)

    allowed_types = set(colors["type"].str.replace("_", " "))
    dataframe["legend"] = dataframe["type"].where(
        dataframe["type"].isin(allowed_types), "misc feature"
    )
    return dataframe


def _new_figure() -> tuple[Any, HoverTool]:
    hover = HoverTool(
        tooltips=(
            '<font size="3"><b>@name</b> — @type   @pi_permatch_int</font> <br> @blurb'
        )
    )
    figure_factory: Any = figure
    range_factory: Any = Range1d
    plot = figure_factory(
        height=PLOT_DIMENSIONS,
        width=PLOT_DIMENSIONS,
        title="",
        toolbar_location=None,
        toolbar_sticky=False,
        match_aspect=True,
        sizing_mode="scale_width",
        tools=["save", "pan"],
        x_range=range_factory(
            -PLOT_SIZE, PLOT_SIZE, bounds=(-0.5, 0.5), min_interval=0.1
        ),
        y_range=range_factory(
            -PLOT_SIZE, PLOT_SIZE, bounds=(-0.5, 0.5), min_interval=0.1
        ),
    )
    plot.add_tools(hover)
    plot.toolbar.logo = None
    zoom = WheelZoomTool(zoom_on_axis=False)
    plot.add_tools(zoom)
    plot.toolbar.active_scroll = zoom
    plot.circle(
        x=0,
        y=0,
        radius=BASE_RADIUS,
        line_color="#000000",
        fill_color=None,
        line_width=2.5,
    )
    plot.axis.visible = False
    plot.grid.grid_line_color = "#EFEFEF"
    plot.outline_line_color = "#DDDDDD"
    return plot, hover


def _add_text_labels(
    plot: Any,
    dataframe: pd.DataFrame,
    *,
    text_field: str,
    alpha: float = 1.0,
    font_size_field: str | None = None,
    y_adjustments: Mapping[str, int] | None = None,
) -> None:
    for position, (alignment, x_offset, default_y_offset) in TEXT_OFFSETS.items():
        y_offset = (y_adjustments or {}).get(position, default_y_offset)
        options: dict[str, Any] = {}
        if font_size_field:
            options["text_font_size"] = font_size_field
        plot.text(
            x="Lx1",
            y="Ly1",
            x_offset=x_offset,
            y_offset=y_offset,
            text_align=alignment,
            text=text_field,
            alpha=alpha,
            level="overlay",
            source=ColumnDataSource(dataframe[dataframe["text_align"] == position]),
            **options,
        )


def _add_features(plot: Any, hover: HoverTool, dataframe: pd.DataFrame) -> None:
    source = ColumnDataSource(dataframe)
    renderer = plot.patches(
        "x",
        "y",
        fill_color="fill_color",
        line_color="line_color",
        name="features",
        line_width=2.5,
        source=source,
        legend_group="legend",
    )
    hover.renderers = [renderer]
    plot.multi_line(
        xs="lineX",
        ys="lineY",
        line_color="annoLineColor",
        line_width=3,
        level="overlay",
        line_cap="round",
        alpha=0.5,
        source=source,
    )
    _add_text_labels(plot, dataframe, text_field="name")


def _add_markers(plot: Any, sequence_length: int) -> None:
    markers = _sequence_markers(sequence_length)
    plot.multi_line(
        xs="lineX",
        ys="lineY",
        line_color="black",
        line_width=2,
        level="underlay",
        line_cap="round",
        alpha=0.5,
        source=ColumnDataSource(markers),
    )
    _add_text_labels(
        plot,
        markers,
        text_field="bp",
        alpha=0.5,
        font_size_field="size",
        y_adjustments={"right": 6, "left": 6, "t_center": -3},
    )
    plot.add_layout(
        Label(
            x=0,
            y=0,
            y_offset=-8,
            text_align="center",
            text=f"{sequence_length} bp",
            text_color="#7b7b7b",
            text_font_size="16px",
            level="overlay",
        )
    )


def get_bokeh(dataframe: pd.DataFrame, linear: bool = False) -> Any:
    """Create an interactive circular feature map without mutating input data."""
    plot, hover = _new_figure()
    if linear:
        line_length = BASE_RADIUS / 5
        plot.line(
            [0, 0],
            [BASE_RADIUS - line_length, BASE_RADIUS + line_length],
            line_width=4,
            level="overlay",
            line_color="black",
        )

    prepared = _prepare_plot_data(dataframe)
    if prepared.empty:
        return plot

    _add_features(plot, hover, prepared)
    sequence_length = int(prepared.iloc[0]["qlen"])
    _add_markers(plot, sequence_length)
    plot.legend.location = "bottom_left"
    plot.legend.border_line_color = "#EFEFEF"
    plot.legend.visible = True
    return plot
