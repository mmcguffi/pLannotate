"""Tests for optional Bokeh plot rendering."""

import pandas as pd

from plannotate import bokeh_plot


def test_get_bokeh_supports_current_bokeh():
    dataframe = pd.read_csv("tests/test_data/pXampl3.csv")
    original = dataframe.copy(deep=True)

    plot = bokeh_plot.get_bokeh(dataframe)

    assert plot.width == 800
    assert plot.height == 800
    pd.testing.assert_frame_equal(dataframe, original)


def test_get_bokeh_handles_empty_annotations():
    columns = ["qstart", "qend", "score", "qlen"]

    assert bokeh_plot.get_bokeh(pd.DataFrame(columns=columns)) is not None
