import os.path as op

import pandas as pd

from plannotate import resources


def test_get_bokeh():
    from plannotate import bokeh_plot

    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)
    bokeh_plot.get_bokeh(df)


def test_get_bokeh_empty_annotations():
    from plannotate import bokeh_plot

    plot = bokeh_plot.get_bokeh(pd.DataFrame(columns=resources.DF_COLS))

    assert plot is not None
