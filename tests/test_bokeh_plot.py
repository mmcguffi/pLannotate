import os.path as op

import pandas as pd


def test_get_bokeh():
    from plannotate import bokeh_plot

    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)
    bokeh_plot.get_bokeh(df)
