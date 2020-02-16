import pandas as pd
from math import pi

from bokeh.io import push_notebook, show, output_notebook
from bokeh.palettes import Category20c
from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.transform import cumsum
from bokeh.models import HoverTool

def get_bokeh():

    X=0
    X2=X-.1
    Y=0

    featDesc=pd.read_csv("/Users/mattmcguffie/Documents/GitHub/pLannotate/feature_notes.csv",sep="\t")

    #TOOLTIPS = [
    #     ("index", "$index"),
    #    ("(x,y)", "($x, $y)"),
    #    ('Value', '$y'),
    #     ("desc", "@desc"),
    #     ('close',  '$@{adj close}{%0.2f}')
    #]
    #TOOLTIPS="@job: @value{0.2f} % insert desc here"

    TOOLTIPS='<font size="3"><b>@Feature</b> — @Type</font> <br> @Description'

    blue="#4e7fff"
    orange="#f6a35e"
    green="#479f71"
    grey="#808080"
    black="#000000"
    white="#ffffff"
    lineColor=white
    lineThick=3
    featThick = .02
    levelUp=featThick*2.25 #hacky -- change this so multi levels are supported when looped

    hover = HoverTool(names=["1","2"])
    plotSize=.4
    plotDimen=800
    p = figure(plot_height=plotDimen,plot_width=plotDimen, title="", toolbar_location=None, match_aspect=True, #sizing_mode='scale_width',
               tools=[hover,], tooltips=TOOLTIPS, x_range=(-plotSize, plotSize), y_range=(-plotSize, plotSize))

    #p.circle(x=X2, y=Y, size=408, line_color=grey, fill_color=None, line_width=5,level="glyph")
    # p.triangle(x=0.2, y=.2, size=35, angle=pi/4, color=orange, line_width=1,line_join='round')
    # p.triangle(x=-0.21, y=.2, size=35, angle=pi/4, color=blue, line_width=1,line_join='round')

    p.annular_wedge(x=X2, y=Y, inner_radius=.205-.004, outer_radius=.205+.004,
                    start_angle=0, end_angle=2*pi,line_color=None,fill_color=grey)

    p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick, outer_radius=.205+featThick, direction="anticlock",
                    start_angle=0, end_angle=pi/4,
            line_color=lineColor, line_width=lineThick,fill_color=orange, legend_group='Type',source=featDesc.iloc[[0]])

    p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick, outer_radius=.205+featThick, direction="anticlock",
                    start_angle=pi/2, end_angle=pi,
            line_color=lineColor, line_width=lineThick,fill_color=blue, legend_group='Type',source=featDesc.iloc[[1]])

    p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick+levelUp, outer_radius=.205+featThick+levelUp, direction="anticlock",
                    start_angle=pi/3, end_angle=pi/1.5,
            line_color=lineColor, line_width=lineThick,fill_color=green, legend_group='Type',source=featDesc.iloc[[2]])

    p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick, outer_radius=.205+featThick, direction="anticlock",
                    start_angle=1.5*pi, end_angle=1.9*pi,
            line_color=white, line_width=lineThick,fill_color=grey, legend_group='Type',source=featDesc.iloc[[3]])

    p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick, outer_radius=.205+featThick, direction="anticlock",
                start_angle=1.1*pi, end_angle=1.4*pi,
        line_color=grey, line_width=2, fill_color=white, legend_group='Type',source=featDesc.iloc[[4]])
    #
    # tSize=35
    # p.triangle(x=X2+0.008, y=.205, size=tSize, angle=pi/6, color=blue, line_width=1,line_join='round')
    # p.triangle(x=X2-0.008, y=-.205, size=tSize, angle=pi/2, color=grey, line_width=1,line_join='round')

    p.axis.axis_label=None
    p.axis.visible=False
    p.grid.grid_line_color = None
    p.outline_line_color = None
    p.legend.location = (230,325)
    p.legend.border_line_color=None

    return(p)
