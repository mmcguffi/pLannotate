import pandas as pd
from math import pi

import pandas as pd
import numpy as np
from Bio import SeqIO

from bokeh.io import push_notebook, show, output_notebook
from bokeh.palettes import Category20c
from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.transform import cumsum
from bokeh.models import HoverTool

import streamlit as st

def get_bokeh(inRecord,plasLen):
    X=0
    X2=X-0
    Y=0

    featDesc=pd.read_csv("/Users/mattmcguffie/Documents/GitHub/pLannotate/feature_notes.csv",sep="\t",index_col=0)
    TOOLTIPS='<font size="3"><b>@Feature</b> — @Type</font> <br> @Description'

    blue="#4e7fff"
    orange="#f6a35e"
    green="#479f71"
    grey="#808080"
    black="#000000"
    white="#ffffff"
    lineColor=white
    lineThick=1
    featThick = .02
    levelUp=featThick*2.25 #hacky -- change this so multi levels are supported when looped

    hover = HoverTool(names=["1","2"])
    plotSize=.4
    plotDimen=800
    p = figure(plot_height=plotDimen,plot_width=plotDimen, title="", toolbar_location=None, match_aspect=True,sizing_mode='scale_width',
               tools=[hover,], tooltips=TOOLTIPS, x_range=(-plotSize, plotSize), y_range=(-plotSize, plotSize))

    #backbone line
    p.annular_wedge(x=X2, y=Y, inner_radius=.205-.004, outer_radius=.205+.004,
            start_angle=0, end_angle=2*pi,line_color=None,fill_color=grey)


    inRecord['rstart']=(pi/2) - ((inRecord["start"]/plasLen)*2*pi)
    inRecord['rend']=(pi/2) - ((inRecord["end"]/plasLen)*2*pi)

    inRecord=inRecord.join(featDesc)
    colorDf=pd.DataFrame(
        (("blue","#4e7fff","origin of replication","#ffffff"),
        ("orange","#f6a35e","promoter","#ffffff"),
        ("green","#479f71","CDS","#ffffff"),
        ("grey","#808080","misc feature","#ffffff"),
        ("white","#ffffff","primer bind","#808080"),
        ("black","#000000",None,"#ffffff")),
        columns=["color","hex","Type","line_color"])
    df=inRecord.merge(colorDf,on=["Type"])

    for index in df.index:
        p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick, outer_radius=.205+featThick, direction="clock",
                start_angle=df.loc[index]['rstart'], end_angle=df.loc[index]['rend'],
                line_color="line_color", line_width=lineThick,fill_color='hex', legend_group='Type',source=df.loc[[index]])

        rstart=df.loc[index]['rstart']
        rend=df.loc[index]['rend']
        theta=(rstart+rend)/2
        if rstart<rend:
            theta+=pi
        normRadius=.205
        longRadius=.31

        Lx0=np.cos(theta)*normRadius
        Ly0=np.sin(theta)*normRadius

        Lx1=np.cos(theta)*longRadius
        Ly1=np.sin(theta)*longRadius

        p.text(x=Lx1, y=Ly1,x_offset=0,y_offset=-Ly1*50,text_align="center",
                         text='Feature', level="annotation", source=df.loc[[index]])
        p.line(x=[Lx0,Lx1], y=[Ly0,Ly1], line_color=grey, line_width=3,level="underlay",line_cap='round',alpha=.5)

    p.axis.axis_label=None
    p.axis.visible=False
    p.grid.grid_line_color = None
    p.outline_line_color = None
    p.legend.location = (230,325)
    p.legend.border_line_color=None
    p.legend.visible=False

    st.write(df)
    return p
