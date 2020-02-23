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
from bokeh.models import HoverTool, ColumnDataSource

import streamlit as st

def draw_acr(inSeries,p):
    r1=inSeries['rend']
    r2=inSeries['rstart']
    frame=inSeries['sframe']
    baseRadius=.205
    thickness=.015
    segLen=r1-r2
    if frame==1: #reverses for direction
        r1,r2=r2,r1
    shift=pi/2
    N=int(25*segLen)+3 #number of lines/sampling size
    theta = np.linspace(shift-r1, shift-r2, N) #regularly samples between space
    x1=(baseRadius+thickness)*np.cos(theta)
    y1=(baseRadius+thickness)*np.sin(theta)
    x1=x1[:-2] #pops last 2 lines so arrow can be drawn
    y1=y1[:-2]

    #adds a single point between both arcs (arrow tip)
    x1=np.append(x1,(baseRadius)*np.cos(shift-r2))
    y1=np.append(y1,(baseRadius)*np.sin(shift-r2))

    #creates second arc (reversed so closed polygon can be drawn)
    x2=(baseRadius-thickness)*np.cos(theta[::-1])
    y2=(baseRadius-thickness)*np.sin(theta[::-1])
    x2=x2[2:]
    y2=y2[2:]

    #honestly I need to learn what does this better
    x = np.hstack((x1, x2))
    y = np.hstack((y1, y2))
    source = ColumnDataSource(dict(x=x, y=y))

    p.patch(x=x, y=y, fill_color=inSeries['hex'],line_color=inSeries['line_color'])

def get_bokeh(inRecord):
    X=0
    X2=X-0
    Y=0

    featDesc=pd.read_csv("/Users/mattmcguffie/Documents/GitHub/pLannotate/feature_notes.csv",sep="\t",index_col=0)
    TOOLTIPS='<font size="3"><b>@Feature</b> — @Type   @pi_permatch_int% </font> <br> @Description'

    blue="#4e7fff"
    orange="#f6a35e"
    green="#479f71"
    grey="#808080"
    black="#000000"
    white="#ffffff"
    lineColor=white
    lineThick=0
    featThick = .015
    levelUp=featThick*2.25 #hacky -- change this so multi levels are supported when looped

    hover = HoverTool(names=["1","2"])
    plotSize=.35
    plotDimen=800
    p = figure(plot_height=plotDimen,plot_width=plotDimen, title="", toolbar_location=None, match_aspect=True,sizing_mode='scale_width',
               tools=[hover,], tooltips=TOOLTIPS, x_range=(-plotSize, plotSize), y_range=(-plotSize, plotSize))

    #backbone line
    p.annular_wedge(x=X2, y=Y, inner_radius=.205-.001, outer_radius=.205+.001,
            start_angle=0, end_angle=2*pi,line_color=None,fill_color=black)

    inRecord['rstart']=((inRecord["qstart"]/inRecord["qlen"])*2*pi)
    inRecord['rend']  =((inRecord["qend"]/inRecord["qlen"])*2*pi)
    inRecord['rstart'] = np.where(inRecord['rstart'] < 0, inRecord['rstart'] + (2*pi), inRecord['rstart'])
    inRecord['rend'] = np.where(inRecord['rend'] < 0, inRecord['rend'] + (2*pi), inRecord['rend'])
    inRecord['rend'] = np.where(inRecord['rend'] < inRecord['rstart'], inRecord['rend'] + (2*pi), inRecord['rend'])
    # shorten=np.arctan(featThick/.205)
    # inRecord['arcLen']=np.abs((inRecord['rstart'] - inRecord['rend']))
    # inRecord['shorten'] = np.where((shorten >= np.abs(inRecord['arcLen'])), 0, shorten)
    # inRecord['rend']=inRecord['rend']+inRecord['shorten']

    st.write(inRecord)

    inRecord=inRecord.join(featDesc)

    fullColorDf=pd.read_csv("./colors.csv",index_col=0)
    fragColorDf=fullColorDf.copy()
    fragColorDf[['hex','line_color']]=fragColorDf[['line_color','hex']]
    fragColorDf["hex"]="#ffffff"

    fullColorDf.to_csv("~/Desktop/colors.csv")

    full=inRecord[inRecord["fragment"]==False]
    full=full.merge(fullColorDf,how = "left",on=["Type"])
    full=full.fillna({"color":"grey","hex":"#808080","line_color":"#000000"})
    st.write("full")
    st.write(full)

    frag=inRecord[inRecord["fragment"]==True]
    frag=frag.merge(fragColorDf,how = "left",on=["Type"])
    frag=frag.fillna({"color":"grey","hex":"#ffffff","line_color":"#808080"})

    df=full.append(frag).reset_index(drop=True)#.set_index('Feature')

    st.write(df)

    for index in df.index:
        #if df.loc[index]['type']=="promoter":continue
        if df.loc[index]['type']=="primer_bind":continue
        # p.annular_wedge(x=X2, y=Y, name="1", inner_radius=.205-featThick,
        #         outer_radius=.205+featThick, direction="clock", start_angle='rstart',
        #         end_angle='rend', line_color='line_color', line_width=lineThick,
        #         fill_color='hex', legend_group='Type',source=df.loc[[index]])

        draw_acr(df.loc[index],p)
        # #draw arrow
        # #############
        # # inRecord['rstart']=(pi/2) - ((inRecord["qstart"]/inRecord["qlen"])*2*pi)
        # #############
        # x=[-.1,-.1,0, 0, featThick]
        # y=[-featThick, featThick,-featThick, featThick, 0]
        # arrowTipLen=featThick
        # x=[-.002,-.002, 0, arrowTipLen,0]
        # y=[-featThick, featThick, featThick, 0,-featThick]
        # l=[x,y]
        # arrow_theta=(pi/2)-rend
        # rotate= [[np.cos(arrow_theta),np.sin(arrow_theta)],
        #         [-np.sin(arrow_theta),np.cos(arrow_theta)]]
        # l=np.dot(rotate,l)
        #
        # x=l[0]
        # y=l[1]
        # xTrans=np.cos(rend)*normRadius
        # yTrans=np.sin(rend)*normRadius
        # x=x+xTrans
        # y=y+yTrans
        #
        # p.patch(x,y, color=df.loc[index]['hex'],alpha=1,line_color=None, line_width=2,level='overlay')
        # ##############

        rstart=df.loc[index]['rstart']
        rend=df.loc[index]['rend']
        theta=(pi/2)-(rstart+rend)/2

        normRadius=.205
        longRadius=.27

        Lx0=np.cos(theta)*normRadius
        Ly0=np.sin(theta)*normRadius

        Lx1=np.cos(theta)*longRadius
        Ly1=np.sin(theta)*longRadius


        lineColor=df.loc[index]['hex']
        if lineColor == "#ffffff":
            lineColor="#808080"

        p.text(x=Lx1, y=Ly1,name="1",x_offset=Lx1*100,y_offset=-Ly1*50,text_align="center",
                         text='Feature', level="annotation", source=df.loc[[index]])
        p.line(x=[Lx0,Lx1], y=[Ly0,Ly1], line_color=lineColor, line_width=3,level="underlay",line_cap='round',alpha=.5)

    p.axis.axis_label=None
    p.axis.visible=False
    p.grid.grid_line_color = None
    p.outline_line_color = None
    #p.legend.location = (230,325)
    p.legend.border_line_color=None
    p.legend.visible=False

    #df.to_csv("~/Desktop/test.csv")
    #st.write(df)
    return p
