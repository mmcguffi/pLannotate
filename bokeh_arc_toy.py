
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid
from bokeh.models.glyphs import Patch
from bokeh.plotting import figure, output_file, show

import numpy as np
from math import pi
plotSize=.35
plotDimen=800

p = figure(
    title=None, plot_width=800, plot_height=800,
    min_border=0, toolbar_location=None,match_aspect=True, x_range=(-plotSize, plotSize), y_range=(-plotSize, plotSize))

def draw_acr(r2,r1,frame):
    baseRadius=.205
    thickness=.015
    segLen=r1-r2
    if frame==-1: #reverses for direction
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

    p.patch(x=x, y=y, fill_color="#a6cee3")
    return p
plot=draw_acr(0.0992,1.016382261,1)
plot=draw_acr(1.0446,1.3412485344,-1)

show(plot)
