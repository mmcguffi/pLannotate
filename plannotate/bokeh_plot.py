from math import pi

from Bio import SeqIO
from bokeh.models import HoverTool, ColumnDataSource, WheelZoomTool, Range1d, Legend, LegendItem
from bokeh.models.annotations import Label
from bokeh.plotting import figure
import numpy as np
import pandas as pd

import plannotate.resources as rsc

global baseRadius
baseRadius=.18

def text_pos(theta,pos="outer"):
    if pos == "inner":
        theta-=pi
    if theta<0: theta+=2*pi
    trQ=pi/3;tlQ=2*pi/3;blQ=4*pi/3;brQ=5*pi/3
    if theta >= blQ and theta <= brQ:
        annoPos="b_center" #bottom
    elif theta >= trQ and theta <= tlQ:
        annoPos="t_center" #top
    elif theta <= trQ or theta >= brQ:
        annoPos="right"
    else:
        annoPos="left"
    return annoPos

def calc_glyphs(inSeries):
    r1=inSeries['rend']
    r2=inSeries['rstart']
    frame=inSeries['sframe']

    level = inSeries['level']

    thickness=.017
    featRadius=float(baseRadius)
    featRadius+=thickness*2.3*level

    segLen=r1-r2
    if frame==1: #reverses for direction
        r1,r2=r2,r1
    shift=pi/2 #corrects for starting at the correct polar space
    N=int(25*segLen)+3 #number of lines/sampling size
    theta = np.linspace(shift-r1, shift-r2, N) #regularly samples between space
    x1=(featRadius+thickness)*np.cos(theta) #x=r*cos(θ), classic polar eq
    y1=(featRadius+thickness)*np.sin(theta)
    
    # calculates the angle between segments
    line_theta_avg = np.mean([r1,r2])

    if inSeries['has_orientation'] == True:
        x1=x1[:-2] #pops last 2 lines so arrow can be drawn
        y1=y1[:-2]
        
        # calculates the angle between the non-arrow segments
        # over-rides previous calcuation
        line_theta_avg = np.arctan2(np.mean(x1), np.mean(y1))

        #adds a single point between both arcs (arrow tip)
        x1=np.append(x1,(featRadius)*np.cos(shift-r2))
        y1=np.append(y1,(featRadius)*np.sin(shift-r2))

    #creates second arc (reversed so closed polygon can be drawn)
    x2=(featRadius-thickness)*np.cos(theta[::-1])
    y2=(featRadius-thickness)*np.sin(theta[::-1])

    #trims bottom part of arrow
    if inSeries['has_orientation'] == True:
        x2=x2[2:]
        y2=y2[2:]

    #concatenates the two -- could change to make more readible
    x = np.hstack((x1, x2))
    y = np.hstack((y1, y2))
    x=list(x)
    y=list(y)

    #calculate text placement/lines
    theta=(pi/2)-line_theta_avg

    Lx0=np.cos(theta)*(featRadius + thickness)
    Ly0=np.sin(theta)*(featRadius + thickness) 
    longRadius=featRadius * 1.3
    Lx1=np.cos(theta)*longRadius
    Ly1=np.sin(theta)*longRadius

    lineX=[Lx0,Lx1]
    lineY=[Ly0,Ly1]

    annoLineColor=inSeries['fill_color']
    if annoLineColor == "#ffffff":
        annoLineColor=inSeries['line_color']

    annoPos = text_pos(theta)

    return pd.Series([x,y,Lx1,Ly1,annoLineColor,lineX,lineY,theta,annoPos])

def calc_num_markers(plasLen):
    #calculate chunk size(s) and positions for drawing lines
    chunkSize = round((plasLen//5)/500)*500
    if chunkSize == 0: chunkSize = 500
    #chunkSize = int(np.ceil((plasLen//5)/500)*500)
    chunks = pd.Series(range(0, plasLen-int(chunkSize/2), int(chunkSize)))
    chunks = chunks[chunks<plasLen]
    chunks = chunks.replace(0,1)
    chunksR = (chunks/plasLen) * 2 * pi

    theta=(pi/2)-chunksR #rotate for canvas

    offset = .155
    Lx0 = np.cos(theta) * offset
    Ly0 = np.sin(theta) * offset
    longRadius = offset / 1.08
    Lx1 = np.cos(theta) * longRadius
    Ly1 = np.sin(theta) * longRadius

    lineX=list(zip(Lx0,Lx1))
    lineY=list(zip(Ly0,Ly1))
    ticks = pd.DataFrame(list(zip(lineX,lineY,theta)),columns=['lineX','lineY','theta'])
    ticks['text_align'] = ticks['theta'].apply(text_pos,pos='inner')
    ticks['bp'] = chunks
    ticks['Lx1'] = Lx1
    ticks['Ly1'] = Ly1
    ticks['size'] = "12px"

    return ticks


def calc_level(inDf):
    # calculates the level to be rendered at
    # highest-scoring hits are priority level 0 (on plasmid "ring")
    # if a level is already occupied, chooses next higher ring
    inDf = inDf.sort_values(by = ['score'], ascending = [False])
    levels = inDf[['qstart','qend','score','qlen']].copy()#.sort_values(by = ['qlen'], ascending = [False])
    levels['qstart'] = np.where(levels['qstart'] >= levels['qend'], levels['qstart'] - levels['qlen'], levels['qstart'])
    
    calculated_levels = pd.DataFrame(columns = ['index', 's', 'e', 'level'])
    for index in levels.index:
        s = levels.loc[index]['qstart']
        e = levels.loc[index]['qend']
        intervals = pd.arrays.IntervalArray.from_arrays(calculated_levels['s'].to_list(), calculated_levels['e'].to_list())
        #overlap = calculated_levels['s'].between(s,e) | calculated_levels['e'].between(s,e)
        overlap = calculated_levels[intervals.overlaps(pd.Interval(s, e))]
        
        if overlap.empty:
            new_level = 0
        else:
            new_level = 0
            #iterates until lowest possible level is available
            while new_level in set(overlap['level']):
                new_level += 1
        
        # calculated_levels = calculated_levels.append({'index' : index, 's' : s, 'e' : e, 'level' : new_level},
        # ignore_index = True)
        calculated_levels = pd.concat([
            calculated_levels, 
            pd.DataFrame([[index, s, e, new_level]], columns = ['index', 's', 'e', 'level']), 
            ],
            ignore_index = True)
            
    inDf = inDf.join(calculated_levels[['level']])
    
    ################################################################
    #weird error where sometimes the level is NaN
    inDf['level'] = inDf['level'].fillna(3)
    ################################################################
    
    # inDf['level']=None
    # for i in inDf.index:
    #     df=inDf[inDf.index<i]
    #     s=inDf.loc[i]['qstart_dup']
    #     e=inDf.loc[i]['qend_dup']
    #     startBound=((df['qstart_dup']<=s) & (df['qend_dup']>=s))
    #     endBound=((df['qstart_dup']<=e) & (df['qend_dup']>=e))
    #     occupied_levels = list(set(df[startBound]['level']) | set(df[endBound]['level']))

    #     new_level = 0
    #     while new_level in occupied_levels:
    #         new_level += 1
    #     inDf.at[i,'level'] = new_level
    return inDf

    # # classic interval scheduling algorithm
    # inDf = inDf.sort_values(by="qend")
    # inDf = inDf.reset_index()
    # levels = {0:0}

    # for i in range(1,len(inDf)):
    #     if inDf.iloc[i]['qstart'] > inDf.iloc[i-1]['qend']:
    #         levels[i] = 0
    #     else: #while loop needed here for deeply nested?
    #         if (levels[i-2] < levels[i-1]) and (inDf.iloc[i]['qstart'] > inDf.iloc[i-2]['qend']):
    #             levels[i] = levels[i-1] - 1
    #         else:
    #             levels[i] = levels[i-1] + 1
    
    # levels = pd.DataFrame(levels.values(),columns = ["level"])
    # inDf = inDf.join(levels)
    # inDf = inDf.set_index("index")
    # inDf.index.name = None
    # inDf = inDf.sort_index()

    # return inDf

def get_bokeh(df, linear = False):
    
    #df = df.fillna("")

    X=0
    Y=0

    TOOLTIPS='<font size="3"><b>@Feature</b> — @Type   @pi_permatch_int</font> <br> @Description'

    hover = HoverTool(name="features")
    plotSize=.45
    plotDimen=600

    x_range = Range1d(-plotSize, plotSize, bounds=(-.5, .5), min_interval=.1)
    y_range = Range1d(-plotSize, plotSize, bounds=(-.5, .5), min_interval=.1)
    toolbar = None
    p = figure(frame_height=plotDimen,frame_width=plotDimen, title="",
                toolbar_location=toolbar, toolbar_sticky=False, match_aspect=True,
                sizing_mode='scale_width', tools=['save',hover,'pan'], tooltips=TOOLTIPS,
                #x_range=(-plotSize, plotSize), y_range=(-plotSize, plotSize))
                x_range=x_range, y_range=y_range)
    p.toolbar.logo = None
    p.add_tools(WheelZoomTool(zoom_on_axis=False))
    p.toolbar.active_scroll = p.select_one(WheelZoomTool)

    #backbone line
    p.circle(x=X, y=Y, radius=baseRadius, line_color="#000000", fill_color=None, line_width=2.5)

    df = calc_level(df)
    
    if linear:
        line_length = baseRadius / 5
        p.line([0,0], [baseRadius - line_length, baseRadius + line_length],
                line_width = 4, level="overlay", line_color = "black")

    df['pi_permatch_int']=df['pi_permatch'].astype('int')
    df['pi_permatch_int'] = df['pi_permatch_int'].astype(str) + "%"
    df.loc[df['db'] == "Rfam", 'pi_permatch_int'] = "" #removes percent from infernal hits

    df['rstart']=((df["qstart"]/df["qlen"])*2*pi)
    df['rend']  =((df["qend"]/df["qlen"])*2*pi)
    df['rstart'] = np.where(df['rstart'] < 0, df['rstart'] + (2*pi), df['rstart'])
    df['rend'] = np.where(df['rend'] < 0, df['rend'] + (2*pi), df['rend'])
    df['rend'] = np.where(df['rend'] < df['rstart'], df['rend'] + (2*pi), df['rend'])

    df['Type'] = df['Type'].str.replace('rep_origin','origin of replication')
    
    #DDE0BD
    #C97064
    #C9E4CA
    fullColorDf=pd.read_csv(rsc.get_resource("data", "colors.csv"),index_col=0)
    fragColorDf=fullColorDf.copy()
    fragColorDf[['fill_color','line_color']]=fragColorDf[['line_color','fill_color']]
    fragColorDf["fill_color"]="#ffffff"

    full=df[df["fragment"]==False]
    full=full.merge(fullColorDf,how = "left",on=["Type"])
    full['legend'] = full['Type']
    full=full.fillna({"color":"grey","fill_color":"#808080","line_color":"#000000"})

    frag=df[df["fragment"]==True]
    frag=frag.merge(fragColorDf,how = "left",on=["Type"])
    frag=frag.fillna({"color":"grey","fill_color":"#ffffff","line_color":"#808080"})

    # df=full.append(frag).reset_index(drop=True) #.set_index('Feature')
    df = pd.concat([full,frag],ignore_index=True).reset_index(drop=True)

    # add orientation column
    orient = pd.read_csv(rsc.get_resource("data", "feature_orientation.csv"),header=None, names = ["Type","has_orientation"])
    orient['Type'] = orient['Type']
    orient['has_orientation'] = orient['has_orientation'].map({"T":True})
    df = df.merge(orient, on="Type", how = "left")
    df['Type'] = df['Type'].str.replace("_"," ")
    df['has_orientation'] = df['has_orientation'].fillna(value=False)

    df[['x','y',"Lx1","Ly1","annoLineColor","lineX","lineY","theta","text_align"]]=df.apply(calc_glyphs,axis=1)

    df['legend'] = df['Type']
    #allowedTypes = ['CDS',"promoter","origin of replication","swissprot"]
    allowedTypes = fullColorDf['Type']
    mask = ~df['legend'].isin(allowedTypes)
    df.loc[mask, 'legend'] = 'misc feature'

    #plot annotations
    source = ColumnDataSource(df)
    p.patches('x', 'y', fill_color='fill_color', line_color='line_color',
            name="features", line_width=2.5, source=source, legend_group = "legend")
    p.multi_line(xs="lineX", ys="lineY", line_color="annoLineColor", line_width=3,
            level="overlay",line_cap='round',alpha=.5, source=source)

    #`text_align` cannot read from `source` -- have to do this workaround
    right=ColumnDataSource(df[df['text_align']=='right'])
    left=ColumnDataSource(df[df['text_align']=='left'])
    bCenter=ColumnDataSource(df[df['text_align']=='b_center'])
    tCenter=ColumnDataSource(df[df['text_align']=='t_center'])

    text_level = 'overlay'
    p.text(x="Lx1", y="Ly1",name="2",x_offset=3,y_offset=8, text_align="left",
            text='Feature', level=text_level, source=right)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=-5,y_offset=8, text_align="right",
            text='Feature', level=text_level, source=left)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=0,y_offset=15, text_align="center",
            text='Feature', level=text_level, source=bCenter)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=0,y_offset=0, text_align="center",
            text='Feature', level=text_level, source=tCenter)

    #calculate chunk size(s) for drawing lines
    plasLen = df.iloc[0]['qlen']
    ticks = calc_num_markers(plasLen)
    ticks_cds = ColumnDataSource(ticks)
    p.multi_line(xs="lineX", ys="lineY", line_color="black", line_width=2,
                level="underlay", line_cap='round', alpha = .5, source = ticks_cds)

    right=ColumnDataSource(ticks[ticks['text_align']=='right'])
    left=ColumnDataSource(ticks[ticks['text_align']=='left'])
    bCenter=ColumnDataSource(ticks[ticks['text_align']=='b_center'])
    tCenter=ColumnDataSource(ticks[ticks['text_align']=='t_center'])
    p.text(x="Lx1", y="Ly1",name="2",x_offset=3,y_offset=6, text_align="left",
            text='bp', alpha = .5, text_font_size = 'size',level=text_level, source=right)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=-5,y_offset=6, text_align="right",
            text='bp', alpha = .5, text_font_size = 'size',level=text_level, source=left)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=0,y_offset=15, text_align="center",
            text='bp', alpha = .5, text_font_size = 'size',level=text_level, source=bCenter)
    p.text(x="Lx1", y="Ly1",name="2",x_offset=0,y_offset=-3, text_align="center",
            text='bp', alpha = .5, text_font_size = 'size', level=text_level, source=tCenter)

    p.add_layout(Label(x=0, y=0,name="2",x_offset=0,y_offset=-8, text_align="center",
            text=f"{plasLen} bp", text_color = "#7b7b7b", text_font_size = '16px', level=text_level))

    p.axis.axis_label=None
    p.axis.visible=False
    p.grid.grid_line_color = "#EFEFEF"
    p.outline_line_color = "#DDDDDD"
    p.legend.location = 'bottom_left'
    p.legend.border_line_color = "#EFEFEF"
    p.legend.visible = True

    #st.write(df)

    return p
