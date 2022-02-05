# %%
import os, sys, io, random
import string
import numpy as np
import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

from IPython.display import HTML

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, CustomJS, Slider, HoverTool, NumeralTickFormatter, Arrow, NormalHead
from bokeh.models import LinearAxis, Range1d, DataRange1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn
import panel.widgets as pnw
pn.extension()

from pybioviz import dashboards, utils, plotters
from importlib import reload


def plot_features(features, start=0, end=None, fontsize="8pt", plot_width=800, plot_height=150,
                  tools="xpan, xwheel_zoom, save", color='#abdda4', rows=3, key='gene'):
    """Bokeh sequence alignment view"""
    
    df = utils.features_to_dataframe(features)#, cds=True)    
    df = df[(df.type!='region') & (df['type']!='source')]   
    df['length'] = df.end-df.start
    #df['level'] = 1
    #df['color'] = random_colors(len(df)) #'green'
    df['x'] = df.start+df.length/2
    df = df.fillna('')
    def get_arrow(x):
        if x.strand == 1:
            return x.end
        else:
            return x.start

    df['arrow_start'] = df.apply(get_arrow,1)
    df['arrow_end'] = df.apply(lambda x: x.arrow_start+50 if x.strand==1 else x.arrow_start-50, 1)
    
    #def get_y(x):
    #    df['col'].shift()
    np.random.seed(8) 
    
    #df['y'] = np.random.randint(1,9, len(df))
    y = list(range(0,rows)) * len(df)
    df['y'] = y[:len(df)]
    
    if end == None:
        end = df.end.max()+100
        if end>10000:
            end=10000
    #print (df[:3])
    text = df.gene
    S = df.start.min()
    N = df.end.max()+10        
    x = list(df.start+df.length/2)
    h = 20

    source = ColumnDataSource(df)
    x_range = Range1d(start,end,min_interval=1)  
    viewlen = end-start
    hover = HoverTool(
        tooltips=[            
            ("gene", "@gene"),     
            ("locus_tag", "@locus_tag"),
            ("product", "@product"), 
            ("strand", "@strand"),
            ("length", "@length"),             
        ],        
    )  
    tools=[hover, tools]    
    #sequence text view with ability to scroll along x axis
    p = figure(title=None, plot_width=plot_width, plot_height=plot_height, x_range=x_range,
                y_range=(-1,rows), tools=tools, min_border=0, toolbar_location='right')#, lod_factor=1)
    #display text only at certain zoom level?
    #print (viewlen)
    if viewlen<20000:
        tags = Text(x="x", y="y", y_offset=-8, text=key, text_align='center',text_color="black", 
                     text_font="monospace",text_font_size=fontsize, name="genetext")
        p.add_glyph(source, tags)
    #rects
    rects = Rect(x="x", y="y", width="length", height=.4, fill_color=color, fill_alpha=0.4, name='rects')
    #arrows
    arr = Arrow(source=source, x_start="arrow_start", x_end="arrow_end", y_start="y", y_end="y", 
                line_color="black", name='arrows', end=NormalHead(size=8))
    p.add_glyph(source, rects)    
    p.add_layout(arr)
    
    p.grid.visible = False
    p.yaxis.visible = False
    p.xaxis.major_label_text_font_style = "bold"
    p.yaxis.minor_tick_line_width = 0
    p.yaxis.major_tick_line_width = 0
    p.toolbar.logo = None
    p.xaxis.formatter = NumeralTickFormatter(format="(0,0)")
    return p

feats = utils.gff_to_features('example.gff')
#feats = utils.genbank_to_features('1765.416.gbk',key=12)
#feats = utils.gff_to_features('Mbovis_AF212297.gff')
p = plot_features(features=feats, start=100)
pn.pane.Bokeh(p)



# %%
import panel as pn
#need to load panel extension first
pn.extension(comms='vscode')
from pybioviz import dashboards
app = dashboards.genome_features_viewer('example.gff')
app
# %%
