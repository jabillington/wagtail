# %%
# %% Load dependencies
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
from pybioviz.utils import genbank_to_features, get_fasta_length, get_fasta_sequence
from pybioviz import dashboards, utils, plotters
from importlib import reload
pn.extension(comms='vscode')

# %%

def genbank_to_sequence(gb_file, key=0):
    """Read genbank record features"""

    if gb_file is None or not os.path.exists(gb_file):
        return
    rec = list(SeqIO.parse(open(gb_file,'r'),'genbank'))[key]
    return rec.seq


def plot_sequence_mod(seq, plot_width=1000, plot_height=20, fontsize='10pt', xaxis=True, tools="", tbc =1):
    """Plot a single sequence.

    Args:
        seq: sequence to plot, a string
        xaxis: display x-axis tick labels or not
        tools: which bokeh tools to display, if needed
    """

    if seq is None or len(seq)==0:
        return plotters.plot_empty('no sequence',plot_width=plot_width,plot_height=plot_height)
    text = list(seq)
    N = len(seq)
    x = np.array(range(N))+1

    colors = utils.get_sequence_colors([seq])
    source = ColumnDataSource(dict(x=x, text=text, colors=colors))
    x_range = Range1d(0,N, bounds='auto')
    p = figure(plot_width=plot_width, plot_height=plot_height, x_range=x_range, y_range=(0,1),
               tools=tools, min_border=0, toolbar_location='below')
    rects = Rect(x="x", y=0,  width=tbc, height=2, fill_color="colors", line_color=None, fill_alpha=0.0)
    p.add_glyph(source, rects)

    glyph = Text(x="x", y=0, text="text", text_align='center', text_color="black",
                 text_font="monospace", text_font_size=fontsize)
    p.add_glyph(source, glyph)
    p.grid.visible = False
    if xaxis == False:
        p.xaxis.visible = False
    else:
        if plot_height<40:
            p.plot_height = 50
        if tools != "":
            p.plot_height = 70
    p.yaxis.visible = False
    p.toolbar.logo = None
    return p
    


def plasmid_features_viewer(gb_file, plot_width=900):
    """Gene feature viewer app"""
    
    if gb_file is None:
        return
    
    features = utils.genbank_to_features(gb_file)
   
    loc_pane = pnw.TextInput(name='location',value='',width=150)
    search_pane = pnw.TextInput(name='find_gene',value='',width=220)
    slider = pnw.IntSlider(name='start',start=0,end=10000,step=500,value=1,width=plot_width)
    xzoom_slider = pnw.IntSlider(name='zoom',start=1,end=500,value=100,step=5,width=100)
    fasta_seq = genbank_to_sequence(gb_file)
    feature_pane = pn.pane.Bokeh(height=100,margin=10)
    seq_pane = pn.pane.Bokeh(height=50, margin = 10)
    
    seqlen = len(fasta_seq)
    slider.end = seqlen
    
    def search_features(event):
        """Find a feature"""
        
        term = search_pane.value        
        feats = utils.genbank_to_features(gb_file)
        df = utils.features_to_dataframe(feats)    
        df['gene'] = df.gene.fillna('')
        f = df[df.gene.str.contains(term)].iloc[0]
        slider.value = int(f.start)-100
        update(event)
        return   
    
    def update(event):      
        print (event.obj.name)
        if event.obj.name in ['start', 'zoom']:
            xzoom = xzoom_slider.value*800
            start = int(slider.value)
            N = xzoom/16
            end = int(start+N)
            loc_pane.value = str(start)+':'+str(end)        
        elif event.obj.name == 'location':            
            vals = loc_pane.value.split(':')
            start = int(vals[0])
            end = int(vals[1])
            slider.value = start        
            

        p = feature_pane.object
        p.x_range.start = start
        p.x_range.end = end
    
        sequence = fasta_seq[start: end]

        if (end-start) <= 300:
            seq_pane.object = plot_sequence_mod(sequence, plot_width, plot_height=50,fontsize='pt',xaxis=False,tbc = 1)            
        else:
            seq_pane.object = plotters.plot_empty()





    slider.param.watch(update,'value',onlychanged=True) 
    xzoom_slider.param.watch(update,'value')    
    search_pane.param.watch(search_features,'value')
    loc_pane.param.watch(update,'value',onlychanged=True)    

    feature_pane.object = plotters.plot_features(features, 0, 10000, plot_width=plot_width, tools="", rows=4)
    seq_pane.object = plotters.plot_empty()
    top = pn.Row(loc_pane,xzoom_slider)
    main = pn.Column(feature_pane, seq_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main, sizing_mode='stretch_width',width_policy='max',margin=20)
    return app
# %%
plasmid_features_viewer('d378_attb-entry.gb')

# %%


def test_app():
    """Test dashboard"""
    
    def refresh(event):
        plot1.object = plotters.test1(cols=col_sl.value,rows=row_sl.value)
        plot2.object = plotters.test2(rows=row_sl.value)
        return
    
    title = pn.pane.Markdown('## pybioviz test plots')
    plot1 = pn.pane.Bokeh()    
    plot2 = pn.pane.Bokeh()
    col_sl = pnw.IntSlider(name='cols',value=30,start=5,end=200,step=1)
    col_sl.param.watch(refresh, 'value')
    row_sl = pnw.IntSlider(name='rows',value=10,start=5,end=100,step=1)
    row_sl.param.watch(refresh, 'value')
    col_sl.param.trigger('value')    
    app = pn.Column(title,col_sl,row_sl,plot1,plot2)
    return app

app = test_app()
app


# %%


def test_app():
    """Test dashboard"""
    
    def refresh(event):
        plot1.object = plotters.test1(cols=col_sl.value,rows=row_sl.value)
        plot2.object = plotters.test2(rows=row_sl.value)
        return
    
    title = pn.pane.Markdown('## pybioviz test plots')
    plot1 = pn.pane.Bokeh()    
    plot2 = pn.pane.Bokeh()
    col_sl = pnw.IntSlider(name='cols',value=30,start=5,end=200,step=1)
    col_sl.param.watch(refresh, 'value')
    row_sl = pnw.IntSlider(name='rows',value=10,start=5,end=100,step=1)
    row_sl.param.watch(refresh, 'value')
    col_sl.param.trigger('value')    
    app = pn.Column(title,col_sl,row_sl,plot1,plot2)
    return app

app = test_app()
app



# %%
