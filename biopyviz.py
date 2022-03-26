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

from pybioviz.pybioviz.utils import genbank_to_features
pn.extension()

from pybioviz import dashboards, utils, plotters
from importlib import reload



# %%
import panel as pn
#need to load panel extension first
pn.extension(comms='vscode')
from pybioviz import dashboards
app = dashboards.genome_features_viewer('example.gff')
app
# %%


# %%
def plasmid_features_viewer(gff_file, ref_file=None, plot_width=900):
    """Plasmid feature viewer app"""
    
    if gff_file is None:
        return
    
    features = utils.genbank_to_features(gff_file)
    df = utils.features_to_dataframe(features)
    
    loc_pane = pnw.TextInput(name='location',value='',width=150)
    slider = pnw.IntSlider(name='start',start=0,end=10000,step=500,value=1,width=plot_width)
    xzoom_slider = pnw.IntSlider(name='zoom',start=1,end=500,value=100,step=5,width=100)
    chrom_select = pnw.Select(name='chrom',width=220)
    feature_pane = pn.pane.Bokeh(height=100,margin=10)
    seq_pane = pn.pane.Bokeh(height=50,margin=10)
    
    if ref_file is not None:
        seqlen = utils.get_fasta_length(ref_file)
        slider.end = seqlen
    else:
        slider.end = int(df.end.max())
    
    def pan(event):
        p = feature_pane.object
        rng = p.x_range.end-p.x_range.start        
        inc = int(rng/10)
        print (event.obj.name)
        if event.obj.name == '<':
            slider.value = int(slider.value) - inc        
        else:
            slider.value = int(slider.value) + inc   
        update(event)
        return
    
    def update(event):      
        print (event.obj.name)
        if event.obj.name in ['start', 'zoom']:
            xzoom = xzoom_slider.value*200
            start = int(slider.value)
            N = xzoom/2
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
        if ref_file:
            sequence = utils.get_fasta_sequence(ref_file, start, end)
            seq_pane.object = plotters.plot_sequence(sequence, plot_width, plot_height=100,fontsize='20pt',xaxis=TRUE)            
        return
        
    slider.param.watch(update,'value',onlychanged=True)
    xzoom_slider.param.watch(update,'value')           
    loc_pane.param.watch(update,'value',onlychanged=True)    
    if ref_file != None:
        chrom_select.options = utils.get_fasta_names(ref_file)
    #plot
    p = feature_pane.object = plotters.plot_features(features, 0, 10000, plot_width=plot_width, tools="", rows=4)
    
    top = pn.Row(xzoom_slider)
    main = pn.Column(feature_pane, seq_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main, sizing_mode='stretch_width',width_policy='max',margin=20)
    return app

# %%
import panel as pn
import panel.widgets as pnw
pn.extension()
from pybioviz import dashboards, utils, plotters
from importlib import reload

#need to load panel extension first
pn.extension(comms='vscode')
from pybioviz import dashboards
app = plasmid_features_viewer('d378_attb-entry.gb')
app

# %%
utils.genbank_to_features('d378_attb-entry.gb')
# %%
