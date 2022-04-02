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
    
def plasmid_features_viewer(gb_file, plot_width=900):
    """Gene feature viewer app"""
    
    if gb_file is None:
        return
    
    features = utils.genbank_to_features(gb_file)
    df = utils.features_to_dataframe(features)
    
    loc_pane = pnw.TextInput(name='location',value='',width=150)
    search_pane = pnw.TextInput(name='find_gene',value='',width=220)
    slider = pnw.IntSlider(name='start',start=0,end=10000,step=500,value=1,width=plot_width)
    xzoom_slider = pnw.IntSlider(name='zoom',start=1,end=500,value=100,step=5,width=100)
    left_button = pnw.Button(name='<',width=40)
    right_button = pnw.Button(name='>',width=40)
    fasta_seq = genbank_to_sequence(gb_file)
    feature_pane = pn.pane.Bokeh(height=100,margin=10)
    seq_pane = pn.pane.Bokeh(height=50, margin = 10)
    debug_pane = pn.pane.Str('debug',width=200,style={'background':'yellow','margin': '4pt'})
    
    seqlen = len(fasta_seq)
    slider.end = seqlen
    
    def search_features(event):
        """Find a feature"""
        
        term = search_pane.value        
        feats = utils.genbank_to_features(gb_file)
        df = utils.features_to_dataframe(feats)    
        df['gene'] = df.gene.fillna('')
        f = df[df.gene.str.contains(term)].iloc[0]
        #debug_pane.object = str(f.start)
        slider.value = int(f.start)-100
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
    
        sequence = fasta_seq[start: end]
        seq_pane.object = plotters.plot_sequence(sequence, plot_width, plot_height=50,fontsize='9pt',xaxis=False)            

        
    slider.param.watch(update,'value',onlychanged=True)
  
    xzoom_slider.param.watch(update,'value')       
    search_pane.param.watch(search_features,'value')    
    loc_pane.param.watch(update,'value',onlychanged=True)    

    feature_pane.object = plotters.plot_features(features, 0, 10000, plot_width=plot_width, tools="", rows=4)
    seq_pane.object = plotters.plot_sequence(fasta_seq, plot_width, plot_height=50,fontsize='9pt',xaxis=False)  
    top = pn.Row(loc_pane,xzoom_slider)
    main = pn.Column(feature_pane, seq_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main, sizing_mode='stretch_width',width_policy='max',margin=20)
    return app
# %%
plasmid_features_viewer('d378_attb-entry.gb')

# %%
