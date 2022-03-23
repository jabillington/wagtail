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
    """Gene feature viewer app"""    
    
    if gff_file is None:
        return
    
    features = utils.gff_to_features(gff_file)
    df = utils.features_to_dataframe(features)
    utils.get_fasta_sequence(ref_file)
    loc_pane = pnw.TextInput(name='location',value='',width=150)
    slider = pnw.IntSlider(name='start',start=1,end=1000,step=500,value=1,width=plot_width)
    xzoom_slider = pnw.IntSlider(name='zoom',start=1,end=1000,value=500,step=5,width=100)  
    feature_pane = pn.pane.Bokeh(height=200,margin=10)
    seq_pane = pn.pane.Bokeh(height=100,margin=10)
    highlight_sel = pnw.Select(name='highlight mode',value='default',options=['default',''],width=100)
    rowheight_sl = pnw.IntSlider(name='row height',value=10,start=5,end=20,width=120)
    
    if ref_file is not None:
        seqlen = utils.get_fasta_length(ref_file)
        slider.end = seqlen
    else:
        slider.end = int(df.end.max())

    def search_features(event):
        """Find a feature"""
        
        term = search_pane.value        
        feats = utils.gff_to_features(gff_file)
        df = utils.features_to_dataframe(feats)    
        df['gene'] = df.gene.fillna('')
        f = df[df.gene.str.contains(term)].iloc[0]
        #debug_pane.object = str(f.start)
        slider.value = int(f.start)-100
        update(event)
        return   
    
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
            seq_pane.object = plotters.plot_sequence(sequence, plot_width, plot_height=50,fontsize='9pt',xaxis=False)            
        return
        
    slider.param.watch(update,'value',onlychanged=True)
    #slider.param.trigger('value')    
    xzoom_slider.param.watch(update,'value')       
    search_pane.param.watch(search_features,'value')    
    loc_pane.param.watch(update,'value',onlychanged=True)    
    left_button.param.watch(pan,'clicks')
    right_button.param.watch(pan,'clicks')

    p = feature_pane.object = plotters.plot_features(features, 0, 10000, plot_width=plot_width, tools="", rows=4)
    
    #side = pn.Column(file_input,css_classes=['form'],width=200,margin=20)
    top = pn.Row(loc_pane,xzoom_slider)
    main = pn.Column(feature_pane, seq_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main, sizing_mode='stretch_width',width_policy='max',margin=20)
    return app


# %%
import panel as pn
#need to load panel extension first
pn.extension(comms='vscode')
from pybioviz import dashboards
app = plasmid_features_viewer('example.gff')
app
# %%
