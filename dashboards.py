#!/usr/bin/env python
"""
    Implements viewers/panel apps for pybioviz
    Created June 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import os,sys,io
import numpy as np
import pandas as pd
from . import utils, plotters
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool)
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn
import panel.widgets as pnw

def sequence_alignment_viewer(filename=None):
    """Sequence alignment viewer"""
    
    title = pn.pane.Markdown()
    aln_btn = pnw.Button(name='align',width=100,button_type='primary')
    file_input = pnw.FileInput(name='load file',width=100,accept='.fa,.fasta,.faa')
    aligner_sel = pnw.Select(name='aligner',value='muscle',options=['muscle','clustal','maaft'],width=100)
    highlight_sel = pnw.Select(name='highlight mode',value='default',options=['default',''],width=100)
    rowheight_sl = pnw.IntSlider(name='row height',value=10,start=5,end=20,width=120)
    seq_pane = pn.pane.HTML(name='sequences',height=200,css_classes=['scrollingArea'])
    bokeh_pane = pn.pane.Bokeh()

    def update_title(filename):
        title.object = '### Sequence aligner: %s' %filename
        
    def update_file(event):        
        nonlocal seqtext
        seqtext = file_input.value.decode('utf-8')
        title.object = file_input.filename
        update_title(file_input.filename)
        return

    def align(event):
        #this function does the alignment
        nonlocal seqtext        
        if seqtext is not None:
            sequences = SeqIO.parse(io.StringIO(seqtext),format='fasta')
        elif filename is not None:    
            sequences = SeqIO.parse(filename, format='fasta')
        else:      
            return
        sequences = list(sequences)
        #print (sequences)
        aligner = aligner_sel.value
        if aligner == 'muscle':
            aln = utils.muscle_alignment(sequences)    
        elif aligner == 'clustal':
            aln = utils.clustal_alignment(sequences)
        elif aligner == 'mafft':
            aln = utils.mafft_alignment(sequences)
        if aln is None:
            bokeh_pane.object = plotters.plot_empty('%s not installed?' %aligner,900)
        else:
            #the bokeh pane is then updated with the new figure
            bokeh_pane.object = plotters.plot_sequence_alignment(aln, row_height=rowheight_sl.value)  
        return 
   
    seqtext = None
    file_input.param.watch(update_file,'value')
    rowheight_sl.param.watch(align,'value')
    aln_btn.param.watch(align, 'clicks')
    aln_btn.param.trigger('clicks')
    update_title(filename)
    side = pn.Column(aln_btn,file_input,aligner_sel,highlight_sel,rowheight_sl,seq_pane,css_classes=['form'],width=200,margin=20)   
    app = pn.Column(title,pn.Row(side, bokeh_pane), sizing_mode='stretch_width',width_policy='max',margin=20)
    return app



def plasmid_features_viewer(gff_file, ref_file=None, plot_width=900):
    """Gene feature viewer app"""
    
    if gff_file is None:
        return
    
    features = utils.gff_to_features(gff_file)
    df = utils.features_to_dataframe(features)
    utils.get_fasta_sequence(ref_file)
    loc_pane = pnw.TextInput(name='location',value='',width=150)
    slider = pnw.IntSlider(name='start',start=0,end=10000,step=500,value=1,width=plot_width)
    xzoom_slider = pnw.IntSlider(name='zoom',start=1,end=1000,value=500,step=5,width=100)  
    feature_pane = pn.pane.Bokeh(height=200,margin=10)
    seq_pane = pn.pane.Bokeh(height=100,margin=10)
    highlight_sel = pnw.Select(name='highlight mode',value='default',options=['default',''],width=100)
    rowheight_sl = pnw.IntSlider(name='row height',value=10,start=5,end=20,width=120)
    # seq_pane = pn.pane.HTML(name='sequences',height=200,css_classes=['scrollingArea'])
    
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
    if ref_file != None:
    #plot
    p = feature_pane.object = plotters.plot_features(features, 0, 10000, plot_width=plot_width, tools="", rows=4)
    
    #side = pn.Column(file_input,css_classes=['form'],width=200,margin=20)
    top = pn.Row(loc_pane,xzoom_slider)
    main = pn.Column(feature_pane, seq_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main, sizing_mode='stretch_width',width_policy='max',margin=20)
    return app

