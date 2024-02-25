# Import all packages
import sys
import pandas as pd
import importlib
import numpy as np

import pickle
from Bio import SeqIO
import gzip
import mappy as mp
import re

import holoviews as hv
import ninetysix as ns
import colorcet as cc
import bokeh.io
hv.extension('bokeh')

# Translate codon to amino acids
def translate(seq):
      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
        'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0: # Return empty if length is not divisible by 3
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

# Get mutated sequence by splitting template sequence
def get_mut(temp_seq, aa_seq):
    i = 0
    mut_ls = []
    for i in range(len(list(zip(temp_seq, aa_seq)))):
        if temp_seq[i] != aa_seq[i]:
            mut_ls.append(temp_seq[i]+str(i+1)+aa_seq[i])
        i = i + 1
    return mut_ls

# Function for making plate maps
def _make_platemap(df, title, cmap=None):
    """Generates a plate heatmap from evSeq data using Holoviews with
    bokeh backend.

    Called via `generate_platemaps`; see docs there.
    """
    # Convert SeqDepth to log for easier visualization.
    df['logseqdepth'] = np.log(
        df['Alignment Count'],
        out=np.zeros_like(df['Alignment Count'], dtype=float),
        where=(df['Alignment Count'] != 0)
    )
    
    # Create necessary Row and Column values and sort
    df = df.sort_values(['Column', 'Row'])
    df['Column'] = df['Column'].astype('str')

    # Set some base opts
    opts = dict(invert_yaxis=True, title=title, show_legend=True)
    
    # logseqdepth heatmap
    seq_depth_cmap = list(reversed(cc.CET_D9))
    
    # Set the center
    center = np.log(10)

    add_min = False
    if df['logseqdepth'].min() >= center:
        add_min = True

    # Adjust if it is greater than max of data (avoids ValueError)
    if df['logseqdepth'].max() <= center:
        # Adjust the center
        center = df['logseqdepth'].median()

    # center colormap
    if not add_min:
        color_levels = ns.viz._center_colormap(df['logseqdepth'], center)
    else:
        color_levels = ns.viz._center_colormap(
            list(df['logseqdepth']) + [np.log(1)],
            center
        )

    
    # Get heights
    n_rows = len(df['Row'].unique())
    n_cols = len(df['Column'].unique())
    height = int(50* n_rows)
    width = height * n_cols // n_rows

    # add tooltips
    tooltips = [
        ('Mutations', '@Mutations'),
        ('Alignment Count', '@Alignment Count'),
        ('Alignment Probability', '@Alignment Probability')
    ]

    def hook(plot, element):
        plot.handles['y_range'].factors = list('HGFEDCBA')
        plot.handles['x_range'].factors = [str(value) for value in range(1,13)]

    # generate the heatmap
    hm = hv.HeatMap(
        df,
        kdims=['Column', 'Row'],
        vdims=[
            'logseqdepth',
            'Mutations',
            'Alignment Count',
            'Alignment Probability',
        ],
    ).redim.values(
        row=np.unique(df["Row"]),
        Column=np.unique(df["Column"])
    ).opts(
        **opts,
        colorbar=True,
        cmap=seq_depth_cmap,
        height=height,
        width=width,
        line_width=4,
        clipping_colors={'NaN': '#DCDCDC'},
        color_levels=color_levels,
        tools=['hover'],
        colorbar_opts=dict(
            title='LogSeqDepth',
            background_fill_alpha=0
        ),
        hooks=[hook]
    )
 # function to bin the alignment frequencies into more relevant groupings
    def bin_align_freq(value):
        if value > 0.95:
            bin_vals = '0.95+'
        if value <= 0.95 and value > 0.9:
            bin_vals = '0.90-0.95'
        if value <= 0.9 and value > 0.8:
            bin_vals = '0.80-0.90'
        # anything below 0.8 should really be discarded
        if value < 0.8:
            bin_vals = '<0.80'
        if value == 0.0:
            bin_vals = '<0.80'
        return bin_vals
    
    # Bin alignment frequencies for easier viz
    bins = ['0.95+', '0.90-0.95', '0.80-0.90','<0.80']
    if cmap is None:
        cmap = [cc.bmy[int((1.1-i)*len(cc.bmy))]
                for i in [0.95, 0.9, 0.8, 0.4]]
    if 'stoplight' in cmap:
        cmap = ['#337D1F', '#94CD35', '#FFC300', '#C62C20']
    else:
        # Validate colormap
        if not isinstance(cmap, (list, tuple)):
            raise ValueError('cmap argument must be a list or tuple')
        if len(cmap) > 4:
            raise ValueError(
                'cmap argument has too many entries; only 4 should be passed'
            )
    cmap = {bin: color for bin, color in zip(bins, cmap)}

    # apply binning function to the AlignmentFrequency
    df['AlignmentProbabilityBinned'] = df['Alignment Probability'].apply(
        bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows*1.2

    # alignment frequency heatmap for edges around wells
    boxes = hv.Points(
        df.sort_values(['Alignment Probability'], ascending=False),
        ['Column', 'Row'],
        'AlignmentProbabilityBinned'
    ).opts(
        **opts,
        marker='square',
        line_color='AlignmentProbabilityBinned',
        line_join='miter',
        cmap=cmap,
        line_width=6,
        fill_alpha=0,
        line_alpha=1,
        legend_position='top',
        size=box_size,
    )
    
    # Use in apply statement for residue labels
    def split_variant_labels(mutation_string):
        
        num_mutations = len(mutation_string.split('_'))

        if  num_mutations > 4:
            return str(num_mutations)+' muts'

        mutation_string = mutation_string.replace('?','')
        new_line_mutations = mutation_string.replace('_','\n')
        
        return new_line_mutations
    
    _df = df.copy()
    _df['Labels'] = _df['Mutations'].apply(split_variant_labels)

    # Set the font size based on if #PARENT# is in a well and num of mutations
    max_num_mutations = _df['Labels'].apply(lambda x: len(x.split('\n'))).max()
    has_parent = ('#PARENT#' in _df['Labels'])
    
    if max_num_mutations > 3 or has_parent:
        label_fontsize = '8pt'
    else:
        label_fontsize = '8pt'

    labels = hv.Labels(
        _df,
        ['Column', 'Row'],
        'Labels',
    ).opts(
        text_font_size=label_fontsize,
        **opts,
        text_color = '#000000'
    )
    # return formatted final plot
    return (hm*boxes*labels).opts(frame_height=550,
                                  frame_width=550 * 3 // 2,
                                  border=50,
                                  show_legend=True)

# Main function to return heatmap
def generate_platemaps(
    max_combo_data,
    cmap=None,
    widget_location='top_left',
):
    """Saves a plate heatmap html generated from from evSeq data.
    
    Input:
    ------
    max_combo_data: path (str) or DartaFrame
        Path to 'Combos_Coupled_Max.csv' from an evSeq experiment or
        a pandas DataFrame of that file.
    cmap: list-like or str, default None
        The colormap to use for the well outline indicating alignment
        frequency. If None, defaults to a Plasma-like (colorcet.bmy)
        colormap. If 'stoplight', uses a green-yellow-red colormap (not 
        the most colorblind friendly, but highly intuitive). Otherwise
        you may pass any list -like object containing four colors (e.g.,
        ['#337D1F', '#94CD35', '#FFC300', '#C62C20'] for 'stoplight').
    widget_location: string, default 'top_left'
        Location of the widget for navigating plots. Must be one of:
        ['left', 'bottom', 'right', 'top', 'top_left', 'top_right',
        'bottom_left', 'bottom_right', 'left_top', 'left_bottom', 
        'right_top', 'right_bottom'].
    
    Returns:
    --------
    hm_holomap: an interactive Platemap
    """

    # Convert to dataframe if necessary
    if isinstance(max_combo_data, str):
        max_combo_df = pd.read_csv(max_combo_data)
    else:
        max_combo_df = max_combo_data.copy()
    
    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()
    
    # dictionary for storing plots
    hm_dict = {}
   
    # make logseqdepth column
    max_combo_df['logseqdepth'] = np.log(
        max_combo_df['Alignment Count'], 
        out=np.zeros_like(
            max_combo_df['Alignment Count'], 
            dtype=float
        ),
        where=max_combo_df['Alignment Count'] != 0
    )
    # Set the center
    center = np.log(10)

    add_min = False
    if max_combo_df['logseqdepth'].min() >= center:
        add_min = True

    # Adjust if it is greater than max of data (avoids ValueError)
    if max_combo_df['logseqdepth'].max() <= center:
        # Adjust the center
        center = max_combo_df['logseqdepth'].median()

    # center colormap
    if not add_min:
        color_levels = ns.viz._center_colormap(
            max_combo_df['logseqdepth'], center
        )
    else:
        color_levels = ns.viz._center_colormap(
            list(max_combo_df['logseqdepth']) + [np.log(1)],
            center
        )

    # Uniform color levels
    for _hm in hm_dict.values():
        _hm.opts({'HeatMap': {'color_levels': color_levels}})
    
    # Generate plots for each plate
    for plate in unique_plates:
        
        # Split to just the information of interest
        df = max_combo_df.loc[max_combo_df.Plate == plate].copy()
         
        # generate a holoviews plot
        hm_dict[plate] = _make_platemap(df, title=plate, cmap=cmap)  

    # plot from the dictionary
    hm_holomap = hv.HoloMap(
        hm_dict, 
        kdims=['Plate']
    )
    # Update widget location
    hv.output(widget_location=widget_location)

    return hm_holomap

