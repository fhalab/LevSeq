#!/usr/bin/env python
# coding: utf-8

# # Every Variant Sequencing with Oxford Nanopore Technologies
# 
# This script is being used after sequencing. The raw pod5 files can be basecalled or the already basecalled files can be used directly (fastq.gz)
# 
# ## Workflow
# 
# ### 1. Basecalling (Optional)
# 
# - The raw reads are stored in the main folder of ONT (e.g /var/lib/minknow/data). Enter the experiment name as input. 
# - Sequences are basecalled based on the model of choice. If enough computational power is available, we recommend "sup" method
# 
# ### 2. Demultiplexing 
# - Each reead is assigned to a well/plate combination. 
# 
# ### 3. Variant Calling
# - Minimap2 for creating Multiple Sequence Alignment (MSA)
# - Base Frequency Caller is being used for variant calling
# 
# 

# ### Packages 

# In[87]:


# Import all packages

import sys
sys.path.append("/home/emre/github_repo/MinION")

from minION.util import IO_processor
from minION.basecaller import Basecaller

from minION.variantcaller import *

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import importlib
importlib.reload(IO_processor)

import pickle
from Bio import SeqIO
import gzip
import subprocess
import mappy as mp
import holoviews as hv
import re

import ninetysix as ns
import colorcet as cc
import warnings

import bokeh.io
import holoviews as hv
from holoviews import opts

hv.extension('bokeh')
bokeh.io.output_notebook()


# In[88]:


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
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein


# In[89]:


def get_mut(temp_seq, aa_seq):
    i = 0
    mut_ls = []
    for i in range(len(list(zip(temp_seq, aa_seq)))):
        if temp_seq[i] != aa_seq[i]:
            mut_ls.append(temp_seq[i]+str(i+1)+aa_seq[i])
        i = i + 1
    return mut_ls


# In[90]:


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


# In[91]:


#### Heatmap ####
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


# In[92]:


def save_platemap_to_file(heatmaps, outputdir, name):
    if not os.path.exists(os.path.join(outputdir, "Platemaps")):
        os.makedirs(os.path.join(outputdir, "Platemaps"))
    file_path = os.path.join(outputdir, "Platemaps", name)
    hv.renderer('bokeh').save(heatmaps, file_path)

def save_csv(df,outputdir,name):
    if not os.path.exists(os.path.join(outputdir, "Results")):
        os.makedirs(os.path.join(outputdir, "Results"))
    file_path = os.path.join(outputdir, "Results", name + ".csv")
    df.to_csv(file_path)


# ### Meta Data 
# 
# - Provide the following arguments:
# 
# - Result Path: Path where the minion result folder will be created. All experiment results are then stored within the folder
# - Experiment Name: The experiment name is assigned when running the sequencer. Use the same name for identification
# 

# In[93]:


result_path = Path("/home/longy/")
experiment_name = "20240208-JR"
basecall_model_type = "sup"
result_folder = IO_processor.create_folder( experiment_name,
                                            basecall_model_type, 
                                            target_path=result_path)




# Create Barcode fasta file 
barcode_path = "../minION/barcoding/minion_barcodes.fasta" # Path to standard barcode file
front_prefix = "NB"
back_prefix = "RB"
bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
barcode_path = result_folder / "minion_barcodes_filtered.fasta"

# Barcode indexes
front_min = 1
front_max = 96
back_min = 1
back_max = 12

# Expected fragment sizes
min_size = 800
max_size = 5000

bp.filter_barcodes(barcode_path, (front_min,front_max), (back_min,back_max))


file_to_experiment= f"/var/lib/minknow/data/{experiment_name}"
template_fasta = "/home/longy/jr.fasta"

# Basecalling
basecall_folder = result_folder / "basecalled"
basecall_folder.mkdir(parents=True, exist_ok=True)
experiment_folder = IO_processor.find_experiment_folder(experiment_name) # Folder where pod5 files are located

# Demultiplexing
experiment_name = experiment_name + "_" + basecall_model_type
result_folder_path = IO_processor.find_folder(result_path, experiment_name)


# In[94]:


# Add conditions to avoid running the script accidentally
skip_basecalling = True
skip_demultiplex = False
skip_variant_calling = False


# ### Step 1 (Optional): Basecall reads
# 
# - Basecall can usually be done while sequencing (if GPU available?)
# - Otherwise, basecall afterwards

# In[95]:


if not skip_basecalling:


    pod5_files = IO_processor.find_folder(experiment_folder, "pod5")
    bc = Basecaller(basecall_model_type, pod5_files, basecall_folder, fastq = True)
    bc.run_basecaller()


# In[96]:


# Find fastq files
file_to_fastq = IO_processor.find_folder(experiment_folder, "fastq_pass")
print(file_to_fastq)


# ### Step 2: Demultiplex with SW
# - Demultiplex with SW 

# In[97]:


if not skip_demultiplex:
    path_to_code = "/home/emre/github_repo/MinION/source/source/demultiplex"
    prompt = f"{path_to_code} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100} -m {min_size} -x {max_size}"
    subprocess.run(prompt, shell=True)


# In[99]:


demultiplex_folder = result_folder 
print(demultiplex_folder)


# ### Step 3: Call Variant with PileUP Analysis
# 
# - Call Variant with min freq of 0.4 & depth min 15

# Read Summary file (Optional):
# 

# In[100]:


demultiplex_folder_name = result_folder


# In[ ]:


if not skip_variant_calling:
    vc = VariantCaller(experiment_folder, 
                   template_fasta, 
                   demultiplex_folder_name=demultiplex_folder_name, 
                   padding_start=0, 
                   padding_end=0)
    
    variant_df = vc.get_variant_df(qualities=True, 
                                threshold=0.2,
                                min_depth=5)
    seq_gen = IO_processor.SequenceGenerator(variant_df, template_fasta)
    variant_df = seq_gen.get_sequences()
    #TODO: Save the variant_df to a file after running. Currently it is not saved.


# In[ ]:


#20 - 30
variant_df.to_csv(result_folder / "variant_df.csv", index=False)  


# In[11]:


variant_df.iloc[1,2]


# In[102]:


variant_df


# In[103]:


df_variants_ = variant_df.copy()


# In[104]:


for seq_record in SeqIO.parse(open(template_fasta),'fasta'):
    temp_seq = str(seq_record.seq).upper()
df_variants_.insert(0, 'template', temp_seq)
df_variants_


# In[107]:


df_variants_['Variant'].tolist()
df_variants_['Variant'].fillna('',inplace = True)


# In[109]:


# Loop through dataframe and replace mutations
mut_ls = []
for i in df_variants_.index:
    if isinstance(df_variants_['Variant'][i], np.ndarray):
        df_variants_['Variant'][i] =  df_variants_['Variant'][i].tolist()
        
    if df_variants_['Variant'][i] == '':
        mut_ls.append('NA')
    elif pd.isnull(df_variants_['Variant'][i]):
        mut_ls.append('NA')
    elif df_variants_['Variant'][i] == '#PARENT#':
        mut_ls.append(df_variants_['template'][i])
    elif 'DEL' in df_variants_['Variant'][i]:
        mut_ls.append('Deletion')
    
    else:
        val_new = [x[-1] for x in df_variants_['Variant'][i].split('_')]
        index = [int(s) for s in re.findall(r'\d+', df_variants_['Variant'][i])]
        index_bp = []
        var_seq = temp_seq
        for m in range(len(index)):
            index_bp.append(index[m]-1)
            var_seq = var_seq[:index_bp[m]] + val_new[m]+ var_seq[index_bp[m] + 1:]
        mut_ls.append(var_seq)


# In[110]:


# Translate mutated sequence to protein
aa_ls = []
for i in range(len(mut_ls)):
    if str(mut_ls[i]).upper() != 'NA':
        aa_ls.append(translate(str(mut_ls[i]).upper()))
    else:
        aa_ls.append('NAN')
df_variants_['Template Sequence'] = aa_ls


# In[111]:


# Compare to template sequence and get mutations
mut = []
temp_aa = translate(temp_seq)
for i in range(len(aa_ls)):
    mut.append('_'.join(get_mut(temp_aa, aa_ls[i])))
df_variants_['Mutations'] = mut
df_variants_['Alignment Probability'] = df_variants_['Probability'].fillna(0.0)
df_variants_['Alignment Count'] = df_variants_['Alignment_count'].fillna(0.0)
df_variants_


# In[112]:


# Fill in parents into mutations column
for i in df_variants_.index:
    if df_variants_['Alignment Probability'].iloc[i] == 0.0 and df_variants_['Mutations'].iloc[i] == '':
        df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace('', '#N.A.#')
    if df_variants_['Mutations'].iloc[i] == '':
        df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace('', '#PARENT#') 
df_variants_


# In[113]:


Well = df_variants_['Well'].tolist()
column = [Well[i].strip('ABCDEFGH') for Well[i] in Well]
row = [Well[i].rstrip('0123456789') for Well[i] in Well]


# In[114]:


df_variants_['Row'] = row
df_variants_['Column'] = column
df_variants_['Plate'] = df_variants_['Plate'].astype(str)
df_variants_.loc[df_variants_['Plate'] == '9', ['Plate']] = '09'
df_variants_


# In[115]:


hm_ = generate_platemaps(df_variants_)

hm_


# In[86]:


save_platemap_to_file(hm_, demultiplex_folder, experiment_name)
save_csv(df_variants_, demultiplex_folder, experiment_name)


# In[ ]:




