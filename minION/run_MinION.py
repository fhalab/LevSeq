# Import MinION objects
from minION.util import IO_processor
from minION.basecaller import Basecaller
from minION.variantcaller import *

# Import external packages
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import sys
import importlib
import tqdm

# Get barcode used
def barcode_user(cl_arg):
    # Set some default values if user did not provide barcodes
    fmin = 1
    fmax = 96
    bmin = 1
    bmax = 12
    bc_df = pd.read_csv(cl_arg['barcodes'])
    fmin = bc_df['NB-min'][0]
    fmax = bc_df['NB-max'][0]
    bmin = bc_df['RB-min'][0]
    bmax = bc_df['RB-max'][0]

    return fmin,fmax,bmin,bmax

# Get output directory
def result_path(path):
    return path

# Get fastq input directory, this is the basecalled folder
def fastq_path(path):
    return path

# Filter barcode
def filter_bc(cl_args):
    front_min,front_max,back_min,back_max = barcode_user(cl_args)
    barcode_path = "minION/barcoding/minion_barcodes.fasta" # Path to standard barcode file
    front_prefix = "NB"
    back_prefix = "RB"
    bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
    bp.filter_barcodes(barcode_path, (front_min,front_max), (back_min,back_max))

# Filter template sequence length
def filter_seq(cl_args)
# Get reference fasta (parent sequence)
def parent_user(csv):
    return fasta

# Demultiplex the basecalled fastq into plate-well folders
def demux_fastq(result_path, ref_csv):
    path_to_code = "source/source/demultiplex"
    prompt = f"{path_to_code} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100}"
    subprocess.run(prompt, shell=True)

# Run MinION
def run_MinION(cl_args, tqdm_fn = tqdm.tqdm):
    
    
