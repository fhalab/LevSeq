# Import MinION objects
from minION.util import IO_processor
from minION.basecaller import Basecaller
from minION.variantcaller import *

# Import external packages
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from importlib import resources
import subprocess
import sys
import importlib
import tqdm
import os

# Get barcode used
def barcode_user(cl_args):
    # Set some default values if user did not provide barcodes
    fmin = 1
    fmax = 96
    bmin = 1
    bmax = 12
    bc_df = pd.read_csv(cl_arg["barcodes"])
    fmin = bc_df["NB-min"][0]
    fmax = bc_df["NB-max"][0]
    bmin = bc_df["RB-min"][0]
    bmax = bc_df["RB-max"][0]

    return fmin,fmax,bmin,bmax

# Get output directory
def get_input_folder(cl_args):
    input_folder = IO_processor.find_experiment_folder(cl_args['name'], clargs['folder'])
    return input_folder

# Get fastq input directory, this is the basecalled folder
def fastq_path(cl_args):
    return IO_processor.find_folder(cl_args['folder'], "fastq_pass")

# Create result folder
def create_result_folder(cl_args):
    result_folder = IO_processor.create_folder(
            cl_args["name"],
            basecall_model,
            target_path = Path(cl_args["output"]))
    return result_folder
# Basecall reads
def basecall_reads(cl_args):

# Filter barcode
def filter_bc(cl_args, result_folder):
    front_min,front_max,back_min,back_max = barcode_user(cl_args)
    barcode_path = "minION/barcoding/minion_barcodes.fasta" # Path to standard barcode file
    front_prefix = "NB"
    back_prefix = "RB"
    bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
    barcode_path = os.join(result_folder, "minion_barcodes_filtered.fasta")
    bp.filter_barcodes(barcode_path, (front_min,front_max), (back_min,back_max))
    
    return barcode_path_filtered
# Filter template sequence length
def filter_seq(cl_args)

# Get reference fasta (parent sequence)
def parent_fasta(cl_args):
    template_fasta = clarg['refseq']
    return template_fasta

# Demultiplex the basecalled fastq into plate-well folders
def demux_fastq(file_to_fastq, result_folder, barcode_path):
    # Obtain path of executable from package
    with resources.path('MinION.source.source', 'demultiplex') as executable_path:
        # Use subprocess to run the executable
        prompt = f"{str(executable_path)} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100}"
        subprocess.run(prompt, shell=True)

# Run MinION
def run_MinION(cl_args, tqdm_fn = tqdm.tqdm):
    # Find input fastq
    file_to_fastq = fastq_path(cl_args)
    
    # Find data folder
    experiment_folder = get_input_folder(cl_args)

    # Create result folder
    result_folder = create_result_folder(cl_args)

    # Basecall if asked
    if cl_args["perform_basecalling"]:
        basecall_reads(cl_args)

    # Filter barcodes and store new barcode file 
    barcode_path = filter_bc(cl_args, result_folder)

    # Demultiplex if not skipped
    if not cl_args['skip_demultiplexing']:
        demux_fastq(file_to_fastq, result_folder, barcode_path)

    # Call Variants if not skipped
