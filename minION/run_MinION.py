# Import MinION objects
from MinION.minION.util import IO_processor
from MinION.minION.basecaller import Basecaller
from MinION.minION.variantcaller import *

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
    bc_df = pd.read_csv(cl_args["barcodes"])
    fmin = bc_df["NB-min"][0]
    fmax = bc_df["NB-max"][0]
    bmin = bc_df["RB-min"][0]
    bmax = bc_df["RB-max"][0]

    return int(fmin),int(fmax),int(bmin),int(bmax)

# Get output directory
def get_input_folder(cl_args):
    input_folder = IO_processor.find_experiment_folder(cl_args['name'], cl_args['folder'])
    return input_folder

# Get fastq input directory, this is the basecalled folder
def fastq_path(folder):
    return IO_processor.find_folder(folder, "fastq_pass")

# Create result folder
def create_result_folder(cl_args):
    basecall_model = 'sup'
    result_folder = IO_processor.create_folder(
            cl_args['name'],
            basecall_model,
            target_path = Path(cl_args['output']))
    return result_folder

# Basecall reads
def basecall_reads(cl_args):
    print('basecalling')

# Filter barcode
def filter_bc(cl_args, result_folder):
    front_min,front_max,back_min,back_max = barcode_user(cl_args)
    # Obtain path of executable from package
    with resources.path('MinION.minION.barcoding', 'minion_barcodes.fasta') as barcode_path:
        front_prefix = "NB"
        back_prefix = "RB"
        bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
    
    barcode_path_filter = os.path.join(result_folder, "minion_barcodes_filtered.fasta")
    bp.filter_barcodes(barcode_path_filter, (1, 96), (9, 12))
    return barcode_path

# Filter template sequence length
def filter_seq(cl_args):
    return seq_min, seq_max

# Get reference fasta (parent sequence)
def parent_fasta(cl_args):
    template_fasta = cl_args['refseq']
    return template_fasta

# Demultiplex the basecalled fastq into plate-well folders
def demux_fastq(file_to_fastq, result_folder, barcode_path):
    # Obtain path of executable from package
    with resources.path('MinION.source.source', 'demultiplex') as executable_path:
        # Get min and max sequence length if user specified, otherwise use default
        seq_min = 800
        seq_max = 5000
        # Use subprocess to run the executable
        prompt = f"{str(executable_path)} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100} -m {seq_min} -x {seq_max}"
        subprocess.run(prompt, shell=True)

# Variant calling using VariantCaller class and generate dataframe
def call_variant(experiment_folder, template_fasta, demultiplex_folder_name):
    vc = VariantCaller(experiment_folder, 
            template_fasta, 
            demultiplex_folder_name = demultiplex_folder_name,
            padding_start = 0,
            padding_end = 0)
    
    variant_df = vc.get_variant_df(qualities = True,
            threshold = 0.2,
            min_depth = 5)
    return variant_df

# Run MinION
def run_MinION(cl_args, tqdm_fn = tqdm.tqdm):
    # Find specific experiment in the upper directory of nanopore data
    experiment_folder = get_input_folder(cl_args)

    # Find fastq from experiment folder
    file_to_fastq = fastq_path(experiment_folder)

    # Create result folder
    result_folder = create_result_folder(cl_args)
    
    # Get template sequence
    template_fasta = parent_fasta(cl_args)

    # Basecall if asked
    if cl_args["perform_basecalling"]:
        basecall_reads(cl_args)

    # Filter barcodes and store new barcode file 
    barcode_path = filter_bc(cl_args, result_folder)

    # Demultiplex if not skipped
    if not cl_args['skip_demultiplexing']:
        demux_fastq(file_to_fastq, result_folder, barcode_path)

    # Call Variants if not skipped
    if not cl_args['skip_variantcalling']:
        variant_df = call_variant(experiment_folder, template_fasta, result_folder)
        variant_df.to_csv(result_folder / "variant_df.csv", index=False)  
