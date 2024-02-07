"""
Contain argument parsers used for command line interface and web interface
"""
# Import packages
import os
import tqdm
import argparse
from time import strftime

from pathlib import Path

# Import local packages
from minION.util import IO_processor
from minION.run_MinION import run_MinION

# Get the working directory
CWD = os.getcwd()

# Set default arguments 
padding_start = 0
padding_end = 0
min_depth = 5
threshold = 0.2

# Build the CLI argparser
def build_cli_parser():

    # Initialize
    parser = argparse.ArgumentParser()

    # Add required arguments
    required_args_group = parser.add_argument_group("Required Arguments","Arguments required for each run")
    required_args_group.add_argument("refseq", 
            help = "fasta file containing reference sequence information.",
            required = True)
    required_args_group.add_argument("folder",
            help = "Folder containing fastq.pass or pod5_pass files. Nanopore experiment saved location",
            required = True)

    # Add optional arguments
    optional_args_group = parser.add_argument_group("Optional Arguments", "Aditional arguments")
    optional_args_group.add_argument("--output",
            help = "Save location of result, default to current workign directory")
    optional_args_group.add_argument("--skip_basecalling",
            action = "store_true",
            help = "Skip the basecalling step, default is false")
    optional_args_group.addargument("--skip_demultiplexing",
            action = "store_true",
            help = "Skip the demultiplexing step, default is false")
    optional_args_group.add_argument("--skip_variantcalling",
            action = "store_true",
            help = "Skip the variant calling step, default is false")
    return parser

# Execute MinION
def execute_MinION():
    # Build parser
    parser = build_cli_parser()

    # Parse the arguments
    CL_ARGS = vars(parser.parse_args())
    
