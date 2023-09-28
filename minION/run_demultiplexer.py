###
# This script is used to run the demultiplexer
# Author: Emre Guersoy
###

import os
import glob
import numpy as np
from pathlib import Path

# MinION modules
from minION.util.parser import create_parser, check_parser
from minION.util.IO_processor import create_folder, find_experiment_folder, find_folder, get_barcode_dict, concatenate_fastq_files
from minION.util.globals import BARCODES
from minION.demultiplexer import run_demultiplexer



def main():
    

    experiment_folder = Path("/var/lib/minknow/data/Masked_3RBC-Minion/1/20230926_1648_MN41105_FAX47927_8a712363/fastq_pass/")

    result_folder = Path("/home/emre/minION_results/MinION_RBC_0902723_sup/demultiplex")

    