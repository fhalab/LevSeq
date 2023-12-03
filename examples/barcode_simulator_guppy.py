
import sys
sys.path.append("/home/emre/github_repo/MinION")
from minION import analyser
from minION import consensus
from minION import demultiplexer
import importlib
from minION.util import IO_processor
from minION.util.globals import BARCODES, MEDAKA_MODELS, DEFAULT_TARGETS
importlib.reload(IO_processor)
importlib.reload(analyser)
importlib.reload(consensus)
importlib.reload(demultiplexer)
from pathlib import Path
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
#import plotly.express as px
import random
from Bio import SeqIO
import subprocess
import re
import shutil
import glob
from tqdm import tqdm
ref_seq = Path("/home/emre/github_repo/MinION/minION/refseq/hetcpiii_padded.fasta")
folders = glob.glob("data/min_read_depth/seq/*")

for var_path in tqdm(folders):
    var_namr = os.path.basename(var_path)
    depths = glob.glob(f"{var_path}/depth*")
    for depth in depths:
        
        folder_path = Path(depth)
        print("Processing", folder_path)
        consensus.get_consensus(folder_path, ref_seq, output_name = "consensus.fastq", qualities = True, consensus_folder = folder_path)