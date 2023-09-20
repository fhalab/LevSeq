import sys
sys.path.append("/home/emre/github_repo/MinION")

import os 
import pandas as pd

from minION.util.IO_processor import get_barcode_dict
from minION import analyse



result_folder = "/home/emre/minION_results"

experiment_folder = "202309_errorprone_3_rtrimmed_f24"

template_fasta = "/home/emre/github_repo/MinION/minION/refseq/hetcpiii.fasta"
 
demultiplex_folder = os.path.join(result_folder, experiment_folder, "demultiplex")

barcode_dicts = get_barcode_dict(demultiplex_folder)

variant_df = analyse.get_variant_df(demultiplex_folder, template_fasta, sequences=True)

variant_df.to_csv("/home/emre/github_repo/MinION/examples/test.csv", index=False)