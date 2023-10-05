import os
from pathlib import Path
from minION.util.IO_processor import concatenate_fastq_files
from minION.util.globals import MEDAKA_MODELS
import subprocess

def process_fastq(fastq_folder, filename = "pre_consensus", prefix = "fastq_runid", delete = False):
    """Processes all runid_fastq files in a folder. It concatenates all reads in one file and returns the path of the file."""
    
    concatenate_fastq_files(fastq_folder, filename = filename, prefix = prefix, delete = delete)

    return os.path.join(fastq_folder, "pre_consensus.fastq")


def consensus_prompt(pre_consensus_file, output_dir, ref, n_threads = 4, model = "default"):
    """Function to get the medaka prompt"""

    model = MEDAKA_MODELS[model]

    prompt = f'medaka_consensus -i {pre_consensus_file} -t {n_threads} -m {model} -o {output_dir} -d {ref} -f'

    return prompt

def medeka_stitch_prompt(barcode_folder : Path, ref, output_name,qualities = True, n_threads = 4):
    """Function to get the medaka stitch prompt
    Input: - Path to Demultiplex folder
           - Path to reference sequence"""


    hdf_file = os.path.abspath(os.path.join(barcode_folder, "medaka", "consensus_probs.hdf"))

    # Check if the hdf file exists
    if not os.path.exists(hdf_file):
        raise Exception("The hdf file does not exist")

    prompt = f'medaka stitch {hdf_file} {ref} {output_name} --threads {n_threads} --quiet'

    if qualities:
        prompt += ' --qualities'

    return prompt

def run_medaka(prompt):
    """Function to run medaka"""
    return subprocess.run(prompt, shell=True)

