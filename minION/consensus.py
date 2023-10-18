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

def run_medaka(prompt):
    """Function to run medaka"""
    return subprocess.run(prompt, shell=True)


def create_consensus_folder(folder_path : Path, folder_name = "consensus"):
    """Function to create the consensus folder"""

    consensus_folder = folder_path / folder_name
    Path(consensus_folder).mkdir(parents=True, exist_ok=True)



def mini_align(folder_path : Path, ref : Path, n_threads = 1, prefix = "" ,output_name = "alignment", consensus_folder = "consensus"):
    """Mini align uses minimap2 and samtools to align the reads to a reference sequence
    Args:
        - folder_path (Path): Path to the folder containing the fastq files
        - ref (Path): Path to the reference sequence
        - n_threads (int, optional): Number of threads to use. Defaults to 1.
        - output_name (str, optional): Name of the output file. Defaults to "alignment.bam". If you want to change
        it, makes sure to add extension and refernce it for down stream tasks.
    
    Returns:
        - subprocess.CompletedProcess: Completed process from subprocess.run
    
    """
        
    prompt = f'mini_align -r {ref} -i {folder_path}/{prefix}*.fastq -t {n_threads} -m -p alignment && mv *.bam *.bam.bai {folder_path}/{consensus_folder}' 

    return subprocess.run(prompt, shell=True)


def medaka_consensus(folder_path : Path, consensus_folder = "consensus"):
    """Run Medaka consensus on aligned bam files"""

    # Check if consensus folder exists
    consensus_path = folder_path / consensus_folder

    # Check if the hdf file exists
    if not os.path.exists(consensus_path):
        raise Exception("Consensus folder does not exist")
    
    prompt = f"medaka consensus {consensus_path}/*.bam {consensus_path}/pre_consensus.hdf --batch 200 --threads 4"

    return subprocess.run(prompt, shell=True)

def medaka_stitch(folder_path : Path, ref : Path, output_name = "consensus.fastq", qualities = True, consensus_folder = "consensus"):
    """Run Medaka stitch on the hdf file"""

    # Check if consensus folder exists
    consensus_path = folder_path / consensus_folder

    # Check if the hdf file exists
    if not os.path.exists(consensus_path):
        raise Exception("Consensus folder does not exist")
    
    prompt = f"medaka stitch {consensus_path}/pre_consensus.hdf {ref} {consensus_path}/{output_name} --threads 4"

    if qualities:
        prompt += " --qualities"

    return subprocess.run(prompt, shell=True)

def medaka_variant(folder_path : Path, ref : Path, output_name = "pre_consensus.hdf"): #TODO, not functional yet
    """Run Medaka variant on the hdf file"""

    # Check if consensus folder exists
    consensus_path = folder_path / "consensus"

    # Check if the hdf file exists
    if not os.path.exists(consensus_path):
        raise Exception("Consensus folder does not exist")
    
    prompt = f"medaka variant {consensus_path}/pre_consensus.hdf {ref} {output_name} --threads 4"

    return subprocess.run(prompt, shell=True)

def get_consensus(folder_path : Path, ref : Path, output_name = "consensus.fastq", qualities = True, consensus_folder = "consensus"):
    """Function to get the consensus sequence from the fastq files"""

    create_consensus_folder(folder_path, consensus_folder)

    # Mini align
    mini_align(folder_path, ref, 4, "concat", "alignment.bam", consensus_folder)

    # Medaka consensus
    medaka_consensus(folder_path, consensus_folder)

    # Medaka stitch
    medaka_stitch(folder_path, ref, output_name, qualities, consensus_folder)
