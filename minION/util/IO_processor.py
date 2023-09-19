# This script contains function to check for input data and process it for downstream analysis.

import os
import glob
from minION.util.globals import MINKNOW_PATH
from Bio import SeqIO
from pathlib import Path


def check_data_folder(path):
    """
    Check if minknown data folder exists. If not, assert an error.
    """
    if not os.path.exists(path):
        raise Exception("MinKNOW data folder does not exist. Please check if you have installed MinKNOW and if the path is correct.")

def find_experiment_folder(experiment_name, minknow_path = MINKNOW_PATH):
    """Find the experiment folder in the minknow data folder. Access the name from Meta File."""

    check_data_folder(minknow_path)

    for (path, dirs, files) in os.walk(minknow_path, topdown=True):
        if experiment_name in dirs:
            return os.path.join(path, experiment_name)
    raise Exception("Experiment folder does not exist. Please check if you have chosen the right experiment name.")


def find_folder(start_path, target_folder_name):
    """
    Find a specific folder within the start_path. Access the name from Meta File.
    
    Parameters:
    - start_path (str): The path to begin the search from.
    - target_folder_name (str): The name of the folder to search for.

    Returns:
    - str: The full path to the found folder.

    Raises:
    - Exception: If the target folder is not found.
    """
    for (path, dirs, _) in os.walk(start_path, topdown=True):
        if target_folder_name in dirs:
            return os.path.join(path, target_folder_name)
    raise Exception(f"{target_folder_name} folder does not exist. Please check if you have chosen the right name.")


def find_file(start_path, prefix, extension):
    """
    Find a file in the given start_path that matches the specified prefix and extension.
    
    Parameters:
    - start_path (str): The directory to begin the search from.
    - prefix (str): The prefix of the filename to search for.
    - extension (str): The extension of the file to search for (e.g., '.txt').

    Returns:
    - str: The full path to the found file.
    
    Raises:
    - Exception: If no matching file is found.
    """
    
    for dirpath, _, filenames in os.walk(start_path):
        for filename in filenames:
            if filename.startswith(prefix) and filename.endswith(extension):
                return os.path.join(dirpath, filename)
    
    raise Exception(f"No file with prefix '{prefix}' and extension '{extension}' was not found in {start_path}.")



def concat_all_fastq(path: Path, filename : str = None, prefix : str = None, delete = True):
    """
    Concatenate all fastq files in a folder
    Input: path, where the fastq files are located
    Output: path of the concatenated fastq file
    """
    files = glob.glob(f"{path}/*.fastq")
    n_files = len(files)
    
    if n_files == 1:
        single_file = files[0]
        if os.path.basename(single_file) != f"{filename}.fastq":
            os.rename(single_file, os.path.join(path, f"{filename}.fastq"))

    
    elif n_files == 0:
        return "No fastq files found"

    else:
        if prefix is not None:
            os.system(f"cat {path}/{prefix}*.fastq > {path}/{filename}.fastq")
        
            if delete:
                os.system(f"rm {path}/{prefix}*.fastq")
        
        else:
            os.system(f"cat {path}/*.fastq > {path}/{filename}.fastq")

            if delete:
                os.system(f"rm {path}/*.fastq")


        
    return "Concatenation successful"


def read_fasta_file(path, score = False):
    """
    Read fasta file from an input path
    Input:  - path, where the fasta file is located
            - score, if True, return sequence and quality scores
    Output: Dictionary with sequence only or sequence and quality scores
    """

    if score:
        sequences_and_scores = {"Sequence" : [], "Quality-Score" : []}

        for record in SeqIO.parse(path, "fastq"):
            sequences_and_scores["Sequence"].append(str(record.seq))
            sequences_and_scores["Quality-Score"].append(record.letter_annotations["phred_quality"])

        return sequences_and_scores

    else:
        sequences = {"Sequence" : []}

        file_extension = os.path.splitext(path)[1][1:] # Get file extension
        
        for record in SeqIO.parse(path, file_extension):
            sequences["Sequence"].append(str(record.seq))

        return sequences



def create_folder(experiment_name: str, target_path: Path = None, output_name: str = None) -> Path:
    """
    When Starting minION, the function checks if a minION result folder exists. If not, it creates one. 
    It also checks for the subfolder of the experiment.
    If no information is given about the folder, it creates a default folder in the current directory.
    
    Input:  - target_path, path of the folder to be created
    Output: - Path object representing the experiment folder
            - Raises Exception if the folder could not be created or path is invalid
    """
    
    if target_path is None:
        # Get current working directory
        curr_dir = Path.cwd()
    else:
        if not target_path.exists():
            raise Exception("Target path does not exist. Please check if you have chosen the right path.")
        curr_dir = target_path
    
    # Create minION_results folder if it doesn't exist
    minION_results_dir = curr_dir / "minION_results"
    minION_results_dir.mkdir(exist_ok=True)
    
    # Create experiment folder

    if output_name is None:
        result_folder = minION_results_dir / experiment_name
    else:
        result_folder = minION_results_dir / output_name

    result_folder.mkdir(exist_ok=True)

    return result_folder

    
def get_rbc_barcode_folders(demultiplex_folder):
    """Get the barcode folder from the demultiplex folder
    Input:  - result_folder, where the demultiplex folder is located
    Output: - barcode_folder, where the barcode folders are located"""

    if not os.path.exists(demultiplex_folder):
        raise Exception("Demultiplex folder does not exist. Run minION to get the demultiplex folder.")
    
    reverse_barcodes = glob.glob(os.path.join(demultiplex_folder, "barcode*"))
    print(f"Reverse Barcdoes: {reverse_barcodes}")
    return reverse_barcodes

def get_fbc_barcode_folders(demultiplex_folder, rbc_barcode_folders):
    """Get the barcode folder from the demultiplex folder
    Input:  - result_folder, where the demultiplex folder is located
            - rbc_barcode_folders, name of reverse barcode folders
    Output: - barcode_folder, where the barcode folders are located"""

    fbc_folders = {}

    if rbc_barcode_folders is None:
        raise Exception("Reverse barcode folders are not given. Please check if you have chosen the right path.")
    
    for folder in rbc_barcode_folders:
        fbc_path = os.path.join(folder, "barcode*")
        fbc_folders[folder] = glob.glob(fbc_path)

    print(fbc_folders)
    return fbc_folders

def get_barcode_dict(demultiplex_folder):
    """Get a dictionary of folder paths, where reverse is the key and forward is stored in a list"""

    rbc_folders = get_rbc_barcode_folders(demultiplex_folder)
    fbc_folders = get_fbc_barcode_folders(demultiplex_folder, rbc_folders)

    return fbc_folders

def trim_fasta(input_file, output_file, trim_length=12):
    with open(input_file, 'r') as in_fasta, open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            trimmed_seq = record.seq[trim_length:]
            record.seq = trimmed_seq
            SeqIO.write(record, out_fasta, 'fasta')

