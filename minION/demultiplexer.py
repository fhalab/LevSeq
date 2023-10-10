import os
import glob
from pathlib import Path
from util.globals import BARCODES
from minION.util.IO_processor import concatenate_fastq_files
import subprocess
import concurrent.futures as futures


class Demultiplexer:
    """This is a custom demultiplexer adapted from guppy barcode. It is adapted for evSeq application and also uses Smith Waterman algorithm to find the best match for the barcode.
    
    """

    def __init__(self) -> None:
        pass
        
    def read_barcodes(self, barcode_file : str) -> dict:
        """Read the barcode file and return a dictionary of barcodes

        Read1 : 5' - Front Barcode ---------------------- Reverse Barcode Complement - 3'

        Read2 : 5' - Reverse Barcode ---------------------- Front Barcode Complement - 3'
        
        Args:
            - barcode_file: Path to the barcode file
        
        Return:
            - Dictionary of barcodes
        """
        pass
    
    def reverse_complement(self, seq):
        """ 
        Get the reverse complement of a sequence/barcode. Depending on the strand, reads will either have front or reverse complement barcode.
         
        """
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        complement = [complement_dict[base] for base in seq]

        return ''.join(complement)[::-1]



def get_prompt(result_folder, output_folder, data_path, barcode_kit, score = 60):
    """Get prompt for guppy_barcoder
    Args: 
        - Experiment folder, Fastq folder after bassecalling
        - output_folder, Output folder for barcoding
        - barcode_kit, Barcode kit used

    Return: 
        - Prompt for guppy_barcoder
    """

    prompt = f'guppy_barcoder --input_path {result_folder} --save_path {output_folder} --data_path {data_path} --barcode_kits {barcode_kit} --min_score_barcode_front {score}'

    return prompt


def run_prompt(prompt):
    """Run os prompts"""
    return subprocess.run(prompt, shell=True)


def run_demultiplexer(result_folder : Path, BARCODES : dict, fbc_score : int = 60, rbc_score : int = 50, output_folder = None, basecall_folder : Path = None) -> bool:
    """Create the prompt to run Guppy Barcoder. The function first checks if the basecalled folder exists. If not, it assumes that the fastq files are in the experiment folder.
    
    Args:
        - result_folder: Path to the experiment folder
        - BARCODES: Dictionary of barcodes
        - fbc_score: Score for the front barcode
        - rbc_score: Score for the reverse barcode
        - output_folder: Output folder for the demultiplexed files
    
    Return:
        - True if the function ran successfully
    """
    
    if output_folder is None:
        output_folder = os.path.join(result_folder, "demultiplex_60_un")
    
    else:
        output_folder = output_folder

    Path(output_folder).mkdir(parents=True, exist_ok=True)

    if basecall_folder is None:
        fastqfolder = os.path.join(result_folder, "basecalled_filtered_un")
        print("Fastqfolder:", fastqfolder)

    else:
        fastqfolder = basecall_folder
    
    #Check if the input path exists
    if not os.path.exists(fastqfolder):
        raise Exception("Basecalled folder does not exist. Please check if you have chosen the right experiment name.")
    
    data_path = os.path.join(os.path.dirname(__file__), "barcoding")
    print("Data_Path:", data_path)

    # Check if the output folder exists
    if not os.path.exists(output_folder):
        raise Exception("Demultiplex folder does not exist. Please check if you have chosen the right experiment name.")

    barcode_rbc = BARCODES["Barcode-kit-RBC-rev"]

    barcode_fbc = BARCODES["Barcode-kit-FBC-rev"]

    print("Barcode_rbc:", barcode_rbc, "Barcode_fbc:", barcode_fbc)

    rbc_prompt = get_prompt(fastqfolder, output_folder, data_path, barcode_rbc, score=rbc_score)

    print("RBC Prompt:", rbc_prompt)

    run_prompt(rbc_prompt) # Generate plate folders

    rbc_files = glob.glob(os.path.join(output_folder, "barcode*"))
    print(output_folder)
    print(rbc_files)
    # Check if the barcode folder exists
    if not rbc_files:
        raise Exception("Barcode folder does not exist. Either no barcodes were found or the barcode score is too high. Rerun the experiment and adapt the barcode score")


    # TODO: Use concurrent function to run the following code in parallel

    for rv_barcodes in rbc_files:
        concatenate_fastq_files(rv_barcodes, "demultiplexed", "fastq_runid", delete = True)
        front_prompt = get_prompt(rv_barcodes, rv_barcodes, data_path, barcode_fbc, score=fbc_score)
        run_prompt(front_prompt)

    return True

