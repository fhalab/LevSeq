import os
import glob
from pathlib import Path
from util.globals import BARCODES
from minION.util.IO_processor import concat_all_fastq


def get_prompt(result_folder, output_folder, data_path, barcode_kit, score = 60):
    """Get prompt for guppy_barcoder
    Input: 
    Experiment folder, Fastq folder after bassecalling
    output_folder, Output folder for barcoding
    barcode_kit, Barcode kit used

    Output: Prompt for guppy_barcoder"""

    prompt = f'guppy_barcoder --input_path {result_folder} --save_path {output_folder} --data_path {data_path} --barcode_kits {barcode_kit} --min_score_barcode_front {score}'

    return prompt


def run_prompt(prompt):
    """Run os prompts"""
    return os.system(prompt)


def run_demultiplexer(result_folder, BARCODES, fbc_score = 60, rbc_score = 50):
    """Run Demultiplexer """
    
    # Create a demultiplex folder if not exists
    output_folder = os.path.join(result_folder, "demultiplex")
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    fastqfolder = os.path.join(result_folder, "basecalled")

    #Check if the input path exists
    if not os.path.exists(fastqfolder):
        raise Exception("Basecalled folder does not exist. Please check if you have chosen the right experiment name.")
    
    data_path = os.path.join(os.path.dirname(__file__), "barcoding")
    print("Data_Path:", data_path)

    # Check if the output folder exists
    if not os.path.exists(output_folder):
        raise Exception("Demultiplex folder does not exist. Please check if you have chosen the right experiment name.")

    barcode_rbc = BARCODES["Barcode-kit-RBC"]

    barcode_fbc = BARCODES["Barcode-kit"]

    rbc_prompt = get_prompt(fastqfolder, output_folder, data_path, barcode_rbc, score=rbc_score)

    run_prompt(rbc_prompt) # Generate plate folders

    rbc_files = glob.glob(os.path.join(output_folder, "barcode*"))

    # Check if the barcode folder exists
    if not rbc_files:
        raise Exception("Barcode folder does not exist. Either no barcodes were found or the barcode score is too high. Rerun the experiment and adapt the barcode score")

    for rv_barcodes in rbc_files:
        concat_all_fastq(rv_barcodes, "demultiplexed", "fastq_runid", delete = True)
        front_prompt = get_prompt(rv_barcodes, rv_barcodes, data_path, barcode_fbc, score=fbc_score)
        run_prompt(front_prompt)

    return True

