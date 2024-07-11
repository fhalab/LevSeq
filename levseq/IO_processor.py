###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

# This script contains function to check for input data and process it for downstream analysis.

import os
import glob
from Bio import SeqIO
from pathlib import Path
import gzip
import subprocess
import re
import numpy as np
import pandas as pd


class SequenceGenerator:
    """
    A class for processing Sequence data. The class uses variant data frame to generate sequences.
    """

    def __init__(self, variant_df, reference_path, padding_start=50, padding_end=50):
        self.variant_df = variant_df
        self.reference = get_template_sequence(reference_path)
        self.padding_start = padding_start
        self.padding_end = padding_end

    def get_sequences(self):
        self.variant_df["Sequence"] = self.variant_df.apply(self.get_sequence, axis=1)
        return self.variant_df

    def get_sequence(self, row):
        variant = row["Variant"]
        if pd.isnull(variant):
            return np.nan
        else:
            # Process the variant and generate the sequence
            seq = self.generate_sequence(str(variant))
            return seq

    def generate_sequence(self, variant: str, ref_only=True):
        """
        Generate sequence from variant string. E.g "A123T
        """

        new_seq = self.reference

        if isinstance(variant, float):
            return float('nan')

        elif variant != "#PARENT#":

            variants = variant.split("_")

            for var in variants:
                match = re.match(r'[A-Za-z]+(\d+)([A-Za-z]+)', var)
                if match:
                    position = int(match.group(1))
                    adj_pos = position - 1 + self.padding_start
                    new_nucleotide = match.group(2)
                else:
                    raise ValueError(f"Invalid variant format: {var}")
                if new_nucleotide == "DEL":
                    new_seq = new_seq[:adj_pos] + "-" + new_seq[adj_pos + 1:]

                else:
                    new_seq = new_seq[:adj_pos] + new_nucleotide + new_seq[adj_pos + 1:]

        # Remove "-" from sequence
        new_seq = new_seq.replace("-", "")

        if ref_only:
            return new_seq[self.padding_start: -self.padding_end]

        return new_seq


class BarcodeProcessor:
    """
    A class for processing barcode fasta file. 
    """

    def __init__(self, barcode_path, front_prefix="NB", reverse_prefix="RB"):
        self.barcode_path = barcode_path
        self.front_prefix = front_prefix
        self.reverse_prefix = reverse_prefix

    def filter_barcodes(self, filtered_fasta_file, front_range, reverse_number):
        """
        Filters barcodes in the given ranges and writes them to a new fasta file.

        Args:
        - filtered_fasta_file (str): The path to the filtered fasta file.
        - front_range (tuple): A tuple of two integers representing the range for front barcodes.
        - reverse_number (int): An integer representing the specific reverse barcode number.
        """
        records = list(SeqIO.parse(self.barcode_path, "fasta"))
        min_front, max_front = front_range

        filtered_records = [
            record for record in records
            if (record.id.startswith(self.front_prefix) and min_front <= int(
                record.id[len(self.front_prefix):]) <= max_front) or
               (record.id.startswith(self.reverse_prefix) and int(
                   record.id[len(self.reverse_prefix):]) == reverse_number)
        ]

        with open(filtered_fasta_file, "w") as output_handle:
            SeqIO.write(filtered_records, output_handle, "fasta")

    def get_barcode_dict(self, demultiplex_folder: Path) -> dict:
        """
        Get a dictionary of folder paths, where reverse is the key and forward is stored in a list
        
        Args:
        - demultiplex_folder, where the demultiplex folder is located
        
        Returns:
        - barcode_dict, dictionary of reverse and front barcode paths
        """

        rbc_folders = get_rbc_barcode_folders(demultiplex_folder, prefix=self.reverse_prefix)
        fbc_folders = get_fbc_barcode_folders(rbc_folders, prefix=self.front_prefix)

        return fbc_folders

    def get_rbc_barcode_folders(self, demultiplex_folder: Path) -> list:
        """
        Extract the reverse barcode folders (rbc) from the demultiplex folder
        Args:
            - demultiplex_folder (Path): Where the demultiplex folder is located.
        Returns:
            - reverse_barcodes (list): Where the barcode folders are located. 
        """

        if not demultiplex_folder.exists():
            raise FileNotFoundError(
                f"Demultiplex folder '{demultiplex_folder}' does not exist. Run levseq to get the demultiplex folder.")

        reverse_barcodes = list(demultiplex_folder.glob(f"{self.reverse_prefix}*"))

        if not reverse_barcodes:
            # Optionally, use logging here instead of raising an exception
            raise FileNotFoundError(
                f"No reverse barcodes found in {demultiplex_folder}. Either no barcodes were found or the barcode score is too high. Rerun the experiment or adapt the barcode score in the TOML file.")

        return reverse_barcodes


def find_folder(start_path, target_folder_name):
    """
    Find a specific folder within the start_path. Access the name from Meta File.
    
    Args:
    - start_path (str): The path to begin the search from.
    - target_folder_name (str): The name of the folder to search for.

    Returns:
    - str: The full path to the found folder.
    - Exception: If the target folder is not found.
    """
    for (path, dirs, _) in os.walk(start_path, topdown=True):
        if target_folder_name in dirs:
            return os.path.join(path, target_folder_name)
    raise Exception(f"{target_folder_name} folder does not exist. Please check if you have chosen the right name.")


def check_data_folder(path: Path) -> None:
    """
    Check if minknown data folder exists. If not, assert an error.

    Args:
    - path (str): The path to the minknow data folder.

    Returns:
    - FileNotFoundError: If the minknow data folder does not exist.
    """
    if not Path(path).exists():
        raise FileNotFoundError(
            f"MinKNOW data folder '{path}' does not exist. Please check if the path is correct.")
    else:
        return path

def find_experiment_folder(experiment_name: str, minknow_path) -> None:
    """Find the experiment folder in the minknow data folder. Access the name from Meta File.
    
    Args:
    - experiment_name (str): The name of the experiment to search for.
    - minknow_path (str): The path to the minknow data folder.

    Returns:
    - str: The full path to the found experiment folder.
    """

    check_data_folder(Path(minknow_path))

    for (path, dirs, files) in os.walk(minknow_path, topdown=True):
        if experiment_name in dirs:
            return os.path.join(path, experiment_name)

    raise Exception("Experiment folder does not exist. Please check if you have chosen the right experiment name.")

def find_folder(start_path, target_folder_name):
    """
    Find a specific folder within the start_path. Access the name from Meta File.
    
    Args:
    - start_path (str): The path to begin the search from.
    - target_folder_name (str): The name of the folder to search for.

    Returns:
    - str: The full path to the found folder.
    - Exception: If the target folder is not found.
    """
    for (path, dirs, _) in os.walk(start_path, topdown=True):
        if target_folder_name in dirs:
            return os.path.join(path, target_folder_name)
    raise Exception(f"{target_folder_name} folder does not exist. Please check if you have chosen the right name.")


def find_experiment_files(start_path: Path, target_folder_name: list) -> Path or None:
    """
    Check if the experimenter has already basecalled before. If pod5_pass or fastq_pass folder exists, the function returns True.

    Args:
        - start_path (Path): The path to begin the search from.
        - target_folder_name (str): The name of the folder to search for.

    Returns:
        - Target folder (Path): The full path to the found folder.
        - False: If the target folder is not found.

    """

    for value in target_folder_name:
        try:
            return find_folder(start_path, value)
        except:
            continue

    return None


def extract_files_from_folder(path: Path) -> list:
    """Extract all files from a folder"""
    return [file for file in path.iterdir() if file.is_file()]


def filter_fastq_by_length(input_fastq: Path, output_fastq: Path, min_length: int, max_length: int):
    """
    Filter a FASTQ file based on read length and write to a new FASTQ file.

    Args:
        - input_fastq: Path to the input FASTQ gzip file.
        - output_fastq: Path to the output FASTQ gzip file.
        - min_length: Minimum length of reads to retain.
        - max_length: Maximum length of reads to retain.
    Returns:
        - Fastq file with reads between min_length and max_length.
        - N_reads
    """
    with gzip.open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:

        while True:

            header = infile.readline().strip()
            sequence = infile.readline().strip()
            plus_sign = infile.readline().strip()
            quality_scores = infile.readline().strip()

            if not header:
                break
            N_reads = 0
            if min_length <= len(sequence) <= max_length:
                outfile.write(f"{header}\n")
                outfile.write(f"{sequence}\n")
                outfile.write(f"{plus_sign}\n")
                outfile.write(f"{quality_scores}\n")
                N_reads += 1

    return N_reads


def find_file(start_path: Path, prefix: str, extension: str) -> Path:
    """
    Find a file in the given start_path that matches the specified prefix and extension.
    
    Args:
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


def concatenate_fastq_files(path: Path, filename: str = "concatenated", prefix: str = "*", delete: bool = True) -> str:
    """
    Concatenate all fastq files in a directory.
    
    Args:
        - path (Path): Directory where fastq files are located.
        - filename (str): Name of the concatenated file. Default is 'concatenated'.
        - prefix (str): Prefix of files to be concatenated. Default is '*' (all files).
        - delete (bool): Whether to delete the original files after concatenation.
    
    Returns:
        - str: Message indicating the result of the operation.
    """
    search_pattern = path / f"{prefix}.fastq"
    fastq_files = [Path(p) for p in glob.glob(str(search_pattern))]

    if not fastq_files:
        return "No fastq files found."

    if len(fastq_files) == 1:
        single_file = fastq_files[0]

        if single_file.name != f"{filename}.fastq":
            target_path = path / f"{filename}.fastq"

            if target_path.exists():
                return f"Target file {target_path} already exists. Rename operation aborted."

            single_file.rename(target_path)
            return f"Single fastq file in {path} found and renamed."

    else:
        # Concatenate multiple files
        target_path = path / f"{filename}.fastq"

        if target_path.exists():
            return f"Target file {target_path} already exists. Concatenation aborted."

        with target_path.open("w") as outfile:
            for fastq_file in fastq_files:
                with fastq_file.open() as infile:
                    outfile.write(infile.read())

        if delete:
            for fastq_file in fastq_files:
                fastq_file.unlink()

        return "Concatenation successful."


def read_fasta_file(path: Path, score=False) -> dict:
    """
    Read fasta file from an input path

    Args:  
        - path, where the fasta file is located
        - score, if True, return sequence and quality scores

    Returns: 
        - Dictionary with sequence only or sequence and quality scores
    """

    if score:
        sequences_and_scores = {"Sequence": [], "Quality-Score": []}

        for record in SeqIO.parse(path, "fastq"):
            sequences_and_scores["Sequence"].append(str(record.seq))
            sequences_and_scores["Quality-Score"].append(record.letter_annotations["phred_quality"])

        return sequences_and_scores

    else:
        sequences = {"Sequence": []}

        file_extension = os.path.splitext(path)[1][1:]  # Get file extension

        for record in SeqIO.parse(path, "fasta"):
            sequences["Sequence"].append(str(record.seq))

        return sequences


def create_folder(experiment_name: str, target_path: Path = None, output_name: str = None) -> Path:
    """
    When Starting levseq, the function checks if a levseq result folder exists. If not, it creates one.
    It also checks for the subfolder of the experiment.
    If no information is given about the folder, it creates a default folder in the current directory.
    
    Args:

    - target_path, path of the folder to be created

    Returns: 

    - Path object representing the experiment folder
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
    minION_results_dir = curr_dir / "EVSeqL_Output"
    minION_results_dir.mkdir(exist_ok=True)

    # Create experiment folder

    experiment_name = f"{experiment_name}"

    if output_name is None:
        result_folder = minION_results_dir / experiment_name
    else:
        output_name = f"{output_name}"
        result_folder = minION_results_dir / output_name

    result_folder.mkdir(exist_ok=True)

    return result_folder


def get_rbc_barcode_folders(demultiplex_folder: Path, prefix="barcode") -> list:
    """Extract the reverse barcode folders (rbc) from the demultiplex folder
    Args:
    - demultiplex_folder (Path): Where the demultiplex folder is located.
    Returns:
    - reverse_barcodes (list): Where the barcode folders are located. 
    """

    if not demultiplex_folder.exists():
        raise FileNotFoundError(
            f"Demultiplex folder '{demultiplex_folder}' does not exist. Run levseq to get the demultiplex folder.")

    reverse_barcodes = list(demultiplex_folder.glob(f"{prefix}*"))

    if not reverse_barcodes:
        # Optionally, use logging here instead of raising an exception
        raise FileNotFoundError(
            f"No reverse barcodes found in {demultiplex_folder}. Either no barcodes were found or the barcode score is too high. Rerun the experiment or adapt the barcode score in the TOML file.")

    return reverse_barcodes


def get_fbc_barcode_folders(rbc_barcode_folders: list, prefix="barcode") -> dict:
    """Extract the forward barcode folders (fbc) within the reverse barcode folders

    Args:  
    - result_folder, where the demultiplex folder is located
    - rbc_barcode_folders, name of reverse barcode folders

    Returns: 
    - Barcode folders (dict), where the forward barcode folders are located    
    """

    fbc_folders = {}

    if rbc_barcode_folders is None:
        raise Exception("Reverse barcode folders are not given. Please check if you have chosen the right path.")

    for folder in rbc_barcode_folders:
        fbc_folders[folder] = list(folder.glob(f"{prefix}*"))

    if not fbc_folders:
        raise Exception(
            f"Forward barcodes in {rbc_barcode_folders} do not exist. Either no barcodes were found or the barcode score is too high. Rerun the experiment or adapt the barcode score in the TOML file. ")

    if not any(fbc_folders.values()):
        raise Exception(f"Forward barcodes in {rbc_barcode_folders} do not exist. Please check the folders.")

    return fbc_folders


def get_barcode_dict(demultiplex_folder: Path, front_prefix="barcode", reverse_prefix="barcode") -> dict:
    """
    Get a dictionary of folder paths, where reverse is the key and forward is stored in a list
    
    Args:
    - demultiplex_folder, where the demultiplex folder is located
    
    Returns:
    - barcode_dict, dictionary of reverse and front barcode paths
    """

    rbc_folders = get_rbc_barcode_folders(demultiplex_folder, prefix=reverse_prefix)
    fbc_folders = get_fbc_barcode_folders(rbc_folders, prefix=front_prefix)

    return fbc_folders


def trim_fasta(input_file, output_file, trim_length=12):
    with open(input_file, 'r') as in_fasta, open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            trimmed_seq = record.seq[trim_length:]
            record.seq = trimmed_seq
            SeqIO.write(record, out_fasta, 'fasta')


def filter_fastq_by_length(result_folder, input_fastq: Path, min_length: int, max_length: int, ind: int):
    """
    Filter a FASTQ file based on read length and write to a new FASTQ file. Currently we use fastq-filter.

    Args:
        - input_fastq: Path to the input FASTQ file or folder.
        - output_fastq: Path to the output FASTQ file.
        - min_length: Minimum length of reads to retain.
        - max_length: Maximum length of reads to retain.
    Returns:
        - Fastq file with reads between min_length and max_length.
        - N_reads
    """

    if input_fastq.is_dir():
        input_files = f'{input_fastq}/*.fastq'

    elif input_fastq.is_file():
        input_files = input_fastq

    else:
        raise Exception("Input file is not a file or a folder. Please check if you have chosen the right path.")

    basecall_folder = os.path.join(result_folder, "basecalled_filtered")
    Path(basecall_folder).mkdir(parents=True, exist_ok=True)

    prompt = f'fastq-filter -o {basecall_folder}/basecalled_filtered{ind}.fastq.gz -l {min_length} -L {max_length} {input_files} --quiet'

    subprocess.run(prompt, shell=True)


def get_template_sequence(path: Path) -> str:
    """
    Read template sequence fasta file
        Args:  
            - path, where the fasta file is located
        Returns: 
            - Template sequence
    """

    template = read_fasta_file(path)

    return template["Sequence"][0]
