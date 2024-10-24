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

# Import MinION objects
from levseq import *

# Import external packages
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from importlib import resources
import subprocess
from Bio import SeqIO
import tqdm
import platform
import subprocess
import os
import re
import gzip
import shutil

import panel as pn
import holoviews as hv
from holoviews.streams import Tap

output_notebook()

pn.extension()
pn.config.comms = "vscode"

hv.extension("bokeh")


# Get barcode used
def barcode_user(cl_args, i):
    try:
        # Set some default values if user did not provide barcodes
        fmin = 1
        fmax = 96
        bc_df = pd.read_csv(cl_args["summary"])
        rbc = bc_df["barcode_plate"][i]
        logging.info(f"Demultiplex executed successfully for index {i}.")

        return int(fmin), int(fmax), int(rbc)

    except Exception as e:
        logging.error("Demultiplex failed to execute for index {i}.", exc_info=True)
        raise

# Split fastq file into 4000 reads per chunk
def split_fastq_file(fastq_file: Path, output_dir: Path, reads_per_file: int):
    """
    Splits a FASTQ file into multiple files, each containing up to a specified number of reads.
    
    Parameters:
    - fastq_file (Path): The input FASTQ file to be split.
    - output_dir (Path): The directory where the split files will be saved.
    - reads_per_file (int): The number of reads per output file.
    """
    try:
        with open(fastq_file, "rt") as handle:
            record_iter = SeqIO.parse(handle, "fastq")
            file_count = 0
            while True:
                chunk = []
                try:
                    for _ in range(reads_per_file):
                        chunk.append(next(record_iter))
                except StopIteration:
                    # End of the file reached
                    if chunk:
                        output_file = output_dir / f"{fastq_file.stem}_part{file_count}.fastq"
                        with open(output_file, "wt") as out_handle:
                            SeqIO.write(chunk, out_handle, "fastq")
                        logging.info(f"Created {output_file} with {len(chunk)} reads")
                    break  # Exit the loop once we reach the end of the records
                output_file = output_dir / f"{fastq_file.stem}_part{file_count}.fastq"
                with open(output_file, "wt") as out_handle:
                    SeqIO.write(chunk, out_handle, "fastq")
                logging.info(f"Created {output_file} with {len(chunk)} reads")
                file_count += 1
        
        logging.info(f"Splitting complete for {fastq_file}. {file_count} parts created.")
    except Exception as e:
        logging.error(f"Failed to split FASTQ file {fastq_file}: {str(e)}", exc_info=True)
        raise

def cat_fastq_files(folder_path: str, output_path: str, reads_per_file: int = 4000):
    """
    Copies all .fastq and .fastq.gz files from the provided folder_path to the output_path.
    If there's only one .fastq file, it will be split into smaller files of a specified number
    of reads each.
    
    Parameters:
    - folder_path (str): The path to the directory containing .fastq and .fastq.gz files.
    - output_path (str): The path to the directory where the files should be copied or split.
    - reads_per_file (int): The number of reads per output file when splitting a FASTQ file.
    
    Returns:
    - str: The path to the output directory where the files were copied.
    
    Raises:
    - ValueError: If the folder_path is not a valid directory or if no .fastq/.fastq.gz files are found.
    - Exception: For any other errors that occur during file copying or splitting.
    """
    try:
        folder_path = Path(folder_path)
        output_path = Path(output_path)

        if not folder_path.is_dir():
            raise ValueError(f"The provided path {folder_path} is not a valid directory")

        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)

        # Find all fastq and fastq.gz files excluding those in fastq_fail folders
        fastq_files = []
        for root, dirs, files in os.walk(folder_path):
            if "fastq_fail" not in root:
                for file in files:
                    if file.endswith((".fastq", ".fastq.gz")):
                        fastq_files.append(Path(root) / file)

        if not fastq_files:
            raise ValueError(f"No FASTQ files found in {folder_path}")

        # If there is only one FASTQ file, split it into smaller files
        if len(fastq_files) == 1 and fastq_files[0].suffix == '.fastq':
            logging.info(f"Splitting single FASTQ file into {reads_per_file} reads per file")
            split_fastq_file(fastq_files[0], output_path, reads_per_file)
        else:
            # Copy each .fastq or .fastq.gz file to the output directory
            for fastq_file in fastq_files:
                destination = output_path / fastq_file.name
                shutil.copy(fastq_file, destination)
                logging.info(f"Copied {fastq_file} to {destination}")

        logging.info(f"All FASTQ files processed successfully to {output_path}")
        return str(output_path)

    except Exception as e:
        logging.error(f"Failed to copy or split fastq files. An error occurred: {str(e)}", exc_info=True)
        raise

# Create result folder
def create_result_folder(cl_args: dict) -> str:
    folder_name = cl_args.get("name")
    if not folder_name:
        raise ValueError("The 'name' key is required in cl_args")
    output_path = cl_args.get("output", os.getcwd())
    result_folder = Path(output_path) / folder_name
    # Create the directory if it doesn't exist
    result_folder.mkdir(parents=True, exist_ok=True)
    return str(result_folder)


# Return and create filtered barcodes
def filter_bc(cl_args: dict, name_folder: Path, i: int) -> Path:
    front_min, front_max, rbc = barcode_user(cl_args, i)
# Use importlib.resources to get the path to the barcode file
    try:
        with resources.path('levseq.barcoding', 'minion_barcodes.fasta') as barcode_path:
            barcode_path = Path(barcode_path)
    except ImportError:
        # Fallback method if the above fails
        package_root = Path(__file__).resolve().parent.parent
        barcode_path = package_root / "levseq" / "barcoding" / "minion_barcodes.fasta"
    # Ensure the barcode file exists
    if not barcode_path.exists():
        raise FileNotFoundError(f"Barcode file not found: {barcode_path}")

    front_prefix = "NB"
    back_prefix = "RB"
    barcode_path_filter = os.path.join(name_folder, "levseq_barcodes_filtered.fasta")
    
    filter_barcodes(
        str(barcode_path),
        str(barcode_path_filter),
        (front_min, front_max),
        rbc,
        front_prefix,
        back_prefix,
    )
    
    return barcode_path_filter

# Filter barcodes
def filter_barcodes(
    input_fasta, output_fasta, barcode_range, rbc, front_prefix, back_prefix
):
    front_min, front_max = barcode_range
    filtered_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        if (
            record.id.startswith(front_prefix)
            and front_min <= int(record.id[len(front_prefix) :]) <= front_max
        ) or (
            record.id.startswith(back_prefix)
            and int(record.id[len(back_prefix) :]) == rbc
        ):
            filtered_records.append(record)

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fasta")

# Demultiplex fastq reads into plate and wells format
def demux_fastq(file_to_fastq, result_folder, barcode_path):
    # Determine the system architecture
    system_architecture = platform.machine().lower()

    # Choose the appropriate executable based on the architecture
    if system_architecture == 'arm64':
        executable_name = "demultiplex-arm64"
    elif system_architecture == 'aarch64':
        executable_name = "demultiplex"
    elif system_architecture == 'x86_64':
        executable_name = "demultiplex-x86"
    else:
        raise ValueError(f"Unsupported architecture: {system_architecture}")

    # Use importlib.resources to get the path to the executable
    try:
        with resources.path('levseq.barcoding', executable_name) as executable_path:
            executable_path = Path(executable_path)
    except ImportError:
        # Fallback method if the above fails
        package_root = Path(__file__).resolve().parent.parent
        executable_path = package_root / "levseq" / "barcoding" / executable_name

    # Ensure the executable exists
    if not executable_path.exists():
        raise FileNotFoundError(f"Executable not found: {executable_path}")

    # Get min and max sequence length if user specified, otherwise use default
    seq_min = 150 
    seq_max = 10000

    # Use subprocess to run the executable
    prompt = f"{executable_path} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w 100 -r 100 -m {seq_min} -x {seq_max}"
    subprocess.run(prompt, shell=True, check=True)

# Variant calling using VariantCaller class and generate dataframe
def call_variant(experiment_name, experiment_folder, template_fasta, filtered_barcodes):
    try:
        vc = VariantCaller(
            experiment_name,
            experiment_folder,
            template_fasta,
            filtered_barcodes,
            padding_start=0,
            padding_end=0,
        )
        variant_df = vc.get_variant_df(threshold=0.5, min_depth=5)
        logging.info("Variant calling to create consensus reads successful")
        return variant_df
    except Exception as e:
        logging.error("Variant calling failed", exc_info=True)
        raise


# Saving heatmaps and csv in the results folder
def save_platemap_to_file(heatmaps, outputdir, name, show_msa):
    if not os.path.exists(os.path.join(outputdir, "Platemaps")):
        os.makedirs(os.path.join(outputdir, "Platemaps"))
    file_path = os.path.join(outputdir, "Platemaps", name)
    if show_msa:
        heatmaps.save(file_path + "_msa.html", embed=True)
    else:
        hv.renderer("bokeh").save(heatmaps, file_path)


def save_csv(df, outputdir, name):
    if not os.path.exists(os.path.join(outputdir, "Results")):
        os.makedirs(os.path.join(outputdir, "Results"))
    file_path = os.path.join(outputdir, "Results", name + ".csv")
    df.to_csv(file_path)


# Generate dataframe for visualization
def create_df_v(variants_df):
    # Make copy of dataframe
    df_variants_ = variants_df.copy()

    # Fill in empty cells
    df_variants_["Variant"].tolist()
    df_variants_["Variant"] = df_variants_["Variant"].replace(np.nan, "", regex=True)

    # Create nc_variant column
    df_variants_["nc_variant"] = df_variants_.apply(
        lambda row: create_nc_variant(row["Variant"], row["refseq"]), axis=1
    )

    # Translate nc_variant to aa_variant
    df_variants_["aa_variant"] = df_variants_["nc_variant"].apply(
        lambda x: "Deletion" if x == "Deletion" else translate(x)
    )
    # Fill in 'Deletion' in 'aa_variant' column
    df_variants_.loc[
        df_variants_["nc_variant"] == "Deletion", "aa_variant"
    ] = "Deletion"

    # Compare aa_variant with translated refseq and generate mutations column
    df_variants_["Mutations"] = df_variants_.apply(get_mutations, axis=1)

    # Fill in empty empty values
    df_variants_["Alignment Probability"] = df_variants_[
        "Average mutation frequency"
    ].fillna(0.0)
    df_variants_["Alignment Count"] = df_variants_["Alignment Count"].fillna(0.0)

    # Fill in Deletion into mutations Column
    for i in df_variants_.index:
        if df_variants_["nc_variant"].iloc[i] == "Deletion":
            df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace(
                "", "-"
            )

    # Add row and columns
    Well = df_variants_["Well"].tolist()
    row = []
    column = []
    for well in Well:
        if len(well) >= 2:
            row.append(well[0])
            if well[1:].isdigit():
                column.append(well[1:])
            else:
                column.append("")
        else:
            row.append("")
            column.append("")

    df_variants_["Row"] = row
    df_variants_["Column"] = column
    df_variants_["Plate"] = df_variants_["name"].astype(str)

    # Update 'Plate' column from '1'-'9' to '01'-'09'
    df_variants_["Plate"] = df_variants_["Plate"].apply(
        lambda x: f"0{x}" if len(x) == 1 else x
    )
    # Select the desired columns in the desired order
    restructured_df = df_variants_[
        [
            "barcode_plate",
            "Plate",
            "Well",
            "Variant",
            "Alignment Count",
            "Average mutation frequency",
            "P value",
            "P adj. value",
            "Mutations",
            "nc_variant",
            "aa_variant",
        ]
    ]
    # Set 'Mutations' and 'Variant' columns to '#N.A.#' if 'Alignment Count' is smaller than 5
    restructured_df.loc[
        restructured_df["Alignment Count"] < 6, ["Mutations", "Variant"]
    ] = "#N.A.#"
    df_variants_.loc[
        df_variants_["Alignment Count"] < 6, ["Mutations", "Variant"]
    ] = "#N.A.#"
    restructured_df.loc[
        restructured_df["Mutations"] == "#PARENT#", ["Alignment Probability"]
    ] = 1.0
    df_variants_.loc[
        df_variants_["Mutations"] == "#PARENT#", ["Alignment Probability"]
    ] = 1.0

    return restructured_df, df_variants_


def create_nc_variant(variant, refseq):
    if isinstance(variant, np.ndarray):
        variant = variant.tolist()
    if variant == "" or pd.isnull(variant):
        return refseq
    elif variant == "#PARENT#":
        return refseq
    elif "DEL" in variant:
        return "Deletion"
    else:
        mutations = variant.split("_")
        nc_variant = list(refseq)
        for mutation in mutations:
            if len(mutation) >= 2:
                position = int(re.findall(r"\d+", mutation)[0]) - 1
                original = mutation[0]
                new = mutation[-1]
            if position < len(nc_variant) and nc_variant[position] == original:
                nc_variant[position] = new
        return "".join(nc_variant)


def is_valid_dna_sequence(sequence):
    return all(nucleotide in 'ATGC' for nucleotide in sequence) and len(sequence) % 3 == 0

def get_mutations(row):
    try:
        refseq = row["refseq"]
        
        if not is_valid_dna_sequence(refseq):
            return "Invalid refseq provided, check template sequence. Only A, T, G, C and sequence dividable by 3 are accepted."

        refseq_aa = translate(refseq)
        variant_aa = row["aa_variant"]
        alignment_count = row["Alignment Count"]

        if variant_aa == "Deletion":
            return ""
        else:
            mutations = []
            if len(refseq_aa) == len(variant_aa):
                for i in range(len(refseq_aa)):
                    if refseq_aa[i] != variant_aa[i]:
                        mutations.append(f"{refseq_aa[i]}{i+1}{variant_aa[i]}")
                if not mutations:
                    if alignment_count < 5:
                        return "#N.A.#"
                    else:
                        return "#PARENT#"
            else:
                return "LEN"
        return "_".join(mutations) if mutations else ""

    except Exception as e:
        logging.error(
            "Translation to amino acids failed, check template sequence. Only A, T, G, C and sequence dividable by 3 are accepted.",
            exc_info=True,
        )   
        raise

# Process the summary file
def process_ref_csv(cl_args):
    ref_df = pd.read_csv(cl_args["summary"])

    result_folder = create_result_folder(cl_args)

    variant_csv_path = os.path.join(result_folder, "variants.csv")
    if os.path.exists(variant_csv_path):
        variant_df = pd.read_csv(variant_csv_path)
    else:
        variant_df = pd.DataFrame(
            columns=["barcode_plate", "name", "refseq", "variant"]
       )
    for i, row in ref_df.iterrows():
        barcode_plate = row["barcode_plate"]
        name = row["name"]
        refseq = row["refseq"].upper()

        # Create a subfolder for the current iteration using the name value
        name_folder = os.path.join(result_folder, name)
        os.makedirs(name_folder, exist_ok=True)

        # Write the refseq to a temporary fasta file
        temp_fasta_path = os.path.join(name_folder, f"temp_{name}.fasta")
        with open(temp_fasta_path, "w") as f:
            f.write(f">{name}\n{refseq}\n")
        # Create filtered barcode path
        barcode_path = filter_bc(cl_args, name_folder, i)
        # Find fastq.gz files
        output_dir = Path(result_folder) / "basecalled_reads"
        output_dir.mkdir(parents=True, exist_ok=True)

        file_to_fastq = cat_fastq_files(cl_args.get("path"), output_dir)

        if not cl_args["skip_demultiplexing"]:
            demux_fastq(output_dir, name_folder, barcode_path)
        if not cl_args["skip_variantcalling"]:
            variant_result = call_variant(
                f"{name}", name_folder, temp_fasta_path, barcode_path
            )
            variant_result["barcode_plate"] = barcode_plate
            variant_result["name"] = name
            variant_result["refseq"] = refseq

            variant_df = pd.concat([variant_df, variant_result])

        # Remove the temporary fasta file
        # os.remove(temp_fasta_path)
    variant_df.to_csv(variant_csv_path, index=False)
    return variant_df


# Run LevSeq
def run_LevSeq(cl_args, tqdm_fn=tqdm.tqdm):
    # Create output folder
    result_folder = create_result_folder(cl_args)

    # Configure logging to save in the output directory
    log_format = "%(asctime)s:%(levelname)s:%(message)s"

    # INFO level logger
    info_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_run.log"))
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(logging.Formatter(log_format))

    # ERROR level logger
    error_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_error.log"))
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(log_format))

    # Configure the root logger
    logging.basicConfig(level=logging.INFO, handlers=[info_handler, error_handler])
    try:
        # Process summary file by row using demux, call_variant function
        variant_df = process_ref_csv(cl_args)

        # Check if variants.csv already exist
        variant_csv_path = os.path.join(result_folder, "variants.csv")
        if os.path.exists(variant_csv_path):
            variant_df = pd.read_csv(variant_csv_path)
            df_variants, df_vis = create_df_v(variant_df)
        # Clean up and prepare dataframe for visualization
        else:
            df_variants, df_vis = create_df_v(variant_df)

        processed_csv = os.path.join(result_folder, "visualization.csv")
        df_vis.to_csv(processed_csv, index=False)

        layout = generate_platemaps(
            max_combo_data=df_vis,
            result_folder=result_folder,
            show_msa=cl_args["show_msa"],
        )

        # Saving heatmap and csv
        save_platemap_to_file(
            heatmaps=layout,
            outputdir=result_folder,
            name=cl_args["name"],
            show_msa=cl_args["show_msa"],
        )
        save_csv(df_variants, result_folder, cl_args["name"])
        logging.info("Run successful, see visualization and results")
    except Exception as e:
        logging.error(
            "An error occured while executing LevSeq, check log file for detail",
            exc_info=True,
        )
        raise
