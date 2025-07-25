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
from levseq.filter_orientation import filter_demultiplexed_folder
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
import matplotlib

import os
import logging
import pandas as pd
from pathlib import Path
import shutil
import subprocess
from Bio import SeqIO
import platform
import numpy as np
import tqdm

import os
import logging
import pandas as pd
from pathlib import Path
import shutil
import subprocess
from Bio import SeqIO
import platform
import numpy as np
import tqdm
import panel as pn
import holoviews as hv
from importlib import resources
from holoviews.streams import Tap

# Utility function to configure logging
def configure_logging(result_folder, cl_args):
    import sys
    from levseq import __version__
    
    # Define a more detailed log format with clean separation
    log_format = "%(asctime)s : %(levelname)s : %(message)s"
    
    # Create log handlers
    info_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_run.log"))
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(logging.Formatter(log_format))

    error_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_error.log"))
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(log_format))

    # Set up basic configuration with both handlers
    logging.basicConfig(level=logging.INFO, handlers=[info_handler, error_handler])
    
    # Log version information and command used to run
    command_used = " ".join(sys.argv)
    logging.info(f"LevSeq Version: {__version__}")
    logging.info(f"Command: {command_used}")
    
    # Log essential run parameters
    logging.info(f"Run name: {cl_args.get('name', 'Not specified')}")
    logging.info(f"Input path: {cl_args.get('path', 'Not specified')}")
    logging.info(f"Summary file: {cl_args.get('summary', 'Not specified')}")
    
    # Log optional parameters if specified
    if cl_args.get('output') and cl_args.get('output') != os.getcwd():
        logging.info(f"Output directory: {cl_args.get('output')}")
    if cl_args.get('oligopool'):
        logging.info("Running in oligopool mode")
    if cl_args.get('skip_demultiplexing'):
        logging.info("Skipping demultiplexing step")
    if cl_args.get('skip_variantcalling'):
        logging.info("Skipping variant calling step")
    if cl_args.get('threshold'):
        logging.info(f"Using variant threshold: {cl_args.get('threshold')}")

# Create result folder
def create_result_folder(cl_args):
    folder_name = cl_args.get("name")
    if not folder_name:
        raise ValueError("The 'name' key is required in cl_args")
    output_path = cl_args.get("output", os.getcwd())
    result_folder = Path(output_path) / folder_name
    result_folder.mkdir(parents=True, exist_ok=True)
    return str(result_folder)

# Split FASTQ file into chunks
def split_fastq_file(fastq_file: Path, output_dir: Path, reads_per_file: int):
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
                    if chunk:
                        output_file = output_dir / f"{fastq_file.stem}_part{file_count}.fastq"
                        with open(output_file, "wt") as out_handle:
                            SeqIO.write(chunk, out_handle, "fastq")
                        logging.info(f"Created {output_file} with {len(chunk)} reads")
                    break
                output_file = output_dir / f"{fastq_file.stem}_part{file_count}.fastq"
                with open(output_file, "wt") as out_handle:
                    SeqIO.write(chunk, out_handle, "fastq")
                logging.info(f"Created {output_file} with {len(chunk)} reads")
                file_count += 1
        logging.info(f"Splitting complete for {fastq_file}. {file_count} parts created.")
    except Exception as e:
        logging.error(f"Failed to split FASTQ file {fastq_file}: {str(e)}", exc_info=True)
        raise

# Concatenate or split FASTQ files
def cat_fastq_files(folder_path: str, output_path: str, reads_per_file: int = 4000):
    try:
        folder_path = Path(folder_path)
        output_path = Path(output_path)
        if not folder_path.is_dir():
            raise ValueError("The provided path is not a valid directory: %s" % folder_path)
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        fastq_files = []
        for root, dirs, files in os.walk(folder_path):
            if "fastq_fail" not in root:
                for file in files:
                    if file.endswith(('.fastq', '.fastq.gz')):
                        fastq_files.append(Path(root) / file)
        if not fastq_files:
            raise ValueError("No FASTQ files found in %s" % folder_path)
        if len(fastq_files) == 1 and fastq_files[0].suffix == '.fastq':
            logging.info("Splitting single FASTQ file into %d reads per file", reads_per_file)
            split_fastq_file(fastq_files[0], output_path, reads_per_file)
        else:
            for fastq_file in fastq_files:
                destination = output_path / fastq_file.name
                # Skip copying if source and destination are identical
                if str(fastq_file) == str(destination):
                    logging.info("Skipping copy of %s (source and destination are identical)", fastq_file)
                    continue
                try:
                    shutil.copy(fastq_file, destination)
                    logging.info("Copied %s to %s", fastq_file, destination)
                except shutil.SameFileError:
                    logging.info("Skipping copy of %s (source and destination are identical files)", fastq_file)
        logging.info("All FASTQ files processed successfully to %s", output_path)
        return str(output_path)
    except Exception as e:
        logging.error("Failed to copy or split fastq files. An error occurred: %s", str(e), exc_info=True)
        raise

# Filter barcodes
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

def filter_bc(cl_args: dict, name_folder: Path, i: int) -> Path:
    front_min, front_max, rbc = barcode_user(cl_args, i)
    try:
        with resources.path('levseq.barcoding', 'minion_barcodes.fasta') as barcode_path:
            barcode_path = Path(barcode_path)
    except ImportError:
        package_root = Path(__file__).resolve().parent.parent
        barcode_path = package_root / "levseq" / "barcoding" / "minion_barcodes.fasta"
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

# Filter barcodes from input fasta
def filter_barcodes(input_fasta, output_fasta, barcode_range, rbc, front_prefix, back_prefix):
    front_min, front_max = barcode_range
    filtered_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if (
            record.id.startswith(front_prefix)
            and front_min <= int(record.id[len(front_prefix):]) <= front_max
        ) or (
            record.id.startswith(back_prefix)
            and int(record.id[len(back_prefix):]) == rbc
        ):
            filtered_records.append(record)
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fasta")

# Demultiplex fastq reads
def demux_fastq(file_to_fastq, result_folder, barcode_path):
    system_architecture = platform.machine().lower()
    if system_architecture == 'arm64':
        executable_name = "demultiplex-arm64"
    elif system_architecture == 'aarch64':
        executable_name = "demultiplex"
    elif system_architecture == 'x86_64':
        executable_name = "demultiplex-x86"
    else:
        raise ValueError(f"Unsupported architecture: {system_architecture}")
    try:
        with resources.path('levseq.barcoding', executable_name) as executable_path:
            executable_path = Path(executable_path)
    except ImportError:
        package_root = Path(__file__).resolve().parent.parent
        executable_path = package_root / "levseq" / "barcoding" / executable_name
    if not executable_path.exists():
        raise FileNotFoundError(f"Executable not found: {executable_path}")
    seq_min = 200
    seq_max = 10000
    prompt = f"{executable_path} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w 100 -r 100 -m {seq_min} -x {seq_max}"
    subprocess.run(prompt, shell=True, check=True)

# Variant calling using VariantCaller class

def call_variant(experiment_name, experiment_folder, template_fasta, filtered_barcodes, threshold=0.5, oligopool=False):
    try:
        vc = VariantCaller(
            experiment_name,
            experiment_folder,
            template_fasta,
            filtered_barcodes,
            padding_start=0,
            padding_end=0,
            oligopool=oligopool
        )
        variant_df = vc.get_variant_df(threshold=threshold, min_depth=5)
        logging.info("Variant calling to create consensus reads successful")
        return variant_df
    except Exception as e:
        logging.error("Variant calling failed", exc_info=True)
        raise

def assign_alignment_probability(row):
    if row["Variant"] == "#PARENT#":
        if row["Alignment Count"] > 20:
            return 1
        elif 10 <= row["Alignment Count"] <= 20:
            return (row["Alignment Count"] - 10) / 10  # Ranges from 0 to 1 linearly
        else:
            return 0
    else:
        return row["Average mutation frequency"]


# Full version of create_df_v function
def create_df_v(variants_df):
    # Make copy of dataframe
    df_variants_ = variants_df.copy()

    # Fill in empty cells
    df_variants_["Variant"] = df_variants_["Variant"].replace(np.nan, "", regex=True)

    # Create nc_variant column
    df_variants_["nc_variant"] = df_variants_.apply(
        lambda row: create_nc_variant(row["Variant"], row["refseq"]), axis=1
    )

    # Translate nc_variant to aa_variant
    df_variants_["aa_variant"] = df_variants_["nc_variant"].apply(
        lambda x: x if x in ["Deletion", "#N.A.#", 'Insertion'] else translate(x)
    )
    # Fill in 'Deletion' in 'aa_variant' column
    df_variants_.loc[
        df_variants_["nc_variant"] == "#DEL#", "aa_variant"
    ] = "#DEL#"
    df_variants_.loc[
        df_variants_["nc_variant"] == "#INS#", "aa_variant"
    ] = "#INS#"

    # Compare aa_variant with translated refseq and generate Substitutions column
    df_variants_["Substitutions"] = df_variants_.apply(get_mutations, axis=1)
    # Adding sequence quality to Alignment Probability before filling in empty values
    df_variants_["Alignment Probability"] = df_variants_.apply(assign_alignment_probability, axis=1)
    df_variants_["Alignment Probability"] = df_variants_["Alignment Probability"].fillna(0.0)
    df_variants_["Alignment Count"] = df_variants_["Alignment Count"].fillna(0.0)

    # Fill in Deletion into Substitutions Column, keep #N.A.# unchanged
    for i in df_variants_.index:
        if df_variants_["nc_variant"].iloc[i] == "Deletion":
            df_variants_.Substitutions.iat[i] = df_variants_.Substitutions.iat[i].replace("", "#DEL#")
        elif df_variants_["nc_variant"].iloc[i] == "#N.A.#":
            df_variants_.Substitutions.iat[i] = "#N.A.#"

    # Low read counts override low mutations
    df_variants_["Substitutions"] = ["#LOW#" if a < 10 and a > 0 else s for a, s in df_variants_[["Alignment Count", "Substitutions"]].values]

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

    # First rename columns as before
    df_variants_.rename(columns={
        "Variant": "nucleotide_mutation",
        "Substitutions": "amino_acid_substitutions",
        "nc_variant": "nt_sequence",
        "aa_variant": "aa_sequence"
        }, inplace=True)

    # Create a copy for restructuring to avoid affecting the original
    restructured_df = df_variants_.copy()
    restructured_df.columns = restructured_df.columns.str.lower().str.replace(r'[\s-]', '_', regex=True)
    # Fix the specific column name
    restructured_df.columns = restructured_df.columns.str.replace('p_adj._value', 'p_adj_value')

    # Select the desired columns in the desired order
    restructured_df = restructured_df[[
        "barcode_plate",
        "plate",
        "well",
        "alignment_count",
        "nucleotide_mutation",
        "amino_acid_substitutions",
        "alignment_probability",
        "average_mutation_frequency",
        "p_value",
        "p_adj_value",
        "nt_sequence",
        "aa_sequence"
    ]]

    return restructured_df, df_variants_

# Helper functions for create_df_v
def create_nc_variant(variant, refseq):
    if isinstance(variant, np.ndarray):
        variant = variant.tolist()
    if variant == "" or pd.isnull(variant):
        return "#N.A.#"  # Return #N.A.# if variant is empty or null
    elif variant == "#PARENT#":
        return refseq
    elif "DEL" in variant:
        return "#DEL#"
    elif "INS" in variant:
        return "#INS#"
    else:
        mutations = variant.split("_")
        nc_variant = list(refseq)
        for mutation in mutations:
            try:
                position = int(re.findall(r"\d+", mutation)[0]) - 1
                original = mutation[0]
                new = mutation[-1]
                if position < len(nc_variant) and nc_variant[position] == original:
                    nc_variant[position] = new
            except:
                print('WARNING! UNABLE TO PROCESS THIS')
                print(mutation)
        return "".join(nc_variant)


def is_valid_dna_sequence(sequence):
    return all(nucleotide in 'ATGC' for nucleotide in sequence) and len(sequence) % 3 == 0

def get_mutations(row):
    try:
        alignment_count = row["Alignment Count"]
        
        # Check if alignment_count is zero and return "#N.A.#" if true
        if alignment_count == 0:
            return "#N.A.#"
        if alignment_count <= 10:
            return "#LOW#"

        refseq = row["refseq"]

        if not is_valid_dna_sequence(refseq):
            return "Invalid refseq provided, check template sequence. Only A, T, G, C and sequence dividable by 3 are accepted."

        refseq_aa = translate(refseq)
        variant_aa = row["aa_variant"]

        if variant_aa == "Deletion":
            return ""
        else:
            mutations = []
            if len(refseq_aa) == len(variant_aa):
                for i in range(len(refseq_aa)):
                    if refseq_aa[i] != variant_aa[i]:
                        mutations.append(f"{refseq_aa[i]}{i+1}{variant_aa[i]}")
                if not mutations:
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


# Save plate maps and CSV
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

# Function to process the reference CSV and generate variants
def process_ref_csv_oligopool(cl_args, tqdm_fn=tqdm.tqdm):
    ref_df = pd.read_csv(cl_args["summary"])
    result_folder = create_result_folder(cl_args)
    variant_csv_path = os.path.join(result_folder, "variants.csv")
    variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "variant"])
    
    # First get the different barcode plates (these will be unique)
    barcode_plates = ref_df["barcode_plate"].unique()
    ref_df["barcode_index"] = [i for i in range(len(ref_df))]
    barcode_to_index = dict(zip(ref_df.barcode_plate, ref_df.barcode_index))
    for barcode_plate in barcode_plates:
        if not cl_args["skip_demultiplexing"]:
            i = barcode_to_index[barcode_plate]
            name_folder = os.path.join(result_folder, f'RB{barcode_plate}')
            os.makedirs(name_folder, exist_ok=True)
            barcode_path = filter_bc(cl_args, name_folder, i)
            output_dir = Path(result_folder) / f"{cl_args['name']}_fastq"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            file_to_fastq = cat_fastq_files(cl_args.get("path"), output_dir)
            try:
                demux_fastq(output_dir, name_folder, barcode_path)
            except Exception as e:
                logging.error("An error occurred during demultiplexing for sample {}. Skipping this sample.".format(barcode_plate), exc_info=True)
                continue
            # Check this - need to see if the code works... ToDo: Ariane
    # Now they are all demultiplexed, we can call variants
    if not cl_args["skip_variantcalling"]:
        for i, row in tqdm_fn(ref_df.iterrows(), total=len(ref_df), desc="Processing Samples"):
            barcode_plate = row["barcode_plate"]
            name = row["name"]
            refseq = row["refseq"].upper()
            # Get the name folder and barcode path
            temp_fasta_path = os.path.join(result_folder, f"temp_{name}.fasta")
            if not os.path.exists(temp_fasta_path):
                with open(temp_fasta_path, "w") as f:
                    f.write(f">{name}\n{refseq}\n")
            else:
                logging.info(f"Fasta file for {name} already exists. Skipping write.")
            try:
                filtered_barcodes = filter_bc(cl_args, result_folder, i)
                variant_result = call_variant(f"{name}", result_folder, temp_fasta_path, filtered_barcodes,
                                              oligopool=True)
                variant_result["barcode_plate"] = barcode_plate
                variant_result["name"] = name
                variant_result["refseq"] = refseq
                variant_df = pd.concat([variant_df, variant_result])
            except Exception as e:
                logging.error("An error occurred during variant calling for sample {}. Skipping this sample.".format(name), exc_info=True)
                continue
        
        variant_df.to_csv(variant_csv_path, index=False)
    # visualize it as well
    return variant_df, ref_df


# Function to process the reference CSV and generate variants
def process_ref_csv(cl_args, tqdm_fn=tqdm.tqdm):
    ref_df = pd.read_csv(cl_args["summary"])
    result_folder = create_result_folder(cl_args)
    variant_csv_path = os.path.join(result_folder, "variants.csv")

    variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "variant"])
    
    for i, row in tqdm_fn(ref_df.iterrows(), total=len(ref_df), desc="Processing Samples"):
        barcode_plate = row["barcode_plate"]
        name = row["name"]
        refseq = row["refseq"].upper()

        name_folder = os.path.join(result_folder, name)
        os.makedirs(name_folder, exist_ok=True)
        
        temp_fasta_path = os.path.join(name_folder, f"temp_{name}.fasta")
        if not os.path.exists(temp_fasta_path):
            with open(temp_fasta_path, "w") as f:
                f.write(f">{name}\n{refseq}\n")
        else:
            logging.info(f"Fasta file for {name} already exists. Skipping write.")
        
        barcode_path = filter_bc(cl_args, name_folder, i)
        output_dir = Path(result_folder) / f"{cl_args['name']}_fastq"
        output_dir.mkdir(parents=True, exist_ok=True)

        if not cl_args["skip_demultiplexing"]:
            file_to_fastq = cat_fastq_files(cl_args.get("path"), output_dir)
            try:
                demux_fastq(output_dir, name_folder, barcode_path)

                # Add filtering step here with multithreading
                filtered_counts = filter_demultiplexed_folder(
                        name_folder, 
                        refseq,
                        num_threads=10
                )
                logging.info(f"Orientation filtering completed for {name}")
                total_reads = sum(counts['total'] for counts in filtered_counts.values())
                kept_reads = sum(counts['kept'] for counts in filtered_counts.values())
                logging.info(f"Total filtering results: {kept_reads}/{total_reads} reads kept ({kept_reads/total_reads*100:.2f}%)")
                for file, counts in filtered_counts.items():
                    logging.info(f"{file}: {counts['kept']}/{counts['total']} reads kept")


            except Exception as e:
                logging.error("An error occurred during demultiplexing/filtering for sample {}. Skipping this sample.".format(name), exc_info=True)
                continue
        
        if not cl_args["skip_variantcalling"]:
            try:
                threshold = cl_args.get("threshold") if cl_args.get("threshold") is not None else 0.5
                variant_result = call_variant(
                    f"{name}", name_folder, temp_fasta_path, barcode_path, threshold=threshold
                )
                variant_result["barcode_plate"] = barcode_plate
                variant_result["name"] = name
                variant_result["refseq"] = refseq

                variant_df = pd.concat([variant_df, variant_result])
            except Exception as e:
                logging.error("An error occurred during variant calling for sample {}. Skipping this sample.".format(name), exc_info=True)
                continue
    
    variant_df.to_csv(variant_csv_path, index=False)

    return variant_df, ref_df


# Main function to run LevSeq and ensure saving of intermediate results if an error occurs
def run_LevSeq(cl_args, tqdm_fn=tqdm.tqdm):
    result_folder = create_result_folder(cl_args)
    # Ref folder for saving ref csv file
    ref_folder = os.path.join(result_folder, "ref")
    os.makedirs(ref_folder, exist_ok=True)
    
    configure_logging(result_folder, cl_args)
    logging.info("Logging configured. Starting analysis...")

    variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "variant"])

    try:
        if cl_args["oligopool"]:
            variant_df, ref_df = process_ref_csv_oligopool(cl_args, tqdm_fn)
        else:
            variant_df, ref_df = process_ref_csv(cl_args, tqdm_fn)
        ref_df_path = os.path.join(ref_folder, cl_args["name"]+".csv")
        ref_df.to_csv(ref_df_path, index=False)

        if variant_df.empty:
            logging.warning("No data found during CSV processing. The CSV is empty.")
    except Exception as e:
        variant_csv_path = os.path.join(result_folder, "variants_partial.csv")
        variant_df.to_csv(variant_csv_path, index=False)
        logging.error("An error occurred during processing summary file. Partial results saved at {}".format(variant_csv_path), exc_info=True)
        raise
    
    try:
        variant_csv_path = os.path.join(result_folder, "variants.csv")
        if os.path.exists(variant_csv_path):
            variant_df = pd.read_csv(variant_csv_path)
        
        if variant_df.empty:
            raise ValueError("The variant DataFrame is empty after processing. Unable to continue.")
        
        df_variants, df_vis = create_df_v(variant_df)
        processed_csv = os.path.join(result_folder, "visualization_partial.csv")
        df_vis.to_csv(processed_csv, index=False)
        if cl_args["oligopool"]:
            make_oligopool_plates(df_vis, result_folder=result_folder, save_files=True)
    except Exception as e:
        processed_csv = os.path.join(result_folder, "visualization_partial.csv")
        if 'df_vis' in locals():
            df_vis.to_csv(processed_csv, index=False)
        logging.error("An error occurred while preparing data for visualization. Partial visualization saved at {}".format(processed_csv), exc_info=True)
        raise
    
    try:
        layout = generate_platemaps(
            max_combo_data=df_vis,
            result_folder=result_folder,
            show_msa=cl_args["show_msa"],
        )
        save_platemap_to_file(
            heatmaps=layout,
            outputdir=result_folder,
            name=cl_args["name"],
            show_msa=cl_args["show_msa"],
        )
        save_csv(df_variants, result_folder, cl_args["name"])
        logging.info("Run successful, see visualization and results")
    except Exception as e:
        partial_csv_path = os.path.join(result_folder, "variants_partial.csv")
        df_variants.to_csv(partial_csv_path, index=False)
        logging.error("An error occurred during visualization. Partial CSV saved at {}".format(partial_csv_path), exc_info=True)
        raise

# This modification saves the results at each critical stage, ensuring that even in the case of failure,
# the user has access to intermediate results and does not lose all the progress.

