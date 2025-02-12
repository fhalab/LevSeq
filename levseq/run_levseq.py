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

# Import LevSeq objects
from levseq import *
from levseq.oligo import OligoVariantCaller
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
def configure_logging(result_folder):
    log_format = "%(asctime)s:%(levelname)s:%(message)s"
    info_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_run.log"))
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(logging.Formatter(log_format))

    error_handler = logging.FileHandler(os.path.join(result_folder, "LevSeq_error.log"))
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(log_format))

    logging.basicConfig(level=logging.INFO, handlers=[info_handler, error_handler])

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
                shutil.copy(fastq_file, destination)
                logging.info("Copied %s to %s", fastq_file, destination)
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
    """
    Filter barcodes - no changes needed for oligo mode as it uses the same barcode filtering
    """
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
    
    return Path(barcode_path_filter)

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
def call_variant(experiment_name, experiment_folder, template_fasta, filtered_barcodes, ref_sequences_df=None, oligo_mode=False):
    """
    Call variants or find best matching sequences depending on mode
    
    Args:
        experiment_name (str): Name of the experiment
        experiment_folder (str): Path to experiment folder
        template_fasta (str): Path to template fasta file (not used in oligo mode)
        filtered_barcodes (str): Path to filtered barcodes file
        ref_sequences_df (pd.DataFrame, optional): DataFrame containing reference sequences for oligo mode
        oligo_mode (bool): Whether to run in oligo pool mode
    
    Returns:
        pd.DataFrame: DataFrame containing variant or match results
    """
    try:
        if oligo_mode and ref_sequences_df is not None:
            # Initialize OligoVariantCaller directly
            oligo_caller = OligoVariantCaller(
                experiment_name,
                experiment_folder,
                ref_sequences_df,
                filtered_barcodes
            )
            # Use the class method to get variants
            variant_df = oligo_caller.get_variant_df(min_depth=5)
            logging.info(f"Oligo sequence matching completed for {experiment_name}")
        else:
            vc = VariantCaller(
                experiment_name,
                experiment_folder,
                template_fasta,
                filtered_barcodes,
                padding_start=0,
                padding_end=0,
            )
            variant_df = vc.get_variant_df(threshold=0.5, min_depth=5)
            logging.info(f"Variant calling completed for {experiment_name}")
        
        return variant_df
    except Exception as e:
        logging.error(f"Error in variant calling for {experiment_name}: {str(e)}", exc_info=True)
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

def create_df_v(variants_df, oligo_mode=False):
    """
    Create visualization and results DataFrames
    
    Args:
        variants_df (pd.DataFrame): Input variants DataFrame
        oligo_mode (bool): Whether to process in oligo pool mode
    
    Returns:
        tuple: (results DataFrame, visualization DataFrame)
    """
    if oligo_mode:
        # For oligo mode, create a simplified DataFrame
        df_variants = variants_df.copy()
        
        # Create visualization DataFrame
        df_vis = pd.DataFrame({
            'Well': df_variants['Well'],
            'Plate': df_variants['barcode_plate'].astype(str),
            'Reference': df_variants['Reference'],
            'Reference_ID': df_variants['Reference_ID'],
            'Alignment Count': df_variants['Alignment Count'],
            'Match_Frequency': df_variants['Match_Frequency'],
            'Total_Reads': df_variants['Total_Reads']
        })
        
        # Extract row and column information from well IDs
        Well = df_vis['Well'].tolist()
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
        
        df_vis['Row'] = row
        df_vis['Column'] = column
        
        # Ensure plate numbers are properly formatted
        df_vis['Plate'] = df_vis['Plate'].apply(
            lambda x: f"0{x}" if len(str(x)) == 1 else str(x)
        )
        
        return df_vis, df_vis
     
    else:
        # Original variant processing logic
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
        
        # Adding sequence quality to Alignment Probability
        df_variants_["Alignment Probability"] = df_variants_.apply(assign_alignment_probability, axis=1)
        df_variants_["Alignment Probability"] = df_variants_["Alignment Probability"].fillna(0.0)
        df_variants_["Alignment Count"] = df_variants_["Alignment Count"].fillna(0.0)

        # Process special cases in Substitutions Column
        for i in df_variants_.index:
            if df_variants_["nc_variant"].iloc[i] == "Deletion":
                df_variants_.Substitutions.iat[i] = df_variants_.Substitutions.iat[i].replace("", "#DEL#")
            elif df_variants_["nc_variant"].iloc[i] == "#N.A.#":
                df_variants_.Substitutions.iat[i] = "#N.A.#"

        # Handle low read counts
        df_variants_["Substitutions"] = ["#LOW#" if a < 10 and a > 0 else s 
                                       for a, s in df_variants_[["Alignment Count", "Substitutions"]].values]

        # Extract well information
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

        # Format plate numbers
        df_variants_["Plate"] = df_variants_["Plate"].apply(
            lambda x: f"0{x}" if len(x) == 1 else x
        )

        # Rename columns
        df_variants_.rename(columns={
            "Variant": "nucleotide_mutation",
            "Substitutions": "amino_acid_substitutions",
            "nc_variant": "nt_sequence",
            "aa_variant": "aa_sequence"
        }, inplace=True)

        # Create restructured DataFrame
        restructured_df = df_variants_.copy()
        restructured_df.columns = restructured_df.columns.str.lower().str.replace(r'[\s-]', '_', regex=True)
        restructured_df.columns = restructured_df.columns.str.replace('p_adj._value', 'p_adj_value')

        # Select and order columns
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


def save_platemap_to_file(heatmaps, outputdir, name, show_msa, oligo_mode=False):
    """Save plate maps to file with support for oligo mode"""
    if not os.path.exists(os.path.join(outputdir, "Platemaps")):
        os.makedirs(os.path.join(outputdir, "Platemaps"))
    
    file_path = os.path.join(outputdir, "Platemaps", name)
    
    if oligo_mode:
        # Save oligo visualization
        hv.save(heatmaps, file_path + "_oligo.html")
    else:
        # Save standard visualization
        if show_msa:
            heatmaps.save(file_path + "_msa.html", embed=True)
        else:
            hv.renderer("bokeh").save(heatmaps, file_path)

def save_csv(df, outputdir, name, oligo_mode=False):
    """Save results to CSV with support for oligo mode"""
    if not os.path.exists(os.path.join(outputdir, "Results")):
        os.makedirs(os.path.join(outputdir, "Results"))
    
    if oligo_mode:
        file_path = os.path.join(outputdir, "Results", f"{name}_oligo_results.csv")
    else:
        file_path = os.path.join(outputdir, "Results", f"{name}.csv")
    
    df.to_csv(file_path, index=False)
    logging.info(f"Saved results to {file_path}")

# Function to process the reference CSV and generate variants
"""
def process_ref_csv(cl_args, tqdm_fn=tqdm.tqdm):    
    Process the reference CSV and generate variants or matches

    Args:
        cl_args (dict): Command line arguments
        tqdm_fn: Progress bar function
    
    Returns:
        tuple: (variant DataFrame, reference DataFrame)
    
    ref_df = pd.read_csv(cl_args["summary"])
    result_folder = create_result_folder(cl_args)
    variant_csv_path = os.path.join(result_folder, "variants.csv")

    # Initialize columns based on mode
    if cl_args.get("oligo"):
        variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "Reference_ID"])
    else:
        variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "variant"])
    
    unique_barcodes = ref_df['barcode_plate'].unique()
    for barcode in tqdm_fn(unique_barcodes, total=len(unique_barcodes), desc="Processing Samples"):
        barcode_plate = row["barcode_plate"]
        name = row["name"]
        refseq = row["refseq"].upper()

        name_folder = os.path.join(result_folder, name)
        os.makedirs(name_folder, exist_ok=True)
        
        # Create template fasta file
        temp_fasta_path = os.path.join(name_folder, f"temp_{name}.fasta")
        if not os.path.exists(temp_fasta_path):
            with open(temp_fasta_path, "w") as f:
                f.write(f">{name}\n{refseq}\n")
            logging.info(f"Created template fasta for {name}")
        
        # Filter barcodes
        barcode_path = filter_bc(cl_args, name_folder, i)
        output_dir = Path(result_folder) / f"{cl_args['name']}_fastq"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Perform demultiplexing if not skipped - same for both modes
        if not cl_args["skip_demultiplexing"]:
            try:

                file_to_fastq = cat_fastq_files(cl_args.get("path"), output_dir)
                demux_fastq(file_to_fastq, name_folder, barcode_path)
                logging.info(f"Demultiplexing completed for {name}")
            except Exception as e:
                logging.error(f"Demultiplexing failed for sample {name}: {str(e)}", exc_info=True)
                continue
        
        # Perform variant calling if not skipped
        if not cl_args["skip_variantcalling"]:
            try:
                # Call variants with appropriate mode
                variant_result = call_variant(
                    name,
                    name_folder,
                    temp_fasta_path,
                    barcode_path,
                    ref_sequences_df=ref_df if cl_args.get("oligo") else None,
                    oligo_mode=cl_args.get("oligo", False)
                )
                
                if variant_result is not None and not variant_result.empty:
                    variant_result["barcode_plate"] = barcode_plate
                    variant_result["name"] = name
                    if not cl_args.get("oligo"):
                        variant_result["refseq"] = refseq
                    variant_df = pd.concat([variant_df, variant_result])
                    logging.info(f"Variant calling completed for {name}")
                
            except Exception as e:
                logging.error(f"Variant calling failed for sample {name}: {str(e)}", exc_info=True)
                continue
    
    # Save complete variant results
    if not variant_df.empty:
        variant_df.to_csv(variant_csv_path, index=False)
        logging.info(f"Saved variant results to {variant_csv_path}")
    else:
        logging.warning("No variant results were generated")
    
    return variant_df, ref_df
"""

def create_oligo_heatmap(data, plate_number):
    """
    Create heatmap visualization for oligo mode data with 96-well plate format
    """
    import holoviews as hv
    from holoviews import opts
    import pandas as pd
    import numpy as np
    
    # Filter data for the specific plate
    plate_data = data[data['barcode_plate'] == plate_number].copy()
    
    # Create complete 96-well plate data frame
    rows = list('ABCDEFGH')
    cols = list(range(1, 13))
    complete_wells = []
    
    for row in rows:
        for col in cols:
            well = f"{row}{col}"
            # Find matching data for this well
            well_data = plate_data[plate_data['Well'] == well]
            
            if not well_data.empty:
                row_data = well_data.iloc[0]
                complete_wells.append({
                    'Row': row,
                    'Column': col,
                    'Alignment_Count': float(row_data['Alignment Count']),
                    'Match_Frequency': float(row_data['Match_Frequency']),
                    'Reference_ID': str(row_data['Reference_ID']),
                    'Total_Reads': int(row_data['Total_Reads']),
                    'Label': f"{row_data['Reference_ID']}\n({int(row_data['Alignment Count'])}/{row_data['Total_Reads']})",
                    'Well': well
                })
            else:
                # Add empty well
                complete_wells.append({
                    'Row': row,
                    'Column': col,
                    'Alignment_Count': 0.0,
                    'Match_Frequency': 0.0,
                    'Reference_ID': '',
                    'Total_Reads': 0,
                    'Label': '',
                    'Well': well
                })
    
    # Convert to DataFrame for HoloViews
    df = pd.DataFrame(complete_wells)
    
    # Create heatmap
    heatmap = hv.HeatMap(
        df,
        kdims=['Column', 'Row'],
        vdims=['Alignment_Count', 'Match_Frequency', 'Reference_ID', 'Total_Reads', 'Label', 'Well']
    )
    
    # Create labels
    labels = hv.Labels(
        {
            'Column': df['Column'],
            'Row': df['Row'],
            'text': df['Label']
        },
        ['Column', 'Row'], 
        'text'
    )
    
    # Style the heatmap
    heatmap = heatmap.opts(
        opts.HeatMap(
            width=800,
            height=600,
            tools=['hover'],
            colorbar=True,
            toolbar='above',
            title=f'Plate {int(plate_number)} - Oligo Matches',
            cmap='Blues',
            xlabel='Column',
            ylabel='Row',
            frame_width=800,
            frame_height=400,
            line_color='black',
            line_width=2,
            invert_yaxis=True,  # Invert y-axis to show A at top
            hover_tooltips=[
                ('Well', '@Well'),
                ('Reference ID', '@Reference_ID'),
                ('Alignment Count', '@Alignment_Count{0}'),
                ('Match Frequency', '@Match_Frequency{0.00}'),
                ('Total Reads', '@Total_Reads')
            ]
        )
    )
    
    # Style the labels
    labels = labels.opts(
        opts.Labels(
            text_font_size='8pt',
            text_color='black',
            text_align='center',
            text_baseline='middle'
        )
    )
    
    # Combine heatmap and labels
    plate_map = heatmap * labels
    
    return plate_map

def visualize_results(variant_df, result_folder, cl_args):
    """
    Create and save visualizations based on mode
    """
    if cl_args.get("oligo"):
        # Create visualization directory
        vis_dir = os.path.join(result_folder, "Platemaps")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Get unique plates
        plates = sorted(variant_df['barcode_plate'].unique())
        
        # Initialize holoviews
        hv.extension('bokeh')
        
        # Create visualization for each plate
        all_plates = []
        for plate in plates:
            try:
                plate_vis = create_oligo_heatmap(variant_df, plate)
                all_plates.append(plate_vis)
                logging.info(f"Created visualization for plate {plate}")
            except Exception as e:
                logging.error(f"Error creating visualization for plate {plate}: {str(e)}", exc_info=True)
                continue
        
        if not all_plates:
            logging.error("No plate visualizations could be created")
            return
        
        try:
            # Arrange plates in a grid
            if len(all_plates) > 1:
                layout = hv.Layout(all_plates).cols(2)
            else:
                layout = all_plates[0]
            
            # Save visualization
            output_path = os.path.join(vis_dir, f"{cl_args['name']}_oligo")
            hv.save(layout, output_path + '.html')
            logging.info(f"Saved visualization to {output_path}.html")
            
            # Save results CSV
            results_dir = os.path.join(result_folder, "Results")
            os.makedirs(results_dir, exist_ok=True)
            csv_path = os.path.join(results_dir, f"{cl_args['name']}_oligo_results.csv")
            variant_df.to_csv(csv_path, index=False)
            logging.info(f"Saved results to {csv_path}")
            
        except Exception as e:
            logging.error(f"Error saving visualizations: {str(e)}", exc_info=True)
            
    else:
        # Use existing non-oligo visualization logic
        df_variants, df_vis = create_df_v(variant_df)
        
        layout = generate_platemaps(
            max_combo_data=df_vis,
            result_folder=result_folder,
            show_msa=cl_args["show_msa"],
        )
        
        save_platemap_to_file(
            heatmaps=layout,
            outputdir=result_folder,
            name=cl_args["name"],
            show_msa=cl_args["show_msa"]
        )
        save_csv(df_variants, result_folder, cl_args["name"])

def run_LevSeq(cl_args, tqdm_fn=tqdm.tqdm):
    """
    Main function to run LevSeq pipeline with support for oligo pool sequencing
    """
    # Create result folder and configure logging
    result_folder = create_result_folder(cl_args)
    configure_logging(result_folder)
    logging.info("Starting LevSeq pipeline")
    
    try:
        # Process reference CSV
        ref_df = pd.read_csv(cl_args["summary"])
        
        # Get unique barcode plates
        unique_barcodes = sorted(ref_df['barcode_plate'].dropna().unique())
        print(f"Number of unique barcodes: {len(unique_barcodes)}")

        # Initialize empty variant DataFrame
        variant_df = pd.DataFrame()
        
        # Process each unique barcode plate
        for barcode_plate in tqdm_fn(unique_barcodes,total=len(unique_barcodes), desc="Processing plates"):
            # Get entries for this barcode plate
            plate_data = ref_df[ref_df['barcode_plate'] == barcode_plate]
            # Get the name (should be same for all entries with this barcode_plate)
            name = plate_data.iloc[0]['name']
            
            logging.info(f"Processing plate {barcode_plate} with name {name}")
            
            name_folder = os.path.join(result_folder, name)
            os.makedirs(name_folder, exist_ok=True)
            
            # Get the first reference sequence for this plate
            refseq = plate_data.iloc[0]['refseq'].upper()
            
            # Create template fasta file
            temp_fasta_path = os.path.join(name_folder, f"temp_{name}.fasta")
            if not os.path.exists(temp_fasta_path):
                with open(temp_fasta_path, "w") as f:
                    f.write(f">{name}\n{refseq}\n")
                logging.info(f"Created template fasta for {name}")
            
            # Filter barcodes
            idx = ref_df[ref_df['barcode_plate'] == barcode_plate].index[0]
            barcode_path = filter_bc(cl_args, name_folder, idx)
            
            # Set up output directory for fastq files
            output_dir = Path(result_folder) / f"{cl_args['name']}_fastq"
            output_dir.mkdir(parents=True, exist_ok=True)

            # Perform demultiplexing if not skipped
            if not cl_args["skip_demultiplexing"]:
                try:
                    file_to_fastq = cat_fastq_files(cl_args.get("path"), output_dir)
                    demux_fastq(file_to_fastq, name_folder, barcode_path)
                    logging.info(f"Demultiplexing completed for {name}")
                except Exception as e:
                    logging.error(f"Demultiplexing failed for {name}: {str(e)}")
                    continue

            # Perform variant calling if not skipped
            if not cl_args["skip_variantcalling"]:
                try:
                    if cl_args.get("oligo"):
                        variant_result = call_variant(
                            str(name),
                            name_folder,
                            temp_fasta_path,
                            barcode_path,
                            ref_sequences_df=ref_df,  # Pass only relevant entries
                            oligo_mode=True
                        )
                    else:
                        variant_result = call_variant(
                            str(name),
                            name_folder,
                            temp_fasta_path,
                            barcode_path
                        )

                    if variant_result is not None and not variant_result.empty:
                        variant_result["barcode_plate"] = barcode_plate
                        variant_result["name"] = name
                        variant_df = pd.concat([variant_df, variant_result])
                        logging.info(f"Variant calling completed for {name}")

                except Exception as e:
                    logging.error(f"Variant calling failed for {name}: {str(e)}")
                    continue

        # Save complete variant results
        if not variant_df.empty:
            variant_csv_path = os.path.join(result_folder, "variants.csv")
            variant_df.to_csv(variant_csv_path, index=False)
            logging.info(f"Saved variant results to {variant_csv_path}")
            
            # Create visualizations
            try:
                visualize_results(variant_df, result_folder, cl_args)
                logging.info("Visualization completed successfully")

            except Exception as e:
                logging.error(f"Error in visualization: {str(e)}", exc_info=True)
                # Continue execution even if visualization fails
        else:
            logging.warning("No variants were generated")
            
    except Exception as e:
        logging.error(f"Error in LevSeq pipeline: {str(e)}", exc_info=True)
        if not variant_df.empty:
            error_csv_path = os.path.join(result_folder, "variants_error.csv")
            variant_df.to_csv(error_csv_path, index=False)
            logging.info(f"Saved partial results to {error_csv_path}")
        raise
    
    logging.info("LevSeq pipeline completed successfully")
    return {
        'result_folder': result_folder,
        'variant_df': variant_df
    }
