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
from minION import *

# Import external packages
from pathlib import Path
import numpy as np
import pandas as pd
from importlib import resources
import subprocess
from Bio import SeqIO
import tqdm
import os
import re

# Get barcode used
def barcode_user(cl_args,i):
    # Set some default values if user did not provide barcodes
    fmin = 1
    fmax = 96
    bc_df = pd.read_csv(cl_args["summary"])
    rbc = bc_df["barcode_plate"][i]

    return int(fmin), int(fmax), int(rbc)


# Get output directory
def get_input_folder(cl_args):
    input_folder = IO_processor.check_data_folder(cl_args['path'])
    return input_folder


# Get fastq input directory, this is the basecalled folder
def fastq_path(folder):
    return IO_processor.find_folder(folder, "fastq_pass")


# Create result folder
def create_result_folder(cl_args):
    basecall_model = 'sup'
    result_folder = IO_processor.create_folder(
        cl_args['name'],
        basecall_model,
        target_path=Path(cl_args['output']))
    return result_folder


# Basecall reads
def basecall_reads(cl_args):
    print('basecalling')


# Filter barcode
def filter_bc(cl_args, result_folder, i):
    front_min, front_max, rbc = barcode_user(cl_args, i)
    # Obtain path of executable from package
    with resources.path('minION.barcoding', 'minion_barcodes.fasta') as barcode_path:
        front_prefix = "NB"
        back_prefix = "RB"
        bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
    barcode_path_filter = os.path.join(result_folder, "minion_barcodes_filtered.fasta")
    bp.filter_barcodes(barcode_path_filter, (front_min, front_max), rbc)
    return barcode_path_filter


# Filter template sequence length
def filter_seq(cl_args):
    return seq_min, seq_max


# Get reference fasta (parent sequence)
def parent_fasta(cl_args):
    template_fasta = cl_args['refseq']
    return template_fasta


# Demultiplex the basecalled fastq into plate-well folders
def demux_fastq(file_to_fastq, result_folder, barcode_path):
    # Plan B to locate using relative path relying on the git folder
    current_file_dir = Path(__file__).parent
    # Obtain path of executable from package
    with resources.path('minION.barcoding', 'demultiplex-x86') as executable_path:
        # Get min and max sequence length if user specified, otherwise use default
        seq_min = 800
        seq_max = 5000
        # Use subprocess to run the executable
        prompt = f"{str(executable_path)} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100} -m {seq_min} -x {seq_max}"
        subprocess.run(prompt, shell=True)


# Variant calling using VariantCaller class and generate dataframe
def call_variant(experiment_folder, template_fasta, demultiplex_folder_name):
    vc = VariantCaller(experiment_folder,
                       template_fasta,
                       demultiplex_folder_name=demultiplex_folder_name,
                       padding_start=0,
                       padding_end=0)

    variant_df = vc.get_variant_df(threshold=0.2,
                                   min_depth=5)
    return variant_df


# Saving heatmaps and csv in the results folder
def save_platemap_to_file(heatmaps, outputdir, name):
    if not os.path.exists(os.path.join(outputdir, "Platemaps")):
        os.makedirs(os.path.join(outputdir, "Platemaps"))
    file_path = os.path.join(outputdir, "Platemaps", name)
    hv.renderer('bokeh').save(heatmaps, file_path)


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
    df_variants_['Variant'].tolist()
    df_variants_['Variant'] = df_variants_['Variant'].replace(np.nan, '', regex=True)

    # Create nc_variant column
    df_variants_['nc_variant'] = df_variants_.apply(lambda row: create_nc_variant(row['Variant'], row['refseq']), axis=1)
    
    # Translate nc_variant to aa_variant
    df_variants_['aa_variant'] = df_variants_['nc_variant'].apply(lambda x: 'Deletion' if x == 'Deletion' else translate(x))
    # Fill in 'Deletion' in 'aa_variant' column
    df_variants_.loc[df_variants_['nc_variant'] == 'Deletion', 'aa_variant'] = 'Deletion'

    # Compare aa_variant with translated refseq and generate mutations column
    df_variants_['Mutations'] = df_variants_.apply(get_mutations, axis=1)

    # Fill in empty empty values
    df_variants_['Alignment Probability'] = df_variants_['Average mutation frequency'].fillna(0.0)
    df_variants_['Alignment Count'] = df_variants_['Alignment_count'].fillna(0.0)

    # Fill in parents into mutations Column
    for i in df_variants_.index:
        if df_variants_['nc_variant'].iloc[i] == 'Deletion':
            df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace('', '-')
        if df_variants_['Average mutation frequency'].iloc[i] == 0.0 and df_variants_['Mutations'].iloc[i] == '':
            df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace('', '#N.A.#')
        if df_variants_['Mutations'].iloc[i] == '':
            df_variants_.Mutations.iat[i] = df_variants_.Mutations.iat[i].replace('', '#PARENT#')

    # Add row and columns
    Well = df_variants_['Well'].tolist()
    column = [Well[i].strip('ABCDEFGH') for Well[i] in Well]
    row = [Well[i].rstrip('0123456789') for Well[i] in Well]
    df_variants_['Row'] = row
    df_variants_['Column'] = column
    df_variants_['Plate'] = df_variants_['name'].astype(str)
    
    # Update 'Plate' column from '1'-'9' to '01'-'09'
    df_variants_['Plate'] = df_variants_['Plate'].apply(lambda x: f'0{x}' if len(x) == 1 else x)

    # Select the desired columns in the desired order
    restructured_df = df_variants_[['barcode_plate', 'Plate', 'Well', 'Variant', 'Alignment Count', 'Average mutation frequency', 'P value', 'P adj. value', 'Mutations', 'nc_variant', 'aa_variant']]
    # Set 'Mutations' and 'Variant' columns to '#N.A.#' if 'Alignment Count' is smaller than 5
    restructured_df.loc[restructured_df['Alignment Count'] < 10, ['Mutations', 'Variant']] = '#N.A.#'

    return restructured_df, df_variants_

def create_nc_variant(variant, refseq):
    if isinstance(variant, np.ndarray):
        variant = variant.tolist()
    if variant == '' or pd.isnull(variant):
        return refseq
    elif variant == '#PARENT#':
        return refseq
    elif 'DEL' in variant:
        return 'Deletion'
    else:
        mutations = variant.split('_')
        nc_variant = list(refseq)
        for mutation in mutations:
            if len(mutation) >= 2:
                position = int(re.findall(r'\d+', mutation)[0]) - 1
                original = mutation[0]
                new = mutation[-1]
            if position < len(nc_variant) and nc_variant[position] == original:
                nc_variant[position] = new
        return ''.join(nc_variant)

def get_mutations(row):
    refseq_aa = translate(row['refseq'])
    variant_aa = row['aa_variant']

    if variant_aa == 'Deletion':
        return ''
    else:
        mutations = []
        if len(refseq_aa) == len(variant_aa):
            for i in range(len(refseq_aa)):
                if refseq_aa[i] != variant_aa[i]:
                    mutations.append(f"{refseq_aa[i]}{i+1}{variant_aa[i]}")
            if not mutations:
                return '#PARENT#'
        else:
            return 'LEN'
    return '_'.join(mutations) if mutations else ''

# Process the summary file
def process_ref_csv(cl_args):
    ref_df = pd.read_csv(cl_args['summary'])
    result_folder = create_result_folder(cl_args)

    variant_csv_path = os.path.join(result_folder, "variants.csv")
    if os.path.exists(variant_csv_path):
        variant_df = pd.read_csv(variant_csv_path)
    else:
        variant_df = pd.DataFrame(columns=["barcode_plate", "name", "refseq", "variant"])
    for i, row in ref_df.iterrows():
        barcode_plate = row["barcode_plate"]
        name = row["name"]
        refseq = row["refseq"]

        # Create a subfolder for the current iteration using the name value
        name_folder = os.path.join(result_folder, name)
        os.makedirs(name_folder, exist_ok=True)

        # Write the refseq to a temporary fasta file
        temp_fasta_path = os.path.join(name_folder, f"temp_{name}.fasta")
        with open(temp_fasta_path, "w") as f:
            f.write(f">{name}\n{refseq}\n")
        barcode_path = filter_bc(cl_args, name_folder, i)
        file_to_fastq = fastq_path(get_input_folder(cl_args))

        if not cl_args['skip_demultiplexing']: 
            demux_fastq(file_to_fastq, name_folder, barcode_path)
        
        if not cl_args['skip_variantcalling']: 
            variant_result = call_variant(result_folder, temp_fasta_path, f"{name}")
            variant_result["barcode_plate"] = barcode_plate
            variant_result["name"] = name
            variant_result["refseq"] = refseq

            variant_df = pd.concat([variant_df, variant_result])
        
        # Remove the temporary fasta file
        os.remove(temp_fasta_path)
    variant_df.to_csv(variant_csv_path, index=False)
    return variant_df

# Run MinION    

def run_MinION(cl_args, tqdm_fn=tqdm.tqdm):
    # Find specific experiment in the upper directory of nanopore data
    experiment_folder = get_input_folder(cl_args)

    # Find fastq from experiment folder
    file_to_fastq = fastq_path(experiment_folder)
    
    # Basecall if asked
    if cl_args["perform_basecalling"]:
        basecall_reads(cl_args)
    # Process summary file by row
    variant_df = process_ref_csv(cl_args)
    
    # Check if variants.csv already exist
    result_folder = create_result_folder(cl_args) 
    variant_csv_path = os.path.join(result_folder, "variants.csv")
    if os.path.exists(variant_csv_path):
        variant_df = pd.read_csv(variant_csv_path)
        df_variants,df_vis = create_df_v(variant_df)
    # Clean up and prepare dataframe for visualization
    else:
        df_variants,df_vis = create_df_v(variant_df)

    processed_csv = os.path.join(result_folder, 'processed.csv')
    df_variants.to_csv(processed_csv, index = False)
    # Generate heatmap
    hm_ = generate_platemaps(df_vis)

    # Saving heatmap and csv
    save_platemap_to_file(hm_, result_folder, cl_args['name'])
    save_csv(df_variants, result_folder, cl_args['name'])
