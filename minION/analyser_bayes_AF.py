# Import
from minION.util import IO_processor
from minION import analyser
from minION import consensus
import importlib
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import gzip
import math
import re
import pickle
import itertools
import pysam
import subprocess



def get_bases_from_pileup_simulation(bam_file, chrom, positions):
    bases_dict = {position: {} for position in positions}
    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, min(positions) - 1, max(positions) + 1,
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            pos = pileup_column.pos + 1
            if pos in positions:
                for pileup_read in pileup_column.pileups:
                    read_name = pileup_read.alignment.query_name

                    # Handle deletions
                    if pileup_read.is_del:
                        base = '-'  # Symbol to represent a deletion
                    elif not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    else:
                        continue

                    # Add base to the dictionary
                    if read_name not in bases_dict[pos]:
                        bases_dict[pos][read_name] = base

    # Get unique read names and sort them
    read_names = sorted(set().union(*[bases_dict[pos].keys() for pos in bases_dict]))

    # Create DataFrame for bases
    df_bases = pd.DataFrame(index=read_names, columns=positions)
    
    # Populate DataFrame
    for pos in positions:
        for read_name in bases_dict[pos]:
            df_bases.at[read_name, pos] = bases_dict[pos][read_name]
    
    # Fill NaN with "-"
    df_bases = df_bases.fillna("-")

    return df_bases



def get_bases_from_pileup(bam_file, chrom, positions):
    bases_dict = {position: {} for position in positions}
    qualities_dict = {position: {} for position in positions}
    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, min(positions) - 1, max(positions) + 1,
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            pos = pileup_column.pos + 1
            if pos in positions:
                for pileup_read in pileup_column.pileups:
                    read_name = pileup_read.alignment.query_name

                    # Handle deletions
                    if pileup_read.is_del:
                        base = '-'  # or any symbol you prefer to represent a deletion
                        quality = 0  # Assign a default quality for deletions
                    elif not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                        quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                    else:
                        continue

                    # Add base and quality to the dictionaries
                    if read_name not in bases_dict[pos]:
                        bases_dict[pos][read_name] = base
                        qualities_dict[pos][read_name] = quality

    # Get unique read names and sort them
    read_names = sorted(set().union(*[bases_dict[pos].keys() for pos in bases_dict]))

    # Create DataFrames
    df_bases = pd.DataFrame(index=read_names, columns=positions)
    df_qualities = pd.DataFrame(index=read_names, columns=positions)
    
    # Populate DataFrames
    for pos in positions:
        for read_name in bases_dict[pos]:
            df_bases.at[read_name, pos] = bases_dict[pos][read_name]
            df_qualities.at[read_name, pos] = qualities_dict[pos][read_name]
    
    # Fill NaN with "-" for bases and 0 for qualities
    df_bases = df_bases.fillna("-")
    df_qualities = df_qualities.fillna(10) # 10 is the lowest quality filter we used for filtering

    return df_bases, df_qualities


def get_soft_pop_frequency(bam_file, template, reference, nb_positions, min_depth = 5):

    # Min depth based on the alphabet size
    
    # Check also for variant by sampling random positions

    bases_df = get_bases_from_pileup(bam_file, reference, nb_positions)
    
    frequency_df = bases_df.apply(get_variant_name, axis=1, args=(template, nb_positions)).value_counts().reset_index()

    frequency_df.columns = ['Population', 'N_reads']
    
    frequency_df["Frequency"] = frequency_df["N_reads"] / frequency_df["N_reads"].sum()


    # Filter for frequency > 0.4 and depth > 15
    frequency_df = frequency_df[(frequency_df["Frequency"] > min_freq) & (frequency_df["N_reads"] > min_depth)]

    return frequency_df


def add_neighbouring_positions(positions, nb_neighbours, max_index):
    if positions is np.nan:
        return np.nan
    
    elif isinstance(positions, int):
        positions = [positions]

    new_positions = []
    for position in positions:
        for new_pos in range(position - nb_neighbours, position + nb_neighbours + 1):
            if 1 <= new_pos <= max_index:  # Check if the new position is within valid index range
                new_positions.append(new_pos)
    return sorted(set(new_positions))

def calculate_mean_quality_for_reads(bases_df, qual_df, nb_positions, nb_neighbours):
    if isinstance(nb_positions, int):
        nb_positions = [nb_positions]  # Convert single integer to a list

    read_mean_qualities = {}

    for nb_position in nb_positions:
        if nb_position not in bases_df.columns or nb_position not in qual_df.columns:
            continue  # Skip positions that are not present in either DataFrame

        neighbor_positions = range(nb_position - nb_neighbours, nb_position + nb_neighbours + 1)

        for read_name in qual_df.index:
            total_qual = 0
            valid_count = 0

            for position in neighbor_positions:
                if position not in bases_df.columns:
                    continue  # Skip positions that are outside the DataFrame's columns

                base = bases_df.at[read_name, position]
                quality = qual_df.at[read_name, position]

                if base != "-" and not pd.isna(quality):
                    total_qual += quality
                    valid_count += 1

            if valid_count == 0:
                continue  # Skip if no valid qualities were found
            else:
                if read_name not in read_mean_qualities:
                    read_mean_qualities[read_name] = {}
                read_mean_qualities[read_name][nb_position] = total_qual / valid_count

    # Convert the dictionary into a DataFrame
    mean_quality_df = pd.DataFrame.from_dict(read_mean_qualities, orient='index')
    updated_base_df = bases_df[nb_positions] 

    return updated_base_df, mean_quality_df

def add_neighbouring_positions(positions, nb_neighbours, max_index):
    if positions is None:
        return None
    
    elif isinstance(positions, int):
        positions = [positions]

    new_positions = []
    for position in positions:
        for new_pos in range(position - nb_neighbours, position + nb_neighbours + 1):
            if 1 <= new_pos <= max_index:  # Check if the new position is within valid index range
                new_positions.append(new_pos)
    return sorted(set(new_positions))

def calculate_mean_quality_for_reads(bases_df, qual_df, nb_positions, nb_neighbours):
    if isinstance(nb_positions, int):
        nb_positions = [nb_positions]  # Convert single integer to a list

    read_mean_qualities = {}

    for nb_position in nb_positions:
        if nb_position not in bases_df.columns or nb_position not in qual_df.columns:
            continue  # Skip positions that are not present in either DataFrame

        neighbor_positions = range(nb_position - nb_neighbours, nb_position + nb_neighbours + 1)

        for read_name in qual_df.index:
            total_qual = 0
            valid_count = 0

            for position in neighbor_positions:
                if position not in bases_df.columns:
                    continue  # Skip positions that are outside the DataFrame's columns

                base = bases_df.at[read_name, position]
                quality = qual_df.at[read_name, position]

                if base != "-" and not pd.isna(quality):
                    total_qual += quality
                    valid_count += 1

            if valid_count == 0:
                continue  # Skip if no valid qualities were found
            else:
                if read_name not in read_mean_qualities:
                    read_mean_qualities[read_name] = {}
                read_mean_qualities[read_name][nb_position] = total_qual / valid_count

    # Convert the dictionary into a DataFrame
    mean_quality_df = pd.DataFrame.from_dict(read_mean_qualities, orient='index')
    updated_base_df = bases_df[nb_positions] 

    return updated_base_df, mean_quality_df

def get_non_error_prop(quality_score):
    """Convert quality score to non-error probability."""
    return 1 - 10 ** (-quality_score / 10)

def get_softmax_count_df(bases_df, qual_df, nb_positions):

    alphabet = "ACTG-"
    softmax_counts = {position: [] for position in nb_positions}
    
    for position in nb_positions:
        for base in alphabet:
            base_mask = bases_df[position] == base
            base_counts = base_mask.sum()
            # Calculate the non-error probability for each base and sum them up
            soft_count = sum(base_mask * qual_df[position].apply(get_non_error_prop))
            softmax_counts[position].append(soft_count)

    softmax_count_df = pd.DataFrame(softmax_counts, columns=nb_positions, index=list(alphabet))

    # Apply softmax to each column (position)
    softmax_count_df = softmax_count_df.apply(lambda x: x / x.sum(), axis=0)

    return softmax_count_df

def get_softmax_count_df_Simulation(bases_df, qual_df, nb_positions):
    
    alphabet = "ACTG-"
    softmax_counts = {position: [] for position in nb_positions}

    for position in nb_positions:
        for base in alphabet:
            base_mask = bases_df[position] == base
            base_counts = base_mask.sum()
            soft_count = sum(base_mask * 0.99)
            softmax_counts[position].append(soft_count)
    
    softmax_count_df = pd.DataFrame(softmax_counts, columns=nb_positions, index=list(alphabet))
    softmax_count_df = softmax_count_df.apply(lambda x: x / x.sum(), axis=0)
    return softmax_count_df
    
def get_softmax(soft_count):
    """Calculate the softmax of a dictionary of soft counts."""
    # Calculate the sum of the non-error probabilities
    total = sum(soft_count.values())
    # Calculate the softmax for each base
    return {base: count / total for base, count in soft_count.items()}

def call_potential_populations(softmax_df, ref_seq):
    positions = softmax_df.columns
    top_combinations = []
    
    # Get the top 2 variants for each position
    for position in positions:
        top_variants = softmax_df[position].nlargest(2)

        if top_variants.iloc[1] < 0.1:
            top_combinations.append([top_variants.index[0]])
        
        else:
            top_combinations.append(top_variants.index.tolist())

        potential_combinations = list(itertools.product(*top_combinations))

    
    variants = {"Variant" : [], "Probability" : []}
    
    for combination in potential_combinations:
        final_variant = []
        for i, pos in enumerate(positions):

            if combination[i] == ref_seq[pos - 1]:
                continue

            elif combination[i] == "-":
                var = f"{ref_seq[pos - 1]}{pos}DEL"
                final_variant.append(var)
            else:
                var = f"{ref_seq[pos - 1]}{pos}{combination[i]}"
                final_variant.append(var)

        final_variant = '_'.join(final_variant)
        if final_variant == "":
            final_variant = "#PARENT#"

        joint_prob = np.prod([softmax_df.at[combination[i], positions[i]] for i in range(len(positions))])
    
        variants["Variant"].append(final_variant)
        variants["Probability"].append(joint_prob)

    return variants

def get_variant_soft(bam_file, template_seq, ref_name, padding = 50):

    variants = {"Variant" : [], "Position" : [], "Alignment Probability" : [], "Alignment Count" : []}

    alignment_count = int(subprocess.run(f"samtools view -c {bam_file}", shell=True, capture_output=True).stdout.decode("utf-8").strip())

    # if alignment_count < 5:
    #     print("Not enough alignments")
    #     return None


    template = analyser.get_template_sequence(template_seq)

    padding_start, padding_end = padding, padding
    range_positions = range(padding_start + 1, len(template) - padding_end + 1) 

    freq_dist = pd.DataFrame(analyser.get_highest_non_ref_base_freq_2(bam_file, ref_name, range_positions, template, qualities=False)[0]).T.rename(columns={0:"Base", 1:"Frequency"})

    nb_positions = analyser.get_nb_positions(freq_dist, 0.3)

    available_positions = [pos for pos in range_positions if pos not in nb_positions]


    if len(nb_positions) == 0:
        # Select random 3 positions
        nb_positions = np.random.choice(available_positions, 3, replace=False)

    elif len(nb_positions) == 1:
        add_pos  = np.random.choice(available_positions, 2, replace=False)
        nb_positions = np.append(nb_positions, add_pos)

    elif len(nb_positions) > 15:
        print("Too many positions, either contaminated or sequencing error")
        nb_positions = np.random.choice(range_positions, 3, replace=False)

    #bases_df, qual_df = get_bases_from_pileup(bam_file, ref_name, add_neighbouring_positions(nb_positions, 2, len(template)))
    bases_df = get_bases_from_pileup_simulation(bam_file, ref_name, add_neighbouring_positions(nb_positions, 2, len(template)))
    #bases_df, qual_df = calculate_mean_quality_for_reads(bases_df, qual_df, nb_positions, 2)
    #qual_df = qual_df.fillna(10)
    #softmax_df = get_softmax_count_df(bases_df, qual_df, nb_positions)
    qual_df = pd.DataFrame() # Filler
    softmax_df = get_softmax_count_df_Simulation(bases_df, qual_df, nb_positions)
    variant_df = pd.DataFrame(call_potential_populations(softmax_df, template)).sort_values(by="Probability", ascending=False)

    # Take top variant
    variants["Variant"] = variant_df["Variant"].iloc[0]
    variants["Position"] = nb_positions
    variants["Alignment Probability"] = variant_df["Probability"].iloc[0]
    variants["Alignment Count"] = alignment_count

    
    # print("Error in getting variant")
    # variants["Variant"] = np.nan
    # variants["Position"] = None
    # variants["Alignment Probability"] = None
    # variants["Alignment Count"] = alignment_count

    return variants


def get_variant_df_soft(demultiplex_folder: Path, ref_seq : Path, ref_name : str, barcode_dicts : dict = None, merge = True, min_depth= 5, padding=50, rowwise = False, alignment_name = "alignment_minimap.bam"):


    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = analyser.template_df(barcode_dicts, rowwise=False)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Alignment Probability": [], "Alignment Count": []}

    template = analyser.get_template_sequence(ref_seq) # Reference sequence

    summary = analyser.read_summary_file(demultiplex_folder)
    n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index() 



    for barcode_id, barcode_dict in barcode_dicts.items():

        rbc = os.path.basename(barcode_id)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)
            print("Processing", rbc, fbc)
            count = n_counts[(n_counts["RBC"] == rbc) & (n_counts["FBC"] == fbc)]["count"].values[0]

            # # If alignment file exist continue
            if not os.path.exists(os.path.join(front_barcode, "alignment_minimap.bam")):
                print(f"Alignment file in {front_barcode} does not exist, running alignment and indexing")
                analyser.run_alignment_and_indexing(ref_seq, front_barcode, site_saturation=True)
            
            else: 
                print("Alignment file already exists, skipping alignment and indexing")


            bam_file = front_barcode / alignment_name


            if not bam_file.exists() or count < min_depth:
                print(f"{bam_file} does not exist.")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(np.nan)
                variants["Variant"].append(np.nan)
                variants["Alignment Count"].append(np.nan)
                variants["Alignment Probability"].append(np.nan)
                print(f"Skipping Variant: {fbc}/{rbc}")
                continue

            # try:
            if padding == 0:
                print("Padding is 0. Implementing soft alignment")
                nn_variants = get_variant_soft(bam_file, ref_seq, ref_name, padding = padding)
            
            else: 
                nn_variants = get_variant_soft(bam_file, ref_seq, ref_name, padding = padding)
                print(nn_variants)

            if nn_variants is None:
                print("Empty variant list")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(np.nan)
                variants["Variant"].append(np.nan)
                variants["Alignment Count"].append(np.nan)
                variants["Alignment Probability"].append(np.nan)
                print(f"Skipping Variant: {fbc}/{rbc}")
                continue

            
            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(nn_variants["Position"])
            variants["Variant"].append(nn_variants["Variant"])
            #variants["Reads"].append(count)
            # Check if Alignment count is a number
            if isinstance(nn_variants["Alignment Count"], int) & (isinstance(nn_variants["Alignment Probability"], float) or nn_variants["Alignment Probability"] == "-"):
                variants["Alignment Count"].append(nn_variants["Alignment Count"])
                variants["Alignment Probability"].append(nn_variants["Alignment Probability"])


            else:
                print(f"Skipping {rbc}/{fbc} due incomplete data")
                variants["Alignment Count"].append(np.nan)
                variants["Alignment Probability"].append(np.nan)

            print(f"Variant: {fbc}/{rbc} {nn_variants['Alignment Count']} {nn_variants['Alignment Probability']}")
        
        # except Exception as e:
        #     # Append 'NA' in case of an exception
        #     print(f"Error processing {rbc}/{fbc}: {e}")
        #     variants["RBC"].append(rbc)
        #     variants["FBC"].append(fbc)
        #     variants["Position"].append(None)
        #     variants["Variant"].append("NA")
        #     variants["Alignment Count"].append("NA")
        #     variants["Alignment Frequency"].append("NA")



    if merge:
        variant_template_df = analyser.template_df(barcode_dicts, rowwise=rowwise)
        variant_df = analyser.rename_barcode(pd.DataFrame(variants).merge(n_counts, on=["RBC","FBC"] , how="left"))
        variant_df["Variant"] = variant_df["Variant"].apply(analyser.format_variant_list)
        variant_df["Variant"] = variant_df["Variant"].apply(lambda x: analyser.adjust_variant(x, padding))

        return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    else:
        return variants

def call_variant_BF(bam_file, chrom, positions, reference_sequence, qualities = False):
    """ Calls variants from a BAM file.
    Args:
        - bam_file (str): Path to the BAM file.
        - chrom (str): Chromosome or contig name.
        - positions (list): List of positions to call variants.
        - reference_sequence (str): Reference sequence of the chromosome or contig.
    Returns:
        - variants (list): List of variants.
    """

    variants = {"Variant" : [], "Position" : [], "Frequency" : []}

    if qualities:
        bases, qualities = get_highest_non_ref_base_freq(bam_file, chrom, positions, reference_sequence)
    
    bases = analyser.get_highest_non_ref_base_freq_2(bam_file, chrom, positions, reference_sequence, qualities=False)

    for position in positions:
        ref_base = reference_sequence[position - 1].upper()
        non_ref_base, freq = bases[0][position]
        if non_ref_base and freq >= 0.35:

            if non_ref_base == "-":
                non_ref_base = "DEL"

            variant = f"{ref_base}{position}{non_ref_base}"

            variants["Variant"].append(variant)
            variants["Position"].append(int(position))
            variants["Frequency"].append(freq)
            #variants["Quality-Score"].append(qualities[position])
    
    if variants["Variant"] == []:
        variants["Variant"].append("#PARENT#")
        variants["Position"].append(np.nan)
        variants["Frequency"].append(np.nan)
        #variants["Quality-Score"].append("-")
            

    return variants
