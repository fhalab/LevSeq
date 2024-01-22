from minION.util.IO_processor import read_fasta_file, find_file, get_barcode_dict
from minION.util.globals import CODONS
from uncertainties import ufloat 
import pandas as pd
from pathlib import Path
import os
from statistics import mean
from Bio import Align
from Bio.Seq import Seq
import pysam 
from collections import Counter
import subprocess
import numpy as np
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging




def get_sequence_distribution(path : Path) -> list:
    """Get the sequence distribution from one or multiple fastq files
    Args: 
        - path, path to fastq file
    Return:
        - List of sequence distribution
    """

    if path is None:
        return [0]

    sequences = read_fasta_file(path, True)

    return sequences["Sequence"]

def get_consensus_sequence(path : Path, score = True) -> dict:
    """Get the consensus sequence with the Phred Quality Score
    Input:  - path, where the fasta file is located
            - score, if True, return sequence and quality scores
    Output: Dictionary with sequence only or sequence and quality scores"""

    if path is None:
        return {"Sequence" : "NA", "Quality-Score" : "NA"}

    sequences = read_fasta_file(path, True)

    return sequences

def get_template_sequence(path : Path) -> str:
    """Read template sequence fasta file
    Input:  - path, where the fasta file is located
    Output: Template sequence"""
    
    template = read_fasta_file(path)
    
    return template["Sequence"][0]

def get_n_reads(path : Path) -> int:
    """Get the number of reads from a fastq file
    Input:  - path, where the fastq file is located
    Output: Number of reads"""

    n_reads = 0

    with open(path, "r") as fastq:
        for line in fastq:
            if line.startswith("@"):
                n_reads += 1

    return n_reads


def check_sequence(sequence : str = "") -> bool:
    """Check if the sequence is a multiple of 3, so that it can be translated into a protein sequence. Otherwise the sequence is not valid. Reason can be low quality, ...
    Input:  - sequence, DNA/RNA sequence
    Output: - True if sequence is a multiple of 3, False otherwise"""
    if len(sequence) % 3 == 0:
        return True
    else:
        return False

def align_sequences(consensus_path : Path, template_path : Path) -> dict:
    """Align consensus sequences to template sequence.
    Args:  
        - template_path, path to template sequence
        - consensus_path, path to consensus sequence
    Return: 
        - Dictionary with aligned sequences
    """

    template = get_template_sequence(template_path)

    consensus = get_consensus_sequence(consensus_path, True)

    return bioalign(template, consensus), template


def swap_NN(template : str, consensus : str , min_score : int = 1):
    """Swap Nucleotides in consensus with template if the score is below min_score

    Input:  
            - template, template sequence
            - consensus, consensus sequence 
            - min score, minimum quality score 

    Output:
            - Dictionary with swapped sequences"""

    seq = list(consensus["Sequence"][0])  # Convert to list for efficient in-place modifications
    quality = consensus["Quality-Score"]


    if len(template) != len(seq):
        return {"Sequence": "NA", "Quality-Score": "NA"}

    for i in range(len(template)):
        if template[i] != seq[i] and quality[i] < min_score:
            print(f"Swapping {seq[i]} with {template[i]} at position {i} because quality score is {quality[i]}")
            seq[i] = template[i]

    # Assuming mean_quality_score returns the average score for the sequence
    aa_scores = mean_quality_score(quality)

    return {"Sequence": "".join(seq), "Quality-Score": aa_scores}

def mean_quality_score(quality_scores):
    """Calculate the mean quality score of a sequence
    Input:  
        - quality_scores, List of Ascii Phred Scores
    Output: 
        - Mean quality score for each codon"""

    if "NA" in quality_scores:
        return ["NA"]
    
    aa_scores = []

    for i in range(0, len(quality_scores), 3):

        triplet = quality_scores[i:i+3]


        mean_score = round(sum(triplet) / len(triplet), 3)

        aa_scores.append(mean_score)

    return aa_scores


def AA_seq(consensus_path, template_fasta):
    """Translate consensus sequences into protein sequences.
    Input:  - consensus_path, path to consensus sequence
            - template_fasta, path to template sequence
    Output: -Dictionary with protein sequences and their corresponding quality scores"""

    consensus_aligned, template = align_sequences(consensus_path, template_fasta)

    consensus_swapped = swap_NN(template, consensus_aligned) # Swap if quality score is higher than min threshold

    consensus_aa = translate_sequence(consensus_swapped["Sequence"])

    aa_quality = consensus_swapped["Quality-Score"]

    return {"Sequence": [consensus_aa["Protein-Sequence"]], "Quality-Score": aa_quality}


def translate_sequence(sequences : list = None):
    """
    Translate a list of sequences into protein sequences.
    Input:  - sequences, a list of DNA sequences
    Output: Dictionary containing DNA sequences and their corresponding protein sequences
    """
    if sequences is None:
        return proteins

    # Check if sequences is a list
    if not isinstance(sequences, list):
        sequences = [sequences]

    proteins = {"DNA-Sequence": [],
                "Protein-Sequence": []}
    
    for sequence in sequences:
        
        proteins["DNA-Sequence"].append(sequence)
        
        if check_sequence(sequence):
            proteins["Protein-Sequence"].append("".join([CODONS[sequence[i:i+3]] for i in range(0, len(sequence), 3) if i+3 <= len(sequence)]))
        else:
            print("Not a valid sequence", proteins["DNA-Sequence"])
            proteins["Protein-Sequence"].append("NA")

    return proteins

def uncertainty(x,y):
    """Boolean outpur if y is within the lower and upper boundry of x
    Input:  - x, e.g Template sequence
            - y, e.g the consensus sequence"""

def bioalign(template, consensus):
    """Align the template and consensus sequence"""


    seq1 = Seq(template)
    seq2 = Seq(consensus["Sequence"][0])
    score = consensus["Quality-Score"][0]

    thres = 5

    n_template = ufloat(len(template), thres)

    if len(seq2) > len(seq1):
        print("Query sequence is longer than target sequence. Cannot align.")

    elif uncertainty(n_template, len(seq2)) == False:
        print("Consensus and template sequence have different lengths. Cannot align.")
        return {"Sequence" : "NA", "Quality-Score" : "NA"}

    aligner = Align.PairwiseAligner()

    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = -100  # Penalize end gaps in the query sequence
    aligner.mismatch_score = 0
    aligner.match_score = 2
    aligner.mode = "global"

    alignments = aligner.align(seq1,seq2)[0]


    
    aligned_seq = "".join([s1 if s2 == '-' else s2 for s1, s2 in zip(alignments.target, alignments.query)])

    score_final = []
    score_index = 0
    for char in alignments.query:
        if char == '-':
            score_final.append(0)
        else:
            score_final.append(score[score_index])
            score_index += 1

    return {"Sequence" : [aligned_seq], "Quality-Score" : score_final}


# def call_variant(template, consensus, quality_score):
#     """Compare the template and consensus sequence and call variants
#     Input:  
#         - template, template sequence
#         - consensus, consensus sequence
#         - quality_score, quality score of the consensus sequence
#     Output: 
#         - Dictionary with Mutaiton and at which Position"""
    
    
#     if len(template) != len(consensus) or consensus == "NA" or quality_score == "NA":
#         return {"Variant" : "NA", "Position" : "NA", "Quality-Score" : "NA"}

#     variants = {"Variant" : [], "Position" : [], "Quality-Score" : []}
    
#     for i in range(min(len(template), len(consensus))):
#         if template[i] != consensus[i]:
#             pos = i + 1
#             variants["Variant"].append(f"{template[i]}{pos}{consensus[i]}") 
#             variants["Position"].append(pos)
#             variants["Quality-Score"].append(quality_score[i])
        
#     if variants["Variant"] == []:
#         variants["Variant"].append("#PARENT")
#         variants["Position"].append("-")
#         variants["Quality-Score"].append("-")
#     return variants

def call_variant_nn(template, consensus, quality_score):

    if len(template) != len(consensus) or consensus == "NA" or quality_score == "NA":
        return {"Variant" : "NA", "Position" : "NA", "Quality-Score" : "NA"}

    variants = {"Variant" : [], "Position" : [], "Quality-Score" : []}

    for i in range(min(len(template), len(consensus))):
        if template[i] != consensus[i]:
            pos = i + 1
            variants["Variant"].append(f"{template[i]}{pos}{consensus[i]}") 
            variants["Position"].append(pos)
            variants["Quality-Score"].append(quality_score[i])
        elif template[i] != consensus[i] and template[i] == "N":
            pass
        
    if variants["Variant"] == []:
        variants["Variant"].append("#PARENT#")
        variants["Position"].append("-")
        variants["Quality-Score"].append("-")

    return variants

def call_variant_cust(bam_file, chrom, positions, reference_sequence):
    """ Calls variants from a BAM file.
    Args:
        - bam_file (str): Path to the BAM file.
        - chrom (str): Chromosome or contig name.
        - positions (list): List of positions to call variants.
        - reference_sequence (str): Reference sequence of the chromosome or contig.
    Returns:
        - variants (list): List of variants.
    """

    variants = {"Variant" : [], "Position" : [], "Quality-Score" : []}

    bases, qualities = get_highest_non_ref_base_freq(bam_file, chrom, positions, reference_sequence)

    for position in positions:
        ref_base = reference_sequence[position - 1].upper()
        non_ref_base, freq = bases[position]
        if non_ref_base and freq >= 0.35:

            if non_ref_base == "-":
                non_ref_base = "DEL"

            variant = f"{ref_base}{position}{non_ref_base}"

            variants["Variant"].append(variant)
            variants["Position"].append(int(position))
            variants["Quality-Score"].append(qualities[position])
    
    if variants["Variant"] == []:
        variants["Variant"].append("#PARENT#")
        variants["Position"].append("-")
        variants["Quality-Score"].append("-")
            

    return variants



def read_summary_file(demultiplex_folder, summary_file = "barcoding_summary", file_type = ".txt"):
    """Read barcoding summary files """
    path = find_file(demultiplex_folder, summary_file, file_type)

    return pd.read_csv(path, sep = "\t")


def convert_ascii_to_phred(ascii):
    """Convert ascII phred scores to numerical Phred Score (Q). To get the numerical score, the ASCII value of the character is subtracted by 33.
    Input:  - ascii, ASCII character
    Output: Numerical Phred Score"""
    return ord(ascii) - 33


def barcode_to_well(barcode,):
    match = re.search(r'\d+', barcode)
    if match:
        number = int(match.group())
        rows = 'ABCDEFGH'
        row = rows[(number-1) // 12]
        col = (number-1) % 12 + 1
        return f"{row}{col}"
    else:
        return "NA"



def rename_barcode_guppy(variant_df):
    variant_df = variant_df.rename(columns={'RBC': 'Plate', 'FBC': 'Well'})

    # Extracting numbers from the 'Plate' column and converting to int
    variant_df["Plate"] = variant_df['Plate'].str.extract('(\d+)').astype(int)

    # Applying the barcode_to_well function to the 'Well' column
    variant_df["Well"] = variant_df["Well"].apply(barcode_to_well)

    return variant_df

def rename_barcode(variant_df):
    variant_df = variant_df.rename(columns={'RBC': 'Plate', 'FBC': 'Well'})

    # Extracting numbers from the 'Plate' column and converting to int
    variant_df["Plate"] = variant_df['Plate'].str.extract('(\d+)').astype(int)

    # Applying the barcode_to_well function to the 'Well' column
    variant_df["Well"] = variant_df["Well"].apply(barcode_to_well)

    return variant_df

def barcode_arrangements(summary_path : Path, column = "RBC") -> pd.DataFrame:
    """
    Get number of reads for each barcode/well

    Input:
        - Path to Summary file (e.g Demultiplex folder)
        - column, column name of the barcode
    
    Output:
        - DataFrame: Contains the number of reads for each barcode/well.

    """

    barcode_counts = (
        read_summary_file(summary_path)[column]
        .value_counts()
        .reset_index()
        .rename(columns={"count": "N_reads"})
    )

    return barcode_counts[barcode_counts[column] != "unclassified"]

# def barcode_arrangements(rbc, column = "RBC"):
#     """
#     Get the number of reads for each barcode/well.
    
#     Input:
#     - rbc (str): reverse barcode.
    
#     Output:
#     - DataFrame: Contains the number of reads for each barcode/well.
#     """
    
#     barcode_counts = (
#         read_summary_file(rbc)[column]
#         .value_counts()
#         .reset_index()
#         .rename(columns={"count": "N_reads"})
#     )

#     return barcode_counts[barcode_counts[column] != "unclassified"]

def count_reads_in_fastq(filename):
    """
    Count the number of reads in a FASTQ file.

    Parameters:
    - filename: The name of the FASTQ file.

    Returns:
    - The number of reads in the FASTQ file.
    """
    with open(filename, 'r') as f:
        return sum(1 for line in f) // 4


def quality_score(qulity_score, position):
    """Return the quality score at the mutated positions
    Input:  
        - quality_score, list of quality scores
        - position, list of positions
    Output:
        - quality score at the mutated positions"""
    
    return [qulity_score[i] for i in position]

def template_df(barcode_dicts : dict = None, rowwise = True):
    """To have coherent df for each experiment, a template df is created. The template also have the desired plates and columns in the desired order
    Input:
        - demultiplex_folder, folder where the demultiplexed files are located
        - rowwise, if True, the reverse barcodes are rows and not plates
        """
    
    if barcode_dicts is None:
        raise ValueError("No barcode dictionary provided")
    
    n_rbc = len(barcode_dicts.items())

    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    columns = [i for i in range(1, 13)]

    if rowwise:
        template = {"Plate": [], "Well": []}

        for row in rows:
            for column in columns:
                template["Plate"].append(1)
                template["Well"].append(f"{row}{column}")
    
    else:

        template = {"Plate": [], "Well": []}

        for i in range(n_rbc):
            for row in rows:
                for column in columns:
                    template["Plate"].append(i+1)
                    template["Well"].append(f"{row}{column}")
        
    return pd.DataFrame(template)

        
def get_variant_df(demultiplex_folder: Path, ref_seq : Path, barcode_dicts : dict = None, consensus_folder_name = "consensus" , sequences = False):
    """Call Variants from consensus sequences and return Variants for each barcode/well

    Input:  
        - Demultiplex folder (Path), folder where the demultiplexed files are located
        - ref_seq (Path), path to reference sequence
        - barcode_dicts (dict), dictionary with barcodes
        - sequences (bool), if True, return the consensus sequence

    Output: 
        - DataFrame of Variants for each reverse Barcode (RBC = Plate) and forward Barcode (FBC = Well). 
            For each Well, the variants are called and the number of reads for each barcode/well is counted
    """

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = template_df(barcode_dicts)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Quality-Score": []}

    if sequences:
        variants["Sequence"] = []

    template = get_template_sequence(ref_seq) # Reference sequence
    template_aa = translate_sequence([template])

    #reads_df = pd.DataFrame(columns=["RBC", "FBC"])

    # summary = read_summary_file(demultiplex_folder)
    
    # n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index()

    for barcode_id, barcode_dict in barcode_dicts.items():

        rbc = os.path.basename(barcode_id)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)

            front_barcode = os.path.join(front_barcode, consensus_folder_name)

            fasta_file = os.path.join(front_barcode, "consensus.fastq")

            # Check if consensus file exists
            if not os.path.exists(fasta_file):
                #print(f"Consensus file in {front_barcode} does not exist, skipping {fbc}")
                continue

            #consensus_aa = AA_seq(fasta_file, ref_seq)

            try:
                consensus = get_consensus_sequence(fasta_file, True)
                quality = mean_quality_score(consensus["Quality-Score"][0])
                consensus_aa = translate_sequence(consensus["Sequence"])
            
            except:
                print(f"Skipping {rbc}/{fbc}")
                print(consensus)
                continue

            if sequences:
                variants["Sequence"].append(consensus_aa["Protein-Sequence"][0])

            aa_variants = call_variant(template_aa["Protein-Sequence"][0], consensus_aa["Protein-Sequence"][0], quality)



            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(aa_variants["Position"])
            variants["Variant"].append(aa_variants["Variant"])
            variants["Quality-Score"].append(aa_variants["Quality-Score"])
    

    #variant_df = rename_barcode(pd.DataFrame(variants).merge(n_counts, on=["RBC","FBC"] , how="left"))

    # return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    return variants

def get_variant_df_nn(demultiplex_folder: Path, ref_seq : Path, barcode_dicts : dict = None, consensus_folder_name = "consensus" , sequences = False, merge = True):
    """Call Variants from consensus sequences and return Variants for each barcode/well

    Input:  
        - Demultiplex folder (Path), folder where the demultiplexed files are located
        - ref_seq (Path), path to reference sequence
        - barcode_dicts (dict), dictionary with barcodes
        - sequences (bool), if True, return the consensus sequence

    Output: 
        - DataFrame of Variants for each reverse Barcode (RBC = Plate) and forward Barcode (FBC = Well). 
            For each Well, the variants are called and the number of reads for each barcode/well is counted
    """

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = template_df(barcode_dicts, rowwise=False)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Quality-Score": []}

    if sequences:
        variants["Sequence"] = []

    template = get_template_sequence(ref_seq) # Reference sequence

    # summary = read_summary_file(demultiplex_folder)
    # n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index()


    for barcode_id, barcode_dict in barcode_dicts.items():

        rbc = os.path.basename(barcode_id)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)

            front_barcode = os.path.join(front_barcode, consensus_folder_name)

            fasta_file = os.path.join(front_barcode, "consensus.fastq")

            # Check if consensus file exists
            if not os.path.exists(fasta_file):
                #print(f"Consensus file in {front_barcode} does not exist, skipping {fbc}")
                continue

            try:
                consensus = get_consensus_sequence(fasta_file, True)
            
            except:
                print(f"Skipping {rbc}/{fbc}")
                print(consensus)
                continue

            if sequences:
                variants["Sequence"].append(consensus["Sequence"][0])

            nn_variants = call_variant_nn(template, consensus["Sequence"][0], consensus["Quality-Score"][0])



            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(nn_variants["Position"])
            variants["Variant"].append(nn_variants["Variant"])
            variants["Quality-Score"].append(nn_variants["Quality-Score"])
    

    variant_df = rename_barcode(pd.DataFrame(variants).merge(n_counts, on=["RBC","FBC"] , how="left"))

    if merge:
        return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    else:
        return variant_df

def get_variant_df_custom(demultiplex_folder: Path, ref_seq : Path, barcode_dicts : dict = None, merge = True, padding = 50, min_depth= 15):

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = template_df(barcode_dicts, rowwise=False)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Quality-Score": []}

    template = get_template_sequence(ref_seq) # Reference sequence

    summary = read_summary_file(demultiplex_folder)
    n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index()


    for barcode_id, barcode_dict in barcode_dicts.items():

        rbc = os.path.basename(barcode_id)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)

            count = n_counts[ (n_counts["RBC"] == rbc) & (n_counts["FBC"] == fbc)]["count"].values[0]

            if count < min_depth:
                print(f"Skipping {rbc}/{fbc} because of low read count")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(["NA"])
                variants["Variant"].append(["NA"])
                variants["Quality-Score"].append(["NA"])
                continue

            
            bam_file = front_barcode / "alignment_minimap.bam"

            print(bam_file)

            if not bam_file.exists():
                print(f"{bam_file} does not exist.")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(["NA"])
                variants["Variant"].append(["NA"])
                variants["Quality-Score"].append(["NA"])
                continue

            try:
                nn_variants = call_variant_cust(bam_file, "HetCPIII", range(padding,len(template)-padding + 1), template) # TODO read header from reference file for name
            
            except:
                print(f"Skipping {rbc}/{fbc}")
                continue

            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(nn_variants["Position"])
            variants["Variant"].append(nn_variants["Variant"])
            variants["Quality-Score"].append(nn_variants["Quality-Score"])
            
    

    variant_df = rename_barcode(pd.DataFrame(variants).merge(n_counts, on=["RBC","FBC"] , how="left"))

    if merge:
        return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    else:
        return variant_df




    
def get_variant_vcf(demultiplex_folder: Path, ref_seq : Path, barcode_dicts : dict = None, consensus_folder_name = "consensus" , sequences : bool = False) -> pd.DataFrame:
    """ Call variants from medaka variant option

    Args:
        - demultiplex_folder, folder where the demultiplexed files are located
        - ref_seq, path to reference sequence
        - barcode_dicts, dictionary with barcodes
        - sequences, if True, return the consensus sequence
        - consensus_folder_name, folder where the consensus sequences are located
        - sequences, if True, return the consensus sequence
    """

    pass

def minION_summarise(experiment_folder):
    """Summarises the the whole ev-Seq-minION"""

    minION_summary = {"N_reads_basecalled": "NA"
                      }

    # 1) minION sequencing run

    # 2) basecalling

    basecalled_folder = os.path.join(experiment_folder, "basecalled")
    file = os.path.join(basecalled_folder, "basecalled.fastq")
    n_reads_basecalled = count_reads_in_fastq(file)

    # 3) Barcode demultiplexing

    demultiplex_folder = os.path.join(experiment_folder, "demultiplex")
    
    # 4) Consensus calling

    # 5) Variant calling



    minION_summary["N_reads_basecalled"] = n_reads_basecalled


    return minION_summary


def extract_query_sequence(text):
    lines = text.split('\n')
    query_sequence = ""

    for line in lines:
        if line.startswith("query"):
            sequence = line.split()[-1]  # assuming the sequence is the last item in the line
            query_sequence += sequence

    return query_sequence




def get_highest_non_ref_base_freq(bam_file, chrom, positions, reference_sequence):

    base_frequencies = {position: Counter() for position in positions}
    base_qualities = {position: [] for position in positions}  # To store base qualities
    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, min(positions) - 1, max(positions),
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            if pileup_column.pos + 1 in positions:
                for pileup_read in pileup_column.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        base_quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                        base_frequencies[pileup_column.pos + 1].update([base])
                        base_qualities[pileup_column.pos + 1].append(base_quality)
                    
                    elif pileup_read.is_del:
                        base_frequencies[pileup_column.pos + 1].update(["-"])
                        base_qualities[pileup_column.pos + 1].append(0)


    # Calculate highest non-reference base frequency
    highest_non_ref_base_freq = {}
    mean_quality_score = {}
    for position in positions:
        counts = base_frequencies[position]
        qualities = base_qualities[position]
        ref_base = reference_sequence[position - 1].upper()  
        
        total_bases = sum(counts.values())
        if total_bases > 0:
            non_ref_bases = {base: count for base, count in counts.items() if base != ref_base}
            if non_ref_bases:
                max_base = max(non_ref_bases, key=non_ref_bases.get)
                max_freq = non_ref_bases[max_base] / total_bases
                
                # Calculate the average quality for the highest non-ref base
                max_base_qualities = [qual for base, qual in zip(counts.elements(), qualities) if base == max_base]
                avg_quality = sum(max_base_qualities) / len(max_base_qualities) if max_base_qualities else 0
                
                highest_non_ref_base_freq[position] = (max_base, max_freq)
                mean_quality_score[position] = int(avg_quality)
            else:
                highest_non_ref_base_freq[position] = (None, 0)
                mean_quality_score[position] = 0
        else:
            highest_non_ref_base_freq[position] = (None, 0)
            mean_quality_score[position] = 0

    return highest_non_ref_base_freq, mean_quality_score



def get_highest_non_ref_base_freq_2(bam_file, chrom, positions, reference_sequence, qualities=True):

    base_frequencies = {position: Counter() for position in positions}
    base_qualities = {position: [] for position in positions} if qualities else None

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, min(positions) - 1, max(positions),
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            if pileup_column.pos + 1 in positions:
                for pileup_read in pileup_column.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        base_frequencies[pileup_column.pos + 1].update([base])
                        if qualities:
                            base_quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                            base_qualities[pileup_column.pos + 1].append(base_quality)
                    elif pileup_read.is_del:
                        base_frequencies[pileup_column.pos + 1].update(["-"])
                        if qualities:
                            base_qualities[pileup_column.pos + 1].append(0)

    highest_non_ref_base_freq = {}
    mean_quality_score = {position: 0 for position in positions} if qualities else None
    for position in positions:
        counts = base_frequencies[position]
        ref_base = reference_sequence[position - 1].upper()  
        total_bases = sum(counts.values())
        if total_bases > 0:
            non_ref_bases = {base: count for base, count in counts.items() if base != ref_base}
            if non_ref_bases:
                max_base = max(non_ref_bases, key=non_ref_bases.get)
                max_freq = non_ref_bases[max_base] / total_bases
                highest_non_ref_base_freq[position] = (max_base, max_freq)

                if qualities:
                    max_base_qualities = [qual for base, qual in zip(counts.elements(), base_qualities[position]) if base == max_base]
                    mean_quality_score[position] = int(sum(max_base_qualities) / len(max_base_qualities)) if max_base_qualities else 0
            else:
                highest_non_ref_base_freq[position] = (None, 0)
        else:
            highest_non_ref_base_freq[position] = (None, 0)

    return highest_non_ref_base_freq, mean_quality_score if qualities else highest_non_ref_base_freq


def get_mean_quality_score(bam_file, chrom, positions, reference_sequence):
    """Get the mean quality score of the highest non-reference base at each position.
    Args:
        - bam_file (str): Path to the BAM file.
        - chrom (str): Chromosome or contig name.
        - positions (list): List of positions to call variants.
        - reference_sequence (str): Reference sequence of the chromosome or contig.
    Returns:
        - mean_quality_score (dict): Dictionary with positions as keys and mean quality scores as values.
    """
    base_qualities = {position: [] for position in positions}  

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, min(positions) - 1, max(positions),
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            if pileup_column.pos + 1 in positions:
                for pileup_read in pileup_column.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        base_quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                        base_frequencies[pileup_column.pos + 1].update([base])
                        base_qualities[pileup_column.pos + 1].append(base_quality)
                    
                    elif pileup_read.is_del:
                        base_frequencies[pileup_column.pos + 1].update(["NaN"])
                        base_qualities[pileup_column.pos + 1].append("NaN")

    mean_quality_score = {}

    return mean_quality_score

def call_and_filter_vcf(input_path, reference, allele_frequency):
    """ Uses bcftools to call variants and filter them based on allele frequency.
    Args:
        - input_path (str): Path to the input files.
        - reference (str): Path to the reference genome.
        - allele_frequency (float): Allele frequency threshold.
    Returns:
        - Subprocess run
    """
    
    
    prompt_pileup = f"bcftools mpileup -d 4000 -Ou -f {reference}  {input_path}/alignment.bam > {input_path}/pileup.bcf"

    prompt_call = f"bcftools call -mv -Ob -o {input_path}/raw_variants.bcf {input_path}/pileup.bcf"

    prompt_view = f"bcftools view -i 'INFO/AF>{allele_frequency}' -Ob -o {input_path}/filtered_variants.vcf {input_path}/raw_variants.bcf"


    subprocess.run(prompt_pileup, shell=True)

    subprocess.run(prompt_call, shell=True)

    subprocess.run(prompt_view, shell=True)

    print(f"Variant calling and filtering completed. Output saved to {input_path}/raw_variants_python.bcf")


def extract_positions_from_vcf(vcf_file : str) -> list:
    """ Extracts the positions of the variants from a VCF file.
    Args:
        - vcf_file (str): Path to the VCF file. 
    Returns:
        - positions (list): List of variant positions.
    """
    positions = []
 
    vcf = pysam.VariantFile(vcf_file)

    for record in vcf:
        positions.append(record.pos)
    
    vcf.close()

    return positions

def extract_mutations_from_vcf(vcf_file):
    """ Extracts the mutations from a VCF file.
    Args:
        - vcf_file (str): Path to the VCF file.
    Returns:
        - formatted_mutations (str): Formatted mutations (e.g A100T_G223C).
    """

    mutations = []

    vcf = pysam.VariantFile(vcf_file)

    for record in vcf:

        ref_allele = record.ref
        alt_alleles = record.alts

        for alt_allele in alt_alleles:
            position = record.pos
            mutation = f"{ref_allele}{position}{alt_allele}"
            mutations.append(mutation)

    vcf.close()

    formatted_mutations = "_".join(mutations)

    return formatted_mutations


def get_base_counts_at_position(bam_file, chrom, position):
    """
    Extract unique bases and gaps (deletions) and their counts from a BAM file 
    at a specific position using pileup.
    
    Args:
    - bam_file (str): path to the BAM file.
    - chrom (str): chromosome or contig name.
    - position (int): 1-based position to extract bases from.

    Returns:
    - dict: unique bases and gaps with their counts at the specified position.
    """
    
    bases = []

    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, position - 1, position, 
                                        min_base_quality=0, 
                                        min_mapping_quality=0, 
                                        truncate=True):
            if pileup_column.pos == position - 1: 
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del:
                        bases.append("-")  
                    elif not pileup_read.is_refskip:
                        bases.append(pileup_read.alignment.query_sequence[pileup_read.query_position])

    base_counts = Counter(bases)

    return base_counts


def generate_heatmap_data(bam_file, chrom, positions):
    data = []
    all_bases = set(['A', 'T', 'C', 'G', '-'])
    for position in positions:
        base_counts = get_base_counts_at_position(bam_file, chrom, position)
        col = [base_counts.get(base, 0) for base in all_bases]
        data.append(col)
    df = pd.DataFrame(data, index=positions, columns=list(all_bases))
    # Order as in all_bases
    # TODO: Specify if a character if it is an actual gap or no coverage
    df = df[['A', 'T', 'C', 'G', '-']]
    return df


def get_most_common_base(position, heatmap_data, ref_base):
    non_ref_data = heatmap_data.loc[position].drop(ref_base)
    return non_ref_data.idxmax(), non_ref_data.max()



def get_bases_from_pileup(bam_file, chrom, positions):
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
                        base = '-'  # or any symbol you prefer to represent a deletion
                    elif not pileup_read.is_refskip:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    else:
                        continue

                    # Add base to the dictionary
                    if read_name not in bases_dict[pos]:
                        bases_dict[pos][read_name] = base

    # Get unique read names and sort them
    read_names = sorted(set().union(*[bases_dict[pos].keys() for pos in bases_dict]))

    # Create DataFrame
    df = pd.DataFrame(index=read_names, columns=positions)
    
    # Populate DataFrame
    for pos in positions:
        for read_name, base in bases_dict[pos].items():
            df.at[read_name, pos] = base
    
    # Fill NaN with "-"
    df = df.fillna("-")

    return df
def get_max_base_at_position(bam_file, chrom, position):
    base_counts = get_base_counts_at_position(bam_file, chrom, position)
    return base_counts.most_common(1)[0][0]


def get_variant_name(entry, reference, nb_positions):
    variant = []
    for i, mut in enumerate(entry):
        
        ref_AA = reference[nb_positions[i] - 1]
        
        if mut == "-":
                v = f'{ref_AA}{nb_positions[i]}DEL'
        
        elif mut == ref_AA:
            continue

        else:
                v = f'{ref_AA}{nb_positions[i]}{mut}'
        variant.append(v)

    return "_".join(variant)

def get_nb_positions(freq_dist, threshold):
    nb_positions = []
    for index, row in freq_dist.iterrows():
        if row["Frequency"] > threshold:
            nb_positions.append(index)
    return nb_positions

def get_pop_frequency(bam_file, template, reference, nb_positions, min_freq=0.4, min_depth= 15):

    # Get PileUP

    # Check also for variant by sampling random positions

    bases_df = get_bases_from_pileup(bam_file, reference, nb_positions)
    
    frequency_df = bases_df.apply(get_variant_name, axis=1, args=(template, nb_positions)).value_counts().reset_index()

    frequency_df.columns = ['Population', 'N_reads']
    
    frequency_df["Frequency"] = frequency_df["N_reads"] / frequency_df["N_reads"].sum()


    # Filter for frequency > 0.4 and depth > 15
    frequency_df = frequency_df[(frequency_df["Frequency"] > min_freq) & (frequency_df["N_reads"] > min_depth)]

    return frequency_df

def call_variant_pop_frequency(bam_file, template, reference, min_freq=0.4, min_depth= 15, padding_start = 50, padding_end = 51):
    """ 
    Call Variant from BAM file based on Basecall Frequency & Alignment Frequency

    Args:
        - bam_file (str): Path to the BAM file.
        - template (str): Template name >NAME from fasta or fastq files.
        - reference (str): ID Name of reference e.g HetCPIII
        - nb_positions (list): List of positions to call variants.
        - min_freq (float): Minimum frequency of the variant.
        - min_depth (int): Minimum depth of the variant.
    
    Returns:
        - variant (str): Variant name.
    """

    variants = {"Variant" : [], "Position" : [], "Alignment Count" : [], "Alignment Frequency" : []}

    try:
        alignment_count = int(subprocess.run(f"samtools view -c {bam_file}", shell=True, capture_output=True).stdout.decode("utf-8").strip())
        range_positions = range(padding_start + 1, len(template) - padding_end + 1) 
        freq_dist = pd.DataFrame(get_highest_non_ref_base_freq_2(bam_file, reference, range_positions, template, qualities=False)[0]).T.rename(columns={0:"Base", 1:"Frequency"})

        nb_positions = get_nb_positions(freq_dist, min_freq)

        print(nb_positions)

        if nb_positions == [] and alignment_count > min_depth:
            variants["Variant"].append("#PARENT#")
            variants["Position"].append("-")
            variants["Alignment Count"].append(alignment_count)
            variants["Alignment Frequency"].append("-")
            return variants
            

        elif nb_positions == [] and alignment_count < min_depth:
            variants["Variant"].append("NA")
            variants["Position"].append("-")
            variants["Alignment Count"].append(alignment_count)
            variants["Alignment Frequency"].append("-")
            return variants

        elif len(nb_positions) == 1 and alignment_count > min_depth:

            ref_base = template[nb_positions[0] - 1]
            pos = nb_positions[0]
            base = freq_dist["Base"][nb_positions[0]]

            if base == "-":
                variants["Variant"].append(f"{ref_base}{pos}DEL")
            
            else:
                variants["Variant"].append(f"{ref_base}{pos}{base}")

            variants["Position"].append(nb_positions[0])
            variants["Alignment Count"].append(alignment_count)
            variants["Alignment Frequency"].append(freq_dist["Frequency"][nb_positions[0]])
            return variants

        # Important, after calling variant substract the position by the padding length  

        elif len(nb_positions) > 1 and alignment_count > min_depth:

            freq_df = get_pop_frequency(bam_file, template, reference, nb_positions, min_freq=min_freq, min_depth= 15)
          
            for index, row in freq_df.iterrows():
               
                variant = row["Population"]
                variants["Variant"].append(variant)
                variants["Position"].append(nb_positions)
                variants["Alignment Count"].append(row["N_reads"])
                variants["Alignment Frequency"].append(row["Frequency"])
        else:
            variants["Variant"].append("NA")
            variants["Position"].append("-")
            variants["Alignment Count"].append(alignment_count)
            variants["Alignment Frequency"].append("-")

    except:
        print("Error in calling variant")
        variants["Variant"].append("NA")
        variants["Position"].append("-")
        variants["Alignment Count"].append("-")
        variants["Alignment Frequency"].append("-")

    return variants


def format_variant_list(variant_list):
    """ Convert a list of integer variants to a string format. """
    if isinstance(variant_list, list):
        return '_'.join(str(v) for v in variant_list)
    return variant_list

def adjust_variant(variant, padding_start):
    """ 
    Adjusts the variant position based on the padding length. After variant calling, the position of the variant is based on the reference sequence.
    """

    if "#PARENT#" in variant:
        return "#PARENT#"
    
    elif "Mixed" in variant:
        return "Mixed"
    
    elif variant == "NA":
        return "NA"
    
    else:
        variants = variant.split('_')
        adjusted_variants = []

        for v in variants:
            # Find the position number using regular expression
            match = re.search(r'([A-Za-z]+)(\d+)([A-Za-z]+)', v)
            if match:
                refAA, pos, newAA = match.groups()
                
                adjusted_pos = max(int(pos) - padding_start, 1)  
                adjusted_variants.append(f"{refAA}{adjusted_pos}{newAA}")

    return '_'.join(adjusted_variants)


def run_alignment_and_indexing(ref, output_dir):
    """
    Aligns sequences using minimap2, converts to BAM, sorts and indexes the BAM file.

    Args:
    ref (str): Path to the reference file.
    fasta_file (str): Path to the FASTA file containing reads.
    output_dir (Path or str): Directory to store output files.

    Returns:
    None
    """
    output_dir = Path(output_dir)

    fastq_files = list(output_dir.glob("*.fastq"))
    if not fastq_files:
        raise FileNotFoundError("No FASTQ files found in the specified output directory.")

    fastq_files_str = " ".join(str(file) for file in fastq_files)

    minimap_cmd = f"minimap2 -ax map-ont {ref} {fastq_files_str} > {output_dir}/alignment_minimap.sam"
    subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    view_cmd = f"samtools view -bS {output_dir}/alignment_minimap.sam > {output_dir}/alignment_minimap.bam"
    subprocess.run(view_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    sort_cmd = f"samtools sort {output_dir}/alignment_minimap.bam -o {output_dir}/alignment_minimap.bam"
    subprocess.run(sort_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    index_cmd = f"samtools index {output_dir}/alignment_minimap.bam"
    subprocess.run(index_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)




def get_variant_df_AF(demultiplex_folder: Path, ref_seq : Path, ref_name : str, barcode_dicts : dict = None, merge = True, min_freq=0.4, min_depth= 15, padding=50):


    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = template_df(barcode_dicts, rowwise=False)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Alignment Count": [], "Alignment Frequency": []}

    template = get_template_sequence(ref_seq) # Reference sequence

    summary = read_summary_file(demultiplex_folder)
    n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index() 



    for barcode_id, barcode_dict in barcode_dicts.items():

        rbc = os.path.basename(barcode_id)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)
            print("Processing", rbc, fbc)
            count = n_counts[(n_counts["RBC"] == rbc) & (n_counts["FBC"] == fbc)]["count"].values[0]

            # # If alignment file exist continue
            # if not os.path.exists(os.path.join(front_barcode, "alignment_minimap.bam")):
            #     print(f"Alignment file in {front_barcode} does not exist, running alignment and indexing")
            #     run_alignment_and_indexing(ref_seq, front_barcode)
            
            # else: 
            #     print("Alignment file already exists, skipping alignment and indexing")
            run_alignment_and_indexing(ref_seq, front_barcode)

            bam_file = front_barcode / "alignment_minimap.bam"


            if not bam_file.exists():
                print(f"{bam_file} does not exist.")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(["NA"])
                variants["Variant"].append(["NA"])
                variants["Alignment Count"].append(["NA"])
                variants["Alignment Frequency"].append(["NA"])
                print(f"Variant: {fbc}/{rbc}")

                continue

            # try:
            if padding == 0:
                nn_variants = call_variant_pop_frequency(bam_file, template, ref_name, min_freq, min_depth, padding_start=0, padding_end=0)
            
            else: 
                nn_variants = call_variant_pop_frequency(bam_file, template, ref_name, min_freq, min_depth, padding_start=50, padding_end=51) # TODO read header from reference file for name

            if nn_variants["Variant"] == []:
                print("Empty variant list")
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(["NA"])
                variants["Variant"].append(["NA"])
                variants["Alignment Count"].append(["NA"])
                variants["Alignment Frequency"].append(["NA"])
                print(f"Variant: {fbc}/{rbc}")
                continue

            for i, variant in enumerate(nn_variants["Variant"]):
                print(variant)
                variants["RBC"].append(rbc)
                variants["FBC"].append(fbc)
                variants["Position"].append(nn_variants["Position"][i])
                variants["Variant"].append(variant)
                #variants["Reads"].append(count)
                # Check if Alignment count is a number
                if isinstance(nn_variants["Alignment Count"][i], int) & (isinstance(nn_variants["Alignment Frequency"][i], float) or nn_variants["Alignment Frequency"][i] == "-"):
                    variants["Alignment Count"].append(nn_variants["Alignment Count"][i])
                    variants["Alignment Frequency"].append(nn_variants["Alignment Frequency"][i])

        
                else:
                    print(f"Skipping {rbc}/{fbc} due incomplete data", nn_variants["Alignment Count"][i], nn_variants["Alignment Frequency"][i])
                    variants["Alignment Count"].append("NA")
                    variants["Alignment Frequency"].append("NA")

                print(f"Variant: {fbc}/{rbc} {variant} {nn_variants['Alignment Count'][i]} {nn_variants['Alignment Frequency'][i]}")
        
        # except Exception as e:
        #     # Append 'NA' in case of an exception
        #     print(f"Error processing {rbc}/{fbc}: {e}")
        #     variants["RBC"].append(rbc)
        #     variants["FBC"].append(fbc)
        #     variants["Position"].append("NA")
        #     variants["Variant"].append("NA")
        #     variants["Alignment Count"].append("NA")
        #     variants["Alignment Frequency"].append("NA")



    if merge:
        variant_df = rename_barcode(pd.DataFrame(variants).merge(n_counts, on=["RBC","FBC"] , how="left"))
        variant_df["Variant"] = variant_df["Variant"].apply(format_variant_list)
        variant_df["Variant"] = variant_df["Variant"].apply(lambda x: adjust_variant(x, padding))

        return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    else:
        return variants


def process_barcode(barcode_id, n_counts, ref_seq, template, ref_name, min_freq, min_depth, padding):

    fbc = os.path.basename(barcode_id)
    rbc = os.path.basename(os.path.dirname(barcode_id))


    print(f"Processing {rbc}/{fbc}")

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Alignment Count": [], "Alignment Frequency": []}
    bam_file = os.path.join(barcode_id, "alignment_minimap.bam")
    count = n_counts[(n_counts["RBC"] == rbc) & (n_counts["FBC"] == fbc)]["count"].values[0]

    if not os.path.exists(bam_file):
        print(f"Alignment file in {barcode_id} does not exist, running alignment and indexing")
        run_alignment_and_indexing(ref_seq, barcode_id)

    if not os.path.exists(bam_file):
        print(f"{bam_file} does not exist.")
        append_na_to_variants(variants, rbc, fbc)
        return variants

    try:

        nn_variants = call_variant_pop_frequency(bam_file, template, ref_name, min_freq, min_depth, padding_start = padding, padding_end = padding + 1)

        for i, variant in enumerate(nn_variants["Variant"]):
            if isinstance(nn_variants["Alignment Count"][i], int) and (isinstance(nn_variants["Alignment Frequency"][i], float) or nn_variants["Alignment Frequency"][i] == "-"):
                append_variant_info(variants, rbc, fbc, nn_variants, i)
            else:
                print(f"Skipping {rbc}/{fbc}", nn_variants["Alignment Count"][i], nn_variants["Alignment Frequency"][i])
                append_na_to_variants(variants, rbc, fbc)
    except Exception as e:
        print(f"Error processing {rbc}/{fbc}: {e}")
        append_na_to_variants(variants, rbc, fbc)

    return variants



def append_na_to_variants(variants, rbc, fbc):
    variants["RBC"].append(rbc)
    variants["FBC"].append(fbc)
    variants["Position"].append("NA")
    variants["Variant"].append("NA")
    variants["Alignment Count"].append("NA")
    variants["Alignment Frequency"].append("NA")

def append_variant_info(variants, rbc, fbc, nn_variants, index):
    variants["RBC"].append(rbc)
    variants["FBC"].append(fbc)
    variants["Position"].append(nn_variants["Position"][index])
    variants["Variant"].append(nn_variants["Variant"][index])
    variants["Alignment Count"].append(nn_variants["Alignment Count"][index])
    variants["Alignment Frequency"].append(nn_variants["Alignment Frequency"][index])



def get_variant_df_AF_parallel(demultiplex_folder, ref_seq, ref_name, barcode_dicts=None, merge=True, min_freq=0.1, min_depth=15, num_jobs=8):
    padding = 50

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)

    variant_template_df = template_df(barcode_dicts, rowwise=False)
    template = get_template_sequence(ref_seq)
    summary = read_summary_file(demultiplex_folder)
    n_counts = summary.groupby(["RBC", "FBC"])["FBC"].value_counts().reset_index(name="count")

    barcode_ids = [os.path.join(demultiplex_folder, rbc, fbc) for rbc, barcodes in barcode_dicts.items() for fbc in barcodes]

    futures = []


    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize the structure to store all variants
    all_variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Alignment Count": [], "Alignment Frequency": []}

    # Use ThreadPoolExecutor to process barcodes in parallel
    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        # Submit tasks
        futures = {executor.submit(process_barcode, barcode_id, n_counts, ref_seq, template, ref_name, min_freq, min_depth, padding): barcode_id for barcode_id in barcode_ids}

        # Process results as they complete
        for future in as_completed(futures):
            barcode_id = futures[future]
            try:
                result = future.result()
                logging.info(f"Result received for {barcode_id}")
                # Extend all_variants with results from this future
                if result:
                    for key in all_variants:
                        all_variants[key].extend(result[key])
            except Exception as e:
                logging.error(f"Error with barcode {barcode_id}: {e}")

    variant_df = pd.DataFrame(all_variants).merge(n_counts, on=["RBC", "FBC"], how="left")
    variant_df["Variant"] = variant_df["Variant"].astype(str)
    variant_df["Variant"] = variant_df["Variant"].apply(lambda x: adjust_variant(x, padding))

    if merge:
        return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    else:
        return variant_df


def get_mixed_population(variant_df):

    # Get row that are not unique
    mixed_pop = variant_df[variant_df.duplicated(subset=["Plate", "Well"], keep=False)]

    return mixed_pop


