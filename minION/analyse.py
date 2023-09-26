from minION.util.IO_processor import read_fasta_file, find_file, get_barcode_dict
from minION.util.globals import CODONS
from uncertainties import ufloat 
import pandas as pd
import os
from statistics import mean
from Bio import Align
from Bio.Seq import Seq

def get_consensus_sequence(path, score = True):
    """Get the consensus sequence with the Phred Quality Score
    Input:  - path, where the fasta file is located
            - score, if True, return sequence and quality scores
    Output: Dictionary with sequence only or sequence and quality scores"""

    if path is None:
        return {"Sequence" : "NA", "Quality-Score" : "NA"}

    sequences = read_fasta_file(path, True)

    return sequences

def get_template_sequence(path):
    """Read template sequence fasta file
    Input:  - path, where the fasta file is located
    Output: Template sequence"""
    
    template = read_fasta_file(path)
    
    return template["Sequence"][0]

def get_n_reads(path):
    """Get the number of reads from a fastq file
    Input:  - path, where the fastq file is located
    Output: Number of reads"""

    n_reads = 0

    with open(path, "r") as fastq:
        for line in fastq:
            if line.startswith("@"):
                n_reads += 1

    return n_reads


def check_sequence(sequence : str = ""):
    """Check if the sequence is a multiple of 3, so that it can be translated into a protein sequence. Otherwise the sequence is not valid. Reason can be low quality, ...
    Input:  - sequence, DNA/RNA sequence
    Output: - True if sequence is a multiple of 3, False otherwise"""
    if len(sequence) % 3 == 0:
        return True
    else:
        return False

def align_sequences(consensus_path, template_path):
    """Align consensus sequences to template sequence.
    Input:  - template_path, path to template sequence
            - consensus_path, path to consensus sequence
    Output: Dictionary with aligned sequences"""

    template = get_template_sequence(template_path)

    consensus = get_consensus_sequence(consensus_path, True)

    return bioalign(template, consensus), template


def swap_NN(template, consensus, min_score=10):
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

    # Retrieve the second sequence from the alignment
    aligned_seq1 = str(alignments).split('\n')[0] # Index 1 is the special character indicating the alignment |
    aligned_seq2 = str(alignments).split('\n')[2] # Does not work for biopython 1.81 anymore, TODO adjust for new version

    
    aligned_seq = "".join([s1 if s2 == '-' else s2 for s1, s2 in zip(aligned_seq1, aligned_seq2)])

    score_final = []
    score_index = 0
    for char in aligned_seq2:
        if char == '-':
            score_final.append(0)
        else:
            score_final.append(score[score_index])
            score_index += 1

    return {"Sequence" : [aligned_seq], "Quality-Score" : score_final}


def call_variant(template, consensus, quality_score):
    """Compare the template and consensus sequence and call variants
    Input:  
        - template, template sequence
        - consensus, consensus sequence
    Output: 
        - Dictionary with Mutaiton and at which Position"""
    
    
    if len(template) != len(consensus) or consensus == "NA" or quality_score == "NA":
        return {"Variant" : "NA", "Position" : "NA", "Quality-Score" : "NA"}

    variants = {"Variant" : [], "Position" : [], "Quality-Score" : []}
    
    for i in range(min(len(template), len(consensus))):
        if template[i] != consensus[i]:
            variants["Variant"].append(f"{template[i]}->{consensus[i]}")
            variants["Position"].append(i + 1)
            variants["Quality-Score"].append(quality_score[i])
    return variants

def read_summary_file(demultiplex_folder):
    """Read barcoding summary files """
    path = find_file(demultiplex_folder, "barcoding_summary", ".txt")

    return pd.read_csv(path, sep = "\t")


def convert_ascii_to_phred(ascii):
    """Convert ascII phred scores to numerical Phred Score (Q). To get the numerical score, the ASCII value of the character is subtracted by 33.
    Input:  - ascii, ASCII character
    Output: Numerical Phred Score"""
    return ord(ascii) - 33


def barcode_to_well(barcode):
    number = int(barcode.split('barcode')[-1])
    rows = 'ABCDEFGH'
    row = rows[(number-1) // 12]
    col = (number-1) % 12 + 1
    return f"{row}{col}"

def rename_barcode(variant_df):
    """Rename the barcode to well names"""
        # Define the shape of a 96-well plate

    if 'RBC' not in variant_df.columns or 'FBC' not in variant_df.columns:
        first_col = variant_df.columns[0]
        second_col = variant_df.columns[1]
        variant_df = variant_df.rename(columns={first_col: 'Plate', second_col: 'Well'})
    
    else:
        variant_df = variant_df.rename(columns={'RBC': 'Plate', 'FBC': 'Well'})

    variant_df["Plate"] = variant_df['Plate'].str.extract('(\d+)').astype(int)
    
    variant_df["Well"] = variant_df["Well"].apply(barcode_to_well)

    return variant_df



def barcode_arrangements(rbc):
    """
    Get the number of reads for each barcode/well.
    
    Input:
    - rbc (str): reverse barcode.
    
    Output:
    - DataFrame: Contains the number of reads for each barcode/well.
    """
    
    barcode_counts = (
        read_summary_file(rbc)["barcode_arrangement"]
        .value_counts()
        .reset_index()
        .rename(columns={"count": "N_reads"})
    )

    return barcode_counts[barcode_counts["barcode_arrangement"] != "unclassified"]

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

def template_df(barcode_dicts : dict = None):
    """To have coherent df for each experiment, a template df is created. The template also have the desired plates and columns in the desired order
    Input:
        - demultiplex_folder, folder where the demultiplexed files are located"""
    
    if barcode_dicts is None:
        raise ValueError("No barcode dictionary provided")
    
    n_rbc = len(barcode_dicts.items())

    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    columns = [i for i in range(1, 13)]

    template = {"Plate": [], "Well": []}

    for i in range(n_rbc):
        for row in rows:
            for column in columns:
                template["Plate"].append(i+1)
                template["Well"].append(f"{row}{column}")
        
    return pd.DataFrame(template)

        
def get_variant_df(demultiplex_folder, ref_seq, barcode_dicts : dict = None, sequences = False):
    """Call Variants from consensus sequences and return Variants for each barcode/well

    Input:  
        - result_folder, folder where the demultiplexed files are located
        - ref_seq, reference sequence
        - barcode_dicts, dictionary with barcodes and their corresponding sequences (Optional)

    Output: 
        - DataFrame of Variants for each reverse Barcode (RBC = Plate) and forward Barcode (FBC = Well). 
            For each Well, the variants are called and the number of reads for each barcode/well is counted"""

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variant_template_df = template_df(barcode_dicts)

    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": [], "Quality-Score": []}

    if sequences:
        variants["Sequence"] = []

    template = get_template_sequence(ref_seq)
    template_aa = translate_sequence([template])

    reads_df = pd.DataFrame(columns=["RBC", "FBC"])

    for barcode_id, barcode_dict in barcode_dicts.items():
        rbc = os.path.basename(barcode_id)

        n_counts = barcode_arrangements(barcode_id) # df with number of reads for each barcode/well

        n_counts["RBC"] = str(rbc)
        n_counts.rename(columns={"barcode_arrangement": "FBC"}, inplace=True)

        reads_df = pd.concat([reads_df, n_counts], ignore_index=True)

        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)

            fasta_file = os.path.join(front_barcode, "final_consensus.fasta")

            consensus_aa = AA_seq(fasta_file, ref_seq)

            if sequences:
                variants["Sequence"].append(consensus_aa["Sequence"][0][0])

            aa_variants = call_variant(template_aa["Protein-Sequence"][0], consensus_aa["Sequence"][0][0], consensus_aa["Quality-Score"])


            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(aa_variants["Position"])
            variants["Variant"].append(aa_variants["Variant"])
            variants["Quality-Score"].append(aa_variants["Quality-Score"])
    

    variant_df = rename_barcode(pd.DataFrame(variants).merge(reads_df, on=["RBC","FBC"] , how="left"))

    return variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")
    


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