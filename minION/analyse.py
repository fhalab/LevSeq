from minION.util.IO_processor import read_fasta_file, find_file, get_barcode_dict
from minION.util.globals import CODONS
import pandas as pd
import os
from statistics import mean
from Bio import Align, AlignIO
from Bio.Seq import Seq, translate

def minION_summarise():
    """Summarises the the whole ev-See-minION"""

    # 1) minION sequencing run

    # 2) basecalling

    # 3) Barcode demultiplexing

    # 4) Consensus calling

    # 5) Variant calling

def get_consensus_sequence(path, score = True):
    """Get the consensus sequence with the Phred Quality Score
    Input:  - path, where the fasta file is located
            - score, if True, return sequence and quality scores
    Output: Dictionary with sequence only or sequence and quality scores"""

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

def align_sequences(template_path, consensus_path):
    """Align consensus sequences to template sequence.
    Input:  - template_path, path to template sequence
            - consensus_path, path to consensus sequence
    Output: Dictionary with aligned sequences"""

    template = get_template_sequence(template_path)
    consensus = get_consensus_sequence(consensus_path, True)

    return bioalign(template, consensus)


def swap_NN(template, consensus, min_score = 30):
    """Swap Nucliotides in consensus with template if the score is below min_score"""
    seq = consensus["Sequence"][0]
    quality = consensus["Quality-Score"]

    if len(template) != len(seq):
        raise ValueError(f"'template' and 'seq' have different lengths: {len(template)} vs {len(seq)}")


    for i in range(len(template)):
        if template[i] != seq[i]:
            if quality[i] < min_score:
                seq = seq[:i] + template[i] + seq[i+1:]

    aa_scores = mean_quality_score(quality)

    return {"Sequence" : seq, "Quality-Score" : aa_scores}

def mean_quality_score(quality_scores):
    """Calculate the mean quality score of a sequence
    Input:  - quality_scores, List of Ascii Phred Scores
    Output: - Mean quality score for each codon"""
    aa_scores = []
    n_quality_scores = len(quality_scores)
    i = 1

    while i <= n_quality_scores:
        if i % 3 == 0:
            mean_score = round(mean(quality_scores[i-3:i]),3)
            aa_scores.append(mean_score)
        i += 1



    return aa_scores


def translate_sequence(sequences=None):
    """
    Translate a list of sequences into protein sequences.
    Input:  - sequences, a list of DNA sequences
    Output: Dictionary containing DNA sequences and their corresponding protein sequences
    """

    proteins = {"DNA-Sequence": [],
                "Protein-Sequence": []}

    if sequences is None:
        return proteins

    for sequence in sequences:
        proteins["DNA-Sequence"].append(sequence)
        
        if check_sequence(sequence):
            proteins["Protein-Sequence"].append("".join([CODONS[sequence[i:i+3]] for i in range(0, len(sequence), 3) if i+3 <= len(sequence)]))
        else:
            proteins["Protein-Sequence"].append("NA")

    return proteins


def bioalign(template, consensus):
    """Align the template and consensus sequence"""


    seq1 = Seq(template)
    seq2 = Seq(consensus["Sequence"][0])
    score = consensus["Quality-Score"][0]
    aligner = Align.PairwiseAligner()
    
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -10
    aligner.target_end_gap_score = -100
    aligner.mismatch_score = 0
    
    alignments = aligner.align(seq1,seq2)[0]


    # Retrieve the second sequence from the alignment
    aligned_seq1 = str(alignments).split('\n')[0] # Index 1 is the special character indicating the alignment |
    aligned_seq2 = str(alignments).split('\n')[2]

    aligned_seq = "".join([seq1[i] if char == '-' else char for i, char in enumerate(aligned_seq2)])

    

    inds = aligned_seq.find('-')
    if inds != -1:
        score_final = score[:inds] + '!' + score[inds + 1:]
    else:
        score_final = score

    return {"Sequence" : [aligned_seq], "Quality-Score" : score_final}


def get_AA_sequence(template_fasta, consensus_path):
    """Get the AA sequence from the consensus sequence
    Input:  - template_fasta, path to template fasta file
            - consensus_path, path to consensus sequence
    Output: Dictionary with AA sequence and quality scores"""


    aligned_seq = align_sequences(template_fasta, consensus_path)

    return translate_sequence([aligned_seq["Sequence"][0]])

def call_variant(template, consensus):
    """Compare the template and consensus sequence and call variants
    Input:  - template, template sequence
            - consensus, consensus sequence
    Output: Dictionary with Mutaiton and at which Position"""
    
    if len(template) != len(consensus):
        return {"Variant" : "NA", "Position" : "NA"}

    variants = {"Variant" : [], "Position" : []}
    
   
    for i in range(min(len(template), len(consensus))):
        if template[i] != consensus[i]:
            variants["Variant"].append(f"{template[i]}->{consensus[i]}")
            variants["Position"].append(i + 1)
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

def get_variant_df(demultiplex_folder, ref_seq, barcode_dicts : dict = None):
    """Call Variants from consensus sequences and return Variants for each barcode/well
    Input:  - result_folder, folder where the demultiplexed files are located
            - ref_seq, reference sequence
            - barcode_dicts, dictionary with barcodes and their corresponding sequences (Optional)
    Output: Dictionary with variants for each barcode/well"""

    if barcode_dicts is None:
        barcode_dicts = get_barcode_dict(demultiplex_folder)
    
    variants = {"RBC": [], "FBC": [], "Position": [], "Variant": []}

    for barcode_id, barcode_dict in barcode_dicts.items():
        rbc = os.path.basename(barcode_id)
        for front_barcode in barcode_dict:

            fbc = os.path.basename(front_barcode)

            fasta_file = os.path.join(front_barcode, "final_consensus.fasta")

            sequences = get_consensus_sequence(fasta_file)

            template = get_template_sequence(ref_seq)
            consensus = sequences["Sequence"][0]
            score = sequences["Quality-Score"][0]

            aligned_seq = bioalign(template, consensus, score)

            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(call_variant(aligned_seq["Sequence"][0], template)["Position"])
            variants["Variant"].append(call_variant(aligned_seq["Sequence"][0], template)["Variant"])

    return variants


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

def annotate_well(variant_df, control):
    """Annotate wells for plates with control wells"""
    
    

