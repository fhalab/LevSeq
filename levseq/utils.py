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
# Import all packages
import pandas as pd
import pysam
import random
import os
import numpy as np
from copy import deepcopy
from collections import defaultdict
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from pathlib import Path
from scipy.stats import combine_pvalues
from Bio import SeqIO
from Bio.PDB.Polypeptide import aa1


amino_acid_to_codon = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
    'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    '*': 'TAA'
}


ALL_AAS = deepcopy(list(aa1))

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


# Get mutated sequence by splitting template sequence
def get_mut(temp_seq, aa_seq):
    i = 0
    mut_ls = []
    for i in range(len(list(zip(temp_seq, aa_seq)))):
        if temp_seq[i] != aa_seq[i]:
            mut_ls.append(temp_seq[i]+str(i+1)+aa_seq[i])
        i = i + 1
    return mut_ls


def check_demultiplexing(demultiplex_folder: Path, reverse_prefix="RB", forward_prefix="NB", verbose=True):
    """
    Check if the demultiplexing was done correctly. If not, return the user that the sequences were not demultiplexed.

    Args:
        - demultiplex_folder: Path to the folder containing the demultiplexed fastq files
        - verbose: If True, print the name of each parent folder and the count of child folders

    Return:
        - Tuple: Number of parent folders and child folders
    """
    demultiplex_path = Path(demultiplex_folder)
    parent_folder_count = 0
    child_folder_count = 0

    for child in demultiplex_path.iterdir():
        if child.is_dir() and (child.name.startswith(reverse_prefix) or child.name.startswith(forward_prefix)):
            parent_folder_count += 1
            child_folder_count += len(list(child.iterdir()))
            if verbose:
                print(f"Parent folder '{child.name}' contains {len(list(child.iterdir()))} folders.")

    return parent_folder_count, child_folder_count


def get_template_df(plate_numbers: list, barcode_dicts: dict = None, rowwise=True):
    """
    To have coherent df for each experiment, a template df is created. The template also have the desired plates and columns in the desired order
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
                    template["Plate"].append(plate_numbers[i])
                    template["Well"].append(f"{row}{column}")

    return pd.DataFrame(template)


def get_dummy_plate_df(plate_name='Plate', well_name='Well', number_of_wells=96):
    """
    Make a dummy plate.
    Plate	Well	Path	Alignment_count	P value	Mixed Well	Variant	Average mutation frequency	P adj. value
    """
    df = pd.DataFrame([i for i in range(0, number_of_wells)], columns=['index'])
    df['Plate'] = plate_name
    df['Well'] = well_name
    df['Path'] = ''
    df['Alignment_count'] = 0
    df['P value'] = 1.0
    df['Mixed Well'] = False
    df['Variant'] = ''
    df['mutation'] = ''
    df['frequency'] = 0
    df['P adj.'] = 0
    df['value'] = 0
    df.set_index('index', inplace=True)
    return df


def make_well_df_from_reads(seqs, read_ids, read_quals):
    """
    Make a dataframe in a specific format taking the reads and read IDs and filtering duplicates based on the
    read quality. Keeps the highest quality scoring read for a given read ID.
    """
    seq_df = pd.DataFrame([list(s) for s in seqs])  # Convert each string to a list so that we get positions nicely
    # Also add in the read_ids and sort by the quality to only take the highest quality one
    seq_df['read_id'] = read_ids
    seq_df['seqs'] = seqs
    seq_df['read_qual'] = [0 if isinstance(r, str) else r for r in read_quals]
    seq_df = seq_df.sort_values(by='read_qual', ascending=False)
    # Should now be sorted by the highest quality
    seq_df = seq_df.drop_duplicates(subset=['read_id'], keep='first')
    return seq_df.drop(columns=['read_id', 'seqs', 'read_qual'])


def calculate_mutation_significance_across_well(seq_df):
    """
    Calculate the background error as just the mean frequency of non-reference sequneces (here we "smooth"
    out the induced mutations.
    """
    mean_error = np.mean(seq_df['freq_non_ref'].values)
    seq_df.reset_index(inplace=True)
    i = 0
    if mean_error > 0.4:
        print('-----------------------------------------')
        print("WARNING!!! Your mean error rate across was too high!!! It was: ", mean_error)
        print('-----------------------------------------')

    # Using this we can calculate the significance of the different errors
    for ref_seq, num_a, num_t, num_g, num_c, num_dels, num_insertions, num_reads, num_total_non_ref_reads in seq_df[
        ['ref', 'A', 'T', 'G', 'C', 'N', 'I', 'total_reads', 'total_other']].values:
        actual_seq, val, p_value, p_a, p_t, p_g, p_c, p_n, p_i = calc_mutation_significance_for_position_in_well(ref_seq, num_a,
                                                                                             num_t, num_g, num_c,
                                                                                             num_dels, num_insertions, num_reads,
                                                                                             num_total_non_ref_reads,
                                                                                             mean_error)
        seq_df.at[i, 'p(a)'] = p_a
        seq_df.at[i, 'p(t)'] = p_t
        seq_df.at[i, 'p(g)'] = p_g
        seq_df.at[i, 'p(c)'] = p_c
        seq_df.at[i, 'p(n)'] = p_n
        seq_df.at[i, 'p(i)'] = p_n
        seq_df.at[i, 'p_value'] = p_value
        seq_df.at[i, 'percent_most_freq_mutation'] = val
        seq_df.at[i, 'most_frequent'] = actual_seq
        seq_df.at[i, 'mean_mutation_error_rate'] = mean_error  # This value should match your ePCR results
        i += 1

    # Do multiple test correction to correct each of the pvalues
    for p in ['p_value', 'p(a)', 'p(t)', 'p(g)', 'p(c)', 'p(n)', 'p(i)']:
        # Do B.H which is the simplest possibly change to have alpha be a variable! ToDo :D
        padjs = multipletests(seq_df[p].values, alpha=0.05, method='fdr_bh')
        seq_df[f'{p} adj.'] = padjs[1]
    return seq_df

def alignment_from_cigar(cigar: str, alignment: str, ref: str, query_qualities: list):
    """
    Generate the alignment from the cigar string.
    Operation	Description	Consumes query	Consumes reference
    0 M	alignment match (can be a sequence match or mismatch)	yes	yes
    1 I	insertion to the reference	yes	no
    2 D	deletion from the reference	no	yes
    3 N	skipped region from the reference	no	yes
    4 S	soft clipping (clipped sequences present in SEQ)	yes	no
    5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
    6 P	padding (silent deletion from padded reference)	no	no
    7 =	sequence match	yes	yes
    8 X	sequence mismatch	yes	yes
    """
    new_seq = ''
    ref_seq = ''
    qual = []
    inserts = {}
    pos = 0
    ref_pos = 0
    for op, op_len in cigar:
        if op == 0:  # alignment match (can be a sequence match or mismatch)
            new_seq += alignment[pos:pos + op_len]
            qual += query_qualities[pos:pos + op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
        elif op == 1:  # insertion to the reference
            inserts[pos] = alignment[pos - 1:pos + op_len]
            # new_seq += alignment[pos:pos + op_len]
            pos += op_len
        elif op == 2:  # deletion from the reference
            new_seq += '-' * op_len
            qual += [-1] * op_len
            ref_seq += ref[ref_pos:ref_pos + op_len]
            ref_pos += op_len
        elif op == 3:  # skipped region from the reference
            new_seq += '*' * op_len
            qual += [-2] * op_len
            ref_pos += op_len
        elif op == 4:  # soft clipping (clipped sequences present in SEQ)
            #inserts[pos] = alignment[pos:pos + op_len]
            pos += op_len
        elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
            continue
        elif op == 6:  # padding (silent deletion from padded reference)
            continue
        elif op == 7:  # sequence mismatch
            new_seq += alignment[pos:pos + op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            qual += query_qualities[pos:pos + op_len]
            pos += op_len
            ref_pos += op_len
    return new_seq, ref_seq, qual, inserts

def get_reads_for_well(parent_name, bam_file_path: str, ref_str: str, msa_path=None):
    """
    Rows are the reads, columns are the columns in the reference. Insertions are ignored.
    """
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    # Ensure the BAM file is indexed
    if not os.path.exists(bam_file_path + ".bai"):
        pysam.index(bam_file_path)

    seqs = []
    read_ids = []
    read_quals = []
    insert_map = defaultdict(list)
    for read in bam.fetch(until_eof=True):
        # Ensure we have at least 75% coverage
        if read.query_sequence is not None and len(read.query_sequence) > 0.75 * len(
                ref_str) and read.cigartuples is not None:
            seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                       read.query_qualities)
            # Make it totally align
            seq = "-" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
            seqs.append(seq)
            for i, insert in ins.items():
                insert_map[i].append(insert)
            read_ids.append(f'{read.query_name}')
            read_quals.append(qual)

    # Check if we want to write a MSA
    if msa_path is not None:
        with open(msa_path, 'w+') as fout:
            # Write the reference first
            fout.write(f'>{parent_name}\n{ref_str}\n')
            for i, seq in enumerate(seqs):
                fout.write(f'>{read_ids[i]}\n{str(seq)}\n')

    # Do this for all wells
    seq_df = make_well_df_from_reads(seqs, read_ids, read_quals)
    alignment_count = len(seq_df.values)
    rows_all = make_row_from_read_pileup_across_well(seq_df, ref_str, parent_name, insert_map)
    bam.close()

    if len(rows_all) > 2:  # Check if we have anything to return
        seq_df = pd.DataFrame(rows_all)
        seq_df.columns = ['gene_name', 'position', 'ref', 'most_frequent', 'freq_non_ref', 'total_other',
                          'total_reads', 'p_value', 'percent_most_freq_mutation', 'A', 'p(a)', 'T', 'p(t)', 'G', 'p(g)',
                          'C', 'p(c)', 'N', 'p(n)', 'I', 'p(i)', 'Warnings']
        return calculate_mutation_significance_across_well(seq_df), alignment_count

def make_row_from_read_pileup_across_well(well_df, ref_str, label, insert_map):
    """
    Given a pileup of reads, we want to get some summary information about that sequence
    """
    rows = []
    for col in well_df:
        vc = well_df[col].values
        ref_seq = ref_str[col]  # Keep track of the reference
        total_reads = len(vc)  # Check if there are at least 25% with a different value compared to the reference.
        total_other = len(vc[vc != ref_seq])
        freq_non_ref = total_other / total_reads
        actual_seq = ref_seq

        # Dummy values that will be filled in later once we calculate the background error rate
        warning = ''
        if total_reads < 20:
            warning = f'WARNING: you had: {total_reads}, we recommend looking at the BAM file or using a second sequencing method on this well.'
        # Check if there was an insert
        if insert_map.get(col) and len(insert_map[col][0]) > total_reads/2:  # i.e. at least half have the insert
            if warning:
                warning += '\nINSERT'
            else:
                warning = f'WARNING: INSERT.'
            rows.append([label, col, ref_seq, actual_seq, freq_non_ref, total_other, total_reads, 1.0, 0.0,
                         len(vc[vc == 'A']), 1.0, len(vc[vc == 'T']), 1.0, len(vc[vc == 'G']), 1.0,
                         len(vc[vc == 'C']), 1.0, len(vc[vc == '-']), 1.0, len(insert_map.get(col)),
                         1.0, warning])
        if ref_seq != '-':
            rows.append([label, col, ref_seq, actual_seq, freq_non_ref, total_other, total_reads, 1.0, 0.0,
                         len(vc[vc == 'A']), 1.0, len(vc[vc == 'T']), 1.0, len(vc[vc == 'G']), 1.0,
                         len(vc[vc == 'C']), 1.0, len(vc[vc == '-']), 1.0, 0,
                         1.0, warning])
    return rows


def calc_mutation_significance_for_position_in_well(ref_seq, num_a, num_t, num_g, num_c, num_dels, num_insertions,
                                                    num_reads, num_total_non_ref_reads, background_error_rate):
    """
    Use the binomial test to check if we have a significant result for any of the observed reads.
    """
    p_a = binomtest(num_a, num_reads, background_error_rate, 'greater').pvalue
    p_t = binomtest(num_t, num_reads, background_error_rate, 'greater').pvalue
    p_g = binomtest(num_g, num_reads, background_error_rate, 'greater').pvalue
    p_c = binomtest(num_c, num_reads, background_error_rate, 'greater').pvalue
    p_n = binomtest(num_dels, num_reads, background_error_rate, 'greater').pvalue
    p_i = binomtest(num_insertions, num_reads, background_error_rate, 'greater').pvalue
    val = 0
    actual_seq = ref_seq
    p_value = float('nan')  # Could also use 0 not sure what is optimal here!
    if num_total_non_ref_reads == 0:
        val = 0.0  # i.e. they were 100% the reference
        p_value = 0.0  # i.e. they are all this
    else:
        if num_a > 0 and 'A' != ref_seq and num_a / num_reads > val:
            val = num_a / num_reads
            actual_seq = 'A'
            p_value = p_a
        if num_t > 0 and 'T' != ref_seq and num_t / num_reads > val:
            val = num_t / num_reads
            actual_seq = 'T'
            p_value = p_t
        if num_g > 0 and 'G' != ref_seq and num_g / num_reads > val:
            val = num_g / num_reads
            actual_seq = 'G'
            p_value = p_g
        if num_c > 0 and 'C' != ref_seq and num_c / num_reads > val:
            val = num_c / num_reads
            actual_seq = 'C'
            p_value = p_c
        if num_dels > 0 and '-' != ref_seq and num_dels / num_reads > val:
            val = num_dels / num_reads
            actual_seq = 'DEL'
            p_value = p_n
        if num_insertions > 0 and num_insertions / num_reads > val:
            val = num_insertions / num_reads
            actual_seq = 'INS'
            p_value = p_i
    return actual_seq, val, p_value, p_a, p_t, p_g, p_c, p_n, p_i


def postprocess_variant_df(df, cutoff=5, output_path=None):
    """
    Postprocess the variant DF to check for any positions that appear to have a higher than expected
    difference to the parent or that occur too many times.
    """
    mutation_map = defaultdict(lambda: defaultdict(int))
    all_plates = pd.DataFrame()
    for plate in set(df['Plate'].values):
        plate_df = df[df['Plate'] == plate]
        positions = []
        for m in plate_df['Variant'].values:
            if str(m) != 'nan':
                m = m.split('_')
                for mutation in m:
                    if 'DEL' not in mutation and 'INS' not in mutation:
                        position = mutation[1:-1]  # i.e. trim off what it was
                        # Get the position and also keep what it was mutated to
                        mutation_map[position][mutation[-1]] += 1  # get what it was mutated too
                    elif 'DEL' in mutation:
                        position = mutation[1:].replace('DEL', '')  # i.e. trim off what it was
                        # Get the position and also keep what it was mutated to
                        mutation_map[position]['DEL'] += 1  # get what it was mutated too
                    elif 'INS' in mutation:
                        position = mutation[1:].replace('INS', '')  # i.e. trim off what it was
                        # Get the position and also keep what it was mutated to
                        mutation_map[position]['INS'] += 1  # get what it was mutated too
                    positions.append(position)
        # Make into a DF that has positions and then the number and types of muattions
        positions = list(set(positions))
        positions.sort()
        rows = []
        for position in positions:
            pos = mutation_map[position]
            a_ = pos['A'] or 0
            t_ = pos['T'] or 0
            c_ = pos['C'] or 0
            g_ = pos['G'] or 0
            del_ = pos['DEL'] or 0
            ins_ = pos['INS'] or 0
            total = a_ + t_ + g_ + c_ + del_ + ins_
            rows.append([position, total, a_, t_, g_, c_, del_, ins_])
            # CHeck if the well has a problem at a specific position
            if total > cutoff:
                print(f"Warning! Position {position} in plate {plate} was mutated: {total} times. "
                      f"This may be an error with your parent.")
        p_df = pd.DataFrame(rows, columns=['Position', 'Total wells mutated in', 'A', 'T', 'G', 'C', 'DEL', 'INS'])
        if output_path:
            # Save QC file if the user specifies a path.
            p_df.to_csv(f'{output_path}Plate_{plate}_QC.csv', index=False)
        p_df['Plate'] = plate
        all_plates = pd.concat([all_plates, p_df])
    all_plates = all_plates.sort_values(by='Total wells mutated in', ascending=False)
    return all_plates  # this gives an idea about what might be mutated.


def get_variant_label_for_well(seq_df, threshold):
    """
    Classify/label the variants and identify whether there is a mixed well at position i.
    """
    # Now use the filter for wells which have a certain threshold of non-reference mutations
    # Filter based on significance to determine whether there is a
    non_refs = seq_df[seq_df['freq_non_ref'] > threshold].sort_values(by='position')
    mixed_well = False
    # Have section for inserts to check if they are > 50% of the reads
    if seq_df['p(i) adj.'].min() < 0.05 and seq_df['I'].max() > len(seq_df) / 2:
        label = '+'
        probability = np.mean([1 - x for x in non_refs['freq_non_ref'].values])
        combined_p_value = float("nan")

    elif len(non_refs) > 0:
        positions = non_refs['position'].values
        refs = non_refs['ref'].values
        label = [f'{refs[i]}{positions[i] + 1}{actual}' for i, actual in enumerate(non_refs['most_frequent'].values)]
        # Check if it is a mixed well i.e. there were multiple with significant greater than 0.05
        padj_vals = non_refs[['p(a) adj.', 'p(t) adj.', 'p(g) adj.', 'p(c) adj.', 'p(n) adj.', 'p(i) adj.']].values

        for p in padj_vals:
            c_sig = 0
            for padj in p:
                if padj < 0.05:  # Have this as a variable
                    c_sig += 1
            if c_sig > 1:  # potential mixed well
                mixed_well = True
        label = '_'.join(label)
        # Only keep the frequency of the most frequent mutation
        probability = np.mean([x for x in non_refs['percent_most_freq_mutation'].values])
        # Combine the values
        chi2_statistic, combined_p_value = combine_pvalues([x for x in non_refs['p_value adj.'].values], method='fisher')
    else:
        label = '#PARENT#'
        probability = np.mean([1 - x for x in non_refs['freq_non_ref'].values])
        combined_p_value = float("nan")
    # Return also the mean mutation rate for the well
    mean_mutation_rate = np.mean([1 - x for x in non_refs['freq_non_ref'].values])
    return label, probability, combined_p_value, mixed_well, mean_mutation_rate