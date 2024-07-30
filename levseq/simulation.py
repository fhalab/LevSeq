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
import random

from levseq.variantcaller import *
import math


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
    df["True Variant"] = ''
    df.set_index('index', inplace=True)
    return df


def mutate_sequence(sequence, mutation_frequency, bases=None):
    """
    Mutates a given nucleotide sequence at a specified mutation frequency.
    """
    bases = bases if bases is not None else ['A', 'T', 'G', 'C', '-']  # Inlucde deletions

    sequence_list = list(sequence)

    # Iterate over the sequence and mutate bases with the given probability
    for i, base in enumerate(sequence_list):
        if random.random() < mutation_frequency:
            # Choose a new base different from the current one (possibly extend this to have a preference?)
            new_base = random.choice([b for b in bases if b != base])
            sequence_list[i] = new_base

    # Convert the list back to a string
    mutated_sequence = ''.join(sequence_list)
    return mutated_sequence


def insert_nt(original_nt_seq, protein_mutations, codon_usage):
    nt_seq_list = list(original_nt_seq)
    for pos, new_aa in protein_mutations:
        # Convert protein position to nucleotide position
        nt_pos = pos * 3  # Assuming the codon starts at this position

        # Select a codon for the new amino acid
        new_codon = codon_usage[new_aa]

        # Replace the original codon in the nucleotide sequence
        nt_seq_list[nt_pos:nt_pos + 3] = list(new_codon)

    # Convert the list back to a string
    mutated_nt_seq = ''.join(nt_seq_list)
    return mutated_nt_seq


def generate_ssm_library(positions, parent_sequence_aa, parent_sequence_nt, codon_usage):
    """
    For each position, generate a SSM library for a given parent and a set of positions.
    """
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
    library = []
    nt_seq_list = list(parent_sequence_nt)
    for position in positions:
        for aa in amino_acids:
            if parent_sequence_aa[position] != aa:  # i.e. don't dup the original
                nt_pos = position * 3  # Assuming the codon starts at this position

                # Select a codon for the new amino acid
                new_codon = codon_usage[aa]
                nt_seq_list[nt_pos:nt_pos + 3] = list(new_codon)

                # Convert the list back to a string
                mutated_nt_seq = ''.join(nt_seq_list)
                library.append(mutated_nt_seq)

    return library


def make_experiment(run_label, read_depth, sequencing_error_rate, parent_sequence, library_number,
                    number_of_wells, epcr_mutation_rate, frequency_cutoff=0.5, number_wells_to_mix=0,
                    mixture_rate=0, qc_files_path=None):
    # Make a full experiment setup
    mixed_wells = None
    if number_wells_to_mix > 0:
        # mixed_wells tells us which wells are truely mixed
        mutated_sequence, mixed_wells = make_mixed_well_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number,
                                                              epcr_mutation_rate, number_wells_to_mix, mixture_rate)
    else:
        mutated_sequence = make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number,
                                                   epcr_mutation_rate)

    variant_df = get_dummy_plate_df(run_label, 'Well', number_of_wells)
    mutant_to_well_df = {}
    current_well = 0
    variant_df['True Mixed Well'] = False
    for mutant in tqdm(mutated_sequence):
        parent_name = 'Parent'
        reads = []
        read_ids = []
        quals = []
        for i, seq in enumerate(mutated_sequence[mutant]):
            read_ids.append(f'read_{i}')
            reads.append(seq)
            quals.append(100)  # Dummy don't need

        well_df = make_well_df_from_reads(reads, read_ids, quals)
        rows_all = make_row_from_read_pileup_across_well(well_df, parent_sequence, parent_name)
        well_df = pd.DataFrame(rows_all)
        well_df.columns = ['gene_name', 'position', 'ref', 'most_frequent', 'freq_non_ref', 'total_other',
                           'total_reads', 'p_value', 'percent_most_freq_mutation', 'A', 'p(a)', 'T', 'p(t)', 'G',
                           'p(g)', 'C', 'p(c)', 'N', 'p(n)']
        well_df = calculate_mutation_significance_across_well(well_df)
        if qc_files_path is not None:
            # Save QC data
            qc_well_df = make_well_df_for_saving(reads, read_ids, quals)
            write_msa_for_df(qc_well_df, well_df, parent_name, parent_sequence,
                             os.path.join(qc_files_path, f'{run_label}_{current_well}.fa'),
                             os.path.join(qc_files_path, f'{run_label}_{current_well}.csv'))
        label, frequency, combined_p_value, mixed_well = get_variant_label_for_well(well_df, frequency_cutoff)
        mutant_to_well_df[f'{mutant}_{current_well}'] = well_df
        variant_df.at[current_well, "Mixed Well"] = mixed_well
        variant_df.at[current_well, "Variant"] = label
        variant_df.at[current_well, "True Variant"] = mutant
        variant_df.at[current_well, "frequency"] = frequency
        variant_df.at[current_well, "P value"] = combined_p_value
        variant_df.at[current_well, "Well"] = f'Well {current_well}'
        variant_df.at[current_well, "Alignment_count"] = read_depth
        if mixed_wells is not None:
            variant_df.at[current_well, "True Mixed Well"] = mixed_wells[mutant] # Save this as a true mixed well
        current_well += 1

    # Before returning adjust the pvalues
    variant_df['P adj.'] = multipletests(list(variant_df["P value"].values), alpha=0.05, method='fdr_bh')[1]
    # Also get the accuracy
    variant_df = check_variants(variant_df, parent_sequence)
    return variant_df


def make_well_df_for_saving(seqs, read_ids, read_quals):
    """
    Make a dataframe in a specific format taking the reads and read IDs and filtering duplicates based on the
    read quality. Keeps the highest quality scoring read for a given read ID.
    """
    seq_df = pd.DataFrame([list(s) for s in seqs]) # Convert each string to a list so that we get positions nicely
    # Also add in the read_ids and sort by the quality to only take the highest quality one
    seq_df['read_id'] = read_ids
    seq_df['read_qual'] = read_quals
    seq_df['seqs'] = seqs
    seq_df = seq_df.sort_values(by='read_qual', ascending=False)
    # Should now be sorted by the highest quality
    seq_df = seq_df.drop_duplicates(subset=['read_id'], keep='first')
    return seq_df


def write_msa_for_df(reads_across_well_df, well_df, parent_name, parent_sequence, msa_path, df_path):
    """ This is for checking that we have the correct data."""
    read_ids = reads_across_well_df['read_id']
    seqs = reads_across_well_df['seqs']
    # Check if we want to write a MSA
    if msa_path is not None:
        with open(msa_path, 'w+') as fout:
            # Write the reference first
            fout.write(f'>{parent_name}\n{parent_sequence}\n')

            for i, seq in enumerate(seqs):
                fout.write(f'>{read_ids[i]}\n{"".join(seq)}\n')
    well_df.to_csv(df_path)


def generate_epcr_library(parent_sequence, mutation_rate, library_number):
    """
    For a parent make a number of sequenes using the error prone PCR.
    """
    return [mutate_sequence(parent_sequence, mutation_rate, ['A', 'T', 'G', 'C']) for c in range(0, library_number)]


def make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number, epcr_mutation_rate):
    """
    library_number would normally be the number of wells for example.
    """
    # Simulate the mutation frequncey
    # First make the library
    library = generate_epcr_library(parent_sequence, epcr_mutation_rate, library_number)
    # For this library for each one, simulate the number of reads with a sequencing error rate
    reads_per_well = {}
    for seq in library:
        reads_per_well[seq] = [mutate_sequence(seq, sequencing_error_rate) for i in range(0, read_depth)]
    return reads_per_well


def make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions, parent_sequence_aa,
                           codon_usage):
    """
    library_number would normally be the number of wells for example.
    """
    # Simulate the mutation freq
    library = generate_ssm_library(positions, parent_sequence_aa, parent_sequence, codon_usage)
    # For this library for each one, simulate the number of reads with a sequencing error rate
    reads_per_well = {}
    for seq in library:
        # i.e. the key is the seq and the value is the mutated reads
        reads_per_well[seq] = [mutate_sequence(seq, sequencing_error_rate) for i in range(0, read_depth)]
    return reads_per_well


def make_mixed_well_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number,
                                       epcr_mutation_rate, number_wells_to_mix, mixture_rate):
    """
    Make a mixed well experiment to test code with.
    """
    # Simulate the mutation frequncey
    # First make the library
    library = generate_epcr_library(parent_sequence, epcr_mutation_rate, library_number)
    # For this library for each one, simulate the number of reads with a sequencing error rate
    reads_per_well = {}
    reads_mutated_label = {}
    for seq in library:
        # Randomly mix some of the wells at the mixture rate. Here we'll just randomly "dope" in some of the randomly
        reads_per_well[seq] = [mutate_sequence(seq, sequencing_error_rate) for i in range(0, read_depth)]
        reads_mutated_label[seq] = False
    wells_to_mix = random.sample(list(reads_per_well.keys()), number_wells_to_mix)
    # Combine them
    for well_seq in wells_to_mix:
        # For each well, randomly select one from the other wells
        dope_in_seq = random.sample(wells_to_mix, 1)[0] # We only want one!
        if dope_in_seq != well_seq:  # Make sure the sequences are different
            # Swap out a percentage of the wells from dope in into the other well
            number_to_add_in = math.floor(read_depth*mixture_rate)
            for read_position in range(0, number_to_add_in):
                # Just make the top X the other doped in seq
                reads_per_well[well_seq][read_position] = reads_per_well[dope_in_seq][read_position]
                reads_mutated_label[well_seq] = True

    return reads_per_well, reads_mutated_label


def check_variants(variant_df, parent_sequence):
    """ This just checks if the variants are actually correct! """

    corrects = []
    incorrects = []
    true_positives = 0
    true_negatives = 0
    false_negatives = 0
    false_positives = 0
    for predicted_variant, true_variant in variant_df[['Variant', 'True Variant']].values:
        count_correct = 0
        count_incorrect = 0
        for mutation in predicted_variant.split('_'):
            try:
                if 'PARENT' in mutation:
                    # Check that the two seqeunces are correct
                    for i in range(0, len(true_variant)):
                        if true_variant[i] == parent_sequence[i]:
                            count_correct += 1
                            true_positives += 1
                        else:
                            count_incorrect += 1
                            false_negatives += 1
                else:
                    # true_variant is a sequence while predicated variant is just the mutations
                    if 'DEL' not in mutation:
                        mut_pos = int(mutation[1:-1])  # A1T
                        mut = mutation[-1]
                    else:
                        mut_pos = int(mutation[1:].replace('DEL', ''))
                        mut = 'DEL'
                    if true_variant[mut_pos - 1] == mut:
                        count_correct += 1
                        true_positives += 1
                    else:
                        count_incorrect += 1
                    try:
                        if parent_sequence[mut_pos - 1] != mutation[0]:
                            print("WARNING!")
                    except:
                        print(mut_pos, len(parent_sequence))
            except:
                print(mutation)
        corrects.append(count_correct)
        incorrects.append(count_incorrect)
    variant_df['correct'] = corrects
    variant_df['incorrect'] = incorrects
    variant_df['accuracy'] = np.array(corrects) / (np.array(corrects) + np.array(incorrects))
    return variant_df
