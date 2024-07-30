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

import shutil
import tempfile
import unittest
from sciutil import *

u = SciUtil()

from levseq.variantcaller import *
from levseq.simulation import *


class TestClass(unittest.TestCase):

    @classmethod
    def setup_class(self):
        local = True
        # Create a base object since it will be the same for all the tests
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))

        self.data_dir = os.path.join(THIS_DIR, 'test_data/')
        if local:
            self.tmp_dir = os.path.join(THIS_DIR, 'test_data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='test_data')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestVariantCalling(TestClass):

    def test_making_ssm(self):
        u.dp(["Testing SSM"])

        parent_sequence = "ATGAGT"
        parent_sequence_aa = 'MS'
        sequencing_error_rate = 0  # Use a 0 error rate so that we don't get confused to start off with
        read_depth = 5
        positions = [0]
        mutated_sequence = make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions,
                                                  parent_sequence_aa, amino_acid_to_codon)
        # Check that the key for a change from M to T is in there i.e. ATG in pos[0] goes to ACT
        assert 'ACTAGT' in mutated_sequence
        # Check that M to F is in there
        assert 'TTTAGT' in mutated_sequence
        # Make sure nothing not supposed to be in there is not in tehre
        assert 'TTTTTT' not in mutated_sequence

        # Check since we have a 0 seq error rate that we don't have any errors so they are all the expected one
        for seq in mutated_sequence['ACTAGT']:
            assert seq == 'ACTAGT'

        # Do some prints
        u.dp(["Parent sequence:", parent_sequence])
        u.dp(["Mutated sequence:", mutated_sequence])

    def test_making_ssm_with_seq_error(self):
        u.dp(["Testing SSM with error rate"])

        parent_sequence = "ATGAGT"
        mutant = 'ACTAGT'
        parent_sequence_aa = 'MS'
        sequencing_error_rate = 0.5  # Use a
        read_depth = 100
        positions = [0]
        mutated_sequence = make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions,
                                                  parent_sequence_aa, amino_acid_to_codon)
        # Check since we have a 0 seq error rate that we don't have any errors so they are all the expected one
        count_errors = 0
        for seq in mutated_sequence[mutant]:
            for i, pos in enumerate(seq):
                if pos != mutant[i]:
                    count_errors +=1
        # Won't be exact so let's just check it's about 50%
        u.dp([f"Expected error rate: {sequencing_error_rate}", f"Actual error rate: {count_errors/(len(mutant)*read_depth)}"])
        assert count_errors/(len(mutant)*read_depth) > 0.4
        assert count_errors/(len(mutant)*read_depth) < 0.6

        # Change the rate to lower and check again
        sequencing_error_rate = 0.1
        mutated_sequence = make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions,
                                                  parent_sequence_aa, amino_acid_to_codon)
        # Check since we have a 0 seq error rate that we don't have any errors so they are all the expected one
        count_errors = 0
        for seq in mutated_sequence[mutant]:
            for i, pos in enumerate(seq):
                if pos != mutant[i]:
                    count_errors +=1
        # Won't be exact so let's just check it's about 50%
        u.dp([f"Expected error rate: {sequencing_error_rate}", f"Actual error rate: {count_errors/(len(mutant)*read_depth)}"])
        assert count_errors/(len(mutant)*read_depth) > 0
        assert count_errors/(len(mutant)*read_depth) < 0.2

    def test_making_epcr(self):
        u.dp(["Testing ePCR no error"])

        parent_sequence = "ATGAGT"
        mutant = 'ACTAGT'
        sequencing_error_rate = 0
        library_number = 100000  # Do it big so we do actually get the mutant lol still may fail by random chance
        epcr_mutation_rate = 0.8
        read_depth = 100

        mutated_sequence = make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number,
                                                   epcr_mutation_rate)
        assert mutant in list(mutated_sequence.keys())
        avg_err = 0
        for seq in mutated_sequence.keys():
            avg_err += len([c for i, c in enumerate(seq) if parent_sequence[i] != c])/len(seq)
        avg_err = avg_err/len(list(mutated_sequence.keys()))
        u.dp([f"Expected ePCR rate: {epcr_mutation_rate}", f"Actual ePCR rate: {avg_err}"])
        assert avg_err > 0.7
        assert avg_err < 0.9

        # Also then do it with very little error rate and see if it appears
        epcr_mutation_rate = 0
        mutated_sequence = make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence, library_number,
                                                   epcr_mutation_rate)
        assert mutant not in list(mutated_sequence.keys())

        # Then also check the average mutation rate is approx what we expect compared to the parent
        avg_err = 0
        for seq in mutated_sequence.keys():
            avg_err += len([c for i, c in enumerate(seq) if parent_sequence[i] != c])/len(seq)
        avg_err = avg_err/len(list(mutated_sequence.keys()))
        u.dp([f"Expected ePCR rate: {epcr_mutation_rate}", f"Actual ePCR rate: {avg_err}"])
        assert avg_err == 0  # Should be exact for 0

    def test_making_epcr_with_error(self):
        u.dp(["Testing ePCR with error"])

        parent_sequence = "ATGAGT"
        mutant = 'ACTAGT'
        sequencing_error_rate = 0.5
        library_number = 100000  # Do it v big so we do actually get the mutant lol
        epcr_mutation_rate = 0.8
        read_depth = 100

        mutated_sequence = make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence,
                                                   library_number,
                                                   epcr_mutation_rate)
        # Check since we have a 0 seq error rate that we don't have any errors so they are all the expected one
        count_errors = 0
        for seq in mutated_sequence[mutant]:
            for i, pos in enumerate(seq):
                if pos != mutant[i]:
                    count_errors += 1
        # Won't be exact so let's just check it's about 50%
        u.dp([f"Expected error rate: {sequencing_error_rate}",
              f"Actual error rate: {count_errors / (len(mutant) * read_depth)}"])
        assert count_errors / (len(mutant) * read_depth) > 0.4
        assert count_errors / (len(mutant) * read_depth) < 0.6

        # Change the rate to lower and check again
        sequencing_error_rate = 0.1
        mutated_sequence = make_epcr_de_experiment(read_depth, sequencing_error_rate, parent_sequence,
                                                   library_number,
                                                   epcr_mutation_rate)
        # Check since we have a 0 seq error rate that we don't have any errors so they are all the expected one
        count_errors = 0
        for seq in mutated_sequence[mutant]:
            for i, pos in enumerate(seq):
                if pos != mutant[i]:
                    count_errors += 1
        # Won't be exact so let's just check it's about 50%
        u.dp([f"Expected error rate: {sequencing_error_rate}",
              f"Actual error rate: {count_errors / (len(mutant) * read_depth)}"])
        assert count_errors / (len(mutant) * read_depth) > 0
        assert count_errors / (len(mutant) * read_depth) < 0.2

    def test_make_epcr_experiment(self):
        experiment_df = pd.DataFrame()
        read_depth = 20
        number_of_wells = 96
        epcr_mutation_rate = 0.02
        frequency_cutoff = 0.5
        library_number = 96  # Usually do a 96 well plate
        parent_sequence = 'ATGACTCCCTCGGACATCCCGGGATATGATTATGGGCGTGTCGAGAAGTCACCCATCACGGACCTTGAGTTTGACCTTCTGAAGAAGACTGTCATGTTAGGTGAAAAGGACGTAATGTACTTGAAAAAGGCGTGTGACGTTCTGAAAGATCAAGTTGATGAGATCCTTGACTTGGCGGGTGGTTGGGTAGCATCAAATGAGCATTTGATTTATTACTTCTCCAATCCGGATACAGGAGAGCCTATTAAGGAATACCTGGAACGTGTACGCGCTCGCTTTGGAGCCTGGATTCTGGACACTACCTGCCGCGACTATAACCGTGAATGGTTAGACTACCAGTACGAAGTTGGGCTTCGTCATCACCGTTCAAAGAAAGGGGTCACAGACGGAGTACGCACCGTGCCCCATATCCCACTTCGTTATCTTATCGCATGGATCTATCCTATCACCGCCACTATCAAGCCATTTTTGGCTAAGAAAGGTGGCTCTCCGGAAGACATCGAAGGGATGTACAACGCTTGGTTCAAGTCTGTAGTTTTACAAGTTGCCATCTGGTCACACCCTTATACTAAGGAGAATGACTGGCTCGAGCACCACCACCACCACCACTGA'
        for sequencing_error in range(0, 100, 50):
            sequencing_error_rate = sequencing_error / 100.0
            run_df = make_experiment(f'SeqError_{sequencing_error}', read_depth, sequencing_error_rate, parent_sequence,
                                     library_number, number_of_wells, epcr_mutation_rate, frequency_cutoff)
            run_df.reset_index(inplace=True)
            experiment_df = pd.concat([experiment_df, run_df])

    def test_making_well_df_from_reads(self):
        u.dp(["Testing calling variants using SSM with no error"])

        parent_sequence = "ATGAGT"
        mutant = 'ACTAGT'
        parent_sequence_aa = 'MS'
        sequencing_error_rate = 0.0
        read_depth = 100
        positions = [0]
        mutated_sequence = make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions,
                                                  parent_sequence_aa, amino_acid_to_codon)
        parent_name = 'parent'
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
                           'p(g)',
                           'C', 'p(c)', 'N', 'p(n)']
        well_df = calculate_mutation_significance_across_well(well_df)
        label, frequency, combined_p_value, mixed_well = get_variant_label_for_well(well_df, 0.5)
        # This should be mutated at 100% - the rate of our sequencing errror
        u.dp([f"Input parent: {parent_sequence}", f"Variant: {mutant}"])
        u.dp(["label", label, f"frequency", frequency, f"combined_p_value", combined_p_value, "mixed_well", mixed_well, ])

        assert label == 'T2C_G3T'  # Second position has been changed to a C from a T and the third from a G to a T
        assert frequency == 1.0
        assert combined_p_value < 0.05
        assert mixed_well is False

    def test_calling_variant_with_error(self):
        u.dp(["Testing calling variants using SSM with error"])

        parent_sequence = "ATGAGT"
        mutant = 'ACTAGT'
        parent_sequence_aa = 'MS'
        sequencing_error_rate = 0.1
        read_depth = 100
        positions = [0]
        mutated_sequence = make_ssm_de_experiment(read_depth, sequencing_error_rate, parent_sequence, positions,
                                                  parent_sequence_aa, amino_acid_to_codon)
        parent_name = 'parent'
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
                           'p(g)',
                           'C', 'p(c)', 'N', 'p(n)']
        well_df = calculate_mutation_significance_across_well(well_df)
        label, frequency, combined_p_value, mixed_well = get_variant_label_for_well(well_df, 0.5)
        # This should be mutated at 100% - the rate of our sequencing errror
        u.dp([f"Input parent: {parent_sequence}", f"Variant: {mutant}"])
        u.dp(["label", label, f"frequency", frequency, f"combined_p_value", combined_p_value, "mixed_well", mixed_well])

        assert label == 'T2C_G3T'  # Second position has been changed to a C from a T and the third from a G to a T
        assert frequency != 1.0
        assert combined_p_value < 0.05
        assert mixed_well is False

    def test_mixed_wells(self):
        # Test whether we're able to call mixed well populations
        u.dp(["Testing ePCR with error"])

        parent_sequence = 'ATGACTCCCTCGGACATCCCGGGATATGATTATGGGCGTGTCGAGAAGTCACCCATCACGGACCTTGAGTTTGACCTTCTGAAGAAGACTGTCATGTTAGGTGAAAAGGACGTAATGTACTTGAAAAAGGCGTGTGACGTTCTGAAAGATCAAGTTGATGAGATCCTTGACTTGGCGGGTGGTTGGGTAGCATCAAATGAGCATTTGATTTATTACTTCTCCAATCCGGATACAGGAGAGCCTATTAAGGAATACCTGGAACGTGTACGCGCTCGCTTTGGAGCCTGGATTCTGGACACTACCTGCCGCGACTATAACCGTGAATGGTTAGACTACCAGTACGAAGTTGGGCTTCGTCATCACCGTTCAAAGAAAGGGGTCACAGACGGAGTACGCACCGTGCCCCATATCCCACTTCGTTATCTTATCGCATGGATCTATCCTATCACCGCCACTATCAAGCCATTTTTGGCTAAGAAAGGTGGCTCTCCGGAAGACATCGAAGGGATGTACAACGCTTGGTTCAAGTCTGTAGTTTTACAAGTTGCCATCTGGTCACACCCTTATACTAAGGAGAATGACTGGCTCGAGCACCACCACCACCACCACTGA'
        sequencing_error_rate = 0.1
        library_number = 10
        epcr_mutation_rate = 0.02
        read_depth = 10
        number_wells_to_mix = 10
        mixture_rate = 0.5  # 50%
        run_label = 'mutations'
        number_of_wells = 10
        frequency_cutoff = 0.3
        df = make_experiment(run_label, read_depth, sequencing_error_rate, parent_sequence, library_number,
                             number_of_wells, epcr_mutation_rate, frequency_cutoff, number_wells_to_mix, mixture_rate,
                             qc_files_path='/Users/arianemora/Documents/code/MinION/tmp/')
        check_variants(df, parent_sequence)
        df.to_csv('TestMixedWells.csv')


