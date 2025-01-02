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
import matplotlib.pyplot as plt
from levseq import *

u = SciUtil()


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

"""
Functions for user interactions.
"""

parent = 'ATGCCGCAAATTCCCGGTTACACTTATGGAGATCCCGCCTTACCTCCCAGTCCCGTTTCTTTAGAGGAATTGGAGCGTC\
TGAAGGCTCGCTTGCTGTGGACGGAGGCGGATGATAAAGCATTAGAGCAGGCCGGTAAAGTATTAGAGGACCAAGTAGA\
AGAAGTGTTAGATTTACTGCAGGGCTTTGTCGGGAGCCATCCTCATTTACTTCACTATTTCACTGATCCTCAGGGGAAC\
CCGATCCCCGACTATCTTGAGCGTGTTCGCCGCCGTTTCGGACAGTGGATTCTGGATACCTGCTTTCGTCCCAAAGACG\
AAACCTGGCTTCGTTATCAGCATGAGATTGGCTTACGTCATCATCACACAAAAAAAAACCAAACTGACGGCGTGACCTC\
CGTCCCACATATTCCGTTGCGCTATTTGATCAGTTCTATTTATCCCATTACAGCCACCATCAAACCCTTCCTGACTAAG\
AAGGGCCACAATCCCGAGGAGGTGGAGCGTATGTATCAGGCATGGTTCAAGGCAGTTGTATTGCAAGTAGCACTTTGGT\
CCTATCCTTACACTCAGCCTGGCGACTTTCTCGAGCACCACCACCACCACCACTGA'


class TestUser(TestClass):

    def test_get_seqs(self):
        ### tests making a MSA from the variant DF.
        variant_df = pd.read_csv('test_ZZ/variant_new_0.5.csv')
        seqs, seq_ids = convert_variant_df_to_seqs(variant_df, parent)
        print(seqs, seq_ids)
        assert 'MPQIPGYTYGDPALPPS' in seqs[0]

    def test_msa(self):
        variant_df = pd.read_csv('test_ZZ/variant_new_0.5.csv')
        seqs, seq_ids = convert_variant_df_to_seqs(variant_df, parent)
        # Using the sequences let's now get the MSA
        alignment = make_msa(seqs, seq_ids, 'aln.fa')

    def test_encode(self):
        variant_df = pd.read_csv('test_ZZ/variant_new_0.5.csv')
        seqs, seq_ids = convert_variant_df_to_seqs(variant_df, parent)
        # Using the sequences let's now get the MSA
        alignment = make_msa(seqs, seq_ids, 'aln.fa')
        # Check how one hot encoding of this works
        # Accessing individual records
        encodings = []
        for record in alignment:
            print(record.id, record.seq)
            encoded = one_hot_encode(record.seq)
            encodings.append(np.array(encoded))
        encodings = np.array(encodings)
        # PCA it
        pca = make_pca(encodings)
        plt.scatter(pca[:, 0], pca[:, 1], color='blue', label='Transformed Data')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('PCA Result')
        plt.legend()
        plt.show()

    def test_epcr_pca(self):
        read_depth = 20
        number_of_wells = 96
        epcr_mutation_rate = 0.02
        frequency_cutoff = 0.5
        library_number = 96  # Usually do a 96 well plate
        parent_sequence = 'ATGACTCCCTCGGACATCCCGGGATATGATTATGGGCGTGTCGAGAAGTCACCCATCACGGACCTTGAGTTTGACCTTCTGAAGAAGACTGTCATGTTAGGTGAAAAGGACGTAATGTACTTGAAAAAGGCGTGTGACGTTCTGAAAGATCAAGTTGATGAGATCCTTGACTTGGCGGGTGGTTGGGTAGCATCAAATGAGCATTTGATTTATTACTTCTCCAATCCGGATACAGGAGAGCCTATTAAGGAATACCTGGAACGTGTACGCGCTCGCTTTGGAGCCTGGATTCTGGACACTACCTGCCGCGACTATAACCGTGAATGGTTAGACTACCAGTACGAAGTTGGGCTTCGTCATCACCGTTCAAAGAAAGGGGTCACAGACGGAGTACGCACCGTGCCCCATATCCCACTTCGTTATCTTATCGCATGGATCTATCCTATCACCGCCACTATCAAGCCATTTTTGGCTAAGAAAGGTGGCTCTCCGGAAGACATCGAAGGGATGTACAACGCTTGGTTCAAGTCTGTAGTTTTACAAGTTGCCATCTGGTCACACCCTTATACTAAGGAGAATGACTGGCTCGAGCACCACCACCACCACCACTGA'
        sequencing_error = 30
        sequencing_error_rate = sequencing_error / 100.0
        run_df = make_experiment(f'SeqError_{sequencing_error}', read_depth, sequencing_error_rate, parent_sequence,
                                 library_number, number_of_wells, epcr_mutation_rate, frequency_cutoff)
        run_df.to_csv('run.csv')
        run_df = pd.read_csv('run.csv')
        seqs, seq_ids = convert_variant_df_to_seqs(run_df, parent)
        # Using the sequences let's now get the MSA
        alignment = make_msa(seqs, seq_ids, 'aln.fa')
        # Check how one hot encoding of this works
        # Accessing individual records
        encodings = []
        for record in alignment:
            print(record.id, record.seq)
            encoded = one_hot_encode(record.seq)
            encodings.append(np.array(encoded))
        encodings = np.array(encodings)
        # PCA it
        pca = make_pca(encodings)
        plt.scatter(pca[:, 0], pca[:, 1], color='blue', label='Transformed Data')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('PCA Result')
        plt.legend()
        plt.show()

    def test_convert(self):
        seq = "ATGACTCCCTCGGACATCCCGGGGTATGATTATGGGCGTGTCGAGAAGTCACCCATCACGGACCTTGAGTTTGACCTTCTGAAGAAGACTGTCATGTTAGGTGAAGAGGACATAATGTACTTGAAAAAGGCGGCTGACGTTCTGAAAGATCAAGTTGATGAGATCCTTGACCTGGTTGGTGGTTGGGTAGCATCAAATGAGCATTTGATTTATTACTTCTCCAATCCGGATACAGGAGAGCCTATTAAAGAATACTGTGAACGTATACGCGCTCGCATTGGAGCCTGGGTTCTGGACACTACCTGCCGCGACTATAACCGTGAATGGTTAGACTACCAGTACGAAGTTGGGCTTCGTCATCACCGTTCAAAGAAAGGGGTCACAGACGGAGTACGCACCGTGCCCAATACCCCACTTCGTTATCTTATCGCAGGTATCTATCCTCTTACCGCCACTATCAAGCCACTTTTAGCTGAGAAAGGTGGCTCTCCGGAAGACATCGAAGGGATGTACAACGCTTGGCTCAAGTCTGTAGTTCTACAAGTTGCCATCTGGTCACACCCTTATACTAAGGAGAATGACTGGCTCGAGCACCACCACCACCACCAC"
        print(translate(seq))

