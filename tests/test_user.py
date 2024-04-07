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

from minION.variantcaller import *
from minION.simulation import *
from minION.user import *


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

    def test_msa(self):
        ### tests making a MSA from the variant DF.
        variant_df = pd.read_csv('test_ZZ/variant_new_0.5_v6.csv')
        seqs, seq_ids = convert_variant_df_to_msa(variant_df, parent)
        print(seqs, seq_ids)
