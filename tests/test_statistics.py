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


class TestStats(TestClass):

    def test_stats(self):
        # Read in once we have a bam file etc
        ref_str = 'TACAATGCAGCCAACCTGAAGAGCATCAACTTCCATCGCCAGATTGAGGAGAAGTCCCTGGCCCAGCTGAAAAGACAGAT\
        GAAGCGCATCCGTGCCAAACAGGAGAAGTACAGGCAGAGCCAGGCAAGTCGTGGCCAACTCCAGTCCAAAGACCCTCAGG\
        ATCCCAGCCAGGAGCCAGGGCCTGACAGCCCAGGGGGCTCCTCCCCGCCACGGAGACAGTGGTGGCGCCCCTGGCTGGAC\
        CACGCCACAGTCATCCACTCTGGCGACTACTTCCTGTTTGAGTCAGATAGCGAGGAGGAAGAGGAGGCCCTACCTGAGGA\
        CCCCAGGCCTGCAGCTCAGAGTGCCTTCCAGATGGCATACCAGGCATGGGTAACCAATGCCCAGACAGTGCTGAGGCAGC\
        GTCGGGAGCGGGCACGGCAGGAGCGGGCAGAGCAGCTGGCTTCTGGAGGTGACTTGAACCCAGATGTGGAACCAGTAGAT\
        GTCCCAGAAGATGAGATGGCAGGCCGTAGCCACATGATGCAGCGTGTGCTAAGCACCATGCAGTTCCTGTGGGTGCTGGG\
        CGAGACCGGTAAACGCTAGCTACGGCACTGGACAGACGTAGGAGCATCTAGACGCCTCGACTGTGCCTTCTAGTTGCCAG\
        CCATCTG'
        well_df = get_reads_for_well('ref', 'test_data/brian/GCACTGGACAGACGTAGG.bam', ref_str,
                                     msa_path=f'msa.fa')
        well_df.to_csv('test_brian.csv')

