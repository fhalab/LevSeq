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
from levseq.run_levseq import process_ref_csv
u = SciUtil()
import math

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

class TestDeploy(TestClass):
    
    def test_deploy(self):
        cmd_list = [
            'docker',  # Needs to be installed as vina.
            'run',
            '--rm',
            '-v',
            f'{os.getcwd()}:/levseq_results',
            'levseq',
            'test_deploy',
            'test_data/laragen_run/levseq-1.2.7/',
            'test_data/laragen_run/20241116-LevSeq-Review-Validation-levseq_ref.csv'
        ]
        # ToDo: add in scoring function for ad4
        cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(cmd_return.stdout, cmd_return)
    
    def test_variant_calling(self):
        # Take as input the demultiplexed fastq files and the reference csv file
        cl_args = {'skip_demultiplexing': True, 'skip_variantcalling': False}
        cl_args["name"] = 'test_deploy'
        cl_args['path'] = 'test_data/laragen_run/levseq-1.2.7/'
        cl_args["summary"] = 'test_data/laragen_run/20241116-LevSeq-Review-Validation-levseq_ref.csv'
        variant_df, ref_df = process_ref_csv(cl_args)
        # Now we want to check all the variants are the same as in the original case:
        checked_variants_df = pd.read_csv('test_data/laragen_run/levseq-1.2.7/variants_gold_standard.csv')
        checked_variants = checked_variants_df['Variant'].values
        checked_sig = checked_variants_df['P adj. value'].values
        i = 0
        for variant, pval in variant_df[['Variant', 'P adj. value']].values:
            print(variant, checked_variants[i])
            if checked_variants[i]:
                if variant:
                    assert variant == checked_variants[i]
            # if pval < 0.05:
            #     assert checked_sig[i] < 0.05
            # elif math.isnan(pval):
            #     assert math.isnan(checked_sig[i])
            # else:
            #     assert checked_sig[i] >= 0.05
            print(pval, checked_sig[i])
            i += 1

