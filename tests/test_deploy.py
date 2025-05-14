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
            f'{os.getcwd()}/test_data/laragen_run:/levseq_results',
            'levseq',
            'test_deploy',
            'levseq_results/levseq-1.2.7/',
            'levseq_results/20241116-LevSeq-Review-Validation-levseq_ref.csv'
        ]
        # docker run --rm -v /Users/arianemora/Documents/code/LevSeq/tests/test_data/laragen_run:/levseq_results levseq test_docker levseq_results/ levseq_results/20241116-LevSeq-Review-Validation-levseq_ref.csv
        print(' '.join(cmd_list))
        # ToDo: add in scoring function for ad4

        # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # print(cmd_return.stdout, cmd_return)

    def test_variant_calling(self):
        # Take as input the demultiplexed fastq files and the reference csv file
        cl_args = {'skip_demultiplexing': True, 'skip_variantcalling': False, 'threshold': 0.5}
        cl_args["name"] = 'test_deploy'
        cl_args['path'] = 'test_data/laragen_run/levseq-1.2.7/'
        cl_args["summary"] = 'test_data/laragen_run/20241116-LevSeq-Review-Validation-levseq_ref.csv'
        variant_df, ref_df = process_ref_csv(cl_args)
        variant_df.to_csv('laragen_test_run.csv')
        # Now we want to check all the variants are the same as in the original case:
        checked_variants_df = pd.read_csv('test_data/laragen_run/levseq-1.2.7/variants_gold_standard.csv')
        checked_variants = checked_variants_df['Variant'].values
        checked_sig = checked_variants_df['Average mutation frequency'].values
        checked_alignments = checked_variants_df['Alignment Count'].values

        i = 0
        for variant, freq, alignment_count, pval in variant_df[['Variant', 'Average mutation frequency',
                                                                'Alignment Count', 'P adj. value']].values:
            print(variant, checked_variants[i])
            if checked_variants[i]:
                if variant:
                    assert variant == checked_variants[i]
                    print(alignment_count, checked_alignments[i])
                    if freq != checked_sig[i]:
                        print(freq, checked_sig[i])
            i += 1


# docker run --rm -v /Users/arianemora/Documents/code/LevSeq/data/degradeo/20250121-JR-IM-HS:/levseq_results levseq 20250121-JR-IM-HS_oligopool levseq_results/ levseq_results/ref_seq_oligopools_single.csv --skip_variantcalling
# levseq oligpool_20250121-JR-IM-HS /Users/arianemora/Documents/code/LevSeq/data/degradeo/20250121-JR-IM-HS/  /Users/arianemora/Documents/code/LevSeq/data/degradeo/20250121-JR-IM-HS/ref_seq_oligopools_all.csv --skip_variantcalling
# levseq results results/  /Users/arianemora/Documents/code/LevSeq/data/degradeo/20250121-JR-IM-HS/ref_seq_oligopools_all.csv --skip_demultiplexing --oligopool