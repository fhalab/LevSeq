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
# Import MinION objects
from levseq.run_levseq import *

# Import external packages
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from importlib import resources
from Bio import SeqIO
import tqdm
import os
import shutil

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

    def test_making_pools(self):
        u.dp(["Testing SSM"])
        cl_args = {'skip_demultiplexing': True, 'skip_variantcalling': False}
        cl_args["name"] = 'oligopools'
        cl_args['path'] = '/Users/arianemora/Documents/projects/LevSeq/oligopools/'
        cl_args["summary"] = '/Users/arianemora/Documents/projects/LevSeq/oligopools/oligopool_seqs.csv'
        variant_df = process_ref_csv(cl_args)

        # Check if variants.csv already exist
        variant_csv_path = os.path.join('oligopools', "variants.csv")
        if os.path.exists(variant_csv_path):
            variant_df = pd.read_csv(variant_csv_path)
            df_variants, df_vis = create_df_v(variant_df)
        # Clean up and prepare dataframe for visualization
        else:
            df_variants, df_vis = create_df_v(variant_df)

    def test_pools(self):
        u.dp(["Testing SSM"])
        cl_args = {'skip_demultiplexing': True, 'skip_variantcalling': False}
        cl_args["name"] = 'oligopools'
        cl_args['path'] = '/Users/arianemora/Documents/projects/LevSeq/oligopools/'
        cl_args["summary"] = '/Users/arianemora/Documents/projects/LevSeq/oligopools/oligopool_seqs.csv'
        variant_df = process_ref_csv(cl_args)

        # Check if variants.csv already exist
        variant_csv_path = os.path.join('oligopools', "variants.csv")
        if os.path.exists(variant_csv_path):
            variant_df = pd.read_csv(variant_csv_path)
            df_variants, df_vis = create_df_v(variant_df)
        # Clean up and prepare dataframe for visualization
        else:
            df_variants, df_vis = create_df_v(variant_df)
