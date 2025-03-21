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

    def test_demultipluxing_pools(self):
        # Take as input the demultiplexed fastq files and the reference csv file
        cl_args = {'skip_demultiplexing': False, 'skip_variantcalling': False, 'threshold': 0.5, 'oligopool': True, 'show_msa': False}
        cl_args["name"] = 'oligotest_21032025'
        cl_args['path'] = 'test_oligopool_2103/'
        cl_args["summary"] = 'test_oligopool_2103/test_oligopool_2103.csv'
        run_LevSeq(cl_args)
