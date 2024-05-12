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
from minION import *
#
# result_path = Path("test_data/")
# experiment_name = "RL-8-70"
# basecall_model_type = "sup"
# result_folder = IO_processor.create_folder(experiment_name,
#                                             basecall_model_type,
#                                             target_path=result_path)
#
# Create Barcode fasta file
barcode_path = "../minION/barcoding/minion_barcodes.fasta" # Path to standard barcode file
front_prefix = "NB"
back_prefix = "RB"
bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
barcode_path = os.path.join("/Users/ariane/Documents/code/MinION/manucript/notebooks/Ape AGW/Data/EVSeqL_Output/YL-EvSeqL1/300-1", "minion_barcodes_filtered.fasta")

# Barcode indexes
front_min = 1
front_max = 96
back_min = 9
back_max = 12

# Expected fragment sizes
min_size = 800
max_size = 5000
#
bp.filter_barcodes(barcode_path, (front_min,front_max), (back_min,back_max))
#
# file_to_experiment = f"test_data/{experiment_name}"
# template_fasta = "data/20220216-ZZ_sup/zz_parent.fasta"
#
# # Basecalling
# basecall_folder = os.path.join(result_folder, "basecalled")
# experiment_folder = ""
#
# # Demultiplexing
# experiment_name = experiment_name + "_" + basecall_model_type
# result_folder_path = IO_processor.find_folder(result_path, experiment_name)

path_to_code = "demultiplex"


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

class TestMinIon(TestClass):

    def test_basecalling(self):
        """
        ToDo test!
        """
        pod5_files = IO_processor.find_folder(experiment_folder, "pod5")
        bc = Basecaller(basecall_model_type, pod5_files, basecall_folder, fastq=True)
        bc.run_basecaller()

    def test_find_fastq(self):
        # Find fastq files
        file_to_fastq = IO_processor.find_folder(experiment_folder, "fastq_pass")
        print(file_to_fastq)

    def test_demultiplex(self):
        # Find fastq files
        file_to_fastq = IO_processor.find_folder(experiment_folder, "fastq_pass")
        print(file_to_fastq)
        prompt = f"{path_to_code} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100} -m {min_size} -x {max_size}"
        subprocess.run(prompt, shell=True)
        demultiplex_folder = result_folder
        print(demultiplex_folder)

    def test_variant_calling(self):
        # To run this will need to update to having stuff not using subprocess...
        vc = VariantCaller(Path('test_data/RL/'),
                           Path(template_fasta),
                           demultiplex_folder_name='',
                           padding_start=0,
                           padding_end=0)

        variant_df = vc.get_variant_df(threshold=0.5,
                                       min_depth=5,
                                       output_dir='test_ZZ/',
                                       num_threads=20)
        variant_df.to_csv('test_ZZ/variant_new_0.5_v6.csv')

    def test_new_minION(self):
        vc = VariantCaller(Path('/Users/ariane/Documents/code/MinION/manucript/notebooks/Ape AGW/Data/EVSeqL_Output/YL-EvSeqL1/300-1/'),
                           Path('/Users/ariane/Documents/code/MinION/manucript/notebooks/Ape AGW/Data/EVSeqL_Output/YL-EvSeqL1/300-1/temp_300-1.fasta'),
                           demultiplex_folder_name='',
                           padding_start=0,
                           padding_end=0)
        variant_df = vc.get_variant_df(threshold=0.5,
                                       min_depth=5,
                                       output_dir='/Users/ariane/Documents/code/MinION/manucript/notebooks/Ape AGW/Data/EVSeqL_Output/YL-EvSeqL1/300-1/msa/',
                                       num_threads=20)
        variant_df.to_csv('/Users/ariane/Documents/code/MinION/manucript/notebooks/Ape AGW/Data/EVSeqL_Output/YL-EvSeqL1/variant_new_0.5_v6.csv')
