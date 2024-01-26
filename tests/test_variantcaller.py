import unittest
import os
import sys
import glob
from pathlib import Path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from minION.variantcaller import check_demultiplexing, VariantCaller
import pandas as pd

class TestVariantCaller(unittest.TestCase):

    def setUp(self):
        self.experiment_folder = Path("/home/emre/minION_results/20240112-RL-8Plates-FLongle-2_sup/")
        self.reference_path = Path("/home/emre/github_repo/MinION/data/refseq/hetcpiii.fasta")
        self.demultiplex_folder_name = "demultiplexed"
        self.demultiplex_path = self.experiment_folder / self.demultiplex_folder_name
        self.front_barcode_prefix = "NB"
        self.reverse_barcode_prefix = "RB"
        self.variant_caller = VariantCaller(self.experiment_folder, 
                            self.reference_path, 
                            barcodes=True, 
                            demultiplex_folder_name=self.demultiplex_folder_name, 
                            front_barcode_prefix=self.front_barcode_prefix, 
                            reverse_barcode_prefix=self.reverse_barcode_prefix)

    def test_check_demultiplex(self):
        """
        Test if the demultiplexing was successful
        """
   
        parent_folder_count, child_folder_count = check_demultiplexing(self.demultiplex_path, verbose=False)

        expected_parent_folders = 8
        expected_child_folders = 752

        self.assertEqual(parent_folder_count, expected_parent_folders, "Number of parent folders does not match expected value.")
        self.assertEqual(child_folder_count, expected_child_folders, "Number of child folders does not match expected value.")

    def test_initialization(self):
        self.assertEqual(self.variant_caller.experiment_folder, self.experiment_folder)
        self.assertEqual(self.variant_caller.demultiplex_folder, self.demultiplex_path)
        #self.assertEqual(variant_caller.alignment_name, "alignment_minimap.bam")
        #self.assertEqual(variant_caller.depth, 4000)
    
    def test_barcode_df(self):

        data = self.variant_caller._get_sample_name()
        self.assertEqual(len(data["Parent"]), 744)
        self.assertEqual(len(data["Child"]), 744)
        self.assertEqual(len(data["Path"]), 744)

    
    def test_variant_df(self):
        df = self.variant_caller.variant_df

        #Check firt and last row
        self.assertEqual(df.iloc[0]["Plate"], 5)
        self.assertEqual(df.iloc[0]["Well"], "A1")

        self.assertEqual(df.iloc[-1]["Plate"], 12)
        self.assertEqual(df.iloc[-1]["Well"], "H12")

        self.assertEqual(df.shape, (768, 3))

    def test_alignment(self):
        
        path_to_test_folder = Path("/home/emre/minION_results/test_minion")

        #Get folder names
        folder_names = [folder for folder in path_to_test_folder.iterdir() if folder.is_dir()]

        
        for folder in folder_names:
            self.variant_caller._align_sequences(self.reference_path, folder, fastq_prefix = "demulipl")
            #Check if the alignment files are created
            alignment_files = list(folder.glob("*.bam"))
            self.assertEqual(len(alignment_files), 1)


        


if __name__ == "__main__":
    unittest.main()