# This test script is used to test the functionality of the data accusition of finished runs from ONT.

import unittest
import os
from levseq.globals import MINKNOW_PATH
from levseq.IO_processor import find_experiment_folder, find_file, get_fbc_barcode_folders, get_rbc_barcode_folders, get_barcode_dict
from levseq.consensus import process_fastq, consensus_prompt, run_medaka, medeka_stitch_prompt
from levseq.util import json_processor import read_json


class TestminION(unittest.TestCase):
    pass

class TestDataAcc(unittest.TestCase):
    def test_find_experiment_folder(self, minknow_path = MINKNOW_PATH):
        """The filename of an experiment """

        data = {"Experiment": "20230905_errorprone-3_test"}

        path = find_experiment_folder(data["Experiment"], minknow_path)

        message = "The path of the experiment folder is not correct."

        self.assertEqual(path, "/var/lib/minknow/data/20230905_errorprone-3_test", message)


    def test_read_json(self):
        """Check if the json file is read correctly and checked correctly"""
        
        json_path = "data/toy_json_files/toy_json.json"
        json_path = os.path.join(os.path.dirname(__file__), json_path)
        data = read_json(json_path)

        self.assertEqual(data["Experiment"]["Name"], "20230905_errorprone-3_test", "The experiment name is not correct")
        self.assertEqual(data["basecaller"]["Barcode-Kit"], "MinION-BARCODES", "The path is not correct")

    def test_find_files(self):
        """Check if the files are found correctly"""
        start_path = os.path.join(os.path.dirname(__file__), "data/toy_ont_folder")
        file_prefix = "final_summary"
        file_extension = ".txt"

        true_path = os.path.join(os.path.dirname(__file__), "data/toy_ont_folder/no_sample/20230607_1644_MN41105_flg114_f6cc1727/final_summary_flg114_f6cc1727_54b00f86.txt")

        file_path = find_file(start_path, file_prefix, file_extension)

        self.assertEqual(true_path, file_path, "The path is not correct")


class TestminION(unittest.TestCase):

    DEMULTIPLEXER_PATH = os.path.join(os.path.dirname(__file__), "data/toy_ont_folder/demultiplex")
    REF_SEQ = os.path.join(os.path.dirname(__file__), "data/refseq/hetcpiii.fasta")

    def setUp(self):
        self.barcode_dict = get_barcode_dict(self.DEMULTIPLEXER_PATH)

class Test_Data_Acc(TestminION):
    """This class test if the data accusition is running correctly. It assumes that levseq is already installed on your computer and the MINKNOW_PATH is set correctly."""

    def test_find_experiment_folder(self):
        """The filename of an experiment """

        data = {"Experiment": "20230905_errorprone-3_test"}

        path = find_experiment_folder(data["Experiment"])

        message = "The path of the experiment folder is not correct."

        self.assertEqual(path, "/var/lib/minknow/data/20230905_errorprone-3_test", message)


class Test_Demultiplexer(TestminION):

    def test_run_demultiplexer(self):
        """Check if the demultiplexer is running correctly"""
        pass

    def test_rbc_barcode_folder(self):
        """Check if the barcode folders are created correctly"""
        barcode_folders = get_rbc_barcode_folders(self.DEMULTIPLEXER_PATH)
        self.assertEqual(len(barcode_folders), 3, "The number of barcode folders is not correct")
    
    def test_fbc_barcode_folder(self):
        """Check if the reverse barcode folders are created correctly"""
        rbc_barcode_folders = get_rbc_barcode_folders(self.DEMULTIPLEXER_PATH)
        fbc_folders = get_fbc_barcode_folders(self.DEMULTIPLEXER_PATH, rbc_barcode_folders)
        self.assertEqual(len(fbc_folders[rbc_barcode_folders[1]]), 96, "The number of barcode folders is not correct")
    
    def test_barcode_dict(self):
        self.assertEqual(len(self.barcode_dict), 3, "The number of barcode folders is not correct")

    def test_process_fastq(self):
        """Check if the fastq files are processed correctly"""
        for reverse in self.barcode_dict.keys():
            for forward in self.barcode_dict[reverse]:
                # Concat all fastq files
                fastq_file = process_fastq(forward)
                self.assertEqual(fastq_file, os.path.join(forward, "pre_consensus.fastq"), "The path of the fastq file is not correct")

    def test_medaka_consensus(self):
        """Check if the medaka consensus is running correctly"""
        test_file = os.path.join(self.DEMULTIPLEXER_PATH, "barcode01/barcode23/pre_consensus.fastq")
        output_dir = os.path.join(os.path.dirname(test_file), "medaka")
        prompt = consensus_prompt(test_file, output_dir, self.REF_SEQ, n_threads=4, model="default")
        run_medaka(prompt)
        
        files_to_check = ["consensus.fasta", "consensus_probs.hdf"]
        self.assertIn(files_to_check[0], os.listdir(output_dir), "The consensus file is not created")
        self.assertIn(files_to_check[1], os.listdir(output_dir), "The consensus probs file is not created")

    def test_medaka_stitch(self):
        """Check if Medaka stitch is running correctly and if the files are created correctly"""
        barcode_folder = os.path.join(self.DEMULTIPLEXER_PATH, "barcode01/barcode23")
        final_consensus = os.path.join(barcode_folder, "final_consensus.fasta")
        prompt = medeka_stitch_prompt(barcode_folder, self.REF_SEQ, final_consensus, qualities=True)
        run_medaka(prompt)
        self.assertTrue(os.path.exists(final_consensus), "The final consensus file is not created")

    def test_run_complete_medaka(self):
        """Check if the whole medaka process is running correctly"""
        for reverse in self.barcode_dict.keys():
            for forward in self.barcode_dict[reverse]:
                # Concat all fastq files
                fastq_file = process_fastq(forward)
                # Output directory
                output_dir = os.path.join(forward, "medaka")
                # Run Consensus
                prompt = consensus_prompt(fastq_file, output_dir, self.REF_SEQ, n_threads=4, model="default")
                run_medaka(prompt)
                # Run Stitch
                final_consensus = os.path.join(forward, "final_consensus.fasta")
                prompt = medeka_stitch_prompt(forward, self.REF_SEQ, final_consensus, qualities=True)
                run_medaka(prompt)
                self.assertTrue(os.path.exists(final_consensus), "The final consensus file is not created")


if __name__ == "__main__":
    unittest.main()