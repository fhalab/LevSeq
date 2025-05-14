import os
import shutil
import unittest
from pathlib import Path
import tempfile
import logging

from levseq.run_levseq import cat_fastq_files

class TestCopyFastq(unittest.TestCase):
    def setUp(self):
        # Create temporary directories for testing
        self.temp_dir = tempfile.mkdtemp()
        self.source_dir = Path(self.temp_dir) / "source"
        self.dest_dir = Path(self.temp_dir) / "dest"
        self.source_dir.mkdir()
        self.dest_dir.mkdir()
        
        # Create a sample fastq file
        self.sample_fastq = self.source_dir / "test.fastq"
        with open(self.sample_fastq, 'w') as f:
            f.write("@SEQ_ID\nGATTACA\n+\n!!!!!!!") 
        
        # Configure logging to only log to console for testing
        logging.basicConfig(level=logging.INFO)
    
    def tearDown(self):
        # Clean up temporary directory
        shutil.rmtree(self.temp_dir)
    
    def test_copy_to_new_location(self):
        """Test copying to a new location works normally"""
        result = cat_fastq_files(self.source_dir, self.dest_dir, reads_per_file=10)
        # The function performs splitting for single .fastq files, so check for split files
        self.assertTrue((self.dest_dir / "test_part0.fastq").exists())
    
    def test_same_location(self):
        """Test copying when source and destination are the same location"""
        # Create a duplicate directory structure to simulate same paths
        same_dir = Path(self.temp_dir) / "same"
        same_dir.mkdir()
        same_fastq = same_dir / "test.fastq"
        with open(same_fastq, 'w') as f:
            f.write("@SEQ_ID\nGATTACA\n+\n!!!!!!!")
        
        # This should not raise an exception
        result = cat_fastq_files(same_dir, same_dir, reads_per_file=10)
        
        # The file should still exist
        self.assertTrue(same_fastq.exists())
        
    def test_compressed_files_same_path(self):
        """Test the case where multiple compressed FASTQ files are in the same path"""
        # Create a duplicate directory with compressed FASTQ files
        compressed_dir = Path(self.temp_dir) / "compressed"
        compressed_dir.mkdir()
        
        # Create two compressed files (we'll fake the compression by just naming them .gz)
        comp_file1 = compressed_dir / "test1.fastq.gz"
        comp_file2 = compressed_dir / "test2.fastq.gz"
        
        with open(comp_file1, 'w') as f:
            f.write("mock compressed file 1") 
        with open(comp_file2, 'w') as f:
            f.write("mock compressed file 2")
            
        # Try to copy to the same location
        result = cat_fastq_files(compressed_dir, compressed_dir, reads_per_file=10)
        
        # Make sure the files still exist and weren't modified
        self.assertTrue(comp_file1.exists())
        self.assertTrue(comp_file2.exists())
        
if __name__ == '__main__':
    unittest.main()