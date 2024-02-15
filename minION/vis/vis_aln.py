"""
A script for visualizing the alignment of reads to a reference genome.
"""

from __future__ import annotations

import os
import pysam

from copy import deepcopy

from Bio import SeqIO

from IPython.display import display

import igv_notebook


class VisAln:
    """
    A class for visualizing the alignment of reads to a reference genome.
    """

    def __init__(self, fasta_path: str, bam_path: str) -> None:

        """
        Args:
        - fasta_path: str, path to the FASTA file.
        - bam_path: str, path to the BAM file.
        """

        self._fasta_path = fasta_path
        self._bam_path = bam_path
        
        igv_notebook.init()
        self._view()

        
    def _process_bam(self) -> tuple[str, str]:

        """
        Process a BAM file to match the sequence names in a FASTA file.

        Args:
        - bam_path: str, path to the BAM file.
        - fasta_seq_name: str, name of the sequence in the FASTA file.

        Returns:
        - str, path to the processed BAM file.
        """

        # Extract the sequence names from the BAM file
        bam_seq_name = extract_bam_seq_name(self._bam_path)

        if len(bam_seq_name) > 1 or bam_seq_name[0] == self.fasta_seq_name:
            return self._bam_path

        else:
            print(f"Processing {self._bam_path} \nto match sequence name with {self._fasta_path}...")
            # Extract the sequence names from the BAM file
            processed_bam_dir = checkNgen_folder(os.path.join(os.path.dirname(self._bam_path), "processed_bam"))
            processed_bam_path = os.path.join(processed_bam_dir, os.path.basename(self._bam_path))
            
            # Open the original BAM file and create a new BAM file with modified header
            with pysam.AlignmentFile(self._bam_path, "rb") as original_bam:

                # Get the original header as a dictionary
                copy_header_dict = deepcopy(original_bam.header.to_dict())
                
                # Modify the header lines to match the sequence names in the FASTA file
                for i, seq_name in enumerate([self.fasta_seq_name]):
                    copy_header_dict["SQ"][i]["SN"] = seq_name

                # Create a new pysam.AlignmentHeader from the modified dictionary
                new_header = pysam.AlignmentHeader.from_dict(copy_header_dict)

                # Create a temporary BAM file with the modified header
                with pysam.AlignmentFile(processed_bam_path, "wb", header=new_header) as processed_bam:
                    # Copy the alignments from the original BAM file to the temporary BAM file
                    for read in original_bam:
                        processed_bam.write(read)
                
            return processed_bam_path

    def _view(self) -> None:

        """
        View the alignment of reads to a reference genome.
        """

        # View the alignment of reads to a reference genome
        b = igv_notebook.Browser(
            {
                "reference": {
                    "id": self.fasta_seq_name,
                    "name": self.fasta_seq_name,
                    "fastaPath": self._fasta_path,
                    "indexPath": self.fasta_fai_path,
                }
            })

        b.load_track(
            {
                "name": self.fasta_seq_name,
                "path": self.processed_bam_path,
                "format": "bam",
                "type": "alignment",
                "sourceType": "file",
            })
    

    @property
    def fasta_seq_name(self) -> str:
        """Return the name of the sequence in the FASTA file."""
        return extract_fasta_seq_name(self._fasta_path)

    @property
    def fasta_fai_path(self) -> str:
        """Return the path to the FASTA index file."""
        return gen_fai(self._fasta_path)
    
    @property
    def bam_seq_name(self) -> list:
        """Return the name of the sequence in the BAM file."""
        return extract_bam_seq_name(self._bam_path)

    @property
    def processed_bam_path(self) -> str:
        """Return the path to the processed BAM file."""
        return self._process_bam()
    


def extract_fasta_seq_name(fasta_path: str) -> str:

    """
    Extract the name of the sequence from a fasta file.

    Args:
    - fasta_path: str, path to the fasta file.

    Returns:
    - str, name of the sequence.
    """

    names = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        names.append(record.id)
    assert len(names) == 1, f"{names} more than one"
    return names[0]


def extract_bam_seq_name(bam_path: str) -> list:

    """
    Extract the name of the sequence from a bam file.

    Args:
    - bam_path: str, path to the bam file.

    Returns:
    - list, name of the sequence.
    """

    sequence_names = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for header_line in bam.header["SQ"]:
            sequence_names.append(header_line.get("SN"))
    
    return sequence_names
    

def gen_fai(fasta_path: str) -> str:

    """
    Generate a FASTA index file.

    Args:
    - fasta_path: str, path to the FASTA file.

    Returns:
    - str, path to the FASTA index file.
    """

    fai_path = fasta_path + ".fai"
    
    if os.path.exists(fai_path):
        os.remove(fai_path)
    
    print(f"Making {fai_path} ...")
    pysam.faidx(fasta_path)
        
    return fai_path


def checkNgen_folder(folder_path: str) -> str:

    """
    Check if the folder and its subfolder exists
    create a new directory if not
    Args:
    - folder_path: str, the folder path
    """

    split_list = os.path.normpath(folder_path).split("/")

    # check if absolute
    if os.path.isabs(folder_path):
        split_list[0] = "/" + split_list[0]

    for p, _ in enumerate(split_list):
        subfolder_path = "/".join(split_list[: p + 1])
        if not os.path.exists(subfolder_path):
            print(f"Making {subfolder_path} ...")
            os.mkdir(subfolder_path)
    return folder_path