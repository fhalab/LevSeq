'''
Script for variant calling
'''


from minION.util.IO_processor import get_barcode_dict
from minION import analyser
from minION import consensus
import subprocess
import pysam
import os
from collections import Counter
from pathlib import Path
from typing import Union
import pandas as pd
import re

'''
The variant caller starts from demultiplexed fastq files. 

1) Before variant calling, check in the demultiplexed folder if the alignment file exists. If not, return the user that the sequences were not demultiplexed
2) If the files exist, create MSA using minimap2
3) Call variant with soft alignment

'''


def check_demultiplexing(demultiplex_folder : Path, reverse_prefix = "RB", forward_prefix = "NB", verbose=True):
    """
    Check if the demultiplexing was done correctly. If not, return the user that the sequences were not demultiplexed.

    Args:
        - demultiplex_folder: Path to the folder containing the demultiplexed fastq files
        - verbose: If True, print the name of each parent folder and the count of child folders

    Return:
        - Tuple: Number of parent folders and child folders
    """
    demultiplex_path = Path(demultiplex_folder)
    parent_folder_count = 0
    child_folder_count = 0

    for child in demultiplex_path.iterdir():
        if child.is_dir() and (child.name.startswith(reverse_prefix) or child.name.startswith(forward_prefix)):
            parent_folder_count += 1
            child_folder_count += len(list(child.iterdir()))
            if verbose:
                print(f"Parent folder '{child.name}' contains {len(list(child.iterdir()))} folders.")

    return parent_folder_count, child_folder_count


    
class VariantCaller:

    """
    Variant caller class.

    """

    def __init__(self, experiment_folder : Path, 
                 reference_path : Path,
                 barcodes = True, 
                 demultiplex_folder_name = "demultiplex", 
                 front_barcode_prefix = "NB", reverse_barcode_prefix = "RB", rowwise=False) -> None:
        
        self.reference = analyser.get_template_sequence(reference_path)
        self.experiment_folder = experiment_folder

        if barcodes:
            self.demultiplex_folder = experiment_folder / demultiplex_folder_name
            self.barcode_dict = get_barcode_dict(self.demultiplex_folder, front_prefix= front_barcode_prefix, reverse_prefix=reverse_barcode_prefix)
            self.variant_df = self._get_sample_name()
            self.variant_df = self._rename_barcode(rowwise = rowwise,merge=True)

        self.alignment_name = "alignment_minimap.bam"
        self.reference = self.experiment_folder / "reference.fasta"
        self.depth = 4000
        self.allele_frequency = 0.1

    def _get_sample_name(self):
        variant_df = {"Parent": [], "Child": [], "Path": []}

        for _, barcode_dict in self.barcode_dict.items():
            for barcode_path in barcode_dict:
                
                child = os.path.basename(barcode_path)
                parent = os.path.basename(os.path.dirname(barcode_path))

                variant_df["Parent"].append(parent)
                variant_df["Child"].append(child)
                variant_df["Path"].append(barcode_path)

        return pd.DataFrame(variant_df)


    def _rename_barcode(self, rowwise = False, parent_name = "Plate", child_name = "Well", merge = True):

        df = self.variant_df

        if rowwise:
            df = df.rename(columns={'Parent': "Row", 'Child': child_name})
            df["Well"] = df["Well"].apply(self._barcode_to_well)
            df["Plate"] = df['Plate'].str.extract('(\d+)').astype(int)

        else:
            df = df.rename(columns={'Parent': parent_name, 'Child': child_name})
            df["Plate"] = df['Plate'].str.extract('(\d+)').astype(int)
            df["Well"] = df["Well"].apply(self._barcode_to_well)

            #Get Plate Numbers
            plate_numbers = df[parent_name].unique()
            plate_numbers.sort()

            if merge:
                template_df = get_template_df(plate_numbers, self.barcode_dict, rowwise=rowwise)
                df = pd.merge(template_df, df, on=[parent_name, child_name], how="left")

        return df

 
    
    @staticmethod
    def _barcode_to_well(barcode):
        match = re.search(r'\d+', barcode)
        if match:
            number = int(match.group())
            rows = 'ABCDEFGH'
            row = rows[(number-1) // 12]
            col = (number-1) % 12 + 1
            return f"{row}{col}"
        else:
            return "NA"


    def _align_sequences(self, ref: Path, output_dir: Path, scores : list = [4,2,10], fastq_prefix = "demultiplexed", site_saturation: bool = False, alignment_name: str = "alignment_minimap") -> None:
        """
        Aligns sequences using minimap2, converts to BAM, sorts and indexes the BAM file.

        Args:
            - ref (Path): Path to the reference file.
            - output_dir (str or Path): Directory to store output files.
            - scores (list, optional): List of match, mismatch and gap opening scores. Defaults to [4,2,10].
            - site_saturation (bool, optional): If True, uses site saturation parameters for minimap2. Defaults to False.
            - alignment_name (str, optional): Name of the alignment file. Defaults to "alignment_minimap".

        Returns:
            - None
        """

        fastq_files = output_dir.glob(f"{fastq_prefix}*.fastq")

        if not fastq_files:
            raise FileNotFoundError("No FASTQ files found in the specified output directory.")
        
        fastq_files_str = " ".join(str(file) for file in fastq_files)

        if site_saturation:
            alignment_name = "alignment_minimap_site_saturation"

            match_score = 4
            mismatch_score = 2
            gap_opening_penalty = 10

            minimap_cmd = f"minimap2 -ax map-ont -A {match_score} -B {mismatch_score} -O {gap_opening_penalty},24 {ref} {fastq_files_str} > {output_dir}/{alignment_name}.sam"
            subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        else:
            minimap_cmd = f"minimap2 -ax map-ont -A {scores[0]} -B {scores[1]} -O {scores[2]},24 {ref} {fastq_files_str} > {output_dir}/{alignment_name}.sam"
            subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        view_cmd = f"samtools view -bS {output_dir}/{alignment_name}.sam > {output_dir}/{alignment_name}.bam"
        subprocess.run(view_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        sort_cmd = f"samtools sort {output_dir}/{alignment_name}.bam -o {output_dir}/{alignment_name}.bam"
        subprocess.run(sort_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        index_cmd = f"samtools index {output_dir}/{alignment_name}.bam"
        subprocess.run(index_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Remove SAM file
        os.remove(f"{output_dir}/{alignment_name}.sam")

    def _get_highest_non_ref_base_freq(self, bam_file : Path, ref_name : str, positions : list, reference_sequence, qualities=True):

        base_frequencies = {position: Counter() for position in positions}
        base_qualities = {position: [] for position in positions} if qualities else None

        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for pileup_column in bam.pileup(ref_name, min(positions) - 1, max(positions),
                                            min_base_quality=0, 
                                            min_mapping_quality=0, 
                                            truncate=True):
                if pileup_column.pos + 1 in positions:
                    for pileup_read in pileup_column.pileups:
                        if not pileup_read.is_del and not pileup_read.is_refskip:
                            base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                            base_frequencies[pileup_column.pos + 1].update([base])
                            if qualities:
                                base_quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                                base_qualities[pileup_column.pos + 1].append(base_quality)
                        elif pileup_read.is_del:
                            base_frequencies[pileup_column.pos + 1].update(["-"])
                            if qualities:
                                base_qualities[pileup_column.pos + 1].append(0)

        highest_non_ref_base_freq = {}
        mean_quality_score = {position: 0 for position in positions} if qualities else None
        for position in positions:
            counts = base_frequencies[position]
            ref_base = reference_sequence[position - 1].upper()  
            total_bases = sum(counts.values())
            if total_bases > 0:
                non_ref_bases = {base: count for base, count in counts.items() if base != ref_base}
                if non_ref_bases:
                    max_base = max(non_ref_bases, key=non_ref_bases.get)
                    max_freq = non_ref_bases[max_base] / total_bases
                    highest_non_ref_base_freq[position] = (max_base, max_freq)

                    if qualities:
                        max_base_qualities = [qual for base, qual in zip(counts.elements(), base_qualities[position]) if base == max_base]
                        mean_quality_score[position] = int(sum(max_base_qualities) / len(max_base_qualities)) if max_base_qualities else 0
                else:
                    highest_non_ref_base_freq[position] = (None, 0)
            else:
                highest_non_ref_base_freq[position] = (None, 0)

        return highest_non_ref_base_freq, mean_quality_score if qualities else highest_non_ref_base_freq

    
    





def get_template_df(plate_numbers : list, barcode_dicts : dict = None, rowwise = True):
    """To have coherent df for each experiment, a template df is created. The template also have the desired plates and columns in the desired order
    Input:
        - demultiplex_folder, folder where the demultiplexed files are located
        - rowwise, if True, the reverse barcodes are rows and not plates
        """
    
    if barcode_dicts is None:
        raise ValueError("No barcode dictionary provided")
    
    n_rbc = len(barcode_dicts.items())

    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    columns = [i for i in range(1, 13)]

    if rowwise:
        template = {"Plate": [], "Well": []}

        for row in rows:
            for column in columns:
                template["Plate"].append(1)
                template["Well"].append(f"{row}{column}")
    
    else:

        template = {"Plate": [], "Well": []}

        for i in range(n_rbc):
            for row in rows:
                for column in columns:
                    template["Plate"].append(plate_numbers[i])
                    template["Well"].append(f"{row}{column}")
        
    return pd.DataFrame(template)
