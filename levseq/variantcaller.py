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
import pandas as pd
import logging
from levseq.utils import *
import subprocess
import os
from collections import defaultdict
import glob
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
from Bio import SeqIO
import re
from tqdm import tqdm
import warnings
'''
Script for variant calling

The variant caller starts from demultiplexed fastq files. 

1) Before variant calling, check in the demultiplexed folder if the alignment file exists. If not, return the user that the sequences were not demultiplexed
2) If the files exist, create MSA using minimap2
3) Call variant with soft alignment

'''
# Set up logging with a default level of WARNING
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
# Suppress numpy warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

class VariantCaller:
    """
    Variant caller class.

    """

    def __init__(self, experiment_name, experiment_folder: Path, template_fasta: Path, barcode_path: Path, padding_start: int = 0, padding_end: int = 0) -> None:
        self.barcode_path = barcode_path
        self.experiment_name = experiment_name
        self.experiment_folder = experiment_folder
        self.padding_start = padding_start
        self.padding_end = padding_end
        self.template_fasta = template_fasta
        self.alignment_name = 'alignment_minimap'
        self.variant_dict = {}
        self.ref_name = experiment_name
        self.ref_str = str(SeqIO.read(template_fasta,'fasta').seq)
        self.variant_df = self.build_variant_df_from_barcodes(barcode_path, experiment_name)

    def build_variant_df_from_barcodes(self, barcode_path, experiment_name) -> pd.DataFrame:
        """
        Build variant dataframe from barcodes, forward and reverse barcodes.
        """
        forward_barcode_ids = []
        reverse_barcode_ids = []
        for record in SeqIO.parse(barcode_path, "fasta"):
            if record.id.startswith('NB'):
                forward_barcode_ids.append(record.id)
            elif record.id.startswith('RB'):
                reverse_barcode_ids.append(record.id)
        # Make the dataframe using these and converting them to something more readable (i.e. the name the user assigned
        # to the plate)
        barcode_ids = []
        renamed_ids = []
        plates = []
        wells = []
        self.variant_dict = defaultdict(dict)
        for reverse_barcode in reverse_barcode_ids:
            for forward_barcode in forward_barcode_ids:
                barcode_ids.append(f'{reverse_barcode}_{forward_barcode}')
                well = self._barcode_to_well(forward_barcode)
                plate = experiment_name
                renamed_ids.append(f'{plate}_{well}')
                plates.append(experiment_name)
                wells.append(well)
                self.variant_dict[f'{plate}_{well}'] = {'Plate': experiment_name, 'Well': well,
                                                        'Barcodes': f'{reverse_barcode}_{forward_barcode}',
                                                        'Path': os.path.join(self.experiment_folder, f'{reverse_barcode}/{forward_barcode}')}
        df = pd.DataFrame()
        df['Plate'] = plates
        df['Well'] = wells
        df['Barcode'] = barcode_ids
        df['ID'] = renamed_ids
        return df

    @staticmethod
    def load_reference(reference_path):
        # The reference enables multiple parents to be used for different
        # WARNING: this assumes all the parents are the same
        ref_seq = str(SeqIO.read(template_fasta,'fasta').seq)
        barcode_to_plate_name = experiment_name
        return 'Parent', ref_seq, barcode_to_plate_name

    @staticmethod
    def _barcode_to_well(barcode):
        match = re.search(r'\d+', barcode)
        if match:
            number = int(match.group())
            rows = 'ABCDEFGH'
            row = rows[(number - 1) // 12]
            col = (number - 1) % 12 + 1
            return f"{row}{col}"
        else:
            return "NA"

    def _align_sequences(self, output_dir, filename, scores=[4, 2, 10], alignment_name="alignment_minimap"):
        try:
            all_fastq = os.path.join(output_dir, '*.fastq')
            fastq_list = glob.glob(all_fastq)
            fastq_files = os.path.join(output_dir, f"demultiplexed_{filename}.fastq")

            if not fastq_list:
                logger.error("No FASTQ files found in the specified output directory.")
                return

            # Combining fastq files into one
            with open(fastq_files, 'w') as outfile:
                for fastq in fastq_list:
                    with open(fastq, 'r') as infile:
                        outfile.write(infile.read())

            # Alignment using minimap2
            minimap_cmd = f"minimap2 -ax map-ont -A {scores[0]} -B {scores[1]} -O {scores[2]},24 '{self.template_fasta}' '{fastq_files}' > '{output_dir}/{alignment_name}.sam'"
            subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Convert SAM to BAM and sort
            view_cmd = f"samtools view -bS '{output_dir}/{alignment_name}.sam' > '{output_dir}/{alignment_name}.bam'"
            subprocess.run(view_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            sort_cmd = f"samtools sort '{output_dir}/{alignment_name}.bam' -o '{output_dir}/{alignment_name}.bam'"
            subprocess.run(sort_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Index the BAM file
            index_cmd = f"samtools index '{output_dir}/{alignment_name}.bam'"
            subprocess.run(index_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Cleanup SAM file to save space
            os.remove(f"{output_dir}/{alignment_name}.sam")
        except Exception as e:
            logger.error(f"Error during alignment for {filename}: {e}")

    def _run_variant_thread(self, args):
        barcode_ids, threshold, min_depth, output_dir = args
        # Overall progress bar for all barcodes in this thread
        with tqdm(barcode_ids, desc="Processing barcodes", leave=False) as pbar:
            for barcode_id in pbar:
                try:
                    row = self.variant_dict.get(barcode_id)
                    bam_file = os.path.join(row["Path"], f'{self.alignment_name}.bam')

                    # Check if alignment file exists, if not, align sequences
                    if not os.path.exists(bam_file):
                        logger.info(f"Aligning sequences for {row['Path']}")
                        self._align_sequences(row["Path"], row['Barcodes'])

                    # Placeholder function calls to demonstrate workflow
                    well_df, alignment_count = get_reads_for_well(self.experiment_name, bam_file,
                                                                  self.ref_str, f'{row["Path"]}/msa.fa')
                    self.variant_dict[barcode_id]['Alignment Count'] = alignment_count
                    if well_df is not None:
                        well_df.to_csv(f"{row['Path']}/seq_{barcode_id}.csv", index=False)
                        label, freq, combined_p_value, mixed_well, avg_error_rate = get_variant_label_for_well(well_df, threshold)
                        self.variant_dict[barcode_id]['Variant'] = label
                        self.variant_dict[barcode_id]['Mixed Well'] = mixed_well
                        self.variant_dict[barcode_id]['Average mutation frequency'] = freq
                        self.variant_dict[barcode_id]['P value'] = combined_p_value
                        self.variant_dict[barcode_id]['Average error rate'] = avg_error_rate
                except Exception as e:
                    logger.error(f"Error processing barcode {barcode_id}: {e}")
                finally:
                    pbar.update(1)

    def get_variant_df(self, threshold: float = 0.5, min_depth: int = 5, output_dir='', num_threads=10):
        """
        Get Variant Data Frame for all samples in the experiment

        Args:
            - alignment_file (Path): Path to the alignment file (.bam).
            - qualities (bool, optional): If True, include base qualities in the analysis. Defaults to True.
            - threshold (float, optional): Threshold for calling a variant. Defaults to 0.5.
        """
        self.variant_df['P value'] = float("nan")
        self.variant_df['Mixed Well'] = False
        pool = ThreadPool(num_threads)
        data = []
        num = int(len(self.variant_df) / num_threads)
        self.variant_df.reset_index(inplace=True)
        for i in range(0, len(self.variant_df), num):
            end_i = i + num if i + num < len(self.variant_df) else len(self.variant_df)
            sub_df = self.variant_df.iloc[i: end_i]['ID'].values
            sub_data = [sub_df, threshold, min_depth, output_dir]
            data.append(sub_data)

        # Thread it
        pool.map(self._run_variant_thread, data)

        self.variant_df['Variant'] = [self.variant_dict[b_id].get('Variant') for b_id in self.variant_df['ID'].values]
        self.variant_df['Mixed Well'] = [self.variant_dict[b_id].get('Mixed Well') for b_id in self.variant_df['ID'].values]
        self.variant_df['Average mutation frequency'] = [self.variant_dict[b_id].get('Average mutation frequency') for b_id in self.variant_df['ID'].values]
        self.variant_df['P value'] = [self.variant_dict[b_id].get('P value') if self.variant_dict[b_id].get('P value') else 1.0 for b_id in self.variant_df['ID'].values]
        self.variant_df['Alignment Count'] = [self.variant_dict[b_id].get('Alignment Count') for b_id in self.variant_df['ID'].values]
        self.variant_df['Average error rate'] = [self.variant_dict[b_id].get('Average error rate') for b_id in self.variant_df['ID'].values]

        # Adjust p-values using bonferroni make it simple
        self.variant_df['P adj. value'] = len(self.variant_df) * self.variant_df["P value"].values
        self.variant_df['P adj. value'] = [1 if x > 1 else x for x in self.variant_df["P adj. value"].values]

        return self.variant_df

    def _get_alignment_count(self, sample_folder_path: Path):
        """
        Get the number of alignments in a BAM file.
        
        Args:
            - sample_folder_path (Path): Path to the folder containing the BAM file.
            
        Returns:
            - int: Number of alignments in the BAM file.
        """
        bam_file = os.path.join(sample_folder_path, self.alignment_name)

        if not os.path.exists(bam_file):
            return 0

        try:
            alignment_count = int(
                subprocess.run(f"samtools view -c {bam_file}", shell=True, capture_output=True).stdout.decode(
                    "utf-8").strip())
        except:
            # ToDo: Return a meaningful error here
            print(f'Warning! your bamfile: {bam_file} had no counts! Check the header manually.')
            return 0
        return alignment_count

    def _apply_alignment_count(self):
        """
        Get alignment count for each sample
        """
        for barcode_id, entry in self.variant_dict.items():
            self.variant_dict[barcode_id]["Alignment_count"] = self._get_alignment_count(entry["Path"])
