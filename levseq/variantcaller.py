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
import shutil
import os
from collections import defaultdict
import glob
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
from Bio import SeqIO
import re
from tqdm import tqdm
import warnings
import math
'''
Script for variant calling

The variant caller starts from demultiplexed fastq files. 

1) Before variant calling, check in the demultiplexed folder if the alignment file exists. If not, return the user that the sequences were not demultiplexed
2) If the files exist, create MSA using minimap2
3) Call variant with soft alignment

'''

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)  # Set default level for this module
# Use the logger in this file.

class VariantCaller:
    """
    Variant caller class.

    """

    def __init__(self, experiment_name, experiment_folder: Path, template_fasta: Path, barcode_path: Path,
                 padding_start: int = 0, padding_end: int = 0, oligopool=True) -> None:
        self.barcode_path = barcode_path
        self.experiment_name = experiment_name
        self.experiment_folder = experiment_folder
        self.padding_start = padding_start
        self.oligopool = oligopool
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
                if self.oligopool:
                    self.variant_dict[f'{plate}_{well}'] = {'Plate': experiment_name, 'Well': well,
                                                            'Barcodes': f'{reverse_barcode}_{forward_barcode}',
                                                            'Path': os.path.join(self.experiment_folder,
                                                                                 f'{reverse_barcode}/{reverse_barcode}/{forward_barcode}')}
                else:
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
            fastq_list = sorted(glob.glob(all_fastq))
            if not fastq_list:
                logger.error("No FASTQ files found in the specified output directory.")
                return

            sam_path = os.path.join(output_dir, f"{alignment_name}.sam")
            # Alignment using minimap2
            minimap_cmd = [
                "minimap2", "-ax", "map-ont",
                "-A", str(scores[0]),
                "-B", str(scores[1]),
                "-O", f"{scores[2]},24",
                str(self.template_fasta),
                *fastq_list,
            ]
            with open(sam_path, "w") as sam_handle:
                minimap_result = subprocess.run(
                    minimap_cmd,
                    stdout=sam_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            if minimap_result.returncode != 0:
                logger.error(
                    "minimap2 failed for %s: %s",
                    filename,
                    minimap_result.stderr.strip(),
                )
                return
            # print(minimap_cmd)
            # Convert SAM to BAM
            unsorted_bam = os.path.join(output_dir, f"{alignment_name}.unsorted.bam")
            with open(unsorted_bam, "wb") as bam_handle:
                view_result = subprocess.run(
                    ["samtools", "view", "-bS", sam_path],
                    stdout=bam_handle,
                    stderr=subprocess.PIPE,
                )
            if view_result.returncode != 0:
                logger.error(
                    "samtools view failed for %s: %s",
                    filename,
                    view_result.stderr.decode().strip(),
                )
                return

            # Sort BAM (support both modern and legacy samtools syntax)
            sorted_bam = os.path.join(output_dir, f"{alignment_name}.bam")
            sort_result = subprocess.run(
                ["samtools", "sort", "-o", sorted_bam, unsorted_bam],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if sort_result.returncode != 0 or not os.path.exists(sorted_bam):
                legacy_prefix = os.path.join(output_dir, f"{alignment_name}.sorted")
                legacy_result = subprocess.run(
                    ["samtools", "sort", unsorted_bam, legacy_prefix],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                )
                if legacy_result.returncode != 0:
                    logger.error(
                        "samtools sort failed for %s: %s",
                        filename,
                        legacy_result.stderr.decode().strip(),
                    )
                    return
                legacy_bam = f"{legacy_prefix}.bam"
                if not os.path.exists(legacy_bam):
                    logger.error("samtools sort did not produce %s", legacy_bam)
                    return
                shutil.move(legacy_bam, sorted_bam)

            # Index the BAM file
            if not os.path.exists(sorted_bam):
                logger.error("samtools sort did not produce %s", sorted_bam)
                return
            index_cmd = ["samtools", "index", sorted_bam]
            index_result = subprocess.run(
                index_cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if index_result.returncode != 0:
                logger.error(
                    "samtools index failed for %s: %s",
                    filename,
                    index_result.stderr.decode().strip(),
                )
                return

            # Cleanup SAM file to save space
            os.remove(sam_path)
            os.remove(unsorted_bam)
        except Exception as e:
            logger.error(f"Error during alignment for {filename}: {e}")

    def _run_variant_thread(self, args):
        barcode_ids, threshold, min_depth, output_dir = args
        logger.info("Variant calling: processing %d barcodes", len(barcode_ids))
        # Overall progress bar for all barcodes in this thread (disabled to reduce console spam)
        with tqdm(barcode_ids, desc="Processing barcodes", leave=False, disable=True) as pbar:
            for barcode_id in pbar:
                try:
                    row = self.variant_dict.get(barcode_id)
                    bam_file = os.path.join(row["Path"], f'{self.alignment_name}_{barcode_id}.bam')

                    # Check if alignment file exists, if not, align sequences
                    if not os.path.exists(bam_file):
                        logger.info(f"Aligning sequences for {row['Path']}")
                        self._align_sequences(
                            row["Path"],
                            row['Barcodes'],
                            alignment_name=f'{self.alignment_name}_{barcode_id}',
                        )
                    elif not os.path.exists(f"{bam_file}.bai"):
                        subprocess.run(
                            ["samtools", "index", bam_file],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                        )

                    # Placeholder function calls to demonstrate workflow
                    well_df, alignment_count = get_reads_for_well(self.experiment_name, bam_file,
                                                                  self.ref_str, f'{row["Path"]}/{self.alignment_name}_{barcode_id}.fa')
                    if well_df is not None:
                        if self.oligopool:
                            if len(well_df.values) < 10:
                                continue
                        self.variant_dict[barcode_id]['Alignment Count'] = alignment_count
                        well_df.to_csv(f"{row['Path']}/seq_{barcode_id}.csv", index=False)
                        # Suppress noisy numerical warnings from downstream stats on sparse wells.
                        with warnings.catch_warnings():
                            warnings.filterwarnings("ignore", category=RuntimeWarning)
                            label, freq, combined_p_value, mixed_well, avg_error_rate = get_variant_label_for_well(
                                well_df, threshold
                            )
                        self.variant_dict[barcode_id]['Variant'] = label
                        self.variant_dict[barcode_id]['Mixed Well'] = mixed_well
                        self.variant_dict[barcode_id]['Average mutation frequency'] = freq
                        self.variant_dict[barcode_id]['P value'] = combined_p_value
                        self.variant_dict[barcode_id]['Average error rate'] = avg_error_rate
                except Exception as e:
                    logger.error(f"Error processing barcode {barcode_id}: {e}")
                finally:
                    pbar.update(1)

    def get_variant_df(self, threshold: float = 0.5, min_depth: int = 5, output_dir='', num_threads=20):
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
        if num_threads > 1:
            for i in range(0, len(self.variant_df), num):
                end_i = i + num if i + num < len(self.variant_df) else len(self.variant_df)
                sub_df = self.variant_df.iloc[i: end_i]['ID'].values
                sub_data = [sub_df, threshold, min_depth, output_dir]
                data.append(sub_data)

            # Thread it
            pool.map(self._run_variant_thread, data)
        else:
            self._run_variant_thread([self.variant_df, threshold, min_depth, output_dir])

        self.variant_df['Variant'] = [self.variant_dict[b_id].get('Variant') for b_id in self.variant_df['ID'].values]
        self.variant_df['Mixed Well'] = [self.variant_dict[b_id].get('Mixed Well') for b_id in
                                         self.variant_df['ID'].values]
        self.variant_df['Average mutation frequency'] = [self.variant_dict[b_id].get('Average mutation frequency') for
                                                         b_id in self.variant_df['ID'].values]
        self.variant_df['P value'] = [self.variant_dict[b_id].get('P value') for b_id in self.variant_df['ID'].values]
        self.variant_df['Alignment Count'] = [self.variant_dict[b_id].get('Alignment Count') for b_id in
                                              self.variant_df['ID'].values]
        self.variant_df['Average error rate'] = [self.variant_dict[b_id].get('Average error rate') for b_id in
                                                 self.variant_df['ID'].values]
        # Adjust p-values using bonferroni make it simple
        self.variant_df['P adj. value'] = [len(self.variant_df) * p if p else None for p in self.variant_df["P value"].values]
        self.variant_df['P adj. value'] = [1 if x and x > 1 else x for x in self.variant_df["P adj. value"].values]
        if self.oligopool:
            # Filter this so we don't get all the junk
            self.variant_df = self.variant_df[self.variant_df['Alignment Count'] > 2]
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
