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
from minION.IO_processor import get_barcode_dict, get_template_sequence
from minION.utils import *
import subprocess
import os
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
import pandas as pd
import re
from tqdm import tqdm

'''
Script for variant calling

The variant caller starts from demultiplexed fastq files. 

1) Before variant calling, check in the demultiplexed folder if the alignment file exists. If not, return the user that the sequences were not demultiplexed
2) If the files exist, create MSA using minimap2
3) Call variant with soft alignment

'''


class VariantCaller:
    """
    Variant caller class.

    """

    def __init__(self, experiment_folder: Path,
                 reference_path: Path,
                 barcodes=True,
                 demultiplex_folder_name="demultiplex",
                 front_barcode_prefix="NB", reverse_barcode_prefix="RB", rowwise=False,
                 padding_start: int = 0, padding_end: int = 0) -> None:
        self.reference_path = reference_path
        self.reference = get_template_sequence(reference_path)
        self.experiment_folder = experiment_folder
        self.ref_name = self._get_ref_name()
        self.padding_start = padding_start
        self.padding_end = padding_end
        self.alignment_name = "alignment_minimap.bam"

        if barcodes:
            self.demultiplex_folder = Path(os.path.join(experiment_folder, demultiplex_folder_name))
            self.barcode_dict = get_barcode_dict(self.demultiplex_folder, front_prefix=front_barcode_prefix,
                                                 reverse_prefix=reverse_barcode_prefix)
            self.variant_df = self._get_sample_name()
            self.variant_df = self._rename_barcodes(rowwise=rowwise, merge=True)
            self.variant_df = self._apply_alignment_count()

    def _get_sample_name(self):
        variant_df = {"Parent": [], "Child": [], "Path": []}

        for _, barcode_dict in self.barcode_dict.items():
            for barcode_path in barcode_dict:
                child = os.path.basename(barcode_path)
                parent = os.path.basename(os.path.dirname(barcode_path))

                variant_df["Parent"].append(parent)
                variant_df["Child"].append(child)
                variant_df["Path"].append(Path(barcode_path))

        return pd.DataFrame(variant_df)

    def _rename_barcodes(self, rowwise=False, parent_name="Plate", child_name="Well", merge=True):

        df = self.variant_df

        if rowwise:
            df = df.rename(columns={'Parent': "Row", 'Child': child_name})
            df["Well"] = df["Well"].apply(self._barcode_to_well)
            df["Plate"] = df['Plate'].str.extract('(\d+)').astype(int)

        else:
            df = df.rename(columns={'Parent': parent_name, 'Child': child_name})
            df["Plate"] = df['Plate'].str.extract('(\d+)').astype(int)
            df["Well"] = df["Well"].apply(self._barcode_to_well)

            # Get Plate Numbers
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
            row = rows[(number - 1) // 12]
            col = (number - 1) % 12 + 1
            return f"{row}{col}"
        else:
            return "NA"

    def _align_sequences(self, output_dir: Path, scores: list = [4, 2, 10], fastq_prefix="demultiplexed",
                         site_saturation: bool = False, alignment_name: str = "alignment_minimap") -> None:
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

            minimap_cmd = f"minimap2 -ax map-ont -A {match_score} -B {mismatch_score} -O {gap_opening_penalty},24 {self.reference_path} {fastq_files_str} > {output_dir}/{alignment_name}.sam"
            subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        else:
            minimap_cmd = f"minimap2 -ax map-ont -A {scores[0]} -B {scores[1]} -O {scores[2]},24 {self.reference_path} {fastq_files_str} > {output_dir}/{alignment_name}.sam"
            subprocess.run(minimap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        view_cmd = f"samtools view -bS {output_dir}/{alignment_name}.sam > {output_dir}/{alignment_name}.bam"
        subprocess.run(view_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        sort_cmd = f"samtools sort {output_dir}/{alignment_name}.bam -o {output_dir}/{alignment_name}.bam"
        subprocess.run(sort_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        index_cmd = f"samtools index {output_dir}/{alignment_name}.bam"
        subprocess.run(index_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Remove SAM file
        os.remove(f"{output_dir}/{alignment_name}.sam")

    def _run_variant_thread(self, args):
        """
        Runs a thread of variant calling.
        """
        data, threshold, min_depth, output_dir = args[0], args[1], args[2], args[3]
        for i, row in tqdm(data.iterrows()):
            if isinstance(row["Path"], os.PathLike):
                bam_file = os.path.join(row["Path"], self.alignment_name)
                # Check if the alignment file exists
                if not os.path.exists(bam_file):
                    # Try aligning the sequences
                    print(f"Aligning sequences for {row['Path']}")
                    self._align_sequences(row["Path"])

                self._apply_alignment_count()

                # Check alignment count
                if self.variant_df["Alignment_count"][i] < min_depth or isinstance(bam_file, float):
                    self.variant_df.at[i, "Variant"] = float("nan")
                    self.variant_df.at[i, "Probability"] = float("nan")
                    print(isinstance(bam_file,float))
                    continue

                fname = '_'.join(bam_file.split("/")[1:3])
                well_df = get_reads_for_well(self.ref_name, bam_file, str(self.reference_path),
                                             msa_path=f'{output_dir}msa_{fname}.fa')
                if well_df is not None:
                    well_df.to_csv(f'{output_dir}seq_{fname}.csv')
                    label, probability, combined_p_value, mixed_well = get_variant_label_for_well(well_df, threshold)
                    self.variant_df.at[i, "Mixed Well"] = mixed_well
                    self.variant_df.at[i, "Variant"] = label
                    self.variant_df.at[i, "Probability"] = probability
                    self.variant_df.at[i, "P value"] = combined_p_value

                else:
                    self.variant_df.at[i, "Variant"] = float("nan")
                    self.variant_df.at[i, "Probability"] = float("nan")

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
            sub_df = self.variant_df.iloc[i: end_i]
            sub_data = [sub_df, threshold, min_depth, output_dir]
            data.append(sub_data)

        # Thread it
        pool.map(self._run_variant_thread, data)

        # Adjust p-values using bonferroni make it simple
        self.variant_df['P adj. value'] = len(self.variant_df) * self.variant_df["P value"].values
        self.variant_df['P adj. value'] = [1 if x > 1 else x for x in self.variant_df["P adj. value"].values]
        self.variant_df.rename(columns={'Probability': "Average mutation frequency"}, inplace=True)

        return self.variant_df

    def _get_ref_name(self):
        with open(self.reference_path, "r") as f:

            # Check if the reference is in fasta format
            ref_name = f.readline().strip()
            if not ref_name.startswith(">"):
                raise ValueError("Reference file is not in fasta format")
            else:
                ref_name = ref_name[1:]
        return ref_name

    def _get_alignment_count(self, sample_folder_path: Path):
        """
        Get the number of alignments in a BAM file.
        
        Args:
            - sample_folder_path (Path): Path to the folder containing the BAM file.
            
        Returns:
            - int: Number of alignments in the BAM file.
        """
        if not isinstance(sample_folder_path, Path):
            return 0

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

        self.variant_df["Path"] = self.variant_df["Path"].apply(lambda x: Path(x) if isinstance(x, str) else x)

        self.variant_df["Alignment_count"] = self.variant_df["Path"].apply(self._get_alignment_count)

        return self.variant_df

