'''
Script for variant calling
'''


from MinION.minION.util.IO_processor import get_barcode_dict, read_fasta_file, get_template_sequence
import subprocess
import pysam
import os
from collections import Counter
from pathlib import Path
import pandas as pd
import re
import itertools
import numpy as np
from tqdm import tqdm  

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
                 front_barcode_prefix = "NB", reverse_barcode_prefix = "RB", rowwise=False,
                 padding_start : int = 0, padding_end : int = 0) -> None:
        self.reference_path = reference_path
        self.reference = get_template_sequence(reference_path)
        self.experiment_folder = experiment_folder
        self.ref_name = self._get_ref_name()
        self.padding_start = padding_start
        self.padding_end = padding_end
        self.alignment_name = "alignment_minimap.bam"

        if barcodes:
            self.demultiplex_folder = Path(os.path.join(experiment_folder, demultiplex_folder_name))
            self.barcode_dict = get_barcode_dict(self.demultiplex_folder, front_prefix= front_barcode_prefix, reverse_prefix=reverse_barcode_prefix)
            self.variant_df = self._get_sample_name()
            self.variant_df = self._rename_barcodes(rowwise = rowwise, merge=True)
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


    def _rename_barcodes(self, rowwise = False, parent_name = "Plate", child_name = "Well", merge = True):

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


    def _align_sequences(self, output_dir: Path, scores : list = [4,2,10], fastq_prefix = "demultiplexed", site_saturation: bool = False, alignment_name: str = "alignment_minimap") -> None:
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

    def _get_highest_non_ref_base_freq(self, bam_file : Path, positions : list, qualities=True, threshold : float = 0.2):
        """
        The aim of this function is to get the highest probability non-reference base call.
        """
        base_frequencies = {position: Counter() for position in positions}
        base_qualities = {position: [] for position in positions} if qualities else None

        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for pileup_column in bam.pileup(self.ref_name, min(positions) - 1, max(positions), min_base_quality=0,
                                            min_mapping_quality=0, truncate=True):
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
            ref_base = self.reference[position - 1].upper()  
            total_bases = sum(counts.values())
            if total_bases > 0:
                non_ref_bases = {base: count for base, count in counts.items() if base != ref_base}
                if non_ref_bases:
                    max_base = max(non_ref_bases, key=non_ref_bases.get)
                    max_freq = non_ref_bases[max_base] / total_bases
                    if max_freq > threshold:  # Add only if frequency is above threshold
                        highest_non_ref_base_freq[position] = (max_base, max_freq)

                        if qualities:
                            max_base_qualities = [qual for base, qual in zip(counts.elements(), base_qualities[position]) if base == max_base]
                            mean_quality_score[position] = int(sum(max_base_qualities) / len(max_base_qualities)) if max_base_qualities else 0

        return highest_non_ref_base_freq, mean_quality_score if qualities else highest_non_ref_base_freq

    def get_variant_df(self, qualities = True, threshold : float = 0.2, min_depth : int = 5):
        """
        Get Variant Data Frame for all samples in the experiment

        Args:
            - alignment_file (Path): Path to the alignment file (.bam).
            - qualities (bool, optional): If True, include base qualities in the analysis. Defaults to True.
            - threshold (float, optional): Threshold for calling a variant. Defaults to 0.2.
        """

        for i, row in tqdm(self.variant_df.iterrows()):

            try:
            
                bam_file = row["Path"] / self.alignment_name

                #Check if the alignment file exists
                if not bam_file.exists():
                    #Try aligning the sequences
                    print(f"Aligning sequences for {row['Path']}")
                    self._align_sequences(row["Path"])

                self._apply_alignment_count()

                #Check alignment count
                if self.variant_df["Alignment_count"][i] < min_depth or isinstance(bam_file, float):
                    self.variant_df.at[i, "Variant"] = float("nan")
                    self.variant_df.at[i, "Probability"] = float("nan")
                    continue

                variant = self.call_variant(bam_file, qualities=qualities, threshold=threshold)
                self.variant_df.at[i, "Variant"] = variant["Variant"].values
                self.variant_df.at[i, "Probability"] = variant["Probability"].values

            except Exception as e:
                print(e)
                self.variant_df.at[i, "Variant"] = float("nan")
                self.variant_df.at[i, "Probability"] = float("nan")
                continue

        return self.variant_df
    
    def call_variant(self, alignment_file : Path, qualities = True, threshold : float = 0.2, top_N : int = 1):
        """
        Call Variant for a given alignment file

        Args:
            - alignment_file (Path): Path to the alignment file (.bam).
            - qualities (bool, optional): If True, include base qualities in the analysis. Defaults to True. Currently only true works
            - threshold (float, optional): Threshold for calling a variant. Defaults to 0.2.
        """
        variants = {"Variant": [], "Probability": []}
        nb_positions = self._get_postion_range(self.padding_start, self.padding_end)
        freq_dist = self._get_highest_non_ref_base_freq(alignment_file, nb_positions, threshold=threshold, qualities=qualities)[0]
        positions = self._get_nb_positions(freq_dist)

        if not positions: #TODO: Generate 3 random positions before calling parent
            # Choose 3 random positions
            positions = np.random.choice(nb_positions, 3, replace=False)

        elif len(positions) > 10:
            print(f"Too many positions: {len(positions)}, Skipping...")
            #Append nan
            variants["Variant"].append(float("nan"))
            variants["Probability"].append(float("nan"))

            return pd.DataFrame(variants)
       
        pileup_df, qual_df = self._get_bases_from_pileup(alignment_file, positions)
        softmax_df = self._get_softmax_count_df(pileup_df, qual_df, positions)
        variants = self._call_potential_populations(softmax_df, call_threshold=0.1)
        variants = pd.DataFrame(variants).nlargest(top_N, "Probability")


        return variants


    def _get_ref_name(self):
        with open(self.reference_path , "r") as f:

            #Check if the reference is in fasta format
            ref_name = f.readline().strip()
            if not ref_name.startswith(">"):
                raise ValueError("Reference file is not in fasta format")
            else:
                ref_name = ref_name[1:]
        return ref_name

    def _get_alignment_count(self, sample_folder_path : Path):
        """
        Get the number of alignments in a BAM file.
        
        Args:
            - sample_folder_path (Path): Path to the folder containing the BAM file.
            
        Returns:
            - int: Number of alignments in the BAM file.
        """
        if not isinstance(sample_folder_path, Path):
            return 0

        bam_file = sample_folder_path / self.alignment_name

        if not bam_file.exists():
            return 0 

        alignment_count = int(subprocess.run(f"samtools view -c {bam_file}", shell=True, capture_output=True).stdout.decode("utf-8").strip())

        return alignment_count

    def _apply_alignment_count(self):
        """
        Get alignment count for each sample
        """

        self.variant_df["Path"] = self.variant_df["Path"].apply(lambda x: Path(x) if isinstance(x, str) else x)


        self.variant_df["Alignment_count"] = self.variant_df["Path"].apply(self._get_alignment_count)

        return self.variant_df

    def _get_postion_range(self, padding_start : int = 0, padding_end : int = 0):
        """
        Get the positions of the non-reference bases in the reference sequence.

        Args:
            - padding_start (int, optional): Number of bases to pad at the start of the reference sequence. Defaults to 0.
            - padding_end (int, optional): Number of bases to pad at the end of the reference sequence. Defaults to 0.
        
        Returns:
            - list: List of positions of the non-reference bases in the reference sequence.
        """
        
        return range(padding_start + 1, len(self.reference) - padding_end + 1)

    def _get_nb_positions(self, base_dict : dict):
        """
        Get the positions of the non-reference bases in the reference sequence.

        Returns:
            - list: List of positions of the non-reference bases in the reference sequence.
        """
        if not base_dict:
            return []
        else:
            return list(base_dict.keys())

    def _get_bases_from_pileup(self, bam_file, positions, n_neighbours : int = 2):
        all_positions = []  # Include neighbouring positions as well

        for position in positions:
            all_positions.extend(self._get_neighbouring_position(position, neighbour_range=n_neighbours))

        bases_dict = {position: {} for position in all_positions}
        qualities_dict = {position: {} for position in all_positions}
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for pileup_column in bam.pileup(self.ref_name, min(all_positions) - 1, max(all_positions) + 1,
                                            min_base_quality=0,
                                            min_mapping_quality=0,
                                            truncate=True,
                                            max_depth=500):
                pos = pileup_column.pos + 1
                if pos in all_positions:
                    for pileup_read in pileup_column.pileups:
                        read_name = pileup_read.alignment.query_name
                        if pileup_read.is_del:
                            base = '-'
                            quality = 0  # Assign a default quality for deletions
                        elif not pileup_read.is_refskip:
                            base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                            quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                        else:
                            continue

                        if read_name not in bases_dict[pos]:
                            bases_dict[pos][read_name] = base
                            qualities_dict[pos][read_name] = quality

        read_names = sorted(set().union(*[bases_dict[pos].keys() for pos in bases_dict]))
        df_bases = pd.DataFrame(index=read_names, columns=positions)  # Only include original positions
        df_qualities = pd.DataFrame(index=read_names, columns=positions)  # Include all positions

        for pos in all_positions:  # Include all positions
            for read_name in bases_dict[pos]:
                df_bases.at[read_name, pos] = bases_dict[pos][read_name]
                df_qualities.at[read_name, pos] = qualities_dict[pos][read_name]

        #Drop NAs
        df_bases = df_bases.dropna()
        df_qualities = df_qualities.dropna()

        # Calculate average quality for each position based on its neighboring positions
        for pos in positions:
            neighbors = self._get_neighbouring_position(pos)
            df_qualities[pos] = df_qualities[neighbors].mean(axis=1)

        df_qualities = df_qualities[positions]
        df_bases = df_bases[positions]  

        return df_bases, df_qualities
    
    def _get_softmax_count_df(self, bases_df, qual_df, nb_positions):

        alphabet = "ACTG-"

        softmax_counts = {position: [] for position in nb_positions}
        
        for position in nb_positions:
            for base in alphabet:
                base_mask = bases_df[position] == base
                base_counts = base_mask.sum()

                soft_count = sum(base_mask * qual_df[position].apply(VariantCaller._get_non_error_prop))
                softmax_counts[position].append(soft_count)

        softmax_count_df = pd.DataFrame(softmax_counts, columns=nb_positions, index=list(alphabet))

        # Apply softmax to each column (position)
        softmax_count_df = softmax_count_df.apply(lambda x: x / x.sum(), axis=0)

        return softmax_count_df

    def _call_potential_populations(self, softmax_df, call_threshold : float = 0.1):
        """

        """
        positions = softmax_df.columns
        top_combinations = []
        
        # Get the top 2 variants for each position
        for position in positions:
            top_variants = softmax_df[position].nlargest(2)

            if top_variants.iloc[1] < call_threshold:
                top_combinations.append([top_variants.index[0]])
            
            else:
                top_combinations.append(top_variants.index.tolist())

            potential_combinations = list(itertools.product(*top_combinations))

        
        variants = {"Variant" : [], "Probability" : []}
        
        for combination in potential_combinations:
            final_variant = []
            for i, pos in enumerate(positions):

                if combination[i] == self.reference[pos - 1]:
                    continue

                elif combination[i] == "-":
                    var = f"{self.reference[pos - 1]}{pos}DEL"
                    final_variant.append(var)
                else:
                    var = f"{self.reference[pos - 1]}{pos}{combination[i]}"
                    final_variant.append(var)

            final_variant = '_'.join(final_variant)
            if final_variant == "":
                final_variant = "#PARENT#"

            
            joint_prob = VariantCaller.joint_probability_score(softmax_df, combination, positions)
            # TODO: Add calculate score function, so different type of scores can be used:
            # Input: Softmax_df
            # Output: Variant name with Score


            variants["Variant"].append(final_variant)
            variants["Probability"].append(joint_prob)

        return variants

    def _get_neighbouring_position(self, position : int, neighbour_range : int = 2):
        """
        Get the neighbouring positions of the non-reference bases in the reference sequence.

        Args:
            - positions (list): List of positions of the non-reference bases in the reference sequence.
            - neighbour_range (int, optional): Range of neighbouring positions to consider. Defaults to 2.
        
        Returns:
            - list: List of neighbouring positions of the non-reference bases in the reference sequence.
        """

        # Get min range and max range


        pos_range = self._get_postion_range(padding_start=self.padding_start, padding_end=self.padding_end)

        min_range = min(pos_range)
        max_range = max(pos_range)
        
        if position < min_range or position > max_range:
            raise ValueError(f"Position {position} is out of range. The position should be between {min_range} and {max_range}.")

        neighbouring_positions = list(range(position - neighbour_range, position + neighbour_range + 1))

        neighbouring_positions = [pos for pos in neighbouring_positions if pos >= min_range and pos <= max_range]

        return neighbouring_positions
    

    @staticmethod
    def _get_non_error_prop(quality_score):
        """
        Convert quality score to non-error probability.
        
        Args:
            - Phred quality_score (int): Quality score.
        
        Returns:
            - float: Non-error probability.
        """
        return 1 - 10 ** (-quality_score / 10)


    @staticmethod
    def joint_probability_score(softmax_df, combination, positions):
        """
        Calculate the joint probability score for a given combination of variants.

        Args:
            - softmax_df (pd.DataFrame): Softmax count dataframe.
        
        Returns:
            - float: Joint probability score.
        """
        return np.prod([softmax_df.at[combination[i], positions[i]] for i in range(len(positions))])
    
    @staticmethod
    def joint_log_probability_score(softmax_df, combination, positions):
        """
        Calculate the joint probability score for a given combination of variants.

        Args:
            - softmax_df (pd.DataFrame): Softmax count dataframe.
        
        Returns:
            - float: Joint probability score.
        """
        return np.sum([np.log(softmax_df.at[combination[i], positions[i]]) for i in range(len(positions))])



def get_template_df(plate_numbers : list, barcode_dicts : dict = None, rowwise = True):
    """
    To have coherent df for each experiment, a template df is created. The template also have the desired plates and columns in the desired order
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


