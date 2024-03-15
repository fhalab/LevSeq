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
import subprocess
import pysam
import os
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
import pandas as pd
import re
from collections import defaultdict
import numpy as np
from tqdm import tqdm
from scipy.stats import binomtest, combine_pvalues
from statsmodels.stats.multitest import multipletests

'''
Script for variant calling

The variant caller starts from demultiplexed fastq files. 

1) Before variant calling, check in the demultiplexed folder if the alignment file exists. If not, return the user that the sequences were not demultiplexed
2) If the files exist, create MSA using minimap2
3) Call variant with soft alignment

'''


def check_demultiplexing(demultiplex_folder: Path, reverse_prefix="RB", forward_prefix="NB", verbose=True):
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
                    continue
                fname = '_'.join(bam_file.split("/")[1:3])
                seq_df = self._get_seq_df([self.ref_name], bam_file, str(self.reference_path),
                                          msa_path=f'{output_dir}msa_{fname}.fa')
                if seq_df is not None:
                    seq_df.to_csv(f'{output_dir}seq_{fname}.csv')
                    # Now use the filter for wells which have a certain threshold of non-reference mutations
                    non_refs = seq_df[seq_df['freq_non_ref'] > threshold].sort_values(by='position')
                    if len(non_refs) > 0:
                        positions = non_refs['position'].values
                        refs = non_refs['ref'].values
                        label = [f'{refs[i]}{positions[i] + 1}{actual}' for i, actual in
                                 enumerate(non_refs['most_frequent'].values)]
                        # Check if it is a mixed well i.e. there were multiple with significant greater than 0.05
                        padj_vals = non_refs[['p(a) adj.', 'p(t) adj.', 'p(g) adj.', 'p(c) adj.', 'p(n) adj.']].values
                        for p in padj_vals:
                            c_sig = 0
                            for padj in p:
                                if padj < 0.05:  # Have this as a variable
                                    c_sig += 1
                            if c_sig > 1:  # potential mixed well
                                self.variant_df.at[i, "Mixed Well"] = True
                        label = '_'.join(label)
                        # Only keep the frequency of the most frequent mutation
                        probability = np.mean([x for x in non_refs['percent_most_freq_mutation'].values])
                        # Combine the values
                        chi2_statistic, combined_p_value = combine_pvalues([x for x in non_refs['p_value adj.'].values],
                                                                           method='fisher')
                    else:
                        label = '#PARENT#'
                        probability = np.mean([1 - x for x in non_refs['freq_non_ref'].values])
                        combined_p_value = float("nan")
                    self.variant_df.at[i, "Variant"] = label
                    self.variant_df.at[i, "Probability"] = probability
                    self.variant_df.at[i, "P value"] = combined_p_value

                else:
                    self.variant_df.at[i, "Variant"] = float("nan")
                    self.variant_df.at[i, "Probability"] = float("nan")

    def get_variant_df(self, qualities=True, threshold: float = 0.5, min_depth: int = 5, output_dir='', num_threads=10):
        """
        Get Variant Data Frame for all samples in the experiment

        Args:
            - alignment_file (Path): Path to the alignment file (.bam).
            - qualities (bool, optional): If True, include base qualities in the analysis. Defaults to True.
            - threshold (float, optional): Threshold for calling a variant. Defaults to 0.2.
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

        # Thread it!
        pool.map(self._run_variant_thread, data)

        # Adjust p-values using bonferroni make it simple
        self.variant_df['P adj. value'] = len(self.variant_df) * self.variant_df["P value"].values
        self.variant_df['P adj. value'] = [1 if x > 1 else x for x in self.variant_df["P adj. value"].values]
        self.variant_df.rename(columns={'Probability': "Average mutation frequency"}, inplace=True)
        # ToDo: make this nicer, but ok for now, postprocess to identify if any mutations appear multiple times
        return self.variant_df

    def postprocess_variant_df(self, df, cutoff=5, output_path=None):
        mutation_map = defaultdict(lambda: defaultdict(int))
        all_plates = pd.DataFrame()
        for plate in set(df['Plate'].values):
            plate_df = df[df['Plate'] == plate]
            positions = []
            for m in plate_df['Variant'].values:
                if str(m) != 'nan':
                    m = m.split('_')
                    for mutation in m:
                        if 'DEL' not in mutation:
                            position = mutation[1:-1] # i.e. trim off what it was
                            # Get the position and also keep what it was mutated to
                            mutation_map[position][mutation[-1]] += 1 # get what it was mutated too
                        else:
                            position = mutation[1:].replace('DEL', '')  # i.e. trim off what it was
                            # Get the position and also keep what it was mutated to
                            mutation_map[position]['DEL'] += 1  # get what it was mutated too
                        positions.append(position)
            # Make into a DF that has positions and then the number and types of muattions
            positions = list(set(positions))
            positions.sort()
            totals, p, a, t, g, c, n = [], [], [], [], [], [], []
            for ps in positions:
                pos = mutation_map[ps]
                a_= pos['A'] or 0
                t_ = pos['T'] or 0
                c_ = pos['C'] or 0
                g_ = pos['G'] or 0
                del_ = pos['DEL'] or 0
                a.append(a_)
                t.append(t_)
                g.append(g_)
                c.append(c_)
                n.append(del_)

                total = a_ + t_ + g_ + c_ + del_
                totals.append(total)
                p.append(ps)
                # CHeck if the well has a problem at a specific position
                if total > cutoff:
                    print(f"Warning! Position {ps} in plate {plate} was mutated: {total} times. This may be an error with your parent.")
            p_df = pd.DataFrame()
            p_df['Position'] = p
            p_df['Total wells mutated in'] = totals
            p_df['A'] = a
            p_df['T'] = t
            p_df['G'] = g
            p_df['C'] = c
            p_df['DEL'] = n
            if output_path:
                # Save QC file if the user specifies a path.
                p_df.to_csv(f'{output_path}Plate_{plate}_QC.csv', index=False)
            p_df['Plate'] = plate
            all_plates = pd.concat([all_plates, p_df])
        all_plates = all_plates.sort_values(by='Total wells mutated in', ascending=False)
        return all_plates # this gives an idea about what might be mutated.

    @staticmethod
    def _alignment_from_cigar(cigar: str, alignment: str, ref: str, query_qualities: list) -> tuple[
        str, str, list, list]:
        """
        Generate the alignment from the cigar string.
        Operation	Description	Consumes query	Consumes reference
        0 M	alignment match (can be a sequence match or mismatch)	yes	yes
        1 I	insertion to the reference	yes	no
        2 D	deletion from the reference	no	yes
        3 N	skipped region from the reference	no	yes
        4 S	soft clipping (clipped sequences present in SEQ)	yes	no
        5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
        6 P	padding (silent deletion from padded reference)	no	no
        7 =	sequence match	yes	yes
        8 X	sequence mismatch	yes	yes

        Parameters
        ----------
        cigar: 49S9M1I24M2D9M2I23M1I3M1D13M1D10M1D13M3D12M1D23M1
        alignment: GCTGATCACAACGAGAGCTCTCGTTGCTCATTACCCCTAAGGAACTCAAATGACGGTTAAAAACTTGTTTTGCT
        ref: reference string (as above but from reference)

        Returns
        -------

        """
        new_seq = ''
        ref_seq = ''
        qual = []
        inserts = []
        pos = 0
        ref_pos = 0
        for op, op_len in cigar:
            if op == 0:  # alignment match (can be a sequence match or mismatch)
                new_seq += alignment[pos:pos + op_len]
                qual += query_qualities[pos:pos + op_len]

                ref_seq += ref[ref_pos:ref_pos + op_len]
                pos += op_len
                ref_pos += op_len
            elif op == 1:  # insertion to the reference
                inserts.append(alignment[pos - 1:pos + op_len])
                pos += op_len
            elif op == 2:  # deletion from the reference
                new_seq += '-' * op_len
                qual += [-1] * op_len
                ref_seq += ref[ref_pos:ref_pos + op_len]
                ref_pos += op_len
            elif op == 3:  # skipped region from the reference
                new_seq += '*' * op_len
                qual += [-2] * op_len
                ref_pos += op_len
            elif op == 4:  # soft clipping (clipped sequences present in SEQ)
                inserts.append(alignment[pos:pos + op_len])
                pos += op_len
            elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
                continue
            elif op == 6:  # padding (silent deletion from padded reference)
                continue
            elif op == 7:  # sequence mismatch
                new_seq += alignment[pos:pos + op_len]
                ref_seq += ref[ref_pos:ref_pos + op_len]
                qual += query_qualities[pos:pos + op_len]
                pos += op_len
                ref_pos += op_len
        return new_seq, ref_seq, qual, inserts

    def _calc_mutation_significance(self, ref_seq, num_a, num_t, num_g, num_c, num_dels, num_reads,
                                    num_total_non_ref_reads, background_error_rate):
        p_a = binomtest(num_a, num_reads, background_error_rate, 'greater').pvalue
        p_t = binomtest(num_t, num_reads, background_error_rate, 'greater').pvalue
        p_g = binomtest(num_g, num_reads, background_error_rate, 'greater').pvalue
        p_c = binomtest(num_c, num_reads, background_error_rate, 'greater').pvalue
        p_n = binomtest(num_dels, num_reads, background_error_rate, 'greater').pvalue
        val = 0
        actual_seq = ref_seq
        p_value = float('nan')  # Could also use 0 not sure what is optimal here!
        if num_total_non_ref_reads == 0:
            val = 0.0  # i.e. they were 100% the reference
            p_value = 1.0  # i.e. they are all this
        else:
            if num_a > 0 and 'A' != ref_seq and num_a / num_reads > val:
                val = num_a / num_reads
                actual_seq = 'A'
                p_value = p_a
            if num_t > 0 and 'T' != ref_seq and num_t / num_reads > val:
                val = num_t / num_reads
                actual_seq = 'T'
                p_value = p_t
            if num_g > 0 and 'G' != ref_seq and num_g / num_reads > val:
                val = num_g / num_reads
                actual_seq = 'G'
                p_value = p_g
            if num_c > 0 and 'C' != ref_seq and num_c / num_reads > val:
                val = num_c / num_reads
                actual_seq = 'C'
                p_value = p_c
            if num_dels > 0 and '-' != ref_seq and num_dels / num_reads > val:
                val = num_dels / num_reads
                actual_seq = 'DEL'
                p_value = p_n
        return actual_seq, val, p_value, p_a, p_t, p_g, p_c, p_n

    def _get_seq_df(self, positions: dict, bam: str, ref: str, min_coverage=20, k=8, msa_path=None):
        """
        Makes a pileup over a gene.

        Rows are the reads, columns are the columns in the reference. Insertions are ignored.
        Parameters
        ----------
        positions: positions
        bam: bam file read in by pysam, pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref: fasta file read in by pysam, pysam.FastaFile(sam/bam file)

        Returns
        -------

        """
        bam = pysam.AlignmentFile(bam, "rb")
        fasta = pysam.FastaFile(ref)
        rows_all = []
        for pos in tqdm(positions):
            reads = []
            try:
                for read in bam.fetch(pos):
                    # Check if we want this read
                    reads.append(read)
            except:
                x = 1

            if len(reads) > min_coverage:
                ref_str = fasta[pos]
                seqs = []
                read_ids = []
                read_quals = []
                for read in reads:
                    if read.query_sequence is not None:
                        seq, ref, qual, ins = self._alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                                         read.query_qualities)
                        # Make it totally align
                        seq = "-" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
                        seqs.append(list(seq))
                        read_ids.append(f'{read.query_name}')
                        read_quals.append(read.qual)

                # Again check that we actually had enough!
                if len(seqs) > min_coverage:
                    seq_df = pd.DataFrame(seqs)
                    # Also add in the read_ids and sort by the quality to only take the highest quality one
                    seq_df['read_id'] = read_ids
                    seq_df['read_qual'] = read_quals
                    seq_df['seqs'] = seqs
                    seq_df = seq_df.sort_values(by='read_qual', ascending=False)
                    # Should now be sorted by the highest quality
                    seq_df = seq_df.drop_duplicates(subset=['read_id'], keep='first')
                    # Reset these guys
                    read_ids = seq_df['read_id'].values
                    seqs = seq_df['seqs'].values

                    seq_df = seq_df.drop(columns=['read_qual', 'read_id', 'seqs'])
                    for col in seq_df:
                        if col - k / 2 > 0:
                            vc = seq_df[col].values
                            ref_seq = ref_str[col]  # Keep track of the reference
                            if ref_seq != '-':
                                # Check if there are at least 25% with a different value compared to the reference.
                                total_reads = len(vc)
                                total_other = len(vc[vc != ref_seq])
                                freq_non_ref = total_other/total_reads
                                a = len(vc[vc == 'A'])
                                t = len(vc[vc == 'T'])
                                g = len(vc[vc == 'G'])
                                c = len(vc[vc == 'C'])
                                nn = len(vc[vc == '-'])
                                val = 0.0
                                actual_seq = ref_seq
                                # Dummy values that will be filled in later once we calculate the background error rate
                                rows_all.append([pos, col, ref_seq, actual_seq, freq_non_ref, total_other, total_reads,
                                                 1.0, val, a, 1.0, t, 1.0, g, 1.0, c, 1.0, nn, 1.0])
                # Check if we want to write a MSA
                if msa_path is not None:
                    with open(msa_path, 'w+') as fout:
                        # Write the reference first
                        fout.write(f'>{self.ref_name}\n{ref_str}\n')

                        for i, seq in enumerate(seqs):
                            fout.write(f'>{read_ids[i]}\n{"".join(seq)}\n')

        bam.close()
        fasta.close()

        if len(rows_all) > 1:  # Check if we have anything to return
            seq_df = pd.DataFrame(rows_all)
            seq_df.columns = ['gene_name', 'position', 'ref', 'most_frequent', 'freq_non_ref', 'total_other',
                              'total_reads', 'p_value', 'percent_most_freq_mutation', 'A', 'p(a)', 'T', 'p(t)', 'G', 'p(g)',
                              'C', 'p(c)', 'N', 'p(n)']
            # Calculate the background error as just the mean frequency of non-reference sequneces (here we "smooth"
            # out the induced mutations.
            mean_error = np.mean(seq_df['freq_non_ref'].values)
            seq_df.reset_index(inplace=True)
            i = 0
            # Using this we can calculate the significance of the different errors
            for ref_seq, num_a, num_t, num_g, num_c, num_dels, num_reads, num_total_non_ref_reads in seq_df[['ref', 'A', 'T', 'G', 'C', 'N', 'total_reads', 'total_other']].values:
                actual_seq, val, p_value, p_a, p_t, p_g, p_c, p_n = self._calc_mutation_significance(ref_seq, num_a,
                                                                                                 num_t, num_g, num_c,
                                                                                                 num_dels, num_reads,
                                                                                                 num_total_non_ref_reads,
                                                                                                 mean_error)
                seq_df.at[i, 'p(a)'] = p_a
                seq_df.at[i, 'p(t)'] = p_t
                seq_df.at[i, 'p(g)'] = p_g
                seq_df.at[i, 'p(c)'] = p_c
                seq_df.at[i, 'p(n)'] = p_n
                seq_df.at[i, 'p_value'] = p_value
                seq_df.at[i, 'percent_most_freq_mutation'] = val
                seq_df.at[i, 'most_frequent'] = actual_seq
                i += 1

            # Do multiple test correction to correct each of the pvalues
            for p in ['p_value', 'p(a)', 'p(t)', 'p(g)', 'p(c)', 'p(n)']:
                # Do B.H which is the simplest
                padjs = multipletests(seq_df[p].values, alpha=0.05, method='fdr_bh')
                seq_df[f'{p} adj.'] = padjs[1]

            return seq_df

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


def get_template_df(plate_numbers: list, barcode_dicts: dict = None, rowwise=True):
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
