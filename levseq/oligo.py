import os
import logging
import subprocess
import sys
from pathlib import Path
import pandas as pd
import re
from Bio import SeqIO

class OligoVariantCaller:
    def __init__(self, experiment_name, experiment_folder, ref_sequences_df, filtered_barcodes, 
                 minimap2_path='minimap2', samtools_path='samtools'):
        """
        Initialize OligoVariantCaller with multiple reference sequences
        """
        self.experiment_name = experiment_name
        self.experiment_folder = Path(experiment_folder)
        self.ref_sequences_df = ref_sequences_df
        self.filtered_barcodes = filtered_barcodes
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        
        # Check dependencies
        self._check_dependencies()

    def _check_dependencies(self):
        """Check if required tools are available"""
        missing_deps = []
        try:
            subprocess.run([self.minimap2_path, '--version'], 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            missing_deps.append(f"minimap2 not found at: {self.minimap2_path}")
        
        try:
            subprocess.run([self.samtools_path, '--version'], 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            missing_deps.append(f"samtools not found at: {self.samtools_path}")
        
        if missing_deps:
            for dep in missing_deps:
                logging.error(dep)
            raise FileNotFoundError("Missing required dependencies")

    def prepare_reference_sequences(self, fastq_name, output_dir):
        """
        Create reference FASTA file containing all sequences from CSV
        
        Args:
            fastq_name (str): Name of the FASTQ file being processed (for file naming)
            output_dir (Path): Directory to store reference files
        """
        # Create references directory if it doesn't exist
        references_dir = output_dir / "references"
        references_dir.mkdir(exist_ok=True)
        
        # Create combined FASTA file with all references from the CSV
        combined_fasta = references_dir / f"{fastq_name}_references.fasta"
        with open(combined_fasta, 'w') as f:
            for _, row in self.ref_sequences_df.iterrows():
                f.write(f">{row['id']}\n{row['refseq'].upper()}\n")
        
        # Create minimap2 index
        mmi_path = references_dir / f"{fastq_name}_references.mmi"
        try:
            subprocess.run([
                self.minimap2_path, '-d',
                str(mmi_path),
                str(combined_fasta)
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logging.info(f"Created minimap2 index for {len(self.ref_sequences_df)} references")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error: minimap2 indexing failed - {e.stderr.decode()}")
            raise
        
        return mmi_path

    def align_reads(self, fastq_file, ref_index, output_dir):
        """
        Align reads using minimap2 and store results in the output directory
        """
        try:
            # Create alignments directory
            alignments_dir = output_dir / "alignments"
            alignments_dir.mkdir(exist_ok=True)
            
            # Set up file paths
            fastq_name = Path(fastq_file).stem
            sam_path = alignments_dir / f"{fastq_name}.sam"
            bam_path = alignments_dir / f"{fastq_name}.bam"
            sorted_bam = alignments_dir / f"{fastq_name}.sorted.bam"
            
            # Run minimap2 alignment
            logging.info(f"Aligning {fastq_file} against reference sequences")
            subprocess.run([
                self.minimap2_path,
                '-ax', 'map-ont',
                '-k', '14',
                '-w', '5',
                '--secondary=no',
                str(ref_index),
                str(fastq_file),
                '-o', str(sam_path)
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Convert SAM to BAM and sort
            subprocess.run([
                self.samtools_path, 'view',
                '-bS', sam_path,
                '-o', bam_path
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            subprocess.run([
                self.samtools_path, 'sort',
                bam_path,
                '-o', sorted_bam
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Get alignment statistics
            stats = subprocess.run([
                self.samtools_path, 'flagstat',
                sorted_bam
            ], capture_output=True, text=True, check=True)
            
            # Get detailed alignment information
            alignment_info = subprocess.run([
                self.samtools_path, 'view',
                sorted_bam
            ], capture_output=True, text=True, check=True)
            
            # Process alignment statistics
            total_reads = 0
            for line in stats.stdout.split('\n'):
                if 'total' in line:
                    total_reads = int(line.split()[0])
                    break
            
            # Count alignments for each reference
            ref_alignments = {}
            for line in alignment_info.stdout.split('\n'):
                if line.strip():
                    fields = line.split('\t')
                    if len(fields) >= 3:
                        ref_id = fields[2]
                        if ref_id != '*':  # Skip unmapped reads
                            ref_alignments[ref_id] = ref_alignments.get(ref_id, 0) + 1
            
            # Save detailed alignment results
            results_file = alignments_dir / f"{fastq_name}_alignments.txt"
            with open(results_file, 'w') as f:
                f.write(f"Total reads: {total_reads}\n\n")
                f.write("Alignments by reference:\n")
                for ref_id, count in sorted(ref_alignments.items(), key=lambda x: (-x[1], x[0])):
                    prob = count / total_reads if total_reads > 0 else 0
                    f.write(f"{ref_id}: {count} reads ({prob:.2%})\n")
            
            # Calculate alignment probabilities for each reference
            alignments = []
            for ref_id, count in ref_alignments.items():
                alignments.append({
                    'Reference_ID': ref_id,
                    'total_reads': total_reads,
                    'mapped_reads': count,
                    'alignment_probability': count / total_reads if total_reads > 0 else 0
                })
            
            return alignments
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {fastq_file}: {e.stderr.decode()}")
            return None

    def get_variant_df(self, min_depth=5):
        """Process all FASTQ files in RB-NB structure"""
        results = []
        
        # Find all RB folders
        rb_folders = [f for f in self.experiment_folder.glob("RB*")]
        logging.info(f"Found {len(rb_folders)} RB folders")
        
        for rb_folder in rb_folders:
            nb_folders = [f for f in rb_folder.glob("NB*")]
            
            for nb_folder in nb_folders:
                fastq_files = list(nb_folder.glob("*.fastq"))
                
                if not fastq_files:
                    continue
                
                for fastq_file in fastq_files:
                    try:
                        fastq_name = fastq_file.stem
                        logging.info(f"Processing {fastq_file}")
                        
                        # Prepare references specifically for this FASTQ
                        ref_index = self.prepare_reference_sequences(fastq_name, nb_folder)
                        
                        # Align against all references
                        alignments = self.align_reads(fastq_file, ref_index, nb_folder)
                        
                        if alignments:
                            # Get best alignment
                            best_alignment = max(alignments, key=lambda x: x['alignment_probability'])
                            
                            if best_alignment['total_reads'] >= min_depth:
                                well = self._barcode_to_well(nb_folder.name)
                                ref_row = self.ref_sequences_df[
                                    self.ref_sequences_df['id'] == best_alignment['Reference_ID']
                                ].iloc[0]
                                
                                results.append({
                                    'Well': well,
                                    'Reference_ID': best_alignment['Reference_ID'],
                                    'Reference': ref_row['refseq'],
                                    'Alignment Count': best_alignment['mapped_reads'],
                                    'Match_Frequency': best_alignment['alignment_probability'],
                                    'Total_Reads': best_alignment['total_reads'],
                                    'barcode_plate': self.experiment_name
                                })
                                
                                # Save individual results file
                                result_file = nb_folder / "alignments" / f"{fastq_name}_result.csv"
                                pd.DataFrame([results[-1]]).to_csv(result_file, index=False)
                        
                    except Exception as e:
                        logging.error(f"Error processing {fastq_file}: {str(e)}")
                        continue
        
        if not results:
            logging.warning("No valid alignments found")
            return pd.DataFrame()
        
        # Save complete results
        results_df = pd.DataFrame(results)
        final_results = nb_folder.parent / "final_results.csv"
        results_df.to_csv(final_results, index=False)
        
        return results_df

    @staticmethod
    def _barcode_to_well(barcode):
        """Convert NB barcode to well position"""
        match = re.search(r'\d+', barcode)
        if match:
            number = int(match.group())
            rows = 'ABCDEFGH'
            row = rows[(number - 1) // 12]
            col = (number - 1) % 12 + 1
            return f"{row}{col}"
        return "NA"
