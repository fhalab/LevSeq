from Bio import SeqIO
from Bio.Seq import Seq
import os
from pathlib import Path
import logging
from Bio.Align import PairwiseAligner
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def calculate_alignment_score(seq1, seq2):
    """Calculate alignment score between two sequences using PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignment = aligner.align(seq1, seq2)[0]
    return alignment.score / max(len(seq1), len(seq2))

def filter_single_file(args):
    """
    Filter a single fastq file. Used for parallel processing.
    
    Args:
        args: tuple containing (input_file, parent_seq, parent_rev_comp)
    Returns:
        tuple: (file_path, total_reads, kept_reads, filtered_records)
    """
    input_file, parent_seq, parent_rev_comp = args
    kept_reads = []
    total_reads = 0
    kept_count = 0
    
    is_forward = "forward" in str(input_file).lower()
    
    for record in SeqIO.parse(input_file, "fastq"):
        total_reads += 1
        seq = str(record.seq)
        
        forward_score = calculate_alignment_score(seq, str(parent_seq))
        reverse_score = calculate_alignment_score(seq, str(parent_rev_comp))
        
        # If it's in forward file (plate barcode was rev comp)
        # Then read should align to reverse complement parent sequence
        if is_forward and reverse_score > forward_score:
            kept_reads.append(record)
            kept_count += 1
        # If it's in reverse file (plate barcode was forward)
        # Then read was already reverse complemented by demultiplexer
        # So it should align to forward parent sequence
        elif not is_forward and forward_score > reverse_score:
            kept_reads.append(record)
            kept_count += 1
            
    return str(input_file), total_reads, kept_count, kept_reads

def filter_demultiplexed_folder(experiment_folder, parent_sequence, num_threads=8):
    """
    Filter demultiplexed files using multiple threads.
    
    Args:
        experiment_folder (str): Path to experiment folder containing RBC/FBC structure
        parent_sequence (str): Parent sequence for alignment checking
        num_threads (int): Number of threads to use
    """
    exp_path = Path(experiment_folder)
    filtered_counts = {}
    
    # Prepare parent sequences once
    parent_seq = Seq(parent_sequence)
    parent_rev_comp = parent_seq.reverse_complement()
    
    # Collect all fastq files
    fastq_files = []
    for rbc_dir in exp_path.glob("RB*"):
        if not rbc_dir.is_dir():
            continue
        for fbc_dir in rbc_dir.glob("NB*"):
            if not fbc_dir.is_dir():
                continue
            fastq_files.extend(list(fbc_dir.glob("*.fastq")))
    
    if not fastq_files:
        logging.warning(f"No fastq files found in {experiment_folder}")
        return filtered_counts
    
    # Prepare arguments for parallel processing
    file_args = [(f, parent_seq, parent_rev_comp) for f in fastq_files]
    
    # Process files in parallel with progress bar
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(filter_single_file, args) for args in file_args]
        
        with tqdm(total=len(fastq_files), desc="Filtering files") as pbar:
            for future in as_completed(futures):
                try:
                    file_path, total, kept, filtered_records = future.result()
                    
                    # Write filtered reads
                    temp_file = Path(file_path).parent / f"temp_{Path(file_path).name}"
                    SeqIO.write(filtered_records, temp_file, "fastq")
                    shutil.move(str(temp_file), file_path)
                    
                    filtered_counts[file_path] = {
                        'total': total,
                        'kept': kept,
                        'filtered': total - kept
                    }
                    
                    logging.info(f"Processed {file_path}: {kept}/{total} reads kept")
                    pbar.update(1)
                    
                except Exception as e:
                    logging.error(f"Error processing file {file_path}: {str(e)}")
                    pbar.update(1)
    
    return filtered_counts
