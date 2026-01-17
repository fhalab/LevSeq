from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import gzip
import logging
import math
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

_VALID_BASES = set("ACGT")


def build_kmer_set(seq, kmer_size):
    """Build a set of k-mers from the parent sequence."""
    if kmer_size <= 0:
        return set()
    seq = seq.upper()
    if len(seq) < kmer_size:
        return set()
    kmers = set()
    for i in range(len(seq) - kmer_size + 1):
        kmer = seq[i:i + kmer_size]
        if set(kmer) <= _VALID_BASES:
            kmers.add(kmer)
    return kmers


def sample_kmer_positions(seq_len, kmer_size, samples, skip_front, skip_back):
    start_min = max(skip_front, 0)
    start_max = seq_len - max(skip_back, 0) - kmer_size
    if start_max < start_min or samples <= 0:
        return []
    span = start_max - start_min
    if span == 0:
        return [start_min]
    if samples == 1:
        return [start_min + span // 2]
    step = span / (samples - 1)
    positions = [int(round(start_min + i * step)) for i in range(samples)]
    seen = set()
    uniq_positions = []
    for pos in positions:
        if pos < start_min or pos > start_max:
            continue
        if pos in seen:
            continue
        seen.add(pos)
        uniq_positions.append(pos)
    return uniq_positions


def count_kmer_hits(seq, positions, kmer_size, parent_kmers, rev_kmers):
    forward_hits = 0
    reverse_hits = 0
    sampled = 0
    for pos in positions:
        kmer = seq[pos:pos + kmer_size]
        if len(kmer) != kmer_size:
            continue
        if set(kmer) <= _VALID_BASES:
            sampled += 1
            if kmer in parent_kmers:
                forward_hits += 1
            if kmer in rev_kmers:
                reverse_hits += 1
    return forward_hits, reverse_hits, sampled


def iter_fastq_records(handle):
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not seq or not plus or not qual:
            break
        yield header, seq, plus, qual

def filter_single_file(args):
    """
    Filter a single fastq file. Used for parallel processing.
    
    Args:
        args: tuple containing (input_file, parent_kmers, rev_kmers, kmer_size, samples,
                                skip_front, skip_back, min_delta, min_ratio)
    Returns:
        tuple: (file_path, total_reads, kept_reads, temp_file)
    """
    (input_file, parent_kmers, rev_kmers, kmer_size, samples,
     skip_front, skip_back, min_delta, min_ratio) = args
    total_reads = 0
    kept_count = 0
    
    is_forward = "forward" in str(input_file).lower()

    input_path = Path(input_file)
    temp_file = input_path.parent / f"temp_{input_path.name}"
    position_cache = {}
    open_fn = gzip.open if input_path.suffix == ".gz" else open
    with open_fn(input_path, "rt") as input_handle, open(temp_file, "w") as output_handle:
        for header, seq_line, plus, qual in iter_fastq_records(input_handle):
            total_reads += 1
            seq = seq_line.strip().upper()
            seq_len = len(seq)
            if seq_len not in position_cache:
                position_cache[seq_len] = sample_kmer_positions(
                    seq_len, kmer_size, samples, skip_front, skip_back
                )
            positions = position_cache[seq_len]
            forward_hits, reverse_hits, sampled = count_kmer_hits(
                seq, positions, kmer_size, parent_kmers, rev_kmers
            )
            if sampled == 0:
                continue
            required_delta = max(min_delta, int(math.ceil(min_ratio * sampled)))
            if required_delta > sampled:
                required_delta = sampled

            # If it's in forward file (plate barcode was rev comp)
            # Then read should align to reverse complement parent sequence
            if is_forward and (reverse_hits - forward_hits) >= required_delta:
                output_handle.write(header)
                output_handle.write(seq_line)
                output_handle.write(plus)
                output_handle.write(qual)
                kept_count += 1
            # If it's in reverse file (plate barcode was forward)
            # Then read was already reverse complemented by demultiplexer
            # So it should align to forward parent sequence
            elif not is_forward and (forward_hits - reverse_hits) >= required_delta:
                output_handle.write(header)
                output_handle.write(seq_line)
                output_handle.write(plus)
                output_handle.write(qual)
                kept_count += 1

    return str(input_file), total_reads, kept_count, str(temp_file)

def filter_demultiplexed_folder(
    experiment_folder,
    parent_sequence,
    num_threads=8,
    kmer_size=6,
    samples=40,
    skip_front=100,
    skip_back=0,
    min_delta=4,
    min_ratio=0.1,
):
    """
    Filter demultiplexed files using a k-mer orientation heuristic.
    
    Args:
        experiment_folder (str): Path to experiment folder containing RBC/FBC structure
        parent_sequence (str): Parent sequence for alignment checking
        num_threads (int): Number of threads to use
        kmer_size (int): Length of k-mer used for orientation checks
        samples (int): Number of k-mers sampled per read
        skip_front (int): Bases to skip from the front of the read
        skip_back (int): Bases to skip from the end of the read
        min_delta (int): Minimum hit difference to keep a read
        min_ratio (float): Minimum hit difference as a ratio of sampled k-mers
    """
    exp_path = Path(experiment_folder)
    filtered_counts = {}
    
    # Prepare parent sequences once
    parent_seq_obj = Seq(parent_sequence)
    parent_seq = str(parent_seq_obj).upper()
    parent_rev_comp = str(parent_seq_obj.reverse_complement()).upper()
    parent_kmers = build_kmer_set(parent_seq, kmer_size)
    rev_kmers = build_kmer_set(parent_rev_comp, kmer_size)
    
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
    file_args = [
        (f, parent_kmers, rev_kmers, kmer_size, samples, skip_front, skip_back, min_delta, min_ratio)
        for f in fastq_files
    ]
    
    # Process files in parallel with progress bar
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(filter_single_file, args) for args in file_args]
        
        with tqdm(total=len(fastq_files), desc="Filtering files") as pbar:
            for future in as_completed(futures):
                try:
                    file_path, total, kept, temp_file = future.result()

                    shutil.move(temp_file, file_path)
                    
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
