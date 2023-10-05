### 
#
# Quality control for ONT NanoPore Sequencing Results
# Author: Emre Gürsoy
#
###

import os
import gzip
from pathlib import Path
import seaborn as sns
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
from Bio import SeqIO
import pickle
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures as futures


def extract_quality_from_file(file : Path) -> dict:
    """Extract single file zipped or unzipped"""
    quality = {}
    # TODO : Check if the file is zipped or not
    with gzip.open(file, "rt") as fastq_file:

        for record in SeqIO.parse(fastq_file, "fastq"):
            quality[record.id] = record.letter_annotations["phred_quality"]
    
    name = os.path.basename(file)


    return quality




def get_position_quality(experiment_folder : Path, qc_path : Path = None): 
    """Plot the quality of the reads based on their position in the read.

    Args:
        basecalled_path (Path): Path to the basecalled folder.
    """
    

    # basecalled folder
    basecalled_path = experiment_folder / "basecalled_filtered"

    # qc folder
    if qc_path is None:
        qc_path = experiment_folder / "qc"
    
    else:
        qc_path = Path(qc_path)
    
    qc_path.mkdir(parents=True, exist_ok=True)


    if (qc_path / "quality_dist.pkl").exists():

        print("Quality distribution file exists. Loading the file...")

        with open(qc_path / "quality_dist.pkl", "rb") as file:
            quality = pickle.load(file)
            return quality


    gzip_files = list(basecalled_path.glob("*.fastq.gz"))
    quality = {}

    with ThreadPoolExecutor() as executor:

        future_to_file = {executor.submit(extract_quality_from_file, file): file for file in gzip_files}

        for future in futures.as_completed(future_to_file):

            file = future_to_file[future]

            try:
                quality_result = future.result() # Gather results and update the main quality dictionary
                quality.update(quality_result)

            except Exception as e:
                print(f"Error processing {file}: {e}")



    # Save the quality dictionary as pkl file 
                
    filename = qc_path / "quality_dist.pkl"

    with open(filename, "wb") as file:
        pickle.dump(quality, file)


def plot_position_quality(quality_dist):
    """Plot the quality of the reads based on their position in the read.
    
    Args:
        - quality_dist (dict): Dictionary of the quality distribution of the reads.
    
    Returns:
        - SeaBorn lineplot
    """

    fig, ax = plt.subplots(figsize=(20,10))

    # Define x axis as the position of the base
    x = range(0,1000)

    sns.lineplot(x=x, y=quality_dist["read_id"])

    plt.xlabel("Position of the base")
    plt.ylabel("Quality Score")


def summarise_quality_plot(quality : dict):
    """Summarise the quality distribution of the reads."""

    max_length = max(len(v) for v in quality.values())
    means = []
    lower_bounds = []
    upper_bounds = []

    for i in range(max_length):
        
        values_at_i = [v[i] for k, v in quality.items() if len(v) > i]

        # Compute mean
        mean = np.mean(values_at_i)
        means.append(mean)

        # Compute 95% CI
        ci_low, ci_high = stats.t.interval(0.95, len(values_at_i)-1, loc=mean, scale=stats.sem(values_at_i))
        lower_bounds.append(ci_low)
        upper_bounds.append(ci_high)

    return {"mean": means, "lower_bound": lower_bounds, "upper_bound": upper_bounds}

def save_summary(summary : dict, experiment_folder : Path):
    """Save the summary of the quality distribution."""
    filename = experiment_folder / "qc" / "quality_summary.pkl"

    with open(filename, "wb") as file:
        pickle.dump(summary, file)


def plot_summary(summary : dict):
    """Plot qc summary"""

    means = summary["mean"]
    lower_bounds = summary["lower_bound"]
    upper_bounds = summary["upper_bound"]

    plt.figure(figsize=(10, 6))
    sns.lineplot(x=range(1000), y=means)
    plt.fill_between(range(1000), lower_bounds, upper_bounds, color='blue', alpha=0.2)
    plt.title("Mean Quality Score with 95% CI")
    plt.xlabel("Position")
    plt.ylabel("Quality Score")
    plt.savefig("qc/quality_summary.png")
    

if __name__ == "__main__":

    experiment_folder = Path('/home/emre/minION_results/MinION_RBC_0902723_sup')

    # quality = get_position_quality(experiment_folder)

    # summary = summarise_quality_plot(quality)

    # save_summary(summary, experiment_folder)
    
    summary = pickle.load(open(experiment_folder / "qc" / "quality_summary.pkl", "rb"))

    plot_summary(summary)



   