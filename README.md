# Variant Sequencing with Nanopore

In directed evolution, sequencing every variant enhances data insight and creates datasets suitable for AI/ML methods. This method is presented as an extension of the original Every Variant Sequencer using Illumina technology. With this approach, sequence variants can be generated within a day at an extremely low cost.

## Prerequisites

Before using this repository, ensure the following preparations:

- Order forward and reverse primers compatible with the desired plasmid.
- Successfully install Oxford Nanopore's software. [Link to installation guide](#).
- Clone this GitHub repository to your local directory.

## How to Use evSeq-Nanopore

The wet lab part is detailed in the method section. Once samples are prepared, the multiplexed sample is used for sequencing, and the sequencing data is stored in the `../data` folder on Nanopore.

After sequencing, variants are identified. For simple applications, we recommend using the notebook `run_minion.ipynb`.

### Steps:

1. **Basecalling**: This step converts Nanopore's FAST5 files to sequences. For basecalling, we use Nanopore's basecaller, Medaka, which can run in parallel with sequencing (recommended) or afterward.

2. **Demultiplexing**: After sequencing, the reads, stored as bulk FASTQ files, are sorted. During demultiplexing, each read is assigned to its correct plate/well combination and stored as a FASTQ file.

3. **Variant Calling**: For each sample, the consensus sequence is compared to the reference sequence. A variant is called if it differs from the reference sequence. The success of variant calling depends on the number of reads sequenced and their quality.
