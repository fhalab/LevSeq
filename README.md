# Variant Sequencing with Nanopore

In directed evolution, sequencing every variant enhances data insight and creates datasets suitable for AI/ML methods. This method is presented as an extension of the original Every Variant Sequencer using Illumina technology. With this approach, sequence variants can be generated within a day at an extremely low cost.

## Prerequisites

Before using this repository, ensure the following preparations:

- Order forward and reverse primers compatible with the desired plasmid.
- Successfully install Oxford Nanopore's software. [Link to installation guide](#).
- Clone this GitHub repository to your local directory.

## How to Use LevSeq.

The wet lab part is detailed in the method section. Once samples are prepared, the multiplexed sample is used for sequencing, and the sequencing data is stored in the `../data` folder on Nanopore.

After sequencing, variants are identified. For simple applications, we recommend using the notebook `run_levseq.ipynb`.

### Steps:

1. **Basecalling**: This step converts Nanopore's FAST5 files to sequences. For basecalling, we use Nanopore's basecaller, Medaka, which can run in parallel with sequencing (recommended) or afterward.


2. **Demultiplexing**: After sequencing, the reads, stored as bulk FASTQ files, are sorted. During demultiplexing, each read is assigned to its correct plate/well combination and stored as a FASTQ file.

3. **Variant Calling**: For each sample, the consensus sequence is compared to the reference sequence. A variant is called if it differs from the reference sequence. The success of variant calling depends on the number of reads sequenced and their quality.


### Installation:

### Other software
LevSeq uses Samtools and Minimap2 so these need to be installed and added to your path.

Check these work by running:

```
samtools
```

```
minimap2
```

Installation via pip:

```
pip install levseq
```

For installing the whole pipeline, you'll need to use the docker image. For this, install docker as required for your 
operating system (https://docs.docker.com/engine/install/).
```

```
### Development

To build the package via pip you need to run:
```
python setup.py sdist bdist_wheel
```

To then install locally:
```
 pip install dist/levseq-0.1.0.tar.gz
```

To build the docker image run (within the main folder that contains the `Dockerfile`):

```
 docker build -t levseq .
```
**Note** when building on a mac with a MX chip, you need to set your docker env variable to: 
`export DOCKER_DEFAULT_PLATFORM=linux/amd64` or add the flag `--platform linux/amd64` when building.

i.e. you need to run on mac MX:
This gives us the access to the minION command line interface via:

```
docker run levseq
```

```
 docker run -v /Users/arianemora/Documents/projects/LevSeq/data:/levseq_results/ levseq 20240502 levseq_results/20240502/ levseq_results/20240502-YL-ParLQ-ep2.csv
```

levseq --skip_demultiplex oligopools /Users/arianemora/Documents/projects/LevSeq/oligopools /Users/arianemora/Documents/projects/LevSeq/oligopools/oligopool_seqs.csv 


### Steps to rebuild the C++ executables

### Mac intel chip
To rebuild on mac move into the `source/source` folder and execute the following command:

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/local/bin/gcc-13 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-13 ../source
```

Note it expects c to be installed in `/usr/local/bin` if this is not the case on your machine you will need to update 
accordingly. 

After building you need to make the make file with the following command:

```
make -j
```

The demultiplex file should now function!
