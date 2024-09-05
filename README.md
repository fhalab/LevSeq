# Variant Sequencing with Nanopore

In directed evolution, sequencing every variant enhances data insight and creates datasets suitable for AI/ML methods. This method is presented as an extension of the original Every Variant Sequencer using Illumina technology. With this approach, sequence variants can be generated within a day at an extremely low cost.

## Prerequisites

Before using this repository, ensure the following preparations:

- Order forward and reverse primers compatible with the desired plasmid, see methods section of [our paper](http://biorxiv.org/cgi/content/short/2024.09.04.611255v1?rss=1).
- Successfully install Oxford Nanopore's software. [Link to installation guide](https://nanoporetech.com/).
- Clone this GitHub repository to your local directory.
- Data to reproduce the results and to test are available on zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13694463.svg)](https://doi.org/10.5281/zenodo.13694463)

## How to Use LevSeq.

The wet lab part is detailed in the method section of the paper. 

Once samples are prepared, the multiplexed sample is used for sequencing, and the sequencing data is stored in the `../data` folder as per the typical Nanopore flow (refer to Nanopore documentation for this).

After sequencing, you can identify variants, demultiplex, and combine with your variant function here! For simple applications, we recommend using the notebook `example/Example.ipynb`.

### Steps of LevSeq:

1. **Basecalling**: This step converts Nanopore's FAST5 files to sequences. For basecalling, we use Nanopore's basecaller, Medaka, which can run in parallel with sequencing (recommended) or afterward.

2. **Demultiplexing**: After sequencing, the reads, stored as bulk FASTQ files, are sorted. During demultiplexing, each read is assigned to its correct plate/well combination and stored as a FASTQ file.

3. **Variant Calling**: For each sample, the consensus sequence is compared to the reference sequence. A variant is called if it differs from the reference sequence. The success of variant calling depends on the number of reads sequenced and their quality.


### Installation:

We aimed to make LevSeq as simple to use as possible, this means you should be able to run it all using pip. However, if you have issues we recomend using the Docker instance!

#### Dependencies 

1. Samtools: https://www.htslib.org/download/ 
2. Minimap2: https://github.com/lh3/minimap2

Once these are installed, you can install levSeq via pip, available as a release tar.gz file.

Installation via pip:

```
conda create --name levseq python=3.10
```
From the LevSeq folder run:
```
pip install release/levseq-0.1.0.tar.gz
```

For installing the whole pipeline, you'll need to use the docker image. For this, install docker as required for your 
operating system (https://docs.docker.com/engine/install/).

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

This gives us the access to the lebSeq command line interface via:

```
docker run levseq
```
Note! The docker image should work with linux, and mac, however, different mac architectures may have issues (owing to the different M1/M3 processers.)

Basically the -v connects a folder on your computer with the output from the minION sequencer with the docker image that will take these results and then perform 
demultiplexing and variant calling.

```
 docker run -v /Users/XXXX/Documents/LevSeq/data:/levseq_results/ levseq 20240502 levseq_results/20240502/ levseq_results/20240502-YL-ParLQ-ep2.csv
```

In this command: `/Users/XXXX/Documents/LevSeq/data` is a folder on your computer, which contains a subfolder `20240502` 

### Steps to rebuild the C++ executables

This is to run the code locally, rather than via the docker instance, note this will be dependent on your computer and may not 
function as above, we highly recomend using the docker version for this reason.

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

### Issues

If you have any issues, please post an issue! We're trying to make this as user friendly as possible but there may still be issues. 

If you solve something code wise, submit a pull request! We would love community input.