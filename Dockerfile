FROM continuumio/anaconda3

# Do the usual things
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y wget
# Need the below to make samtools
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libncurses-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y libcurl4-openssl-dev
# For matplotlib
RUN apt-get install -y --no-install-recommends -f pkg-config
RUN apt-get install -y libfreetype-dev
RUN apt-get install -y libfreetype6
RUN apt-get install -y libfreetype6-dev


# System packages
RUN apt-get update && apt-get install -yq curl wget jq vim

RUN apt-get install bash

# forgot about this one needed for Drummer https://bedtools.readthedocs.io/en/latest/content/installation.html
RUN apt-get install bedtools
# Clean up mess to make things smaller
RUN apt-get clean
# -c conda-forge r-base
COPY requirements.txt requirements.txt
RUN conda create --name levseq python=3.9.18 --y
#RUN echo "source activate levseq" > ~/.bashrc
RUN conda init

RUN activate levseq

RUN pip install -r requirements.txt
# Install all the software
COPY software /software

# This will just make it easy to debug since then we can run stuff from the notebook
RUN pip install notebook

# install pycoQC which is a quality control for basecalled files
# https://a-slide.github.io/pycoQC/installation/
RUN pip install pycoQC

# Install samtools
RUN tar -xvjf /software/htslib-1.20.tar.bz2
RUN tar -xvjf /software/bcftools-1.20.tar.bz2
RUN tar -xvjf /software/samtools-1.20.tar.bz2
WORKDIR /htslib-1.20/
RUN make install
WORKDIR /bcftools-1.20/
RUN make install
WORKDIR /samtools-1.20
RUN make install

# Change back
WORKDIR /

# Add folder that we'll output data to
# COPY docker_data /docker_data
RUN mkdir /levseq_results

# Add to paths
RUN export PATH="/htslib-1.20:/bcftools-1.20:/samtools-1.20:/software/minimap2:$PATH"

# Install for demultiplexing
RUN conda install -y conda-forge::gcc=13.1
RUN conda install -y conda-forge::gxx
RUN apt install -y gcc
# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install tzdata package
RUN apt-get update && apt-get install -y tzdata

# Set your timezone
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install -y libstdc++6

# For some reason it's not wanting to play nice so gonna just do it the ugly way...
RUN cp -r /software/minimap2/* /usr/local/bin

# Install levseq via pip and remove these two steps
COPY dist/levseq-0.1.0.tar.gz /
COPY dist/levseq-0.1.0-py3-none-any.whl /

# Add in some sample data ToDo.!
RUN pip install levseq-0.1.0.tar.gz

# Set an entry point to CLI for pipeline
COPY levseq /levseq
COPY setup.py /
COPY README.md /
COPY LICENSE /
RUN mkdir /source
COPY source /source
RUN apt install g++ build-essential
WORKDIR /
RUN python setup.py install
ENTRYPOINT ["levseq"]