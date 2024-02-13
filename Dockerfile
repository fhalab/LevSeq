FROM ubuntu:latest
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
# forgot about this one needed for Drummer https://bedtools.readthedocs.io/en/latest/content/installation.html
RUN apt-get install bedtools
# Clean up mess to make things smaller
RUN apt-get clean

# install conda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
# -c conda-forge r-base
RUN conda create -n rnamod -c bioconda -c conda-forge -c anaconda python=3.6 pysam=0.13 \
    pandas=0.23.4 pybedtools=0.8.0 bedtools=2.25 rpy2=2.8.5 r-base=3.4.1 tqdm=4.40.2 numpy=1.11.3
RUN conda install -c conda-forge r-base

## Activate ELIGOS environment
RUN conda activate rnamod

## Install samplesizeCMH module (Eligos2)
RUN Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")'
## Install Epinano packages
RUN Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("outliers", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("reshape2", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("car", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("ggrepel", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'

# Create new conda env (need to use 3.6 because of Eligos2 and Epinano)
RUN conda create --name rnamod python=3.6 scikit-learn==0.20.2 cloudpickle==1.6.0 fsspec==0.3.3 toolz==0.11.1 pandas==0.25.1 dask==2.5.2
RUN activate rnamod
RUN conda install -c bioconda gffread

# Install all the software
COPY software /software

# This will just make it easy to debug since then we can run stuff from the notebook
RUN pip install notebook

# install pycoQC which is a quality control for basecalled files
# https://a-slide.github.io/pycoQC/installation/
RUN pip install pycoQC

# Install samtools
RUN tar -xvjf /software/htslib-1.15.1.tar.bz2
RUN tar -xvjf /software/bcftools-1.15.1.tar.bz2
RUN tar -xvjf /software/samtools-1.15.1.tar.bz2
WORKDIR /htslib-1.15.1/
RUN make install
WORKDIR /bcftools-1.15.1/
RUN make install
WORKDIR /samtools-1.15.1
RUN make install

# Install differr
RUN pip install /software/differr_nanopore_DRS/dist/differr-0.2.tar.gz

# Change back
WORKDIR /

# Add folder that we'll output data to
# COPY docker_data /docker_data
RUN mkdir /examples

# Add to paths
RUN export PATH="/htslib-1.15.1:/bcftools-1.15.1:/samtools-1.15.1:/software/minimap2:/software/eligos2/Scripts:/software/eligos2:$PATH"

# For some reason it's not wanting to play nice so gonna just do it the ugly way...
RUN cp -r /software/minimap2/* /usr/local/bin
RUN activate rnamod

# Install modrunner remove these two steps
COPY dist/modrunner-1.0.6.tar.gz /
COPY dist/modrunner-1.0.6-py3-none-any.whl /

RUN activate rnamod

RUN pip install typer

RUN conda install -c bioconda gffread
RUN conda install -c bioconda bedtools
RUN pip install -r /software/eligos2_requirements.txt

# Need to install Modrunner
# Set an entry point to CLI for pipeline
COPY modrunner /modrunner
COPY setup.py /
COPY README.md /
COPY LICENSE /
WORKDIR /
RUN python setup.py install
ENTRYPOINT ["modrunner"]