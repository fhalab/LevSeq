FROM ubuntu:latest AS build-demultiplex

# Install system dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    gcc-13 \
    g++-13 \
    git \
    zlib1g-dev \
    build-essential \
    libc6-dev \
    linux-libc-dev \
    && rm -rf /var/lib/apt/lists/*

# Set GCC 13 as the default compiler
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100

WORKDIR /demultiplex

COPY executable/source .

# Use CMake with Release flag and specify the C and C++ compilers
RUN find . -name "CMakeCache.txt" -delete \
    && cmake -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=gcc-13 \
             -DCMAKE_CXX_COMPILER=g++-13 \
             . \
    && make -j


FROM ubuntu:latest AS dependencies
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

# Define the build argument
ARG TARGETARCH

# Install Miniconda depending on the architecture
ENV CONDA_DIR=/opt/conda
RUN if [ "$TARGETARCH" = "arm64" ]; then \
        wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O ~/miniconda.sh; \
    elif [ "$TARGETARCH" = "amd64" ]; then \
        wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    else \
        echo "Unsupported architecture: $TARGETARCH" && exit 1; \
    fi && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH
# -c conda-forge r-base
COPY requirements.txt requirements.txt

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN source /opt/conda/bin/activate
RUN conda init bash
RUN exec bash
RUN conda init bash
RUN source ~/.bashrc
RUN conda create --name levseq python=3.12
# Add levseq to the path
ENV PATH="/opt/conda/envs/levseq/bin:$PATH"
RUN echo "source activate levseq" > ~/.bashrc
RUN source activate levseq && conda install -c conda-forge h5py
RUN source activate levseq && pip install -r requirements.txt
# Install all the software
COPY software /software

# This will just make it easy to debug since then we can run stuff from the notebook
RUN source activate levseq && pip install notebook

# install pycoQC which is a quality control for basecalled files
# https://a-slide.github.io/pycoQC/installation/
RUN source activate levseq && pip install pycoQC

# Download the required software packages using wget
RUN wget -O /software/htslib-1.15.1.tar.bz2 https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
RUN wget -O /software/bcftools-1.15.1.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2
RUN wget -O /software/samtools-1.15.1.tar.bz2 https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
RUN wget -O /software/minimap2-2.28.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2

# Install samtools
RUN tar -xvjf /software/htslib-1.15.1.tar.bz2
RUN tar -xvjf /software/bcftools-1.15.1.tar.bz2
RUN tar -xvjf /software/samtools-1.15.1.tar.bz2
RUN tar -xvjf /software/minimap2-2.28.tar.bz2

# Build and Install
WORKDIR /htslib-1.15.1/
RUN make install
WORKDIR /bcftools-1.15.1/
RUN make install
WORKDIR /samtools-1.15.1
RUN make install
WORKDIR /minimap2-2.28
RUN make || make arm_neon=1 aarch64=1 || make arm_neon=11
# Move the binary to PATH
RUN cp minimap2 /usr/local/bin
RUN chmod +x /usr/local/bin/minimap2

# Verify the installation
RUN minimap2 --version

# Change back
WORKDIR /

# Add folder that we'll output data to
# COPY docker_data /docker_data
RUN mkdir /levseq_results

# Add to paths
RUN export PATH="/htslib-1.15.1:/bcftools-1.15.1:/samtools-1.15.1:/software/minimap2:$PATH"

# Install for demultiplexing
RUN source activate levseq && conda install conda-forge::gcc=13.1
RUN source activate levseq && conda install conda-forge::gxx
RUN source activate levseq && apt install gcc
# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install tzdata package
RUN apt-get update && apt-get install -y tzdata

# Set your timezone
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y libstdc++6


# Set an entry point to CLI for pipeline
COPY levseq /levseq
COPY setup.py /
COPY MANIFEST.in /
COPY README.md /
COPY LICENSE /
RUN mkdir /executable
COPY executable /executable
WORKDIR /

# Copy the binary from build-demultiplex stage
COPY --from=build-demultiplex /demultiplex/bin/demultiplex /levseq/barcoding/demultiplex

# Create the wheel
RUN source activate levseq
RUN python setup.py sdist bdist_wheel

# Install --> should update this to the latest pip version
RUN source activate levseq && pip install dist/levseq-1.2.1.tar.gz

# Add entry point script
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh", "levseq"]
